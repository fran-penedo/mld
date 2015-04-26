import gurobipy as g
from stl import *

def label(name, i, j):
    return name + "_" + str(i) + "_" + str(j)


def delta_label(l, i):
    return l + "_" + str(i)


def add_max_constr(m, label, args, K, nnegative=True):
    return add_minmax_constr(m, label, args, K, 'max', nnegative)


def add_min_constr(m, label, args, K, nnegative=True):
    return add_minmax_constr(m, label, args, K, 'min', nnegative)


# TODO handle len(args) == 2 differently
def add_minmax_constr(m, label, args, K, op='min', nnegative=True):
    if op not in ['min', 'max']:
        raise ValueError('Expected one of [min, max]')

    y = {}
    y[label] = m.addVar(lb=0 if nnegative else -K, ub=K, name=label)

    if len(args) == 0:
        m.update()
        m.addConstr(y[label] == K)

    else:
        for i in range(len(args)):
            l = delta_label(label, i)
            y[l] = m.addVar(vtype=g.GRB.BINARY, name=l)

        m.update()

        for i in range(len(args)):
            x = delta_label(label, i)
            if op == 'min':
                m.addConstr(y[label] <= args[i])
                m.addConstr(y[label] >= args[i] - y[x] * K)
            else:
                m.addConstr(y[label] >= args[i])
                m.addConstr(y[label] <= args[i] + y[x] * K)

        m.addConstr(g.quicksum([y[delta_label(label, i)]
                                for i in range(len(args))]) == len(args) - 1)

    return y


def add_abs_var(m, label, var, obj):
    if m.modelSense == g.GRB.MINIMIZE:
        z = m.addVar(name="abs_"+label, obj=obj)
        m.update()
        m.addConstr(z >= var)
        m.addConstr(z >= -var)
        return z

    else:
        return None


class Signal(object):

    def __init__(self, signal, bounds):
        self.signal = signal
        self.bounds = bounds


def _stl_expr(m, label, f, t):
    bounds = f.args[0].bounds
    y = m.addVar(name=label, lb=bounds[0], ub=bounds[1])
    m.update()
    m.addConstr(y == f.args[0].signal[t])
    return y, bounds


def _stl_not(m, label, f, t):
    x, bounds = _stl_expr(m, label + "_not", f.args[0], t)
    y = m.addVar(name=label, lb=bounds[0], ub=bounds[1])
    m.update()
    m.addConstr(y == -x)


def _stl_and_or(m, label, f, t, op):
    xx = []
    boundss = []
    for i, ff in enumerate(f.args):
        x, bounds = add_stl_constr(m, label + "_" + op + str(i), ff, t)
        xx.append(x)
        boundss.append(bounds)

    # I'm not gonna bother using the best bounds
    bounds = map(max, zip(*boundss))
    K = max(map(abs, bounds))
    add = add_min_constr if op == "min" else add_max_constr
    y = add(m, label, xx, K, nnegative=False)[label]
    return y, bounds


def _stl_and(m, label, f, t):
    return _stl_and_or(m, label, f, t, "and")


def _stl_or(m, label, f, t):
    return _stl_and_or(m, label, f, t, "or")


def _stl_next(m, label, f, t):
    return add_stl_constr(m, label, f.args[0], t+1)


def _stl_always_eventually(m, label, f, t, op):
    xx = []
    boundss = []
    for i in range(f.bnd[0], f.bnd[1] + 1):
        x, bounds = add_stl_constr(m, label + "_" + op + str(i), f.args[0],
                                   t + i)
        xx.append(x)
        boundss.append(bounds)

    # I'm not gonna bother using the best bounds
    bounds = map(max, zip(*boundss))
    K = max(map(abs, bounds))
    add = add_min_constr if op == "alw" else add_max_constr
    y = add(m, label, xx, K, nnegative=False)[label]
    return y, bounds


def _stl_always(m, label, f, t):
    return _stl_always_eventually(m, label, f, t, "alw")


def _stl_eventually(m, label, f, t):
    return _stl_always_eventually(m, label, f, t, "eve")


def add_stl_constr(m, label, f, t=0):
    return {
        EXPR: _stl_expr,
        NOT: _stl_not,
        AND: _stl_and,
        OR: _stl_or,
        ALWAYS: _stl_always,
        NEXT: _stl_next,
        EVENTUALLY: _stl_eventually
    }[f.op](m, label, f, t)


def add_penalty(m, label, var, obj):
    m.setAttr("OBJ", [var], [-obj])
    y = add_abs_var(m, label, var, obj)
    return y


def add_always_constr(m, label, a, b, rho, K, t=0):
    y = add_min_constr(m, label, rho[(t + a):(t + b + 1)], K)[label]
    return y


def add_always_penalized(m, label, a, b, rho, K, obj, t=0):
    y = add_always_constr(m, label, a, b, rho, K, t)
    add_penalty(m, label, y, obj)
    return y


def create_model(A, b, B, E, d, c, xcap, x0, cost, N):
    # Big number for min alpha/beta (xcap - x)
    K_out = max(xcap) / min([min([1] + [abs(i) for i in row if i > 0])
                             for row in B])
    # Big number required for y = min(x, c, min)
    # No need to add bound on x because it's included in K_out
    K = max(max(c), K_out) * 10
    # Bounds for x
    mm = 0
    M = max(xcap)

    m = g.Model("traffic")
    m.modelSense = g.GRB.MINIMIZE

    x = {}
    y = {}
    z = {}
    u = {}
    v = {}
    for i in range(len(B)):
        for j in range(N):
            # actual decision variables
            # system state
            labelx = label("x", i, j)
            # auxiliary binary variable z = y * u
            labelz = label("z", i, j)
            # control variable
            labelu = label("u", i, j)

            x[labelx] = m.addVar(
                obj=cost[i], lb=0, ub=g.GRB.INFINITY, name=labelx)
            if j < N - 1:
                z[labelz] = m.addVar(lb=0, ub=g.GRB.INFINITY, name=labelz)
                u[labelu] = m.addVar(vtype=g.GRB.BINARY, name=labelu)

    m.update()

    for i in range(len(B)):
        for j in range(N):
            # v = min(alpha/beta(xcap - x))
            v.update(
                add_min_constr(
                    m, label("v", i, j),
                    [E[k][i] / B[k][i] * (xcap[i] - x[label("x", i, j)])
                     for k in range(len(B)) if B[k][i] > 0], K_out))
            # y = min(x, c, v)
            y.update(
                add_min_constr(
                    m, label("y", i, j),
                    [x[label("x", i, j)], c[i], v[label("v", i, j)]], K))

    # Dynamics constraints
    for i in range(len(B)):
        for j in range(N):
            if j == 0:
                m.addConstr(x[label("x", i, j)] == x0[i])
            else:
                m.addConstr(x[label("x", i, j)] == x[label("x", i, j - 1)] +
                            d[i][j - 1] +
                            g.quicksum(
                                B[i][k] * z[label("z", k, j - 1)]
                                for k in range(len(B))))

    # Transformation constraints
    for i in range(len(B)):
        for j in range(N - 1):
            m.addConstr(z[label("z", i, j)] <= M * u[label("u", i, j)])
            m.addConstr(z[label("z", i, j)] >= mm * u[label("u", i, j)])
            m.addConstr(
                z[label("z", i, j)] <=
                y[label("y", i, j)] - mm * (1 - u[label("u", i, j)]))
            m.addConstr(
                z[label("z", i, j)] >=
                y[label("y", i, j)] - M * (1 - u[label("u", i, j)]))

    # Control constraints
    for j in range(N - 1):
        for k in range(len(A)):
            m.addConstr(
                g.quicksum(A[k][i] * u[label("u", i, j)]
                           for i in range(len(B))) <= b[k])

    # return m, [u, x, z, y, ymin]
    return m, [u, x, z, y]


def run():

    A = [[1, 1, 0, 0],
         [0, 0, 1, 1]]
    b = [1, 1]
    B = [[-1, 0, 0, 0],
         [0, -1, 0, 0],
         [0.6, 0.5, -1, 0],
         [0, 0, 0, -1]]
    E = [[1, 1, 1, 1],
         [1, 1, 1, 1],
         [1, 1, 1, 1],
         [1, 1, 1, 1]]
    d = [[10, 0, 10, 0],
         [5, 5, 5, 5],
         [0, 0, 0, 0],
         [5, 5, 5, 5]]
    c = [20, 10, 20, 10]
    xcap = [300, 300, 300, 300]
    x0 = [0, 0, 0, 0]
    cost = [2, 1, 2, 1]
    N = 5

    m, var = create_model(A, b, B, E, d, c, xcap, x0, cost, N)
    u = var[0]
    x = var[1]
    z = var[2]
    y = var[3]
    #ymin = var[4]

#    add_always_penalized(
#        m, "alw", 0, 3, [-x[label("x", 3, j)] + 9 for j in range(N - 1)], 100, 100)
    formula = Formula(ALWAYS, bounds=[0, 3],
                      args=[Formula(EXPR,
                        args=[Signal(
                            [-x[label("x", 3, j)] + 9 for j in range(N - 1)],
                            [-xcap[3], xcap[3]])])])
    alw, bounds = add_stl_constr(m, "alw", formula)
    add_penalty(m, "alw", alw, 100)

    m.update()

    m.optimize()

    m.write('foo.lp')

    if m.status == g.GRB.status.OPTIMAL:
        print('t '),
        for i in range(len(B)):
            print('x' + str(i) + ' '),
        for i in range(len(B)):
            print('u' + str(i) + ' '),
        for i in range(len(B)):
            print('z' + str(i) + ' '),
        for i in range(len(B)):
            print('y' + str(i) + ' '),
        # for i in range(len(B)):
        #    print('ymin' + str(i) + ' '),
        print('')

        ux = m.getAttr('x', u)
        xx = m.getAttr('x', x)
        zx = m.getAttr('x', z)
        yx = m.getAttr('x', y)
        #yminx = m.getAttr('x', ymin)
        for j in range(N - 1):
            print(str(j) + ' '),
            for i in range(len(B)):
                print('%d ' % round(xx[label('x', i, j)])),
            for i in range(len(B)):
                print('%d ' % round(ux[label('u', i, j)])),
            for i in range(len(B)):
                print('%d ' % round(zx[label('z', i, j)])),
            for i in range(len(B)):
                print('%d ' % round(yx[label('y', i, j)])),
        #    for i in range(len(B)):
        #        print('%d ' % round(yminx[label('ymin', i, j)])),
            print('')

        print(str(N - 1) + ' '),
        for i in range(len(B)):
            print('%d ' % round(xx[label('x', i, N - 1)])),


if __name__ == '__main__':
    run()
