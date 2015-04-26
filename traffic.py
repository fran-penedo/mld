from mld import *


class TrafficModel(object):

    def __init__(self, alpha, beta, d, c, xcap, A, b):
        """TODO: to be defined1.

        :alpha: Capacity available to upstream matrix. Diagonal should be 1.
        :beta: Turning ratio matrix. Diagonal should be -1.
        :d: External inputs for each link.
        :c: Flow capacity for each link.
        :xcap: Capacity of each link.
        :A: Canonical form coefficient matrix for control variable constraints.
        :b: Right hand side vector for control variable constraints.

        """
        self._alpha = alpha
        self._beta = beta
        self._d = d
        self._c = c
        self._xcap = xcap
        self._A = A
        self._b = b


def create_milp(model, x0, cost, N):
    # Some aliases to reduce clutter
    E = model._alpha
    B = model._beta
    d = model._d
    c = model._c
    xcap = model._xcap
    A = model._A
    b = model._b

    # Big number for min alpha/beta (xcap - x)
    K_out = max(xcap) / min([min([1] + [abs(i) for i in row if i > 0])
                             for row in B])
    # Big number required for y = min(x, c, min)
    # No need to add bound on x because it's included in K_out
    K = max(max(c), K_out) * 10

    m = g.Model("traffic")
    m.modelSense = g.GRB.MINIMIZE

    x = {}
    y = {}
    z = {}
    u = {}
    v = {}
    # actual decision variables
    for i in range(len(B)):
        for j in range(N):
            # system state
            labelx = label("x", i, j)
            # control variable
            labelu = label("u", i, j)

            x[labelx] = m.addVar(
                obj=cost[i], lb=0, ub=g.GRB.INFINITY, name=labelx)
            if j < N - 1:
                u[labelu] = m.addVar(vtype=g.GRB.BINARY, name=labelu)

    m.update()

    # y, v and z variables
    for i in range(len(B)):
        for j in range(N - 1):
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
            # z = u y
            labelz = label("z", i, j)
            z[labelz] = add_milp_var(m, labelz, u[label("u", i, j)],
                                     y[label("y", i, j)],
                                     max(xcap), 0)

    # Dynamics constraints
    for i in range(len(B)):
        for j in range(N):
            if j == 0:
                m.addConstr(x[label("x", i, j)] == x0[i])
            else:
                # x = x + d + u*beta*y = x + d + beta*z
                m.addConstr(x[label("x", i, j)] == x[label("x", i, j - 1)] +
                            d[i][j - 1] +
                            g.quicksum(
                                B[i][k] * z[label("z", k, j - 1)]
                                for k in range(len(B))))

    # Control constraints: Au <= b
    for j in range(N - 1):
        for k in range(len(A)):
            m.addConstr(
                g.quicksum(A[k][i] * u[label("u", i, j)]
                           for i in range(len(B))) <= b[k])

    return m, [u, x, z, y]


def label(name, i, j):
    return name + "_" + str(i) + "_" + str(j)


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

    model = TrafficModel(E, B, d, c, xcap, A, b)
    x0 = [0, 0, 0, 0]
    cost = [2, 1, 2, 1]
    N = 5

    m, var = create_milp(model, x0, cost, N)
    u = var[0]
    x = var[1]
    z = var[2]
    y = var[3]
    #ymin = var[4]

#    add_always_penalized(
# m, "alw", 0, 3, [-x[label("x", 3, j)] + 9 for j in range(N - 1)], 100,
# 100)
    formula = Formula(ALWAYS, bounds=[0, 3],
                      args=[Formula(EXPR,
                                    args=[Signal(
                                        [-x[label("x", 3, j)] +
                                         9 for j in range(N - 1)],
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
