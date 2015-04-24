import gurobipy as g


def label(name, i, j):
    return name + "_" + str(i) + "_" + str(j)


def delta_label(l, i):
    return l + "_" + str(i)


def add_max_constr(m, label, args, K):
    y = {}
    y[label] = m.addVar(lb=0, ub=K, name=label)

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
            m.addConstr(y[label] >= args[i])
            m.addConstr(y[label] <= args[i] + y[x] * K)

        m.addConstr(g.quicksum([y[delta_label(label, i)]
                                for i in range(len(args))] == len(args) - 1))

    return y


# TODO handle len(args) == 2 differently
def add_min_constr(m, label, args, K):
    y = {}
    y[label] = m.addVar(lb=0, ub=K, name=label)

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
            m.addConstr(y[label] <= args[i])
            m.addConstr(y[label] >= args[i] - y[x] * K)

        m.addConstr(g.quicksum([y[delta_label(label, i)]
                                for i in range(len(args))]) == len(args) - 1)

    return y


def add_always_constr(m, label, a, b, rho, K, t=0):
    y = add_min_constr(m, label, rho[(t + a):(t + b + 1)], K)[label]
    # Just y>=0 (lower bound is part of the definition of the variable)
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

    m.modelSense = g.GRB.MINIMIZE
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
                                E[i][k] * B[i][k] * z[label("z", k, j - 1)]
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

    add_always_constr(
        m, "alw", 0, 3, [-x[label("x", 3, j)] + 9 for j in range(N - 1)], 100)
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
