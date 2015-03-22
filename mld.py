import gurobipy as g


def label(name, i, j):
    return name + "_" + str(i) + "_" + str(j)


def add_max_constr(m, label, args, K):
    def delta_label(l, i):
        return l + "_" + str(i)

    y = {}
    y[label] = m.addVar(lb=0, ub=K, name=label)
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
    def delta_label(l, i):
        return l + "_" + str(i)

    y = {}
    y[label] = m.addVar(lb=0, ub=K, name=label)
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


def always(m, label, a, b, rho, K, t=0):
    y = add_min_constr(m, label, rho[(t + a):(t + b + 1)], K)[label]
    # Just y>=0 (lower bound is part of the definition of the variable)
    return y


def create_model(A, b, B, d, c, x0, cost, N):
    # Big number required for y = min(x, c)
    K = max(c) * 100
    # Bounds for x
    mm = 0
    M = 10 * max(c)

    m = g.Model("traffic")

    x = {}
    y = {}
    #ymin = {}
    z = {}
    u = {}
    for i in range(len(B)):
        for j in range(N):
            # actual decision variables
            # system state
            labelx = label("x", i, j)
            # auxiliary variable y = min(x, c)
            # labely = label("y", i, j)
            # auxiliary binary variable needed to describe y
            # labelymin = label("ymin", i, j)
            # auxiliary binary variable z = y * u
            labelz = label("z", i, j)
            # control variable
            labelu = label("u", i, j)

            x[labelx] = m.addVar(
                obj=cost[i], lb=0, ub=g.GRB.INFINITY, name=labelx)
            if j < N - 1:
                # y[labely] = m.addVar(lb=0, ub=c[i], name=labely)
                # ymin[labelymin] = m.addVar(vtype=g.GRB.BINARY, name=labelymin)
                z[labelz] = m.addVar(lb=0, ub=g.GRB.INFINITY, name=labelz)
                u[labelu] = m.addVar(vtype=g.GRB.BINARY, name=labelu)

    m.modelSense = g.GRB.MINIMIZE
    m.update()

    for i in range(len(B)):
        for j in range(N):
            y.update(
                add_min_constr(
                    m, label("y", i, j), [x[label("x", i, j)], c[i]], K))

    # Dynamics constraints
    for i in range(len(B)):
        for j in range(N):
            if j == 0:
                m.addConstr(x[label("x", i, j)] == x0[i])
            else:
                m.addConstr(x[label("x", i, j)] == x[label("x", i, j - 1)] +
                            d[i][j - 1] +
                            g.quicksum(B[i][k] * z[label("z", k, j - 1)]
                                       for k in range(len(B))))

    # Transformation constraints
    for i in range(len(B)):
        for j in range(N - 1):
            #m.addConstr(y[label("y", i, j)] <= x[label("x", i, j)])
            # m.addConstr(
            #    y[label("y", i, j)] >= x[label("x", i, j)] - K * ymin[label("ymin", i, j)])
            # m.addConstr(
            # y[label("y", i, j)] >= c[i] - K * (1 - ymin[label("ymin", i, j)]))
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
    d = [[10, 0, 10, 0],
         [5, 5, 5, 5],
         [0, 0, 0, 0],
         [5, 5, 5, 5]]
    c = [20, 10, 20, 10]
    x0 = [0, 0, 0, 0]
    cost = [2, 1, 2, 1]
    N = 5

    m, var = create_model(A, b, B, d, c, x0, cost, N)
    u = var[0]
    x = var[1]
    z = var[2]
    y = var[3]
    #ymin = var[4]

    always(
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


def foo():
    m = g.Model("mip1")

    # Create variables
    x = m.addVar(vtype=g.GRB.BINARY, name="x")
    y = m.addVar(vtype=g.GRB.BINARY, name="y")
    z = m.addVar(vtype=g.GRB.BINARY, name="z")

    # Integrate new variables
    m.update()

    # Set objective
    m.setObjective(x + y + 2 * z, g.GRB.MAXIMIZE)

    # Add constraint: x + 2 y + 3 z <= 4
    m.addConstr(x + 2 * y + 3 * z <= 4, "c0")

    m.optimize()

    print(len(m.getConstrs()))


if __name__ == '__main__':
    run()
