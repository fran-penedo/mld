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
        self._for = []


    def add_formula(self, f, cost=None, label=""):
        self._for.append((f, cost, label))


    def run_once(self, u, x0):
        v = [min([self._xcap] +
                 [self._alpha[k][i] / self._beta[k][i] *
                     (self._xcap[i] - x0[i]) for k in range(len(self._beta))
                     if self._beta[k][i] > 0])
             for i in range(len(self._beta))]

        z = [u[i] * min([x0[i], self._c[i], v[i]])
             for i in range(len(self._beta))]

        x = [x0[i] + self._d[i] +
             sum([self._beta[i][k] * z[k]
                  for k in range(len(self._beta))])
             for i in range(len(self._beta))]

        return x


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

    for f, cost, lbl in model._for:
        var, bounds = add_stl_constr(m, lbl, f)
        if cost is None:
            m.setAttr("LB", [var], [0])
        else:
            add_penalty(m, lbl, var, cost)

    return m, [u, x, z, y]


def label(name, i, j):
    return name + "_" + str(i) + "_" + str(j)


def print_solution(m, var, n, N):
    (u, x, z, y) = var

    if m.status == g.GRB.status.OPTIMAL:
        print('t '),
        for i in range(n):
            print('x' + str(i) + ' '),
        for i in range(n):
            print('u' + str(i) + ' '),
        for i in range(n):
            print('z' + str(i) + ' '),
        for i in range(n):
            print('y' + str(i) + ' '),
        print('')

        ux = m.getAttr('x', u)
        xx = m.getAttr('x', x)
        zx = m.getAttr('x', z)
        yx = m.getAttr('x', y)
        for j in range(N - 1):
            print(str(j) + ' '),
            for i in range(n):
                print('%d ' % round(xx[label('x', i, j)])),
            for i in range(n):
                print('%d ' % round(ux[label('u', i, j)])),
            for i in range(n):
                print('%d ' % round(zx[label('z', i, j)])),
            for i in range(n):
                print('%d ' % round(yx[label('y', i, j)])),
            print('')

        print(str(N - 1) + ' '),
        for i in range(n):
            print('%d ' % round(xx[label('x', i, N - 1)])),


def rhc_traffic(model, x0, cost, N, H):
    pass


