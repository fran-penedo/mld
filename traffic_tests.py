from traffic import *
from nose.tools import assert_almost_equals

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
d = [8, 5, 0, 5]
c = [20, 10, 20, 10]
xcap = [300, 300, 300, 300]

x0 = [0, 0, 0, 0]
cost = [2, 1, 2, 1]
N = 5


def test_rhc():
    model = TrafficModel(E, B, d, c, xcap, A, b)

    formula = Formula(EVENTUALLY, bounds=[0, 3],
                      args=[Formula(EXPR,
                                    args=[Signal(
                                        [lambda j: label("x", 1, j)],
                                        lambda x: -x[0] + 9,
                                        [-xcap[1], xcap[1]])])])
    model.add_formula(formula, 100, "eve")
    rhc_traffic(model, x0, cost, 10, 2)


def test_milp():
    model = TrafficModel(E, B, d, c, xcap, A, b)
    formula = Formula(ALWAYS, bounds=[0, 3],
                      args=[Formula(EXPR,
                                    args=[Signal(
                                        [lambda j: label("x", 3, j)],
                                        lambda x: -x[0] + 9,
                                        [-xcap[3], xcap[3]])])])
    model.add_formula(formula, 100, "alw")

    m, var = create_milp(model, x0, cost, N)
    (u, x, z, y) = var

    m.update()
    m.optimize()

    m.write('foo.lp')

    print_solution(m, var, len(B), N)

    xcur = x0
    uu = m.getAttr("x", u)
    xx = m.getAttr("x", x)
    for j in range(N - 1):
        ucur = [uu[label("u", i, j)] for i in range(len(B))]
        xnext, bcur = model.run_once(ucur, xcur)
        for nom, sim in zip(xcur,
                            [xx[label("x", i, j)] for i in range(len(B))]):
            assert_almost_equals(nom, sim)
        xcur = xnext

    for nom, sim in zip(xcur,
                        [xx[label("x", i, N - 1)] for i in range(len(B))]):
        assert_almost_equals(nom, sim)
