from traffic import *
from nose.tools import assert_list_equal


def test():
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

    model = TrafficModel(E, B, d, c, xcap, A, b)
    x0 = [0, 0, 0, 0]
    cost = [2, 1, 2, 1]
    N = 5

    def xsignal(i, j):
        return lambda milp: -milp.getVarByName(label("x", i, j)) + 9

    formula = Formula(ALWAYS, bounds=[0, 3],
                      args=[Formula(EXPR,
                                    args=[Signal(
                                        [xsignal(3, j) for j in range(N - 1)],
                                        [-xcap[3], xcap[3]])])])
    model.add_formula(formula, 100, "alw")

    m, var = create_milp(model, x0, cost, N)
    (u, x, z, y) = var

    m.update()
    m.optimize()

    m.write('foo.lp')

    print_solution(m, var, len(B), N)

    xcur = x0
    for j in range(N-1):
        ucur = [u[label("u", i, j)] for i in range(len(B))]
        xcur = model.run_once(ucur, xcur)
        assert_list_equal(xcur, [x[label("x", i, j)] for i in range(len(B))])

