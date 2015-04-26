from traffic import *


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
    (u, x, z, y) = var

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

    print_solution(m, var)
