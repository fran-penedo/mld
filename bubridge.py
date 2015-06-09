import traffic as t
from stl import *
import matplotlib.pyplot as plt

data = {
    'cost': [],
    'robustness': {}
}
cost = [1 for i in range(13)]


def log(model, x, u, r):
    data['cost'].append([cost[i] * x[i] for i in range(len(x))])
    for phi, rho in r.items():
        if phi in data['robustness']:
            data['robustness'][phi].append(rho)


def plot():
    f = plt.figure()
    f1 = f.add_subplot(111)
    f1.plot(data['robustness']['ncarsf'])
    plt.show()



def run():
    beta = [[0 if i != j else -1 for i in range(13)] for j in range(13)]
    beta[9][1] = 1
    beta[8][9] = 0.8
    beta[11][5] = 0.2
    beta[12][5] = 0.8
    beta[10][7] = 1
    beta[8][2] = 0.2
    beta[10][2] = 0.2
    beta[0][2] = 0.6
    beta[6][3] = 0.7
    beta[7][3] = 0.2
    beta[6][4] = 0.3
    beta[7][4] = 0.7

    alpha = [[1 for i in range(13)] for j in range(13)]

    xcap = [30 for i in range(13)]
    xcap[0] = 10
    xcap[11] = 10
    xcap[12] = 10
    xcap[1] = xcap[2] = xcap[3] = xcap[4] = xcap[5] = 10000

    c = [6, 6, 6, 6, 6, 2, 6, 6, 6, 6, 3, 6, 6]

    d = [0 for i in range(13)]
    d[1] = 6
    d[2] = 4
    d[3] = 6
    d[4] = 6
    d[5] = 6

    cons = [
        ([(1, -1)], -1),
        ([(5, -1)], -1),
        ([(11, -1)], -1),
        ([(0, -1)], -1),
        ([(10, -1)], -1),
        ([(9, 2), (12, 1), (2, 1)], 2),
        ([(7, 2), (12, 1), (2, 1)], 2),
        ([(3, 2), (8, 1), (4, 1)], 2)
    ]

    A = [[0 for i in range(13)] for j in range(len(cons))]
    for i, lhs in enumerate(zip(*cons)[0]):
        for j, a in lhs:
            A[i][j] = a
    b = zip(*cons)[1]

    x0 = [0 for i in range(13)]
    cost = [1 for i in range(13)]

    model = t.TrafficModel(alpha, beta, d, c, xcap, A, b)

    def slabel(s, i):
        return lambda j: t.label(s, i, j)

    lightf = Formula(OR, [
        Formula(NOT, [
            Formula(AND, [
                Formula(EXPR, [Signal(
                    [slabel("u", 9)], lambda x: -x[0], [-1, 1])]),
                Formula(NEXT, [
                    Formula(EXPR, [Signal(
                        [slabel("u", 9)], lambda x: x[0] - 1, [-1, 1])]),
                ])

            ])
        ]),
        Formula(NEXT, [
            Formula(NEXT, [
                Formula(EXPR, [Signal(
                    [slabel("u", 8)], lambda x: x[0] - 1, [-1, 1])]),
            ])
        ])
    ])

    blockf = Formula(OR, [
        Formula(EXPR, [Signal(
            [slabel("x", 6)], lambda x: (xcap[6] - x[0]) - 5, [-5, xcap[6]]
        )]),
        Formula(EXPR, [Signal(
            [slabel("u", 4)], lambda x: -x[0], [-1, 1]
        )])
    ])

    flowf = Formula(OR, [
        Formula(NOT, [
            Formula(EXPR, [Signal(
                [slabel("x", 4)], lambda x: x[0] - c[4], [-c[4], xcap[4]]
            )])
        ]),
        Formula(EXPR, [Signal(
            [slabel("y", 4)], lambda x: x[0] - c[4], [-c[4], xcap[4]]
        )])
    ])

    ncarsf = Formula(EVENTUALLY, bounds=[0, 2], args=[
        Formula(EXPR, [Signal(
            [slabel("x", 9), slabel("x", 10)], lambda x: -x[0] - x[1] + 10,
            [- c[9] - c[10], 10]
        )])
    ])

    model.add_formula(flowf, 100, "flowf")
    model.add_formula(blockf, 100, "blockf")
    model.add_formula(lightf, 500, "lightf")
    model.add_formula(ncarsf, 100, "ncarsf")

    data['robustness']['flowf'] = []
    data['robustness']['blockf'] = []
    data['robustness']['lightf'] = []
    data['robustness']['ncarsf'] = []

    t.rhc_traffic(model, x0, cost, 10, 2, log)
    plot()


if __name__ == "__main__":
    run()
