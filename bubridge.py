import traffic as t


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
    for lhs in zip(*cons)[0]:
        for i, a in lhs:
            A[i] = a
    b = zip(*cons)[1]

    x0 = [0 for i in range(13)]
    cost = [1 for i in range(13)]

    model = t.TrafficModel(alpha, beta, d, c, xcap, A, b)

    t.rhc_traffic(model, x0, cost, 10, 2)


if __name__ == "__main__":
    run()
