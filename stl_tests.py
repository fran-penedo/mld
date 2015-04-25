from stl import *


def test_horizon():
    f = Formula(NEXT, [
        Formula(ALWAYS, bounds=[3, 5], args=[
            Formula(EVENTUALLY, bounds=[2, 4], args=[
                Formula(AND, args=[
                    Formula(OR, args=[
                        Formula(EXPR, []),
                        Formula(ALWAYS, bounds=[1, 3],
                                args=[Formula(EXPR, [])]),
                        Formula(EVENTUALLY, bounds=[0, 2],
                                args=[Formula(EXPR, [])])
                    ]),
                    Formula(NOT, [Formula(EXPR, [])])
                ])
            ])
        ])
    ])

    assert f.horizon() == 13
