"""Independent exact checks for the n=3 collinear calculation."""

from itertools import combinations

import sympy as sp
from sympy.polys.matrices import DomainMatrix

import check_n3 as q


def test_chart_and_nonadjacent_orderings():
    lam, til = q.spinors()
    total = sum((lam[i] * til[i].T for i in ("1", "2", "3", "P")), sp.zeros(2))
    assert total == sp.zeros(2)
    a, b = q.spinor_brackets()
    assert sp.factor(a("1", "2") * b("2", "1")) == 1 - q.z
    assert sp.factor(a("2", "3") * b("3", "2")) == q.z
    assert sp.factor(a("1", "3") * b("3", "1")) == -1
    assert len(set(q.ORDERS)) == 6
    for order in q.ORDERS:
        assert order.index("a") - order.index("b") not in (-1, 1)
        assert set(order) == {"1", "2", "3", "a", "b"}
        assert tuple(i for i in order if i in {"1", "2", "3"}) == ("1", "2", "3")


def test_five_gluon_single_minus_against_independent_1993_formula():
    a, b = q.spinor_brackets()
    for order in q.ORDERS:
        for negative in ("1", "2", "3"):
            assert sp.factor(q.single_minus(order, negative, a, b)
                             - q.single_minus_1993(order, negative, a, b)) == 0


def test_all_plus_standard_order_against_2018_formula():
    a, b = q.spinor_brackets()
    o = q.ORDERS[1]  # (1,b,2,a,3), 2018 Eq. (11.9)
    s = lambda i, j: a(i, j) * b(j, i)
    direct_numerator = (
        -s("1", "b") * s("b", "2")
        -s("1", "3") * s("3", "a")
        +a("b", "2") * a("a", "3") * b("2", "a") * b("3", "b")
    )
    trace_numerator = sum(
        a(o[i], o[j]) * b(o[j], o[k]) * a(o[k], o[l]) * b(o[l], o[i])
        for i, j, k, l in combinations(range(5), 4)
    )
    assert sp.factor(direct_numerator - trace_numerator) == 0


def test_complete_tree_kernel_space_and_calibration():
    tree, _, _, _ = q.data()
    negative_graviton = q.negative_graviton_tree()
    standard = q.standard_kernel()
    assert sp.factor((sp.Matrix(tree[("1", "2")]).T * standard)[0]
                     - 1 / (q.z * (q.z - 1))) == 0
    assert all(sp.factor(tree[neg][j] - tree[("1", "2")][j]) == 0
               for neg in tree for j in range(6))
    for neg in negative_graviton:
        ratio = sp.factor(negative_graviton[neg][0] / tree[("1", "2")][0])
        assert all(sp.factor(negative_graviton[neg][j]
                             - ratio * tree[("1", "2")][j]) == 0 for j in range(6))

    U = sp.symbols("U1:4")
    V = sp.symbols("V1:4")
    pair = [U[i] * (1 - q.z) + V[i] * q.z for i in range(3)]
    polynomial = sp.Poly(q.z * (q.z - 1) * pair[0]
                         - (q.z - 1) * pair[1] + q.z * pair[2], q.z)
    assert sp.linsolve(polynomial.all_coeffs(), (*U, *V)) == sp.FiniteSet(
        (V[0], 0, V[0] - V[1], V[0], V[1], 0)
    )


def test_analytic_no_go_and_function_field_rank():
    ym_residual, target_residual = q.analytic_obstruction()
    assert ym_residual == [0] * 6
    assert sp.factor(target_residual - (2*q.x**2 - 2*q.x - 3)/(6*q.z)) == 0

    A, rhs, _ = q.build_system([sp.Rational(2, 7), sp.Rational(3, 7), sp.Rational(4, 7)])
    field = sp.QQ.frac_field(q.x)
    assert len(DomainMatrix.from_Matrix(A).convert_to(field).rref()[1]) == 9
    assert len(DomainMatrix.from_Matrix(A.row_join(rhs)).convert_to(field).rref()[1]) == 10

    # This point was not used in the matrix rank test.
    xx, zz = sp.Rational(2, 5), sp.Integer(3)
    assert q.analytic_obstruction(xx, zz)[1] == -sp.Rational(29, 150)
