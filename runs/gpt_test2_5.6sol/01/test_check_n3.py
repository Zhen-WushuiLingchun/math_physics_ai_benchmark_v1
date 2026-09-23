import sympy as sp

import check_n3 as c


def test_rational_chart_and_collinear_split():
    total = sum((c.lam[k] * c.tlam[k].T for k in ("1", "2", "3", "P")), sp.zeros(2))
    assert total == sp.zeros(2)
    assert c.sij("1", "2") == 1 - c.z
    assert c.sij("2", "3") == c.z
    assert c.sij("1", "3") == -1
    assert sp.simplify(c.sij("1", "2") + c.sij("2", "3") + c.sij("1", "3")) == 0
    assert sp.simplify(c.lam["a"] * c.tlam["a"].T - c.x * c.lam["P"] * c.tlam["P"].T) == sp.zeros(2)
    assert sp.simplify(c.lam["b"] * c.tlam["b"].T - (1 - c.x) * c.lam["P"] * c.tlam["P"].T) == sp.zeros(2)


def test_ordering_set_is_exactly_six_nonadjacent_insertions():
    assert len(c.ORDERS) == 6
    assert len(set(c.ORDERS)) == 6
    for order in c.ORDERS:
        hard = tuple(k for k in order if k in ("1", "2", "3"))
        assert hard == ("1", "2", "3")
        ia, ib = order.index("a"), order.index("b")
        assert (ia - ib) % 5 not in (1, 4)


def test_standard_kernel_tree_targets_are_x_independent():
    kernel = c.standard_tree_kernel()
    mhv = sp.factor(sum(k * c.tree_mhv_reduced(o) for k, o in zip(kernel, c.ORDERS)))
    anti = sp.factor(sum(k * c.tree_anti_mhv_reduced(o) for k, o in zip(kernel, c.ORDERS)))
    assert mhv == 1 / (c.z**2 * (c.z - 1))
    assert anti == -1 / c.z
    assert c.x not in mhv.free_symbols
    assert c.x not in anti.free_symbols


def test_standard_kernel_has_nonzero_all_plus_defect():
    kernel = c.standard_tree_kernel()
    defect = sp.factor(-sum(k * c.all_plus_reduced(o) for k, o in zip(kernel, c.ORDERS)))
    assert defect != 0
    assert sp.factor(defect.subs({c.x: sp.Rational(2, 5), c.z: sp.Rational(3, 7)})) == sp.Rational(76, 75)


def test_exact_function_field_inconsistency_certificate():
    matrix, rhs, _unknowns, labels = c.build_linear_system()
    chosen, certificate, contradiction = c.inconsistency_certificate(matrix, rhs, labels)
    a_sub = matrix[chosen, :]
    b_sub = rhs[chosen, :]
    assert (certificate.T * a_sub).applyfunc(sp.factor) == sp.zeros(1, matrix.cols)
    assert sp.factor((certificate.T * b_sub)[0]) == -c.x * (c.x - 1)
    assert contradiction == -c.x * (c.x - 1)
    # Two exact holdout points in the physical splitting interval.
    assert contradiction.subs(c.x, sp.Rational(2, 5)) == sp.Rational(6, 25)
    assert contradiction.subs(c.x, sp.Rational(3, 8)) == sp.Rational(15, 64)
