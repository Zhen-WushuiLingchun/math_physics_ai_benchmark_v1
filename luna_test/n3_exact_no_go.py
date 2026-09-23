# -*- coding: utf-8 -*-
"""Exact n=3 collinear EYM -> YM one-loop no-go audit.

All files/results are generated within the task folder. Arithmetic in the
kernel certificate is exact SymPy rational arithmetic, not floating point.
"""
from itertools import combinations
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ

x, t = sp.symbols("x t", nonzero=True)
zeta = 1 - 1 / t
sqrt_x = sp.sqrt(x)
sqrt_1mx = sp.sqrt(1 - x)

# Four-point massless kinematics (1,2,3,P), followed by the collinear pair
# a=xP, b=(1-x)P. This gives s12=-1, s23=t, s13=1-t.
lam = {
    "1": (1, 0),
    "2": (0, 1),
    "3": (1, 1),
    "P": (zeta, 1),
    "a": (sqrt_x * zeta, sqrt_x),
    "b": (sqrt_1mx * zeta, sqrt_1mx),
}
tilde = {
    "1": (1, 0),
    "2": (0, 1),
    "3": (1 / (zeta - 1), -zeta / (zeta - 1)),
    "P": (-1 / (zeta - 1), 1 / (zeta - 1)),
    "a": (-sqrt_x / (zeta - 1), sqrt_x / (zeta - 1)),
    "b": (-sqrt_1mx / (zeta - 1), sqrt_1mx / (zeta - 1)),
}


def det2(v, w):
    return v[0] * w[1] - v[1] * w[0]


def angle(i, j):
    return sp.factor(det2(lam[i], lam[j]))


def square(i, j):
    return sp.factor(det2(tilde[i], tilde[j]))


def sij(i, j):
    # s_ij = <ij>[ji]
    return sp.factor(angle(i, j) * square(j, i))


# Pi_3: hard cyclic order (1,2,3), with a and b in distinct gaps.
ORDERS = [
    ("1", "a", "2", "b", "3"),
    ("1", "b", "2", "a", "3"),
    ("1", "a", "2", "3", "b"),
    ("1", "b", "2", "3", "a"),
    ("1", "2", "a", "3", "b"),
    ("1", "2", "b", "3", "a"),
]


def cyclic_denominator(order, bracket):
    return sp.prod(bracket(order[k], order[(k + 1) % len(order)])
                   for k in range(len(order)))


def ym_allplus_function(order):
    """F_sigma such that A5^(1)(all+) = i/(48 pi^2) F_sigma.

    This is the compact five-gluon all-plus formula. Its numerator is also
    checked below against the sum of chiral traces.
    """
    o = order
    numerator = (
        -sij(o[0], o[1]) * sij(o[1], o[2])
        -sij(o[4], o[0]) * sij(o[3], o[4])
        + angle(o[1], o[2]) * angle(o[3], o[4])
        * square(o[2], o[3]) * square(o[4], o[1])
    )
    return sp.factor(numerator / cyclic_denominator(o, angle))


def ym_allplus_trace_numerator(order):
    """Return sum tr_-(ijkl), with tr_-=2<ij>[jk]<kl>[li]."""
    total = 0
    for inds in combinations(range(5), 4):
        i, j, k, ell = [order[m] for m in inds]
        total += 2 * angle(i, j) * square(j, k) * angle(k, ell) * square(ell, i)
    return sp.factor(total)


def ym_singleminus_function(order, negative_leg):
    """G_{sigma,h} such that A5^(1)(h-,others+) = i/(48 pi^2) G.

    Uses the BDK five-gluon single-minus expression, after cyclically rotating
    the color ordering so that the negative-helicity gluon is first.
    """
    o = list(order)
    start = o.index(negative_leg)
    a, b, c, d, e = o[start:] + o[:start]
    value = (
        square(b, e) ** 3 / (square(a, b) * square(e, a))
        + angle(a, d) ** 3 * square(d, e) * angle(c, e)
        / (angle(a, b) * angle(b, c) * angle(d, e) ** 2)
        - angle(a, c) ** 3 * square(c, b) * angle(d, b)
        / (angle(a, e) * angle(e, d) * angle(c, b) ** 2)
    ) / angle(c, d) ** 2
    return sp.factor(value)


def eym_singleminus_function(h, j, k):
    """E_h such that M_EYM^(1) = i/(4 pi)^2 E_h for cyclic (h,j,k)."""
    return sp.factor(
        square(j, "P") * square(k, "P")
        / (angle(j, "P") * angle(k, "P"))
        * (sij(h, j) ** 2 + sij(h, k) ** 2)
        / (6 * angle(j, k) * square(j, h) * square(k, h))
    )


# Basic kinematic checks and the independent Mandelstam basis.
assert sp.simplify(sij("1", "2") + 1) == 0
assert sp.simplify(sij("2", "3") - t) == 0
assert sp.simplify(sij("1", "3") - (1 - t)) == 0
assert sp.simplify(sij("1", "P") - t) == 0
assert sp.simplify(sij("2", "P") - (1 - t)) == 0
assert sp.simplify(sij("3", "P") + 1) == 0

# Exact tree calibration: the representative has only sigma=(1,b,2,a,3).
# The sign is fixed by the CHY single-graviton tree expansion in the stated
# F=...+g f convention.
K_tree = -x * (1 - x) * sij("2", "P")
order_tree = ORDERS[1]
A5_tree = sp.I * angle("1", "2") ** 4 / cyclic_denominator(order_tree, angle)
M_collinear_tree = sp.factor(K_tree * A5_tree)

# Independent insertion-form check at reference spinor q=1.
A4_1P23 = sp.I * angle("1", "2") ** 4 / (
    angle("1", "P") * angle("P", "2") * angle("2", "3") * angle("3", "1")
)
A4_12P3 = sp.I * angle("1", "2") ** 4 / (
    angle("1", "2") * angle("2", "P") * angle("P", "3") * angle("3", "1")
)
ref1_dot = lambda i: sp.factor(angle("1", i) * square(i, "P") / angle("1", "P"))
M_insertion_tree = sp.factor(
    -ref1_dot("1") * A4_1P23
    - (ref1_dot("1") + ref1_dot("2")) * A4_12P3
)
assert sp.simplify(M_collinear_tree - M_insertion_tree) == 0
assert sp.simplify(M_collinear_tree + sp.I * t ** 3 / (t - 1)) == 0

# Representative all-plus defect. The EYM all-plus target is zero at eps^0.
F_tree_order = ym_allplus_function(order_tree)
delta_allplus = sp.factor(-K_tree * sp.I / (48 * sp.pi ** 2) * F_tree_order)
delta_closed = sp.I / (48 * sp.pi ** 2) * t ** 3 * (t * (1 - x) ** 2 - x ** 2)
assert sp.simplify(delta_allplus - delta_closed) == 0

# Independent all-plus check: sum of chiral traces equals twice the compact
# numerator (the factor two is the standard tr_- convention).
compact_numerator = sp.factor(F_tree_order * cyclic_denominator(order_tree, angle))
assert sp.simplify(
    ym_allplus_trace_numerator(order_tree) - 2 * compact_numerator
) == 0
assert sp.simplify(delta_closed.subs({x: sp.Rational(1, 2), t: 2})
                   - sp.I / (24 * sp.pi ** 2)) == 0

# To test the entire allowed kernel space, let every ordering have an
# independent linear kernel K_sigma = U_sigma*s12 + V_sigma*s23.
U = sp.symbols("U0:6")
V = sp.symbols("V0:6")
K = [-U[i] + t * V[i] for i in range(6)]  # s12=-1, s23=t

loop_equations = [
    sum(K[i] * ym_allplus_function(o) for i, o in enumerate(ORDERS)),
]
# Use two required one-minus positions; the third would add constraints.
for h, j, k in [("1", "2", "3"), ("2", "3", "1")]:
    loop_equations.append(
        sum(K[i] * ym_singleminus_function(o, h)
            for i, o in enumerate(ORDERS))
        - 3 * eym_singleminus_function(h, j, k)
    )


def exact_constraint_matrix(x0):
    rows, rhs, tags = [], [], []
    for eq_index, equation in enumerate(loop_equations):
        at_x = sp.factor(equation.subs(x, x0))
        numerator = sp.factor(sp.radsimp(sp.cancel(at_x).as_numer_denom()[0]))
        poly = sp.Poly(sp.expand(numerator), t)
        for degree in range(poly.degree(), -1, -1):
            coefficient = sp.expand(poly.coeff_monomial(t ** degree))
            if coefficient == 0:
                continue
            rows.append([sp.factor(sp.diff(coefficient, var)) for var in U + V])
            rhs.append(sp.factor(-coefficient.subs({var: 0 for var in U + V})))
            tags.append((eq_index, degree))
    return sp.Matrix(rows), sp.Matrix(rhs), tags


def function_field_constraint_matrix():
    """Polynomial coefficient system over Q(x), before x specialization.

    Clearing each equation's rational denominator in (x,t) gives an
    equivalent polynomial identity.  Its coefficients are rational functions
    of x and define the function-field matrix used in the rank argument.
    """
    rows, rhs, tags = [], [], []
    for eq_index, equation in enumerate(loop_equations):
        numerator = sp.together(sp.cancel(equation)).as_numer_denom()[0]
        poly = sp.Poly(sp.expand(numerator), t)
        for degree in range(poly.degree(), -1, -1):
            coefficient = sp.cancel(poly.coeff_monomial(t ** degree))
            if coefficient == 0:
                continue
            rows.append([sp.cancel(sp.diff(coefficient, var)) for var in U + V])
            rhs.append(sp.cancel(-coefficient.subs({var: 0 for var in U + V})))
            tags.append((eq_index, degree))
    return sp.Matrix(rows), sp.Matrix(rhs), tags


M, rhs, tags = exact_constraint_matrix(sp.Rational(1, 3))
assert tags == [
    (0, 5), (0, 4), (0, 3),
    (1, 8), (1, 7), (1, 6), (1, 5), (1, 4), (1, 3),
    (2, 6), (2, 5), (2, 4), (2, 3), (2, 2), (2, 1),
]
assert all(entry.is_Rational for entry in M)
assert all(entry.is_Rational for entry in rhs)
rank_M = DomainMatrix.from_Matrix(M).convert_to(QQ).rank()
rank_aug = DomainMatrix.from_Matrix(M.row_join(rhs)).convert_to(QQ).rank()
assert M.shape == (15, 12)
assert (rank_M, rank_aug) == (12, 13)

# Explicit left-null certificate for this row order: lambda^T M=0 but
# lambda^T rhs=-253920. Thus the finite exact system is inconsistent.
certificate = sp.Matrix([[
    169680, 162288, -296352,
    -122163, -402012, -21876, 57624, -77616, -500192,
    915417, 225276, 233556, 212856, 0, 0,
]])
assert certificate * M == sp.zeros(1, 12)
assert (certificate * rhs)[0] == -253920

# Independent exact rational specialization.
M_quarter, rhs_quarter, _ = exact_constraint_matrix(sp.Rational(1, 4))
assert DomainMatrix.from_Matrix(M_quarter).convert_to(QQ).rank() == 12
assert DomainMatrix.from_Matrix(M_quarter.row_join(rhs_quarter)).convert_to(QQ).rank() == 13

# Function-field no-go: specialize the coefficients of the unspecialized
# rational-function system. A nonzero augmented minor at x=1/3 proves that
# the same minor is not the zero element of Q(x), even if a candidate kernel
# itself were allowed to have a pole there.
M_x, rhs_x, tags_x = function_field_constraint_matrix()
x_third = sp.Rational(1, 3)
for entry in M_x.row_join(rhs_x):
    assert sp.denom(sp.cancel(entry)).subs(x, x_third) != 0
M_x_third = M_x.subs(x, x_third)
rhs_x_third = rhs_x.subs(x, x_third)
rank_M_x_third = DomainMatrix.from_Matrix(M_x_third).convert_to(QQ).rank()
rank_aug_x_third = DomainMatrix.from_Matrix(
    M_x_third.row_join(rhs_x_third)
).convert_to(QQ).rank()
assert (rank_M_x_third, rank_aug_x_third) == (12, 13)

print("Pi_3 order count: 6")
print("tree calibration: passed; stripped MHV value = -I*t^3/(t-1)")
print("all-plus representative defect: I/(48*pi^2)*t^3*(t*(1-x)^2-x^2)")
print("exact x=1/3 loop constraints: 15 equations, 12 unknowns")
print("rank(M), rank([M|b]) over QQ: 12, 13")
print("left-null certificate contraction: -253920")
print("independent exact x=1/4 rank check: 12, 13")
print("function-field matrix: regular at x=1/3; specialized ranks 12, 13")
