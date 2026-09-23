"""Exact n=3 collinear-kernel experiment for the YM -> EYM problem.

This file is deliberately self-contained.  It works on a rational four-point
spinor-helicity chart, evaluates the six allowed non-adjacent orderings, and
builds the tree/all-plus/single-minus linear constraints on a general local
kernel K_sigma = a_sigma s + b_sigma t.

The overall factors i/(48 pi^2) of the one-loop Yang--Mills amplitudes are
removed.  With the conventions used below, the one-minus EYM target therefore
comes with an overall factor 1/2.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


x, z = sp.symbols("x z")
sqrt_x = sp.sqrt(x)
sqrt_y = sp.sqrt(1 - x)


# A rational chart for p1+p2+p3+P=0.  It has
# s=s12=1-z, t=s23=z, u=s13=-1.
lam = {
    "1": sp.Matrix([1, 0]),
    "2": sp.Matrix([0, 1]),
    "3": sp.Matrix([1, 1]),
    "P": sp.Matrix([1, z]),
}
tlam = {
    "1": sp.Matrix([-1, -1]),
    "2": sp.Matrix([-1, -z]),
    "3": sp.Matrix([1, 0]),
    "P": sp.Matrix([0, 1]),
}
lam["a"], tlam["a"] = sqrt_x * lam["P"], sqrt_x * tlam["P"]
lam["b"], tlam["b"] = sqrt_y * lam["P"], sqrt_y * tlam["P"]


def det2(v: sp.Matrix, w: sp.Matrix) -> sp.Expr:
    return sp.expand(v[0] * w[1] - v[1] * w[0])


def angle(i: str, j: str) -> sp.Expr:
    return det2(lam[i], lam[j])


def square(i: str, j: str) -> sp.Expr:
    return det2(tlam[i], tlam[j])


def sij(i: str, j: str) -> sp.Expr:
    # Convention s_ij=<ij>[ji].
    return sp.expand(angle(i, j) * square(j, i))


def pt(order: tuple[str, ...]) -> sp.Expr:
    return sp.prod(angle(order[i], order[(i + 1) % len(order)]) for i in range(len(order)))


def pt_square(order: tuple[str, ...]) -> sp.Expr:
    return sp.prod(square(order[i], order[(i + 1) % len(order)]) for i in range(len(order)))


def tree_mhv_reduced(order: tuple[str, ...]) -> sp.Expr:
    """Parke--Taylor amplitude with its ordering-independent numerator removed."""
    return sp.factor(1 / pt(order))


def tree_anti_mhv_reduced(order: tuple[str, ...]) -> sp.Expr:
    """Parity-conjugate Parke--Taylor amplitude, numerator removed."""
    return sp.factor(1 / pt_square(order))


def tr_minus(i: str, j: str, k: str, ell: str) -> sp.Expr:
    return angle(i, j) * square(j, k) * angle(k, ell) * square(ell, i)


def all_plus_reduced(order: tuple[str, ...]) -> sp.Expr:
    """One-loop all-plus YM amplitude with i/(48 pi^2) removed.

    In this normalization the numerator is the sum of chiral traces over all
    ordered four-subsets.
    """
    numerator = sum(tr_minus(*(order[i] for i in inds)) for inds in combinations(range(5), 4))
    return sp.factor(numerator / pt(order))


def single_minus_reduced(order: tuple[str, ...], negative: str) -> sp.Expr:
    """Five-point one-minus YM amplitude with i/(48 pi^2) removed.

    Rotate the cyclic ordering so that the negative-helicity leg is first and
    use the compact five-point rational form.  All other legs are positive.
    """
    p = order.index(negative)
    q = order[p:] + order[:p]
    one, two, three, four, five = q
    assert one == negative
    term1 = -square(two, five) ** 3 / (square(one, two) * square(five, one))
    term2 = (
        angle(one, four) ** 3
        * square(four, five)
        * angle(three, five)
        / (angle(one, two) * angle(two, three) * angle(four, five) ** 2)
    )
    term3 = -(
        angle(one, three) ** 3
        * square(three, two)
        * angle(four, two)
        / (angle(one, five) * angle(five, four) * angle(three, two) ** 2)
    )
    return sp.factor((term1 + term2 + term3) / angle(three, four) ** 2)


ORDERS: tuple[tuple[str, ...], ...] = (
    ("1", "a", "2", "b", "3"),
    ("1", "a", "2", "3", "b"),
    ("1", "b", "2", "a", "3"),
    ("1", "2", "a", "3", "b"),
    ("1", "b", "2", "3", "a"),
    ("1", "2", "b", "3", "a"),
)


def eym_single_minus_target(negative: str) -> sp.Expr:
    """Four-point EYM target divided by i/(48 pi^2)."""
    positives = [leg for leg in ("1", "2", "3") if leg != negative]
    j, k = positives
    spin = square(j, "P") * square(k, "P") / (angle(j, "P") * angle(k, "P"))
    denom = angle(j, k) * square(j, negative) * square(k, negative)
    invariant = sij(negative, j) ** 2 + sij(negative, k) ** 2
    return sp.factor(sp.Rational(1, 2) * spin * invariant / denom)


def standard_tree_kernel() -> list[sp.Expr]:
    # K_1=x(1-x)s_{2P}, K_2=x(1-x)s_{3P}; all other entries vanish.
    # The x(1-x) factor cancels the homogeneous rescaling of the two
    # non-adjacent split legs and is fixed by the x-independent tree target.
    split_weight = x * (1 - x)
    return [split_weight * sij("2", "P"), split_weight * sij("3", "P"), 0, 0, 0, 0]


def numerator_coefficients(
    expr: sp.Expr, variable: sp.Symbol, prefix: str
) -> tuple[list[sp.Expr], list[str]]:
    """Coefficients and labels of the cleared numerator of an identity."""
    num, _den = sp.fraction(sp.cancel(expr))
    poly = sp.Poly(sp.expand(num), variable)
    powers = range(poly.degree(), -1, -1)
    return [sp.factor(c) for c in poly.all_coeffs()], [f"{prefix}:z^{p}" for p in powers]


def build_linear_system() -> tuple[sp.Matrix, sp.Matrix, tuple[sp.Symbol, ...], list[str]]:
    a = sp.symbols("a1:7")
    b = sp.symbols("b1:7")
    unknowns = a + b
    s, t = 1 - z, z
    n_kernel = [a[i] * s + b[i] * t for i in range(6)]
    full_kernel = [sp.expand(k0 + n) for k0, n in zip(standard_tree_kernel(), n_kernel)]

    equations: list[sp.Expr] = []
    labels: list[str] = []

    def add_identity(expr: sp.Expr, prefix: str) -> None:
        coeffs, names = numerator_coefficients(expr, z, prefix)
        equations.extend(coeffs)
        labels.extend(names)

    tree_values = [tree_mhv_reduced(o) for o in ORDERS]
    add_identity(sum(n * v for n, v in zip(n_kernel, tree_values)), "tree-MHV")
    anti_tree_values = [tree_anti_mhv_reduced(o) for o in ORDERS]
    add_identity(sum(n * v for n, v in zip(n_kernel, anti_tree_values)), "tree-anti-MHV")

    all_plus_values = [all_plus_reduced(o) for o in ORDERS]
    add_identity(sum(k * v for k, v in zip(full_kernel, all_plus_values)), "loop-all-plus")

    for negative in ("1", "2", "3"):
        single_values = [single_minus_reduced(o, negative) for o in ORDERS]
        residual = sum(k * v for k, v in zip(full_kernel, single_values))
        residual -= eym_single_minus_target(negative)
        add_identity(residual, f"loop-one-minus-{negative}")

    matrix, rhs = sp.linear_eq_to_matrix(equations, unknowns)
    return matrix.applyfunc(sp.factor), rhs.applyfunc(sp.factor), unknowns, labels


def inconsistency_certificate(
    matrix: sp.Matrix, rhs: sp.Matrix, labels: list[str]
) -> tuple[list[int], sp.Matrix, sp.Expr]:
    """Return a small exact left-null certificate y^T A=0, y^T b!=0."""
    chosen: list[int] = []
    aug = matrix.row_join(rhs)
    current_rank = 0
    for i in range(aug.rows):
        trial = chosen + [i]
        trial_rank = aug[trial, :].rank()
        if trial_rank > current_rank:
            chosen = trial
            current_rank = trial_rank
        if current_rank == aug.rank():
            break

    a_sub = matrix[chosen, :]
    b_sub = rhs[chosen, :]
    candidates = a_sub.T.nullspace()
    certificate = next(v for v in candidates if sp.factor((v.T * b_sub)[0]) != 0)

    denominators = [sp.denom(sp.cancel(v)) for v in certificate]
    common = sp.lcm(denominators)
    certificate = certificate.applyfunc(lambda v: sp.factor(sp.cancel(common * v)))
    contradiction = sp.factor((certificate.T * b_sub)[0])

    assert all(sp.factor(v) == 0 for v in certificate.T * a_sub)
    assert contradiction != 0
    assert len(chosen) == len([labels[i] for i in chosen])
    return chosen, certificate, contradiction


def main() -> None:
    print("MOMENTUM INVARIANTS")
    print("s12 =", sij("1", "2"), "s23 =", sij("2", "3"), "s13 =", sij("1", "3"))
    print()
    print("ORDERED DATA (overall loop constant removed)")
    for i, order in enumerate(ORDERS, 1):
        print(f"sigma_{i} = {order}")
        print("  tree/PT  =", sp.factor(tree_mhv_reduced(order)))
        print("  anti/PT  =", sp.factor(tree_anti_mhv_reduced(order)))
        print("  all-plus =", sp.factor(all_plus_reduced(order)))
        for negative in ("1", "2", "3"):
            print(f"  1minus({negative}) =", sp.factor(single_minus_reduced(order, negative)))
    print()

    k0 = standard_tree_kernel()
    print("STANDARD TREE KERNEL")
    print([sp.factor(k) for k in k0])
    all_plus_defect = -sum(k * all_plus_reduced(o) for k, o in zip(k0, ORDERS))
    print("Delta_all_plus / [i/(48 pi^2)] =", sp.factor(all_plus_defect))
    for negative in ("1", "2", "3"):
        reconstructed = sum(k * single_minus_reduced(o, negative) for k, o in zip(k0, ORDERS))
        defect = eym_single_minus_target(negative) - reconstructed
        print(f"Delta_1minus({negative}) / [i/(48 pi^2)] =", sp.factor(defect))
    print()

    matrix, rhs, unknowns, labels = build_linear_system()
    augmented = matrix.row_join(rhs)
    print("LINEAR SYSTEM")
    print("shape =", matrix.shape)
    print("rank(A) =", matrix.rank())
    print("rank([A|b]) =", augmented.rank())
    solution = sp.linsolve((matrix, rhs), unknowns)
    print("solution =", solution)
    chosen, certificate, contradiction = inconsistency_certificate(matrix, rhs, labels)
    print("certificate rows and coefficients:")
    for row, coefficient in zip(chosen, certificate):
        print(" ", labels[row], ":", sp.factor(coefficient))
    print("y^T b =", contradiction)


if __name__ == "__main__":
    main()
