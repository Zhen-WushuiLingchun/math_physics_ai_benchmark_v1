"""Exact n=3 spinor-helicity calculation for question 1.

This file is a reproducible calculation, not an independent derivation of the
published one-loop input amplitudes.  Its source formulas are documented in
solution_n3.tex.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


x, z = sp.symbols("x z")


def spinors(xx=x, zz=z):
    """Four-point momentum-conserving chart with rational collinear spinors."""
    lam = {
        "1": sp.Matrix([1, 0]),
        "2": sp.Matrix([0, 1]),
        "3": sp.Matrix([1, 1]),
        "P": sp.Matrix([1, zz]),
    }
    til = {
        "1": sp.Matrix([-1, -1]),
        "2": sp.Matrix([-1, -zz]),
        "3": sp.Matrix([1, 0]),
        "P": sp.Matrix([0, 1]),
    }
    lam["a"], til["a"] = xx * lam["P"], til["P"]
    lam["b"], til["b"] = (1 - xx) * lam["P"], til["P"]
    return lam, til


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def spinor_brackets(xx=x, zz=z):
    lam, til = spinors(xx, zz)
    angle = lambda i, j: det(lam[i], lam[j])
    square = lambda i, j: det(til[i], til[j])
    return angle, square


ORDERS = (
    ("1", "a", "2", "b", "3"),
    ("1", "b", "2", "a", "3"),
    ("1", "a", "2", "3", "b"),
    ("1", "b", "2", "3", "a"),
    ("1", "2", "a", "3", "b"),
    ("1", "2", "b", "3", "a"),
)


def parke_taylor(order, negative, angle):
    i, j = negative
    denominator = sp.prod(angle(order[k], order[(k + 1) % 5]) for k in range(5))
    return sp.cancel(angle(i, j) ** 4 / denominator)


def anti_parke_taylor(order, positive, square):
    i, j = positive
    denominator = sp.prod(square(order[k], order[(k + 1) % 5]) for k in range(5))
    return sp.cancel(square(i, j) ** 4 / denominator)


def all_plus(order, angle, square):
    """The i/(4pi)^2 factor is removed; gluon-loop Np=2 is included."""
    numerator = sum(
        angle(order[i], order[j])
        * square(order[j], order[k])
        * angle(order[k], order[l])
        * square(order[l], order[i])
        for i, j, k, l in combinations(range(5), 4)
    )
    denominator = sp.prod(angle(order[k], order[(k + 1) % 5]) for k in range(5))
    return sp.cancel(numerator / (3 * denominator))


def single_minus(order, negative, angle, square):
    """BDK hep-ph/0505055, Eq. (4.9), with the negative leg rotated to 1."""
    j = order.index(negative)
    p = order[j:] + order[:j]
    a = lambda i, j: angle(p[i - 1], p[j - 1])
    b = lambda i, j: square(p[i - 1], p[j - 1])
    result = (
        -b(2, 5) ** 3 / (b(1, 2) * b(5, 1))
        + a(1, 4) ** 3 * b(4, 5) * a(3, 5)
        / (a(1, 2) * a(2, 3) * a(4, 5) ** 2)
        - a(1, 3) ** 3 * b(3, 2) * a(4, 2)
        / (a(1, 5) * a(5, 4) * a(3, 2) ** 2)
    ) / (3 * a(3, 4) ** 2)
    return sp.cancel(result)


def single_minus_1993(order, negative, angle, square):
    """Independent five-point form: BDK hep-ph/9302280, Eq. (4)."""
    j = order.index(negative)
    p = order[j:] + order[:j]
    a = lambda i, j: angle(p[i - 1], p[j - 1])
    b = lambda i, j: square(p[i - 1], p[j - 1])
    sij = lambda i, j: a(i, j) * b(j, i)
    numerator = (
        (sij(2, 3) + sij(3, 4) + sij(4, 5)) * b(2, 5) ** 2
        - b(2, 4) * a(4, 3) * b(3, 5) * b(2, 5)
        - b(1, 2) * b(1, 5) / (a(1, 2) * a(1, 5))
        * (
            a(1, 2) ** 2 * a(1, 3) ** 2 * b(2, 3) / a(2, 3)
            + a(1, 3) ** 2 * a(1, 4) ** 2 * b(3, 4) / a(3, 4)
            + a(1, 4) ** 2 * a(1, 5) ** 2 * b(4, 5) / a(4, 5)
        )
    )
    return sp.cancel(numerator / (3 * b(1, 2) * a(2, 3) * a(3, 4) * a(4, 5) * b(5, 1)))


def eym_single_minus(negative, angle, square):
    """Nandan-Plefka-Travaglini 1803.08497, Eq. (4.31), cyclic relabeling."""
    hard = ("1", "2", "3")
    j = hard.index(negative)
    p1, p2, p3 = hard[j:] + hard[:j]
    s = angle(p1, p2) * square(p2, p1)
    u = angle(p1, p3) * square(p3, p1)
    return sp.cancel(
        square(p2, "P") * square(p3, "P")
        / (angle(p2, "P") * angle(p3, "P"))
        * (s**2 + u**2)
        / (6 * angle(p2, p3) * square(p2, p1) * square(p3, p1))
    )


def data(xx=x, zz=z):
    angle, square = spinor_brackets(xx, zz)
    # The rational spinor gauge differs from the symmetric sqrt(x) split by
    # x(1-x) for the two positive-helicity collinear gluons.
    weight = xx * (1 - xx)
    tree = {
        neg: [sp.cancel(weight * parke_taylor(o, neg, angle)) for o in ORDERS]
        for neg in combinations(("1", "2", "3"), 2)
    }
    plus = [sp.cancel(weight * all_plus(o, angle, square)) for o in ORDERS]
    minus = {
        neg: [sp.cancel(weight * single_minus(o, neg, angle, square)) for o in ORDERS]
        for neg in ("1", "2", "3")
    }
    target_minus = {neg: eym_single_minus(neg, angle, square) for neg in ("1", "2", "3")}
    return tree, plus, minus, target_minus


def negative_graviton_tree(xx=x, zz=z):
    """The three tree calibrations with a negative-helicity graviton."""
    angle, square = spinor_brackets(xx, zz)
    hard = ("1", "2", "3")
    return {
        neg: [
            sp.cancel(anti_parke_taylor(o, tuple(j for j in hard if j != neg), square) / (xx * (1 - xx)))
            for o in ORDERS
        ]
        for neg in hard
    }


def analytic_obstruction(xx=x, zz=z):
    """Pointwise YM relation and its incompatible EYM target."""
    tree, _, minus, target = data(xx, zz)
    coefficients = ((zz - 1) * xx * (xx - 1) / 3,
                    -(zz - 1) / zz, zz - 1, sp.Integer(1))
    ym_residual = [
        sp.factor(sum(c * a for c, a in zip(coefficients, amplitudes)))
        for amplitudes in zip(tree[("1", "2")], minus["1"], minus["2"], minus["3"])
    ]
    tree_target = sp.factor((sp.Matrix(tree[("1", "2")]).T * standard_kernel(xx, zz))[0])
    eym_residual = sp.factor(sum(c * a for c, a in zip(
        coefficients, (tree_target, target["1"], target["2"], target["3"]))))
    return ym_residual, eym_residual


def standard_kernel(xx=x, zz=z):
    """One nonadjacent ordering from Stieberger-Taylor's n=3 tree formula."""
    angle, square = spinor_brackets(xx, zz)
    s2p = angle("2", "P") * square("P", "2")
    return sp.Matrix([0, xx * (1 - xx) * s2p, 0, 0, 0, 0])


def row_for_amplitudes(amplitudes, zz):
    s = 1 - zz  # s12
    t = zz      # s23
    return [sp.factor(s * a) for a in amplitudes] + [sp.factor(t * a) for a in amplitudes]


def build_system(z_values, xx=x, include_tree=True):
    rows, rhs, labels = [], [], []
    for zz in z_values:
        tree, plus, minus, target_minus = data(xx, zz)
        k = standard_kernel(xx, zz)
        if include_tree:
            for neg, amps in tree.items():
                rows.append(row_for_amplitudes(amps, zz))
                rhs.append(sp.factor((sp.Matrix(amps).T * k)[0]))
                labels.append((zz, "tree", neg))
        rows.append(row_for_amplitudes(plus, zz))
        rhs.append(sp.Integer(0))
        labels.append((zz, "all-plus", None))
        for neg, amps in minus.items():
            rows.append(row_for_amplitudes(amps, zz))
            rhs.append(sp.factor(target_minus[neg]))
            labels.append((zz, "single-minus", neg))
    return sp.Matrix(rows), sp.Matrix(rhs), labels


if __name__ == "__main__":
    print("orderings:", ORDERS)
    tree, plus, minus, target = data()
    k = standard_kernel()
    print("tree 12:", sp.factor((sp.Matrix(tree[("1", "2")]).T * k)[0]))
    print("all-plus standard defect:", sp.factor(-(sp.Matrix(plus).T * k)[0]))
    for neg in ("1", "2", "3"):
        print("single-minus", neg, "target:", sp.factor(target[neg]))
        print("single-minus", neg, "standard defect:", sp.factor(target[neg] - (sp.Matrix(minus[neg]).T * k)[0]))
