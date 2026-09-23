"""Exact low-point checks for the collinear YM/EYM question.

Uses a rational four-point spinor chart and no files outside this directory.
"""

from itertools import combinations
import sympy as sp

z, x = sp.symbols("z x")
ra, rb = sp.symbols("ra rb", nonzero=True)
lam = {"1": (1, 0), "2": (0, 1), "3": (1, 1), "P": (1, z)}
til = {"1": (-1, -1), "2": (-1, -z), "3": (1, 0), "P": (0, 1)}
lam["a"] = lam["b"] = lam["P"]
til["a"] = til["b"] = til["P"]
orders = [
    ("1", "a", "2", "b", "3"),
    ("1", "b", "2", "a", "3"),
    ("1", "a", "2", "3", "b"),
    ("1", "b", "2", "3", "a"),
    ("1", "2", "a", "3", "b"),
    ("1", "2", "b", "3", "a"),
]


def det(v, w):
    return sp.expand(v[0] * w[1] - v[1] * w[0])


def angle(i, j):
    return det(lam[i], lam[j])


def square(i, j):
    return det(til[i], til[j])


def scale(i):
    return ra if i == "a" else (rb if i == "b" else 1)


def aang(i, j):
    return angle(i, j) * scale(i) * scale(j)


def ssqr(i, j):
    return square(i, j) * scale(i) * scale(j)


def scale2(i):
    return x if i == "a" else (1 - x if i == "b" else 1)


def pt_den(order, bracket):
    out = sp.Integer(1)
    for i, j in zip(order, order[1:] + order[:1]):
        out *= bracket(i, j)
    return sp.expand(out)


def allplus_numerator(order):
    # The conventional all-plus numerator is sum tr_-(i,j,k,l).
    out = 0
    for inds in combinations(range(5), 4):
        i, j, k, l = (order[q] for q in inds)
        out += (angle(i, j) * square(j, k) * angle(k, l)
                * square(l, i) * scale2(i) * scale2(j)
                * scale2(k) * scale2(l))
    return sp.expand(out)


def singleminus(order, negative):
    """BDK 1992, Eq. (7), with the negative leg rotated to slot 1."""
    p = order.index(negative)
    q = order[p:] + order[:p]
    one, two, three, four, five = q
    term1 = ssqr(two, five)**3 / (ssqr(one, two) * ssqr(five, one))
    term2 = (aang(one, four)**3 * ssqr(four, five) * aang(three, five)
             / (aang(one, two) * aang(two, three) * aang(four, five)**2))
    term3 = -(aang(one, three)**3 * ssqr(three, two) * aang(four, two)
              / (aang(one, five) * aang(five, four) * aang(three, two)**2))
    value = sp.cancel((term1 + term2 + term3) / aang(three, four)**2)
    # Every collinear leg has even total power in this helicity sector.
    value = sp.cancel(value.subs({ra: sp.sqrt(x), rb: sp.sqrt(1-x)}))
    return value


def rows_at_x(xval):
    # Strip common i, powers of x, and loop normalization; only relative
    # values among orderings enter a linear identity.
    plus = [sp.cancel(1 / pt_den(o, angle)) for o in orders]
    minus = [sp.cancel(1 / pt_den(o, square)) for o in orders]
    loop = [sp.cancel(allplus_numerator(o) / pt_den(o, angle)) for o in orders]
    singles = [[singleminus(o, negative) for o in orders]
               for negative in ("1", "2", "3")]
    expressions = [plus, minus, loop] + singles
    expressions = [[sp.cancel(v.subs(x, xval)) for v in row] for row in expressions]
    # K_sigma = u_sigma*(1-z) + v_sigma*z.
    columns = [[sp.cancel((1-z)*f), sp.cancel(z*f)] for row in expressions for f in row]
    columns = [columns[6*r:6*(r+1)] for r in range(len(expressions))]
    rhs = [sp.cancel(-xval*(1-xval)*plus[1].subs(x, xval)),
           sp.cancel(-xval*(1-xval)*minus[1].subs(x, xval)), sp.Integer(0)]
    rhs += [eym_singleminus(hard) for hard in ("1", "2", "3")]
    rhs = [sp.cancel(v.subs(x, xval)) for v in rhs]
    equations = []
    labels = []
    sector_names = ("tree_plus", "tree_minus", "loop_allplus",
                    "loop_minus1", "loop_minus2", "loop_minus3")
    for sector, row, target in zip(sector_names, columns, rhs):
        exprs = [v for pair in row for v in pair] + [target]
        common = sp.lcm([sp.denom(v) for v in exprs])
        polys = [sp.Poly(sp.cancel(v*common), z) for v in exprs]
        for power in range(max(p.degree() for p in polys)+1):
            equations.append([p.nth(power) for p in polys])
            labels.append((sector, power))
    augmented = sp.Matrix(equations)
    return augmented[:, :-1], augmented[:, -1], expressions, labels


def invariant(i, j):
    return sp.expand(-angle(i, j) * square(i, j))


def eym_singleminus(negative):
    # Nandan--Plefka--Travaglini, Eq. (1.6), normalized by i/(48 pi^2).
    cyclic = ("1", "2", "3")
    position = cyclic.index(negative)
    positive = [cyclic[(position+1) % 3], cyclic[(position+2) % 3]]
    a, b = positive
    s = invariant(negative, a)
    u = invariant(negative, b)
    return sp.cancel(sp.Rational(1, 2) * square(a, "P") * square(b, "P")
                     / (angle(a, "P") * angle(b, "P"))
                     * (s*s + u*u)
                     / (angle(a, b) * square(a, negative) * square(b, negative)))


if __name__ == "__main__":
    import sys
    # Independent check against Nandan--Plefka--Travaglini eq. (4.9):
    # their numerator for (1,b,2,a,3) equals the trace-sum numerator.
    bdk_numerator = (
        -(1-x)**2 * invariant("1", "P") * invariant("2", "P")
        - x * invariant("1", "3") * invariant("3", "P")
        + x*(1-x) * angle("P", "2") * angle("P", "3")
        * square("2", "P") * square("3", "P")
    )
    assert sp.cancel(allplus_numerator(orders[1]) - bdk_numerator) == 0
    if "--certificate" in sys.argv:
        a, b, expr, labels = rows_at_x(sp.Rational(1, 3))
        aug = a.row_join(b)
        pivots = aug.T.rref()[1]
        minor = aug[list(pivots), :]
        print("pivot rows 0-based:", pivots)
        print("pivot labels:", [labels[p] for p in pivots])
        print("13x13 determinant:", minor.det(method="domain-ge"))
        print("rank(A),rank(A|b):", a.rank(), aug.rank())
        cvec = sp.Matrix([b[i] if labels[i][0].startswith("loop_minus") else 0
                          for i in range(len(labels))])
        tvec = b-cvec
        free = a.row_join(-cvec)
        free_aug = free.row_join(tvec)
        pivots2 = free_aug.T.rref()[1]
        print("free-normalization ranks:", free.rank(), free_aug.rank())
        print("14x14 rows 0-based:", pivots2)
        print("14x14 labels:", [labels[p] for p in pivots2])
        print("14x14 determinant:",
              free_aug[list(pivots2), :].det(method="domain-ge"))
        assert free.rank() == 13 and free_aug.rank() == 14
        assert free_aug[list(pivots2), :].det(method="domain-ge") == -2720977920
        sys.exit(0)
    for xval in [sp.Rational(1, 3), sp.Rational(2, 5), sp.Rational(1, 2)]:
        a, b, expr, labels = rows_at_x(xval)
        print("x", xval, "shape", a.shape, "rank", a.rank(),
              "aug_rank", a.row_join(b).rank())
        print("allplus representative", sp.factor(expr[2][1]))
        print("tree plus representative", sp.factor(expr[0][1]))
    print("symbolic rows")
    a, b, _, labels = rows_at_x(x)
    print("shape", a.shape, "rank", a.rank(), "aug_rank", a.row_join(b).rank())
