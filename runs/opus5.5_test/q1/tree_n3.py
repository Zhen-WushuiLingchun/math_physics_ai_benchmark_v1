from collinear import *
import sympy as sp
rng = random.Random(2024)
n = 3
sig = orderings(n)
print(len(sig), sig)
basis = [(1, 2), (2, 3)]          # s12, s23 independent at 4 points
x = Fr(1, 3)
hels = [(h, hP) for hP in (2, -2) for h in itertools.product([1, -1], repeat=n)]
rows, rhs = [], []
for k in range(12):
    pt = Point(n, rng)
    for hel, hP in hels:
        M = eym_tree(pt, hel, hP)
        row = []
        for s_ in sig:
            A = ym_collinear_tree(pt, s_, hel, 1 if hP > 0 else -1, x)
            for (i, j) in basis:
                row.append(pt.s(i, j)*A)
        rows.append(row); rhs.append(M)
Am = sp.Matrix(rows); b = sp.Matrix(rhs)
print('rank A =', Am.rank(), ' rank [A|b] =', Am.row_join(b).rank(), ' unknowns', Am.shape[1])
