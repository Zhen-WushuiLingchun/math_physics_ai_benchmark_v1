import sys, pickle, itertools
from collinear import *
from multiprocessing import Pool
n = 4
sig = orderings(n)
basis = [(1, 2), (2, 3), (3, 4), (1, 3), (2, 4)]   # candidate independent set (checked below)
hels = [(h, 2) for h in itertools.product([1, -1], repeat=n) if 2 <= h.count(-1) <= 3] + \
       [(h, -2) for h in itertools.product([1, -1], repeat=n) if 1 <= h.count(-1) <= 2]
def work(args):
    seed, x = args
    rng = random.Random(seed)
    pt = Point(n, rng)
    rows, rhs = [], []
    for hel, hP in hels:
        M = eym_tree(pt, hel, hP)
        row = []
        for s_ in sig:
            A = ym_collinear_tree(pt, s_, hel, 1 if hP > 0 else -1, x)
            row += [pt.s(i, j)*A for (i, j) in basis]
        rows.append(row); rhs.append(M)
    return rows, rhs
if __name__ == '__main__':
    x = Fr(1, 3)
    with Pool(12) as pool:
        res = pool.map(work, [(1000+k, x) for k in range(12)])
    rows = sum([r[0] for r in res], []); rhs = sum([r[1] for r in res], [])
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix
    A = DomainMatrix([[sp.Rational(v.numerator, v.denominator) for v in r] for r in rows], (len(rows), len(rows[0])), sp.QQ)
    b = DomainMatrix([[sp.Rational(v.numerator, v.denominator)] for v in rhs], (len(rhs), 1), sp.QQ)
    print('x=1/3: rows', len(rows), 'rank A', A.rank(), 'rank [A|b]', A.hstack(b).rank(), 'unknowns', len(rows[0]))
    # check Mandelstam basis independence: rank of s-vectors over points
    rng = random.Random(5)
    S = [[Point(n, rng).s(i, j) for (i, j) in basis] for _ in range(8)]
    print('basis rank', sp.Matrix(S).rank())
