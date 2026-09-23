import pickle, glob, itertools, sympy as sp
from collinear import *
from ym_formulas import *
from multiprocessing import Pool
from sympy.polys.matrices import DomainMatrix
n = 4
sig = orderings(n)
basis = [(1, 2), (2, 3), (3, 4), (1, 3), (2, 4)]
hels = [(h, 2) for h in itertools.product([1, -1], repeat=n) if 2 <= h.count(-1) <= 3] + \
       [(h, -2) for h in itertools.product([1, -1], repeat=n) if 1 <= h.count(-1) <= 2]
def tree_rows(args):
    seed, x = args
    rng = random.Random(seed); pt = Point(n, rng)
    rows, rhs = [], []
    for hel, hP in hels:
        M = eym_tree(pt, hel, hP)
        row = []
        for s_ in sig:
            A = ym_collinear_tree(pt, s_, hel, 1 if hP > 0 else -1, x)
            row += [pt.s(i, j)*A for (i, j) in basis]
        rows.append(row); rhs.append(M)
    return rows, rhs
def sval(lam, lt, i, j):
    return 2*dot(bisp(lam[i-1], lt[i-1]), bisp(lam[j-1], lt[j-1]))
def loop_rows_allplus(recs, x):
    rows, rhs = [], []
    for rec in recs:
        if not rec['ok'][(1, 1, 1, 1)]: continue
        lam, lt = rec['lam'], rec['lt']
        row = []
        for s_ in sig:
            A = CxA1(lam, lt, s_, x, None)
            row += [sval(lam, lt, i, j)*A for (i, j) in basis]
        rows.append(row); rhs.append(rec['M'][(1, 1, 1, 1)])
    return rows, rhs
def Q(M): return DomainMatrix([[sp.Rational(v.numerator, v.denominator) for v in r] for r in M], (len(M), len(M[0])), sp.QQ)
if __name__ == '__main__':
    recs = sum([pickle.load(open(f, 'rb')) for f in glob.glob('exact_n4_*.pkl')], [])
    print('EYM n=4 records', len(recs))
    for x in [Fr(1, 3), Fr(2, 7)]:
        with Pool(12) as pool:
            res = pool.map(tree_rows, [(1000+k, x) for k in range(12)])
        T = sum([r[0] for r in res], []); tb = sum([r[1] for r in res], [])
        L, lb = loop_rows_allplus(recs, x)
        rT = Q(T).rank(); rTb = Q(T).hstack(Q([[v] for v in tb])).rank()
        A = Q(T + L); Ab = A.hstack(Q([[v] for v in tb + lb]))
        print('x=%s: tree rank %d/%d ; tree+allplus rank %d, augmented %d (loop rows %d)' % (x, rT, rTb, A.rank(), Ab.rank(), len(L)))
