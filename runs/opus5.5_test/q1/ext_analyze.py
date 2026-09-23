import pickle, sympy as sp
from exactpt import *
from ym_formulas import *
from eym_formulas import *
mp.mp.dps = 50
X0 = Fr(2, 7); x0 = sp.Rational(2, 7)
sig, basis, part, null = pickle.load(open('tree_n3.pkl', 'rb'))
Kp = [mpcf(Fr(str(sp.sympify(v).subs(sp.Symbol('x'), x0)))) for v in part]
Nv = [[mpcf(Fr(str(sp.sympify(t).subs(sp.Symbol('x'), x0)))) for t in vec] for vec in null]
blk = {}
import glob
for a, b, c, d, re, im in sum([pickle.load(open(f, 'rb')) for f in glob.glob('ext_blocks*.pkl')], []):
    blk[(a, b, c, d)] = mp.mpc(mp.mpf(re), mp.mpf(im))
seeds = sorted(set(k[0] for k in blk))
rows = {'mu': [], 'l': [], 'both': [], 'none': []}; rhs = []
for sd in seeds:
    pt = Point(3, random.Random(sd))
    L = [mpv(l) for l in pt.lam]; T = [mpv(t) for t in pt.lt]
    sv = lambda i, j: 2*dot(bisp(L[i-1], T[i-1]), bisp(L[j-1], T[j-1]))
    for m in (None, 1, 2, 3):
        col = []
        for s_ in sig:
            A = CxA1(L, T, s_, mpcf(X0), m)
            col += [sv(i, j)*A for (i, j) in basis]
        M = mp.mpc(0) if m is None else M1_sm_n3(L, T, m-1)
        rhs.append(M - sum(Kp[i]*col[i] for i in range(12)))
        V = [sum(v[i]*col[i] for i in range(12)) for v in Nv]
        Bmu = [sv(i, j)*blk[(sd, si, m, 'mu')] for si in range(6) for (i, j) in basis]
        Bl = [sv(i, j)*blk[(sd, si, m, 'l')] for si in range(6) for (i, j) in basis]
        rows['none'].append(V); rows['mu'].append(V + Bmu); rows['l'].append(V + Bl); rows['both'].append(V + Bmu + Bl)
def svals(R):
    M = mp.matrix(R)
    for i in range(M.rows):
        nr = max(abs(M[i, j]) for j in range(M.cols))
        for j in range(M.cols): M[i, j] /= nr
    return mp.svd_c(M, compute_uv=False)
print('equations:', len(rhs))
for key in ['none', 'mu', 'l', 'both']:
    R = rows[key]; Rb = [r + [b] for r, b in zip(R, rhs)]
    s1 = svals(R); s2 = svals(Rb)
    rk = lambda s: sum(1 for v in s if v > mp.mpf(10)**-30*s[0])
    print('%-5s unknowns %2d : rank %2d, augmented rank %2d   (smallest sv of augmented: %s)' % (key, len(R[0]), rk(s1), rk(s2), mp.nstr(s2[len(s2)-1], 3)))
