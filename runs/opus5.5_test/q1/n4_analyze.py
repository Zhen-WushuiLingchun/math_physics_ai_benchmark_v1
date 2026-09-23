import pickle, glob, sympy as sp
from exactpt import *
from fixedx_n4 import tree_rows, sig, basis, n
from multiprocessing import Pool
from sympy.polys.matrices import DomainMatrix
X0 = Fr(2, 7)
if __name__ == '__main__':
    mp.mp.dps = 40
    with Pool(8) as pool:
        res = pool.map(tree_rows, [(1000+k, X0) for k in range(12)])
    T = sum([r[0] for r in res], []); tb = sum([r[1] for r in res], [])
    A = DomainMatrix([[sp.Rational(v.numerator, v.denominator) for v in r] + [sp.Rational(b.numerator, b.denominator)]
                      for r, b in zip(T, tb)], (len(T), 61), sp.QQ)
    R, piv = A.rref(); Rl = R.to_Matrix()
    assert 60 not in piv
    part = [0]*60
    for r, p in enumerate(piv): part[p] = Rl[r, 60]
    free = [c for c in range(60) if c not in piv]
    null = []
    for f in free:
        v = [0]*60; v[f] = 1
        for r, p in enumerate(piv): v[p] = -Rl[r, f]
        null.append(v)
    print('tree: rank', len(piv), 'null dim', len(null))
    tomp = lambda q: mp.mpf(q.p)/q.q
    Kp = [tomp(sp.Rational(v)) for v in part]; Nv = [[tomp(sp.Rational(t)) for t in v] for v in null]
    data = [pickle.load(open(f, 'rb')) for f in glob.glob('n4_direct_*.pkl')]
    eym = {}; ym = {}
    for d in data:
        for a, b, re, im in d['eym']: eym[(a, b)] = mp.mpc(mp.mpf(re), mp.mpf(im))
        for a, b, c, re, im in d['ym']: ym[(a, b, c)] = mp.mpc(mp.mpf(re), mp.mpf(im))
    seeds = sorted(set(k[0] for k in eym))
    rows, rhs = [], []
    for sd in seeds:
        pt = Point(n, random.Random(sd))
        L = [mpv(l) for l in pt.lam]; Tt = [mpv(t) for t in pt.lt]
        sv = lambda i, j: 2*dot(bisp(L[i-1], Tt[i-1]), bisp(L[j-1], Tt[j-1]))
        for m in [None, 1, 2, 3, 4]:
            col = [sv(i, j)*ym[(sd, si, m)] for si in range(len(sig)) for (i, j) in basis]
            rhs.append(eym[(sd, m)] - sum(Kp[i]*col[i] for i in range(60)))
            rows.append([sum(v[i]*col[i] for i in range(60)) for v in Nv])
    def svals(Rr):
        M = mp.matrix(Rr)
        for i in range(M.rows):
            nr = max(abs(M[i, j]) for j in range(M.cols))
            for j in range(M.cols): M[i, j] /= nr
        return mp.svd_c(M, compute_uv=False)
    s1 = svals(rows); s2 = svals([r + [b] for r, b in zip(rows, rhs)])
    print('one-loop eqs', len(rows), ' unknowns', len(Nv))
    print('V  :', [mp.nstr(s1[i]/s1[0], 2) for i in range(len(s1))])
    print('aug:', [mp.nstr(s2[i]/s2[0], 2) for i in range(len(s2))])
    # Delta of the RREF representative
    print('max |Delta| (rep.):', mp.nstr(max(abs(b) for b in rhs), 5))
    import json
    json.dump({'V': [float(s1[i]/s1[0]) for i in range(len(s1))], 'aug': [float(s2[i]/s2[0]) for i in range(len(s2))], 'neq': len(rows)}, open('paper/spectra_n4.json', 'w'), indent=1)
