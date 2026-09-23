import pickle, sympy as sp
from colK import *
from ym_formulas import *
from eym_formulas import *
from sympy.polys.matrices import DomainMatrix
sig, basis, part, null = pickle.load(open('tree_n3.pkl', 'rb'))
labels = [(s_, b) for s_ in sig for b in basis]
Kp = [KX.from_sympy(sp.sympify(v)) for v in part]
Nv = [[KX.from_sympy(sp.sympify(v)) for v in vec] for vec in null]
rng = random.Random(2025)
from exactpt import small_point
rows, rhs, info = [], [], []
for k in range(6):
    p = small_point(3, rng)
    lam = [tuple(toK(c) for c in l) for l in p.lam]; lt = [tuple(toK(c) for c in l) for l in p.lt]
    sv = {}
    for (i, j) in basis:
        sv[(i, j)] = 2*dot(bisp(lam[i-1], lt[i-1]), bisp(lam[j-1], lt[j-1]))
    for minus in [None, 1, 2, 3]:
        CA = {s_: CxA1(lam, lt, s_, xK, minus) for s_ in sig}
        col = [sv[b]*CA[s_] for (s_, b) in labels]          # coefficient of each kernel coefficient
        M1 = KX(0) if minus is None else M1_sm_n3(lam, lt, minus-1)
        tree_part = sum(Kp[i]*col[i] for i in range(len(col)))
        Delta = M1 - tree_part
        rows.append([sum(v[i]*col[i] for i in range(len(col))) for v in Nv])
        rhs.append(Delta)
        info.append((k, minus, Delta))
pickle.dump(([[str(fromK(v)) for v in r] for r in rows], [str(fromK(v)) for v in rhs]), open('oneloop_n3_sys.pkl', 'wb'))
for k, minus, D in info[:8]:
    print(k, minus, sp.factor(fromK(D)))
A = DomainMatrix([list(r) for r in rows], (len(rows), len(Nv)), KX)
b = DomainMatrix([[v] for v in rhs], (len(rhs), 1), KX)
print('rank V =', A.rank(), ' rank [V|Delta] =', A.hstack(b).rank(), ' (#null dirs = %d, #eqs = %d)' % (len(Nv), len(rows)))
