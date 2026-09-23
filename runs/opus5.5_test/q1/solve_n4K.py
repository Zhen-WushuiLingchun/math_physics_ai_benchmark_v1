import pickle, time, sympy as sp
from symx import *
from sympy.polys.matrices import DomainMatrix
res = pickle.load(open('tree_n4_rows.pkl', 'rb'))
rows = sum([r[0] for r in res], []); rhs = sum([r[1] for r in res], [])
x = X
# pick independent rows at x = 1/3 (exact, fast)
num = [[sp.sympify(v) for v in r] for r in rows]; nb = [sp.sympify(v) for v in rhs]
at = lambda e: sp.Rational(e.subs(x, sp.Rational(2, 7)))
A0 = DomainMatrix([[at(v) for v in r] for r in num], (len(num), 60), sp.QQ)
_, piv = A0.transpose().rref()
sel = list(piv)
print('rank at x=2/7:', len(sel))
t0 = time.time()
A = DomainMatrix([[KX.from_sympy(num[i][j]) for j in range(60)] for i in sel], (len(sel), 60), KX)
b = DomainMatrix([[KX.from_sympy(nb[i])] for i in sel], (len(sel), 1), KX)
R, pv = A.hstack(b).rref()
print('rref time', time.time()-t0, 'consistent', 60 not in pv)
Rl = R.to_Matrix()
part = [0]*60
for r, p in enumerate(pv):
    if p < 60: part[p] = Rl[r, 60]
free = [c for c in range(60) if c not in pv]
null = []
for f in free:
    v = [0]*60; v[f] = 1
    for r, p in enumerate(pv):
        if p < 60: v[p] = -Rl[r, f]
    null.append(v)
print('null dim', len(null))
# verify all rows (incl. unused) at several x values exactly
ok = True
for xv in [sp.Rational(1, 5), sp.Rational(3, 11), sp.Rational(5, 8)]:
    P = [sp.Rational(sp.sympify(p).subs(x, xv)) for p in part]
    Ns = [[sp.Rational(sp.sympify(t).subs(x, xv)) for t in v] for v in null]
    for i in range(len(num)):
        ri = [v.subs(x, xv) for v in num[i]]
        if sum(ri[j]*P[j] for j in range(60)) != nb[i].subs(x, xv): ok = False
        for v in Ns:
            if sum(ri[j]*v[j] for j in range(60)) != 0: ok = False
print('all 240 rows satisfied at 3 x-values:', ok)
pickle.dump(([str(p) for p in part], [[str(t) for t in v] for v in null]), open('tree_n4.pkl', 'wb'))
from collinear import orderings
sig = orderings(4); basis = [(1, 2), (2, 3), (3, 4), (1, 3), (2, 4)]
labels = [(s_, bb) for s_ in sig for bb in basis]
print('particular solution nonzero entries:')
for l, v in zip(labels, part):
    if v != 0: print('  ', l, sp.factor(v))
