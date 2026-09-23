"""symbolic proof check: K^tree and the 8 null vectors for n=3 hold as rational-function identities
   in free spinor variables (momentum conservation solved) and x."""
import pickle, itertools, time, sympy as sp
from trees import *
from collinear import orderings
# Lorentz + little-group frame: lam1=(1,0), lam2=(0,1), lam3=(1,1), lamP=(1,w) (Zariski-open, covariant identities)
syms = sp.symbols('w e1 e2 f1 f2 x')
F = sp.QQ.frac_field(*syms)
g = F.gens
lam = [(F(1), F(0)), (F(0), F(1)), (F(1), F(1)), (F(1), g[0])]
lt12 = [(g[1], g[2]), (g[3], g[4])]; xK = g[5]
# solve lt3, ltP from sum_i lam_i lt_i = 0
a, b = 2, 3
S = [sum(ang(lam[b], lam[i])*lt12[i][c] for i in range(2)) for c in range(2)]
T = [sum(ang(lam[a], lam[i])*lt12[i][c] for i in range(2)) for c in range(2)]
lt = lt12 + [tuple(-S[c]/ang(lam[b], lam[a]) for c in range(2)), tuple(-T[c]/ang(lam[a], lam[b]) for c in range(2))]
assert all(v == 0 for v in vsum([bisp(lam[i], lt[i]) for i in range(4)], F(0)))
refs = [(F(3), F(-7)), (F(5), F(2)), (F(-4), F(9)), (F(7), F(11)), (F(-2), F(5)), (F(13), F(-3))]
def pol(h, la, t, r):
    return eps_plus(la, t, r) if h > 0 else eps_minus(la, t, r)
sig, basis, part, null = pickle.load(open('tree_n3.pkl', 'rb'))
X = sp.Symbol('x')
conv = lambda e: F.from_sympy(sp.sympify(e).subs(X, syms[5]))
Kp = [conv(v) for v in part]; Nv = [[conv(t) for t in v] for v in null]
mom = [bisp(lam[i], lt[i]) for i in range(4)]
sv = {(i, j): 2*dot(mom[i-1], mom[j-1]) for (i, j) in basis}
t0 = time.time(); bad = 0
for hP in (2, -2):
    for hel in itertools.product([1, -1], repeat=3):
        legs = [('g', mom[i], pol(hel[i], lam[i], lt[i], refs[i])) for i in range(3)]
        eP = pol(1 if hP > 0 else -1, lam[3], lt[3], refs[3])
        M = BG(legs, grav=(eP, mom[3]), zero=F(0)).amp()
        hab = 1 if hP > 0 else -1
        col = []
        for s_ in sig:
            L = []
            for lab in s_:
                if lab == 'a': la, t, h, r = lam[3], (xK*lt[3][0], xK*lt[3][1]), hab, refs[4]
                elif lab == 'b': la, t, h, r = lam[3], ((1-xK)*lt[3][0], (1-xK)*lt[3][1]), hab, refs[5]
                else: la, t, h, r = lam[lab-1], lt[lab-1], hel[lab-1], refs[lab-1]
                L.append(('g', bisp(la, t), pol(h, la, t, r)))
            A = BG(L, zero=F(0)).amp()*xK**(-hab)*(1-xK)**(-hab)
            col += [sv[bb]*A for bb in basis]
        ok1 = sum(Kp[i]*col[i] for i in range(12)) == M
        ok2 = all(sum(v[i]*col[i] for i in range(12)) == 0 for v in Nv)
        bad += (not ok1) + (not ok2)
        print(hP, hel, 'kernel identity:', ok1, ' null identities:', ok2, '(%.0fs)' % (time.time()-t0), flush=True)
print('FAILURES', bad)
