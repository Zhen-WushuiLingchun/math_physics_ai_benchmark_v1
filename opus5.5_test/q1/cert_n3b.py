import pickle, sympy as sp
from symx import *
from exactpt import *
from eym_formulas import *
from sympy.polys.matrices import DomainMatrix
rows, rhs = pickle.load(open('oneloop_n3_sys.pkl', 'rb'))
x = X
V = sp.Matrix([[sp.sympify(v) for v in rows[i]] for i in (1, 2, 3)])
D = sp.Matrix([sp.sympify(rhs[i]) for i in (1, 2, 3)])
print('rank V(3x8) =', V.rank(simplify=True))
ns = V.T.nullspace(simplify=True)
y = sp.simplify(ns[0])
y = sp.simplify(y/y[0])
print('y =', [sp.factor(t) for t in y])
print('y.V =', sp.simplify((y.T*V)))
cert = sp.factor(sp.simplify((y.T*D)[0]))
print('y.Delta =', cert)
# point used
rng = random.Random(2025)
p = small_point(3, rng)
print('point lam =', p.lam, '\n      lt =', p.lt)
L = [tuple(sp.Rational(c.numerator, c.denominator) for c in l) for l in p.lam]
T = [tuple(sp.Rational(c.numerator, c.denominator) for c in l) for l in p.lt]
print('formula M1 (minus on 1,2,3):', [M1_sm_n3(L, T, m) for m in range(3)])
# engine check
mp.mp.dps = 90
eng = []
for m in range(3):
    hel = [1, 1, 1]; hel[m] = -1
    gl = [mp_gluon(p.lam[i], p.lt[i], hel[i], rng) for i in range(3)]
    pP, eP = mp_gluon(p.lam[3], p.lt[3], 1, rng)
    tot = eym_one_loop_split(gl, eP, pP)[0]
    eng.append(ratrec(tot, 10**35))
print('engine M1 (rational reconstruction):', eng)
