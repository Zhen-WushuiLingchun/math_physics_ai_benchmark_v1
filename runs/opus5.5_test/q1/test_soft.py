import time, sys
from eymloop import *
mp.mp.dps = 60
rng = random.Random(5)
n = 4
# base spinors
lam = [(rnd(rng), rnd(rng)) for _ in range(n)]
lt12 = [(rnd(rng), rnd(rng)) for _ in range(2)]
lamPh = (rnd(rng), rnd(rng)); ltP = (rnd(rng), rnd(rng))
def point(delta):
    lamP = (lamPh[0]*delta, lamPh[1]*delta)
    L = lam + [lamP]
    known = [(0, lt12[0]), (1, lt12[1]), (4, ltP)]
    a, b = 2, 3
    S = [sum(ang(L[b], L[i])*t[c] for i, t in known) for c in range(2)]
    T = [sum(ang(L[a], L[i])*t[c] for i, t in known) for c in range(2)]
    lta = tuple(-S[c]/ang(L[b], L[a]) for c in range(2)); ltb = tuple(-T[c]/ang(L[a], L[b]) for c in range(2))
    LT = [lt12[0], lt12[1], lta, ltb, ltP]
    assert all(x == 0 for x in vsum([bisp(L[i], LT[i]) for i in range(5)], Fr(0)))
    return L, LT
hel = (-1, 1, 1, 1)
refs = [(mp.mpc(rng.uniform(-2,2)), mp.mpc(rng.uniform(-2,2))) for _ in range(5)]
def legs_at(L, LT):
    gl = []
    for i in range(n):
        l, t = mpv(L[i]), mpv(LT[i])
        gl.append((bisp(l, t), eps_plus(l, t, refs[i]) if hel[i] > 0 else eps_minus(l, t, refs[i])))
    l, t = mpv(L[4]), mpv(LT[4])
    if L[4][0] == 0 and L[4][1] == 0:
        return gl, None, None
    return gl, bisp(l, t), eps_plus(l, t, refs[4])
L0, LT0 = point(Fr(0)) if False else (None, None)
for dl in [Fr(1, 10**4), Fr(1, 10**6)]:
    L, LT = point(dl)
    gl, pP, eP = legs_at(L, LT)
    Sft = sum(dot(eP, p)**2/dot(p, pP) for (p, e) in gl)
    # tree
    Mt = BG([('g',) + x for x in gl], grav=(eP, pP), zero=mp.mpc(0)).amp()
    # 4-pt gluon kinematics at delta->0: recompute with delta=0 exactly
    L0, LT0 = point(Fr(0))
    gl0, _, _ = legs_at(L0, LT0)
    At = BG([('g',) + x for x in gl0], zero=mp.mpc(0)).amp()
    print('delta', dl, ' tree ratio M/(S A) =', mp.nstr(Mt/(Sft*At), 12))
    t0 = time.time()
    M1, dg = eym_one_loop_split(gl, eP, pP)
    A1, parts, d1 = one_loop(gl0)
    print('   one-loop ratio M1/(S A1) =', mp.nstr(M1/(Sft*A1), 12), '(%.0fs)' % (time.time()-t0))
    print('   diag', [mp.nstr(max(d.values()), 2) for d in dg])
