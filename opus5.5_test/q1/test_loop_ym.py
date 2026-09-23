import time, sys
from opp import *
from kin import *
mp.mp.dps = 40

def to_mp_legs(lam, lt, hel, rng):
    out = []
    for i, h in enumerate(hel):
        l = tuple(mpcf(c) for c in lam[i]); t = tuple(mpcf(c) for c in lt[i])
        r = (mp.mpc(rng.uniform(-2, 2)), mp.mpc(rng.uniform(-2, 2)))
        e = eps_plus(l, t, r) if h > 0 else eps_minus(l, t, r)
        out.append((bisp(l, t), e))
    return out

def allplus_formula(lam, lt, n):
    L = [tuple(mpcf(c) for c in x) for x in lam]; T = [tuple(mpcf(c) for c in x) for x in lt]
    num = 0
    for a, b, c, d in itertools.combinations(range(n), 4):
        num += ang(L[a], L[b])*sqb(T[b], T[c])*ang(L[c], L[d])*sqb(T[d], T[a])
    den = 1
    for i in range(n):
        den *= ang(L[i], L[(i+1) % n])
    return num/den

rng = random.Random(3)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
for trial in range(2):
    lam, lt = massless_kin(N, rng)
    for rr in range(2):
        t0 = time.time()
        legs = to_mp_legs(lam, lt, [1]*N, rng)
        tot, parts, diag = one_loop(legs)
        F = allplus_formula(lam, lt, N)
        print('N=%d trial %d ref %d: A=%s  A/formula=%s  (%.1fs)' % (N, trial, rr, mp.nstr(tot, 15), mp.nstr(tot/F, 15), time.time()-t0))
        print('   parts', {k: mp.nstr(v, 8) for k, v in parts.items()})
        print('   diag', {k: mp.nstr(v, 3) for k, v in diag.items()})
