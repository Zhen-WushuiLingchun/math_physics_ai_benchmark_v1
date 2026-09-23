import time, sys
from test_loop_ym import *
mp.mp.dps = 50
def br(L, T):
    A = lambda i, j: ang(L[i-1], L[j-1]); B = lambda i, j: sqb(T[i-1], T[j-1])
    return A, B
def f4(lam, lt):
    L = [tuple(mpcf(c) for c in x) for x in lam]; T = [tuple(mpcf(c) for c in x) for x in lt]
    A, B = br(L, T)
    return A(2,4)*B(2,4)**3/(B(1,2)*A(2,3)*A(3,4)*B(4,1))
def f5(lam, lt):
    L = [tuple(mpcf(c) for c in x) for x in lam]; T = [tuple(mpcf(c) for c in x) for x in lt]
    A, B = br(L, T)
    return (1/A(3,4)**2)*(-B(2,5)**3/(B(1,2)*B(5,1)) + A(1,4)**3*B(4,5)*A(3,5)/(A(1,2)*A(2,3)*A(4,5)**2)
                          - A(1,3)**3*B(3,2)*A(4,2)/(A(1,5)*A(5,4)*A(3,2)**2))
rng = random.Random(8)
for N, f in [(4, f4), (5, f5)]:
    for trial in range(2):
        lam, lt = massless_kin(N, rng)
        t0 = time.time()
        legs = to_mp_legs(lam, lt, [-1] + [1]*(N-1), rng)
        tot, parts, diag = one_loop(legs)
        print('N=%d: A=%s  A/formula=%s (%.1fs)' % (N, mp.nstr(tot, 12), mp.nstr(tot/f(lam, lt), 15), time.time()-t0))
        print('   parts', {k: mp.nstr(v, 8) for k, v in parts.items()}, '\n   diag', {k: mp.nstr(v, 3) for k, v in diag.items()})
