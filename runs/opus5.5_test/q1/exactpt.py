from eymloop import *
from fractions import Fraction as Fr
def small_point(n, rng, R=4):
    """n gluons + P, small integer lambdas, integer-ish lt"""
    while True:
        N = n+1
        lam = [(Fr(rng.randint(-R, R)), Fr(rng.randint(-R, R))) for _ in range(N)]
        lt = [(Fr(rng.randint(-R, R)), Fr(rng.randint(-R, R))) for _ in range(N-2)]
        try:
            a, b = N-2, N-1
            S = [sum(ang(lam[b], lam[i])*lt[i][c] for i in range(N-2)) for c in range(2)]
            T = [sum(ang(lam[a], lam[i])*lt[i][c] for i in range(N-2)) for c in range(2)]
            lt = lt + [tuple(-S[c]/ang(lam[b], lam[a]) for c in range(2)), tuple(-T[c]/ang(lam[a], lam[b]) for c in range(2))]
        except ZeroDivisionError:
            continue
        if all(ang(lam[i], lam[j]) != 0 and sqb(lt[i], lt[j]) != 0 for i in range(N) for j in range(i+1, N)):
            p = Point(n, rng); p.lam, p.lt = lam, lt
            p.mom = [bisp(lam[i], lt[i]) for i in range(N)]
            return p
def ratrec(x, maxden=10**40):
    f = Fr(str(mp.nstr(x.real, mp.mp.dps-5))).limit_denominator(maxden)
    return f
