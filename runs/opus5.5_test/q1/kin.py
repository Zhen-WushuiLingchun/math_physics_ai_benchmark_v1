from fractions import Fraction as Fr
import random
from trees import *

def rnd(rng, lo=-29, hi=29):
    while True:
        v = rng.randint(lo, hi)
        if v != 0:
            return Fr(v, rng.randint(1, 7))

def massless_kin(N, rng):
    """N massless momenta (all incoming, sum=0), real rational spinors (split signature)."""
    while True:
        try:
            lam, lt = _massless_kin(N, rng)
        except ZeroDivisionError:
            continue
        if all(ang(lam[i], lam[j]) != 0 and sqb(lt[i], lt[j]) != 0 for i in range(N) for j in range(i+1, N)):
            return lam, lt

def _massless_kin(N, rng):
    lam = [(rnd(rng), rnd(rng)) for _ in range(N)]
    lt = [(rnd(rng), rnd(rng)) for _ in range(N-2)]
    # solve for lt[N-2], lt[N-1]
    a, b = N-2, N-1
    S = [sum(ang(lam[b], lam[i])*lt[i][c] for i in range(N-2)) for c in range(2)]
    T = [sum(ang(lam[a], lam[i])*lt[i][c] for i in range(N-2)) for c in range(2)]
    lta = tuple(-S[c]/ang(lam[b], lam[a]) for c in range(2))
    ltb = tuple(-T[c]/ang(lam[a], lam[b]) for c in range(2))
    lt = lt + [lta, ltb]
    moms = [bisp(lam[i], lt[i]) for i in range(N)]
    tot = vsum(moms, Fr(0))
    assert all(x == 0 for x in tot), tot
    return lam, lt

def pol(h, lam, lt, ref):
    """h=+1: ref is a lambda;  h=-1: ref is a lambda-tilde"""
    if (ang(ref, lam) if h > 0 else sqb(lt, ref)) == 0:
        ref = (ref[0] + 1, ref[1] + Fr(1, 3))
    return eps_plus(lam, lt, ref) if h > 0 else eps_minus(lam, lt, ref)

def gluon_legs(lam, lt, hel, rng, refs=None):
    legs = []
    for i, h in enumerate(hel):
        r = refs[i] if refs else (rnd(rng), rnd(rng))
        legs.append(('g', bisp(lam[i], lt[i]), pol(h, lam[i], lt[i], r)))
    return legs
