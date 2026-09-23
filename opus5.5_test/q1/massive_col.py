"""tree-level collinear relation with two massive scalars (D-dim cut trees)."""
from colK import *
from sympy.polys.matrices import DomainMatrix
def massive_pointK(k, rng):
    """k gluons + P + scalar s (l) + sb, exact, P spinors given. returns spinors (k gluons + P), l, lb, m2 as KX"""
    while True:
        p = Point(max(k+2, 3), rng)   # generic massless set, only first k+1 used (k gluons + P)
        lam, lt = p.lam[:k+1], p.lt[:k+1]
        moms = [bisp(lam[i], lt[i]) for i in range(k+1)]
        K = vsum(moms, Fr(0))
        l11, l12, l21 = rnd(rng), rnd(rng), rnd(rng)
        if K[0] == 0: continue
        l22 = (-dot(K, K) - l11*K[3] + l12*K[2] + l21*K[1])/K[0]
        l = (l11, l12, l21, l22)
        m2 = dot(l, l)
        if m2 == 0: continue
        cv = lambda v: tuple(toK(c) for c in v)
        return [cv(a) for a in lam], [cv(a) for a in lt], cv(l), cv(neg(add(l, K))), toK(m2)

def eym_massive(lam, lt, l, lb, m2, hel, hP, rng):
    k = len(hel)
    legs = [('s', l, 1)] + [('g', bisp(lam[i], lt[i]), polK(hel[i], lam[i], lt[i], (toK(rnd(rng)), toK(rnd(rng))))) for i in range(k)] + [('sb', lb, 1)]
    eP = polK(1 if hP > 0 else -1, lam[k], lt[k], (toK(rnd(rng)), toK(rnd(rng))))
    return BG(legs, grav=(eP, bisp(lam[k], lt[k])), m2=m2, zero=KX(0)).amp()

def ym_massive_col(lam, lt, l, lb, m2, seq, hel, hab, rng):
    """seq: list over 's','sb', gluon labels 1..k, 'a','b' -- must start with 's' and end with 'sb'"""
    k = len(hel)
    lamP, ltP = lam[k], lt[k]
    legs = []
    for lab in seq:
        if lab == 's': legs.append(('s', l, 1)); continue
        if lab == 'sb': legs.append(('sb', lb, 1)); continue
        if lab == 'a': la, t, h = lamP, (xK*ltP[0], xK*ltP[1]), hab
        elif lab == 'b': la, t, h = lamP, ((1-xK)*ltP[0], (1-xK)*ltP[1]), hab
        else: la, t, h = lam[lab-1], lt[lab-1], hel[lab-1]
        legs.append(('g', bisp(la, t), polK(h, la, t, (toK(rnd(rng)), toK(rnd(rng))))))
    return BG(legs, m2=m2, zero=KX(0)).amp() * xK**(-hab) * (1-xK)**(-hab)
