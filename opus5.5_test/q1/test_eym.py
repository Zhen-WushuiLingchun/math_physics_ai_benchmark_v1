from kin import *
rng = random.Random(7)

def eym(lam, lt, hel, hP, refs=None, refP=None):
    """gluons 0..n-1 with helicities hel, graviton = last spinor with helicity hP (+-2)"""
    n = len(hel)
    legs = gluon_legs(lam[:n], lt[:n], hel, rng, refs)
    rP = refP if refP else (rnd(rng), rnd(rng))
    eP = pol(1 if hP > 0 else -1, lam[n], lt[n], rP)
    return BG(legs, grav=(eP, bisp(lam[n], lt[n])), zero=Fr(0)).amp(), legs, eP

def insertion(lam, lt, hel, legs, eP):
    n = len(hel)
    pP = bisp(lam[n], lt[n])
    tot = Fr(0)
    X = (Fr(0),)*4
    for l in range(1, n):
        X = add(X, legs[l-1][1])
        L = legs[:l] + [('g', pP, eP)] + legs[l:]
        tot += dot(eP, X)*BG(L, zero=Fr(0)).amp()
    return tot

for n in [3, 4, 5]:
    lam, lt = massless_kin(n+1, rng)
    for hel, hP in [([-1, -1] + [1]*(n-2), 2), ([-1] + [1]*(n-1), -2), ([1]*n, 2), ([-1]+[1]*(n-1), 2), ([-1, 1, -1] + [1]*(n-3), 2)]:
        vals = []
        for t in range(3):
            M, legs, eP = eym(lam, lt, hel, hP)
            I = insertion(lam, lt, hel, legs, eP)
            vals.append((M, M/I if I != 0 else None))
        print(n, hel, hP, vals)
    # gluon gauge invariance
    hel = [-1, -1] + [1]*(n-2)
    legs = gluon_legs(lam[:n], lt[:n], hel, rng)
    legs[1] = ('g', legs[1][1], legs[1][1])
    eP = pol(1, lam[n], lt[n], (rnd(rng), rnd(rng)))
    print('  gluon gauge:', BG(legs, grav=(eP, bisp(lam[n], lt[n])), zero=Fr(0)).amp())
