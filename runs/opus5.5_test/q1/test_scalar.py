from kin import *
rng = random.Random(11)

def massive_setup(k, grav):
    # k gluons (+graviton), scalar s with momentum l, sb with -(l+K); need (l+K)^2 = l^2
    N = k + (1 if grav else 0)
    lam, lt = massless_kin(N + 3, rng)   # use N massless momenta from an N+2 point set: K arbitrary
    moms = [bisp(lam[i], lt[i]) for i in range(N)]
    K = vsum(moms, Fr(0))
    # l = (l11,l12,l21,l22), choose three random, solve 2 l.K + K^2 = 0 for l22
    l11, l12, l21 = rnd(rng), rnd(rng), rnd(rng)
    # 2 l.K = l11 K22 + l22 K11 - l12 K21 - l21 K12
    l22 = (-dot(K, K) - l11*K[3] + l12*K[2] + l21*K[1])/K[0]
    l = (l11, l12, l21, l22)
    m2 = dot(l, l)
    lb = neg(add(l, K))
    assert dot(lb, lb) == m2
    return lam, lt, l, lb, m2

for k in [1, 2, 3, 4]:
    for grav in [False, True]:
        lam, lt, l, lb, m2 = massive_setup(k, grav)
        hel = [1, -1, 1, 1][:k]
        vals = []
        for t in range(3):
            legs = gluon_legs(lam[:k], lt[:k], hel, rng)
            L = [('s', l, 1)] + legs + [('sb', lb, 1)]
            gr = None
            if grav:
                eP = pol(1, lam[k], lt[k], (rnd(rng), rnd(rng)))
                gr = (eP, bisp(lam[k], lt[k]))
            A = BG(L, grav=gr, m2=m2, zero=Fr(0)).amp()
            ins = None
            if grav:
                tot = Fr(0); X = l
                for pos in range(k+1):
                    if pos > 0:
                        X = add(X, legs[pos-1][1])
                    LL = [('s', l, 1)] + legs[:pos] + [('g', gr[1], gr[0])] + legs[pos:] + [('sb', lb, 1)]
                    tot += dot(gr[0], X)*BG(LL, m2=m2, zero=Fr(0)).amp()
                ins = A/tot if tot != 0 else None
            vals.append((A, ins))
        # gauge
        legs = gluon_legs(lam[:k], lt[:k], hel, rng)
        legs[0] = ('g', legs[0][1], legs[0][1])
        L = [('s', l, 1)] + legs + [('sb', lb, 1)]
        g0 = BG(L, grav=gr, m2=m2, zero=Fr(0)).amp()
        print(k, grav, [v[0] == vals[0][0] for v in vals], [v[1] for v in vals], 'gauge:', g0)
