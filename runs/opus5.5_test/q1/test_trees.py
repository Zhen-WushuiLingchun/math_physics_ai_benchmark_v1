from kin import *
rng = random.Random(1)

def PT(lam, i, j, N):
    den = 1
    for k in range(N):
        den *= ang(lam[k], lam[(k+1) % N])
    return ang(lam[i], lam[j])**4/den

for N in [4, 5, 6]:
    lam, lt = massless_kin(N, rng)
    hel = [-1, -1] + [1]*(N-2)
    rs = []
    for trial in range(3):
        legs = gluon_legs(lam, lt, hel, rng)
        A = BG(legs, zero=Fr(0)).amp()
        rs.append(A/PT(lam, 0, 1, N))
    print('N=%d MHV A/PT ratios:' % N, rs)
    # gauge invariance
    legs = gluon_legs(lam, lt, hel, rng)
    legs[2] = ('g', legs[2][1], legs[2][1])
    print('   gauge check:', BG(legs, zero=Fr(0)).amp())
    # all plus & one minus vanish
    for hh in ([1]*N, [-1]+[1]*(N-1)):
        legs = gluon_legs(lam, lt, hh, rng)
        print('   hel', hh, BG(legs, zero=Fr(0)).amp())
