from exactpt import *
from opp import Topology, Reducer, MU2S
from eym_formulas import *
import sys
mp.mp.dps = 60
def ym_col_split(pt, sigma, hel, x, rng):
    """C_x A^(1)(sigma) with a^+b^+, split into (a,b same tree) and (a,b separated) masters"""
    n = pt.n; lamP, ltP = pt.lam[n], pt.lt[n]; xm = mpcf(x)
    gl = []
    for lab in sigma:
        if lab == 'a': gl.append(mp_gluon(lamP, (ltP[0]*x, ltP[1]*x), 1, rng))
        elif lab == 'b': gl.append(mp_gluon(lamP, (ltP[0]*(1-x), ltP[1]*(1-x)), 1, rng))
        else: gl.append(mp_gluon(pt.lam[lab-1], pt.lt[lab-1], hel[lab-1], rng))
    ia, ib = sigma.index('a'), sigma.index('b')
    topo = Topology(gl); red = Reducer(topo, MU2S); red.run()
    same = sep = mp.mpc(0)
    for S, v in red.contributions().items():
        together = any(ia in g and ib in g for g, hp in topo.segments(S))
        if together: same += v
        else: sep += v
    f = 1/(xm*(1-xm))
    return same*f, sep*f
if __name__ == '__main__':
    rng = random.Random(2025)
    p = small_point(3, rng)
    x = Fr(2, 7); xm = mpcf(x)
    s13 = 2*dot(bisp(mpv(p.lam[0]), mpv(p.lt[0])), bisp(mpv(p.lam[2]), mpv(p.lt[2])))
    K = -(xm*(1-xm)/4)*s13
    L = [mpv(l) for l in p.lam]; T = [mpv(t) for t in p.lt]
    for m in range(3):
        hel = [1, 1, 1]; hel[m] = -1
        same, sep = ym_col_split(p, (1, 'a', 2, 'b', 3), hel, x, rng)
        M = M1_sm_n3(L, T, m)
        print('minus', m+1, ' M1 =', mp.nstr(M, 15), ' K*same =', mp.nstr(K*same, 15), ' K*sep =', mp.nstr(K*sep, 15), ' M - K*same =', mp.nstr(M - K*same, 10))
