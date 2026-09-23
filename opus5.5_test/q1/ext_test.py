import pickle, sys, sympy as sp
from exactpt import *
from opp import Topology, Reducer, MU2S
from ym_formulas import *
from eym_formulas import *
from multiprocessing import Pool
mp.mp.dps = 50
X0 = Fr(2, 7)
sig = orderings(3); basis = [(1, 2), (2, 3)]
def block(args):
    seed, s_idx, m, kind = args
    rng = random.Random(seed); pt = Point(3, rng)
    sigma = sig[s_idx]; n = 3; lamP, ltP = pt.lam[n], pt.lt[n]; x = X0; xm = mpcf(x)
    hel = [1, 1, 1]
    if m is not None: hel[m-1] = -1
    gl = []
    for lab in sigma:
        if lab == 'a': gl.append(mp_gluon(lamP, (ltP[0]*x, ltP[1]*x), 1, rng))
        elif lab == 'b': gl.append(mp_gluon(lamP, (ltP[0]*(1-x), ltP[1]*(1-x)), 1, rng))
        else: gl.append(mp_gluon(pt.lam[lab-1], pt.lt[lab-1], hel[lab-1], rng))
    if kind == 'mu':
        topo = Topology(gl, mu2pow=1)
    else:
        ia = sigma.index('a')
        q = (mp.mpc(0),)*4
        for g in gl[:ia]: q = add(q, g[0])
        P = bisp(mpv(lamP), mpv(ltP))
        topo = Topology(gl, weight=(mp.mpf(2), P, q))
    red = Reducer(topo, MU2S); red.run()
    v = sum(red.contributions_general().values())/(xm*(1-xm))
    return (seed, s_idx, m, kind, v)
if __name__ == '__main__':
    import sys
    seeds = [int(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else [500+k for k in range(7)]
    jobs = [(sd, si, m, kd) for sd in seeds for si in range(6) for m in (None, 1, 2, 3) for kd in ('mu', 'l')]
    with Pool(8) as pool:
        res = pool.map(block, jobs)
    pickle.dump([(a, b, c, d, str(e.real), str(e.imag)) for a, b, c, d, e in res], open('ext_blocks%s.pkl' % ('_more' if len(sys.argv) > 1 else ''), 'wb'))
    print('done', len(res))
