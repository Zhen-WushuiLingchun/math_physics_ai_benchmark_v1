"""Direct n=4 cross-check at fixed x0 (numerical, 40 digits): tree affine space (exact at x0) +
one-loop equations (EYM via split representation, YM collinear via engine / all-plus formula)."""
import sys, pickle, itertools
from exactpt import *
from opp import Topology, Reducer, MU2S, eym_one_loop_split
from ym_formulas import *
from multiprocessing import Pool
X0 = Fr(2, 7)
n = 4
sig = orderings(n)
basis = [(1, 2), (2, 3), (3, 4), (1, 3), (2, 4)]
SECT = [None, 1, 2, 3, 4]

def ym_block(args):
    seed, si, m = args
    mp.mp.dps = 40
    rng = random.Random(seed); pt = Point(n, rng)
    sigma = sig[si]; lamP, ltP = pt.lam[n], pt.lt[n]; x = X0; xm = mpcf(x)
    hel = [1]*n
    if m is not None:
        hel[m-1] = -1
    if m is None:
        L = [mpv(l) for l in pt.lam]; T = [mpv(t) for t in pt.lt]
        return (seed, si, m, CxA1(L, T, sigma, xm, None))
    gl = []
    for lab in sigma:
        if lab == 'a': gl.append(mp_gluon(lamP, (ltP[0]*x, ltP[1]*x), 1, rng))
        elif lab == 'b': gl.append(mp_gluon(lamP, (ltP[0]*(1-x), ltP[1]*(1-x)), 1, rng))
        else: gl.append(mp_gluon(pt.lam[lab-1], pt.lt[lab-1], hel[lab-1], rng))
    red = Reducer(Topology(gl), MU2S); red.run()
    return (seed, si, m, sum(red.contributions_general().values())/(xm*(1-xm)))

def eym_block(args):
    seed, m = args
    mp.mp.dps = 40
    rng = random.Random(seed); pt = Point(n, rng)
    hel = [1]*n
    if m is not None:
        hel[m-1] = -1
    gl = [mp_gluon(pt.lam[i], pt.lt[i], hel[i], rng) for i in range(n)]
    pP, eP = mp_gluon(pt.lam[n], pt.lt[n], 1, rng)
    return (seed, m, eym_one_loop_split(gl, eP, pP)[0])

if __name__ == '__main__':
    seeds = [int(a) for a in sys.argv[1:]]
    jobs_y = [(sd, si, m) for sd in seeds for si in range(len(sig)) for m in SECT]
    jobs_e = [(sd, m) for sd in seeds for m in SECT]
    with Pool(8) as pool:
        ye = pool.map(eym_block, jobs_e)
        yy = pool.map(ym_block, jobs_y)
    out = {'eym': [(a, b, str(c.real), str(c.imag)) for a, b, c in ye],
           'ym': [(a, b, c, str(d.real), str(d.imag)) for a, b, c, d in yy]}
    pickle.dump(out, open('n4_direct_%d.pkl' % seeds[0], 'wb'))
    print('done', len(ye), len(yy))
