import pickle, sys, time
from exactpt import *
mp.mp.dps = 90
seed = int(sys.argv[1]); npts = int(sys.argv[2])
rng = random.Random(seed)
hels = [(1,1,1,1), (-1,1,1,1), (1,-1,1,1), (1,1,-1,1), (1,1,1,-1)]
out = []
for k in range(npts):
    pt = small_point(4, rng, R=3)
    rec = {'lam': pt.lam, 'lt': pt.lt, 'M': {}, 'ok': {}}
    for hel in hels:
        gl = [mp_gluon(pt.lam[i], pt.lt[i], h, rng) for i, h in enumerate(hel)]
        pP, eP = mp_gluon(pt.lam[4], pt.lt[4], 1, rng)
        tot = eym_one_loop_split(gl, eP, pP)[0]
        r = ratrec(tot, 10**38)
        rec['M'][hel] = r
        rec['ok'][hel] = abs(tot - mp.mpf(r.numerator)/r.denominator) < mp.mpf(10)**-75*max(1, abs(tot))
    out.append(rec)
    print(k, rec['M'], rec['ok'], flush=True)
    pickle.dump(out, open('exact_n4_%d.pkl' % seed, 'wb'))
