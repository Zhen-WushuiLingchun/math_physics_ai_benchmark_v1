import pickle, sys
from exactpt import *
mp.mp.dps = 90
seed = int(sys.argv[1]); npts = int(sys.argv[2])
hels = eval(sys.argv[3]) if len(sys.argv) > 3 else [(-1,1,1)]
rng = random.Random(seed)
out = []
for k in range(npts):
    pt = small_point(3, rng)
    rec = {'lam': pt.lam, 'lt': pt.lt, 'M': {}}
    for hel in hels:
        gl = [mp_gluon(pt.lam[i], pt.lt[i], h, rng) for i, h in enumerate(hel)]
        pP, eP = mp_gluon(pt.lam[3], pt.lt[3], 1, rng)
        tot, dg = eym_one_loop_split(gl, eP, pP)
        r = ratrec(tot, 10**35)
        ok = abs(tot - mp.mpf(r.numerator)/r.denominator) < mp.mpf(10)**(-70)*max(1, abs(tot))
        rec['M'][hel] = r if ok else None
        rec.setdefault('num', {})[hel] = tot
    out.append(rec)
    print(k, rec['M'], flush=True)
pickle.dump(out, open('exact_n3_%d.pkl' % seed, 'wb'))
