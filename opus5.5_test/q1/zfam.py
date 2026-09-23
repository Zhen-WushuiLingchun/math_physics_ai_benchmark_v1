import pickle, sys
from exactpt import *
mp.mp.dps = 90
rng = random.Random(int(sys.argv[1]))
zs = eval(sys.argv[2])
hel = eval(sys.argv[3]) if len(sys.argv) > 3 else (-1,1,1)
pt0 = small_point(3, random.Random(3))
out = []
for z in zs:
    z = Fr(z)
    lam = list(pt0.lam); lt = list(pt0.lt)
    lt[0] = (lt[0][0] + z*lt[1][0], lt[0][1] + z*lt[1][1])
    lam[1] = (lam[1][0] - z*lam[0][0], lam[1][1] - z*lam[0][1])
    gl = [mp_gluon(lam[i], lt[i], h, rng) for i, h in enumerate(hel)]
    pP, eP = mp_gluon(lam[3], lt[3], 1, rng)
    tot = eym_one_loop_split(gl, eP, pP)[0]
    r = ratrec(tot, 10**35)
    out.append((z, lam, lt, r, abs(tot - mp.mpf(r.numerator)/r.denominator) < mp.mpf(10)**-70*max(1,abs(tot))))
    print(z, r, out[-1][-1], flush=True)
pickle.dump(out, open('zfam_%s.pkl' % sys.argv[1], 'wb'))
