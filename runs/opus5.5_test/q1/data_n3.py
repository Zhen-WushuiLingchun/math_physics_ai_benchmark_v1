import time, sys, pickle
from eymloop import *
mp.mp.dps = 50
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 100)
npts = int(sys.argv[2]) if len(sys.argv) > 2 else 6
out = []
for k in range(npts):
    pt = Point(3, rng)
    rec = {'lam': pt.lam, 'lt': pt.lt, 'M': {}}
    for hel in [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]:
        gl = [mp_gluon(pt.lam[i], pt.lt[i], hel[i], rng) for i in range(3)]
        pP, eP = mp_gluon(pt.lam[3], pt.lt[3], 1, rng)
        tot, dg = eym_one_loop_split(gl, eP, pP)
        rec['M'][hel] = tot
        worst = max(max(d['tri_hiL'], d['bub_hiz'], d['high_mu']) for d in dg)
        print(k, hel, mp.nstr(tot, 30), 'hi-mode diag', mp.nstr(worst, 2), flush=True)
    out.append(rec)
pickle.dump(out, open('data_n3_%s.pkl' % (sys.argv[1] if len(sys.argv) > 1 else '100'), 'wb'))
