import time, sys
from eymloop import *
mp.mp.dps = 50
rng = random.Random(21)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
pt = Point(n, rng)
for hel in [(1,)*n, (-1,) + (1,)*(n-1)]:
    for rr in range(2):
        t0 = time.time()
        gl = [mp_gluon(pt.lam[i], pt.lt[i], hel[i], rng) for i in range(n)]
        pP, eP = mp_gluon(pt.lam[n], pt.lt[n], 1, rng)
        tot, terms, diags = eym_one_loop_split(gl, eP, pP, seed=rr+1, return_terms=True)
        print(hel, 'M1 =', mp.nstr(tot, 25), '(%.1fs)' % (time.time()-t0))
        print('   terms', [mp.nstr(t, 10) for t in terms])
        print('   diag', [{k: mp.nstr(v, 2) for k, v in d.items()} for d in diags])
