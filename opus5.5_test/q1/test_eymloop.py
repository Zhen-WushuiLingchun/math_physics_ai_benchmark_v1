import time, sys
from eymloop import *
mp.mp.dps = 50
rng = random.Random(21)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
pt = Point(n, rng)
for hel in [(1,)*n, (-1,) + (1,)*(n-1)]:
    for rr in range(2):
        t0 = time.time()
        tot, parts, diag = eym_one_loop(pt, hel, 2, rng, seed=rr+1)
        print(hel, 'M1 =', mp.nstr(tot, 20), '(%.1fs)' % (time.time()-t0))
        print('   parts', {k: mp.nstr(v, 8) for k, v in parts.items()}, '\n   diag', {k: mp.nstr(v, 3) for k, v in diag.items()})
