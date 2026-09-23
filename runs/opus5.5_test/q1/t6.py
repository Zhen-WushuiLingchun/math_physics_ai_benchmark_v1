import time, sys
from exactpt import *
from ym_formulas import *
mp.mp.dps = 90
x = Fr(2, 7)
sig = orderings(4)[0]
for seed in range(79, 90):
    rng = random.Random(seed)
    pt = small_point(4, rng, R=6)
    t0 = time.time()
    try:
        v, parts, diag = ym_col_one_loop(pt, sig, (-1, 1, 1, 1), 1, x, rng)
    except (ZeroDivisionError, ValueError) as e:
        print('seed', seed, 'degenerate:', type(e).__name__, flush=True); continue
    r = ratrec(v, 10**40)
    print(seed, sig, r, mp.nstr(abs(v - mp.mpf(r.numerator)/r.denominator)/abs(v), 3), '%.0fs' % (time.time()-t0), flush=True)
    break
