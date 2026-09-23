import mpmath as mp
import qnm_mp as Q
from qnm_mp import *
import time

wlit = mp.mpc('0.3736716844', '-0.0889623158')
res = {}
for dx in ['0.025', '0.0125', '0.00625']:
    Q.DX = mp.mpf(dx)
    Q._N = int((Q.XMAX - Q.XMIN) / Q.DX)
    Vs = V_arrays()
    t0 = time.time()
    w0, it = find_qnm_mp(wlit, Vs, h=Q.DX*mp.mpf('0.01'))
    res[dx] = w0
    print(f"dx={dx}: omega0 = {mp.nstr(w0, 22)}  ({it} iters, {time.time()-t0:.0f}s)")

# Richardson (RK4: error ~ dx^4)
d1 = res['0.025'] - res['0.0125']
d2 = res['0.0125'] - res['0.00625']
print("diffs:", mp.nstr(d1, 6), " / ", mp.nstr(d2, 6), " ratio Re:", mp.nstr(d1.real/d2.real, 6))
w_ext = res['0.00625'] + d2/15
print("Richardson: omega0 =", mp.nstr(w_ext, 25))
print("literature:         ", mp.nstr(wlit, 25))
print("literature - rich.  =", mp.nstr(wlit - w_ext, 8))
