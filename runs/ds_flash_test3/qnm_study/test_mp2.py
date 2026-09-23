import time
import mpmath as mp
from qnm_mp import *

t0 = time.time()
Vs = V_arrays()
print(f"V arrays: {time.time()-t0:.2f}s")

wlit = mp.mpc('0.3736716844', '-0.0889623158')
t0 = time.time()
W = rk4_full(wlit, Vs, record=False)[0]
print(f"W(lit) = {mp.nstr(W, 15)}  |W|={mp.nstr(abs(W),4)}   ({time.time()-t0:.2f}s)")

t0 = time.time()
w0, it = find_qnm_mp(wlit, Vs)
print(f"omega0 = {mp.nstr(w0, 30)}   ({it} iters, {time.time()-t0:.1f}s)")
print(f"lit    = {mp.nstr(wlit, 30)}")
print(f"diff   = {mp.nstr(w0 - wlit, 8)}")

h = mp.mpf('0.001')
Wp = (rk4_full(w0+h, Vs, record=False)[0] - rk4_full(w0-h, Vs, record=False)[0])/(2*h)
print(f"W'(omega0) = {mp.nstr(Wp, 20)}")
print(f"|W'| = {mp.nstr(abs(Wp), 10)}")
