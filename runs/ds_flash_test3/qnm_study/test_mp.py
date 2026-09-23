import time
import mpmath as mp
from qnm_mp import *

t0 = time.time()
Vs = V_grid()
print(f"V grid built: {len(Vs)} points, {time.time()-t0:.2f}s")

wlit = mp.mpc('0.3736716844', '-0.0889623158')
t0 = time.time()
W = jost_mp(wlit, Vs)[0]
print(f"W(lit) = {mp.nstr(W, 20)}  |W|={mp.nstr(abs(W),6)}   ({time.time()-t0:.1f}s)")

t0 = time.time()
w0, it = find_qnm_mp(wlit, Vs)
print(f"omega0 = {mp.nstr(w0, 25)}   ({it} iters, {time.time()-t0:.1f}s)")
print(f"lit    = {mp.nstr(wlit, 25)}")
print(f"diff   = {mp.nstr(w0-wlit, 8)}")
print(f"1/(2|Im w0|) = {mp.nstr(1/(2*abs(mp.im(w0))), 12)}")

# Wronskian derivative at root
h = mp.mpf('1e-8')
Wp = (jost_mp(w0+h, Vs)[0] - jost_mp(w0-h, Vs)[0])/(2*h)
print(f"W'(omega0) = {mp.nstr(Wp, 20)}")
# also at a larger h for stability check
h2 = mp.mpf('1e-6')
Wp2 = (jost_mp(w0+h2, Vs)[0] - jost_mp(w0-h2, Vs)[0])/(2*h2)
print(f"W'(omega0) h=1e-6: {mp.nstr(Wp2, 20)}")
