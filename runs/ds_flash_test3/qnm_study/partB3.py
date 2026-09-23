import mpmath as mp, json
import qnm_mp as Q
from qnm_mp import *

mp.mp.dps = 40
load_asymp(2, 8)
Vs0 = V_arrays()
w0 = mp.mpc('0.373672834073936005', '-0.0889633453601956463')
w0, _ = find_qnm_mp(w0, Vs0)
h = mp.mpf('1e-4')
Wp = (rk4_full(w0+h, Vs0, record=False)[0] - rk4_full(w0-h, Vs0, record=False)[0])/(2*h)

with open('partB_overlap.json') as f:
    dat = json.load(f)
def I_of_L(L):
    d = dat['overlap'][str(L)]
    return mp.mpc(d['I_re'], d['I_im'])

for L in [16, 24]:
    I = I_of_L(L)
    print(f"=== L={L}: pred slope -I/W' = {mp.nstr(-I/Wp, 10)}")
    for eps_s in ['1e-7','1e-6','1e-5','1e-4','3e-4']:
        eps = mp.mpf(eps_s)
        Vs = V_arrays(eps=eps, L=mp.mpf(L))
        pred = -eps*I/Wp
        w_exact, it = find_qnm_mp(w0+pred, Vs, tol=mp.mpf('1e-28'))
        d = w_exact - w0
        print(f"  eps={eps_s:6s}: domega/eps = {mp.nstr(d/eps, 12)}  (resid |W|={mp.nstr(abs(rk4_full(w_exact, Vs, record=False)[0]),4)})")
