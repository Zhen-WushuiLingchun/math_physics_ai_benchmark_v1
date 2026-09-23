"""
Verify  domega = -eps I(L)/W'(omega0)  against exact shooting with the smooth bump.
Also delta-model: exact pole equation  W(w) + eps_d psi_L(L) psi_R(L) = 0.
"""
import mpmath as mp
import json
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

print("  eps     L     pred domega                 exact domega                ratio")
rows = []
for eps_s, L in [('1e-3', 16), ('1e-3', 24), ('3e-3', 16), ('3e-3', 24), ('1e-2', 16), ('3e-2', 16)]:
    eps = mp.mpf(eps_s)
    I = I_of_L(L)
    pred = -eps*I/Wp
    # exact: perturbed potential
    Vs = V_arrays(eps=eps, L=mp.mpf(L))
    w_exact, it = find_qnm_mp(w0 + pred, Vs, tol=mp.mpf('1e-25'))
    d_exact = w_exact - w0
    ratio = d_exact/pred
    print(f" {eps_s:6s} {L:4d}  {mp.nstr(pred,10):26s} {mp.nstr(d_exact,10):26s} {mp.nstr(ratio,8)}")
    rows.append({'eps': eps_s, 'L': L, 'pred_re': mp.nstr(pred.real,12), 'pred_im': mp.nstr(pred.imag,12),
                 'exact_re': mp.nstr(d_exact.real,12), 'exact_im': mp.nstr(d_exact.imag,12)})

with open('partB_shift_compare.json','w') as f:
    json.dump(rows, f, indent=1)
print("\nsaved partB_shift_compare.json")
