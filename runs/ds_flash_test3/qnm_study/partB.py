"""
Part B quantities:
 - omega_0 (validated), W'(omega_0)
 - overlap I(L) = int W(y) psi_L(L+y) psi_R(L+y) dy
 - first-order pole shift  domega1 = -eps I / W'
 - exact perturbed pole (smooth bump) by shooting  [small eps]
 - exact delta-model poles  W(w) + eps_d psi_L(L) psi_R(L) = 0  [all eps]
Saves a JSON summary.
"""
import json
import mpmath as mp
import qnm_mp as Q
from qnm_mp import *

mp.mp.dps = 40
load_asymp(2, 8)
Vs0 = V_arrays()          # unperturbed, on current DX grid

w0 = mp.mpc('0.373672834073936005', '-0.0889633453601956463')  # our converged value
# refine with a couple of Newton steps at current grid
w0, _ = find_qnm_mp(w0, Vs0)

h = mp.mpf('1e-4')
Wp = (rk4_full(w0+h, Vs0, record=False)[0] - rk4_full(w0-h, Vs0, record=False)[0])/(2*h)
print("omega0 =", mp.nstr(w0, 20))
print("W'(omega0) =", mp.nstr(Wp, 16), "  |W'| =", mp.nstr(abs(Wp), 10))

def overlap(L, record):
    idx = idx_of(L)
    i1, i2 = idx_of(L-1), idx_of(L+1)
    s = mp.mpf(0)
    for j in range(i1, i2+1):
        x = XMIN + j*DX
        y = x - L
        wgt = 1 if (j == i1 or j == i2) else 2
        s += wgt*Wb_mp(y)*record[0][j][0]*record[1][j][0]
    return s*DX/2

Ls = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 56]
_, _, _, recL, recR = rk4_full(w0, Vs0, record=True)
res = {}
print()
print("  L      I(L)                     |I(L)|        |I| e^{-2|Im w|L}")
for L in Ls:
    I = overlap(mp.mpf(L), (recL, recR))
    amp = abs(I)*mp.e**(2*abs(mp.im(w0))*L)
    res[L] = {'I_re': mp.nstr(I.real, 20), 'I_im': mp.nstr(I.imag, 20),
              'amp': mp.nstr(amp, 12)}
    print(f"{L:4d}  {mp.nstr(I, 12):26s} {mp.nstr(abs(I), 10):12s} {mp.nstr(amp,10)}")

with open('partB_overlap.json', 'w') as f:
    json.dump({'w0': [mp.nstr(w0.real, 25), mp.nstr(w0.imag, 25)],
               'Wp': [mp.nstr(Wp.real, 20), mp.nstr(Wp.imag, 20)],
               'overlap': res, 'DX': str(DX)}, f, indent=1)
print("\nsaved partB_overlap.json")
