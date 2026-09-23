"""
Double-limit table: L = c ln(1/eps), eps = 1e-3 -> L = 20.72, 38.83, 55.26 for c=3,5.62,8.
For each: first-order shift, exact smooth-bump pole shift (shooting), and E/eps (TD).
"""
import json
import numpy as np
import mpmath as mp
import qnm_mp as Q
from qnm_mp import *
from td_solve import run_td

mp.mp.dps = 40
load_asymp(2, 8)
Vs0 = V_arrays()
w0 = mp.mpc('0.373672834073936005', '-0.0889633453601956463')
w0, _ = find_qnm_mp(w0, Vs0)
h = mp.mpf('1e-4')
Wp = (rk4_full(w0+h, Vs0, record=False)[0] - rk4_full(w0-h, Vs0, record=False)[0])/(2*h)

with open('partB_overlap.json') as f:
    dat = json.load(f)
def I_of(L):
    Ls = sorted(int(k) for k in dat['overlap'])
    if L in Ls:
        d = dat['overlap'][str(L)]
        return mp.mpc(d['I_re'], d['I_im'])
    return None

eps = mp.mpf('1e-3')
out = []
print(" c      L        mu=eps e^{2|Im w|L}   pred |dw|/|w0|    exact |dw|/|w0|   E/eps (TD)")
for c, L in [(3.0, 20.717), (5.62, 38.830), (8.0, 55.265)]:
    # overlap: need I(L) at this L; compute directly
    _, _, _, recL, recR = rk4_full(w0, Vs0, record=True)
    Lm = mp.mpf(str(L))
    # integrate over bump support using recorded profiles
    i1, i2 = idx_of(Lm-1), idx_of(Lm+1)
    s = mp.mpf(0)
    for j in range(i1, i2+1):
        x = XMIN + j*DX
        wgt = 1 if (j == i1 or j == i2) else 2
        s += wgt*Wb_mp(x-Lm)*recL[j][0]*recR[j][0]
    I = s*DX/2
    pred = eps*I/Wp
    mu = eps*mp.e**(2*abs(mp.im(w0))*Lm)
    # exact pole
    Vs = V_arrays(eps=eps, L=Lm)
    w_exact, it = find_qnm_mp(w0+pred, Vs, tol=mp.mpf('1e-28'))
    d_exact = w_exact - w0
    # TD
    o = run_td(float(eps), float(L), dx=0.025, dt=0.0125, Tmax=2*float(L)+40,
               xobs=10.0, xmin=-70.0, xmax=102.0)
    o0 = run_td(0.0, float(L), dx=0.025, dt=0.0125, Tmax=2*float(L)+40,
                xobs=10.0, xmin=-70.0, xmax=102.0)
    tau, hh = o['tau'], o['h']
    h0 = o0['h']
    n = min(len(hh), len(h0))
    dh = hh[:n]-h0[:n]
    norm2 = np.trapezoid(h0**2, o0['tau'])
    E = np.sqrt(np.trapezoid(dh**2, tau[:n])/norm2)
    row = dict(c=c, L=L, mu=mp.nstr(mu, 8),
               pred_rel=mp.nstr(abs(pred)/abs(w0), 6),
               exact_rel=mp.nstr(abs(d_exact)/abs(w0), 6),
               pred=[mp.nstr(pred.real, 10), mp.nstr(pred.imag, 10)],
               exact=[mp.nstr(d_exact.real, 10), mp.nstr(d_exact.imag, 10)],
               E_over_eps=float(E/float(eps)))
    out.append(row)
    print(f" {c:5.2f}  {L:7.3f}  {mp.nstr(mu,6):12s}    {row['pred_rel']:12s}     {row['exact_rel']:12s}     {row['E_over_eps']:.4f}")

with open('double_limit.json', 'w') as f:
    json.dump(out, f, indent=1)
print("saved double_limit.json")
