"""
Check the first-order echo formula
   dh(tau) = -(eps/2) int_{-1}^{1} W(y) F(tau - 2L - 2y + 10) dy,
with F(u) = int_0^{u+10} h0(s) ds  (outgoing profile, F' = h0(.-10)).
Also compute window-local mismatch M and the echo/background ratio.
"""
import json
import numpy as np
from scipy.interpolate import interp1d
from td_solve import run_td, initial_C

XOBS, XMIN, XMAX, DX, DT = 10.0, -70.0, 102.0, 0.025, 0.0125

def run(eps, L, Tmax):
    return run_td(eps, L, dx=DX, dt=DT, Tmax=Tmax, xobs=XOBS,
                  xmin=XMIN, xmax=XMAX)

out0 = run(0.0, 16.0, 130.0)
tau0, h0 = out0['tau'], out0['h']
norm2 = np.trapezoid(h0**2, tau0)

# cumulative profile F(u) = int_0^{u+10} h0 ds
cum = np.concatenate([[0.0], np.cumsum(0.5*(h0[1:]+h0[:-1])*DT)])
def Fprof(u):
    i = (u + 10.0)/DT
    i = np.clip(i, 0, len(cum)-1)
    return np.interp(i, np.arange(len(cum)), cum)

def Wb(y):
    y = np.asarray(y)
    out = np.zeros_like(y)
    m = np.abs(y) < 1
    out[m] = np.exp(1 - 1/(1-y[m]**2))
    return out

ys = np.linspace(-1, 1, 401)
Wy = Wb(ys)

out = {}
for eps, L in [(1e-3, 20.0), (1e-3, 28.0)]:
    o = run(eps, L, 2*L+40)
    tau, h = o['tau'], o['h']
    n = min(len(tau), len(tau0))
    dh = h[:n] - h0[:n]
    t = tau[:n]
    # predicted echo
    dh_pred = np.zeros_like(t)
    for i, tt in enumerate(t):
        arg = tt - 2*L + 10.0 - 2*ys
        dh_pred[i] = -0.5*eps*np.trapezoid(Wy*Fprof(arg), ys)
    m1 = (t > 2*L-14) & (t < 2*L+16)
    err = np.max(np.abs(dh[m1]-dh_pred[m1]))/np.max(np.abs(dh[m1]))
    # window-local mismatch over the echo window
    W1, W2 = 2*L-15, 2*L-4
    mw = (t >= W1) & (t <= W2)
    a = np.trapezoid(h[:n][mw]*h0[:n][mw], t[mw])
    ne = np.sqrt(np.trapezoid(h[:n][mw]**2, t[mw]))
    n0 = np.sqrt(np.trapezoid(h0[:n][mw]**2, t[mw]))
    M_win = 1 - abs(a)/(ne*n0)
    ratio = np.sqrt(np.trapezoid(dh[mw]**2, t[mw]))/n0
    # global M at T_end
    M_end = 1 - abs(np.trapezoid(h[:n]*h0[:n], t))/(np.sqrt(np.trapezoid(h[:n]**2,t))*np.sqrt(np.trapezoid(h0[:n]**2,t)))
    print(f"eps={eps:g} L={L:g}: echo formula max-rel-err={err:.3f}  "
          f"M_window[{W1:.0f},{W2:.0f}]={M_win:.4f}  ||dh||/||h0||_window={ratio:.2f}  M_end={M_end:.3e}")
    out[f'{eps}_{L}'] = dict(tau=t, dh=dh, dh_pred=dh_pred, err=float(err),
                             M_win=float(M_win), ratio_win=float(ratio), M_end=float(M_end))

np.savez('echo_check.npz', **{f'{k}_{kk}': vv for k, v in out.items() for kk, vv in v.items()})
print("saved echo_check.npz")
