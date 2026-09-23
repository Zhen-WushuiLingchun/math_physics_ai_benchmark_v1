"""
Exact light-cone term of the first-order Duhamel formula, using the measured
unperturbed field at the bump:
  dh_lc(tau) = -(eps/2) int_{-1}^{1} W(y) psi_0(tau-(L+y-10), L+y) dy
Compare with measured dh. (eps=1e-3)
"""
import numpy as np
from td_solve import run_td, W_bump

XOBS, XMIN, XMAX, DX, DT = 10.0, -70.0, 102.0, 0.025, 0.0125
EPS = 1e-3
res = {}
for L in [20.0, 28.0]:
    T = 2*L + 40
    o = run_td(0.0, L, dx=DX, dt=DT, Tmax=T, xobs=XOBS, xmin=XMIN, xmax=XMAX,
               history_region=(L-1.0, L+1.0))
    p = run_td(EPS, L, dx=DX, dt=DT, Tmax=T, xobs=XOBS, xmin=XMIN, xmax=XMAX)
    tau = o['tau']; dh = p['h'] - o['h']
    hist = o['hist']; xb = o['hist_x']
    wb = W_bump(xb - L)
    nt = len(tau)
    dh_lc = np.zeros(nt)
    i0 = int(round((L-10.0)/DT))
    for n in range(nt):
        # integrand: Wb(y) psi0(tau - (L+y-10), L+y)
        iarg = n - (xb - XOBS)/DT          # tau - (x - xobs)  in index units
        iarg = np.clip(iarg, 0, nt-1)
        ii = np.floor(iarg).astype(int); ff = iarg - ii
        val = hist[ii, np.arange(len(xb))]*(1-ff) + hist[np.minimum(ii+1, nt-1), np.arange(len(xb))]*ff
        dh_lc[n] = -0.5*EPS*np.trapezoid(wb*val, xb)
    m = (tau > 2*L-15) & (tau < 2*L+5)
    err = np.max(np.abs(dh[m]-dh_lc[m]))/np.max(np.abs(dh[m]))
    nrm = np.sqrt(np.trapezoid(dh_lc[m]**2, tau[m])/np.trapezoid(dh[m]**2, tau[m]))
    print(f"L={L}: light-cone first-order: max rel err={err:.3f}, norm ratio={nrm:.3f}")
    res[f'L{L}'] = dict(tau=tau, dh=dh, dh_lc=dh_lc)
    np.savez(f'duhamel2_L{L:.0f}.npz', tau=tau, dh=dh, dh_lc=dh_lc)
print("saved duhamel2_L*.npz")
