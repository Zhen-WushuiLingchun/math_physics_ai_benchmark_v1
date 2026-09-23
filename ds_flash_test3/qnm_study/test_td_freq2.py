import numpy as np, time
from scipy.optimize import least_squares
from td_solve import *

def fit2(tau, h, t1, t2, w0=(0.37, -0.09)):
    m = (tau >= t1) & (tau <= t2)
    t, y = tau[m], h[m]
    def model(p):
        wr, wi, a, ph, c = p
        return a*np.exp(wi*t)*np.cos(wr*t + ph) + c
    def resid(p):
        return model(p) - y
    p0 = [w0[0], w0[1], np.max(np.abs(y)), 0.0, 0.0]
    r = least_squares(resid, p0, method='lm', max_nfev=20000)
    return r.x, np.sqrt(np.mean(r.fun**2))/np.max(np.abs(y))

out = run_td(0.0, 12.0, dx=0.025, dt=0.0125, Tmax=150.0)
tau, h = out['tau'], out['h']
print(f"C={out['C']:.8f}")
for wnd in [(15,40),(20,60),(25,70),(30,80),(20,100)]:
    p, rel = fit2(tau, h, *wnd)
    print(f"window {wnd}: Re w={p[0]:.6f}  Im w={p[1]:+.6f}  A={p[2]:.4e} phi={p[3]:.3f} c={p[4]:.1e}  relRMS={rel:.2e}")

print()
print("literature l=2:  0.3736717 - 0.0889623 i")
print("mp shooting   :  0.3892840 - 0.0441101 i")
