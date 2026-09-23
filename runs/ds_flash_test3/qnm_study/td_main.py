"""
Part A/C: time-domain computations.
 - h_0(tau), h_eps(tau) at x_o = 10
 - causality check: h_eps = h_0 for tau < 2L-13
 - E(T) and M(T)
 - echo first-order formula check
Saves td_results.json and arrays in npz.
"""
import json
import numpy as np
from td_solve import run_td, initial_C, fit_frequency

XMIN, XMAX = -70.0, 102.0
DX, DT = 0.025, 0.0125
XOBS = 10.0

def run(eps, L, dx=DX, dt=DT, Tmax=None):
    return run_td(eps, L, dx=dx, dt=dt, Tmax=Tmax, xobs=XOBS,
                  xmin=XMIN, xmax=XMAX)

# unperturbed reference (large Tmax)
TMAX0 = 130.0
out0 = run(0.0, 16.0, Tmax=TMAX0)
tau0, h0 = out0['tau'], out0['h']
C = out0['C']
norm2_h0 = np.trapezoid(h0**2, tau0)
print(f"C = {C:.10f}   ||h0||^2 = {norm2_h0:.10f}")

# rebuild run_td signature check: xmin/xmax are parameters now
def E_of_T(tau, dh, T):
    m = tau <= T
    return np.sqrt(np.trapezoid(dh[m]**2, tau[m])/norm2_h0)

def M_of_T(tau, he, T):
    m = tau <= T
    a = np.trapezoid(he[m]*h0[:len(tau)][m], tau[m])
    ne = np.sqrt(np.trapezoid(he[m]**2, tau[m]))
    n0 = np.sqrt(np.trapezoid(h0[:len(tau)][m]**2, tau[m]))
    if ne == 0 or n0 == 0:
        return np.nan
    return 1.0 - abs(a)/(ne*n0)

cases = [(1e-3, 12.0), (1e-3, 16.0), (1e-3, 20.0), (1e-3, 24.0),
         (1e-3, 28.0), (1e-3, 32.0), (1e-2, 24.0), (1e-4, 32.0)]
results = {}
store = {}
for eps, L in cases:
    Tmax = 2*L + 45.0
    out = run(eps, L, Tmax=Tmax)
    tau, h = out['tau'], out['h']
    # common time base
    n = min(len(tau), len(tau0))
    dh = h[:n] - h0[:n]
    t = tau[:n]
    # causality
    ic = int(round((2*L-13.0)/DT))
    causal = np.max(np.abs(dh[:max(ic,1)]))
    # E at various T
    E_end = E_of_T(t, dh, t[-1])
    E_echo = E_of_T(t, dh, 2*L-6.0)
    E_early = E_of_T(t, dh, max(2*L-13.0, 0.0))
    M_end = M_of_T(t, h[:n], t[-1])
    M_echo = M_of_T(t, h[:n], 2*L-6.0)
    # echo peak amplitude
    i1 = int(round((2*L-15.0)/DT)); i2 = int(round((2*L-3.0)/DT))
    echo_amp = np.max(np.abs(dh[i1:i2]))
    key = f"eps{eps:g}_L{L:g}"
    results[key] = dict(eps=eps, L=L, Tmax=float(t[-1]), causal_max=float(causal),
                        E_end=float(E_end), E_echo=float(E_echo), E_early=float(E_early),
                        M_end=float(M_end), M_echo=float(M_echo),
                        echo_amp=float(echo_amp), E_over_eps=float(E_end/eps))
    store[key] = dict(tau=t, h=h, dh=dh)
    print(f"eps={eps:g} L={L:g}: causal_max={causal:.2e}  E_end={E_end:.4e} "
          f"(E/eps={E_end/eps:.3f})  E_echo={E_echo:.3e}  M_end={M_end:.2e} M_echo={M_echo:.3f} echo_amp={echo_amp:.3e}")

np.savez('td_arrays.npz', tau0=tau0, h0=h0,
         **{f"{k}_h": v['h'] for k, v in store.items()},
         **{f"{k}_dh": v['dh'] for k, v in store.items()})
with open('td_results.json', 'w') as f:
    json.dump(results, f, indent=1)
print("\nsaved td_results.json / td_arrays.npz")
