"""Analyse collapse of 1-F_rec against candidate critical variables Xi_alpha = (NR-NB) - alpha log_phi k."""
import csv, math
import numpy as np

phi = (1+math.sqrt(5))/2; lphi = math.log(phi)
rows = list(csv.DictReader(open("work/results_recovery.csv")))
for r in rows:
    for kk in r:
        try: r[kk] = float(r[kk])
        except: pass

print("data:")
print(f"{'NR':>3}{'NB':>3}{'k':>3}{'delta':>7}{'xi':>7}{'1-F':>8}{'phi^-xi':>9}{'1-F/phi^-xi':>12}")
for r in sorted(rows, key=lambda r: r['xi']):
    x = 1-r['Fsdp']; p = phi**(-r['xi'])
    print(f"{int(r['NR']):>3}{int(r['NB']):>3}{int(r['k']):>3}{r['NR']-r['NB']:>7.0f}{r['xi']:>7.2f}{x:>8.4f}{p:>9.4f}{x/p:>12.4f}")

print()
print("collapse test: 1-F ~ c * phi^{-Xi_alpha}; fit log(1-F) = log c - (ln phi) * Xi_alpha")
mask = np.array([(1-r['Fsdp']) > 0.015 and (1-r['Fsdp']) < 0.75 for r in rows])
for alpha in [0.0, 0.5, 1.0, 1.5, 2.0]:
    xi = np.array([r['NR']-r['NB'] - alpha*math.log(r['k'])/lphi for r in rows])
    y = np.log(1-np.array([r['Fsdp'] for r in rows]))
    X = np.vstack([np.ones(mask.sum()), xi[mask]]).T
    coef, res, *_ = np.linalg.lstsq(X, y[mask], rcond=None)
    pred = X@coef
    rms = float(np.sqrt(np.mean((y[mask]-pred)**2)))
    # also unconstrained quadratic fit to allow curvature
    X2 = np.vstack([np.ones(mask.sum()), xi[mask], xi[mask]**2]).T
    coef2, *_ = np.linalg.lstsq(X2, y[mask], rcond=None)
    rms2 = float(np.sqrt(np.mean((y[mask]-X2@coef2)**2)))
    print(f"alpha={alpha:.1f}: slope={coef[1]:+.4f} (ideal {-lphi:+.4f}), c={math.exp(coef[0]):.4f}, rms(log)={rms:.4f}, rms(quad)={rms2:.4f}, N={int(mask.sum())}")
