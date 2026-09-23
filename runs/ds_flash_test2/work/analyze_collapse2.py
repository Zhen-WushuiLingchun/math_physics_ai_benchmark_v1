"""Refined collapse tests: candidate vs alternatives; exact-dimension variable; F vs ||C||_1^2 relation."""
import csv, math
import numpy as np

phi = (1+math.sqrt(5))/2; lphi = math.log(phi)
def fib(n):
    a, b = 0, 1
    for _ in range(n): a, b = b, a+b
    return a
rows = list(csv.DictReader(open("work/results_recovery.csv")))
for r in rows:
    for kk in r:
        try: r[kk] = float(r[kk])
        except: pass

def fit(name, var, lo=0.015, hi=0.75):
    mask = np.array([lo < (1-r['Fsdp']) < hi for r in rows])
    xi = np.array([var(r) for r in rows])
    y = np.log(1-np.array([r['Fsdp'] for r in rows]))
    X = np.vstack([np.ones(mask.sum()), xi[mask]]).T
    coef, *_ = np.linalg.lstsq(X, y[mask], rcond=None)
    rms = float(np.sqrt(np.mean((y[mask]-X@coef)**2)))
    print(f"{name:38s}: slope={coef[1]:+.4f} (ideal {-lphi:+.4f})  c={math.exp(coef[0]):.4f}  rms={rms:.4f}  N={int(mask.sum())}")
    return rms

print("=== collapse: log(1-F) vs candidate variables (one-parameter family) ===")
for alpha in [0.0, 0.5, 1.0, 1.25, 1.5, 2.0]:
    fit(f"Xi_alpha, alpha={alpha}", lambda r, a=alpha: r['NR']-r['NB']-a*math.log(r['k'])/lphi)
print()
print("=== alternatives without N_B dependence ===")
fit("NR - log_phi k", lambda r: r['NR']-math.log(r['k'])/lphi)
fit("NR - 2log_phi k", lambda r: r['NR']-2*math.log(r['k'])/lphi)
fit("NR - (3/2)log_phi k", lambda r: r['NR']-1.5*math.log(r['k'])/lphi)
print()
print("=== finite-size exact-dimension variable:  log_phi( F_{NR+1} / (k F_{NB+1}) ) ===")
fit("log_phi(DR/(k DB))", lambda r: math.log(fib(int(r['NR'])+1)/(r['k']*fib(int(r['NB'])+1)))/lphi)

print()
print("=== empirical relation between recovery infidelity and ||C_QB||_1 ===")
x2 = np.array([r['C1'] for r in rows]); y = 1-np.array([r['Fsdp'] for r in rows])
mask = y > 0.01
A = np.vstack([x2[mask]**2]).T
coef, *_ = np.linalg.lstsq(A, y[mask], rcond=None)
pred = A@coef
print(f"1-F vs ||C||_1^2: slope={coef[0]:.4f}, rms={np.sqrt(np.mean((y[mask]-pred)**2)):.4f}, N={mask.sum()}")
A2 = np.vstack([x2[mask]**2, x2[mask]**3]).T
coef2, *_ = np.linalg.lstsq(A2, y[mask], rcond=None)
pred2 = A2@coef2
print(f"1-F vs (||C||_1^2, ||C||_1^3): {coef2}, rms={np.sqrt(np.mean((y[mask]-pred2)**2)):.4f}")

print()
print("=== check ||C||_1^2 vs k*DB*TrC2 ===")
ratios = []
for r in rows:
    DB = fib(int(r['NB'])+1)
    pred = r['k']*DB*r['TrC2']
    ratios.append(r['C1']**2/pred)
print(f"mean={np.mean(ratios):.3f}  std={np.std(ratios):.3f}  min={np.min(ratios):.3f} max={np.max(ratios):.3f}")

print()
print("=== asymptotic constant c: 1-F * phi^Xi, restricted to N>=11 ===")
cs = []
for r in rows:
    N = r['NR']+r['NB']
    if N >= 11 and 1-r['Fsdp'] > 0.01:
        c = (1-r['Fsdp'])*phi**r['xi']
        cs.append((N, int(r['NR']), int(r['NB']), int(r['k']), round(r['xi'],2), round(c,3)))
for row in cs: print("   N=%d (NR,NB,k)=(%d,%d,%d) xi=%.2f :  c=%.3f" % row)
vals = [c for *_, c in cs]
print(f"mean c = {np.mean(vals):.3f} +- {np.std(vals):.3f}")
