"""
Part A: exact entropy formulas, asymptotics, replica conventions, basis invariance.
"""
import numpy as np, math
from scipy.special import digamma
from sympy import fibonacci, sqrt, simplify, nsimplify, Rational, pi, exp, log, N as SN, Symbol, oo, digamma as sdigamma

phi = (1+math.sqrt(5))/2; lphi = math.log(phi)
def fib(n):
    a, b = 0, 1
    for _ in range(n): a, b = b, a+b
    return a

def specs(NR, NB):
    return [fib(NR-1), fib(NR)], [fib(NB-1), fib(NB)]

def page_S(m, n):
    if m > n: m, n = n, m
    return sum(1.0/r for r in range(n+1, m*n+1)) - (m-1)/(2.0*n)

def E_S_alg_exact(NR, NB):
    m, n = specs(NR, NB)
    D = sum(a*b for a, b in zip(m, n))
    val = digamma(D+1)
    for ma, na in zip(m, n):
        al = ma*na
        val -= (al/D)*digamma(al+1)
        val += (al/D)*page_S(ma, na)
    return val

def E_S_qtr_exact(NR, NB):
    m, n = specs(NR, NB)
    D = sum(a*b for a, b in zip(m, n))
    return E_S_alg_exact(NR, NB) + (m[1]*n[1]/D)*lphi

# asymptotic constants
p1 = 1/(phi**2+1); ptau = phi**2/(phi**2+1)
H = -p1*math.log(p1) - ptau*math.log(ptau)
C1 = H - p1*lphi
print("p_1 = %.10f, p_tau = %.10f" % (p1, ptau))
print("H(p*) = %.10f, C1 = H - p_1 ln phi = %.10f" % (H, C1))
print("2phi/sqrt5 = %.10f" % (2*phi/math.sqrt(5)))
print()

print("asymptotic check for E S_alg:  predicted = N_B logphi - 0.5 log5 + C1 - phi^{-|delta|}/2 + 1/(2D)")
print(f"{'(NR,NB)':>10} {'delta':>6} {'exact':>10} {'asympt':>10} {'diff':>12} {'phi^-2NB':>10}")
for (NR, NB) in [(8,8),(9,8),(10,10),(11,10),(12,12),(13,12),(14,12),(15,14),(16,14),(16,12)]:
    D = fib(NR+NB-1)
    ex = E_S_alg_exact(NR, NB)
    delta = NR-NB
    asym = NB*lphi - 0.5*math.log(5) + C1 - phi**(-abs(delta))/2 + 1/(2*D)
    print(f"({NR:>3},{NB:>3}) {delta:>6} {ex:>10.6f} {asym:>10.6f} {ex-asym:>12.3e} {phi**(-2*NB):>10.3e}")

print()
print("same check for E S_qtr (adds p_tau ln phi = %.6f):" % (ptau*lphi))
for (NR, NB) in [(10,10),(12,12),(14,12),(16,14)]:
    D = fib(NR+NB-1); delta = NR-NB
    ex = E_S_qtr_exact(NR, NB)
    asym = NB*lphi - 0.5*math.log(5) + C1 - phi**(-abs(delta))/2 + 1/(2*D) + ptau*lphi
    print(f"({NR:>3},{NB:>3}) {delta:>6} {ex:>10.6f} {asym:>10.6f} {ex-asym:>12.3e}")

print()
print("=== replica conventions for a small case (fib (3,3)), MC ===")
rng = np.random.default_rng(3)
m, n = specs(3,3); D = sum(a*b for a,b in zip(m,n))
Ns = 200000
logZ1 = np.zeros(Ns); logZ2 = np.zeros(Ns); Ss = np.zeros(Ns)
for s in range(Ns):
    psi = rng.normal(size=D) + 1j*rng.normal(size=D); psi /= np.linalg.norm(psi)
    off = 0; pur1 = 0.0; S1 = 0.0
    blocks = []
    for ma, na in zip(m, n):
        X = psi[off:off+ma*na].reshape(ma, na); off += ma*na
        blocks.append(X)
    # S_alg = -Tr rho_R^alg log rho_R^alg = -sum_a Tr(X X^dag log(X X^dag))
    rhoR = np.zeros((sum(m), sum(m)), dtype=complex)
    o = 0
    for X, ma in zip(blocks, m):
        rhoR[o:o+ma, o:o+ma] = X @ X.conj().T; o += ma
    ev = np.linalg.eigvalsh(rhoR); ev = ev[ev > 1e-15]
    Ss[s] = -np.sum(ev*np.log(ev))
    # Z_n for n=2 : Tr (rho_R^alg)^2 and Tr (rho_R^alg)^1 = 1
    Z2 = np.sum(ev**2)
    logZ2[s] = math.log(Z2)
print("E[S_alg]          = %.6f   (exact: %.6f)" % (Ss.mean(), E_S_alg_exact(3,3)))
print("E[-log Z_2]       = %.6f   [quenched 2nd Renyi]" % (-logZ2).mean())
print("log E[Z_2]: E[Z_2] = %.6f -> -log E[Z_2] = %.6f  (exact ETr^2 = %.6f)"
      % (np.exp(logZ2).mean(), -math.log(np.exp(logZ2).mean()),
         sum(a*b*(a+b) for a,b in zip(m,n))/(D*(D+1))))
print("S(E rho_R) = log D = %.6f  [entropy of the averaged state]" % math.log(D))

print()
print("=== basis invariance check (random blockwise unitaries) ===")
def entropy_alg(blocks):
    rhoR = np.zeros((sum(b.shape[0] for b in blocks),)*2, dtype=complex); o = 0
    for X in blocks:
        ma = X.shape[0]; rhoR[o:o+ma, o:o+ma] = X@X.conj().T; o += ma
    ev = np.linalg.eigvalsh(rhoR); ev = ev[ev > 1e-15]
    return -np.sum(ev*np.log(ev)), rhoR
psi = rng.normal(size=D) + 1j*rng.normal(size=D); psi /= np.linalg.norm(psi)
off = 0; blocks = []
for ma, na in zip(m, n):
    blocks.append(psi[off:off+ma*na].reshape(ma, na)); off += ma*na
S0, rhoR0 = entropy_alg(blocks)
errs = []
for trial in range(20):
    bl2 = []
    for X, na in zip(blocks, n):
        Qr, _ = np.linalg.qr(rng.normal(size=(X.shape[0],)*2)+1j*rng.normal(size=(X.shape[0],)*2))
        Qb, _ = np.linalg.qr(rng.normal(size=(na, na))+1j*rng.normal(size=(na, na)))
        bl2.append(Qr @ X @ Qb.conj().T)
    S1, _ = entropy_alg(bl2)
    errs.append(abs(S0-S1))
print("S_alg invariant under X_a -> U_R X_a U_B^dag:  max|dS| = %.2e" % max(errs))
print()
print("F-matrix unitarity:", np.round(np.array([[1/phi, 1/math.sqrt(phi)],[1/math.sqrt(phi), -1/phi]]) @ np.array([[1/phi, 1/math.sqrt(phi)],[1/math.sqrt(phi), -1/phi]]), 12))
