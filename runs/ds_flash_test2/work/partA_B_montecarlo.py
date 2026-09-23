"""
Monte-Carlo checks of Part A (entropy calibration) and Part B (E Tr C_QB^2) formulas.

Part B exact claim:
  E Tr C_QB^2 = (k^2-1)(D*Bsum - Asum) / (k^2 * D * (D^2-1))
  Asum = sum_a m_a^2 n_a, Bsum = sum_a m_a n_a^2, D = sum_a m_a n_a.
Part A exact claims (k=1):
  E Tr (rho_R^alg)^2 = sum_a m_a n_a (m_a+n_a) / (D(D+1))
  E S_alg = psi(D+1) - sum_a (a_a/D) psi(a_a+1)
            + sum_a (a_a/D) [ sum_{r=max+1}^{m n} 1/r - (min-1)/(2 max) ]   (a_a=m_a n_a)
  E S_qtr = E S_alg + sum_a (m_a n_a/D) log d_a
"""
import numpy as np
from scipy.special import digamma
import math

rng = np.random.default_rng(20260923)

def rand_isometry(D, k, rng):
    A = (rng.normal(size=(D, k)) + 1j*rng.normal(size=(D, k)))/np.sqrt(2)
    Q, R = np.linalg.qr(A)
    d = np.diagonal(R).copy()
    d /= np.abs(d)
    Q = Q * d.conj()          # phase-fix; Q is Haar(Stiefel)
    return Q

def blocks_from_V(V, Aspec, Bspec):
    """V: D x k isometry, D = sum_a m_a n_a; return dict of blocks X[(a)] shape (k, m_a, n_a)."""
    X = {}
    off = 0
    for i, (m, n) in enumerate(zip(Aspec, Bspec)):
        X[i] = V[off:off+m*n, :].reshape(m, n, -1).transpose(2, 0, 1)  # (k,m,n)
        off += m*n
    assert off == V.shape[0]
    return X

def C_trace2(V, Aspec, Bspec):
    """Tr C_QB^2 for given random isometry."""
    X = blocks_from_V(V, Aspec, Bspec)
    k = V.shape[1]
    # O_{jl} blocks: (n_a x n_a) matrices  (O_{jl})_{nu nu'} = sum_mu X[j,mu,nu] conj(X[l,mu,nu'])
    # Tr rho_QB^2 = (1/k^2) sum_{jl} sum_a || X_a^{(j)dag} X_a^{(l)} ||_F^2
    tr_qb2 = 0.0
    for a in X:
        Xa = X[a]  # (k,m,n)
        # Gram[j,l,nu,nu'] = sum_mu conj(Xa[j,mu,nu]) Xa[l,mu,nu']
        G = np.einsum('jmn,lmq->jlnq', Xa.conj(), Xa)  # G_{jl} = X_a^{(j)dag} X_a^{(l)}
        # ||G_{jl}||_F^2 summed over j,l
        tr_qb2 += np.sum(np.abs(G)**2)
    tr_qb2 /= k**2
    # Tr rho_B^2 = (1/k^2) sum_{jl} sum_a || X_a^{(j)dag}X_a^{(j)}||??  -- compute directly
    # rho_B = 1/k sum_j O_jj ; Tr rho_B^2 = 1/k^2 sum_{jl} sum_a tr(O_jj O_ll)
    tr_b2 = 0.0
    for a in X:
        Xa = X[a]
        # P_j = Xa[j]^dag Xa[j]  (n x n)  (unnormalized); tr(P_j P_l)
        P = np.einsum('jmn,jmq->jnq', Xa.conj(), Xa)  # (k,n,n) : (X^dag X)_j
        # tr(P_j P_l) = sum_{nq} P[j,n,q] P[l,q,n]
        tr_b2 += np.einsum('jnq,lqn->jl', P, P).sum()
    tr_b2 /= k**2
    return tr_qb2 - tr_b2/k, tr_qb2, tr_b2

def ETrC2_formula(Aspec, Bspec, k):
    D = sum(m*n for m, n in zip(Aspec, Bspec))
    Asum = sum(m*m*n for m, n in zip(Aspec, Bspec))
    Bsum = sum(m*n*n for m, n in zip(Aspec, Bspec))
    return (k**2-1)*(D*Bsum - Asum)/(k**2*D*(D**2-1)), Asum, Bsum, D

print("=== Part B: E Tr C_QB^2 checks ===")
configs = [
    # (name, Aspec, Bspec, k)
    ("single 2x2",      [2], [2], 2),
    ("single 3x2",      [3], [2], 2),
    ("single 2x3",      [2], [3], 3),
    ("single 4x4",      [4], [4], 3),
    ("fib (2,2)",       [1,1], [1,1], 2),
    ("fib (3,3)",       [1,2], [1,2], 2),
    ("fib (3,3) k=4",   [1,2], [1,2], 4),
    ("fib (4,3)",       [2,3], [1,2], 2),
    ("fib (4,3) k=5",   [2,3], [1,2], 5),
    ("fib (5,4)",       [3,5], [2,3], 3),
]
for name, A_, B_, k in configs:
    D = sum(m*n for m, n in zip(A_, B_))
    assert k <= D
    Ns = 4000
    vals = np.zeros(Ns)
    for s in range(Ns):
        V = rand_isometry(D, k, rng)
        vals[s], _, _ = C_trace2(V, A_, B_)
    mc = vals.mean(); se = vals.std(ddof=1)/np.sqrt(Ns)
    exact, Asum, Bsum, D_ = ETrC2_formula(A_, B_, k)
    print(f"{name:18s} k={k} D={D_}: MC={mc:.6f} +- {se:.6f}   exact={exact:.6f}   diff={(mc-exact)/max(se,1e-12):+.1f} se")

print()
print("=== Part A: k=1 entropy formulas (Monte Carlo) ===")
def page_S(m, n):
    """E S of rho_m for Haar random state on C^m x C^n."""
    if m > n: m, n = n, m
    return sum(1.0/r for r in range(n+1, m*n+1)) - (m-1)/(2.0*n)

def E_S_alg_formula(Aspec, Bspec):
    D = sum(m*n for m, n in zip(Aspec, Bspec))
    val = digamma(D+1)
    for m, n in zip(Aspec, Bspec):
        a = m*n
        val -= (a/D)*digamma(a+1)
        val += (a/D)*page_S(m, n)
    return val

def E_S_qtr_formula(Aspec, Bspec, ds):
    return E_S_alg_formula(Aspec, Bspec) + sum((m*n/D)*math.log(d) for (m, n), d in zip(zip(Aspec, Bspec), ds))

def sample_state_blocks(D, Aspec, Bspec, rng):
    """Haar-random pure state in C^D written as blocks; returns p_a, rho_{R,a} eigvals."""
    psi = rng.normal(size=D) + 1j*rng.normal(size=D)
    psi /= np.linalg.norm(psi)
    out = []
    off = 0
    for m, n in zip(Aspec, Bspec):
        X = psi[off:off+m*n].reshape(m, n)   # m x n
        off += m*n
        p = np.sum(np.abs(X)**2)
        rhoR = X @ X.conj().T / p
        out.append((p, np.linalg.eigvalsh(rhoR)))
    return out

def ent(eigs):
    eigs = eigs[eigs > 1e-300]
    return -np.sum(eigs*np.log(eigs))

for name, A_, B_, dspec in [
    ("fib (2,2)", [1,1], [1,1], [1, (1+math.sqrt(5))/2]),
    ("fib (3,3)", [1,2], [1,2], [1, (1+math.sqrt(5))/2]),
    ("fib (4,3)", [2,3], [1,2], [1, (1+math.sqrt(5))/2]),
    ("single 4x4", [4], [4], [1]),
]:
    D = sum(m*n for m, n in zip(A_, B_))
    Ns = 4000
    Salg = np.zeros(Ns); Sqtr = np.zeros(Ns); purity = np.zeros(Ns)
    for s in range(Ns):
        bl = sample_state_blocks(D, A_, B_, rng)
        # S_alg
        Salg[s] = sum(p*ent(e) for p, e in bl) - sum(p*math.log(p) for p, e in bl)
        Sqtr[s] = Salg[s] + sum(p*math.log(d) for (p, e), d in zip(bl, dspec))
        purity[s] = sum((p**2)*np.sum(e**2) for p, e in bl)
    mc1 = Salg.mean(); se1 = Salg.std(ddof=1)/np.sqrt(Ns)
    ex1 = E_S_alg_formula(A_, B_)
    mc2 = Sqtr.mean(); se2 = Sqtr.std(ddof=1)/np.sqrt(Ns)
    ex2 = E_S_qtr_formula(A_, B_, dspec)
    D_ = D
    mns = sum(m*n*(m+n) for m, n in zip(A_, B_))/(D_*(D_+1))
    print(f"{name:12s} D={D}: E S_alg MC={mc1:.5f}+-{se1:.5f} exact={ex1:.5f} | E S_qtr MC={mc2:.5f}+-{se2:.5f} exact={ex2:.5f}")
    print(f"{'':12s}       E Tr(rho_R^alg)^2 MC={purity.mean():.6f}+-{purity.std(ddof=1)/np.sqrt(Ns):.6f} exact={mns:.6f}")
