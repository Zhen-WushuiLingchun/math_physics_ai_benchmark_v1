"""
Optimal entanglement recovery fidelity F_rec for the sector-resolved code.
F_rec = sum_a H_a^2 / k,   H_a = max_{sigma density on B_a} Tr sqrt( L_a^dag (I_k (x) sigma) L_a ),
rho_QB^(a) = L_a L_a^dag (unnormalised, trace q_a).
Certificate (concavity + Euler homogeneity):  h(s) <= H* <= h(s)/2 + lambda_max(grad h(s)).
"""
import numpy as np
from common import *

def sector_states(V, m, n, k):
    """returns list of (rho_QB^(a) factor L_a [(k n_a) x m_a], rho_QR^(a) [(k m_a) x (k m_a)])"""
    A = blocks(V, m, n)
    out = []
    for a in range(2):
        G = np.concatenate([A[a][j].T for j in range(k)], axis=0) / np.sqrt(k)   # rows (j,b), cols r
        GR = np.concatenate([A[a][j] for j in range(k)], axis=0) / np.sqrt(k)    # rows (j,r), cols b
        out.append((G, GR))
    return out

def _psd_sqrt_inv(Y, tol=1e-14):
    w, U = np.linalg.eigh(Y)
    w = np.clip(w, 0, None)
    keep = w > tol * max(w.max(), 1e-300)
    s = np.zeros_like(w); s[keep] = 1 / np.sqrt(w[keep])
    return (U * s) @ U.conj().T, np.sqrt(w).sum()

def reduce_factor(G):
    """replace G (p x q) by L (p x r), r = min(p,q), with LL^dag = GG^dag"""
    p, q = G.shape
    if q <= p:
        return G
    w, U = np.linalg.eigh(G @ G.conj().T)
    w = np.clip(w, 0, None)
    return U * np.sqrt(w)

def maxH(G, k, n, iters=5000, tol=1e-10, sigma0=None):
    """maximise h(sigma) = Tr sqrt(L^dag (I_k x sigma) L). returns (h_lower, h_upper, sigma)"""
    L = reduce_factor(G)
    if n == 1:
        h = np.linalg.svd(L, compute_uv=False).sum()
        return h, h, np.ones((1, 1))
    sig = np.eye(n) / n if sigma0 is None else sigma0
    best = (0, np.inf, sig)
    for it in range(iters):
        Lt = L.reshape(k, n, -1)                          # [j, b, r]
        Y = sum(Lt[j].conj().T @ sig @ Lt[j] for j in range(k))
        Yih, h = _psd_sqrt_inv(Y)
        # gradient = 1/2 Tr_Q[ L Y^{-1/2} L^dag ]
        g = 0.5 * sum(Lt[j] @ Yih @ Lt[j].conj().T for j in range(k))
        g = (g + g.conj().T) / 2
        up = h / 2 + np.linalg.eigvalsh(g)[-1]
        if h > best[0]: best = (h, min(best[1], up), sig)
        else: best = (best[0], min(best[1], up), best[2])
        if best[1] - best[0] < tol * best[0]:
            break
        # multiplicative fixed-point step  sigma <- (2/h) sigma^{1/2} g sigma^{1/2}
        ws, Us = np.linalg.eigh(sig)
        sq = (Us * np.sqrt(np.clip(ws, 0, None))) @ Us.conj().T
        new = sq @ g @ sq * (2 / h)
        new = (new + new.conj().T) / 2
        new /= np.trace(new).real
        sig = new
    return best[0], best[1], best[2]

from scipy.optimize import minimize

def _h_and_grad(sig, Lt):
    k = Lt.shape[0]
    Y = sum(Lt[j].conj().T @ sig @ Lt[j] for j in range(k))
    Yih, h = _psd_sqrt_inv(Y)
    g = 0.5 * sum(Lt[j] @ Yih @ Lt[j].conj().T for j in range(k))
    return h, (g + g.conj().T) / 2

def maxH_lbfgs(G, k, n, tol=1e-11, pre=30, maxiter=5000):
    """L-BFGS on sigma = X X^dag / Tr(X X^dag), warm-started by the multiplicative iteration.
    Returns certified (h_lower, h_upper, sigma)."""
    L = reduce_factor(G)
    if n == 1:
        h = np.linalg.svd(L, compute_uv=False).sum()
        return h, h, np.ones((1, 1))
    Lt = L.reshape(k, n, -1)
    hl, hu, sig = maxH(G, k, n, iters=pre, tol=tol)
    w, U = np.linalg.eigh(sig)
    X0 = (U * np.sqrt(np.clip(w, 1e-300, None))) @ U.conj().T
    def f(x):
        X = (x[:n * n] + 1j * x[n * n:]).reshape(n, n)
        t = np.trace(X @ X.conj().T).real
        s = X @ X.conj().T / t
        h, g = _h_and_grad(s, Lt)
        gr = (2 / t) * (g - (h / 2) * np.eye(n)) @ X          # d h / d conj(X)  (x2 for real/imag split)
        return -h, -np.concatenate([gr.real.ravel(), gr.imag.ravel()])
    x0 = np.concatenate([X0.real.ravel(), X0.imag.ravel()])
    res = minimize(f, x0, jac=True, method='L-BFGS-B',
                   options=dict(maxiter=maxiter, maxcor=30, ftol=1e-16, gtol=1e-13))
    X = (res.x[:n * n] + 1j * res.x[n * n:]).reshape(n, n)
    s = X @ X.conj().T; s /= np.trace(s).real
    h, g = _h_and_grad(s, Lt)
    up = h / 2 + np.linalg.eigvalsh(g)[-1]
    return max(h, hl), min(up, hu), s

def Frec(V, m, n, k, **kw):
    """returns (F_lower, F_upper, per-sector info)"""
    st = sector_states(V, m, n, k)
    lo = up = 0.0; info = []
    for a in range(2):
        G, _ = st[a]
        hl, hu, sig = maxH_lbfgs(G, k, n[a], **kw)
        lo += hl**2 / k; up += hu**2 / k
        q = np.linalg.norm(G) ** 2
        info.append(dict(q=q, F_a=hl**2 / k / q, sigma=sig))
    return lo, up, info
