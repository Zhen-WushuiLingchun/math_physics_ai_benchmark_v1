"""Verify: (a) Petz fidelity closed formula, (b) SDP dual upper bound, (c) dephasing vs full trace."""
import numpy as np, math
import cvxpy as cp

rng = np.random.default_rng(5)

def fib(n):
    a, b = 0, 1
    for _ in range(n): a, b = b, a+b
    return a

def specs(NR, NB):
    return [fib(NR-1), fib(NR)], [fib(NB-1), fib(NB)]

def rand_isometry(D, k, rng):
    A = (rng.normal(size=(D, k)) + 1j*rng.normal(size=(D, k)))/np.sqrt(2)
    Q, R = np.linalg.qr(A); d = np.diagonal(R).copy(); d /= np.abs(d)
    return Q * d.conj()

def blocks(V, m, n):
    X = []; off = 0
    for ma, na in zip(m, n):
        X.append(V[off:off+ma*na, :].reshape(ma, na, -1).transpose(2, 0, 1)); off += ma*na
    return X

def rho_from_blocks(X, k, which='R'):
    dims = [x.shape[1] if which == 'R' else x.shape[2] for x in X]
    d = sum(dims)
    rho = np.zeros((k*d, k*d), dtype=complex)
    off = 0
    for x, dx in zip(X, dims):
        for q in range(k):
            for qp in range(k):
                if which == 'R': blk = x[q] @ x[qp].conj().T
                else: blk = x[q].T @ x[qp].conj()
                rho[q*d+off:q*d+off+dx, qp*d+off:qp*d+off+dx] = blk/k
        off += dx
    return rho, d

def W_from_rho(rho, dR, k):
    W = np.zeros((dR*k, dR*k), dtype=complex)
    for l in range(k):
        for a in range(dR):
            for lp in range(k):
                for ap in range(dR):
                    W[a*k+l, ap*k+lp] = rho[lp*dR+ap, l*dR+a]
    return W

def freq_sdp(rho, dR, k):
    W = W_from_rho(rho, dR, k)
    J = cp.Variable((dR*k, dR*k), hermitian=True)
    Al = [np.kron(np.eye(dR), np.eye(k)[l]) for l in range(k)]
    psum = sum(Al[l] @ J @ Al[l].conj().T for l in range(k))
    prob = cp.Problem(cp.Maximize(cp.real(cp.trace(W @ J))/k), [J >> 0, psum == np.eye(dR)])
    prob.solve(solver=cp.SCS, eps=1e-9, max_iters=60000)
    # dual upper bound: min Tr Y s.t. Y x I_L >= W/k
    Y = cp.Variable((dR, dR), hermitian=True)
    prob2 = cp.Problem(cp.Minimize(cp.real(cp.trace(Y))),
                       [cp.kron(Y, np.eye(k)) - W/k >> 0])
    prob2.solve(solver=cp.SCS, eps=1e-9, max_iters=60000)
    return float(prob.value), float(prob2.value)

def petz_fidelity_formula(X, k):
    tot = 0.0
    for x in X:  # x: (k,m,n)
        S = sum(x[q] @ x[q].conj().T for q in range(x.shape[0]))
        w, U = np.linalg.eigh(S)
        Sinv = (U*(1/np.sqrt(np.maximum(w, 1e-15)))) @ U.conj().T
        G = np.zeros((x.shape[2], x.shape[2]), dtype=complex)
        for j in range(k):
            G += x[j].conj().T @ Sinv @ x[j]
        tot += np.real(np.trace(G @ G))
    return tot/k**2

def petz_fidelity_direct(X, rhoQR, k):
    m = [x.shape[1] for x in X]; dR = sum(m)
    Slist = []
    for x, ma in zip(X, m):
        S = sum(x[q] @ x[q].conj().T for q in range(k))
        w, U = np.linalg.eigh(S)
        Slist.append((U*(1/np.sqrt(np.maximum(w, 1e-15)))) @ U.conj().T)
    def Dfun(Y):
        out = np.zeros((k, k), dtype=complex); off = 0
        for x, ma, Sinv in zip(X, m, Slist):
            Ya = Y[off:off+ma, off:off+ma]; Z = Sinv @ Ya @ Sinv
            for nu in range(x.shape[2]):
                K = x[:, :, nu].T
                out += K.conj().T @ Z @ K
            off += ma
        return out
    T = rhoQR.reshape(k, dR, k, dR); F = 0.0
    for r in range(dR):
        for rp in range(dR):
            Xm = np.zeros((dR, dR), dtype=complex); Xm[r, rp] = 1
            Dr = Dfun(Xm)
            for q in range(k):
                for qp in range(k):
                    F += T[q, r, qp, rp]*Dr[q, qp]
    return float(np.real(F))/k

def rhoQR_full(X, k):
    """Un-dephased R-channel output (includes cross-sector coherences)."""
    dR = sum(x.shape[1] for x in X)
    offs = [0]
    for x in X: offs.append(offs[-1]+x.shape[1])
    rho = np.zeros((k*dR, k*dR), dtype=complex)
    for q in range(k):
        for qp in range(k):
            for a, xa in enumerate(X):
                for ap, xap in enumerate(X):
                    blk = xa[q] @ xap[qp].conj().T   # (m_a, m_ap)
                    rho[q*dR+offs[a]:q*dR+offs[a]+xa.shape[1],
                        qp*dR+offs[ap]:qp*dR+offs[ap]+xap.shape[1]] = blk/k
    return rho, dR

def decoupling_distance(X, k):
    rhoQB, dB = rho_from_blocks(X, k, 'B')
    rhoB = rhoQB.reshape(k, dB, k, dB).trace(axis1=0, axis2=2)
    C = rhoQB - np.kron(np.eye(k)/k, rhoB)
    return 0.5*float(np.sum(np.abs(np.linalg.eigvalsh(C)))), float(np.real(np.trace(C@C)))

print("config          F_sdp      F_dual     F_petz(form)  F_petz(direct)  dec1   TrC2      F_sdp_full")
for (NR, NB, k) in [(3,3,2),(4,3,2),(5,3,2),(6,3,2),(5,4,3),(3,2,2),(4,4,2)]:
    m, n = specs(NR, NB); D = fib(NR+NB-1)
    V = rand_isometry(D, k, rng); X = blocks(V, m, n)
    rhoQR, dR = rho_from_blocks(X, k, 'R')
    fs, fu = freq_sdp(rhoQR, dR, k)
    fp1 = petz_fidelity_formula(X, k)
    fp2 = petz_fidelity_direct(X, rhoQR, k)
    dec, trC2 = decoupling_distance(X, k)
    rhoQRf, dRf = rhoQR_full(X, k)
    fsf, _ = freq_sdp(rhoQRf, dRf, k)
    print(f"({NR},{NB},k={k})     {fs:.6f}  {fu:.6f}   {fp1:.6f}      {fp2:.6f}      {dec:.4f} {trC2:.5f}   {fsf:.6f}")
