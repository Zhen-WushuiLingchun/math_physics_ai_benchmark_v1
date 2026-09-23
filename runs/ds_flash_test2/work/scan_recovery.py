"""
Main numerical experiment:
  For Fibonacci splits (NR,NB) and logical dim k:
    - sample Haar-random isometry V of dim D=F_{N-1}, k columns
    - build rho_QR (input to SDP), rho_QB, C_QB, compute ||C_QB||_1, Tr C_QB^2, ||C_QR||_1
    - solve SDP for optimal F_rec
    - compute Petz (transpose-channel) fidelity
  Save results to CSV.
"""
import numpy as np, math, csv, sys, time
import cvxpy as cp
from itertools import product

rng = np.random.default_rng(20260923)
phi = (1+math.sqrt(5))/2
logphi = math.log(phi)

def fib(n):
    a, b = 0, 1
    for _ in range(n): a, b = b, a+b
    return a

def specs(NR, NB):
    m = [fib(NR-1), fib(NR)]      # a = 1, tau  (index 0,1)
    n = [fib(NB-1), fib(NB)]
    return m, n

def rand_isometry(D, k, rng):
    A = (rng.normal(size=(D, k)) + 1j*rng.normal(size=(D, k)))/np.sqrt(2)
    Q, R = np.linalg.qr(A)
    d = np.diagonal(R).copy(); d /= np.abs(d)
    return Q * d.conj()

def blocks(V, m, n):
    X = []; off = 0
    for ma, na in zip(m, n):
        X.append(V[off:off+ma*na, :].reshape(ma, na, -1).transpose(2, 0, 1))  # (k,ma,na)
        off += ma*na
    return X

def rho_QR(X, k):
    m = [x.shape[1] for x in X]; dR = sum(m)
    rho = np.zeros((k*dR, k*dR), dtype=complex)
    off = 0
    for x, ma in zip(X, m):
        blk = np.zeros((k, k, ma, ma), dtype=complex)
        for q in range(k):
            for qp in range(k):
                blk[q, qp] = x[q] @ x[qp].conj().T
        for q in range(k):
            for qp in range(k):
                rho[q*dR+off:q*dR+off+ma, qp*dR+off:qp*dR+off+ma] = blk[q, qp]/k
        off += ma
    return rho, dR

def rho_QB(X, k):
    n = [x.shape[2] for x in X]; dB = sum(n)
    rho = np.zeros((k*dB, k*dB), dtype=complex)
    off = 0
    for x, na in zip(X, n):
        blk = np.zeros((k, k, na, na), dtype=complex)
        for q in range(k):
            for qp in range(k):
                blk[q, qp] = x[q].T @ x[qp].conj()
        for q in range(k):
            for qp in range(k):
                rho[q*dB+off:q*dB+off+na, qp*dB+off:qp*dB+off+na] = blk[q, qp]/k
        off += na
    return rho, dB

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
    return float(prob.value)

def petz_fidelity(X, rhoQR, k):
    """Petz/transpose-channel decoder for Lambda with respect to I/k; fidelity on rhoQR."""
    m = [x.shape[1] for x in X]; dR = sum(m)
    # precompute S_a^{-1/2} and K^{(a)}_nu : K_nu[:, q] = X[a][q,:,nu] (m x k)
    Slist = []
    for x, ma in zip(X, m):
        S = sum(x[q] @ x[q].conj().T for q in range(k))
        w, U = np.linalg.eigh(S)
        Sinv = (U * (1/np.sqrt(np.maximum(w, 1e-14)))) @ U.conj().T
        Slist.append(Sinv)
    def Dfun(Y):
        out = np.zeros((k, k), dtype=complex)
        off = 0
        for x, ma, Sinv in zip(X, m, Slist):
            Ya = Y[off:off+ma, off:off+ma]
            Z = Sinv @ Ya @ Sinv
            for nu in range(x.shape[2]):
                K = x[:, :, nu].T  # (m,k)
                out += K.conj().T @ Z @ K
            off += ma
        return out
    T = rhoQR.reshape(k, dR, k, dR)
    F = 0.0
    for r in range(dR):
        for rp in range(dR):
            Xm = np.zeros((dR, dR), dtype=complex); Xm[r, rp] = 1
            Dr = Dfun(Xm)
            for q in range(k):
                for qp in range(k):
                    F += T[q, r, qp, rp]*Dr[q, qp]
    return float(np.real(F))/k

def diagnostics(X, rhoQB, dB, k):
    rhoB = rhoQB.reshape(k, dB, k, dB).trace(axis1=0, axis2=2)  # Tr_Q
    C = rhoQB - np.kron(np.eye(k)/k, rhoB)
    trC2 = float(np.real(np.trace(C @ C)))
    n1 = float(np.sum(np.abs(np.linalg.eigvalsh(C))))
    return trC2, n1

# ---------------- run grid ----------------
Nsamp = int(sys.argv[1]) if len(sys.argv) > 1 else 6
out = []
configs = []
for NR, NB, kmax in [(2,2,2),(3,2,3),(3,3,4),(4,2,4),(4,3,5),(4,4,3),(5,2,5),(5,3,5),(5,4,5),(5,5,4),
                     (6,2,5),(6,3,5),(6,4,4),(6,5,4),(6,6,3),(7,2,4),(7,3,4),(7,4,4),(7,5,3),(7,6,2),
                     (8,2,4),(8,3,3),(8,4,3),(8,5,2),(9,2,3),(9,3,3),(9,4,2)]:
    D = fib(NR+NB-1)
    for k in range(2, min(D, kmax)+1):
        m, n = specs(NR, NB)
        dR = sum(m)
        if dR*k > 56: continue
        configs.append((NR,NB,k,D,dR))

print(f"{'NR':>3} {'NB':>3} {'k':>2} {'D':>4} {'xi':>7} {'TrC2':>11} {'||C||1':>8} {'F_sdp':>8} {'F_petz':>8} {'TrC2_qr':>11} {'||Cqr||1':>8}")
for (NR,NB,k,D,dR) in configs:
    xi = NR - NB - math.log(k)/logphi
    tr2 = n1 = fq = fp = tr2r = n1r = 0.0
    for s in range(Nsamp):
        V = rand_isometry(D, k, rng)
        X = blocks(V, *specs(NR,NB))
        rhoQR, _ = rho_QR(X, k)
        rhoQB, dB = rho_QB(X, k)
        t2, nn = diagnostics(X, rhoQB, dB, k)
        # complementary (Q vs R) diagnostics: swap roles by transposing blocks
        Xt = [x.transpose(0,2,1) for x in X]
        rhoQR2, dR2 = rho_QB(Xt, k)   # rho_QB of transposed = rho_QR of original
        t2r, nnr = diagnostics(Xt, rhoQR2, dR2, k)
        f = freq_sdp(rhoQR, dR, k)
        p = petz_fidelity(X, rhoQR, k)
        tr2 += t2; n1 += nn; fq += f; fp += p; tr2r += t2r; n1r += nnr
    tr2/=Nsamp; n1/=Nsamp; fq/=Nsamp; fp/=Nsamp; tr2r/=Nsamp; n1r/=Nsamp
    print(f"{NR:>3} {NB:>3} {k:>2} {D:>4} {xi:>7.3f} {tr2:>11.6f} {n1:>8.4f} {fq:>8.4f} {fp:>8.4f} {tr2r:>11.6f} {n1r:>8.4f}", flush=True)
    out.append(dict(NR=NR,NB=NB,k=k,D=D,dR=dR,xi=xi,TrC2=tr2,C1=n1,Fsdp=fq,Fpetz=fp,TrC2qr=tr2r,C1qr=n1r))

with open("work/results_recovery.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
print("saved work/results_recovery.csv")
