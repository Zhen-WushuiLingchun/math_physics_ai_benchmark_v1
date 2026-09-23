"""
SDP solver for optimal entanglement-recovery fidelity F_rec(V).

Conventions:
  Channel Lambda: L(=C^k) -> R(=C^dR), given as a linear map on matrices.
  rhoQR[(q r),(q' r')] = (1/k) [Lambda(|q><q'|)]_{r r'}.
  Choi J on R x L: J[(r,l),(r',l')] = <l| D(|r><r'|) |l'>;  PSD; Tr_L J = I_R.
  F = (1/k) Tr[J W],  W[(a,q),(a',q')] = rho[(q,a'),(q',a)].
"""
import numpy as np
import cvxpy as cp

def W_from_rho(rho, dR, k):
    """W[a*k+l, a'*k+l'] = rho[(l',a'),(l,a)]  gives (1/k)Tr[W J] = <Phi|(idxD)(rho)|Phi>."""
    W = np.zeros((dR*k, dR*k), dtype=complex)
    for l in range(k):
        for a in range(dR):
            for lp in range(k):
                for ap in range(dR):
                    W[a*k+l, ap*k+lp] = rho[lp*dR+ap, l*dR+a]
    return W

def rho_from_channel(chan, k, dR):
    rho = np.zeros((k*dR, k*dR), dtype=complex)
    for q in range(k):
        for qp in range(k):
            X = np.zeros((k, k), dtype=complex); X[q, qp] = 1
            rho.reshape(k, dR, k, dR)[q, :, qp, :] = chan(X)/k
    return rho

def freq_sdp(rho, dR, k, verbose=False):
    W = W_from_rho(rho, dR, k)
    J = cp.Variable((dR*k, dR*k), hermitian=True)
    Al = [np.kron(np.eye(dR), np.eye(k)[l]) for l in range(k)]
    psum = sum(Al[l] @ J @ Al[l].conj().T for l in range(k))
    prob = cp.Problem(cp.Maximize(cp.real(cp.trace(W @ J))/k),
                      [J >> 0, psum == np.eye(dR)])
    prob.solve(solver=cp.SCS, eps=1e-9, max_iters=60000, verbose=verbose)
    return float(prob.value), J.value

def fidelity_of_decoder(rho, dR, k, Dfun):
    T = rho.reshape(k, dR, k, dR)
    F = 0.0
    for r in range(dR):
        for rp in range(dR):
            X = np.zeros((dR, dR), dtype=complex); X[r, rp] = 1
            Dr = Dfun(X)
            for q in range(k):
                for qp in range(k):
                    F += T[q, r, qp, rp]*Dr[q, qp]
    return float(np.real(F))/k

def check(name, chan, dR, k, decoder=None, expect=None):
    rho = rho_from_channel(chan, k, dR)
    f, J = freq_sdp(rho, dR, k)
    msg = f"{name:38s}: F_sdp = {f:.6f}"
    if decoder is not None:
        fd = fidelity_of_decoder(rho, dR, k, decoder)
        msg += f"  (given: {fd:.6f})"
    if expect is not None:
        msg += f"  expect {expect:.6f}"
    print(msg)
    return f, rho

rng = np.random.default_rng(11)

# 1. identity channel R=L
k = 3; dR = 3
check("identity channel", lambda X: X, dR, k, decoder=lambda X: X, expect=1.0)

# 2. depolarizing
check("fully depolarizing", lambda X: np.eye(dR)*np.trace(X)/dR, dR, k,
      decoder=lambda X: np.eye(k)*np.trace(X)/k, expect=1/9)

# 3. random unitary channel
U, _ = np.linalg.qr(rng.normal(size=(dR, dR)) + 1j*rng.normal(size=(dR, dR)))
check("random unitary channel", lambda X: U @ X @ U.conj().T, dR, k,
      decoder=lambda X: U.conj().T @ X @ U, expect=1.0)

# 4. random isometry + partial trace of B
k, dRA, dRB = 2, 3, 2
V, _ = np.linalg.qr(rng.normal(size=(dRA*dRB, k)) + 1j*rng.normal(size=(dRA*dRB, k)))
def chan_iso(X):
    Y = V @ X @ V.conj().T
    return Y.reshape(dRA, dRB, dRA, dRB).trace(axis1=1, axis2=3)
check("random isometry, Tr_B, dR=3, k=2", chan_iso, dRA, k)

# 5. random small channel: full CPTP via isometry into R x E, trace out E
k, dR, dE = 2, 4, 3
Wm, _ = np.linalg.qr(rng.normal(size=(dR*dE, k)) + 1j*rng.normal(size=(dR*dE, k)))
def chan_gen(X):
    Y = Wm @ X @ Wm.conj().T
    return Y.reshape(dR, dE, dR, dE).trace(axis1=1, axis2=3)
check("random chan (4x3 env), k=2", chan_gen, dR, k)

# 6. exact noiseless: R contains a copy of L + junk
k = 2; dR = 4
Wm, _ = np.linalg.qr(rng.normal(size=(k, k)) + 1j*rng.normal(size=(k, k)))
# isometry L -> R:  |q> -> U|q> (x) |junk>
Junk = np.zeros((k, 1)); Junk[0, 0] = 1
Vemb = np.kron(Wm, Junk)
check("exact embedded unitary, k=2, dR=4", lambda X: Vemb @ X @ Vemb.conj().T, dR, k,
      decoder=lambda Y: Wm.conj().T @ Y[:k, :k] @ Wm, expect=1.0)
