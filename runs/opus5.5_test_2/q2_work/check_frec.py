import numpy as np, cvxpy as cp
from common import *
from frec import *

def F_sdp_hmin(V, m, n, k):
    """F_rec = (1/k) sum_a min{Tr s_a : I_k (x) s_a >= rho_QR^(a)}  (Konig-Renner-Schaffner form, CPTP decoders)"""
    st = sector_states(V, m, n, k); tot = 0.0
    for a in range(2):
        _, GR = st[a]
        rho = GR @ GR.conj().T                               # rows/cols (j, r) : Q (x) R_a
        ma = m[a]
        S = cp.Variable((ma, ma), hermitian=True)
        cons = [cp.kron(np.eye(k), S) - rho >> 0]
        pr = cp.Problem(cp.Minimize(cp.real(cp.trace(S))), cons)
        pr.solve(solver=cp.CLARABEL)
        tot += pr.value
    return tot / k

def F_uhlmann_decoder(V, m, n, k):
    """construct explicit sector-wise Uhlmann decoders from optimal sigma_B, apply them as CPTP maps, evaluate (Q2)"""
    st = sector_states(V, m, n, k); A = blocks(V, m, n)
    Phi = np.eye(k).reshape(k * k) / np.sqrt(k)              # |Phi>_{QL}, index j*k + l
    total = 0.0
    for a in range(2):
        G, GR = st[a]; ma, na = m[a], n[a]
        _, _, sig = maxH_lbfgs(G, k, na)
        # |Psi>_{Q R B} (sector a, unnormalised): coefficient [j, r, b] = A_j[r,b]/sqrt(k)
        Psi = np.stack([A[a][j] for j in range(k)]) / np.sqrt(k)
        # target |Phi>_{QL} (x) |sqrt(sigma)>_{B E}, coefficient [j, l, b, e]
        w, U = np.linalg.eigh(sig); ss = (U * np.sqrt(np.clip(w, 0, None))) @ U.conj().T
        Tgt = np.einsum('jl,be->jlbe', np.eye(k) / np.sqrt(k), ss)
        # overlap operator  M[(l,e), r] = sum_{j,b} conj(Tgt[j,l,b,e]) Psi[j,r,b]
        Mo = np.einsum('jlbe,jrb->ler', Tgt.conj(), Psi).reshape(k * na, ma)
        Uu, s, Wh = np.linalg.svd(Mo, full_matrices=False)
        W = (Uu @ Wh).conj()                                  # optimal (partial) isometry R -> L E  (k na x ma)
        # CPTP completion: D(X) = Tr_E(W X W^dag) + Tr((I - W^dag W) X) I/k
        comp = np.eye(ma) - W.conj().T @ W
        # apply id_Q (x) D to rho_QR^(a)
        rhoQR = GR @ GR.conj().T                              # [(j,r),(j',r')]
        R4 = rhoQR.reshape(k, ma, k, ma)
        W3 = W.reshape(k, na, ma)                             # [l, e, r]
        out = np.einsum('ler,jrts,mes->jlmt', W3, R4.transpose(0, 1, 2, 3), W3.conj()) if False else None
        # (id (x) D)(rho)[(j,l),(j',l')] = sum_{r,r',e} W[l,e,r] rho[(j,r),(j',r')] conj(W[l',e,r'])
        out = np.einsum('ler,jrJs,Les->jlJL', W3, R4, W3.conj()).reshape(k * k, k * k)
        out += np.kron(np.einsum('jrJs,sr->jJ', R4, comp), np.eye(k) / k)
        # check trace preservation of the decoder on this block
        assert abs(np.trace(out).real - np.trace(rhoQR).real) < 1e-10
        total += (Phi.conj() @ out @ Phi).real
    return total

if __name__ == "__main__":
    rng = np.random.default_rng(11)
    for (NR, NB, k) in [(2, 2, 2), (3, 2, 2), (2, 3, 2), (3, 3, 2), (4, 3, 2), (3, 4, 3), (4, 4, 2), (5, 3, 3), (4, 5, 2), (5, 5, 4)]:
        m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
        if k > D: continue
        V = haar_stiefel(D, k, rng)
        lo, up, info = Frec(V, m, n, k)
        fs = F_sdp_hmin(V, m, n, k)
        fu = F_uhlmann_decoder(V, m, n, k)
        print(f"({NR},{NB}) k={k} D={D}: max_sigma form [{lo:.8f}, {up:.8f}] | CPTP SDP (H_min) {fs:.8f} | explicit Uhlmann decoder {fu:.8f}")
