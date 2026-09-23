"""Exact and numerical checks for the Fibonacci recovery problem."""
import numpy as np
from numpy.linalg import svd, norm, eigvalsh

PHI = (1 + np.sqrt(5)) / 2

def fib(n):
    """Fib_0=0, Fib_1=1, Fib_2=1, ..."""
    if n < 0:
        raise ValueError(n)
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

def dims(NR, NB):
    m1, mt = fib(NR - 1), fib(NR)
    n1, nt = fib(NB - 1), fib(NB)
    D = m1 * n1 + mt * nt
    return m1, mt, n1, nt, D

def E_trC2(NR, NB, k):
    m1, mt, n1, nt, D = dims(NR, NB)
    P = m1**2 * n1 + mt**2 * nt
    Q = m1 * n1**2 + mt * nt**2
    num = D * Q - P
    return (k**2 - 1) / (k**2 * D * (D**2 - 1)) * num

def haar_isometry(D, k, rng):
    G = rng.normal(size=(D, k)) + 1j * rng.normal(size=(D, k))
    Q, _ = np.linalg.qr(G)
    # qr of complex may flip phases; columns are Haar
    return Q[:, :k]

def pack_indices(m1, mt, n1, nt):
    """Map flat index -> (sector, alpha, beta)."""
    idx = []
    for a, m, n in ((0, m1, n1), (1, mt, nt)):
        for alpha in range(m):
            for beta in range(n):
                idx.append((a, alpha, beta, m, n))
    return idx

def trC2_from_V(V, m1, mt, n1, nt, k):
    """Direct Tr C^2 from an isometry V (D,k), for validation."""
    D = V.shape[0]
    # Build M_a^{jj'} as n_a x n_a matrices
    # c_j(a,alpha,beta) = V[flat, j]
    secs = []
    offset = 0
    for m, n in ((m1, n1), (mt, nt)):
        # reshape this sector's coefficients to (k, m, n)
        block = V[offset:offset + m * n, :]
        # block[alpha*n+beta, j]
        C = block.reshape(m, n, k)  # careful: numpy reshape is C-order alpha slow?
        # flat = alpha*n + beta if alpha is outer. Our pack uses alpha outer then beta.
        secs.append(C)  # shape (m, n, k)
        offset += m * n
    # M^{jj'}_{beta beta'} = sum_alpha C[alpha,beta,j] * conj(C[alpha,beta',j'])
    def tr_MM(j, jp):
        s = 0.0
        for C in secs:
            # (n,k) after contracting alpha: sum_alpha C[a,b,j] conj(C[a,bp,jp])
            # einsum
            M = np.einsum("abj,abp->bp", C[:, :, j:j+1], np.conj(C[:, :, jp:jp+1]), optimize=True)
            # M shape (n, n) wait einsum "abj,acp->bc"
        return s

    # redo cleaner
    tot_cross = 0.0 + 0.0j
    tot_diag = 0.0 + 0.0j
    Ms = []  # list over a of M[j,jp]
    for C in secs:
        # C[alpha,beta,j]
        # M[j,jp,beta,betap]
        M = np.einsum("abj,acp->jpcb", C, np.conj(C), optimize=True)
        Ms.append(M)
    for j in range(k):
        for jp in range(k):
            tr = 0.0 + 0.0j
            for M in Ms:
                tr += np.trace(M[j, jp] @ M[jp, j])
            tot_cross += tr
    for j in range(k):
        for l in range(k):
            tr = 0.0 + 0.0j
            for M in Ms:
                tr += np.trace(M[j, j] @ M[l, l])
            tot_diag += tr
    val = tot_cross / k**2 - tot_diag / k**3
    return float(np.real(val))

def monte_carlo_trC2(NR, NB, k, trials, seed=0):
    m1, mt, n1, nt, D = dims(NR, NB)
    rng = np.random.default_rng(seed)
    acc = 0.0
    for _ in range(trials):
        V = haar_isometry(D, k, rng)
        acc += trC2_from_V(V, m1, mt, n1, nt, k)
    return acc / trials

def petz_fidelity(V, m1, mt, n1, nt, k):
    """Achievable entanglement fidelity of the Petz decoder at I/k.
    Returns a number in [0,1]. Lower bound on F_rec.
    """
    # Build Kraus-like action of N on basis operators |i><j|
    # N(|i><j|) = ⊕_a sum_beta E_beta |i><j| E_beta^dagger, an operator on ⊕ R^a
    # We compute the k x k matrix of the decoded channel's entanglement fidelity directly.
    # R(X) = sum_a V_a^dagger (Gamma_a^{-1/2} X_a Gamma_a^{-1/2} ⊗ I) V_a
    # F = k^{-2} sum_{ij} <i| R(N(|i><j|)) |j>
    secs = []
    offset = 0
    for m, n in ((m1, n1), (mt, nt)):
        if m == 0 or n == 0:
            secs.append(None)
            continue
        block = V[offset:offset + m * n, :]  # (m*n, k), flat=alpha*n+beta
        C = block.reshape(m, n, k)
        secs.append(C)
        offset += m * n
    # Gamma_a = sum_j sum_beta |phi><phi| = einsum C C*
    gammas = []
    for C in secs:
        if C is None:
            gammas.append(None)
            continue
        # Gamma[alpha,alphap] = sum_{beta,j} C[a,b,j] conj(C[ap,b,j])
        G = np.einsum("abj,cbj->ac", C, np.conj(C), optimize=True)
        # inverse square root on support
        ev, U = np.linalg.eigh((G + G.conj().T) / 2)
        ev = np.clip(ev.real, 0, None)
        inv_sqrt = np.zeros_like(ev)
        mask = ev > 1e-12
        inv_sqrt[mask] = 1 / np.sqrt(ev[mask])
        Gmh = (U * inv_sqrt) @ U.conj().T
        gammas.append(Gmh)
    # For each i,j compute <i|R(N(|i><j|))|j>
    # N(|i><j|)_a [alpha,alphap] = sum_beta C[a,b,i] conj(C[ap,b,j])
    acc = 0.0 + 0.0j
    for i in range(k):
        for j in range(k):
            # R output matrix element <p| R |q> summed for p=i, q=j?
            # <i| R(X) |j> = sum_a <psi_i^a| (Gmh X Gmh ⊗ I) |psi_j^a>
            # = sum_a sum_{alpha,alphap,beta} conj(C[a,b,i]) (Gmh X Gmh)[a,ap] C[ap,b,j]
            val = 0.0 + 0.0j
            for C, Gmh in zip(secs, gammas):
                if C is None:
                    continue
                X = np.einsum("ab,cb->ac", C[:, :, i], np.conj(C[:, :, j]), optimize=True)
                Y = Gmh @ X @ Gmh
                val += np.einsum("ab,ac,cb->", np.conj(C[:, :, i]), Y, C[:, :, j], optimize=True)
            acc += val
    F = acc / k**2
    return float(np.real(F))

def gamma_sqrt_bound(V, m1, mt, n1, nt, k):
    """(Tr sqrt(Gamma))^2 / k^2, candidate upper bound, Gamma=Tr_B Pi."""
    secs = []
    offset = 0
    s = 0.0
    for m, n in ((m1, n1), (mt, nt)):
        if m == 0 or n == 0:
            offset += m * n
            continue
        block = V[offset:offset + m * n, :]
        C = block.reshape(m, n, k)
        G = np.einsum("abj,cbj->ac", C, np.conj(C), optimize=True)
        ev = eigvalsh((G + G.conj().T) / 2).real
        ev = np.clip(ev, 0, None)
        s += np.sum(np.sqrt(ev))
        offset += m * n
    return (s / k) ** 2

def main():
    print("=== dimension checks ===")
    for N in range(2, 12):
        for NR in range(1, N):
            NB = N - NR
            m1, mt, n1, nt, D = dims(NR, NB)
            assert D == fib(N - 1), (N, NR, D, fib(N - 1))
    print("D_N identity OK")

    print("=== E Tr C^2 vs Monte Carlo ===")
    rng_cases = [(2, 2, 2), (2, 3, 2), (3, 3, 2), (3, 4, 2), (4, 4, 2)]
    for NR, NB, k in rng_cases:
        exact = E_trC2(NR, NB, k)
        mc = monte_carlo_trC2(NR, NB, k, trials=400, seed=1)
        print(f"({NR},{NB}) k={k}: exact={exact:.6f} mc={mc:.6f}")

    print("=== Petz fidelity vs (Tr sqrt Gamma)^2/k^2 ===")
    rng = np.random.default_rng(0)
    for NR, NB, k, trials in [(2, 2, 2, 30), (2, 3, 2, 30), (3, 3, 2, 20), (3, 4, 2, 15), (4, 3, 2, 15), (4, 5, 2, 8)]:
        m1, mt, n1, nt, D = dims(NR, NB)
        if D < k:
            continue
        Fs, Bs = [], []
        for _ in range(trials):
            V = haar_isometry(D, k, rng)
            Fs.append(petz_fidelity(V, m1, mt, n1, nt, k))
            Bs.append(gamma_sqrt_bound(V, m1, mt, n1, nt, k))
        Fs, Bs = np.array(Fs), np.array(Bs)
        xi = NR - NB - np.log(k) / np.log(PHI)
        print(
            f"({NR},{NB}) k={k} D={D} dR={m1+mt} Xi={xi:.3f} "
            f"Petz mean={Fs.mean():.4f} min={Fs.min():.4f} "
            f"bound mean={Bs.mean():.4f}  Petz>bound {np.mean(Fs>Bs+1e-6):.2f} "
            f"EtrC2={E_trC2(NR,NB,k):.4e}"
        )

    print("=== scaling of HS trace-distance proxy vs Xi ===")
    # proxy = 0.5 * sqrt(k * sum n) * sqrt(E Tr C^2)
    for k in (2, 3, 5):
        for delta in range(-4, 6):
            # N_B=40, N_R=N_B+delta if positive dims
            NB = 30
            NR = NB + delta
            if NR < 2:
                continue
            m1, mt, n1, nt, D = dims(NR, NB)
            if k > D:
                continue
            etc = E_trC2(NR, NB, k)
            proxy = 0.5 * np.sqrt(k * (n1 + nt) * etc)
            Xi = NR - NB - np.log(k) / np.log(PHI)
            print(f"k={k} delta={delta} Xi={Xi:.3f} proxy={proxy:.6f} EtrC2={etc:.3e} dR={m1+mt}")

if __name__ == "__main__":
    main()
