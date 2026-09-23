"""Tail constant s_k: F_rec - 1/k^2 ~ (s_k/k) c^{-1/2} as c -> inf (fixed k).
s_k = lim_m (1/m) min{ Tr S : I_k (x) S >= W },  W = (km x km) GUE normalised to semicircle on [-2,2]."""
import numpy as np, cvxpy as cp, warnings
warnings.filterwarnings("ignore")
rng = np.random.default_rng(5)
def gue(N):
    A = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2)
    return (A + A.conj().T) / np.sqrt(2 * N)          # semicircle radius 2
for k in [2, 3]:
    for m in [8, 16, 32, 64]:
        vals = []
        S_ = 12 if m <= 32 else 4
        for _ in range(S_):
            W = gue(k * m)
            S = cp.Variable((m, m), hermitian=True)
            pr = cp.Problem(cp.Minimize(cp.real(cp.trace(S))), [cp.kron(np.eye(k), S) - W >> 0])
            pr.solve(solver=cp.CLARABEL)
            vals.append(pr.value / m)
        print(f"k={k} m={m:3d}: (1/m) min Tr S = {np.mean(vals):.4f} +- {np.std(vals)/np.sqrt(len(vals)):.4f}", flush=True)
