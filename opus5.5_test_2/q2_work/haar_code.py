"""Ordinary (sector-free) Haar random code C^k -> C^m (x) C^n, recovery from the m-factor."""
import numpy as np
from common import *
from frec import maxH_lbfgs

def F_haar(k, m, n, rng):
    D = m * n
    V = haar_stiefel(D, k, rng)
    A = V.T.reshape(k, m, n)
    G = np.concatenate([A[j].T for j in range(k)], axis=0) / np.sqrt(k)
    hl, hu, sig = maxH_lbfgs(G, k, n)
    return hl**2 / k, hu**2 / k, G

def bounds_from_G(G, k, n):
    """per-code: MP-type lower bound (sigma = I/n) and lambda_max converse upper bound"""
    rho = G @ G.conj().T; q = np.trace(rho).real
    ev = np.clip(np.linalg.eigvalsh(rho), 0, None)
    low = (np.sqrt(ev).sum()) ** 2 / (k * n)                  # = q * F(rho_hat, I/(kn))
    rB = sum(G.reshape(k, n, -1)[j] @ G.reshape(k, n, -1)[j].conj().T for j in range(k))
    rank = min(G.shape)
    up = min(q, rank * np.linalg.eigvalsh(rB)[-1] / k)        # F <= rank(rho_QB) * lambda_max(rho_B) / k  (unnormalised form)
    return low, up
