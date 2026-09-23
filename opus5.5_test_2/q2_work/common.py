import numpy as np
phi = (1 + 5**0.5) / 2
def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a
def dims(NR, NB):
    """m_a, n_a for a in (1, tau)"""
    m = (fib(NR - 1), fib(NR)); n = (fib(NB - 1), fib(NB))
    return m, n
d = (1.0, phi)
def haar_stiefel(D, k, rng):
    G = (rng.standard_normal((D, k)) + 1j * rng.standard_normal((D, k))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    return Q * (np.diag(R) / np.abs(np.diag(R)))
def blocks(V, m, n):
    """split rows of V (D x k) into sector blocks; returns list of arrays A[a] of shape (k, m_a, n_a)"""
    out, s = [], 0
    for ma, na in zip(m, n):
        B = V[s:s + ma * na, :]; s += ma * na
        out.append(B.T.reshape(V.shape[1], ma, na))
    assert s == V.shape[0]
    return out
