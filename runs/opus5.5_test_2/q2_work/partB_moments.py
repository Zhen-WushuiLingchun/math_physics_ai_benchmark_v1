"""
Exact second moments for the Stiefel code, sector by sector.
Hand-derived (Weingarten, n=2) closed forms, checked against
 (1) a brute-force Weingarten evaluation of E[V (x) V (x) conj V (x) conj V] contracted numerically,
 (2) Monte Carlo over Haar isometries.
"""
import numpy as np
from fractions import Fraction as Fr
from itertools import product
from common import *

def formulas(m, n, k):
    """closed forms (exact rationals): per sector a -> dict"""
    M = [m[i] * n[i] for i in range(2)]; D = sum(M)
    We, Ws = Fr(1, D * D - 1), Fr(-1, D * (D * D - 1))
    out = []
    for a in range(2):
        ma, na, Ma = m[a], n[a], M[a]
        al, be = ma * ma * na, ma * na * na
        rQB = Fr(1, k * k) * ((al * k + be * k * k) * We + (al * k * k + be * k) * Ws)
        rB = Fr(1, k * k) * ((al * k * k + be * k) * We + (al * k + be * k * k) * Ws)
        C2 = (1 - Fr(1, k * k)) * Fr(ma * na, D * D - 1) * (na - Fr(ma, D))
        leak = Fr((k * k - 1) * Ma * (D - Ma), D * (D * D - 1))       # E||Q_a - (TrQ_a/k) I||_2^2
        out.append(dict(TrrhoQB2=rQB, TrrhoB2=rB, TrC2=C2, leak=leak, q=Fr(Ma, D)))
    return out

def weingarten_tensor(D, k):
    """E[V_{x1 j1} V_{x2 j2} conj(V_{y1 l1}) conj(V_{y2 l2})] as array [x1,x2,y1,y2,j1,j2,l1,l2]"""
    We, Ws = 1.0 / (D * D - 1), -1.0 / (D * (D * D - 1))
    I = np.eye(D); J = np.eye(k)
    rows = {0: np.einsum('ac,bd->abcd', I, I), 1: np.einsum('ad,bc->abcd', I, I)}   # sigma=e: x1=y1,x2=y2 ; sigma=s: x1=y2,x2=y1
    cols = {0: np.einsum('ac,bd->abcd', J, J), 1: np.einsum('ad,bc->abcd', J, J)}
    T = 0
    for s in (0, 1):
        for t in (0, 1):
            w = We if s == t else Ws
            T = T + w * np.einsum('abcd,efgh->abcdefgh', rows[s], cols[t])
    return T

def brute(m, n, k):
    """evaluate E Tr rho_QB^2, E Tr rho_B^2 per sector with the brute-force tensor"""
    M = [m[i] * n[i] for i in range(2)]; D = sum(M)
    T = weingarten_tensor(D, k)
    res = []
    off = 0
    for a in range(2):
        ma, na = m[a], n[a]
        # row index of (a, r, b) = off + r*na + b
        idx = lambda r, b: off + r * na + b
        sQB = sB = 0.0
        for r, rp, b, bp in product(range(ma), range(ma), range(na), range(na)):
            x1, x2, y1, y2 = idx(r, b), idx(rp, bp), idx(r, bp), idx(rp, b)
            sub = T[x1, x2, y1, y2]       # [j1,j2,l1,l2]
            # Tr rho_QB^2 = 1/k^2 sum V_{(rb)j} V*_{(rb')j'} V*_{(r'b)j} V_{(r'b')j'}:  j1=j, j2=j', l1=j', l2=j
            sQB += np.einsum('ijji->', sub)
            # Tr rho_B^2 = 1/k^2 sum V_{(rb)j} V*_{(rb')j} V*_{(r'b)j'} V_{(r'b')j'}:   j1=j, j2=j', l1=j, l2=j'
            sB += np.einsum('ijij->', sub)
        res.append((sQB / k**2, sB / k**2))
        off += ma * na
    return res

def mc(m, n, k, S, rng):
    M = [m[i] * n[i] for i in range(2)]; D = sum(M)
    acc = np.zeros((S, 2, 4))
    for s in range(S):
        V = haar_stiefel(D, k, rng)
        A = blocks(V, m, n)
        for a in range(2):
            G = np.concatenate([A[a][j].T for j in range(k)], axis=0)       # (k n_a) x m_a
            rho = G @ G.conj().T / k                                        # rho_QB^(a)
            rB = sum(A[a][j].T @ A[a][j].conj() for j in range(k)) / k     # rho_B^(a)
            C = rho - np.kron(np.eye(k), rB) / k
            Qa = V[sum(M[:a]):sum(M[:a + 1])].conj().T @ V[sum(M[:a]):sum(M[:a + 1])]
            leak = np.linalg.norm(Qa - np.trace(Qa).real / k * np.eye(k)) ** 2
            acc[s, a] = [np.trace(rho @ rho).real, np.trace(rB @ rB).real, np.linalg.norm(C) ** 2, leak]
    return acc.mean(0), acc.std(0) / np.sqrt(S)

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    print("--- brute-force Weingarten vs closed form ---")
    for (NR, NB, k) in [(2, 3, 2), (3, 2, 3), (3, 3, 2), (3, 3, 4), (2, 4, 3), (4, 3, 2)]:
        m, n = dims(NR, NB); f = formulas(m, n, k); b = brute(m, n, k)
        for a in range(2):
            print(f"(NR,NB,k)=({NR},{NB},{k}) a={'1t'[a]}: TrrhoQB2 {float(f[a]['TrrhoQB2']):.10f} vs {b[a][0]:.10f} ; "
                  f"TrrhoB2 {float(f[a]['TrrhoB2']):.10f} vs {b[a][1]:.10f} ; TrC2 {float(f[a]['TrC2']):.10f} vs {b[a][0]-b[a][1]/k:.10f}")
    print("--- Monte Carlo (4000 Haar isometries each) ---")
    for (NR, NB, k) in [(3, 3, 2), (4, 3, 3), (3, 5, 2), (5, 4, 4), (6, 4, 2)]:
        m, n = dims(NR, NB); f = formulas(m, n, k); mu, se = mc(m, n, k, 4000, rng)
        for a in range(2):
            print(f"({NR},{NB},k={k}) a={'1t'[a]}: TrC2 exact {float(f[a]['TrC2']):.6f} MC {mu[a,2]:.6f}+-{se[a,2]:.6f} | "
                  f"TrrhoQB2 {float(f[a]['TrrhoQB2']):.6f} MC {mu[a,0]:.6f} | leak exact {float(f[a]['leak']):.3e} MC {mu[a,3]:.3e}+-{se[a,3]:.1e}")
