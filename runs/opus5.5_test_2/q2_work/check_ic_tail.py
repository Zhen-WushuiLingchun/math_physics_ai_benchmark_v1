"""Tail check (ordinary Haar code, k=2): I_c + ln k  vs (k^2-1)/(2c), and F - 1/k^2 vs s_k/(k sqrt c)."""
import numpy as np
from common import haar_stiefel
from haar_code import F_haar
rng = np.random.default_rng(21)
def ent(ev):
    ev = ev[ev > 1e-15]; return -(ev * np.log(ev)).sum()
k = 2
for m, n in [(8, 64), (8, 128), (12, 384), (8, 400)]:
    c = k * n / m
    Ic, F = [], []
    for _ in range(12):
        V = haar_stiefel(m * n, k, rng)
        A = V.T.reshape(k, m, n) / np.sqrt(k)
        rR = sum(A[j] @ A[j].conj().T for j in range(k)); rB = sum(A[j].T @ A[j].conj() for j in range(k))
        Ic.append(ent(np.linalg.eigvalsh(rR)) - ent(np.linalg.eigvalsh(rB)))
        if n <= 128: F.append(F_haar(k, m, n, rng)[0])
    Ic = np.mean(Ic)
    s = f"  F-1/4={np.mean(F)-0.25:.4f} (s_k/k)/sqrt(c)|s=1: {0.5/np.sqrt(c):.4f}" if F else ""
    print(f"m={m} n={n} c={c:.1f}: I_c+ln k = {Ic+np.log(k):.5f}   (k^2-1)/(2c) = {(k*k-1)/(2*c):.5f}{s}")
