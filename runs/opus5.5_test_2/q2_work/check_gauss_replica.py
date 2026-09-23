import numpy as np, mpmath as mp
from common import *
from partA import exact_A
rng = np.random.default_rng(3)
for NR, NB in [(3, 3), (4, 4), (5, 3)]:
    m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
    S = 200000; acc_num = acc_den = 0.0
    for _ in range(S):
        g = (rng.standard_normal(D) + 1j * rng.standard_normal(D)) / np.sqrt(2 * D)
        X = blocks(g[:, None], m, n)
        for i in range(2):
            ev = np.linalg.eigvalsh(X[i][0] @ X[i][0].conj().T); ev = ev[ev > 1e-300]
            acc_num += (ev * np.log(ev)).sum(); acc_den += ev.sum()
    lhs = -acc_num / acc_den                               # -d/dn ln E Tr rho_g^n at n=1
    ES = float(exact_A(NR, NB)['Salg'])
    rhs = ES - float(mp.digamma(D + 1) - mp.log(D))
    print(f"({NR},{NB}) D={D}: Gaussian replica (MC) {lhs:.5f}  vs  E S_alg - [psi(D+1)-ln D] = {rhs:.5f}   (E S_alg = {ES:.5f})")
