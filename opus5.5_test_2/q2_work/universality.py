import numpy as np, sys, json
from common import *
from frec import Frec, sector_states
from haar_code import F_haar, bounds_from_G
from multiprocessing import Pool

def fib_sample(args):
    NR, NB, k, seed = args
    rng = np.random.default_rng(seed)
    m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
    V = haar_stiefel(D, k, rng)
    lo, up, info = Frec(V, m, n, k)
    st = sector_states(V, m, n, k)
    lw = uu = 0.0
    for a in range(2):
        l, u = bounds_from_G(st[a][0], k, n[a]); lw += l; uu += u
    # sector-wise ordinary Haar codes with the same (m_a, n_a), weighted by p_a = M_a/D
    mix = sum(m[a] * n[a] / D * F_haar(k, m[a], n[a], rng)[0] for a in range(2))
    return lo, lw, uu, mix

def haar_sample(args):
    k, m, n, seed = args
    rng = np.random.default_rng(seed)
    F, _, G = F_haar(k, m, n, rng)
    l, u = bounds_from_G(G, k, n)
    return F, l, u

if __name__ == "__main__":
    k = 2
    print("=== (1) Fibonacci code vs p_a-mixture of ordinary Haar codes with the same (m_a,n_a); MP lower / lambda_max upper bounds ===")
    with Pool(20) as p:
        for NR, NB in [(4, 3), (5, 5), (6, 4), (7, 7), (8, 6), (9, 9), (10, 8), (7, 9), (11, 10)]:
            r = np.array(p.map(fib_sample, [(NR, NB, k, 31 * s + NR * 1000 + NB) for s in range(60)]))
            mu, se = r.mean(0), r.std(0) / np.sqrt(len(r))
            Xi = NR - NB - np.log(k) / np.log(phi)
            print(f"({NR},{NB}) Xi={Xi:+.3f}: F_Fib={mu[0]:.4f}+-{se[0]:.4f} | sum_a p_a F_Haar(m_a,n_a)={mu[3]:.4f}+-{se[3]:.4f} | "
                  f"MP-lower={mu[1]:.4f}  lmax-upper={mu[2]:.4f}", flush=True)
    print("=== (2) ordinary Haar code at fixed c = k n / m, growing n  (k=2) ===")
    res = {}
    with Pool(20) as p:
        for Dl in [-3, -2, -1, 0, 1, 2, 3, 4, 5]:
            c = k * phi ** (-Dl)
            row = []
            for n in [4, 8, 16, 32, 64, 128]:
                m = int(round(k * n / c))
                if m < 1 or m * n > 3e5 or m * n < k: continue
                S = 100 if n <= 32 else 40
                r = np.array(p.map(haar_sample, [(k, m, n, 7 * s + n * 100 + Dl) for s in range(S)]))
                mu, se = r.mean(0), r.std(0) / np.sqrt(S)
                row.append((n, m, k * n / m, mu[0], se[0], mu[1], mu[2]))
                print(f"Delta={Dl:+d} c_target={c:.4f} n={n:4d} m={m:5d} c={k*n/m:.4f}: F={mu[0]:.5f}+-{se[0]:.5f}  MP-lower={mu[1]:.5f} lmax-upper={mu[2]:.5f}", flush=True)
            res[Dl] = row
    json.dump({str(k_): v for k_, v in res.items()}, open("haar_fixed_c_k2.json", "w"))
