"""Growing k: E F_rec vs Xi for k in {3,5,8,13,21}; compare with f_inf(c)=(E_MP sqrt(lambda))^2 and bounds."""
import numpy as np, sys, json, time
from common import *
from frec import Frec, sector_states
from multiprocessing import Pool
from scipy.special import ellipe, ellipk

def f_inf(c):
    al, be = abs(1 - np.sqrt(c)), 1 + np.sqrt(c)
    mm = 1 - al**2 / be**2
    val = be / 3 * ((al**2 + be**2) * ellipe(mm) - 2 * al**2 * (ellipk(mm) if mm < 1 else 0))
    return (val / (np.pi * c)) ** 2

def sample(args):
    NR, NB, k, seed = args
    rng = np.random.default_rng(seed)
    m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
    V = haar_stiefel(D, k, rng)
    lo, up, info = Frec(V, m, n, k)
    st = sector_states(V, m, n, k)
    mpl = 0.0; cert = 0.0
    for a in range(2):
        G = st[a][0]; rho = G @ G.conj().T
        w, U = np.linalg.eigh(rho); w = np.clip(w, 0, None)
        sq = (U * np.sqrt(w)) @ U.conj().T
        mpl += np.sqrt(w).sum() ** 2 / (k * n[a])
        T = sum(sq.reshape(k, n[a], k, n[a])[j, :, j, :] for j in range(k))   # Tr_Q sqrt(rho)
        lam = np.linalg.eigvalsh((T + T.conj().T) / 2)[-1]
        cert += (0.5 * (np.sqrt(w).sum() / np.sqrt(n[a]) + np.sqrt(n[a]) * lam)) ** 2 / k   # first-order (sigma=I/n) upper bound
    return lo, mpl, cert

if __name__ == "__main__":
    out = {}
    t0 = time.time()
    with Pool(22) as p:
        for k in [3, 5, 8, 13, 21]:
            L = np.log(k) / np.log(phi)
            for NB in [7, 9]:
                for Dl in range(int(np.floor(L - 3)), int(np.ceil(L + 4)) + 1):
                    NR = NB + Dl
                    if NR < 2: continue
                    m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
                    if k > D or D > 6e5: continue
                    S = 24
                    r = np.array(p.map(sample, [(NR, NB, k, 97 * s + 13 * NR + NB + 100000 * k) for s in range(S)]))
                    mu, se = r.mean(0), r.std(0) / np.sqrt(S)
                    Xi = Dl - L
                    out[f"{k},{NR},{NB}"] = dict(k=k, NR=NR, NB=NB, Xi=Xi, F=mu[0], Fse=se[0], MPlow=mu[1], cert=mu[2], finf=f_inf(phi ** (-Xi)))
                    print(f"k={k:2d} NR={NR:2d} NB={NB:2d} Xi={Xi:+.3f}: F={mu[0]:.4f}+-{se[0]:.4f}  MP-low={mu[1]:.4f}  1st-order-up={mu[2]:.4f}  "
                          f"f_inf(phi^-Xi)={f_inf(phi**(-Xi)):.4f}  [{time.time()-t0:.0f}s]", flush=True)
    json.dump(out, open("scan_kgrow.json", "w"), indent=1)
