"""Part C scan: E F_rec for k=2 over many (N_R, N_B) incl. unequal cuts; compare with Xi-collapse."""
import numpy as np, sys, time, json
from common import *
from frec import Frec
from bounds import decoupling_bound, collision_bound
from multiprocessing import Pool

def one(args):
    NR, NB, k, seed = args
    rng = np.random.default_rng(seed)
    m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
    V = haar_stiefel(D, k, rng)
    lo, up, info = Frec(V, m, n, k)
    return lo, up, info[0]['F_a'], info[1]['F_a'], info[0]['q']

def run(NR, NB, k, S, base=0):
    with Pool(20) as p:
        r = np.array(p.map(one, [(NR, NB, k, base + 1000 * NR + 17 * NB + 7919 * s + k) for s in range(S)]))
    return r

if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    cases = []
    for Dl in range(-4, 8):
        for NB in range(2, 14):
            NR = NB + Dl
            if NR < 2: continue
            m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
            if k > D or max(n) > 150 or D > 4e5: continue
            cases.append((NR, NB))
    res = {}
    t0 = time.time()
    for NR, NB in cases:
        m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
        S = 48 if max(n) > 50 else 120
        r = run(NR, NB, k, S)
        Xi = NR - NB - np.log(k) / np.log(phi)
        Ed, dec = decoupling_bound(m, n, k); col = collision_bound(m, n, k)
        res[f"{NR},{NB}"] = dict(NR=NR, NB=NB, Xi=Xi, D=D, S=S, F=r[:, 0].mean(), Fse=r[:, 0].std() / np.sqrt(S),
                                 Fsd=r[:, 0].std(), gap=(r[:, 1] - r[:, 0]).max(), F1=r[:, 2].mean(), Ft=r[:, 3].mean(),
                                 dec=dec, col=col, Edelta=Ed)
        e = res[f"{NR},{NB}"]
        print(f"NR={NR:2d} NB={NB:2d} Xi={Xi:+.3f} D={D:7d}  E F={e['F']:.5f}+-{e['Fse']:.5f} (sd {e['Fsd']:.4f}, cert.gap {e['gap']:.1e})"
              f"  F_1={e['F1']:.4f} F_tau={e['Ft']:.4f} | rigorous: collision>={col:.4f} decoupling>={dec:.4f}  [{time.time()-t0:.0f}s]", flush=True)
    json.dump(res, open(f"scan_k{k}.json", "w"), indent=1)
