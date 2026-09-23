"""Do the Page turns of S_alg / S_qtr determine the recovery transition?
For code states (k>=1) compute E S(R), E S(B) in both conventions, I_c = S(R)-S(B) (convention independent),
H(Q|R) = -I_c and H_min(Q|R) = -ln(k F_rec)."""
import numpy as np, json
from common import *
from frec import Frec, sector_states
from multiprocessing import Pool

def ent(ev):
    ev = ev[ev > 1e-15]; return -(ev * np.log(ev)).sum()

def sample(args):
    NR, NB, k, seed = args
    rng = np.random.default_rng(seed)
    m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
    V = haar_stiefel(D, k, rng)
    st = sector_states(V, m, n, k)
    SR = SB = cen = 0.0
    for a in range(2):
        G, GR = st[a]
        rQB = G @ G.conj().T; rQR = GR @ GR.conj().T
        # rho_R^(a) = Tr_Q rho_QR^(a), rho_B^(a) = Tr_Q rho_QB^(a)
        rR = sum(rQR.reshape(k, m[a], k, m[a])[j, :, j, :] for j in range(k))
        rB = sum(rQB.reshape(k, n[a], k, n[a])[j, :, j, :] for j in range(k))
        SR += ent(np.linalg.eigvalsh(rR)); SB += ent(np.linalg.eigvalsh(rB))
        cen += np.trace(rR).real * np.log(d[a])
    F = Frec(V, m, n, k)[0] if k > 1 else 1.0
    return SR, SB, SR + cen, SB + cen, SR - SB, -np.log(k * F), F

if __name__ == "__main__":
    out = {}
    with Pool(22) as p:
        N = 20
        print(f"--- Page curves at fixed N={N}: E S(R) vs N_R, k=1 and k=2 (MC, 24 samples) ---")
        for k in [1, 2]:
            for NR in range(6, 15):
                NB = N - NR
                r = np.array(p.map(sample, [(NR, NB, k, 11 * s + NR + 1000 * k) for s in range(24)])).mean(0)
                out[f"page,{k},{NR}"] = r.tolist()
                print(f"k={k} NR={NR:2d} NB={NB:2d} Delta={NR-NB:+d} Xi={NR-NB-np.log(k)/np.log(phi):+.2f}: "
                      f"S_alg(R)={r[0]:.4f} S_qtr(R)={r[2]:.4f} | S_alg(B)={r[1]:.4f} S_qtr(B)={r[3]:.4f} | "
                      f"I_c={r[4]:+.4f} H(Q|R)={-r[4]:+.4f} H_min(Q|R)={r[5]:+.4f} F={r[6]:.4f}", flush=True)
        print("--- coherent information vs recovery fidelity at several k (N_B=8) ---")
        for k in [2, 5, 13]:
            L = np.log(k) / np.log(phi)
            for Dl in range(int(np.floor(L - 2)), int(np.ceil(L + 3)) + 1):
                NB = 8; NR = NB + Dl
                r = np.array(p.map(sample, [(NR, NB, k, 13 * s + NR + 1000 * k) for s in range(20)])).mean(0)
                out[f"ic,{k},{NR},{NB}"] = r.tolist()
                print(f"k={k:2d} Xi={Dl-L:+.3f}: I_c={r[4]:.4f}  ln k - I_c={np.log(k)-r[4]:.4f}  H_min(Q|R)={r[5]:+.4f}  F={r[6]:.4f}", flush=True)
    json.dump(out, open("entropy_vs_recovery.json", "w"), indent=1)
