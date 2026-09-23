import numpy as np, mpmath as mp
from common import *
mp.mp.dps = 40
H = lambda n: mp.harmonic(n)
def page(m, n):                       # Page's exact mean entropy for C^m (x) C^n
    a, b = min(m, n), max(m, n)
    return H(a * b) - H(b) - mp.mpf(a - 1) / (2 * b)
def exact_A(NR, NB):
    m, n = dims(NR, NB); M = [m[i] * n[i] for i in range(2)]; D = sum(M)
    EH = H(D) - sum(mp.mpf(M[i]) / D * H(M[i]) for i in range(2))       # E H(p), p ~ Dirichlet(M_1, M_tau)
    Salg = EH + sum(mp.mpf(M[i]) / D * page(m[i], n[i]) for i in range(2))
    Sqtr = Salg + mp.mpf(M[1]) / D * mp.log(mp.phi)
    Z2alg = sum(mp.mpf(M[i]) * (m[i] + n[i]) for i in range(2)) / (D * (D + 1))
    Z2qtr = sum(mp.mpf(M[i]) * (m[i] + n[i]) / [1, mp.phi][i] for i in range(2)) / (D * (D + 1))
    return dict(D=D, EH=EH, Salg=Salg, Sqtr=Sqtr, Z2alg=Z2alg, Z2qtr=Z2qtr)

def mc_A(NR, NB, S=20000, seed=1):
    rng = np.random.default_rng(seed)
    m, n = dims(NR, NB); D = sum(m[i] * n[i] for i in range(2))
    out = {k: [] for k in ["Salg", "Sqtr", "SalgB", "SqtrB", "Z2alg", "Z2qtr", "p1", "Z2g", "g4"]}
    for _ in range(S):
        g = (rng.standard_normal(D) + 1j * rng.standard_normal(D)) / np.sqrt(2 * D)   # E|g|^2 = 1
        nrm2 = np.vdot(g, g).real
        psi = g / np.sqrt(nrm2)
        X = blocks(psi[:, None], m, n)
        sa = sq = sb = z2a = z2q = z2g = 0.0
        for i in range(2):
            Xa = X[i][0]
            ev = np.linalg.eigvalsh(Xa @ Xa.conj().T); ev = ev[ev > 1e-300]
            evB = np.linalg.eigvalsh(Xa.T @ Xa.conj()); evB = evB[evB > 1e-300]
            sa += -(ev * np.log(ev)).sum(); sb += -(evB * np.log(evB)).sum()
            sq += -(ev * np.log(ev / d[i])).sum()
            z2a += (ev**2).sum(); z2q += (ev**2).sum() / d[i]
            z2g += (ev**2).sum() * nrm2**2
        pt = np.linalg.norm(X[1]) ** 2
        out["Salg"].append(sa); out["Sqtr"].append(sq); out["SalgB"].append(sb)
        out["SqtrB"].append(sb + pt * np.log(phi))
        out["Z2alg"].append(z2a); out["Z2qtr"].append(z2q); out["p1"].append(1 - pt)
        out["Z2g"].append(z2g); out["g4"].append(nrm2**2)
    return {k: np.array(v) for k, v in out.items()}

if __name__ == "__main__":
    print("=== exact vs Monte Carlo (20000 Haar states each) ===")
    for NR, NB in [(2, 2), (3, 2), (2, 4), (4, 4), (5, 3), (3, 6), (6, 6)]:
        e = exact_A(NR, NB); r = mc_A(NR, NB)
        m, n = dims(NR, NB); M = [m[i] * n[i] for i in range(2)]; D = sum(M)
        se = lambda x: x.std() / np.sqrt(len(x))
        print(f"(NR,NB)=({NR},{NB}) m={m} n={n} D={D}")
        print(f"   E S_alg : exact {float(e['Salg']):.5f}  MC {r['Salg'].mean():.5f} +- {se(r['Salg']):.5f} | S_alg(B) MC {r['SalgB'].mean():.5f}")
        print(f"   E S_qtr : exact {float(e['Sqtr']):.5f}  MC {r['Sqtr'].mean():.5f} +- {se(r['Sqtr']):.5f} | S_qtr(B) MC {r['SqtrB'].mean():.5f}")
        print(f"   E p_1   : exact {M[0]/D:.5f}  MC {r['p1'].mean():.5f};  Var p_1 exact {M[0]*M[1]/(D**2*(D+1)):.3e} MC {r['p1'].var():.3e}")
        print(f"   E Tr(rho_alg^2): exact {float(e['Z2alg']):.5f} MC {r['Z2alg'].mean():.5f} +- {se(r['Z2alg']):.5f}")
        print(f"   E TrQ(rho~^2) : exact {float(e['Z2qtr']):.5f} MC {r['Z2qtr'].mean():.5f} +- {se(r['Z2qtr']):.5f}")
        # three conventions for the second Renyi entropy
        q = -np.log(r['Z2alg']).mean(); a = -float(mp.log(e['Z2alg']))
        gauss_naive = -np.log(r['Z2g'].mean())          # -log E Z2(g) with E Z1(g)=1 (unnormalized Gaussian)
        print(f"   S2_alg: quenched E[-log Z]={q:.5f} (MC) | annealed -log E Z={a:.5f} (exact) | "
              f"unnormalized Gaussian -log E Z2(g)={gauss_naive:.5f} (MC); exact shift log(1+1/D)={np.log(1+1/D):.5f}")
