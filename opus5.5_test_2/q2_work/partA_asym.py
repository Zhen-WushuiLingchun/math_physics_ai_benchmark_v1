import mpmath as mp
from partA import exact_A
from common import dims
mp.mp.dps = 60
ph = mp.phi; L = mp.log(ph)
pt_inf = ph**2 / (1 + ph**2)
def leading(NR, NB):
    Dl = NR - NB; b = min(NR, NB)
    Sq = b * L - mp.mpf(1) / 2 * ph**(-abs(Dl))
    return Sq - pt_inf * L, Sq
print("Leading: E S_qtr ~ min(NR,NB) ln(phi) - phi^{-|Delta|}/2 ;  E S_alg ~ E S_qtr - p_tau ln(phi),  p_tau=phi^2/(1+phi^2)=%.6f" % float(pt_inf))
print(" NR NB  Delta   E S_alg exact      err_alg*phi^N     E S_qtr exact     err_qtr*phi^N")
for Dl in [0, 1, 2, 3, -1]:
    for NB in [6, 10, 14, 18, 22, 26, 27]:
        NR = NB + Dl
        e = exact_A(NR, NB); la, lq = leading(NR, NB)
        N = NR + NB
        print(f"{NR:3d}{NB:3d}{Dl:6d}   {float(e['Salg']):.12f}  {float((e['Salg']-la)*ph**N):+.6f}   {float(e['Sqtr']):.12f}  {float((e['Sqtr']-lq)*ph**N):+.6f}")
# annealed Renyi-2 in both conventions vs Page-like closed forms
print("\nAnnealed purity: E TrQ(rho~^2) vs phi^-NR + phi^-NB ; E Tr(rho_alg^2) vs (2phi/sqrt5)(phi^-NR+phi^-NB)")
for NR, NB in [(10, 10), (14, 12), (12, 16), (20, 20)]:
    e = exact_A(NR, NB)
    t = ph**(-NR) + ph**(-NB)
    print(NR, NB, float(e['Z2qtr'] / t), float(e['Z2alg'] / (2 * ph / mp.sqrt(5) * t)))
