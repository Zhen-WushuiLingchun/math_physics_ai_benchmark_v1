import mpmath as mp
from qnm_mp import *

wlit = mp.mpc('0.3736716844', '-0.0889623158')

# 1) free case
Vs_free = ([mp.mpf(0)]*(_N+1), [mp.mpf(0)]*_N)
W = rk4_full(wlit, Vs_free, record=False)[0]
print("free W =", mp.nstr(W, 15), "   2iw =", mp.nstr(2j*wlit, 15))

# 2) Wronskian x-independence for RW
Vs = V_arrays()
for xm in [-5.0, -2.0, 0.0, 2.0, 5.0]:
    W = rk4_full(wlit, Vs, xm_idx=idx_of(xm), record=False)[0]
    print(f"RW W(xm={xm:5.1f}) = {mp.nstr(W, 12)}")

# 3) Poschl-Teller with the new integrator
def VPT(x):
    c = mp.cosh(x)
    return 2/c**2
Vi = [VPT(XMIN + j*DX) for j in range(_N+1)]
Vh = [VPT(XMIN + (j+mp.mpf('0.5'))*DX) for j in range(_N)]
VsPT = (Vi, Vh)
wPT, it = find_qnm_mp(mp.mpc('1.3228756555', '-0.5'), VsPT)
print("PT found:", mp.nstr(wPT, 20), " exact sqrt(7)/2 - i/2 =", mp.nstr(mp.sqrt(7)/2 - mp.mpf('0.5')*1j, 20))
