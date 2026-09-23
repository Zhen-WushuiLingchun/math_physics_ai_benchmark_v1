import mpmath as mp
import qnm_mp as Q
from qnm_mp import *

wlit = mp.mpc('0.3736716844', '-0.0889623158')
Vs = V_arrays()

W, yL, yR, recL, recR = rk4_full(wlit, Vs, record=True)
print("psi_L(0) =", mp.nstr(yL[0], 12), " |.|=", mp.nstr(abs(yL[0]),6))
print("psi_R(0) =", mp.nstr(yR[0], 12), " |.|=", mp.nstr(abs(yR[0]),6))
print("W =", mp.nstr(W, 12))

# test: run right solution from an inner start point with the SAME init form
# (if solution were a pure e^{i w x}-normalized Jost solution, starting at
#  different XMAX should give proportional results at x=0)
for xmax_try in [30.0, 45.0, 62.0]:
    # build an integration from index of xmax_try to index of 0
    xm_idx = idx_of(xmax_try)
    psi = mp.exp(1j*wlit*xmax_try)
    phi = 1j*wlit*psi
    w = wlit; w2 = w*w
    Vi, Vh = Vs
    for j in range(xm_idx, idx_of(0), -1):
        v0, vh, v1 = Vi[j], Vh[j-1], Vi[j-1]
        h = -DX
        k1p, k1q = phi, (v0-w2)*psi
        s2p = psi+h/2*k1p; s2q = phi+h/2*k1q
        k2p, k2q = s2q, (vh-w2)*s2p
        s3p = psi+h/2*k2p; s3q = phi+h/2*k2q
        k3p, k3q = s3q, (vh-w2)*s3p
        s4p = psi+h*k3p; s4q = phi+h*k3q
        k4p, k4q = s4q, (v1-w2)*s4p
        psi = psi + h/6*(k1p+2*k2p+2*k3p+k4p)
        phi = phi + h/6*(k1q+2*k2q+2*k3q+k4q)
    # ratio to the value obtained from the full run
    print(f"XMAX={xmax_try:5.1f}: psi_R(0) = {mp.nstr(psi,10)}   ratio to full run = {mp.nstr(psi/yR[0],8)}")
