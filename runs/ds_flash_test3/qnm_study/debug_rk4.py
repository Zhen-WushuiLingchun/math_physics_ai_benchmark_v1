import mpmath as mp
mp.mp.dps = 50
from qnm_mp import rk4, index_of, XMIN, XMAX, DX
import qnm_mp

w = mp.mpc('0.3736716844', '-0.0889623158')
Vs = [mp.mpf(0)] * (int((XMAX - XMIN)/DX) + 1)
print("len Vs", len(Vs), "i_match", index_of(0.0))

# left
y0 = (mp.exp(-1j*w*XMIN), -1j*w*mp.exp(-1j*w*XMIN))
yL = rk4(y0, Vs, 0, index_of(0.0), w)
print("psi_L(0) =", mp.nstr(yL[0], 15), "  exact 1")
print("phi_L(0) =", mp.nstr(yL[1], 15), "  exact", mp.nstr(-1j*w, 15))

# right
y0r = (mp.exp(1j*w*XMAX), 1j*w*mp.exp(1j*w*XMAX))
yR = rk4(y0r, Vs, index_of(XMAX), index_of(0.0), w)
print("psi_R(0) =", mp.nstr(yR[0], 15), "  exact 1")
print("phi_R(0) =", mp.nstr(yR[1], 15), "  exact", mp.nstr(1j*w, 15))

W = yL[0]*yR[1] - yL[1]*yR[0]
print("W =", mp.nstr(W, 15), "  exact 2iw =", mp.nstr(2j*w, 15))
print()
print("check step size consistency: DX =", DX, " XMIN =", XMIN)
print("i(XMAX) =", index_of(XMAX))
