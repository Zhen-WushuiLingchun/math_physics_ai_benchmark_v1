import mpmath as mp
mp.mp.dps = 50
from qnm_mp import rk4, jost_mp, index_of, XMIN, XMAX, DX

w = mp.mpc('0.3736716844', '-0.0889623158')
Vs = [mp.mpf(0)] * (int((XMAX - XMIN)/DX) + 1)

W, yL, yR, _, _ = jost_mp(w, Vs)
print("jost_mp: psi_L(0) =", mp.nstr(yL[0], 15))
print("jost_mp: psi_R(0) =", mp.nstr(yR[0], 15))
print("jost_mp: W =", mp.nstr(W, 15))

y0 = (mp.exp(-1j*w*XMIN), -1j*w*mp.exp(-1j*w*XMIN))
yl2 = rk4(y0, Vs, 0, index_of(0.0), w)
print("rk4    : psi_L(0) =", mp.nstr(yl2[0], 15))

# also check the right side with rk4
y0r = (mp.exp(1j*w*XMAX), 1j*w*mp.exp(1j*w*XMAX))
yr2 = rk4(y0r, Vs, index_of(XMAX), index_of(0.0), w)
print("rk4    : psi_R(0) =", mp.nstr(yr2[0], 15))
W2 = yl2[0]*yr2[1]-yl2[1]*yr2[0]
print("rk4    : W =", mp.nstr(W2, 15))
