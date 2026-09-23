import mpmath as mp
mp.mp.dps = 50
from qnm_mp import rho_mp, V0_mp, V_grid, jost_mp, XMIN, DX
import numpy as np
import sys
sys.path.insert(0, '.')
from qnm_lib import rho_of_x_scalar, V0_of_x, jost as jost_dp, potential_spline

print(" x      rho_mp                    rho_dp                  V_mp            V_dp")
for x in [-76.0, -40, -10, -2, 0, 2, 5, 10, 12, 20, 30, 31, 50, 76]:
    rm = rho_mp(x); rd = rho_of_x_scalar(x)
    vm = V0_mp(x); vd = V0_of_x(x)
    print(f"{x:6.1f} {mp.nstr(rm,18):24s} {rd:24.16f} {float(vm):12.8f} {float(vd):12.8f}")

Vs = V_grid()
print()
print("V_grid[0..3] =", [mp.nstr(v, 8) for v in Vs[:4]])
print("V_grid at x=0 (index", int(-XMIN/DX), ") =", mp.nstr(Vs[int(-XMIN/DX)], 20))

wl = mp.mpc('0.3736716844', '-0.0889623158')
W, yL, yR, _, _ = jost_mp(wl, Vs)
print("W(omega_lit) =", mp.nstr(W, 15))
print("psi_L(0) =", mp.nstr(yL[0], 12), " psi_R(0) =", mp.nstr(yR[0], 12))
