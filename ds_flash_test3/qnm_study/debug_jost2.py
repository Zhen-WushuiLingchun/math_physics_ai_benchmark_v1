import numpy as np
from qnm_lib import *
from scipy.interpolate import CubicSpline

Vspl, _ = potential_spline(eps=0.0)
wlit = 0.3736716844 - 0.0889623158j

print("W at literature omega vs integration window:")
for xm_abs in [30, 45, 60, 80, 100, 150]:
    W,_,_,_,_ = jost(wlit, Vspl, xmin=-xm_abs, xmax=xm_abs, xm=0.0, rtol=1e-13, atol=1e-16)
    print(f"  xmax={xm_abs:4d}   W={W.real:+.6e}{W.imag:+.6e}i   |W|={abs(W):.6e}")

print()
w0 = find_qnm(wlit, Vspl, xmin=-60, xmax=60, rtol=1e-13, atol=1e-16)
print(f"Newton from literature: omega0 = {w0.real:.12f} {w0.imag:+.12f}i")
print(f"literature            : omega0 = {wlit.real:.12f} {wlit.imag:+.12f}i")
W0,_,_,_,_ = jost(w0, Vspl, xmin=-60, xmax=60, rtol=1e-13, atol=1e-16)
print(f"|W(omega0)| = {abs(W0):.3e}")
