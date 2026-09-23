import numpy as np
from qnm_lib import *
from scipy.interpolate import CubicSpline

# test 1: free potential V=0 -> W should be exactly 2i*omega
def spline_V(Vfun, xmin=-150, xmax=150, dx=0.05):
    xg = np.arange(xmin, xmax+dx/2, dx)
    return CubicSpline(xg, Vfun(xg))

Vfree = spline_V(lambda x: np.zeros_like(x))
for w in [0.3737-0.0890j, 0.4+0.1j, 0.5-0.2j]:
    W,_,_,_,_ = jost(w, Vfree, rtol=1e-12, atol=1e-15)
    print(f"V=0: omega={w}  W={W:.12g}  2i*omega={2j*w:.12g}")

print()
Vspl, _ = potential_spline(eps=0.0)
wlit = 0.3736716844 - 0.0889623158j
for w in [wlit, 0.3736716844+0j, 0.3730-0.0248j]:
    W,_,_,_,_ = jost(w, Vspl, rtol=1e-12, atol=1e-15)
    print(f"RW: omega={w}  W={W:.10g}  |W|={abs(W):.6g}")

print()
print("scan |W| on grid:")
for re in np.arange(0.25, 0.50, 0.05):
    row = []
    for im in np.arange(-0.15, 0.01, 0.03):
        W,_,_,_,_ = jost(complex(re,im), Vspl, rtol=1e-10, atol=1e-14)
        row.append(f"{abs(W):9.2e}")
    print(f"Re={re:5.3f}  " + " ".join(row))
print("Im values:  " + " ".join([f"{im:9.2f}" for im in np.arange(-0.15,0.01,0.03)]))
