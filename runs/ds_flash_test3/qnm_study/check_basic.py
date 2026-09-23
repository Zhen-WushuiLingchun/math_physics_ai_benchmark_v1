import numpy as np
from qnm_lib import *

print("=== geometry checks ===")
for x in [-70.0, -10.0, -2.0, 0.0, 2.0, 10.0, 12.0, 15.0, 30.0]:
    r = rho_of_x_scalar(x)
    # verify mapping
    xr = r + 2*np.log(r/2 - 1)
    print(f"x={x:8.3f}  rho={r:14.10f}  check x(rho)={xr:14.10f}  V0={6*(r-2)*(r-1)/r**4:.10f}")

print()
print("=== bump ===")
for y in [0.0, 0.25, 0.5, 0.75, 0.9, 0.99]:
    print(f"y={y:5.2f} W={W_bump(y):.10f} W'={dW_dy(y):.10f}")

print()
print("=== initial data normalization ===")
C, I = initial_C()
print(f"C = {C:.12f}   (integral I = {I:.12f})")
# energy check
y = np.arange(-1, 1+1e-5, 1e-5)
w = W_bump(y); dw = dW_dy(y); v0 = V0_of_x(y)
E0 = 0.5*C**2*np.trapezoid(2*dw**2 + v0*w**2, y)
print(f"E0 = {E0:.12f}")

print()
print("=== QNM l=2 (unperturbed) ===")
Vspl, xg = potential_spline(eps=0.0)
w0_guess = 0.3736716844 - 0.0889623158j
w0 = find_qnm(w0_guess, Vspl)
print(f"omega0 = {w0.real:.12f} {w0.imag:+.12f}i")
print(f"literature guess: {w0_guess}")
print(f"|Im| = {-w0.imag:.12f},  c_* = 1/(2|Im|) = {1/(2*(-w0.imag)):.6f}")
