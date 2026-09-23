# Symbolic checks (sympy): static solutions, asymptotic recurrence, horizon series ODE.
import sympy as sp
r, w = sp.symbols('rho omega')
f = 1 - 2/r
Vf = 6/r**2 - 6/r**3          # V0 / f
# 1) static solution rho^3 : d/drho ( f dpsi/drho ) = (V0/f) psi
psi = r**3
print("static rho^3 residual:", sp.simplify(sp.diff(f*sp.diff(psi, r), r) - Vf*psi))
# second static solution v0 = rho^3 * int_rho^oo drho'/(rho'^5 (rho'-2))
s = sp.symbols('s', positive=True)
Iint = sp.integrate(1/(s**5*(s-2)), (s, r, sp.oo))
v0 = sp.simplify(r**3*Iint)
print("v0 =", v0)
print("v0 residual:", sp.simplify(sp.diff(f*sp.diff(v0, r), r) - Vf*v0))
print("v0 large-rho series:", sp.series(v0, r, sp.oo, 5))
# Wronskian in x: W_x(u,v) = f (u v' - u' v)
Wx = sp.simplify(f*(psi*sp.diff(v0, r) - sp.diff(psi, r)*v0))
print("W_x(rho^3, v0) =", Wx)
G00 = sp.simplify(psi*v0/(-Wx))
print("G0(0;y,y) = u0 v0 / |W| =", G00, "  large-rho:", sp.series(G00, r, sp.oo, 2))
# 2) asymptotic recurrence: psi_+ = e^{i w x} sum a_n rho^{-n}
N = 8
a = sp.symbols('a0:%d' % (N+2))
g = sum(a[n]*r**(-n) for n in range(N+2))
expr = sp.diff(f*sp.diff(g, r), r) + 2*sp.I*w*sp.diff(g, r) - Vf*g
expr = sp.expand(expr*r**(N+3))
# check recurrence 2 i w (k+1) a_{k+1} = [k(k+1)-6] a_k + (8-2k^2) a_{k-1}
ok = True
for k in range(0, N):
    coeff = sp.expand(expr).coeff(r, N+3-(k+2))
    am1 = a[k-1] if k >= 1 else 0
    rec = -2*sp.I*w*(k+1)*a[k+1] + (k*(k+1)-6)*a[k] + (8-2*k**2)*am1
    if sp.simplify(coeff - rec) != 0:
        ok = False; print("mismatch at k", k, sp.simplify(coeff-rec))
print("asymptotic recurrence verified:", ok)
