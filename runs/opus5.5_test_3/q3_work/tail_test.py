"""Late-time tail: barrier modifies the tail amplitude by the static factor (8 A_right)^{-2}."""
import numpy as np
from scipy.integrate import solve_ivp
from rwcore import td_leapfrog, W, rho_of_x, x_of_rho, V0_of_rho
from sympy import lambdify, symbols, log
r = symbols('r')
v0expr = (3*r**4*(log(r) - log(r - 2)) - 6*r**3 - 6*r**2 - 8*r - 12)/(96*r)
v0 = lambdify(r, v0expr, 'numpy')
def static_A(eps, L):
    # integrate static solution u = rho^3/8 through the barrier (in rho), compute A = -W_x(u, v0) right of barrier
    rl = float(rho_of_x(np.array([L-1.5]))[0]); rr = float(rho_of_x(np.array([L+1.5]))[0])
    def rhs(rh, y):
        f = 1-2/rh; xx = rh+2*np.log(rh/2-1)
        V = V0_of_rho(rh) + eps*W(np.array([xx-L]))[0]
        return [y[1]/f, V*y[0]/f]
    s = solve_ivp(rhs, (rl, rr), [rl**3/8, (1-2/rl)*3*rl**2/8], method='DOP853', rtol=1e-12, atol=1e-14)
    u, ux = s.y[0, -1], s.y[1, -1]
    f = 1-2/rr; dv0 = (v0(rr+1e-6)-v0(rr-1e-6))/2e-6
    Wx = u*f*dv0 - ux*v0(rr)
    return -Wx
for L, eps in [(20., 0.1), (20., 0.3), (40., 0.1)]:
    A = static_A(eps, L)
    pred = (8*A)**-2
    T = 1600.0
    res = {}
    for h in [1/16, 1/32]:
        tau, h0, dh, rec = td_leapfrog(eps, L, T, h)
        psi0 = rec[0][:-1]; psie = psi0 + rec[1][:-1]
        m = (tau > 800) & (tau < 1600)
        rat = psie[m]/psi0[m]
        res[h] = (tau[m], rat, psi0[m])
    t, rat, p0 = res[1/32]
    idx = [np.argmin(abs(t-tt)) for tt in (800, 1000, 1200, 1400, 1590)]
    print("L=%g eps=%g  8A=%.6f  predicted tail ratio=%.5f   psi_eps/psi_0 at tau=800..1590: %s  (h=1/16: %s)  tau^7 psi0 at 1590: %.4e" % (
        L, eps, 8*A, pred, np.round(rat[idx], 5), np.round(res[1/16][1][[np.argmin(abs(res[1/16][0]-tt)) for tt in (800,1000,1200,1400,1590)]], 5), p0[idx[-1]]*t[idx[-1]]**7))
