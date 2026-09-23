"""QNMs of RW l=2 (s=2) and resonances of RW + eps W(x-L), complex omega.

- Leaver continued fraction (units 2M=1 internally) for the unperturbed QNMs.
- Jost functions A_in(w), A_out(w) for complex w via complex-rho-path integration from
  optimally truncated asymptotic series (started in the recessive direction).
- Barrier coefficients a(w), b(w) (psi_+^eps = a psi_+ + b psi~_+ left of the bump) for complex w.
- Resonance function F(w) = A_in(w) - beta(w) A_out(w), beta=b/a.
"""
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from rwcore import W, x_of_rho, rho_of_x, V0_of_rho, XO
from freqdomain import psi_minus_series, psi_plus_series, RHO_O


# ------------------------------------------------------------------ Leaver
def leaver_cf(w2M, l=2, s=2, N=400):
    """Continued fraction residual; w2M is omega in units 2M=1."""
    w = mp.mpc(w2M)
    def al(n): return n * n + (2 - 2j * w) * n + 1 - 2j * w
    def be(n): return -(2 * n * n + (2 - 8j * w) * n - 8 * w * w - 4j * w + l * (l + 1) + 1 - s * s)
    def ga(n): return n * n - 4j * w * n - 4 * w * w - s * s
    t = mp.mpc(0)
    for n in range(N, 0, -1):
        t = al(n - 1) * ga(n) / (be(n) - t)
    return (be(0) - t) / be(0)


def leaver_cf_inv(w2M, k, l=2, s=2, N=400):
    """k-th inversion of the CF (better for the k-th overtone)."""
    w = mp.mpc(w2M)
    def al(n): return n * n + (2 - 2j * w) * n + 1 - 2j * w
    def be(n): return -(2 * n * n + (2 - 8j * w) * n - 8 * w * w - 4j * w + l * (l + 1) + 1 - s * s)
    def ga(n): return n * n - 4j * w * n - 4 * w * w - s * s
    # tail
    t = mp.mpc(0)
    for n in range(N, k, -1):
        t = al(n - 1) * ga(n) / (be(n) - t)
    up = t
    # head (finite, downward to 0):  beta_k - alpha_{k-1} gamma_k /(beta_{k-1} - ... )
    hd = mp.mpc(0)
    for n in range(0, k):
        hd = al(n) * ga(n + 1) / (be(n) - hd) if n > 0 else al(0) * ga(1) / be(0)
    return (be(k) - hd - up) / be(k)


def qnm_leaver(guess_M, k=0, N=600):
    mp.mp.dps = 30
    f = (lambda z: leaver_cf(z, N=N)) if k == 0 else (lambda z: leaver_cf_inv(z, k, N=N))
    z = mp.findroot(f, mp.mpc(2 * guess_M))
    return complex(z) / 2.0


# ------------------------------------------------------------------ complex path ODE
def _integrate_path(w, y0, r0, r1, rtol=1e-12, atol=1e-30):
    """Integrate (psi, dpsi/dx) along the straight segment r0 -> r1 (complex rho), unperturbed RW."""
    d = r1 - r0
    w2 = w * w

    def rhs(s, y):
        r = r0 + s * d
        f = 1.0 - 2.0 / r
        V = (1.0 - 2.0 / r) * (6.0 / r**2 - 6.0 / r**3)
        return [d * y[1] / f, d * (V - w2) * y[0] / f]
    sol = solve_ivp(rhs, (0.0, 1.0), np.asarray(y0, complex), method='DOP853', rtol=rtol,
                    atol=atol)
    return sol.y[0, -1], sol.y[1, -1]


def psi_plus_c(w, rho_t, sign=+1, target=40.0):
    """psi_+(sign*w, rho_t) (so sign=-1 gives psi~_+ = psi_+(-w)) with its x-derivative,
    at real or complex rho_t, by starting in the recessive direction of psi_+(sign*w)."""
    ws = sign * w
    theta = np.pi / 2 - np.angle(ws)
    D = max(target / abs(ws), 30.0)
    rs = rho_t + D * np.exp(1j * theta)
    q, qp, err = psi_plus_series(ws, np.array([rs]))
    return _integrate_path(w, [q[0], qp[0]], rs, rho_t), err[0]


def psi_minus_c(w, rho_t, r0=3.2):
    p0, pp0 = psi_minus_series(w, np.array([r0]))
    return _integrate_path(w, [p0[0], pp0[0]], r0, rho_t)


def wr(a, ap, b, bp):
    return a * bp - ap * b


def jost(w, rho_m=RHO_O):
    """A_in, A_out (psi_- = A_out psi_+ + A_in psi~_+), plus psi values at rho_m."""
    pm = psi_minus_c(w, rho_m)
    (pp, ppp), e1 = psi_plus_c(w, rho_m, +1)
    (pt, ptp), e2 = psi_plus_c(w, rho_m, -1)
    Wtp = wr(pt, ptp, pp, ppp)       # should be 2 i w
    Ain = wr(pm[0], pm[1], pp, ppp) / (2j * w)
    Aout = wr(pm[0], pm[1], pt, ptp) / (-2j * w)
    return dict(Ain=Ain, Aout=Aout, check=Wtp / (2j * w) - 1, pm=pm, pp=(pp, ppp), pt=(pt, ptp),
                err=max(e1, e2))


# ------------------------------------------------------------------ barrier at complex w
def barrier_c(w, eps, L, nbar=401):
    rl = float(rho_of_x(np.array([L - 1.0]))[0]) - 1e-9
    rr = float(rho_of_x(np.array([L + 1.0]))[0]) + 1e-9
    (P_r, Pp_r), _ = psi_plus_c(w, rr, +1)
    (P_l, Pp_l), _ = psi_plus_c(w, rl, +1)
    (T_l, Tp_l), _ = psi_plus_c(w, rl, -1)
    w2 = w * w

    def rhs(r, y, e):
        f = 1.0 - 2.0 / r
        V = V0_of_rho(r)
        if e != 0.0:
            xx = r + 2.0 * np.log(r / 2.0 - 1.0)
            V = V + e * W(np.array([xx - L]))[0]
        return [y[1] / f, (V - w2) * y[0] / f]
    rgrid = np.linspace(rr, rl, nbar)
    su = solve_ivp(rhs, (rr, rl), [P_r, Pp_r], method='DOP853', t_eval=rgrid, rtol=1e-12,
                   atol=1e-30, args=(0.0,))
    u = su.y[0]
    xg = x_of_rho(rgrid)
    b_born = -np.trapezoid(W(xg - L) * u**2 / (1 - 2 / rgrid), rgrid) / (2j * w)
    out = dict(b_born=b_born, unpert_check=abs(su.y[0, -1] - P_l) / abs(P_l))
    epss = np.atleast_1d(eps)
    res = []
    for e in epss:
        sv = solve_ivp(rhs, (rr, rl), [P_r, Pp_r], method='DOP853', rtol=1e-12, atol=1e-30,
                       args=(float(e),))
        v, vp = sv.y[0, -1], sv.y[1, -1]
        a = wr(T_l, Tp_l, v, vp) / (2j * w)
        b = wr(P_l, Pp_l, v, vp) / (-2j * w)
        res.append((a, b))
    out['a'] = np.array([r[0] for r in res])
    out['b'] = np.array([r[1] for r in res])
    return out


def F_res(w, eps, L, J=None):
    """Resonance function (A_in - beta A_out); zero <=> resonance of V0 + eps W_L."""
    if J is None:
        J = jost(w)
    B = barrier_c(w, eps, L)
    beta = B['b'][0] / B['a'][0]
    return J['Ain'] - beta * J['Aout']
