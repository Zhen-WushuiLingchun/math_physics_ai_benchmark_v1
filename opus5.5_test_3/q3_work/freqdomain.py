"""Frequency-domain (Jost solution) machinery for RW l=2 + bump, real or complex omega.

Conventions: time dependence e^{-i w tau};  Laplace transform  hat f(w) = int_0^inf e^{i w tau} f dtau.
psi_-(w,x) ~ e^{-i w x} (x -> -inf),  psi_+(w,x) ~ e^{+i w x} (x -> +inf).
Derivatives d/dx are denoted by 'p' (e.g. psip = d psi/dx).
"""
import numpy as np
from scipy.integrate import solve_ivp
from rwcore import W, dW, C_NORM, x_of_rho, rho_of_x, V0_of_rho, XO

RHO_O = float(rho_of_x(np.array([XO]))[0])


# ----------------------------------------------------------------- horizon series
def psi_minus_series(w, rho, nmax=4000, tol=1e-18):
    """psi_- = e^{-i w x} g(rho), g = sum c_n (rho-2)^n, c_0 = 1.  Valid for |rho-2|<2.
    Returns psi, dpsi/dx."""
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    z = rho - 2.0
    c = [1.0 + 0j]
    # n = 0:  c1 * 1*(4 - 16 i w) + c0*(-6) = 0
    g = np.ones_like(z, dtype=complex)
    gp = np.zeros_like(z, dtype=complex)
    zp = np.ones_like(z)
    cm2, cm1, cn = 0j, 0j, 1.0 + 0j
    zmax = z.max()
    for n in range(0, nmax):
        num = (cn * (4 * n * (n - 1) + 2 * n - 24j * w * n - 6)
               + cm1 * ((n - 1) * (n - 2) - 12j * w * (n - 1) - 6)
               + cm2 * (-2j * w * (n - 2)))
        cnp1 = -num / ((n + 1) * (4 * n + 4 - 16j * w))
        # accumulate term n+1
        gp += (n + 1) * cnp1 * zp            # derivative uses z^n
        zp = zp * z
        g += cnp1 * zp
        cm2, cm1, cn = cm1, cn, cnp1
        if n > 10 and abs(cnp1) * zmax ** (n + 1) * (n + 2) < tol * np.max(np.abs(g)) and \
                abs(cm1) * zmax ** n * (n + 1) < tol * np.max(np.abs(g)):
            break
    x = x_of_rho(rho)
    f = 1.0 - 2.0 / rho
    e = np.exp(-1j * w * x)
    psi = e * g
    psip = e * (-1j * w * g + f * gp)
    return psi, psip


# ----------------------------------------------------------------- infinity series
def asym_coeffs(w, kmax=400):
    a = np.zeros(kmax + 1, dtype=complex)
    a[0] = 1.0
    a[1] = -6.0 / (2j * w)
    for k in range(1, kmax):
        a[k + 1] = ((k * (k + 1) - 6) * a[k] + (8 - 2 * k * k) * a[k - 1]) / (2j * w * (k + 1))
    return a


def psi_plus_series(w, rho, kmax=2000, tol=1e-18):
    """psi_+ = e^{i w x} sum a_k rho^{-k}, optimally truncated (scaled recurrence for t_k = a_k rho^-k).
    Works for complex w and complex rho. Returns psi, dpsi/dx (x-derivative), est. rel. error."""
    rho = np.atleast_1d(np.asarray(rho, dtype=complex))
    out_psi = np.empty(rho.shape, dtype=complex)
    out_psip = np.empty(rho.shape, dtype=complex)
    err = np.empty(rho.shape)
    for i, r in enumerate(rho):
        # NOTE: for s=2,l=2 one has a_3 = 0 exactly, so terms are not monotone at the start;
        # compute terms, then truncate at the smallest term (k>=4) = optimal truncation.
        ts = [1.0 + 0j]
        tm1, t = 0j, 1.0 + 0j
        mn = np.inf
        for k in range(0, kmax):
            tn = ((k * (k + 1) - 6) * t + (8 - 2 * k * k) * tm1 / r) / (2j * w * (k + 1) * r)
            ts.append(tn)
            tm1, t = t, tn
            if k >= 4:
                mn = min(mn, abs(tn))
                if abs(tn) < tol or abs(tn) > 1e8 * mn:
                    break
        ts = np.array(ts)
        mags = np.abs(ts)
        kk = np.arange(len(ts))
        if mags[-1] < tol:
            kstop = len(ts)
            e_est = mags[-1]
        else:
            kstop = int(np.argmin(np.where(kk >= 4, mags, np.inf)))
            e_est = mags[kstop]
        g = ts[:kstop].sum()
        gp = (-kk[:kstop] * ts[:kstop]).sum() / r
        x = r + 2.0 * np.log(r / 2.0 - 1.0)
        f = 1.0 - 2.0 / r
        e = np.exp(1j * w * x)
        out_psi[i] = e * g
        out_psip[i] = e * (1j * w * g + f * gp)
        err[i] = e_est / abs(g)
    return out_psi, out_psip, err


# ----------------------------------------------------------------- ODE in rho
def _rhs_factory(w, eps, L):
    w2 = w * w

    def rhs(r, y):
        f = 1.0 - 2.0 / r
        V = V0_of_rho(r)
        if eps != 0.0:
            xx = r + 2.0 * np.log(r / 2.0 - 1.0)
            V = V + eps * W(np.array([xx - L]))[0]
        return [y[1] / f, (V - w2) * y[0] / f]
    return rhs


def integrate(w, y0, r0, r_eval, eps=0.0, L=0.0, rtol=1e-12, atol=1e-14):
    """Integrate (psi, dpsi/dx) from rho=r0 to the points r_eval (monotone, same direction)."""
    r_eval = np.atleast_1d(np.asarray(r_eval, float))
    sol = solve_ivp(_rhs_factory(w, eps, L), (r0, r_eval[-1]), np.asarray(y0, complex),
                    method='DOP853', t_eval=r_eval, rtol=rtol, atol=atol * max(1.0, abs(y0[0])))
    if sol.status != 0:
        raise RuntimeError(sol.message)
    return sol.y[0], sol.y[1]


def rho_far(w, L=None, target=45.0):
    r = max(target / abs(w), 60.0)
    if L is not None:
        r = max(r, float(rho_of_x(np.array([L + 1.5]))[0]) + 5.0)
    return r


# ----------------------------------------------------------------- source integral
_GL_N = 400
_gl_x, _gl_w = np.polynomial.legendre.leggauss(_GL_N)


def source_integral(w):
    """I(w) = int_{-1}^{1} psi_-(w,y) S(w,y) dy,  S = -C W' - i w C W."""
    y = _gl_x
    rho = rho_of_x(y)
    psi, _ = psi_minus_series(w, rho)
    S = -C_NORM * dW(y) - 1j * w * C_NORM * W(y)
    return np.sum(_gl_w * psi * S)


# ----------------------------------------------------------------- main per-frequency object
class Freq:
    """All unperturbed Jost data at one real frequency w>0, plus methods for barrier at L."""

    def __init__(self, w, rho_match=None):
        self.w = w
        self.I = source_integral(w)
        rm = RHO_O if rho_match is None else rho_match
        # psi_- at x_o: series valid up to rho<4; integrate from rho=3.2
        r0 = 3.2
        p0, pp0 = psi_minus_series(w, np.array([r0]))
        pm, pmp = integrate(w, [p0[0], pp0[0]], r0, [rm])
        self.pm, self.pmp = pm[0], pmp[0]
        # psi_+ at x_o: from far series, inward
        rf = rho_far(w)
        q, qp, _ = psi_plus_series(w, np.array([rf]))
        pp, ppp = integrate(w, [q[0], qp[0]], rf, [rm])
        self.pp, self.ppp = pp[0], ppp[0]
        # conj solution (real w): psi~_+ = conj(psi_+)
        self.pt, self.ptp = np.conj(self.pp), np.conj(self.ppp)
        wr = lambda a, ap, b, bp: a * bp - ap * b
        self.W0 = wr(self.pm, self.pmp, self.pp, self.ppp)            # W(psi_-, psi_+) = 2 i w A_in
        self.Ain = self.W0 / (2j * w)
        self.Aout = wr(self.pm, self.pmp, self.pt, self.ptp) / (-2j * w)
        self.R = self.Aout / self.Ain
        self.phiR = self.pm / self.Ain                                # psi_-(x_o)/A_in
        self.Aamp = -self.I / (2j * w * self.Ain)                    # outgoing amplitude
        self.psi0_xo = -self.pp * self.I / self.W0                   # hat psi_0(w, x_o)
        self.h0 = -1j * w * self.psi0_xo

    def barrier(self, eps, L, nbar=801, series_ok=25.0):
        """eps: scalar or list. Returns dict of arrays over eps:
        a, b (psi_+^eps = a psi_+ + b psi~_+ left of barrier), beta, Born b1 (per unit eps),
        dh = hat h_eps - hat h_0 (exact), h1 = Born hat h_1 (per unit eps)."""
        w = self.w
        epss = np.atleast_1d(np.asarray(eps, float))
        rl = float(rho_of_x(np.array([L - 1.0]))[0]) - 1e-9
        rr = float(rho_of_x(np.array([L + 1.0]))[0]) + 1e-9
        if w * rr >= series_ok:
            q, qp, err = psi_plus_series(w, np.array([rr]))
            a1, a1p = q, qp
        else:
            rf = rho_far(w, L)
            q, qp, _ = psi_plus_series(w, np.array([rf]))
            a1, a1p = integrate(w, [q[0], qp[0]], rf, [rr])
        rgrid = np.linspace(rr, rl, nbar)
        u, up = integrate(w, [a1[0], a1p[0]], rr, rgrid)            # unperturbed psi_+ on barrier
        P, Pp = u[-1], up[-1]
        Pt, Ptp = np.conj(P), np.conj(Pp)
        wr = lambda a_, ap, b_, bp: a_ * bp - ap * b_
        xg = x_of_rho(rgrid)
        integrand = W(xg - L) * u ** 2 / (1.0 - 2.0 / rgrid)      # dx = drho / f
        b_born = -np.trapezoid(integrand, rgrid) / (2j * w)          # rgrid descending
        out = dict(a=[], b=[], beta=[], dh=[])
        for e in epss:
            v, vp = integrate(w, [a1[0], a1p[0]], rr, [rl], eps=e, L=L)
            a = wr(Pt, Ptp, v[-1], vp[-1]) / (2j * w)
            b = wr(P, Pp, v[-1], vp[-1]) / (-2j * w)
            beta = b / a
            dpsi = self.Aamp * self.phiR * beta / (1.0 - beta * self.R)
            out['a'].append(a); out['b'].append(b); out['beta'].append(beta)
            out['dh'].append(-1j * w * dpsi)
        out = {k: np.array(v_) for k, v_ in out.items()}
        out['b_born'] = b_born
        out['h1'] = -1j * w * self.Aamp * self.phiR * b_born
        out['absPsiPlusBar'] = np.max(np.abs(u))
        return out
