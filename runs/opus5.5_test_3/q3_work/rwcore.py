"""Core definitions for Problem 3 (Regge-Wheeler l=2 odd parity + small bump at x=L).

Units: M = 1.  x = tortoise coordinate, rho = r/M.
"""
import numpy as np
from scipy.special import lambertw
from scipy.integrate import quad
import numba as nb

XO = 10.0  # observer


def rho_of_x(x):
    """Invert x = rho + 2 log(rho/2 - 1):  rho = 2 (1 + W0(exp(x/2 - 1)))."""
    x = np.asarray(x, dtype=float)
    t = x / 2.0 - 1.0
    out = np.empty_like(t)
    small = t < 600
    out[small] = lambertw(np.exp(t[small])).real
    tb = t[~small]
    if tb.size:
        w = tb - np.log(tb)
        for _ in range(8):  # Newton on w + log w - t = 0
            w = w - (w + np.log(w) - tb) / (1.0 + 1.0 / w)
        out[~small] = w
    return 2.0 * (1.0 + out)


def x_of_rho(rho):
    return rho + 2.0 * np.log(rho / 2.0 - 1.0)


def V0_of_rho(rho):
    return (1.0 - 2.0 / rho) * (6.0 / rho**2 - 6.0 / rho**3)


def V0(x):
    return V0_of_rho(rho_of_x(x))


def W(y):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    m = np.abs(y) < 1
    out[m] = np.exp(1.0 - 1.0 / (1.0 - y[m] ** 2))
    return out


def dW(y):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    m = np.abs(y) < 1
    ym = y[m]
    out[m] = np.exp(1.0 - 1.0 / (1.0 - ym**2)) * (-2.0 * ym / (1.0 - ym**2) ** 2)
    return out


def d2W(y):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    m = np.abs(y) < 1
    ym = y[m]
    q = 1.0 - ym**2
    g = -2.0 * ym / q**2
    gp = -2.0 / q**2 - 8.0 * ym**2 / q**3
    out[m] = np.exp(1.0 - 1.0 / q) * (g**2 + gp)
    return out


def _q(fun, a=-1.0, b=1.0):
    return quad(fun, a, b, epsabs=1e-15, epsrel=1e-13, limit=400)[0]


def norm_constant():
    """C such that E0 = (C^2/2) int (2 W'^2 + V0 W^2) dx = 1."""
    I = _q(lambda y: 2.0 * dW(y) ** 2 + V0(y) * W(y) ** 2)
    return np.sqrt(2.0 / I)


C_NORM = norm_constant()


def W_norms():
    L1 = _q(lambda y: W(y))
    L2 = np.sqrt(_q(lambda y: W(y) ** 2))
    Linf = 1.0
    return L1, L2, Linf


def E1_energy():
    """Energy of d_tau psi_0 at tau=0 (conserved), for the Sobolev/energy bound."""
    C = C_NORM
    # d_tau psi(0) = -C W',  d_tau^2 psi(0) = C W'' - V0 C W
    f = lambda y: 0.5 * ((C * d2W(y) - V0(y) * C * W(y)) ** 2 + (C * d2W(y)) ** 2
                         + V0(y) * (C * dW(y)) ** 2)
    return _q(f)


# ---------------------------------------------------------------------------
# Time-domain solver A: leapfrog with Courant number 1 (exact numerical light cone)
# ---------------------------------------------------------------------------
@nb.njit(cache=True)
def _leapfrog(p0, p1, V, S0, S1, src_coef, nsteps, jo, h):
    """Evolve  psi^{n+1} = psi_{j+1}+psi_{j-1}-psi^{n-1} - h^2 V psi^n - h^2 src^n
    Two coupled fields: A (with potential V[0]) and B (with potential V[1], source -src_coef*S*A)
    Here we evolve: psi0 with V0 ; u with V_eps and source  -(src) * psi0.
    p0,p1: arrays shape (2,N) at levels n-1, n.  S: source profile (eps*W_L) acting on field 0.
    Returns recorded field values at jo for both fields at all levels 0..nsteps+1.
    """
    N = p0.shape[1]
    rec = np.zeros((2, nsteps + 2))
    rec[0, 0] = p0[0, jo]; rec[1, 0] = p0[1, jo]
    rec[0, 1] = p1[0, jo]; rec[1, 1] = p1[1, jo]
    pm = p0.copy(); pc = p1.copy(); pn = np.zeros_like(p0)
    h2 = h * h
    for n in range(1, nsteps + 1):
        for j in range(1, N - 1):
            # Gundlach-Price-Pullin diamond form (stable at Courant number 1):
            # psi_N = psi_E + psi_W - psi_S - (h^2/2) V_C (psi_E + psi_W)
            aE = pc[0, j + 1] + pc[0, j - 1]
            bE = pc[1, j + 1] + pc[1, j - 1]
            pn[0, j] = aE - pm[0, j] - 0.5 * h2 * V[0, j] * aE
            pn[1, j] = bE - pm[1, j] - 0.5 * h2 * V[1, j] * bE - 0.5 * h2 * src_coef * S0[j] * aE
        pn[0, 0] = 0.0; pn[0, N - 1] = 0.0; pn[1, 0] = 0.0; pn[1, N - 1] = 0.0
        tmp = pm; pm = pc; pc = pn; pn = tmp
        rec[0, n + 1] = pc[0, jo]; rec[1, n + 1] = pc[1, jo]
    return rec


def td_leapfrog(eps, L, T, h, born=False, xo=XO):
    """Return tau grid, h0(tau), dh(tau) = h_eps - h0 (or Born u1 if born=True, i.e. d/deps at eps=0)
    computed with the CFL=1 leapfrog scheme; h = d psi/d tau by 2nd-order central differences."""
    # domain large enough that the numerical domain of dependence of (T, xo) is inside
    xl = xo - T - 4 * h
    xr = xo + T + 4 * h
    N = int(np.ceil((xr - xl) / h)) + 1
    jo = int(round((xo - xl) / h))
    x = xl + (np.arange(N) - 0) * h
    x = x - x[jo] + xo  # make xo exactly a grid point
    V0x = V0(x)
    WL = W(x - L)
    C = C_NORM
    # initial data: psi(0)=C W(x), psi(h) via Taylor to O(h^5)
    p0 = C * W(x)
    d1 = -C * dW(x)
    # psi_tt = psi'' - V psi ; psi_ttt = psi1'' - V psi1
    tt = C * d2W(x) - V0x * p0
    # psi(h) = C W(x-h) - h^2/2 V psi0 - h^3/6 V psi1 + h^4/24 [ (psi_tt)'' - V psi_tt - psi0'''' ]
    # the free part C W(x-h) contains all pure-derivative terms exactly.
    # h^4 term: psi_tttt = (psi_tt)'' - V psi_tt ; free part of it = C W''''.
    # (psi_tt)'' = C W'''' - (V psi0)'' ;  so correction = -(V psi0)'' - V psi_tt
    Vp0 = V0x * p0
    d2Vp0 = np.gradient(np.gradient(Vp0, h), h)
    p1_0 = C * W(x - h) - 0.5 * h**2 * Vp0 - h**3 / 6.0 * V0x * d1 + h**4 / 24.0 * (-d2Vp0 - V0x * tt)
    nsteps = int(np.ceil(T / h)) + 1
    Vs = np.vstack([V0x, V0x + (0.0 if born else eps) * WL])
    P0 = np.vstack([p0, np.zeros(N)])
    P1 = np.vstack([p1_0, np.zeros(N)])
    coef = 1.0 if born else eps
    # u has zero data at tau=0 and tau=h up to O(h^?): u(h) = 0 exactly? The source acts only where
    # psi0 overlaps W_L, which happens at tau >= L-2 > 0, so u = 0 identically for early times.
    rec = _leapfrog(P0, P1, Vs, WL, WL, coef, nsteps, jo, h)
    tau = h * np.arange(rec.shape[1])
    # central differences for d/dtau
    hp0 = np.zeros_like(rec[0]); hdu = np.zeros_like(rec[1])
    hp0[1:-1] = (rec[0, 2:] - rec[0, :-2]) / (2 * h)
    hdu[1:-1] = (rec[1, 2:] - rec[1, :-2]) / (2 * h)
    hp0[0] = d1[jo]
    return tau[:-1], hp0[:-1], hdu[:-1], rec


# ---------------------------------------------------------------------------
# Time-domain solver B: method of lines, 6th-order FD in x, RK4 in time
# ---------------------------------------------------------------------------
@nb.njit(cache=True)
def _rhs(psi, pi, V, S, eps_src, psi0, out_psi, out_pi, ih2):
    N = psi.shape[0]
    c0 = -49.0 / 18.0; c1 = 3.0 / 2.0; c2 = -3.0 / 20.0; c3 = 1.0 / 90.0
    for j in range(N):
        out_psi[j] = pi[j]
    for j in range(3, N - 3):
        lap = (c0 * psi[j] + c1 * (psi[j + 1] + psi[j - 1]) + c2 * (psi[j + 2] + psi[j - 2])
               + c3 * (psi[j + 3] + psi[j - 3])) * ih2
        out_pi[j] = lap - V[j] * psi[j] - eps_src * S[j] * psi0[j]
    for j in (0, 1, 2, N - 3, N - 2, N - 1):
        out_pi[j] = 0.0


@nb.njit(cache=True)
def _mol(psiA, piA, psiB, piB, VA, VB, S, eps_src, dt, nsteps, jo, h, rec_every):
    N = psiA.shape[0]
    ih2 = 1.0 / (h * h)
    nrec = nsteps // rec_every + 1
    rec = np.zeros((4, nrec))
    k = 0
    zero = np.zeros(N)
    kA1 = np.zeros((2, N)); kA2 = np.zeros((2, N)); kA3 = np.zeros((2, N)); kA4 = np.zeros((2, N))
    kB1 = np.zeros((2, N)); kB2 = np.zeros((2, N)); kB3 = np.zeros((2, N)); kB4 = np.zeros((2, N))
    tA = np.zeros((2, N)); tB = np.zeros((2, N))
    for n in range(nsteps + 1):
        if n % rec_every == 0:
            rec[0, k] = psiA[jo]; rec[1, k] = piA[jo]; rec[2, k] = psiB[jo]; rec[3, k] = piB[jo]
            k += 1
        if n == nsteps:
            break
        _rhs(psiA, piA, VA, S, 0.0, zero, kA1[0], kA1[1], ih2)
        _rhs(psiB, piB, VB, S, eps_src, psiA, kB1[0], kB1[1], ih2)
        for j in range(N):
            tA[0, j] = psiA[j] + 0.5 * dt * kA1[0, j]; tA[1, j] = piA[j] + 0.5 * dt * kA1[1, j]
            tB[0, j] = psiB[j] + 0.5 * dt * kB1[0, j]; tB[1, j] = piB[j] + 0.5 * dt * kB1[1, j]
        _rhs(tA[0], tA[1], VA, S, 0.0, zero, kA2[0], kA2[1], ih2)
        _rhs(tB[0], tB[1], VB, S, eps_src, tA[0], kB2[0], kB2[1], ih2)
        for j in range(N):
            tA[0, j] = psiA[j] + 0.5 * dt * kA2[0, j]; tA[1, j] = piA[j] + 0.5 * dt * kA2[1, j]
            tB[0, j] = psiB[j] + 0.5 * dt * kB2[0, j]; tB[1, j] = piB[j] + 0.5 * dt * kB2[1, j]
        _rhs(tA[0], tA[1], VA, S, 0.0, zero, kA3[0], kA3[1], ih2)
        _rhs(tB[0], tB[1], VB, S, eps_src, tA[0], kB3[0], kB3[1], ih2)
        for j in range(N):
            tA[0, j] = psiA[j] + dt * kA3[0, j]; tA[1, j] = piA[j] + dt * kA3[1, j]
            tB[0, j] = psiB[j] + dt * kB3[0, j]; tB[1, j] = piB[j] + dt * kB3[1, j]
        _rhs(tA[0], tA[1], VA, S, 0.0, zero, kA4[0], kA4[1], ih2)
        _rhs(tB[0], tB[1], VB, S, eps_src, tA[0], kB4[0], kB4[1], ih2)
        for j in range(N):
            psiA[j] += dt / 6.0 * (kA1[0, j] + 2 * kA2[0, j] + 2 * kA3[0, j] + kA4[0, j])
            piA[j] += dt / 6.0 * (kA1[1, j] + 2 * kA2[1, j] + 2 * kA3[1, j] + kA4[1, j])
            psiB[j] += dt / 6.0 * (kB1[0, j] + 2 * kB2[0, j] + 2 * kB3[0, j] + kB4[0, j])
            piB[j] += dt / 6.0 * (kB1[1, j] + 2 * kB2[1, j] + 2 * kB3[1, j] + kB4[1, j])
    return rec


def td_mol(eps, L, T, h, cfl=0.5, born=False, xo=XO, margin=40.0, rec_dt=None):
    xl = xo - T - margin
    xr = xo + T + margin
    N = int(np.ceil((xr - xl) / h)) + 1
    x = xl + np.arange(N) * h
    jo = int(round((xo - xl) / h))
    x = x - x[jo] + xo
    V0x = V0(x)
    WL = W(x - L)
    C = C_NORM
    psiA = C * W(x); piA = -C * dW(x)
    psiB = np.zeros(N); piB = np.zeros(N)
    dt = cfl * h
    if rec_dt is None:
        rec_dt = dt
    rec_every = max(1, int(round(rec_dt / dt)))
    dt = rec_dt / rec_every
    nsteps = int(np.ceil(T / dt))
    nsteps -= nsteps % rec_every
    VB = V0x + (0.0 if born else eps) * WL
    rec = _mol(psiA, piA, psiB, piB, V0x, VB, WL, 1.0 if born else eps, dt, nsteps, jo, h, rec_every)
    tau = np.arange(rec.shape[1]) * dt * rec_every
    return tau, rec[1], rec[3], rec
