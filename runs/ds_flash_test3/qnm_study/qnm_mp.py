"""
High-precision (mpmath) QNM shooting for the RW problem + smooth bump.

Conventions: e^{-i omega tau};  psi_L ~ e^{-i w x} (x->-inf), psi_R ~ e^{+i w x} (x->+inf).
QNM condition: W(w) = psi_L psi_R' - psi_L' psi_R = 0.

RK4 with exact potential evaluations at both integer and half-integer fine-grid
points (no interpolation error, true 4th order).
"""
import mpmath as mp
import pickle
import sympy as _sp

mp.mp.dps = 40

XMIN = mp.mpf(-62)
XMAX = mp.mpf(62)
DX = mp.mpf('0.025')          # integration step = fine grid spacing

# ---- asymptotic coefficients for the right Jost solution
# psi_R = e^{i w x} F(x),  F = 1 + sum_k G_k(log x) x^{-k}
_Lx_s, _w_s = _sp.symbols('Lx omega')
_Nas = 0
_Gfun = {}
_dGfun = {}

def load_asymp(l=2, N=8):
    global _Nas, _Gfun, _dGfun
    with open(f'asymp_coeffs_l{l}_N{N}.pkl', 'rb') as f:
        d = pickle.load(f)
    _Nas = d['N']
    _Gfun = {}
    _dGfun = {}
    for k in range(1, _Nas + 1):
        _Gfun[k] = _sp.lambdify((_Lx_s, _w_s), _sp.expand(d['G'][k]), 'mpmath')
        _dGfun[k] = _sp.lambdify((_Lx_s, _w_s), _sp.expand(_sp.diff(d['G'][k], _Lx_s)),
                                 'mpmath')

load_asymp(2, 8)

def asymp_F(w, x):
    """F(x) and dF/dx with F from the asymptotic series, all in mpmath."""
    x = mp.mpf(x)
    Lx = mp.log(x)
    w = mp.mpc(w)
    xn = 1 / x
    F = mp.mpf(1)
    Fp = mp.mpf(0)
    for k in range(1, _Nas + 1):
        Gk = _Gfun[k](Lx, w)
        dGk = _dGfun[k](Lx, w)
        F += Gk * xn**k
        Fp += (dGk - k * Gk) * xn**(k + 1)
    return F, Fp

# ---------------------------------------------------------------- geometry
def rho_mp(x):
    x = mp.mpf(x)
    if x > 30:
        base = (x - 2 * mp.log(x / 2)) / 2 - 1
        u = mp.log(base)
    else:
        u = (x - 2) / 2
    for _ in range(80):
        eu = mp.exp(u)
        g = 2 * eu + 2 * u - (x - 2)
        gp = 2 * eu + 2
        step = g / gp
        u -= step
        if abs(step) < mp.mpf(10)**(-mp.mp.dps + 4) * max(1, abs(u)):
            break
    return 2 + 2 * mp.exp(u)

def V0_mp(x):
    r = rho_mp(x)
    return 6 * (r - 2) * (r - 1) / r**4

def Wb_mp(y):
    y = mp.mpf(y)
    if abs(y) >= 1:
        return mp.mpf(0)
    return mp.exp(1 - 1 / (1 - y * y))

# ------------------------------------------------- potential on two grids
_N = int((XMAX - XMIN) / DX)            # number of fine steps
def V_arrays(eps=mp.mpf(0), L=mp.mpf(0)):
    """V on fine grid (integer points) and half-shifted grid."""
    Vi = []
    Vh = []
    for j in range(_N + 1):
        x = XMIN + j * DX
        v = V0_mp(x)
        if eps != 0:
            v += eps * Wb_mp(x - L)
        Vi.append(v)
    for j in range(_N):
        x = XMIN + (j + mp.mpf('0.5')) * DX
        v = V0_mp(x)
        if eps != 0:
            v += eps * Wb_mp(x - L)
        Vh.append(v)
    return Vi, Vh

# ----------------------------------------------------------------- RK4
def rk4_full(w, Vs, xm_idx=None, record=True):
    """Integrate both Jost solutions; return W at xm_idx, plus (if record)
    full recorded solutions on the fine grid."""
    if xm_idx is None:
        xm_idx = idx_of(0)
    w = mp.mpc(w)
    w2 = w * w
    Vi, Vh = Vs
    # ---- left solution: from index 0 to _N (full domain)
    psi = mp.exp(-1j * w * XMIN)
    phi = -1j * w * psi
    recL = {0: (psi, phi)} if record else None
    yL_at_xm = None
    for j in range(0, _N):
        v0, vh, v1 = Vi[j], Vh[j], Vi[j + 1]
        k1p, k1q = phi, (v0 - w2) * psi
        s2p = psi + DX / 2 * k1p; s2q = phi + DX / 2 * k1q
        k2p, k2q = s2q, (vh - w2) * s2p
        s3p = psi + DX / 2 * k2p; s3q = phi + DX / 2 * k2q
        k3p, k3q = s3q, (vh - w2) * s3p
        s4p = psi + DX * k3p; s4q = phi + DX * k3q
        k4p, k4q = s4q, (v1 - w2) * s4p
        psi = psi + DX / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        phi = phi + DX / 6 * (k1q + 2 * k2q + 2 * k3q + k4q)
        if record:
            recL[j + 1] = (psi, phi)
        if j + 1 == xm_idx:
            yL_at_xm = (psi, phi)
    yL = (psi, phi)
    # ---- right solution: from index _N down to 0 (backward)
    F0, Fp0 = asymp_F(w, XMAX)
    e0 = mp.exp(1j * w * XMAX)
    psi = e0 * F0
    phi = 1j * w * psi + e0 * Fp0
    recR = {_N: (psi, phi)} if record else None
    yR_at_xm = None
    for j in range(_N, 0, -1):
        v0, vh, v1 = Vi[j], Vh[j - 1], Vi[j - 1]
        h = -DX
        k1p, k1q = phi, (v0 - w2) * psi
        s2p = psi + h / 2 * k1p; s2q = phi + h / 2 * k1q
        k2p, k2q = s2q, (vh - w2) * s2p
        s3p = psi + h / 2 * k2p; s3q = phi + h / 2 * k2q
        k3p, k3q = s3q, (vh - w2) * s3p
        s4p = psi + h * k3p; s4q = phi + h * k3q
        k4p, k4q = s4q, (v1 - w2) * s4p
        psi = psi + h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        phi = phi + h / 6 * (k1q + 2 * k2q + 2 * k3q + k4q)
        if record:
            recR[j - 1] = (psi, phi)
        if j - 1 == xm_idx:
            yR_at_xm = (psi, phi)
    yR = (psi, phi)
    W = yL_at_xm[0] * yR_at_xm[1] - yL_at_xm[1] * yR_at_xm[0]
    return W, yL_at_xm, yR_at_xm, recL, recR

def find_qnm_mp(w0, Vs, xm_idx=None, tol=mp.mpf('1e-28'), itmax=40, h=None):
    if xm_idx is None:
        xm_idx = idx_of(0)
    if h is None:
        h = DX * mp.mpf('0.01')     # ~2.5e-4
    w = mp.mpc(w0)
    for it in range(itmax):
        W0 = rk4_full(w, Vs, xm_idx, record=False)[0]
        dW = (rk4_full(w + h, Vs, xm_idx, record=False)[0]
              - rk4_full(w - h, Vs, xm_idx, record=False)[0]) / (2 * h)
        step = W0 / dW
        w = w - step
        if abs(step) < tol * max(1, abs(w)):
            break
    return w, it

def idx_of(x):
    return int(round((mp.mpf(x) - XMIN) / DX))
