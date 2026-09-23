"""
QNM / ringdown study for problem 3.

Units: G=c=1, M=1.  Tortoise coordinate x = rho + 2 ln(rho/2 - 1).
l=2 odd-parity Regge-Wheeler potential:
    V0(x) = (1 - 2/rho) (6/rho^2 - 6/rho^3) = 6 (rho-2)(rho-1)/rho^4
Bump:  W(y) = exp(1 - 1/(1-y^2)) for |y|<1, 0 otherwise.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

# ---------------------------------------------------------------- geometry
def rho_of_x_scalar(xs):
    # Solve 2 e^u + 2 u = xs - 2  with  u = ln(rho/2 - 1), rho = 2 + 2 e^u.
    if xs > 30.0:
        base = (xs - 2.0 * np.log(xs / 2.0)) / 2.0 - 1.0
        u = np.log(base)
    else:
        u = (xs - 2.0) / 2.0
    for _ in range(300):
        eu = np.exp(u)
        g = 2.0 * eu + 2.0 * u - (xs - 2.0)
        gp = 2.0 * eu + 2.0
        step = g / gp
        u -= step
        if abs(step) < 1e-15 * max(1.0, abs(u)):
            break
    return 2.0 + 2.0 * np.exp(u)

rho_of_x = np.vectorize(rho_of_x_scalar, otypes=[float])

def V0_of_x(x):
    rho = rho_of_x(np.asarray(x, dtype=float))
    return 6.0 * (rho - 2.0) * (rho - 1.0) / rho**4

def W_bump(y):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    m = np.abs(y) < 1.0
    out[m] = np.exp(1.0 - 1.0 / (1.0 - y[m]**2))
    return out

def dW_dy(y):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    m = (np.abs(y) < 1.0) & (np.abs(y) > 1e-14)
    out[m] = -2.0 * y[m] / (1.0 - y[m]**2)**2 * W_bump(y[m])
    return out

# ------------------------------------------------------- initial-data norm
def initial_C(dx=2e-4):
    y = np.arange(-1.0, 1.0 + dx / 2, dx)
    w = W_bump(y)
    dw = dW_dy(y)
    v0 = V0_of_x(y)
    integrand = 2.0 * dw**2 + v0 * w**2
    I = np.trapezoid(integrand, y)
    return 1.0 / np.sqrt(0.5 * I), I

# ------------------------------------------------------------- potentials
def potential_spline(L=0.0, eps=0.0, xmin=-150.0, xmax=150.0, dx=0.05):
    xg = np.arange(xmin, xmax + dx / 2, dx)
    V = V0_of_x(xg)
    if eps != 0.0:
        V = V + eps * W_bump(xg - L)
    return CubicSpline(xg, V), xg

# ------------------------------------------------------------- Jost solver
def _rhs(omega, Vspl):
    w2 = omega * omega
    def rhs(x, y):
        V = float(Vspl(x))
        return [y[1], (V - w2) * y[0]]
    return rhs

def jost(omega, Vspl, xmin=-60.0, xmax=60.0, xm=0.0,
         rtol=1e-12, atol=1e-16, tevalL=None, tevalR=None):
    """Return psi_L, psi_L', psi_R, psi_R' at xm (and profiles if teval given)."""
    w = complex(omega)
    rhs = _rhs(w, Vspl)
    y0L = [np.exp(-1j * w * xmin), -1j * w * np.exp(-1j * w * xmin)]
    solL = solve_ivp(rhs, [xmin, xm], y0L, rtol=rtol, atol=atol,
                     method='DOP853', t_eval=tevalL)
    y0R = [np.exp(1j * w * xmax), 1j * w * np.exp(1j * w * xmax)]
    solR = solve_ivp(rhs, [xmax, xm], y0R, rtol=rtol, atol=atol,
                     method='DOP853', t_eval=tevalR)
    yL = solL.y[:, -1]
    yR = solR.y[:, -1]
    W = yL[0] * yR[1] - yL[1] * yR[0]
    return W, yL, yR, solL, solR

def wronskian(omega, Vspl, **kw):
    return jost(omega, Vspl, **kw)[0]

def find_qnm(omega0, Vspl, h=1e-5, itmax=40, tol=1e-13, **kw):
    w = complex(omega0)
    for it in range(itmax):
        W0 = wronskian(w, Vspl, **kw)
        dW = (wronskian(w + h, Vspl, **kw) - wronskian(w - h, Vspl, **kw)) / (2 * h)
        step = W0 / dW
        w -= step
        if abs(step) < tol * max(1.0, abs(w)):
            break
    return w

def bump_overlap(omega, Vspl, L, dy=0.002, **kw):
    """I(L) = int W(y) psi_L(L+y) psi_R(L+y) dy  at frequency omega."""
    ys = np.arange(-1.0, 1.0 + dy / 2, dy)
    tt = L + ys
    W, yL, yR, solL, solR = jost(omega, Vspl, xm=0.0, tevalL=tt, tevalR=tt, **kw)
    # careful: t_eval ordering for backward integration
    pl = solL.y[0]
    pr = solR.y[0]
    return np.trapezoid(W_bump(ys) * pl * pr, ys)

# -------------------------------------------------------- fast fixed-step
def fast_jost(omega, Vspl, xmin=-150.0, xmax=150.0, xm=0.0, dx=0.02):
    """Fixed-step RK4 complex integrator; returns W and the 4 Jost values at xm."""
    w = complex(omega)
    N = int(round((xm - xmin) / dx))
    def rk4_scan(y0, x_start, n, step):
        y = np.array(y0, dtype=complex)
        for _ in range(n):
            x = x_start + step * _
            k1 = _f(x, y, w, Vspl)
            k2 = _f(x + step / 2, y + step / 2 * k1, w, Vspl)
            k3 = _f(x + step / 2, y + step / 2 * k2, w, Vspl)
            k4 = _f(x + step, y + step * k3, w, Vspl)
            y = y + step / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return y
    def _f(x, y, w, Vspl):
        V = float(Vspl(x))
        return np.array([y[1], (V - w * w) * y[0]], dtype=complex)
    yL = rk4_scan([np.exp(-1j * w * xmin), -1j * w * np.exp(-1j * w * xmin)], xmin, N, dx)
    M = int(round((xmax - xm) / dx))
    yR = rk4_scan([np.exp(1j * w * xmax), 1j * w * np.exp(1j * w * xmax)], xmax, M, -dx)
    W = yL[0] * yR[1] - yL[1] * yR[0]
    return W, yL, yR
