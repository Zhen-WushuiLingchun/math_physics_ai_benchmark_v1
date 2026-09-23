"""
Time-domain 1+1D solver for  psi_tt = psi_xx - V(x) psi,
V = V0(x) + eps*W(x-L);  initial data psi=C W(x), psi_t=-C W'(x).
RK4 in time, 4th-order central differences in space, sponge layers.
"""
import numpy as np
from qnm_lib import rho_of_x, W_bump, dW_dy

def make_grid(L, dx=0.025, xmin=-70.0, pad=70.0):
    xmax = L + pad
    N = int(round((xmax - xmin) / dx))
    x = xmin + dx * np.arange(N + 1)
    rho = rho_of_x(x)
    V0 = 6.0 * (rho - 2.0) * (rho - 1.0) / rho**4
    return x, V0

def sponge_profile(x, width=20.0, sigma0=0.8, xmin=None, xmax=None):
    s = np.zeros_like(x)
    dL = (x - xmin) / width
    dR = (xmax - x) / width
    m = dL < 1.0
    s[m] = sigma0 * (1.0 - dL[m])**2
    m = dR < 1.0
    s[m] = sigma0 * (1.0 - dR[m])**2
    return s

def initial_C(dx=1e-4):
    y = np.arange(-1.0, 1.0 + dx / 2, dx)
    w = W_bump(y); dw = dW_dy(y)
    rho = rho_of_x(y)
    v0 = 6.0 * (rho - 2.0) * (rho - 1.0) / rho**4
    I = np.trapezoid(2.0 * dw**2 + v0 * w**2, y)
    return 1.0 / np.sqrt(0.5 * I)

def run_td(eps, L, dx=0.025, dt=0.0125, Tmax=None, xobs=10.0,
           sponge=True, sponge_width=20.0, sigma0=0.8, store_x=None,
           xmin=-70.0, xmax=None, history_region=None):
    if xmax is None:
        xmax = L + 70.0
    N = int(round((xmax - xmin) / dx))
    x = xmin + dx * np.arange(N + 1)
    rho = rho_of_x(x)
    V0 = 6.0 * (rho - 2.0) * (rho - 1.0) / rho**4
    V = V0 + eps * W_bump(x - L)
    C = initial_C()
    psi = C * W_bump(x)
    psi_t = -C * dW_dy(x)
    # sponge
    sig = sponge_profile(x, width=sponge_width, sigma0=sigma0, xmin=xmin, xmax=xmax) if sponge \
        else np.zeros_like(x)
    if Tmax is None:
        Tmax = 2 * L + 30.0
    nt = int(round(Tmax / dt))
    iobs = int(round((xobs - xmin) / dx))
    h = np.zeros(nt + 1)
    h[0] = psi_t[iobs]
    if history_region is not None:
        ih1 = int(round((history_region[0]-xmin)/dx))
        ih2 = int(round((history_region[1]-xmin)/dx))
        hist = np.zeros((nt+1, ih2-ih1+1))
        hist[0] = psi[ih1:ih2+1]
    # 4th-order Laplacian, one-sided at edges
    def lap(p):
        d2 = np.zeros_like(p)
        d2[2:-2] = (-p[4:] + 16*p[3:-1] - 30*p[2:-2] + 16*p[1:-3] - p[:-4]) / (12*dx*dx)
        # second-order at the 4 edge points
        d2[0] = (p[2] - 2*p[1] + p[0]) / dx**2
        d2[1] = (p[2] - 2*p[1] + p[0]) / dx**2
        d2[-1] = (p[-1] - 2*p[-2] + p[-3]) / dx**2
        d2[-2] = (p[-1] - 2*p[-2] + p[-3]) / dx**2
        return d2
    def rhs(psi, psi_t):
        return psi_t, lap(psi) - V * psi - 2.0 * sig * psi_t
    for n in range(nt):
        k1p, k1v = rhs(psi, psi_t)
        k2p, k2v = rhs(psi + dt/2*k1p, psi_t + dt/2*k1v)
        k3p, k3v = rhs(psi + dt/2*k2p, psi_t + dt/2*k2v)
        k4p, k4v = rhs(psi + dt*k3p, psi_t + dt*k3v)
        psi = psi + dt/6*(k1p + 2*k2p + 2*k3p + k4p)
        psi_t = psi_t + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
        h[n+1] = psi_t[iobs]
        if history_region is not None:
            hist[n+1] = psi[ih1:ih2+1]
    tau = np.arange(nt + 1) * dt
    out = {'tau': tau, 'h': h, 'C': C, 'x': x, 'V': V, 'dx': dx, 'dt': dt}
    if history_region is not None:
        out['hist'] = hist
        out['hist_x'] = x[
            int(round((history_region[0]-xmin)/dx)):int(round((history_region[1]-xmin)/dx))+1]
    return out
    if store_x is not None:
        out['psi'] = psi
        out['psi_t'] = psi_t
    return out

def energy(x, psi, psi_t, V):
    d = np.gradient(psi, x)
    dens = 0.5 * (psi_t**2 + d**2 + V * psi**2)
    return np.trapezoid(dens, x)

def fit_frequency(tau, h, t1, t2, lag=0.5, nmax=400):
    """estimate complex frequency from ratio h(t+lag)/h(t) averaged."""
    i1 = int(round(t1 / (tau[1]-tau[0])))
    i2 = int(round(t2 / (tau[1]-tau[0])))
    il = int(round(lag / (tau[1]-tau[0])))
    ws = []
    for i in range(i1, min(i2, len(h) - il)):
        r = h[i + il] / h[i]
        if abs(h[i]) > 1e-14 and abs(h[i+il]) > 1e-14:
            ws.append(1j * np.log(r) / (il * (tau[1]-tau[0])))
    ws = np.array(ws)
    # robust: median after removing outliers by 2-sigma
    wm = np.median(ws)
    keep = np.abs(ws - wm) < 0.05
    return ws[keep].mean(), ws[keep].std(), len(ws[keep])
