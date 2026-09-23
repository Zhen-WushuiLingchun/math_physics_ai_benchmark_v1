"""
First-order Duhamel with the measured unperturbed Green function:
   dh(tau) = -eps int_0^tau K(tau-tau') S(tau') dtau',
   S(tau') = int_{bump} W(z-L) psi_0(tau',z) dz,
   K(s) = d/dtau G(s; x_o, L)  (impulse response at x_o from x=L).
Compare with the measured dh from the perturbed run.  (eps=1e-3, L=20)
"""
import numpy as np
from numpy.fft import rfft, irfft
from td_solve import initial_C, W_bump, dW_dy, rho_of_x

XOBS, XMIN, XMAX, DX, DT = 10.0, -70.0, 102.0, 0.025, 0.0125
EPS, L = 1e-3, 20.0
TMAX = 2*L + 40

def make(eps, L, ic_kind, Tmax):
    from td_solve import sponge_profile
    xmin, xmax = XMIN, XMAX
    N = int(round((xmax-xmin)/DX))
    x = xmin + DX*np.arange(N+1)
    rho = rho_of_x(x)
    V0 = 6.0*(rho-2.0)*(rho-1.0)/rho**4
    V = V0 + eps*W_bump(x-L)
    C = initial_C()
    if ic_kind == 'data':
        psi = C*W_bump(x); psi_t = -C*dW_dy(x)
    else:  # smooth impulse at x=L (Gaussian, width 0.25)
        psi = np.zeros_like(x)
        sig0 = 0.25
        g = np.exp(-0.5*((x-L)/sig0)**2)
        g = g/(np.trapezoid(g, x))
        psi_t = g
    sig = sponge_profile(x, width=20.0, sigma0=0.8, xmin=xmin, xmax=xmax)
    nt = int(round(Tmax/DT))
    iobs = int(round((XOBS-xmin)/DX))
    def lap(p):
        d2 = np.zeros_like(p)
        d2[2:-2] = (-p[4:]+16*p[3:-1]-30*p[2:-2]+16*p[1:-3]-p[:-4])/(12*DX*DX)
        d2[0] = d2[1] = (p[2]-2*p[1]+p[0])/DX**2
        d2[-1] = d2[-2] = (p[-1]-2*p[-2]+p[-3])/DX**2
        return d2
    def rhs(psi, psi_t):
        return psi_t, lap(psi) - V*psi - 2.0*sig*psi_t
    i1 = int(round((L-1.0-xmin)/DX)); i2 = int(round((L+1.0-xmin)/DX))
    h = np.zeros(nt+1)
    S = np.zeros(nt+1)     # bump-weighted field (or its cumsum use)
    hist = np.zeros((nt+1, i2-i1+1))
    h[0] = psi_t[iobs]
    for n in range(nt):
        k1p,k1v = rhs(psi, psi_t)
        k2p,k2v = rhs(psi+DT/2*k1p, psi_t+DT/2*k1v)
        k3p,k3v = rhs(psi+DT/2*k2p, psi_t+DT/2*k2v)
        k4p,k4v = rhs(psi+DT*k3p, psi_t+DT*k3v)
        psi = psi + DT/6*(k1p+2*k2p+2*k3p+k4p)
        psi_t = psi_t + DT/6*(k1v+2*k2v+2*k3v+k4v)
        h[n+1] = psi_t[iobs]
        hist[n+1] = psi[i1:i2+1]
    return dict(x=x, dx=DX, dt=DT, h=h, hist=hist, i1=i1, i2=i2, nt=nt,
                tau=np.arange(nt+1)*DT)

def Wb_arr(xs, L):
    return W_bump(xs-L)

# runs
unp = make(0.0, L, 'data', TMAX)
imp = make(0.0, L, 'impulse', TMAX)
per = make(EPS, L, 'data', TMAX)

x = unp['x']; i1, i2 = unp['i1'], unp['i2']
xb = x[i1:i2+1]
wb = Wb_arr(xb, L)
# S(tau) = int Wb psi_0 dz
S = np.trapezoid(per['hist']*0 + unp['hist']*wb[None, :], xb, axis=1)
# K = h from impulse
K = imp['h']
# dh_pred(tau) = -eps * (K * S)(tau)  (convolution, trapezoid)
nt = unp['nt']
np_ = 2*nt+1
Kf = rfft(K, np_); Sf = rfft(S, np_)
conv = irfft(Kf*Sf, np_)[:nt+1]*DT
dh_pred = -EPS*conv
dh_meas = per['h'] - unp['h']
tau = unp['tau']
m = (tau > 2*L-15) & (tau < 2*L+5)
err = np.max(np.abs(dh_meas[m]-dh_pred[m]))/np.max(np.abs(dh_meas[m]))
nrm = np.sqrt(np.trapezoid(dh_pred[m]**2, tau[m])/np.trapezoid(dh_meas[m]**2, tau[m]))
print(f"eps={EPS} L={L}: first-echo window: max rel err = {err:.4f},  norm ratio = {nrm:.4f}")
print(f"  measured peak = {np.max(np.abs(dh_meas[m])):.5e},  predicted peak = {np.max(np.abs(dh_pred[m])):.5e}")
np.savez('duhamel_check.npz', tau=tau, dh_meas=dh_meas, dh_pred=dh_pred, h0=unp['h'], S=S, K=K)
print("saved duhamel_check.npz")
