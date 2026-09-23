"""Constants of the source-detector resolvent bound (Theorem B):
  k2(w,L)  = || w G0(w; x_o, .) W_L^{1/2} ||_{L^2}  = (1/2)|phi_R(x_o)| (int W_L |psi_+|^2)^{1/2}
  k3(w,L)  = || W_L^{1/2} R0(w) W_L^{1/2} ||_{L^2 -> L^2}      (Birman-Schwinger operator)
  n(w,L)   = |A(w)|^2 int W_L |psi_+|^2   ->  N_L^2 = (1/pi) int_0^inf n dw  = int int W_L psi0^2 dx dtau
"""
import numpy as np, sys, time, json
from multiprocessing import Pool
from freqdomain import Freq, psi_plus_series, integrate, rho_far, psi_minus_series
from rwcore import W, rho_of_x, x_of_rho

LS = [13.0, 15.0, 20.0, 30.0, 40.0, 60.0, 100.0, 200.0]
NG = 64
gx, gw = np.polynomial.legendre.leggauss(NG)


def work(w):
    F = Freq(w)
    out = []
    for L in LS:
        xs = L + gx                      # Gauss nodes on [L-1, L+1]
        rs = rho_of_x(xs)
        rr = float(rho_of_x(np.array([L + 1.0]))[0]) + 1e-9
        if w * rr >= 25:
            q, qp, _ = psi_plus_series(w, np.array([rr]))
            a1, a1p = q, qp
        else:
            rf = rho_far(w, L)
            q, qp, _ = psi_plus_series(w, np.array([rf]))
            a1, a1p = integrate(w, [q[0], qp[0]], rf, [rr])
        order = np.argsort(-rs)
        u, _ = integrate(w, [a1[0], a1p[0]], rr, rs[order])
        pp = np.empty(NG, complex); pp[order] = u
        pt = np.conj(pp)
        # psi_- on the barrier: integrate outward from the horizon series (stable; avoids the
        # catastrophic cancellation of A_out psi_+ + A_in psi~_+ at low frequency)
        p0, pp0 = psi_minus_series(w, np.array([3.2]))
        asc = np.argsort(rs)
        vm, _ = integrate(w, [p0[0], pp0[0]], 3.2, rs[asc])
        pm = np.empty(NG, complex); pm[asc] = vm
        Wx = W(xs - L)
        J = np.sum(gw * Wx * np.abs(pp) ** 2)
        k2 = 0.5 * abs(F.phiR) * np.sqrt(J)
        # kernel G0(y,y') = -psi_-(y<) psi_+(y>)/W0
        Y1, Y2 = np.meshgrid(xs, xs, indexing='ij')
        PM1, PP2 = np.meshgrid(pm, pp, indexing='ij')
        PP1, PM2 = np.meshgrid(pp, pm, indexing='ij')
        G = np.where(Y1 <= Y2, -PM1 * PP2, -PP1 * PM2) / F.W0
        s = np.sqrt(gw * Wx)
        K = s[:, None] * G * s[None, :]
        k3 = np.linalg.norm(K, 2)
        n = abs(F.Aamp) ** 2 * J
        out.append((k2, k3, n, J))
    return w, out


if __name__ == '__main__':
    ws = np.unique(np.concatenate([np.geomspace(1e-4, 0.05, 150), np.arange(0.05, 12.0, 0.002)]))
    t0 = time.time()
    with Pool(12) as p:
        res = p.map(work, ws, chunksize=8)
    arr = np.array([[list(o) for o in r[1]] for r in res])   # (nw, nL, 4)
    np.savez('kconst.npz', w=ws, L=np.array(LS), k2=arr[:, :, 0], k3=arr[:, :, 1], n=arr[:, :, 2],
             J=arr[:, :, 3])
    print('done', time.time() - t0)
