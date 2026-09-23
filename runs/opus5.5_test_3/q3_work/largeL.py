"""Uniformity in L: exact E_{eps,L}(inf) for large L (up to 1e4) via the exact echo formula
    delta hat h = -i w A phiR beta/(1 - beta R),
with the L-independent factors interpolated from fd_scan (w>=0.01) or computed directly (w<0.01),
and beta(w; eps, L) computed exactly on a grid with spacing pi/(8L).
Also the phase-averaged value  |beta|^2/(1-|beta R|^2)  (the L->inf limit form)."""
import numpy as np, sys, time, json
from multiprocessing import Pool
from scipy.interpolate import CubicSpline
from freqdomain import Freq, psi_plus_series, integrate, rho_far
from rwcore import W, rho_of_x, x_of_rho, V0_of_rho
from scipy.integrate import solve_ivp

EPSS = [0.3, 0.1, 0.01, 0.001]
WMAX = 4.0


def beta_at(w, L, epss, nbar=201):
    rl = float(rho_of_x(np.array([L - 1.0]))[0]) - 1e-9
    rr = float(rho_of_x(np.array([L + 1.0]))[0]) + 1e-9
    if w * rr >= 25:
        q, qp, _ = psi_plus_series(w, np.array([rr]))
        a1, a1p = q, qp
    else:
        rf = rho_far(w, L)
        q, qp, _ = psi_plus_series(w, np.array([rf]))
        a1, a1p = integrate(w, [q[0], qp[0]], rf, [rr])
    rgrid = np.linspace(rr, rl, nbar)
    u, up = integrate(w, [a1[0], a1p[0]], rr, rgrid)
    P, Pp = u[-1], up[-1]
    Pt, Ptp = np.conj(P), np.conj(Pp)
    wr = lambda a_, ap, b_, bp: a_ * bp - ap * b_
    xg = x_of_rho(rgrid)
    bb = -np.trapezoid(W(xg - L) * u ** 2 / (1.0 - 2.0 / rgrid), rgrid) / (2j * w)
    betas = []
    for e in epss:
        v, vp = integrate(w, [a1[0], a1p[0]], rr, [rl], eps=e, L=L)
        a = wr(Pt, Ptp, v[-1], vp[-1]) / (2j * w)
        b = wr(P, Pp, v[-1], vp[-1]) / (-2j * w)
        betas.append(b / a)
    return np.array(betas), bb


def job(args):
    w, L = args
    betas, bb = beta_at(w, L, EPSS)
    if w < 0.01:
        F = Freq(w)
        return w, betas, bb, F.Aamp * F.phiR, F.R
    return w, betas, bb, None, None


if __name__ == '__main__':
    Ls = [float(v) for v in sys.argv[1].split(',')]
    d = np.load('fd_scan.npz')
    H0 = json.load(open('fd_results.json'))['H0']
    ws0 = d['w']
    m = ws0 >= 0.009
    APr = CubicSpline(ws0[m], (d['Aamp'] * d['phiR'])[m].real)
    APi = CubicSpline(ws0[m], (d['Aamp'] * d['phiR'])[m].imag)
    Rr = CubicSpline(ws0[m], d['R'][m].real)
    Ri = CubicSpline(ws0[m], d['R'][m].imag)
    results = {}
    for L in Ls:
        dw = np.pi / (8 * L)
        ws = np.arange(dw / 2, WMAX, dw)
        t0 = time.time()
        with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 12) as p:
            res = p.map(job, [(w, L) for w in ws], chunksize=32)
        betas = np.array([r[1] for r in res])            # (nw, neps)
        bb = np.array([r[2] for r in res])
        AP = np.where(ws < 0.01, 0, APr(ws) + 1j * APi(ws)).astype(complex)
        R = np.where(ws < 0.01, 0, Rr(ws) + 1j * Ri(ws)).astype(complex)
        for k, r in enumerate(res):
            if r[3] is not None:
                AP[k] = r[3]; R[k] = r[4]
        wAP2 = np.abs(ws * AP) ** 2
        row = dict(L=L, born=float(np.sqrt(np.trapezoid(wAP2 * np.abs(bb) ** 2, ws) / np.pi / H0)))
        for j, e in enumerate(EPSS):
            be = betas[:, j]
            ex = np.trapezoid(wAP2 * np.abs(be / (1 - be * R)) ** 2, ws) / np.pi / H0
            av = np.trapezoid(wAP2 * np.abs(be) ** 2 / (1 - np.abs(be * R) ** 2), ws) / np.pi / H0
            low = np.trapezoid((wAP2 * np.abs(be / (1 - be * R)) ** 2)[ws < 0.1], ws[ws < 0.1]) / np.pi / H0
            row[str(e)] = dict(E_over_eps=float(np.sqrt(ex) / e), Eavg_over_eps=float(np.sqrt(av) / e),
                               lowfreq_frac=float(low / ex), max_abs_beta=float(np.max(np.abs(be))),
                               max_enh=float(np.max(np.abs(1 / (1 - be * R)))))
        results[str(L)] = row
        print(json.dumps(row), 'secs', round(time.time() - t0, 1), flush=True)
        np.savez('largeL_%g.npz' % L, w=ws, betas=betas, bb=bb, AP=AP, R=R)
    json.dump(results, open('largeL_%s.json' % sys.argv[1].replace(',', '_'), 'w'), indent=1)
