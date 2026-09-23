import numpy as np, json
from scipy.optimize import least_squares

d = np.load('td_batch.npz')
fd = json.load(open('fd_results.json'))
w0 = 0.3736716844180418 - 0.0889623156889357j
w1 = 0.34671099687916346 - 0.2739148752912348j
g0 = -w0.imag
CASES = [(1e-1, 20.), (1e-2, 20.), (1e-3, 20.), (1e-3, 30.), (1e-4, 50.), (1e-6, 80.), (1e-2, 40.),
         (1e-2, 13.), (1e-3, 60.)]


def get(e, L, n):
    k = 'e%g_L%g_h%d' % (e, L, n)
    return d[k + '_tau'], d[k + '_h0'], d[k + '_dh'], d[k + '_h1'], float(d[k + '_onset'])


def cum(f, tau):
    c = np.concatenate([[0], np.cumsum(0.5 * (f[1:] ** 2 + f[:-1] ** 2) * np.diff(tau))])
    return c


def fit_single(t, y):
    """fit y ~ 2 Re[A e^{-i w (t-t1)}] with complex A, w (normalized least squares, multistart)."""
    t1 = t[0]; sc = np.max(np.abs(y)); tt = t - t1; yy = y / sc
    def res(p):
        A = p[0] + 1j * p[1]; w = p[2] + 1j * p[3]
        return (2 * (A * np.exp(-1j * w * tt)).real - yy)
    best = None
    for wr in (0.2, 0.3, 0.37, 0.45):
        for wi in (-0.03, -0.09, -0.2):
            for ph in (0, 1, 2, 3, 4, 5):
                s = least_squares(res, [0.5 * np.cos(ph), 0.5 * np.sin(ph), wr, wi], xtol=1e-15, ftol=1e-15, gtol=1e-15)
                if best is None or s.cost < best.cost:
                    best = s
    p = best.x
    A = (p[0] + 1j * p[1]) * sc
    w = p[2] + 1j * p[3]
    return w, A * np.exp(1j * w * t1), np.sqrt(2 * best.cost / np.sum(yy ** 2))


out = {}
for (e, L) in CASES:
    tau, h0, dh, h1, on = get(e, L, 256)
    tau2, h0b, dhb, h1b, on2 = get(e, L, 128)
    ts = 2 * L - 13
    H0_256 = np.trapezoid(h0 ** 2, tau); H0_128 = np.trapezoid(h0b ** 2, tau2)
    H0 = (4 * H0_256 - H0_128) / 3
    Ecum = np.sqrt(cum(dh, tau) / H0)
    Ecum_b = np.sqrt(cum(dhb, tau2) / H0)
    Bcum = np.sqrt(cum(h1, tau) / H0)
    Tlist = [ts + 1, ts + 2, ts + 5, ts + 10, ts + 20, ts + 50, tau[-1]]
    row = dict(eps=e, L=L, onset_discrete=on, onset_theory=ts, onset_128=on2, H0=H0)
    row['E'] = {('%.1f' % T): float(np.interp(T, tau, Ecum)) for T in Tlist}
    row['E_over_eps_Tmax'] = float(Ecum[-1] / e)
    row['E_over_eps_Tmax_h128'] = float(Ecum_b[-1] / e)
    row['Born_Tmax'] = float(Bcum[-1])
    # mismatch
    Ms = {}
    for T in Tlist + [9.5, 10.0, ts - 1]:
        m = tau <= T
        he = (h0 + dh)[m]; hh = h0[m]; tt = tau[m]
        n0 = np.sqrt(np.trapezoid(hh ** 2, tt)); ne = np.sqrt(np.trapezoid(he ** 2, tt))
        if n0 == 0 or ne == 0:
            Ms['%.1f' % T] = None; continue
        ip = np.trapezoid(he * hh, tt)
        M = 1 - abs(ip) / (n0 * ne)
        # Born prediction m2 = ||P_perp h1||^2 / (2 ||h0||^2)
        h1m = h1[m]
        par = np.trapezoid(h1m * hh, tt) / n0 ** 2
        perp = h1m - par * hh
        m2 = np.trapezoid(perp ** 2, tt) / (2 * n0 ** 2)
        Ms['%.1f' % T] = dict(M=float(M), M_over_eps2=float(M / e ** 2), m2_born=float(m2),
                              bound=float(np.trapezoid(dh[m] ** 2, tt) / n0 ** 2))
    row['M'] = Ms
    # relative pointwise error in the echo window
    win = (tau > ts) & (tau < ts + 25)
    row['max_dh_echo'] = float(np.max(abs(dh[win])))
    row['max_h0_echo'] = float(np.max(abs(h0[win])))
    row['rel_err_echo'] = row['max_dh_echo'] / row['max_h0_echo']
    row['eps_e2gL'] = e * np.exp(2 * g0 * L)
    # QNM fit of h_eps just before the echo, and of the echo tail
    if ts - 25 > 30:
        m = (tau > ts - 25) & (tau < ts - 1)
        wf, Af, rel = fit_single(tau[m], (h0 + dh)[m])
        row['fit_before_echo'] = dict(w=[wf.real, wf.imag], A=[Af.real, Af.imag], relres=rel)
    m = (tau > ts + 15) & (tau < ts + 45)
    wf, Af, rel = fit_single(tau[m], dh[m])
    row['fit_echo_ringdown'] = dict(w=[wf.real, wf.imag], relres=rel)
    wf2, Af2, rel2 = fit_single(tau[m], (h0 + dh)[m])
    row['fit_total_after_echo'] = dict(w=[wf2.real, wf2.imag], relres=rel2)
    out['e%g_L%g' % (e, L)] = row
    print('eps=%g L=%g onset=%.4f (theory %g; h128 %.4f) E(Tmax)/eps=%.6f [h128 %.6f] FD=%s Born=%.6f | rel_err_echo=%.3g  eps e^{2gL}=%.3g' % (
        e, L, on, ts, on2, Ecum[-1] / e, Ecum_b[-1] / e,
        fd['rows'][[r['L'] for r in fd['rows']].index(L)]['full'].get(str(e)) if L in [r['L'] for r in fd['rows']] else None,
        Bcum[-1], row['rel_err_echo'], row['eps_e2gL']))
    print('    E(T):', row['E'])
    print('    M(T):', {k: (None if v is None else (round(v['M_over_eps2'], 5), round(v['m2_born'], 5))) for k, v in Ms.items()})
    if 'fit_before_echo' in row:
        print('    fit before echo:', np.round(row['fit_before_echo']['w'], 6), 'res', '%.1e' % row['fit_before_echo']['relres'], ' | dh after echo:', np.round(row['fit_echo_ringdown']['w'], 5), 'res %.1e' % row['fit_echo_ringdown']['relres'], ' | h_eps after echo:', np.round(row['fit_total_after_echo']['w'], 5), 'res %.1e' % row['fit_total_after_echo']['relres'])
json.dump(out, open('td_results.json', 'w'), indent=1, default=float)
# QNM residue check on h0: fit amplitudes with fixed w0, w1 on [25, 60]
tau, h0, dh, h1, on = get(1e-3, 60., 256)
m = (tau > 25) & (tau < 70)
t = tau[m]
Mx = np.column_stack([np.cos(w0.real * t) * np.exp(w0.imag * t), np.sin(w0.real * t) * np.exp(w0.imag * t),
                      np.cos(w1.real * t) * np.exp(w1.imag * t), np.sin(w1.real * t) * np.exp(w1.imag * t)])
c, *_ = np.linalg.lstsq(Mx, h0[m], rcond=None)
# 2Re[A e^{-i w t}] = 2 e^{Im w t}(Re A cos - ... ) : 2Re[(a+ib)(cos - i sin)] = 2(a cos + b sin)
A0 = (c[0] + 1j * c[1]) / 2; A1 = (c[2] + 1j * c[3]) / 2
print('fitted amplitude of w0 mode A0 =', A0, ' (predicted -i Res =', -1j * (0.04692739601703409 + 0.014382262417433542j), ')')
print('fitted amplitude of w1 mode A1 =', A1, ' (predicted -i Res =', -1j * (0.3636462808507978 + 0.09329920215952356j), ')')
print('residual rel:', np.linalg.norm(Mx @ c - h0[m]) / np.linalg.norm(h0[m]))
