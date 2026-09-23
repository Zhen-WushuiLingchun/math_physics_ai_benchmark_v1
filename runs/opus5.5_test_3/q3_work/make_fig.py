import numpy as np, json, glob
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
d = np.load('td_batch.npz')
fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
# (a) waveforms
k = 'e0.001_L60_h256'; tau = d[k+'_tau']; h0 = d[k+'_h0']; dh = d[k+'_dh']
a = ax[0, 0]
a.semilogy(tau, np.abs(h0), lw=0.8, label=r'$|h_0|$')
a.semilogy(tau, np.abs(dh), lw=0.8, label=r'$|h_\varepsilon-h_0|$, $\varepsilon=10^{-3}$, $L=60$')
a.axvline(2*60-13, color='k', ls=':', lw=0.8); a.text(2*60-13+2, 1e-1, r'$2L-13$')
a.set_ylim(1e-14, 3); a.set_xlim(0, 320); a.set_xlabel(r'$\tau$'); a.legend(fontsize=8); a.set_title('(a) waveform and exact difference')
# (b) E(T)/eps
a = ax[0, 1]
H0 = 1.006738
for (e, L) in [(1e-1, 20.), (1e-3, 20.), (1e-3, 60.), (1e-6, 80.), (1e-2, 13.)]:
    k = 'e%g_L%g_h256' % (e, L); tau = d[k+'_tau']; dh = d[k+'_dh']
    c = np.concatenate([[0], np.cumsum(0.5*(dh[1:]**2+dh[:-1]**2)*np.diff(tau))])
    a.plot(tau, np.sqrt(c/H0)/e, label=r'$\varepsilon=%g, L=%g$' % (e, L))
a.set_xlabel('T'); a.set_ylabel(r'$\mathcal{E}_{\varepsilon,L}(T)/\varepsilon$'); a.legend(fontsize=8); a.set_title('(b) normalized waveform error (Q3a)')
# (c) resonances
a = ax[1, 0]
w0 = 0.3736716844180418-0.0889623156889357j
r = json.load(open('roots_tdcases.json'))
for key, mk in [('0.001_60', 'o'), ('0.01_40', 's')]:
    pts = np.array([[t[1], t[2]] for t in r[key]['roots']])
    a.plot(pts[:, 0], pts[:, 1], mk, mfc='none', label=r'resonances $\varepsilon,L$=%s' % key.replace('_', ', '))
b = json.load(open('roots_branches.json'))['branches_1e-06_80']
pts = np.array([t['exact'] for t in b]); a.plot(pts[:, 0], pts[:, 1], '^', mfc='none', label=r'resonances $\varepsilon,L=10^{-6},80$')
a.plot([w0.real], [w0.imag], 'k*', ms=12, label=r'unperturbed $\omega_0$')
for c in [8.69, 5.79]:
    a.axhline(-1/(2*c), ls=':', lw=0.7, color='gray')
a.set_xlabel(r'Re $\omega$'); a.set_ylabel(r'Im $\omega$'); a.legend(fontsize=7); a.set_title('(c) pole migration (same systems as (a),(b))')
# (d) E(inf)/eps vs L
a = ax[1, 1]
fd = json.load(open('fd_results.json'))
ll = json.load(open('largeL_300_1000_3000_10000.json'))
eps = ['0.3', '0.1', '0.01', '0.001']
Ls = [13, 60, 100, 300, 1000, 3000, 10000]
small = {}
for f in ['largeL_13_60_100.json']:
    try:
        small = json.load(open(f))
    except Exception:
        pass
for e in eps:
    ys = []
    for L in Ls:
        key = str(float(L))
        row = ll.get(key) or small.get(key)
        ys.append(row[e]['E_over_eps'] if row else np.nan)
    a.semilogx(Ls, ys, 'o-', label=r'$\varepsilon=%s$' % e)
a.axhline(fd['kappa_inf'], color='k', ls='--', lw=0.8, label=r'$\kappa_\infty$ (Born, $L\to\infty$)')
a.set_xlabel('L'); a.set_ylabel(r'$\mathcal{E}_{\varepsilon,L}(\infty)/\varepsilon$'); a.legend(fontsize=8); a.set_title('(d) uniformity in L (frequency domain, exact)')
plt.tight_layout(); plt.savefig('../Q3_figure.png', dpi=130)
print('saved')
