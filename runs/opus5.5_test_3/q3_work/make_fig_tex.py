"""Vector figures for the XeLaTeX report."""
import numpy as np, json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9.5, 'legend.fontsize': 7.5,
                     'mathtext.fontset': 'cm', 'axes.linewidth': 0.7})
OUT = '../latex/figs/'
w0 = 0.3736716844180418 - 0.0889623156889357j
g0 = -w0.imag

# ---------------- Fig 1: summary (same content as Q3_figure.png) ----------------
d = np.load('td_batch.npz')
fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))
k = 'e0.001_L60_h256'; tau = d[k+'_tau']; h0 = d[k+'_h0']; dh = d[k+'_dh']
a = ax[0, 0]
a.semilogy(tau, np.abs(h0), lw=0.6, color='C0', label=r'$|h_0|$')
a.semilogy(tau, np.abs(dh), lw=0.6, color='C3', label=r'$|h_\varepsilon-h_0|$')
a.axvline(2*60-13, color='k', ls=':', lw=0.7); a.text(2*60-13+3, 2e-1, r'$\tau_*=2L-13$', fontsize=8)
a.set_ylim(1e-14, 3); a.set_xlim(0, 320); a.set_xlabel(r'$\tau$')
a.legend(loc='upper right'); a.set_title(r'(a) waveforms, $\varepsilon=10^{-3}$, $L=60$')
a = ax[0, 1]
H0 = 1.006738
for (e, L) in [(1e-2, 13.), (1e-1, 20.), (1e-3, 20.), (1e-3, 60.), (1e-6, 80.)]:
    k = 'e%g_L%g_h256' % (e, L); tau = d[k+'_tau']; dh = d[k+'_dh']
    c = np.concatenate([[0], np.cumsum(0.5*(dh[1:]**2+dh[:-1]**2)*np.diff(tau))])
    a.plot(tau, np.sqrt(c/H0)/e, lw=1, label=r'$\varepsilon=10^{%d}$, $L=%g$' % (round(np.log10(e)), L))
a.set_xlabel(r'$T$'); a.set_ylabel(r'$\mathcal{E}_{\varepsilon,L}(T)/\varepsilon$'); a.legend(loc='lower right')
a.set_title('(b) normalized waveform error (Q3a)')
a = ax[1, 0]
r = json.load(open('roots_tdcases.json'))
pts = np.array([[t[1], t[2]] for t in r['0.001_60']['roots']])
a.plot(pts[:, 0], pts[:, 1], 'o', mfc='none', color='C3', ms=5, label=r'$\varepsilon=10^{-3},L=60$')
pts = np.array([[t[1], t[2]] for t in r['0.01_40']['roots']])
a.plot(pts[:, 0], pts[:, 1], 's', mfc='none', color='C1', ms=5, label=r'$\varepsilon=10^{-2},L=40$')
b = json.load(open('roots_branches.json'))['branches_1e-06_80']
pts = np.array([t['exact'] for t in b])
a.plot(pts[:, 0], pts[:, 1], '^', mfc='none', color='C2', ms=5, label=r'$\varepsilon=10^{-6},L=80$')
a.plot([w0.real], [w0.imag], 'k*', ms=10, label=r'unperturbed $\omega_0$')
a.set_xlabel(r'Re $\omega$'); a.set_ylabel(r'Im $\omega$'); a.legend(loc='lower center', ncol=2)
a.set_title('(c) perturbed resonances (same systems)'); a.set_ylim(-0.135, -0.035)
a = ax[1, 1]
ll = json.load(open('largeL_300_1000_3000_10000.json')); sm = json.load(open('largeL_13_60_100.json'))
fd = json.load(open('fd_results.json'))
Ls = [13, 60, 100, 300, 1000, 3000, 10000]
for e in ['0.3', '0.1', '0.01', '0.001']:
    ys = [(ll.get(str(float(L))) or sm.get(str(float(L))))[e]['E_over_eps'] for L in Ls]
    a.semilogx(Ls, ys, 'o-', ms=3.5, lw=1, label=r'$\varepsilon=%s$' % e)
a.axhline(fd['kappa_inf'], color='k', ls='--', lw=0.8, label=r'$\kappa_\infty$ (Born, $L\to\infty$)')
a.set_xlabel(r'$L$'); a.set_ylabel(r'$\mathcal{E}_{\varepsilon,L}(\infty)/\varepsilon$'); a.legend(loc='upper right')
a.set_title(r'(d) uniformity in $L$ (exact, frequency domain)')
plt.tight_layout(); plt.savefig(OUT + 'fig_summary.pdf'); plt.close()

# ---------------- Fig 2: double limit L = c log(1/eps) ----------------
rb = json.load(open('roots_branches.json'))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
cols = {'3': 'C0', '4': 'C1', '5': 'C2', '6.5': 'C3', '8': 'C4'}
for c in ['3', '4', '5', '6.5', '8']:
    rows = rb['c=' + c]
    x = [np.log10(1/r_['eps']) for r_ in rows]
    y = [r_['root'][1] for r_ in rows]
    ax[0].plot(x, y, 'o-', ms=3.5, lw=1, color=cols[c], label=r'$c=%s$' % c)
    if float(c) > 1/(2*g0):
        ax[0].axhline(-1/(2*float(c)), color=cols[c], ls=':', lw=0.8)
    dist = [abs(complex(*r_['root']) - w0) for r_ in rows]
    ax[1].loglog([r_['eps'] for r_ in rows], dist, 'o-', ms=3.5, lw=1, color=cols[c], label=r'$c=%s$' % c)
ax[0].axhline(w0.imag, color='k', ls='--', lw=0.8, label=r'Im $\omega_0$')
ax[0].set_xlabel(r'$\log_{10}(1/\varepsilon)$'); ax[0].set_ylabel(r'Im $\omega$ (continuation of $\omega_0$)')
ax[0].set_title(r'(a) $L=c\log(1/\varepsilon)$; dotted: $-1/(2c)$'); ax[0].set_ylim(-0.1, -0.025); ax[0].legend(ncol=3, loc='upper center')
ee = np.array([1e-10, 1e-6])
for c, cc in [('3', 3.0), ('4', 4.0)]:
    rows = rb['c=' + c]
    y10 = [abs(complex(*r_['root']) - w0) for r_ in rows if r_['eps'] == 1e-10][0]
    ax[1].loglog(ee, y10 * (ee / 1e-10) ** (1 - 2 * g0 * cc), ':', color='k', lw=1.0)
ax[1].invert_xaxis(); ax[1].set_xlabel(r'$\varepsilon$'); ax[1].set_ylabel(r'$|\omega-\omega_0|$')
ax[1].set_title(r'(b) shift; dotted: slope $1-2\gamma_0c$'); ax[1].legend(ncol=2, loc='lower left')
plt.tight_layout(); plt.savefig(OUT + 'fig_double_limit.pdf'); plt.close()

# ---------------- Fig 3: first order vs Lambert-W ----------------
rn = json.load(open('roots_near_w0.json'))
Lam = np.array([r_['Lambda'] for r_ in rn])
e1 = np.array([abs(complex(*r_['d1']) - (complex(*r_['root']) - w0)) for r_ in rn])
elw = np.array([abs(complex(*r_['dLW']) - (complex(*r_['root']) - w0)) for r_ in rn])
d1 = np.array([abs(complex(*r_['d1'])) for r_ in rn])
o = np.argsort(Lam)
fig, a = plt.subplots(figsize=(3.6, 2.8))
a.loglog(Lam[o], (e1/d1)[o], 'o-', ms=4, lw=1, label='first order')
a.loglog(Lam[o], (elw/d1)[o], 's-', ms=4, lw=1, label='Lambert-$W$ ($k=0$)')
a.axvline(1/np.e, color='k', ls=':', lw=0.8); a.text(1/np.e*1.08, 1.6, r'$1/e$', fontsize=8)
a.set_xlabel(r'$\Lambda=2L|\delta\omega^{(1)}|$'); a.set_ylabel(r'error $/\,|\delta\omega^{(1)}|$')
a.set_title('relative error vs exact resonance'); a.set_ylim(2e-3, 3); a.legend(loc='lower left')
plt.tight_layout(); plt.savefig(OUT + 'fig_lambertw.pdf'); plt.close()
print('ok')
