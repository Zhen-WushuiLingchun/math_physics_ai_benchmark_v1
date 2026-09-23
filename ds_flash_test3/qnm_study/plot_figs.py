import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': 0.3,
                     'figure.dpi': 150})
import mpmath as mp

# ---------------- Fig 1: waveforms + echo check
d = np.load('td_arrays.npz')
tau0 = d['tau0']; h0 = d['h0']
DT = 0.0125
hL20 = d['eps0.001_L20_h']
t = np.arange(len(hL20))*DT
du = np.load('duhamel2_L20.npz')
fig, ax = plt.subplots(2, 1, figsize=(6.4, 5.2), sharex=False)
ax[0].plot(tau0, h0, lw=0.8, color='C0', label=r'$h_0(\tau)$')
ax[0].plot(t, d['eps0.001_L20_h'], lw=0.8, color='C3', alpha=0.7,
           label=r'$h_{\varepsilon,L}(\tau)$   ($\varepsilon=10^{-3},\,L=20$)')
ax[0].set_xlim(0, 60); ax[0].set_xlabel(r'$\tau$'); ax[0].set_ylabel(r'$h$')
ax[0].legend(fontsize=8)
ax[0].annotate(r'$\tau_*=2L-13$', xy=(27, 0.15), xytext=(33, 0.22),
               arrowprops=dict(arrowstyle='->'))
tu = du['tau']; dhu = du['dh']
m = (tu > 2*20-16) & (tu < 2*20+6)
ax[1].plot(tu[m], dhu[m]/1e-3, lw=1.0, color='C3', label=r'$[h_{\varepsilon,L}-h_0]/\varepsilon$ (PDE)')
ax[1].plot(tu[m], du['dh_lc'][m]/1e-3, lw=1.0, ls='--', color='C0',
           label='light-cone Duhamel (measured $\\psi_0$)')
ax[1].axvline(2*20-13, color='k', lw=0.8, ls=':')
ax[1].set_xlabel(r'$\tau$'); ax[1].set_ylabel(r'$[h_{\varepsilon,L}-h_0]/\varepsilon$')
ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig1_waveforms.pdf'); plt.close()

# ---------------- Fig 2: log|dh| onset
fig, ax = plt.subplots(figsize=(6.4, 3.4))
for key, col in [('eps0.001_L12', 'C0'), ('eps0.001_L20', 'C1'), ('eps0.001_L28', 'C2')]:
    h = d[f'{key}_h']; dh_ = d[f'{key}_dh']
    n = min(len(h), len(tau0))
    tt = tau0[:n]
    L = float(key.split('_L')[1])
    mm = tt > 5
    ax.semilogy(tt[mm], np.abs(dh_[mm]), lw=0.7, color=col, label=f'$L={L:.0f}$')
    ax.axvline(2*L-13, color=col, lw=0.7, ls=':')
ax.axhline(1e-3, color='k', lw=0.5, ls='--')
ax.set_xlim(5, 75); ax.set_ylim(1e-15, 3e-3)
ax.set_xlabel(r'$\tau$'); ax.set_ylabel(r'$|h_{\varepsilon,L}-h_0|$')
ax.legend(fontsize=8, ncol=3)
plt.tight_layout(); plt.savefig('fig2_causality.pdf'); plt.close()

# ---------------- Fig 3: E(T)
with open('td_results.json') as f:
    res = json.load(f)
fig, ax = plt.subplots(figsize=(6.4, 3.6))
for key, col in [('eps0.001_L12', 'C0'), ('eps0.001_L16', 'C1'), ('eps0.001_L20', 'C2'),
                 ('eps0.001_L24', 'C3'), ('eps0.001_L28', 'C4'), ('eps0.001_L32', 'C5')]:
    h = d[f'{key}_h']; dh_ = d[f'{key}_dh']
    n = min(len(h), len(tau0))
    tt = tau0[:n]
    norm2 = np.trapezoid(h0**2, tau0)
    integ = np.concatenate([[0.0], np.cumsum(0.5*(dh_[1:]**2+dh_[:-1]**2)*np.diff(tt))])
    E = np.sqrt(integ/norm2)
    L = float(key.split('_L')[1])
    ax.plot(tt[::20], E[::20]/1e-3, lw=0.9, color=col, label=f'$L={L:.0f}$')
ax.axhline(0.22, color='k', ls='--', lw=0.8)
ax.text(2, 0.235, r'$0.22\,\varepsilon$', fontsize=8)
ax.set_xlim(0, 90); ax.set_ylim(0, 0.3)
ax.set_xlabel(r'$\tau$'); ax.set_ylabel(r'$\mathcal{E}_{\varepsilon,L}(T)$  ($\varepsilon=10^{-3}$)')
ax.legend(fontsize=8, ncol=3)
plt.tight_layout(); plt.savefig('fig3_E_of_T.pdf'); plt.close()

# ---------------- Fig 4: pole shift scaling
with open('partB_overlap.json') as f:
    ov = json.load(f)
Wp = mp.mpc(ov['Wp'][0], ov['Wp'][1])
Ls = sorted(int(k) for k in ov['overlap'])
amp = [abs(mp.mpc(ov['overlap'][str(L)]['I_re'], ov['overlap'][str(L)]['I_im'])/Wp) for L in Ls]
Imw = 0.0889633453601956463
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.semilogy(Ls, amp, 'o-', color='C0', ms=4, label=r'$|\delta\omega_1|/\varepsilon=|I(L)/W\'(\omega_0)|$')
LL = np.linspace(12, 60, 100)
ax.semilogy(LL, 0.130*np.exp(2*Imw*LL), 'k--', lw=0.8,
            label=r'$0.130\,e^{2|\mathrm{Im}\,\omega_0|L}$')
ax.axvline(np.log(1e3)/(2*Imw), color='C3', lw=0.8, ls=':')
ax.text(np.log(1e3)/(2*Imw)+0.5, 3e-4, r'$L_*=38.8$  ($c=c_*$, $\varepsilon=10^{-3}$)', fontsize=8)
ax.set_xlabel(r'$L$'); ax.set_ylabel(r'$|\delta\omega_1|/\varepsilon$')
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig4_pole_shift.pdf'); plt.close()

# ---------------- Fig 5: pole proliferation
fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.plot([0.373672834], [-0.088963345], 'k*', ms=12, label=r'unperturbed fundamental')
ax.plot([0.3467109961], [-0.2739158157], 'k^', ms=8, label=r'unperturbed $n=1$')
try:
    rs = json.load(open('roots_smooth_eps0.03_L24.0.json'))
    ax.plot([float(r[0]) for r in rs], [float(r[1]) for r in rs], 'o', ms=6,
            mfc='none', color='C3', label=r'smooth bump $\varepsilon=0.03,\,L=24$')
except Exception as e:
    print(e)
try:
    rd = json.load(open('roots_delta_0.03.json'))
    ax.plot([float(r[0]) for r in rd], [float(r[1]) for r in rd], 's', ms=5,
            mfc='none', color='C0', label=r'delta model $\varepsilon_\delta=0.03,\,L=24$')
except Exception as e:
    print(e)
ax.set_xlim(0.15, 0.55); ax.set_ylim(-0.32, 0.02)
ax.set_xlabel(r'$\mathrm{Re}\,\omega$'); ax.set_ylabel(r'$\mathrm{Im}\,\omega$')
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig5_branches.pdf'); plt.close()

# ---------------- Fig 6: delta model root path
with open('partB_results.json') as f:
    pr = json.load(f)
eps_list = sorted(pr['delta_roots'], key=lambda s: float(s))
xs = [float(pr['delta_roots'][e][0]) for e in eps_list]
ys = [float(pr['delta_roots'][e][1]) for e in eps_list]
mu = [float(e)*0.93/8.02 for e in eps_list]
fig, ax = plt.subplots(1, 2, figsize=(6.4, 3.0))
ax[0].plot(xs, ys, 'o-', ms=3, lw=0.8, color='C0')
ax[0].plot([0.373672834], [-0.088963345], 'k*', ms=12)
ax[0].set_xlabel(r'$\mathrm{Re}\,\omega$'); ax[0].set_ylabel(r'$\mathrm{Im}\,\omega$')
ax[0].set_title(r'delta model: followed root, $L=24$', fontsize=9)
ax[1].loglog(mu, np.abs(np.array(xs)+1j*np.array(ys) - (0.373672834-0.088963345j)),
             'o-', ms=3, lw=0.8, color='C0', label='exact')
ax[1].loglog(mu, 0.116*np.array(mu), 'k--', lw=0.8, label=r'$0.116\,\mu$ (first order)')
ax[1].set_xlabel(r'$\mu=\varepsilon_\delta|\psi_L\psi_R|/|W\'|$')
ax[1].set_ylabel(r'$|\delta\omega|$')
ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig('fig6_delta_path.pdf'); plt.close()
print("figures done")
