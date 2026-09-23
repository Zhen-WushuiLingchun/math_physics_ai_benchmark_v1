"""Singular-value spectra: coefficient matrix V vs augmented [V|Delta] (obstruction = one extra
singular value above the noise floor). Data: spectra_n3.json (ext_analyze.py), spectra_n4.json (n4_analyze.py)."""
import json, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

font_manager.fontManager.addfont(r'C:\Windows\Fonts\Deng.ttf')      # DengXian (static)
plt.rcParams.update({
    'font.family': ['DengXian', 'DejaVu Sans'], 'mathtext.fontset': 'dejavusans',
    'pdf.fonttype': 42, 'font.size': 7.5, 'axes.linewidth': 0.6,
    'axes.edgecolor': '#52514e', 'xtick.color': '#52514e', 'ytick.color': '#52514e',
    'axes.labelcolor': '#0b0b0b', 'text.color': '#0b0b0b',
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
})
C_V, C_AUG = '#2a78d6', '#eb6834'          # validated categorical slots 1, 2 (light)
INK2, GRID = '#52514e', '#e4e3df'

n3 = json.load(open('spectra_n3.json')); n4 = json.load(open('spectra_n4.json'))
panels = [
    (n3['none'], r'$n=3$：仅 $\mathcal{N}_3$', '8 个未知数，48 个方程'),
    (n3['mu'],   r'$n=3$：$\mathcal{N}_3+E_\mu$', '20 个未知数，48 个方程'),
    (n3['l'],    r'$n=3$：$\mathcal{N}_3+E_\ell$', '20 个未知数，48 个方程'),
    (n3['both'], r'$n=3$：$\mathcal{N}_3+E_\mu+E_\ell$', '32 个未知数，48 个方程'),
    (n4,         r'$n=4$：$\mathcal{N}_4$（直接检验）', '21 个未知数，40 个方程'),
]
NOISE = -15   # log10 threshold separating signal from numerical zero (all data: gap >= 1e-18 around it)

fig, axes = plt.subplots(2, 3, figsize=(6.3, 4.1), sharey=True)
axes = axes.ravel()
for ax, (d, title, sub) in zip(axes, panels):
    V = [math.log10(max(v, 1e-60)) for v in d['V']]
    A = [math.log10(max(v, 1e-60)) for v in d['aug']]
    rV = sum(1 for v in V if v > NOISE); rA = sum(1 for v in A if v > NOISE)
    ax.axhspan(-60, NOISE, color='#f1f0ec', zorder=0, lw=0)
    ax.axhline(NOISE, color=INK2, lw=0.5, ls=(0, (3, 2)), zorder=1)
    ax.grid(axis='y', color=GRID, lw=0.4, zorder=0)
    ax.scatter([k + 0.78 for k in range(len(V))], V, s=9, color=C_V, lw=0, zorder=3)
    ax.scatter([k + 1.22 for k in range(len(A))], A, s=9, color=C_AUG, lw=0, zorder=3)
    # the extra augmented singular value = first one missing a partner above the noise line
    kx = rA - 1
    ax.annotate(r'新增 $\sigma=%.1f\times10^{%d}$' % (10**(A[kx] - math.floor(A[kx])), math.floor(A[kx])),
                xy=(kx + 1.22, A[kx]), xytext=(max(0.4, kx + 1.22 - 0.62*len(A)), -9.6),
                fontsize=6.6, color='#0b0b0b', va='center',
                arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2, shrinkA=0, shrinkB=2.5))
    ax.set_title(title, fontsize=8, pad=3)
    ax.text(0.03, 0.035, '%s\n秩 %d → 增广 %d' % (sub, rV, rA), transform=ax.transAxes,
            fontsize=6.3, color=INK2, va='bottom')
    ax.set_xlim(0, len(A) + 1)
    ax.set_ylim(-56, 4)
    ax.set_yticks([0, -10, -20, -30, -40, -50])
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
for ax in axes[3:5]:
    ax.set_xlabel('奇异值序号 $k$')
for ax in (axes[0], axes[3]):
    ax.set_ylabel(r'$\log_{10}(\sigma_k/\sigma_1)$')
# legend panel
lg = axes[5]
lg.axis('off')
lg.scatter([0.08], [0.80], s=16, color=C_V, transform=lg.transAxes, clip_on=False)
lg.text(0.16, 0.80, r'$V$：一圈方程在零核方向上的系数', transform=lg.transAxes, va='center', fontsize=7)
lg.scatter([0.08], [0.66], s=16, color=C_AUG, transform=lg.transAxes, clip_on=False)
lg.text(0.16, 0.66, r'$[V\,|\,\Delta]$：加上缺陷列（增广）', transform=lg.transAxes, va='center', fontsize=7)
lg.text(0.04, 0.53, '灰区：数值零（噪声），虚线 $10^{-15}$', transform=lg.transAxes, va='center', fontsize=6.6, color=INK2)
lg.text(0.04, 0.42, '增广后多出一个远高于噪声的奇异值，\n说明缺陷 $\\Delta$ 不在 $V$ 的列空间中，\n即该扩张下方程组不相容。\n所有面板均取 $x=2/7$。',
        transform=lg.transAxes, va='top', fontsize=6.6, color='#0b0b0b', linespacing=1.45)
fig.tight_layout(pad=0.4, w_pad=0.6, h_pad=1.0)
fig.savefig('fig_spectra.pdf')
fig.savefig('fig_spectra.png', dpi=200)
print('saved')
