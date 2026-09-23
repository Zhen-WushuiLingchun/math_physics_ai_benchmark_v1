import numpy as np, json
d = np.load('fd_scan.npz')
w = d['w']; Ls = d['L']; epss = d['eps']
h0 = d['h0']
# Plancherel: int_0^inf |f|^2 dtau = (1/pi) int_0^inf |fhat|^2 dw  (f real).  Add high-w part of ||h0||^2 from fd_base.
b = np.load('fd_base.npz'); m = b['w'] >= 12
H0_hi = np.trapezoid(np.abs(b['h0'][m])**2, b['w'][m]) / np.pi
H0_lo = np.trapezoid(np.abs(h0)**2, w) / np.pi
H0 = H0_lo + H0_hi
print("||h0||^2 (FD) = %.6f  (w<12: %.6f, w>=12: %.6f)" % (H0, H0_lo, H0_hi))
out = dict(H0=H0, rows=[])
for i, L in enumerate(Ls):
    h1 = d['h1'][:, i]
    E1 = np.sqrt(np.trapezoid(np.abs(h1)**2, w) / np.pi / H0)
    line = "L=%5.1f  Born ||h1||/||h0|| = %.6f |" % (L, E1)
    row = dict(L=L, born=E1, full={})
    for j, e in enumerate(epss):
        dh = d['dh'][:, i, j]
        E = np.sqrt(np.trapezoid(np.abs(dh)**2, w) / np.pi / H0)
        row['full'][str(e)] = E
        line += "  eps=%.0e: E/eps=%.6f (E/eps-born)/eps=%.3f" % (e, E/e, (E/e-E1)/e)
        mb = np.max(np.abs(d['beta'][:, i, j]))
        line += " max|beta|=%.3g" % mb
    print(line)
    out['rows'].append(row)
# Large-L limit formula (Born, eps->0): kappa_inf^2 = (1/pi) int |w A phiR|^2 |What(2w)|^2/(4 w^2) dw / ||h0||^2
from rwcore import W
gx, gw = np.polynomial.legendre.leggauss(200)
What = np.array([np.sum(gw*W(gx)*np.cos(2*ww*gx)) for ww in w])
g = np.abs(w*d['Aamp']*d['phiR'])**2 * What**2/(4*w**2)
kinf = np.sqrt(np.trapezoid(g, w)/np.pi/H0)
print("L->inf Born limit kappa_inf = %.6f" % kinf)
out['kappa_inf'] = kinf
json.dump(out, open('fd_results.json', 'w'), indent=1, default=float)
