import numpy as np, json, time
import qnm
from scipy.special import lambertw
from q_roots import F, muller, born_b, w0, Aout0, dA0
g0 = -w0.imag
res = {}
# (1) cavity branches from Lambert W_k for one strongly perturbed case
for eps, L in [(1e-6, 80.), (1e-4, 50.)]:
    d1 = eps*born_b(w0, L)*Aout0/dA0
    z = -2j*L*d1
    lst = []
    for k in range(-4, 5):
        dk = 1j*lambertw(z, k)/(2*L)
        g = w0 + dk
        r, it = muller(lambda w: F(w, eps, L), g, g*(1+1e-4), g*(1-1e-4)+1e-4j)
        lst.append(dict(k=k, lw=[g.real, g.imag], exact=[r.real, r.imag], it=it))
        print("eps=%g L=%g k=%+d  LambertW %s  exact %s it=%d" % (eps, L, k, np.round(g, 5), np.round(r, 5), it), flush=True)
    res['branches_%g_%g' % (eps, L)] = lst
# (2) double limit L = c log(1/eps): track k=0 root
for c in [3.0, 4.0, 5.0, 6.5, 8.0]:
    lst = []
    for p in [2, 3, 4, 6, 8, 10]:
        eps = 10.0**(-p); L = c*np.log(1/eps)
        if L < 13: continue
        d1 = eps*born_b(w0, L)*Aout0/dA0
        z = -2j*L*d1; Lam = abs(2*L*d1)
        g = w0 + 1j*lambertw(z, 0)/(2*L)
        r, it = muller(lambda w: F(w, eps, L), g, g*(1+1e-4), g*(1-1e-4)+1e-4j)
        lst.append(dict(eps=eps, L=L, Lambda=Lam, d1=[d1.real, d1.imag], root=[r.real, r.imag], it=it))
        print("c=%.1f eps=1e-%d L=%.1f Lambda=%.3g |d1|=%.3e  root=%s  |root-w0|=%.3e  Im shift=%.4f  (-1/2c=%.4f) it=%d" % (
            c, p, L, Lam, abs(d1), np.round(r, 5), abs(r-w0), r.imag-w0.imag, -1/(2*c), it), flush=True)
    res['c=%g' % c] = lst
json.dump(res, open('roots_branches.json', 'w'), indent=1)
