"""
Scan the complex omega plane for roots of the exact perturbed resonance
function, for the delta-model and the smooth bump.
"""
import sys, json, time
import mpmath as mp
import qnm_mp as Q
from qnm_mp import *

mp.mp.dps = 40
load_asymp(2, 8)
w0 = mp.mpc('0.373672834073936005', '-0.0889633453601956463')
Vs0 = V_arrays()
w0, _ = find_qnm_mp(w0, Vs0)

Ldel = mp.mpf(24)

def F_delta(w, eps_d):
    W, yL, yR, recL, recR = rk4_full(w, V_arrays(), record=True)
    j = idx_of(Ldel)
    return W - eps_d*recL[j][0]*recR[j][0]

def F_smooth(w, eps, L):
    W = rk4_full(w, V_arrays(eps=eps, L=L), record=False)[0]
    return W

def at_X(eps):   # crude axis crossing count along Re at a few Im lines
    pass

def scan(Ffun, re_rng, im_rng, nre, nim, dx_scan='0.1'):
    Q.DX = mp.mpf(dx_scan); Q._N = int((Q.XMAX-Q.XMIN)/Q.DX)
    res = []
    t0 = time.time()
    n = 0
    tot = nre*nim
    for i in range(nre):
        re = re_rng[0] + (re_rng[1]-re_rng[0])*i/(nre-1)
        for j in range(nim):
            im = im_rng[0] + (im_rng[1]-im_rng[0])*j/(nim-1)
            w = mp.mpc(re, im)
            try:
                f = abs(Ffun(w))
            except Exception:
                f = mp.inf
            res.append((re, im, f))
            n += 1
            if n % 50 == 0:
                print(f"  ...{n}/{tot}  {time.time()-t0:.0f}s", flush=True)
    return res

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'delta'
    if which == 'delta':
        eps_d = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('3e-2')
        F = lambda w: F_delta(w, eps_d)
        tag = f'delta_{eps_d}'
    else:
        eps = mp.mpf(sys.argv[2])
        L = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf(24)
        F = lambda w: F_smooth(w, eps, L)
        tag = f'smooth_eps{eps}_L{L}'
    data = scan(F, (0.18, 0.52), (-0.30, -0.01), 18, 15)
    with open(f'scan_{tag}.json', 'w') as f:
        json.dump([[a, b, mp.nstr(c, 6)] for a, b, c in data], f)
    print("saved", f'scan_{tag}.json')
