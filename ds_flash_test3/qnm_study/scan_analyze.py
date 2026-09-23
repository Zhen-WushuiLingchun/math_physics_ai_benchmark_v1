"""Find local minima in a scan grid and Newton-refine to exact roots."""
import sys, json
import numpy as np
import mpmath as mp
import qnm_mp as Q
from qnm_mp import *

mp.mp.dps = 40
load_asymp(2, 8)
w0 = mp.mpc('0.373672834073936005', '-0.0889633453601956463')
Vs0 = V_arrays()
w0, _ = find_qnm_mp(w0, Vs0)
Ldel = mp.mpf(24)

def F_delta(w, eps_d, dx=None):
    if dx is not None:
        Q.DX = mp.mpf(dx); Q._N = int((Q.XMAX-Q.XMIN)/Q.DX)
    W, yL, yR, recL, recR = rk4_full(w, V_arrays(), record=True)
    j = idx_of(Ldel)
    return W - eps_d*recL[j][0]*recR[j][0]

def F_smooth(w, eps, L, dx=None):
    if dx is not None:
        Q.DX = mp.mpf(dx); Q._N = int((Q.XMAX-Q.XMIN)/Q.DX)
    return rk4_full(w, V_arrays(eps=eps, L=L), record=False)[0]

def newton(F, w, itmax=25):
    w = mp.mpc(w)
    h = mp.mpf('1e-6')
    for _ in range(itmax):
        F0 = F(w)
        dF = (F(w+h)-F(w-h))/(2*h)
        step = F0/dF
        w -= step
        if abs(step) < mp.mpf('1e-22'):
            break
    return w

if __name__ == '__main__':
    tag = sys.argv[1]
    if tag.startswith('delta'):
        eps_d = mp.mpf(tag.split('_')[1])
        F = lambda w: F_delta(w, eps_d)
        F05 = lambda w: F_delta(w, eps_d, dx='0.05')
    else:
        # smooth_eps0.03_L24
        parts = tag.split('_')
        eps = mp.mpf(parts[1][3:]); L = mp.mpf(parts[2][1:])
        F = lambda w: F_smooth(w, eps, L)
        F05 = lambda w: F_smooth(w, eps, L, dx='0.05')
    data = json.load(open(f'scan_{tag}.json'))
    re = np.array([d[0] for d in data]); im = np.array([d[1] for d in data])
    f = np.array([float(mp.mpf(d[2])) for d in data])
    nre = len(set(re)); nim = len(set(im))
    Fg = f.reshape(nre, nim)
    roots = []
    for i in range(1, nre-1):
        for j in range(1, nim-1):
            v = Fg[i, j]
            if v < Fg[i-1, j] and v < Fg[i+1, j] and v < Fg[i, j-1] and v < Fg[i, j+1]:
                w_guess = mp.mpc(float(re[i*nim+j]), float(im[i*nim+j]))
                try:
                    wr = newton(F05, w_guess)
                    wr = newton(F, wr)
                    # keep only true zeros
                    if abs(F(wr)) < mp.mpf('1e-8'):
                        roots.append(wr)
                except Exception as e:
                    print("  fail", w_guess, e)
    # dedupe
    uniq = []
    for r in roots:
        if all(abs(r-u) > mp.mpf('1e-3') for u in uniq):
            uniq.append(r)
    print(f"tag={tag}: {len(uniq)} roots")
    for r in sorted(uniq, key=lambda z: (z.real, z.imag)):
        print(f"   {mp.nstr(r, 12)}   |F|={mp.nstr(abs(F(r)),3)}")
    with open(f'roots_{tag}.json','w') as fo:
        json.dump([[mp.nstr(r.real,12), mp.nstr(r.imag,12)] for r in uniq], fo)
