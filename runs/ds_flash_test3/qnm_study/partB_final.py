"""
Part B final:
 (1) first-order vs exact smooth-bump pole shift (small eps)
 (2) delta-model exact resonance equation W(w) = eps_d psi_L(L) psi_R(L): root scan
     -> resummation, new branches
 (3) double-limit table
Saves partB_results.json
"""
import mpmath as mp, json, time
import qnm_mp as Q
from qnm_mp import *

mp.mp.dps = 40
load_asymp(2, 8)

def get_I(L):
    with open('partB_overlap.json') as f:
        dat = json.load(f)
    d = dat['overlap'][str(L)]
    return mp.mpc(d['I_re'], d['I_im'])

Vs0 = V_arrays()
w0 = mp.mpc('0.373672834073936005', '-0.0889633453601956463')
w0, _ = find_qnm_mp(w0, Vs0)
h = mp.mpf('1e-4')
Wp = (rk4_full(w0+h, Vs0, record=False)[0] - rk4_full(w0-h, Vs0, record=False)[0])/(2*h)
print("omega0 =", mp.nstr(w0, 18), " W' =", mp.nstr(Wp, 10))

out = {'w0': [mp.nstr(w0.real, 20), mp.nstr(w0.imag, 20)],
       'Wp': [mp.nstr(Wp.real, 16), mp.nstr(Wp.imag, 16)], 'shifts': [], 'delta_roots': {}}

# ---------------- (1) first order vs exact (smooth bump)
print("\n(1) smooth bump: first-order vs exact")
for eps_s, L in [('1e-4', 16), ('1e-4', 24), ('1e-4', 32), ('3e-4', 16), ('3e-4', 24),
                 ('1e-3', 16), ('1e-3', 24)]:
    eps = mp.mpf(eps_s)
    I = get_I(L)
    pred = eps*I/Wp
    Vs = V_arrays(eps=eps, L=mp.mpf(L))
    w_exact, it = find_qnm_mp(w0+pred, Vs, tol=mp.mpf('1e-28'))
    d = w_exact - w0
    row = {'eps': eps_s, 'L': L,
           'pred_re': mp.nstr(pred.real, 10), 'pred_im': mp.nstr(pred.imag, 10),
           'exact_re': mp.nstr(d.real, 10), 'exact_im': mp.nstr(d.imag, 10),
           'ratio_re': mp.nstr((d/pred).real, 6), 'ratio_im': mp.nstr((d/pred).imag, 6)}
    out['shifts'].append(row)
    print(f"  eps={eps_s:6s} L={L:3d}: pred={mp.nstr(pred,8):24s} exact={mp.nstr(d,8):24s} ratio={mp.nstr(d/pred,6)}")

# ---------------- (2) delta model
print("\n(2) delta model: exact roots W(w) = eps_d psi_L(L) psi_R(L)")
Ldel = mp.mpf(24)
def delta_F(w, eps_d, rec=None):
    if rec is None:
        W, yL, yR, recL, recR = rk4_full(w, V_arrays(), record=True)
    else:
        recL, recR = rec
        W = rk4_full(w, V_arrays(), record=False)[0]
    j = idx_of(Ldel)
    K = recL[j][0]*recR[j][0]
    return W - eps_d*K, K

# prepare the recorded profiles at various w is expensive; do it per eval instead:
def delta_F2(w, eps_d):
    W, yL, yR, recL, recR = rk4_full(w, V_arrays(), record=True)
    j = idx_of(Ldel)
    K = recL[j][0]*recR[j][0]
    return W - eps_d*K

# roots along the eps_d family (follow the fundamental continuously)
eps_list = ['1e-3', '1e-2', '3e-2', '1e-1', '3e-1', '1e0', '3e0', '1e1']
w_cur = w0
print("  eps_d      root (followed)                 |eps_d K/W'| estimate")
for es in eps_list:
    eps_d = mp.mpf(es)
    w0_guess = w_cur
    # Newton
    w_cur_try = w0_guess
    for it in range(60):
        F0 = delta_F2(w_cur_try, eps_d)
        hh = mp.mpf('1e-6')
        dF = (delta_F2(w_cur_try+hh, eps_d) - delta_F2(w_cur_try-hh, eps_d))/(2*hh)
        step = F0/dF
        w_cur_try -= step
        if abs(step) < mp.mpf('1e-20'):
            break
    w_cur = w_cur_try
    # strength estimate at the unperturbed w0:
    _, K0 = delta_F(w0, mp.mpf(0))
    mu = eps_d*abs(K0)/abs(Wp)
    out['delta_roots'][es] = [mp.nstr(w_cur.real, 12), mp.nstr(w_cur.imag, 12)]
    print(f"  {es:8s}  {mp.nstr(w_cur, 14):30s}  mu={mp.nstr(mu, 6)}")

with open('partB_results.json', 'w') as f:
    json.dump(out, f, indent=1)
print("\nsaved partB_results.json")
