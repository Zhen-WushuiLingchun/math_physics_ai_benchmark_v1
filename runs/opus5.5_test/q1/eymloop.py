from opp import *
from collinear import *

def mpv(v): return tuple(mpcf(c) for c in v)

def mp_gluon(lam, lt, h, rng):
    l, t = mpv(lam), mpv(lt)
    r = (mp.mpc(rng.uniform(-2, 2)), mp.mpc(rng.uniform(-2, 2)))
    return (bisp(l, t), eps_plus(l, t, r) if h > 0 else eps_minus(l, t, r))

def eym_one_loop(pt, hel, hP, rng, seed=1):
    n = pt.n
    gl = [mp_gluon(pt.lam[i], pt.lt[i], hel[i], rng) for i in range(n)]
    pP, eP = mp_gluon(pt.lam[n], pt.lt[n], 1 if hP > 0 else -1, rng)
    return one_loop(gl, grav=(eP, pP), seed=seed)

def ym_col_one_loop(pt, sigma, hel, hab, x, rng, seed=1):
    """C_x A^(1)(sigma), x an mp number or Fraction"""
    n = pt.n
    lamP, ltP = pt.lam[n], pt.lt[n]
    xm = mpcf(x) if hasattr(x, 'numerator') else x
    gl = []
    for lab in sigma:
        if lab == 'a':
            gl.append(mp_gluon(lamP, (ltP[0]*x, ltP[1]*x), hab, rng))
        elif lab == 'b':
            gl.append(mp_gluon(lamP, (ltP[0]*(1-x), ltP[1]*(1-x)), hab, rng))
        else:
            gl.append(mp_gluon(pt.lam[lab-1], pt.lt[lab-1], hel[lab-1], rng))
    tot, parts, diag = one_loop(gl, seed=seed)
    return tot * xm**(-hab) * (1-xm)**(-hab), parts, diag
