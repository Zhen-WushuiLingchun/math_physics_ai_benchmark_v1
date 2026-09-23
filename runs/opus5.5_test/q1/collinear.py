"""collinear-limit machinery: EYM(1..n;P) vs YM(n+2) with P -> a,b"""
from kin import *
import itertools

def orderings(n):
    """Pi_n: cyclic orderings of (1..n) with a,b inserted in different gaps. gap i = after gluon i."""
    res = []
    for ga in range(1, n+1):
        for gb in range(1, n+1):
            if ga == gb:
                continue
            seq = []
            for i in range(1, n+1):
                seq.append(i)
                if ga == i: seq.append('a')
                if gb == i: seq.append('b')
            res.append(tuple(seq))
    return res

class Point:
    """kinematic point: gluons 1..n + P, real rational spinors. index n = P."""
    def __init__(self, n, rng):
        self.n = n
        self.lam, self.lt = massless_kin(n+1, rng)
        self.mom = [bisp(self.lam[i], self.lt[i]) for i in range(n+1)]
        self.rng = rng
    def s(self, i, j):   # labels 1..n, 'P'
        ii = self.n if i == 'P' else i-1
        jj = self.n if j == 'P' else j-1
        return 2*dot(self.mom[ii], self.mom[jj])

def eym_tree(pt, hel, hP):
    n = pt.n
    legs = [('g', pt.mom[i], pol(hel[i], pt.lam[i], pt.lt[i], (rnd(pt.rng), rnd(pt.rng)))) for i in range(n)]
    eP = pol(1 if hP > 0 else -1, pt.lam[n], pt.lt[n], (rnd(pt.rng), rnd(pt.rng)))
    return BG(legs, grav=(eP, pt.mom[n]), zero=Fr(0)).amp()

def ym_collinear_tree(pt, sigma, hel, hab, x):
    """C_x A(sigma) with a,b helicity hab (+1/-1). exact."""
    n = pt.n
    lamP, ltP = pt.lam[n], pt.lt[n]
    legs = []
    for lab in sigma:
        if lab == 'a':
            lt = (x*ltP[0], x*ltP[1]); lam = lamP; h = hab
        elif lab == 'b':
            lt = ((1-x)*ltP[0], (1-x)*ltP[1]); lam = lamP; h = hab
        else:
            lam, lt, h = pt.lam[lab-1], pt.lt[lab-1], hel[lab-1]
        legs.append(('g', bisp(lam, lt), pol(h, lam, lt, (rnd(pt.rng), rnd(pt.rng)))))
    A = BG(legs, zero=Fr(0)).amp()
    return A * (x**(-hab)) * ((1-x)**(-hab))

def mandelstam_basis(n):
    """independent s_ij for n gluons + P (n+1 massless): pairs (i,j) with i<j among 1..n, excluding those eliminated.
       we use s_{i,j}, 1<=i<j<=n, minus relations: for n+1 points, # independent = (n+1)(n-2)/2."""
    pairs = [(i, j) for i in range(1, n+1) for j in range(i+1, n+1)]
    return pairs  # caller will check rank; for n=3: need 2 of 3
