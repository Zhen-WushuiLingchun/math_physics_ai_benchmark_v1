"""collinear machinery over Q(x) (symbolic x)"""
from collinear import *
from symx import *
from sympy.polys.matrices import DomainMatrix

class PointK:
    def __init__(self, n, rng):
        p = Point(n, rng)
        self.n, self.rng, self.pt = n, rng, p
        self.lam = [tuple(toK(c) for c in l) for l in p.lam]
        self.lt = [tuple(toK(c) for c in l) for l in p.lt]
        self.mom = [bisp(self.lam[i], self.lt[i]) for i in range(n+1)]
    def s(self, i, j):
        ii = self.n if i == 'P' else i-1
        jj = self.n if j == 'P' else j-1
        return 2*dot(self.mom[ii], self.mom[jj])
    def ref(self):
        return (toK(rnd(self.rng)), toK(rnd(self.rng)))

def polK(h, lam, lt, ref):
    if (ang(ref, lam) if h > 0 else sqb(lt, ref)) == 0:
        ref = (ref[0] + 1, ref[1] + KX(sp.Rational(1, 3)))
    return eps_plus(lam, lt, ref) if h > 0 else eps_minus(lam, lt, ref)

def eym_treeK(pk, hel, hP):
    n = pk.n
    legs = [('g', pk.mom[i], polK(hel[i], pk.lam[i], pk.lt[i], pk.ref())) for i in range(n)]
    eP = polK(1 if hP > 0 else -1, pk.lam[n], pk.lt[n], pk.ref())
    return BG(legs, grav=(eP, pk.mom[n]), zero=KX(0)).amp()

def ym_col_treeK(pk, sigma, hel, hab):
    n = pk.n
    lamP, ltP = pk.lam[n], pk.lt[n]
    legs = []
    for lab in sigma:
        if lab == 'a':
            lt = (xK*ltP[0], xK*ltP[1]); lam = lamP; h = hab
        elif lab == 'b':
            lt = ((1-xK)*ltP[0], (1-xK)*ltP[1]); lam = lamP; h = hab
        else:
            lam, lt, h = pk.lam[lab-1], pk.lt[lab-1], hel[lab-1]
        legs.append(('g', bisp(lam, lt), polK(h, lam, lt, pk.ref())))
    A = BG(legs, zero=KX(0)).amp()
    return A * xK**(-hab) * (1-xK)**(-hab)

def solve_affine(rows, rhs):
    """over KX: returns (rank, consistent, particular solution, nullspace basis) using DomainMatrix"""
    A = DomainMatrix([list(r) for r in rows], (len(rows), len(rows[0])), KX)
    b = DomainMatrix([[v] for v in rhs], (len(rhs), 1), KX)
    Ab = A.hstack(b)
    R, piv = Ab.rref()
    ncol = len(rows[0])
    consistent = ncol not in piv
    rank = len([p for p in piv if p < ncol])
    Rl = R.to_Matrix()
    part = [0]*ncol
    for r, p in enumerate(piv):
        if p < ncol:
            part[p] = Rl[r, ncol]
    free = [c for c in range(ncol) if c not in piv]
    null = []
    for f in free:
        v = [0]*ncol; v[f] = 1
        for r, p in enumerate(piv):
            if p < ncol:
                v[p] = -Rl[r, f]
        null.append(v)
    return rank, consistent, part, null, piv
