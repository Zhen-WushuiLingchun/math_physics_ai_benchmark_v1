"""
Color-ordered tree amplitudes for YM + adjoint complex scalar (mass m) + one external
graviton coupled at O(kappa), derived from

  L = tr[-1/4 F^2] + tr[(D Phi)^+ D Phi] - m^2 tr[Phi^+ Phi]
      + kappa tr[ 1/2 h^{mu a} F_{mu nu} F_a^nu - h^{mu nu} (D_mu Phi)^+ D_nu Phi ]
  F = dA - dA - i[A,A],  D Phi = d Phi - i[A,Phi],  tr(T^a T^b) = delta^{ab}  (g' = 1)

Metric (+,-,-,-).  All momenta incoming.  All i's removed ("reduced" rules):
  vertex V = i V',  gluon prop = i(-eta/P^2),  scalar prop = i/(P^2-m^2)
  current J = -Pi' V'  ;  tree amplitude  A = i A'   (functions below return A').
Vectors are bispinors (a11,a12,a21,a22) = p_{alpha alphadot}; dot(a,b) = a.b.
Number type is generic (Fraction for exact real split-signature kinematics, mpc for loops).
"""
from fractions import Fraction
from functools import lru_cache

# ---------------- vectors -----------------
def dot(a, b):
    return (a[0]*b[3] + a[3]*b[0] - a[1]*b[2] - a[2]*b[1]) / 2

def add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3])
def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3])
def scal(c, a): return (c*a[0], c*a[1], c*a[2], c*a[3])
def neg(a): return (-a[0], -a[1], -a[2], -a[3])

def vsum(vs, zero):
    s = (zero, zero, zero, zero)
    for v in vs:
        s = add(s, v)
    return s

def bisp(l, lt):  # |l> [lt|  -> p_{alpha alphadot}
    return (l[0]*lt[0], l[0]*lt[1], l[1]*lt[0], l[1]*lt[1])

def ang(a, b): return a[0]*b[1] - a[1]*b[0]
def sqb(a, b): return a[0]*b[1] - a[1]*b[0]

def inv(x):
    return Fraction(1, x) if isinstance(x, int) else 1/x

def eps_plus(lam, lamt, q):      # q = reference lambda;  eps+ = |q>[k| / <q k>
    return scal(inv(ang(q, lam)), bisp(q, lamt))

def eps_minus(lam, lamt, qt):    # qt = reference lambda-tilde;  eps- = |k>[q| / [k q]
    return scal(inv(sqb(lamt, qt)), bisp(lam, qt))

# dual basis: V(e) linear -> vector V with V(e) = V.e
E = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
def dual_vec(f):
    return (2*f(E[3]), -2*f(E[2]), -2*f(E[1]), 2*f(E[0]))

# ---------------- lines -----------------
# a line is (kind, mom, pol); kind in 'g','s','sb'; pol: vector for 'g', scalar 1 for s/sb

def rot(lst, k):
    return lst[k:] + lst[:k]

def V3g(L):
    tot = 0
    for r in range(3):
        a, b, c = rot(L, r)
        tot += dot(a[1], b[2])*dot(a[2], c[2]) - dot(a[2], b[2])*dot(a[1], c[2])
    return tot

def V4g(L):
    e1, e2, e3, e4 = [l[2] for l in L]
    return 2*dot(e1, e3)*dot(e2, e4) - dot(e1, e4)*dot(e2, e3) - dot(e1, e2)*dot(e3, e4)

def _norm_sb(L):
    k = [l[0] for l in L].index('sb')
    return rot(L, k)

def V3s(L):
    L = _norm_sb(L)
    sb, x, y = L
    if x[0] == 'g':   # (sb, g, s)
        g, s = x, y
        return dot(sub(s[1], sb[1]), g[2])
    else:             # (sb, s, g)
        s, g = x, y
        return dot(sub(sb[1], s[1]), g[2])

def V4s(L):
    L = _norm_sb(L)
    kinds = ''.join(l[0][0] for l in L[1:])   # pattern of remaining three: e.g. 'sgg','ggs','gsg'
    gl = [l for l in L if l[0] == 'g']
    ee = dot(gl[0][2], gl[1][2])
    if kinds == 'sgg' or kinds == 'ggs':
        return ee
    elif kinds == 'gsg':
        return -2*ee
    raise ValueError(kinds)

# graviton vertices; eps = graviton polarization vector (h^{mu nu} = eps^mu eps^nu)
def fvec(eps, l):
    return sub(scal(dot(eps, l[1]), l[2]), scal(dot(eps, l[2]), l[1]))

def Vh2g(eps, L):
    return -dot(fvec(eps, L[0]), fvec(eps, L[1]))

def Vh3g(eps, L):
    tot = 0
    for r in range(3):
        a, b, c = rot(L, r)
        fa = fvec(eps, a)
        tot += -dot(eps, b[2])*dot(fa, c[2]) + dot(fa, b[2])*dot(eps, c[2])
    return tot

def Vh4g(eps, L):
    tot = 0
    for r in range(4):
        a, b, c, d = [l[2] for l in rot(L, r)]
        m1 = dot(eps, a)*dot(b, d)*dot(eps, c)
        m2 = dot(eps, a)*dot(b, c)*dot(eps, d)
        m3 = dot(a, d)*dot(eps, b)*dot(eps, c)
        m4 = dot(a, c)*dot(eps, b)*dot(eps, d)
        tot += -(m1 - m2 - m3 + m4)/2
    return tot

def Vh2s(eps, L):
    L = _norm_sb(L)
    sb, s = L
    return dot(eps, sb[1])*dot(eps, s[1])

def Vh3s(eps, L):
    L = _norm_sb(L)
    sb, x, y = L
    if x[0] == 'g':
        g, s = x, y
        return dot(eps, g[2])*dot(eps, sub(sb[1], s[1]))
    else:
        s, g = x, y
        return dot(eps, g[2])*dot(eps, sub(s[1], sb[1]))

def Vh4s(eps, L):
    L = _norm_sb(L)
    kinds = ''.join(l[0][0] for l in L[1:])
    gl = [l for l in L if l[0] == 'g']
    ee = dot(eps, gl[0][2])*dot(eps, gl[1][2])
    if kinds == 'gsg':
        return 2*ee
    elif kinds in ('sgg', 'ggs'):
        return -ee
    raise ValueError(kinds)

def vertex(L, eps=None):
    """reduced color-ordered vertex V' for cyclically ordered lines L; eps!=None -> graviton attached"""
    kinds = [l[0] for l in L]
    ns = kinds.count('s') + kinds.count('sb')
    n = len(L)
    if eps is None:
        if ns == 0:
            return V3g(L) if n == 3 else V4g(L)
        return V3s(L) if n == 3 else V4s(L)
    else:
        if ns == 0:
            return {2: Vh2g, 3: Vh3g, 4: Vh4g}[n](eps, L)
        return {2: Vh2s, 3: Vh3s, 4: Vh4s}[n](eps, L)


# ---------------- Berends-Giele -----------------
class BG:
    """legs: list of (kind, mom, pol) in color order.  kinds: all 'g', or legs[0]='s', legs[-1]='sb'.
       grav: None or (eps, mom).  m2: scalar mass^2.  zero: number-type zero."""
    def __init__(self, legs, grav=None, m2=0, zero=0):
        self.legs, self.grav, self.m2, self.zero = legs, grav, m2, zero
        self.memo = {}
        self.scalar = legs[0][0] == 's'

    def mom(self, i, j, h):
        P = vsum([self.legs[k][1] for k in range(i, j+1)], self.zero)
        if h:
            P = add(P, self.grav[1])
        return P

    def cur(self, i, j, h):
        """current for legs i..j (+graviton if h). returns (kind, P, J)"""
        key = (i, j, h)
        if key in self.memo:
            return self.memo[key]
        P = self.mom(i, j, h)
        isS = self.scalar and i == 0
        if i == j and not h:
            res = self.legs[i]
        else:
            if isS:
                V = self.vsum(i, j, h, ('sb', neg(P), 1))
                J = -V/(dot(P, P) - self.m2)
                res = ('s', P, J)
            else:
                f = lambda e: self.vsum(i, j, h, ('g', neg(P), e))
                V = dual_vec(f)
                res = ('g', P, scal(1/dot(P, P), V))
        self.memo[key] = res
        return res

    def vsum(self, i, j, h, out):
        """sum of all vertices joining sub-currents of legs i..j (and graviton if h) to line 'out'"""
        tot = self.zero
        eps = self.grav[0] if self.grav else None
        n = j - i + 1
        # splits into 2 and 3 consecutive nonempty pieces
        splits2 = [((i, k), (k+1, j)) for k in range(i, j)]
        splits3 = [((i, k), (k+1, l), (l+1, j)) for k in range(i, j) for l in range(k+1, j)]
        def lines(parts, hs):
            out_l = []
            for (a, b), hh in zip(parts, hs):
                out_l.append(self.cur(a, b, hh))
            return out_l
        def val(L, e=None):
            v = vertex(L, e)
            for l in L:
                if l[0] in ('s', 'sb'):
                    v = v*l[2]
            return v
        for parts in splits2 + splits3:
            np_ = len(parts)
            if not h:
                tot += val(lines(parts, [False]*np_) + [out])
            else:
                for q in range(np_):
                    hs = [False]*np_; hs[q] = True
                    tot += val(lines(parts, hs) + [out])
                tot += val(lines(parts, [False]*np_) + [out], eps)
        if h:
            # graviton attaches to the line carrying all of i..j
            tot += val([self.cur(i, j, False), out], eps)
        return tot

    def amp(self):
        n = len(self.legs)
        h = self.grav is not None
        last = self.legs[-1]
        if self.scalar:
            return self.vsum(0, n-2, h, last)
        else:
            f = lambda e: self.vsum(0, n-2, h, ('g', last[1], e))
            V = dual_vec(f)
            return dot(V, last[2])
