"""
D-dimensional generalized unitarity for color-ordered one-loop amplitudes with an adjoint
(complex) scalar in the loop, D = 4 - 2 eps, loop momentum l = (l4, mu), mu^2 -> 4D mass.
Integrand reduction (OPP, harmonic basis in the transverse space), residues fitted on cuts at a
global set of mu^2 samples; rational part assembled from

  J4[mu^4] = -1/6,  J3[mu^2] = 1/2,  J3[mu^4] = (K1^2+K2^2+K3^2)/24,
  J2[mu^2] = -K^2/6, J2[mu^4] = -K^4/60      (J = int d^D l /(i pi^{D/2}), eps -> 0)

Returned value A = sum a*J  ( = 16 pi^2 * (reduced amplitude A'), one charge-flow orientation ).
Cut function: F_S = (-1)^{|S|} prod(reduced trees).
"""
import mpmath as mp
from trees import *
import itertools, random

mp.mp.dps = 40
I = mp.mpc(0, 1)

def cvec(v):
    return tuple(mp.mpc(c) for c in v)

def mpcf(q):  # Fraction -> mpc
    return mp.mpc(mp.mpf(q.numerator)/q.denominator) if hasattr(q, 'numerator') else mp.mpc(q)

# ------------------------------------------------------------------ topology
class Topology:
    """objects: gluons (list of (mom, pol)) in color order; optional graviton (eps, mom).
       propagators: labels ('A',j) offset K_{1..j}, ('B',j) offset K_{1..j}+P."""
    def __init__(self, gluons, grav=None, weight=None, mu2pow=0):
        self.gl, self.grav, self.weight, self.mu2pow = gluons, grav, weight, mu2pow
        n = len(gluons)
        self.n = n
        zero = mp.mpc(0)
        self.K = [(zero,)*4]
        for (p, e) in gluons:
            self.K.append(add(self.K[-1], p))
        if grav is None:
            self.props = [('A', j) for j in range(n)]
        else:
            self.props = [('A', j) for j in range(n+1)] + [('B', j) for j in range(1, n)]

    def offset(self, pr):
        t, j = pr
        q = self.K[j]
        if t == 'B':
            q = add(q, self.grav[1])
        return q

    def key(self, pr):
        return (pr[1], 0 if pr[0] == 'A' else 1)

    def valid(self, S):
        if self.grav is None:
            return True
        As = [j for (t, j) in S if t == 'A' and j >= 1]
        Bs = [j for (t, j) in S if t == 'B']
        if As and Bs and max(As) > min(Bs):
            return False
        return True

    def cuts(self, k):
        out = []
        for S in itertools.combinations(self.props, k):
            S = tuple(sorted(S, key=self.key))
            if self.valid(S):
                out.append(S)
        return out

    def segments(self, S):
        """list of (gluon indices (0-based), hasP) for consecutive propagators S[i] -> S[i+1] (cyclic)"""
        k, n = len(S), self.n
        segs = []
        for i in range(k-1):
            X, Y = S[i], S[i+1]
            segs.append((list(range(X[1], Y[1])), X[0] == 'A' and Y[0] == 'B'))
        X, Y = S[k-1], S[0]
        glu = list(range(X[1], n)) + list(range(0, Y[1]))
        hasP = (self.grav is not None) and (X[0] == 'A' or Y[0] == 'B')
        segs.append((glu, hasP))
        return segs

    def cut_trees(self, S, l, mu2):
        """product of reduced trees on cut S at loop momentum l (4-vector, l^2 etc not checked)"""
        k = len(S)
        segs = self.segments(S)
        prod = mp.mpc(1)
        for i in range(k):
            qa = add(l, self.offset(S[i]))
            qb = add(l, self.offset(S[(i+1) % k]))
            glu, hasP = segs[i]
            legs = [('s', qa, 1)] + [('g', self.gl[g][0], self.gl[g][1]) for g in glu] + [('sb', neg(qb), 1)]
            grav = self.grav if hasP else None
            prod *= BG(legs, grav=grav, m2=mu2, zero=mp.mpc(0)).amp()
        if self.mu2pow:
            prod *= mu2**self.mu2pow
        if self.weight is not None:      # extra numerator  c * (w . (l + q))
            c, w, q = self.weight
            prod *= c*dot(w, add(l, q))
        return (-1)**k * prod

    def seg_momenta(self, S):
        k = len(S)
        return [sub(self.offset(S[(i+1) % k]), self.offset(S[i])) for i in range(k)]
