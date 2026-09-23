import sympy as sp
from fractions import Fraction as Fr
def exact_pade(pts, maxtot=20):
    """pts [(r,g)] Fractions. minimal P/Q exact."""
    r_ = sp.Symbol('r')
    for tot in range(0, maxtot+1):
        for dp in range(0, tot+1):
            dq = tot - dp
            if dp + dq + 2 > len(pts):
                continue
            M = sp.Matrix([[sp.Rational(r.numerator, r.denominator)**k for k in range(dp+1)] +
                           [-sp.Rational(g.numerator, g.denominator)*sp.Rational(r.numerator, r.denominator)**k for k in range(dq+1)] for r, g in pts])
            ns = M.nullspace()
            if len(ns) >= 1:
                v = ns[0]
                P = sum(v[k]*r_**k for k in range(dp+1)); Q = sum(v[dp+1+k]*r_**k for k in range(dq+1))
                return sp.factor(sp.cancel(P/Q)), dp, dq, len(ns)
    return None
