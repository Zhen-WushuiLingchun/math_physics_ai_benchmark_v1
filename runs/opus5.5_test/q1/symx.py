import sympy as sp
X = sp.Symbol('x')
KX = sp.QQ.frac_field(X)
xK = KX.gens[0]
def toK(q):  # Fraction -> KX
    return KX(sp.Rational(q.numerator, q.denominator)) if hasattr(q, 'numerator') else KX(q)
def fromK(e):
    return KX.to_sympy(e)
