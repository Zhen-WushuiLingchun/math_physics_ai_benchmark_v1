"""analytic one-loop EYM (n=3), our normalization: A = sum a*J (1/16pi^2 stripped, one orientation)"""
from trees import ang, sqb
def M1_sm_n3(L, T, m):
    """single-minus on gluon m (0-based), P^{++}; L,T spinors of (g0,g1,g2,P). cyclic relabel of formula"""
    i, j, k, P = m, (m+1) % 3, (m+2) % 3, 3
    A = lambda a, b: ang(L[a], L[b]); B = lambda a, b: sqb(T[a], T[b])
    s = A(i, j)*B(i, j); u = A(i, k)*B(i, k)
    return -B(j, P)**2*B(k, P)**2*A(i, j)*A(i, k)*(s*s + u*u)/(48*A(j, k)*s*s*u*u)
