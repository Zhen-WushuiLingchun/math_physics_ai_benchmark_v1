"""analytic one-loop YM (scalar loop) in OUR normalization (verified by the D-dim engine):
   all-plus N-pt:  A = -(-1)^N/6 * sum tr_-(i1 i2 i3 i4)/PT
   single-minus 5-pt (minus first): A = (-1)^N/6 * F5_BDK  (N=5 -> -F5/6)"""
import itertools
from trees import ang, sqb

def allplus(L, T):
    N = len(L)
    num = 0
    for a, b, c, d in itertools.combinations(range(N), 4):
        num += ang(L[a], L[b])*sqb(T[b], T[c])*ang(L[c], L[d])*sqb(T[d], T[a])
    den = 1
    for i in range(N):
        den *= ang(L[i], L[(i+1) % N])
    return -(-1)**N*num/den/6

def F5(L, T):
    A = lambda i, j: ang(L[i-1], L[j-1]); B = lambda i, j: sqb(T[i-1], T[j-1])
    return (1/A(3, 4)**2 if not isinstance(A(3,4), int) else 1)*(-B(2, 5)**3/(B(1, 2)*B(5, 1)) + A(1, 4)**3*B(4, 5)*A(3, 5)/(A(1, 2)*A(2, 3)*A(4, 5)**2)
                          - A(1, 3)**3*B(3, 2)*A(4, 2)/(A(1, 5)*A(5, 4)*A(3, 2)**2))

def singleminus5(L, T, m):
    """minus at position m (0-based) in the ordered list"""
    L2 = L[m:] + L[:m]; T2 = T[m:] + T[:m]
    return -F5(L2, T2)/6

def col_spinors(lam, lt, sigma, x):
    """sigma labels: ints 1..n (gluons) or 'a','b'; lam/lt: gluons then P. returns lists + helicity factor for a^+b^+"""
    n = len(lam) - 1
    L, T = [], []
    for lab in sigma:
        if lab == 'a':
            L.append(lam[n]); T.append((x*lt[n][0], x*lt[n][1]))
        elif lab == 'b':
            L.append(lam[n]); T.append(((1-x)*lt[n][0], (1-x)*lt[n][1]))
        else:
            L.append(lam[lab-1]); T.append(lt[lab-1])
    return L, T

def CxA1(lam, lt, sigma, x, minus=None):
    """C_x A^(1)(sigma) with a^+ b^+; minus = gluon label with - helicity or None (all plus)."""
    L, T = col_spinors(lam, lt, sigma, x)
    if minus is None:
        A = allplus(L, T)
    else:
        A = singleminus5(L, T, sigma.index(minus))
    return A/(x*(1-x))
