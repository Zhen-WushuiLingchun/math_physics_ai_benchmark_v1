"""
Fibonacci fusion-category setup: verify Hom-space dimensions.
Conventions:
  A[n] = dim Hom(1, tau^{otimes n}),  B[n] = dim Hom(tau, tau^{otimes n})
  A[0]=1, B[0]=0;  A[n] = B[n-1];  B[n] = A[n-1]+B[n-1]
  Fibonacci numbers F_0=0, F_1=1, F_n = F_{n-1}+F_{n-2}
  Claim: A[n] = F_{n-1} (n>=1), B[n] = F_n, and dim Hom(1,tau^N)=D_N=F_{N-1}.
  Split identity: A[N] = A[NR]*A[NB] + B[NR]*B[NB]  (N = NR + NB)
"""
from sympy import fibonacci, simplify, Rational, sqrt, nsimplify, Symbol

def hom_dims(N):
    A = [0]*(N+1); B = [0]*(N+1)
    A[0] = 1; B[0] = 0
    for n in range(1, N+1):
        A[n] = B[n-1]
        B[n] = A[n-1] + B[n-1]
    return A, B

A, B = hom_dims(20)
ok = True
for n in range(0, 19):
    fa = fibonacci(n-1) if n >= 1 else None
    fb = fibonacci(n)
    if n >= 1 and A[n] != fa: ok = False; print("A mismatch", n)
    if B[n] != fb: ok = False; print("B mismatch", n)
print("Hom dims match Fibonacci:", ok)
print("n      :", list(range(0, 13)))
print("A[n]   :", A[:13])
print("B[n]   :", B[:13])

# split identity
print("\nSplit identity A[NR+NB] = A[NR]A[NB] + B[NR]B[NB]:")
for NR in range(0, 9):
    for NB in range(0, 9):
        assert A[NR+NB] == A[NR]*A[NB] + B[NR]*B[NB], (NR, NB)
print("  verified for all NR,NB <= 8")

# explicit sector multiplicities for the problem (NR,NB >= 2)
print("\nSector data m_a, n_a, D_N for small (NR,NB):")
for NR in range(2, 8):
    for NB in range(2, 8):
        N = NR + NB
        m1, mt = A[NR], B[NR]          # a=1: A[NR] ; a=tau: B[NR]
        n1, nt = A[NB], B[NB]
        D = m1*n1 + mt*nt
        DR = m1 + mt
        DB = n1 + nt
        assert D == A[N] == fibonacci(N-1)
        assert DR == fibonacci(NR+1) and DB == fibonacci(NB+1)
print("  verified: D_N = F_{N-1}, DR = F_{NR+1}, DB = F_{NB+1}")
for (NR, NB) in [(2,2),(3,2),(3,3),(4,3),(5,4),(6,5)]:
    N = NR + NB
    m1, mt = A[NR], B[NR]; n1, nt = A[NB], B[NB]
    D = m1*n1 + mt*nt
    print(f"  (NR,NB)=({NR},{NB}): m=(1:{m1}, tau:{mt}), n=(1:{n1}, tau:{nt}), D={D}, DR={m1+mt}, DB={n1+nt}")

# golden ratio exact checks with Binet: check A[n] = (phi^n - psi^n)/sqrt5 etc.
phi = (1 + sqrt(5))/2
psi = (1 - sqrt(5))/2
print("\nBinet check A[n]=F_{n-1}:")
for n in [5, 9]:
    val = simplify((phi**(n-1) - psi**(n-1))/sqrt(5))
    print(f"  n={n}: {val}  == {fibonacci(n-1)}")
print("phi =", float(phi))
print("log_phi(2) =", float(1/(phi*0+1) * 0) if False else float(__import__('math').log(2)/__import__('math').log(float(phi))))
