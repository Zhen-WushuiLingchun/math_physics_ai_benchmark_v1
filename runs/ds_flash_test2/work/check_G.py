import math
phi=(1+math.sqrt(5))/2
def fib(n):
    p,q=0,1
    for _ in range(n): p,q=q,p+q
    return p
def G(NR,NB):
    m1,mt=fib(NR-1),fib(NR); n1,nt=fib(NB-1),fib(NB)
    D=m1*n1+mt*nt; A=m1*m1*n1+mt*mt*nt; B=m1*n1*n1+mt*nt*nt
    return (D*B-A)/(D*(D*D-1))*phi**NR
print("G = [E Tr C^2]_k / ((k^2-1)/k^2) * phi^{NR},  limit 2phi/sqrt5 =",2*phi/math.sqrt(5))
print(f"{'(NR,NB)':>10} {'delta':>6} {'G':>14} {'G-limit':>12} {'phi^-NB':>10} {'phi^-2NB':>10}")
for NB in [4,6,8,10,12,14]:
    for NR in [NB, NB+1, NB+2]:
        g=G(NR,NB)
        print(f"({NR:>3},{NB:>3}) {NR-NB:>6} {g:>14.9f} {g-2*phi/math.sqrt(5):>12.3e} {phi**-NB:>10.3e} {phi**(-2*NB):>10.3e}")
