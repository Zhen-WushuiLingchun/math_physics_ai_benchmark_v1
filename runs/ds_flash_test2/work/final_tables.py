"""Final verification tables with high sample counts (for the report)."""
import numpy as np, math
from scipy.special import digamma

rng = np.random.default_rng(20260924)
phi=(1+math.sqrt(5))/2; lphi=math.log(phi)
def fib(n):
    p,q=0,1
    for _ in range(n): p,q=q,p+q
    return p
def specs(NR,NB): return [fib(NR-1),fib(NR)],[fib(NB-1),fib(NB)]
def rand_isometry(D,k):
    A=rng.normal(size=(D,k))+1j*rng.normal(size=(D,k))
    Q,R=np.linalg.qr(A); d=np.diagonal(R).copy(); d/=np.abs(d); return Q*d.conj()
def blocks(V,m,n):
    X=[];off=0
    for ma,na in zip(m,n):
        X.append(V[off:off+ma*na,:].reshape(ma,na,-1).transpose(2,0,1)); off+=ma*na
    return X
def C_trace2(X,k):
    tq=0.0; tb=0.0
    for x in X:
        G=np.einsum('jmn,lmq->jlnq',x.conj(),x)
        tq+=np.sum(np.abs(G)**2)
        P=np.einsum('jmn,jmq->jnq',x.conj(),x)
        tb+=np.einsum('jnq,lqn->jl',P,P).sum()
    return float(np.real((tq-tb/k)/k**2))
def etrC2(NR,NB,k):
    m,n=specs(NR,NB); D=sum(a*b for a,b in zip(m,n))
    A=sum(a*a*b for a,b in zip(m,n)); B=sum(a*b*b for a,b in zip(m,n))
    return (k*k-1)*(D*B-A)/(k*k*D*(D*D-1))
def page_S(m,n):
    if m>n: m,n=n,m
    return sum(1.0/r for r in range(n+1,m*n+1))-(m-1)/(2.0*n)
def ES_alg(NR,NB):
    m,n=specs(NR,NB); D=sum(a*b for a,b in zip(m,n)); v=digamma(D+1)
    for ma,na in zip(m,n):
        al=ma*na; v-=(al/D)*digamma(al+1); v+=(al/D)*page_S(ma,na)
    return v
def ES_qtr(NR,NB):
    m,n=specs(NR,NB); D=sum(a*b for a,b in zip(m,n))
    return ES_alg(NR,NB)+(m[1]*n[1]/D)*lphi
def ETr2_exact(NR,NB):
    m,n=specs(NR,NB); D=sum(a*b for a,b in zip(m,n))
    return sum(a*b*(a+b) for a,b in zip(m,n))/(D*(D+1))

print("=== Table B1: E Tr C_QB^2 : Monte Carlo vs exact ===")
Ns=20000
for (NR,NB,k) in [(2,2,2),(3,3,2),(3,3,4),(4,3,2),(5,4,2),(5,4,3),(6,5,3),(4,2,2),(6,3,3)]:
    m,n=specs(NR,NB); D=fib(NR+NB-1)
    vals=np.empty(Ns)
    for s in range(Ns):
        V=rand_isometry(D,k); X=blocks(V,m,n); vals[s]=C_trace2(X,k)
    mc=vals.mean(); se=vals.std(ddof=1)/np.sqrt(Ns); ex=etrC2(NR,NB,k)
    print(f"({NR},{NB},k={k}) D={D}: MC={mc:.6f}+-{se:.6f}  exact={ex:.6f}  ({(mc-ex)/se:+.1f} se)")

print()
print("=== Table A1: entropies and second moment (k=1) ===")
for (NR,NB) in [(2,2),(3,3),(4,3),(4,4),(5,4)]:
    m,n=specs(NR,NB); D=sum(a*b for a,b in zip(m,n))
    vals_S=np.empty(Ns); vals_Sq=np.empty(Ns); vals_P=np.empty(Ns)
    ds=[1.0,phi]
    for s in range(Ns):
        psi=rng.normal(size=D)+1j*rng.normal(size=D); psi/=np.linalg.norm(psi)
        off=0; Salg=0.0; Sqtr=0.0; pur=0.0
        for (ma,na),da in zip(zip(m,n),ds):
            X=psi[off:off+ma*na].reshape(ma,na); off+=ma*na
            p=np.sum(np.abs(X)**2); rho=X@X.conj().T/p
            ev=np.linalg.eigvalsh(rho); ev=ev[ev>1e-15]
            Salg+= -p*np.sum(ev*np.log(ev)) -p*math.log(p)
            Sqtr+= -p*np.sum(ev*np.log(ev)) -p*math.log(p) + p*math.log(da)
            pur += p*p*np.sum(ev**2)
        vals_S[s]=Salg; vals_Sq[s]=Sqtr; vals_P[s]=pur
    mS,seS=vals_S.mean(),vals_S.std(ddof=1)/np.sqrt(Ns)
    mQ,seQ=vals_Sq.mean(),vals_Sq.std(ddof=1)/np.sqrt(Ns)
    mP,seP=vals_P.mean(),vals_P.std(ddof=1)/np.sqrt(Ns)
    print(f"({NR},{NB}) D={D}: ES_alg={mS:.5f}+-{seS:.5f} (ex {ES_alg(NR,NB):.5f}) | "
          f"ES_qtr={mQ:.5f}+-{seQ:.5f} (ex {ES_qtr(NR,NB):.5f}) | ETr2={mP:.6f}+-{seP:.6f} (ex {ETr2_exact(NR,NB):.6f})")

print()
print("=== Table A2: replica conventions, fib(3,3) D=5, Ns=200000 ===")
m,n=specs(3,3); D=sum(a*b for a,b in zip(m,n)); Ns2=200000
logZ2=np.empty(Ns2); Ss=np.empty(Ns2)
for s in range(Ns2):
    psi=rng.normal(size=D)+1j*rng.normal(size=D); psi/=np.linalg.norm(psi)
    off=0; rhoR=np.zeros((sum(m),sum(m)),dtype=complex)
    o=0
    for ma,na in zip(m,n):
        X=psi[off:off+ma*na].reshape(ma,na); off+=ma*na
        rhoR[o:o+ma,o:o+ma]=X@X.conj().T; o+=ma
    ev=np.linalg.eigvalsh(rhoR); ev=ev[ev>1e-15]
    logZ2[s]=math.log(np.sum(ev**2)); Ss[s]=-np.sum(ev*np.log(ev))
print(f"E[S_alg]     = {Ss.mean():.6f}  (exact {ES_alg(3,3):.6f})")
print(f"E[-log Z_2]  = {(-logZ2).mean():.6f}  [quenched]   std/√N = {(-logZ2).std(ddof=1)/math.sqrt(Ns2):.6f}")
print(f"-log E[Z_2]  = {-math.log(np.exp(logZ2).mean()):.6f}  [annealed]   (exact -log {ETr2_exact(3,3):.6f} = {-math.log(ETr2_exact(3,3)):.6f})")
print(f"S(E rho_R)   = log D = {math.log(D):.6f}  [entropy of averaged state]")
