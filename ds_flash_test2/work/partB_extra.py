""" (C) naive Gaussian-columns model vs Stiefel; (D) classical label leakage. """
import numpy as np, math
rng = np.random.default_rng(1)
phi=(1+math.sqrt(5))/2
def fib(n):
    p,q=0,1
    for _ in range(n): p,q=q,p+q
    return p
def etrC2_formula(NR, NB, k):
    m1, mt = fib(NR-1), fib(NR); n1, nt = fib(NB-1), fib(NB)
    D = m1*n1+mt*nt; A = m1*m1*n1+mt*mt*nt; B = m1*n1*n1+mt*nt*nt
    return (k*k-1)*(D*B-A)/(k*k*D*(D*D-1))

print("=== (C) naive model: k independent normalized (non-orthogonal) Gaussian columns ===")
def mc_naive(NR,NB,k,Ns=4000):
    m1, mt = fib(NR-1), fib(NR); n1, nt = fib(NB-1), fib(NB); D = m1*n1+mt*nt
    vals=[]
    for _ in range(Ns):
        X=[]
        for ma,na in [(m1,n1),(mt,nt)]:
            G = rng.normal(size=(k,ma,na))+1j*rng.normal(size=(k,ma,na))
            for j in range(k): G[j]/=np.linalg.norm(G[j])
            X.append(G/math.sqrt(k))
        tq=0.0; tb=0.0
        for x in X:
            Gj = np.einsum('jmn,lmq->jlnq', x.conj(), x)
            tq += np.sum(np.abs(Gj)**2)
            P = np.einsum('jmn,jmq->jnq', x.conj(), x)
            tb += np.einsum('jnq,lqn->jl', P, P).sum()
        vals.append((tq-tb/k)/k**2)
    return np.mean(vals)
for (NR,NB,k) in [(3,3,2),(3,3,4),(5,4,3)]:
    print(f"({NR},{NB},k={k}): naive-column MC = {mc_naive(NR,NB,k):.5f},  Stiefel exact = {etrC2_formula(NR,NB,k):.5f}")

print()
print("=== (D) classical sector-label leakage:  (1/2)|| P_a^{(j)} - pbar_a ||_1 summed/k ===")
def label_leakage(NR,NB,k,Ns=3000):
    m1, mt = fib(NR-1), fib(NR); n1, nt = fib(NB-1), fib(NB); D = m1*n1+mt*nt
    vals=[]
    for _ in range(Ns):
        v = rng.normal(size=(D,k))+1j*rng.normal(size=(D,k))
        Q,_ = np.linalg.qr(v); dg=np.diagonal(_).copy(); dg/=np.abs(dg); Q*=dg.conj()
        off=0; P=np.zeros((k,2))
        for ia,(ma,na) in enumerate([(m1,n1),(mt,nt)]):
            Xb = Q[off:off+ma*na,:].reshape(ma,na,k).transpose(2,0,1); off+=ma*na
            for j in range(k): P[j,ia]=np.sum(np.abs(Xb[j])**2)
        pbar=P.mean(axis=0)
        vals.append(0.5*np.sum(np.abs(P-pbar)))
    return np.mean(vals)
for (NR,NB,k) in [(3,3,2),(5,4,2),(6,5,3),(8,8,2)]:
    N=NR+NB
    print(f"({NR},{NB},k={k}): leakage={label_leakage(NR,NB,k):.5f},  phi^(-N/2)={phi**(-N/2):.5f}")

print()
print("=== sector weights statistics (mean and fluctuation) ===")
for (NR,NB) in [(3,3),(5,4),(8,8)]:
    m1, mt = fib(NR-1), fib(NR); n1, nt = fib(NB-1), fib(NB); D = m1*n1+mt*nt
    a1, at = m1*n1, mt*nt
    # beta moments
    s2 = a1*at/(D*D*(D+1))
    print(f"({NR},{NB}): E p_tau = {at/D:.6f} (limit {phi**2/(phi**2+1):.6f}), std = {math.sqrt(s2):.6f}, 1/sqrt(alpha_tau) = {1/math.sqrt(at):.6f}")
