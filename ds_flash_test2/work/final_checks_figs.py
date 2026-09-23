"""Final checks: (i) S_qtr definitional identity; (ii) concentration of F_sdp across samples; (iii) figures."""
import numpy as np, math, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cvxpy as cp

rng = np.random.default_rng(42)
phi=(1+math.sqrt(5))/2; lphi=math.log(phi)
def fib(n):
    p,q=0,1
    for _ in range(n): p,q=q,p+q
    return p
def specs(NR,NB): return [fib(NR-1),fib(NR)], [fib(NB-1),fib(NB)]

# (i) S_qtr identity: S_qtr = -Tr~(rho~ log rho~) computed from definition vs S_alg + p_tau log phi
NR,NB = 4,3
m,n = specs(NR,NB); D=sum(a*b for a,b in zip(m,n))
psi = rng.normal(size=D)+1j*rng.normal(size=D); psi/=np.linalg.norm(psi)
off=0; blocks=[]
for ma,na in zip(m,n):
    blocks.append(psi[off:off+ma*na].reshape(ma,na)); off+=ma*na
ds = [1.0, phi]
# S_alg
mR=sum(m); rhoR=np.zeros((mR,mR),dtype=complex); o=0
for X,ma in zip(blocks,m):
    rhoR[o:o+ma,o:o+ma]=X@X.conj().T; o+=ma
ev=np.linalg.eigvalsh(rhoR); ev=ev[ev>1e-15]; S_alg=-np.sum(ev*np.log(ev))
# S_qtr from definition: -sum_a d_a Tr(Y_a log Y_a), Y_a = X_a X_a^dag/d_a
S_qtr = 0.0
for X,da in zip(blocks,ds):
    Y = X@X.conj().T/da
    ev2 = np.linalg.eigvalsh(Y); ev2=ev2[ev2>1e-15]
    S_qtr += -da*np.sum(ev2*np.log(ev2))
p = [np.sum(np.abs(X)**2) for X in blocks]
print("(i) S_qtr(def) = %.12f ;  S_alg + p_tau ln phi = %.12f ; diff = %.2e"
      % (S_qtr, S_alg + p[1]*math.log(phi), S_qtr-(S_alg+p[1]*math.log(phi))))

# (ii) concentration
def W_from_rho(rho,dR,k):
    W=np.zeros((dR*k,dR*k),dtype=complex)
    for l in range(k):
        for a in range(dR):
            for lp in range(k):
                for ap in range(dR):
                    W[a*k+l,ap*k+lp]=rho[lp*dR+ap,l*dR+a]
    return W
def freq_sdp(rho,dR,k):
    W=W_from_rho(rho,dR,k)
    J=cp.Variable((dR*k,dR*k),hermitian=True)
    Al=[np.kron(np.eye(dR),np.eye(k)[l]) for l in range(k)]
    psum=sum(Al[l]@J@Al[l].conj().T for l in range(k))
    prob=cp.Problem(cp.Maximize(cp.real(cp.trace(W@J))/k),[J>>0,psum==np.eye(dR)])
    prob.solve(solver=cp.SCS,eps=1e-9,max_iters=60000)
    return float(prob.value)
def rand_isometry(D,k):
    A=rng.normal(size=(D,k))+1j*rng.normal(size=(D,k))
    Q,R=np.linalg.qr(A); d=np.diagonal(R).copy(); d/=np.abs(d); return Q*d.conj()
def blocks_of(V,m,n):
    X=[]; off=0
    for ma,na in zip(m,n):
        X.append(V[off:off+ma*na,:].reshape(ma,na,-1).transpose(2,0,1)); off+=ma*na
    return X
def rhoR_of(X,k):
    m=[x.shape[1] for x in X]; dR=sum(m)
    rho=np.zeros((k*dR,k*dR),dtype=complex); off=0
    for x,ma in zip(X,m):
        for q in range(k):
            for qp in range(k):
                rho[q*dR+off:q*dR+off+ma,qp*dR+off:qp*dR+off+ma]=(x[q]@x[qp].conj().T)/k
        off+=ma
    return rho,dR
def trC2(X,k):
    tq=0.0; tb=0.0
    for x in X:
        G=np.einsum('jmn,lmq->jlnq',x.conj(),x)
        tq+=np.sum(np.abs(G)**2)
        P=np.einsum('jmn,jmq->jnq',x.conj(),x)
        tb+=np.einsum('jnq,lqn->jl',P,P).sum()
    return float(np.real((tq-tb/k)/k**2))

print("(ii) concentration over 40 samples:")
for (NR,NB,k) in [(3,3,2),(5,4,2),(6,5,2),(7,4,2)]:
    m,n=specs(NR,NB); D=fib(NR+NB-1)
    Fs=[]; Ts=[]
    for _ in range(40):
        V=rand_isometry(D,k); X=blocks_of(V,m,n)
        rho,dR=rhoR_of(X,k)
        Fs.append(freq_sdp(rho,dR,k)); Ts.append(trC2(X,k))
    Fs=np.array(Fs); Ts=np.array(Ts)
    xi = NR-NB-math.log(k)/lphi
    print(f"({NR},{NB},k={k}) xi={xi:.2f}: F_sdp mean={Fs.mean():.4f} std={Fs.std(ddof=1):.4f} | TrC2 rel-std={Ts.std(ddof=1)/Ts.mean():.3f}")

# (iii) figures
rows=list(csv.DictReader(open("work/results_recovery.csv")))
for r in rows:
    for kk in r:
        try: r[kk]=float(r[kk])
        except: pass
xi=np.array([r['xi'] for r in rows]); y=1-np.array([r['Fsdp'] for r in rows])
delta=np.array([r['NR']-r['NB'] for r in rows]); kk=np.array([r['k'] for r in rows])

fig,ax=plt.subplots(1,2,figsize=(9.5,3.6))
sc=ax[0].scatter(xi,1-y,c=delta,cmap='viridis',s=28,vmin=0,vmax=5)
ax[0].set_xlabel(r'$\Xi=N_R-N_B-\log_\varphi k$'); ax[0].set_ylabel(r'$F_{\rm rec}$ (SDP optimum)')
ax[0].set_ylim(0.15,1.0); ax[0].axvline(0,color='gray',lw=0.7,ls=':')
plt.colorbar(sc,ax=ax[0],label=r'$\delta=N_R-N_B$')
ax[0].set_title('recovery fidelity vs candidate variable')
edges=np.linspace(-3.1,3.7,30)
means=[]; centers=[]
for i in range(len(edges)-1):
    sel=(xi>=edges[i])&(xi<edges[i+1])
    if sel.sum()>=2:
        means.append(y[sel].mean()); centers.append(0.5*(edges[i]+edges[i+1]))
ax[0].plot(centers,1-np.array(means),'r-o',ms=3,lw=1.2,label='bin mean')
ax[0].legend(fontsize=8)
ax2=ax[1]
xs=np.linspace(-3.1,3.7,200)
ax2.semilogy(xi,y,'o',ms=4,color='C0')
ax2.semilogy(xs,0.20*phi**(-xs),'r-',lw=1.2,label=r'$0.20\,\varphi^{-\Xi}$')
ax2.set_xlabel(r'$\Xi$'); ax2.set_ylabel(r'$1-F_{\rm rec}$')
ax2.legend(fontsize=9); ax2.set_title('numerical crossover law')
plt.tight_layout()
plt.savefig('work/fig_collapse.pdf')
print("saved work/fig_collapse.pdf")

# Page-curve figure: exact E S_alg, E S_qtr vs delta at N=24 (even)
from scipy.special import digamma
def page_S(m,n):
    if m>n: m,n=n,m
    return sum(1.0/r for r in range(n+1,m*n+1))-(m-1)/(2.0*n)
def ES(NR,NB,qtr=False):
    m,n=specs(NR,NB); D=sum(a*b for a,b in zip(m,n))
    v=digamma(D+1)
    for ma,na in zip(m,n):
        al=ma*na; v-=(al/D)*digamma(al+1); v+=(al/D)*page_S(ma,na)
    return v + ((m[1]*n[1]/D)*lphi if qtr else 0.0)
N=24
deltas=np.arange(-16,17,2)
Salg=[ES((N+d)//2,(N-d)//2,False) for d in deltas]
Sqtr=[ES((N+d)//2,(N-d)//2,True) for d in deltas]
plt.figure(figsize=(4.6,3.4))
plt.plot(deltas,Salg,'-o',ms=3,label=r'$\mathbb{E}S_{\rm alg}$')
plt.plot(deltas,Sqtr,'-s',ms=3,label=r'$\mathbb{E}S_{\rm qtr}$')
plt.axvline(0,color='gray',ls=':',lw=0.8)
plt.xlabel(r'$\delta=N_R-N_B$'); plt.ylabel('entropy (nats)')
plt.title(r'Page curves ($N=24$)')
plt.legend(fontsize=9); plt.tight_layout(); plt.savefig('work/fig_page.pdf')
print("saved work/fig_page.pdf")
