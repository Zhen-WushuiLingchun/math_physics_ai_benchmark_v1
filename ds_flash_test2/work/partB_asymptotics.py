"""
(A) symbolic asymptotic expansion of  E Tr C^2 = (k^2-1)(D Bsum - Asum)/(k^2 D(D^2-1))
    with Fibonacci data, in powers of phi^{-NR}, phi^{-NB} (parity-resolved).
(C) illustration: naive "independent non-orthogonal Gaussian columns" model.
(D) sector-label leakage: trace distance of classical label cq-state.
"""
import sympy as sp
import numpy as np, math

# ---------- (A) symbolic expansion ----------
phi = (1+sp.sqrt(5))/2
a, b = sp.symbols('a b', positive=True)     # a = phi^NR, b = phi^NB
sR, sB = sp.symbols('s_R s_B')              # (-1)^{NR}, (-1)^{NB}
sq5 = sp.sqrt(5)
# F_n = (phi^n - psi^n)/sq5 ; psi^n = (-1)^n phi^{-n}
F_R   = (a - sR/a)/sq5                      # F_{NR}
F_Rm1 = (a/phi + sR*phi/a)/sq5              # F_{NR-1}
F_B   = (b - sB/b)/sq5
F_Bm1 = (b/phi + sB*phi/b)/sq5
D = sp.simplify(F_R*F_Rm1*0 + (F_Rm1*F_Bm1 + F_R*F_B)*0 + (a*b - sR*sB/(a*b))/sq5*0)  # placeholder
# exact: D = F_{N-1} with N=NR+NB : phi^N = a b
Fn = (a*b + sp.Rational(1,1)/phi*(0) ) # not used
D = sp.simplify((a*b/phi + sp.Rational(1,1)*(a*b)**(-1)*phi*sp.Symbol('s_N'))/sq5*0)  # skip
# Use direct definitions:
Dexpr = sp.simplify(F_Rm1*F_Bm1 + F_R*F_B)
Asum  = sp.simplify(F_Rm1**2*F_Bm1 + F_R**2*F_B)
Bsum  = sp.simplify(F_Rm1*F_Bm1**2 + F_R*F_B**2)
expr = sp.simplify((Dexpr*Bsum - Asum)/(Dexpr*(Dexpr**2-1)))
# expand in eps = 1/(a b) : substitute a = 1/x, b = 1/y and series in x,y
x, y = sp.symbols('x y', positive=True)
expr_xy = sp.simplify(expr.subs({a: 1/x, b: 1/y}))
ser = sp.series(expr_xy, x, 0, 4).removeO()
ser = sp.series(ser, y, 0, 4).removeO()
ser = sp.expand(ser)
print("=== asymptotic expansion of (D Bsum - Asum)/(D(D^2-1)) in phi^{-NR}, phi^{-NB} ===")
print("raised to common denominator form; collecting powers of x=phi^{-NR}, y=phi^{-NB}:")
poly = sp.Poly(sp.expand(sp.simplify(ser * x**1)), x, y)
# collect terms manually
terms = sp.expand(ser)
for monom in sorted(terms.as_ordered_terms(), key=lambda t: sp.degree(t, x)+sp.degree(t, y)):
    print("  ", sp.simplify(monom))
print()
print("leading term coefficient (multiply by (k^2-1)/k^2):")
lead = sp.simplify(sp.expand(ser).coeff(x,1).coeff(y,0))
print("   coefficient of phi^{-NR}:  ", sp.nsimplify(sp.simplify(lead)))
print("   numeric:", float(lead), " vs 2 phi/sqrt5 =", 2*float(phi)/math.sqrt(5))
print("   second leading coefficient (phi^{NR} phi^{-2NB}):",
      sp.simplify(sp.expand(ser).coeff(x,1).coeff(y,2)), "=", float(sp.simplify(sp.expand(ser).coeff(x,1).coeff(y,2))))
print("   coefficient of phi^{-3NR}*phi^{2NB}?:",
      sp.simplify(sp.expand(ser).coeff(x,3).coeff(y,2)), "=", float(sp.simplify(sp.expand(ser).coeff(x,3).coeff(y,2))))

# numerical check of the collected expansion at moderate sizes
def etrC2_formula(NR, NB, k):
    def fib(n):
        A=[1,0]
        if n==0: return 0
        p,q=0,1
        for _ in range(n): p,q=q,p+q
        return p
    m1, mt = fib(NR-1), fib(NR); n1, nt = fib(NB-1), fib(NB)
    D = m1*n1+mt*nt; Asum = m1*m1*n1+mt*mt*nt; Bsum = m1*n1*n1+mt*nt*nt
    return (k*k-1)*(D*Bsum-Asum)/(k*k*D*(D*D-1))

print()
print("=== numerical residual against leading 2phi/sqrt5 * phi^{-NR} ===")
for (NR,NB,k) in [(8,8,2),(10,10,3),(12,12,4),(12,10,2),(14,12,3),(16,16,5)]:
    ex = etrC2_formula(NR,NB,k); pred = (k*k-1)/(k*k)*2*float(phi)/math.sqrt(5)*float(phi)**(-NR)
    print(f"({NR},{NB},k={k}): exact={ex:.3e} lead={pred:.3e} ratio={ex/pred:.6f}")

# ---------- (C) naive Gaussian model (non-orthogonal columns) ----------
print()
print("=== naive model: k independent normalized Gaussian columns (NOT an isometry) ===")
rng = np.random.default_rng(1)
def mc_naive(D, m1, mt, n1, nt, k, Ns=4000):
    vals = []
    for _ in range(Ns):
        X = []
        for ma, na in [(m1,n1),(mt,nt)]:
            G = rng.normal(size=(k,ma,na)) + 1j*rng.normal(size=(k,ma,na))
            for j in range(k): G[j] /= np.linalg.norm(G[j])   # normalize each column
            X.append(G/math.sqrt(k))                          # total trace 1
        # Tr C^2 (same formula as before)
        tr_qb2 = 0.0; tr_b2 = 0.0
        for x in X:
            Gj = np.einsum('jmn,lmq->jlnq', x.conj(), x)
            tr_qb2 += np.sum(np.abs(Gj)**2)
            P = np.einsum('jmn,jmq->jnq', x.conj(), x)
            tr_b2 += np.einsum('jnq,lqn->jl', P, P).sum()
        vals.append((tr_qb2 - tr_b2/k)/k**2)
    return np.mean(vals)
for (NR,NB,k) in [(3,3,2),(3,3,4),(5,4,3)]:
    def fib(n):
        p,q=0,1
        for _ in range(n): p,q=q,p+q
        return p
    m1,mt = fib(NR-1),fib(NR); n1,nt = fib(NB-1),fib(NB)
    D = m1*n1+mt*nt
    print(f"({NR},{NB},k={k}): naive-model MC = {mc_naive(D,m1,mt,n1,nt,k):.5f},   exact (Stiefel) = {etrC2_formula(NR,NB,k):.5f}")

# ---------- (D) classical label leakage ----------
print()
print("=== classical sector label leakage: (1/2)||rho_QA - I/k (x) p||_1 ===")
def label_leakage(NR,NB,k,Ns=2000):
    def fib(n):
        p,q=0,1
        for _ in range(n): p,q=q,p+q
        return p
    m1,mt = fib(NR-1),fib(NR); n1,nt = fib(NB-1),fib(NB); D = m1*n1+mt*nt
    vals = []
    for _ in range(Ns):
        v = rng.normal(size=(D,k)) + 1j*rng.normal(size=(D,k))
        Q,_ = np.linalg.qr(v); dg = np.diagonal(_).copy(); dg/=np.abs(dg); Q*=dg.conj()
        off=0; P = np.zeros((k,2))
        for ma,na in [(m1,n1),(mt,nt)]:
            Xb = Q[off:off+ma*na,:].reshape(ma,na,k).transpose(2,0,1); off+=ma*na
            for j in range(k): P[j, [(m1,n1),(mt,nt)].index((ma,na))] = np.sum(np.abs(Xb[j])**2)
        pbar = P.mean(axis=0)
        vals.append(0.5*np.sum(np.abs(P - pbar)))
    return np.mean(vals)
for (NR,NB,k) in [(3,3,2),(5,4,2),(6,5,3),(8,8,2)]:
    def fib(n):
        p,q=0,1
        for _ in range(n): p,q=q,p+q
        return p
    D = fib(NR+NB-1)
    print(f"({NR},{NB},k={k}): leakage={label_leakage(NR,NB,k):.5f},  phi^(-N/2)={float(phi)**(-(NR+NB)/2):.5f}")
