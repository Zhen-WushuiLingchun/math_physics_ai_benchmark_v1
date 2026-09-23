import sympy as sp
z,w=sp.symbols('z omega'); rho=z+2
N=8; c=sp.symbols('c0:%d'%(N+3))
g=sum(c[n]*z**n for n in range(N+3))
ode=sp.expand(rho**2*(rho-2)*sp.diff(g,z,2)+2*rho*sp.diff(g,z)-2*sp.I*w*rho**3*sp.diff(g,z)-(6*rho-6)*g)
ok=True
for n in range(0,N):
    rec=c[n+1]*(n+1)*(4*n+4-16*sp.I*w)+c[n]*(4*n*(n-1)+2*n-24*sp.I*w*n-6)+(c[n-1]*((n-1)*(n-2)-12*sp.I*w*(n-1)-6) if n>=1 else 0)+(c[n-2]*(-2*sp.I*w*(n-2)) if n>=2 else 0)
    if sp.simplify(ode.coeff(z,n)-rec)!=0: ok=False; print(n, sp.simplify(ode.coeff(z,n)-rec))
print("horizon recurrence verified:",ok)
# check the ODE for g itself: psi=e^{-i w x} g, x=rho+2log(rho/2-1)
r=sp.symbols('rho',positive=True); G=sp.Function('G')
x=r+2*sp.log(r/2-1); f=1-2/r; V=f*(6/r**2-6/r**3)
psi=sp.exp(-sp.I*w*x)*G(r)
dx=lambda e: f*sp.diff(e,r)
res=sp.simplify((dx(dx(psi))+(w**2-V)*psi)*sp.exp(sp.I*w*x)*r**3/f)
target=r**2*(r-2)*sp.diff(G(r),r,2)+2*r*sp.diff(G(r),r)-2*sp.I*w*r**3*sp.diff(G(r),r)-(6*r-6)*G(r)
print("g-ODE check:", sp.simplify(sp.expand(res*r - target*1)) , "|", sp.simplify(res/target))
