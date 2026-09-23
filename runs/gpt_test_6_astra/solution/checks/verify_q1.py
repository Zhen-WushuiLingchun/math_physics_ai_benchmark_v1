"""Exact verification of the n=3 obstruction; no floating-point rank tests.

Run: python solution/checks/verify_q1.py
Conventions and mathematical proof are in solution/q1_solution.tex.
The physical EYM loop amplitudes are literature inputs, not rederived here
from a complete set of Feynman graphs. Their integral limit is checked.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import random
import sys
import sympy as S
from sympy.polys.matrices import DomainMatrix

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'build'
x, z, k = S.symbols('x z k')
y = 1-x
h = x*x-x+1
ORDERS = [(1,'b',2,'a',3),(1,'a',2,'b',3),
          (1,2,'b',3,'a'),(1,2,'a',3,'b'),
          (1,'a',2,3,'b'),(1,'b',2,3,'a')]
checks = []

def check(label, expression):
    assert S.cancel(expression) == 0, (label, S.factor(expression))
    checks.append(label)

class Spinors:
    def __init__(self,L,T): self.L,self.T=L,T
    def a(self,i,j): return S.det(S.Matrix.hstack(self.L[i],self.L[j]))
    def b(self,i,j): return S.det(S.Matrix.hstack(self.T[i],self.T[j]))
    def s(self,i,j): return self.a(i,j)*self.b(j,i)
    def pt(self,o,square=False):
        f=self.b if square else self.a
        return S.prod(f(o[j],o[(j+1)%len(o)]) for j in range(len(o)))
    def plus(self,o):
        a,b=self.a,self.b
        return sum(a(i,j)*b(j,l)*a(l,m)*b(m,i) for i,j,l,m in combinations(o,4))/self.pt(o)
    def plus_alternative(self,o):
        i,j,l,m,n=o; a,b,s=self.a,self.b,self.s
        return (-s(i,j)*s(j,l)-s(i,n)*s(m,n)+a(j,l)*a(m,n)*b(l,m)*b(n,j))/self.pt(o)
    def minus(self,o,neg,alternative=False):
        pos=o.index(neg); i,j,l,m,n=o[pos:]+o[:pos]
        a,b,s=self.a,self.b,self.s
        if not alternative:
            return (-b(j,n)**3/(b(i,j)*b(n,i))
             +a(i,m)**3*b(m,n)*a(l,n)/(a(i,j)*a(j,l)*a(m,n)**2)
             -a(i,l)**3*b(l,j)*a(m,j)/(a(i,n)*a(n,m)*a(l,j)**2))/a(l,m)**2
        num=((s(j,l)+s(l,m)+s(m,n))*b(j,n)**2-b(j,m)*a(m,l)*b(l,n)*b(j,n)
         -b(i,j)*b(i,n)/(a(i,j)*a(i,n))*(a(i,j)**2*a(i,l)**2*b(j,l)/a(j,l)
            +a(i,l)**2*a(i,m)**2*b(l,m)/a(l,m)+a(i,m)**2*a(i,n)**2*b(m,n)/a(m,n)))
        return num/(b(i,j)*a(j,l)*a(l,m)*a(m,n)*b(n,i))
    def eym(self,neg):
        i,j,l={1:(1,2,3),2:(2,3,1),3:(3,1,2)}[neg]
        a,b,s=self.a,self.b,self.s
        return b(j,4)*b(l,4)/(a(j,4)*a(l,4))*(s(i,j)**2+s(i,l)**2)/(2*a(j,l)*b(j,i)*b(l,i))

L={1:S.Matrix([1,0]),2:S.Matrix([0,1]),3:S.Matrix([1,1]),4:S.Matrix([1,z])}
T={1:S.Matrix([-1,-k]),2:S.Matrix([-1,-k*z]),3:S.Matrix([1,0]),4:S.Matrix([0,k])}
# Asymmetric auxiliary spinors. Multiply plus-pair amplitudes by 1/(xy),
# and minus-pair amplitudes by xy, to recover the question's symmetric split.
L['a']=L['b']=L[4]; T['a']=x*T[4]; T['b']=y*T[4]
sp=Spinors(L,T)
tree=[S.factor(1/sp.pt(o)/x/y) for o in ORDERS]
plus=[S.factor(sp.plus(o)/x/y) for o in ORDERS]
minus=[[S.factor(sp.minus(o,n)/x/y) for o in ORDERS] for n in [1,2,3]]
target=[S.factor(sp.eym(n)) for n in [1,2,3]]
for i in range(2):
    for j in range(2): check(f'momentum_{i}{j}',sum(L[n][i]*T[n][j] for n in [1,2,3,4]))
for n,o in enumerate(ORDERS):
    check(f'plus_two_forms_{n}',plus[n]-sp.plus_alternative(o)/x/y)
    for neg in [1,2,3]: check(f'minus_two_forms_{neg}_{n}',minus[neg-1][n]-sp.minus(o,neg,True)/x/y)

q=[x*x*(1-z)+y*y*z,y*y*(1-z)+x*x*z,x*x*z-y*y,y*y*z-x*x,
   y*y*(1-z)-x*x,x*x*(1-z)-y*y]
D=x*z*(x-1)*(z-1)
for j in range(6): check(f'plus_table_{j}',plus[j]-k*k*q[j]/D)
cs=S.symbols('c0:12')
weights=[cs[2*j]*(1-z)+cs[2*j+1]*z for j in range(6)]
tree_target=S.factor(x*y*sp.s(2,4)*tree[0]).subs(k,1)
tree_pol=S.Poly(S.cancel(sum(w*t for w,t in zip(weights,tree))-tree_target).as_numer_denom()[0],z)
At,bt=S.linear_eq_to_matrix(tree_pol.all_coeffs(),cs)
null=At.nullspace()
assert len(null)==8
checks.append('tree_nullity_8')
for ni,v in enumerate(null):
    W=[v[2*j]*sp.s(1,2)+v[2*j+1]*sp.s(2,3) for j in range(6)]
    check(f'nulltree_plus_{ni}',sum(w*t for w,t in zip(W,tree)))
    check(f'nulltree_parity_{ni}',sum(w*x*y/sp.pt(o,True) for w,o in zip(W,ORDERS)))

# Exact coefficient constraints. The denominators are fixed in the manuscript.
den=[D,2*x*z*(x-1)*(z-1)**3,
     2*x*z**3*(x-1)*(z-1)**3,2*x*z**3*(x-1)*(z-1)]
expr=[sum(w*a.subs(k,1) for w,a in zip(weights,plus))]
expr += [sum(w*a.subs(k,1) for w,a in zip(weights,aa))-e.subs(k,1) for aa,e in zip(minus,target)]
polys=[S.Poly(S.cancel(d*f),z) for d,f in zip(den,expr)]
labels=[]; equations=[]
for sec,poly in zip(['plus','minus1','minus2','minus3'],polys):
    for (power,),coeff in poly.terms(): labels.append([sec,power]);equations.append(coeff)
A,bvec=S.linear_eq_to_matrix(equations,cs)
d=x*(x-1)
terms={('plus',2):-2*d*h,('minus1',5):d*(3*h-1),
       ('minus1',4):h*(2*h-1),('minus1',3):h*h,
       ('minus2',4):-h*h,('minus3',3):-h*h}
v=S.Matrix([[terms.get(tuple(label),0) for label in labels]])
for j,e in enumerate(v*A): check(f'left_null_column_{j}',e)
obstruction=d*h*(2*h+1)
check('certificate_nonzero_rhs',(v*bvec)[0]-obstruction)
assert S.factor(obstruction)!=0
checks.append('function_field_inconsistent')
def rank(M): return len(DomainMatrix.from_Matrix(M).convert_to(S.QQ.frac_field(x)).rref()[1])
assert (rank(A),rank(A.row_join(bvec)))==(9,10)
checks.append('loop_ranks_9_10')
Atp=At.col_join(A[:3,:]);btp=bt.col_join(bvec[:3,:])
assert (rank(Atp),rank(Atp.row_join(btp)))==(7,7)
checks.append('tree_plus_ranks_7_7')

# Cyclic average of the standard tree representative.
W=[x*y*sp.s(2,4)/3,0,x*y*sp.s(3,4)/3,0,x*y*sp.s(1,4)/3,0]
delta_plus=S.factor(-sum(w*a for w,a in zip(W,plus)))
delta_minus=[S.factor(e-sum(w*a for w,a in zip(W,aa))) for e,aa in zip(target,minus)]
check('cyclic_tree_calibration',sum(w*t for w,t in zip(W,tree))-tree_target*k)
check('delta_plus_closed',delta_plus+k**3*(2*h-1)*(z*z-z+1)/(3*z*(z-1)))
expected_dm=[k*k*((2*h-1)*(z*z-2*z+2)+4*(z-1))/(6*(z-1)**2),
 k*k*((2*h-1)*(2*z*z-2*z+1)-4*z*(z-1))/(6*z*z*(z-1)**2),
 k*k*((2*h-1)*(z*z+1)-4*z)/(6*z*z)]
for j in range(3): check(f'delta_minus_closed_{j}',delta_minus[j]-expected_dm[j])

# Independent exact checks on generic four-point spinors, not the gauge used
# to construct the coefficient matrix. No fit is performed on these points.
rng=random.Random(20260917)
accepted=0; samples=[]
while accepted<12:
    ll={i:S.Matrix([rng.randint(-7,7),rng.randint(-7,7)]) for i in [1,2,3,4]}
    if any(S.det(S.Matrix.hstack(ll[i],ll[j]))==0 for i,j in combinations([1,2,3,4],2)): continue
    tt={i:S.Matrix([rng.randint(-7,7),rng.randint(-7,7)]) for i in [1,2]}
    rhs=-(ll[1]*tt[1].T+ll[2]*tt[2].T)
    sol=S.Matrix.hstack(ll[3],ll[4]).inv()*rhs
    tt[3]=sol[0,:].T; tt[4]=sol[1,:].T
    if any(S.det(S.Matrix.hstack(tt[i],tt[j]))==0 for i,j in combinations([1,2,3,4],2)): continue
    r,t=[(S.Rational(3,5),S.Rational(4,5)),(S.Rational(5,13),S.Rational(12,13)),(S.Rational(8,17),S.Rational(15,17))][accepted%3]
    xx=r*r; yy=t*t
    ll['a']=r*ll[4];ll['b']=t*ll[4];tt['a']=r*tt[4];tt['b']=t*tt[4]
    ss=Spinors(ll,tt)
    for j,o in enumerate(ORDERS):
        check(f'holdout{accepted}_plus{j}',ss.plus(o)-ss.plus_alternative(o))
        for neg in [1,2,3]:check(f'holdout{accepted}_minus{neg}_{j}',ss.minus(o,neg)-ss.minus(o,neg,True))
    for ni,nv in enumerate(null):
        ww=[nv[2*j]*ss.s(1,2)+nv[2*j+1]*ss.s(2,3) for j in range(6)]
        check(f'holdout{accepted}_tree{ni}',sum(w/ss.pt(o) for w,o in zip(ww,ORDERS)))
        check(f'holdout{accepted}_parity{ni}',sum(w/ss.pt(o,True) for w,o in zip(ww,ORDERS)))
    samples.append({'x':str(xx),'s':str(ss.s(1,2)),'t':str(ss.s(2,3)),
                    'lambda':{str(i):list(map(str,ll[i])) for i in [1,2,3,4]},
                    'tilde':{str(i):list(map(str,tt[i])) for i in [1,2,3,4]}})
    accepted+=1

# Independent evaluation of the dimension-shifted integral constants from
# Feynman-parameter simplex moments; then the EYM Eq. (4.30) finite part.
def moment(powers):
    return S.prod(S.factorial(p) for p in powers)/S.factorial(sum(powers)+len(powers)-1)
s,t,u=S.symbols('s t u')
check('I4_mu4',-moment([0,0,0,0])+S.Rational(1,6))
check('I3_mu4',s*moment([1,0,1])-s/24)
check('I3_mu2',moment([0,0,0])-S.Rational(1,2))
check('I2_mu2',-s*moment([1,1])+s/6)
box=(-S.Rational(3,2)*s*t+S.Rational(1,2)*s*u-S.Rational(3,2)*t*u)*(-S.Rational(1,6))
tri4=-(-t*t*u-2*t*u*u+u**3)/(s*u)*s/24-(t*t*u+t*u*u)/(s*u)*t/24-(2*t*t+4*t*u+u*u)/u*u/24
tri2=-(t**4+3*t**3*u+3*t*t*u*u+t*u**3)/(s*u)/2-(-t**4-2*t**3*u-t*t*u*u+2*t*u**3+u**4)/(s*u)/2+t*u*u/s/2
bubble=t*(t+2*u)/s*(-s/6)+(t*t+2*t*u+2*u*u)/t*(-t/6)-t*(t+2*u)/u*(-u/6)
integral_groups={name:S.factor(value.subs(t,-s-u)) for name,value in [('box_mu4',box),('triangle_mu4',tri4),('triangle_mu2',tri2),('bubble_mu2',bubble)]}
check('eym_minus_integrated',sum(integral_groups.values())-(s*s+u*u)/12)
check('eym_plus_integrated',-S.Rational(1,12)+S.Rational(1,24)+S.Rational(1,24))

OUT.mkdir(exist_ok=True)
data={'field':'QQ(x)','unknowns':list(map(str,cs)),'orders':ORDERS,'row_labels':labels,
      'A':[[str(e) for e in A.row(i)] for i in range(A.rows)],'b':list(map(str,bvec)),
      'certificate':list(map(str,v)),'certificate_rhs':str(S.factor(obstruction)),
      'rank':9,'augmented_rank':10,'tree_nullspace':[list(map(str,n)) for n in null],
      'delta_plus':str(delta_plus),'delta_minus':list(map(str,delta_minus)),
      'integral_groups':{n:str(e) for n,e in integral_groups.items()},'holdouts':samples}
(OUT/'certificate.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
report={'status':'PASS','exact_assertions':len(checks),'holdout_points':accepted,
        'python':sys.version,'sympy':S.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'checks':checks,'limits':['No numerical rank arguments used.','Not an independent full EYM Feynman-diagram calculation.','n=4 and all-multiplicity extension not established.']}
(OUT/'verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='checks'},ensure_ascii=False,indent=2))
