from explore import *
from sympy.polys.matrices import DomainMatrix
cs=S.symbols('c0:12')
w=[cs[2*i]*(1-z)+cs[2*i+1]*z for i in range(6)]
rows=[]; labels=[]
for label, vals, target in [('plus',[allplus(o).subs(k,1) for o in orders],0)]+[(f'minus{n}',[minus(o,n).subs(k,1) for o in orders],eym(n).subs(k,1)) for n in [1,2,3]]:
    num=S.cancel(sum(a*b for a,b in zip(w,vals))-target).as_numer_denom()[0]
    for (p,),c in S.Poly(num,z).terms(): rows.append(c);labels.append((label,p))
A,bb=S.linear_eq_to_matrix(rows,cs)
dm=DomainMatrix.from_Matrix(A.T).convert_to(S.QQ.frac_field(x))
ns=dm.nullspace().to_Matrix()
for ri in range(ns.rows):
    v=ns[ri,:]
    obs=S.factor((v*bb)[0])
    if obs!=0:
        div=obs
        cert=[S.factor(q/div) for q in list(v)]
        print('cert with rhs 1',[(labels[i],q) for i,q in enumerate(cert) if q!=0])
        print('A proof zero',all(S.cancel(q)==0 for q in (S.Matrix([cert])*A)))
        break
print('done')
