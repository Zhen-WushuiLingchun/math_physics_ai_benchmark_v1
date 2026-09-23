from explore import *
import json
cs=S.symbols('c0:12')
weights=[cs[2*i]*(1-z)+cs[2*i+1]*z for i in range(6)]
tau,rho=S.symbols('tau rho')
tree=[S.factor(1/pt(o)/x/y) for o in orders]
target=S.factor(x*y*sij(2,4)*tree[0]).subs(k,1)
sectors=[('tree',sum(w*t for w,t in zip(weights,tree))-tau*target),
 ('plus',sum(w*allplus(o).subs(k,1) for w,o in zip(weights,orders)))]
for neg in [1,2,3]:
    sectors.append((f'minus{neg}',sum(w*minus(o,neg).subs(k,1) for w,o in zip(weights,orders))-rho*eym(neg).subs(k,1)))
rows=[]; labels=[]
for label,expr in sectors:
    num,den=S.cancel(expr).as_numer_denom()
    print(label,'denominator',S.factor(den))
    poly=S.Poly(num,z)
    for (power,),coeff in poly.terms():
        rows.append(coeff); labels.append((label,power))
A,bvec=S.linear_eq_to_matrix(rows,cs)
print('shape',A.shape)
print('rows',labels)
for start in [0,4,7]:
    print('start',start,'rank',A[start:,:].rank(),'aug',A[start:,:].row_join(bvec[start:,:]).rank())
for j,v in enumerate(A.T.nullspace()):
    obstruction=S.factor((v.T*bvec)[0])
    if obstruction!=0:
        print('certificate',j,[(labels[i],S.factor(v[i])) for i in range(len(v)) if v[i]!=0], 'rhs',obstruction)
        break
print('tree rows', [str(S.factor(row)) for row in rows[:4]])
