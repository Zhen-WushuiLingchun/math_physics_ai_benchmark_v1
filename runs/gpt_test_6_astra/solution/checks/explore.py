import sympy as S
from itertools import combinations
x,z,k=S.symbols('x z k', nonzero=True)
y=1-x
L={1:S.Matrix([1,0]),2:S.Matrix([0,1]),3:S.Matrix([1,1]),4:S.Matrix([1,z])}
T={1:S.Matrix([-1,-k]),2:S.Matrix([-1,-z*k]),3:S.Matrix([1,0]),4:S.Matrix([0,k])}
L['a']=L['b']=L[4]
T['a']=x*T[4]; T['b']=y*T[4]
def a(i,j): return S.det(S.Matrix.hstack(L[i],L[j]))
def b(i,j): return S.det(S.Matrix.hstack(T[i],T[j]))
def sij(i,j): return S.factor(a(i,j)*b(j,i))
def pt(order, square=False):
    f=b if square else a
    return S.prod(f(order[i],order[(i+1)%len(order)]) for i in range(len(order)))
orders=[(1,'b',2,'a',3),(1,'a',2,'b',3),(1,2,'b',3,'a'),(1,2,'a',3,'b'),(1,'a',2,3,'b'),(1,'b',2,3,'a')]
def allplus(order):
    num=sum(a(i,j)*b(j,l)*a(l,m)*b(m,i) for i,j,l,m in combinations(order,4))
    return S.factor(num/pt(order)/x/y)
def minus(order,neg=1):
    j=order.index(neg); o=order[j:]+order[:j]
    i,j,l,m,n=o
    val=(-b(j,n)**3/(b(i,j)*b(n,i))
         +a(i,m)**3*b(m,n)*a(l,n)/(a(i,j)*a(j,l)*a(m,n)**2)
         -a(i,l)**3*b(l,j)*a(m,j)/(a(i,n)*a(n,m)*a(l,j)**2))/a(l,m)**2
    return S.factor(val/x/y)
def minus_bdk(order,neg=1):
    j=order.index(neg); o=order[j:]+order[:j]
    i,j,l,m,n=o
    val=((sij(j,l)+sij(l,m)+sij(m,n))*b(j,n)**2-b(j,m)*a(m,l)*b(l,n)*b(j,n)
         -b(i,j)*b(i,n)/(a(i,j)*a(i,n))*(a(i,j)**2*a(i,l)**2*b(j,l)/a(j,l)
            +a(i,l)**2*a(i,m)**2*b(l,m)/a(l,m)+a(i,m)**2*a(i,n)**2*b(m,n)/a(m,n)))
    return S.factor(val/(b(i,j)*a(j,l)*a(l,m)*a(m,n)*b(n,i))/x/y)
def eym(neg=1):
    i,j,l={1:(1,2,3),2:(2,3,1),3:(3,1,2)}[neg]
    return S.factor(b(j,4)*b(l,4)/(a(j,4)*a(l,4))*
                    (sij(i,j)**2+sij(i,l)**2)/(2*a(j,l)*b(j,i)*b(l,i)))
if __name__=='__main__':
    print('invariants',sij(1,2),sij(2,3),sij(1,3))
    tree=[S.factor(1/pt(o)/x/y) for o in orders]
    plus=[allplus(o) for o in orders]
    print('tree',tree)
    print('allplus (coefficient i/(48 pi^2))',plus)
    print('allplus / tree',[S.factor(q/t) for q,t in zip(plus,tree)])
    cs=S.symbols('c0:12')
    weights=[cs[2*i]*sij(1,2)+cs[2*i+1]*sij(2,3) for i in range(6)]
    target=S.factor(x*y*sij(2,4)*tree[0])
    eq=[]
    for expression in [sum(w*t for w,t in zip(weights,tree))-target,sum(w*p for w,p in zip(weights,plus))]:
        num=S.cancel(expression).as_numer_denom()[0]
        eq+=S.Poly(num,z,k).coeffs()
    mat,rhs=S.linear_eq_to_matrix(eq,cs)
    print('tree plus matrix',mat.shape,'rank',mat.rank(),'aug',mat.row_join(rhs).rank())
    print('solution',S.linsolve((mat,rhs),cs))
    for neg in [1,2,3]:
        vals=[minus(o,neg) for o in orders]
        print('minus',neg,vals,'EYM',eym(neg))
        print('BDK crosscheck',[S.factor(v-minus_bdk(o,neg)) for v,o in zip(vals,orders)])
        num=S.cancel(sum(w*v for w,v in zip(weights,vals))-eym(neg)).as_numer_denom()[0]
        eq+=S.Poly(num,z,k).coeffs()
    mat,rhs=S.linear_eq_to_matrix(eq,cs)
    print('full matrix',mat.shape,'rank',mat.rank(),'aug',mat.row_join(rhs).rank())
    print('full solution',S.linsolve((mat,rhs),cs))
