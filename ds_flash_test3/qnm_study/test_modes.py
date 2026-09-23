import mpmath as mp
import qnm_mp as Q
from qnm_mp import *

def V_l(l):
    # V = (1-2/rho)(l(l+1)/rho^2 - 6/rho^3)
    def f(x):
        r = rho_mp(x)
        return (1-2/r)*(l*(l+1)/r**2 - 6/r**3)
    return f

def V_arrays_generic(Vfun, eps=mp.mpf(0), L=mp.mpf(0)):
    Vi = []; Vh = []
    for j in range(Q._N+1):
        x = Q.XMIN + j*Q.DX
        Vi.append(Vfun(x))
    for j in range(Q._N):
        x = Q.XMIN + (j+mp.mpf('0.5'))*Q.DX
        Vh.append(Vfun(x))
    return Vi, Vh

print("l=2 fundamental (ours): 0.3735791423434 - 0.0890491411376 i")
print("l=2 literature memory  : 0.3736716844180 - 0.0889623152074 i")
print()
for name, l, gr, gi in [("l=2 n=1", 2, '0.3467109961', '-0.2739158157'),
                        ("l=2 n=2", 2, '0.3010531000', '-0.4782760000'),
                        ("l=3 n=0", 3, '0.5994300000', '-0.0927000000'),
                        ("l=4 n=0", 4, '0.8091700000', '-0.0941600000')]:
    Vs = V_arrays_generic(V_l(l))
    w0, it = find_qnm_mp(mp.mpc(gr, gi), Vs)
    print(f"{name}: ours = {mp.nstr(w0, 18)}")
print()
print("memory: l=3 n=0: 0.5994432912 - 0.0927030542 i")
print("memory: l=4 n=0: 0.8091784379 - 0.0941642205 i")
