import mpmath as mp
import qnm_mp as Q
from qnm_mp import *

def V_l(l):
    def f(x):
        r = rho_mp(x)
        return (1-2/r)*(l*(l+1)/r**2 - 6/r**3)
    return f

def V_arr(l):
    Vi = [V_l(l)(Q.XMIN + j*Q.DX) for j in range(Q._N+1)]
    Vh = [V_l(l)(Q.XMIN + (j+mp.mpf('0.5'))*Q.DX) for j in range(Q._N)]
    return (Vi, Vh)

cases = [("l=2 n=0", 2, '0.3736716844', '-0.0889623158'),
         ("l=3 n=0", 3, '0.5994432912', '-0.0927030542'),
         ("l=4 n=0", 4, '0.8091784379', '-0.0941642205')]
for name, l, gr, gi in cases:
    load_asymp(l, 8)
    Vs = V_arr(l)
    w0, it = find_qnm_mp(mp.mpc(gr, gi), Vs)
    Wp = (rk4_full(w0+mp.mpf('0.001'), Vs, record=False)[0]
          - rk4_full(w0-mp.mpf('0.001'), Vs, record=False)[0])/mp.mpf('0.002')
    print(f"{name}: ours = {mp.nstr(w0, 18)}   |W'| = {mp.nstr(abs(Wp),8)}   ({it} iters)")
print()
print("memory: l=2  0.3736716844180418 - 0.0889623152074178 i")
print("memory: l=3  0.5994432912132190 - 0.0927030542069889 i")
print("memory: l=4  0.8091784379324507 - 0.0941642205054141 i")
