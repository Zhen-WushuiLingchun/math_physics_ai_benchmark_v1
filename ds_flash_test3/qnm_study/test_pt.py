import mpmath as mp
from qnm_mp import *
import qnm_mp

# override grid to +-40 for PT test (Im w = -0.5)
qnm_mp.XMIN, qnm_mp.XMAX, qnm_mp.DX = mp.mpf(-40), mp.mpf(40), mp.mpf('0.02')

def run(Vfun, wguess, label):
    n = int(mp.floor((qnm_mp.XMAX - qnm_mp.XMIN)/qnm_mp.DX)) + 1
    Vs = [Vfun(qnm_mp.XMIN + j*qnm_mp.DX) for j in range(n)]
    w, it = find_qnm_mp(wguess, Vs)
    return w

# free case: expect W = 2i w exactly, identically zero shift
Vs0 = [mp.mpf(0)]*(int((qnm_mp.XMAX-qnm_mp.XMIN)/qnm_mp.DX)+1)
wtest = mp.mpc('0.3736716844','-0.0889623158')
W = jost_mp(wtest, Vs0)[0]
print("free: W =", mp.nstr(W, 20), "  2i w =", mp.nstr(2j*wtest, 20))

# Poschl-Teller V = 2 sech^2 x : known QNM  w0 = sqrt(1.75) - 0.5 i
def VPT(x):
    c = mp.cosh(x)
    return 2/c**2
wPT = run(VPT, mp.mpc('1.3228756555','-0.5'), 'PT')
print("PT  : found", mp.nstr(wPT, 25))
print("PT  : exact", mp.nstr(mp.sqrt(mp.mpf(7)/4) - mp.mpf('0.5')*1j, 25))
