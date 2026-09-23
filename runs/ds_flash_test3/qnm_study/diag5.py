import mpmath as mp
import qnm_mp as Q
from qnm_mp import *

Vs = V_arrays()
Vi, Vh = Vs

def march_L(w, i0, i1):
    psi = mp.exp(-1j*w*XMIN)
    phi = -1j*w*psi
    w2 = w*w
    for j in range(i0, i1):
        v0, vh, v1 = Vi[j], Vh[j], Vi[j+1]
        k1p, k1q = phi, (v0-w2)*psi
        s2p = psi+DX/2*k1p; s2q = phi+DX/2*k1q
        k2p, k2q = s2q, (vh-w2)*s2p
        s3p = psi+DX/2*k2p; s3q = phi+DX/2*k2q
        k3p, k3q = s3q, (vh-w2)*s3p
        s4p = psi+DX*k3p; s4q = phi+DX*k3q
        k4p, k4q = s4q, (v1-w2)*s4p
        psi = psi + DX/6*(k1p+2*k2p+2*k3p+k4p)
        phi = phi + DX/6*(k1q+2*k2q+2*k3q+k4q)
    return psi, phi

# real frequency: flux conservation |A|^2-|B|^2 = 1
for wr in ['0.4', '0.2', '0.1', '0.05']:
    w = mp.mpf(wr)
    psi, phi = march_L(mp.mpc(w), 0, Q._N)
    flux = 2*mp.im(mp.conj(psi)*phi)   # Im(psi* psi')  -> |A|^2-|B|^2 up to convention
    print(f"omega={wr}: psi(+62)={mp.nstr(psi,8)}  Im(psi* psi')={mp.nstr(flux,8)}")

# same test but with quadrature accuracy check vs dx
print()
print("dx convergence of flux at omega=0.4:")
import copy
for dxv in ['0.05','0.025','0.0125']:
    Q.DX = mp.mpf(dxv); Q._N = int((Q.XMAX-Q.XMIN)/Q.DX)
    Vs2 = V_arrays()
    Vi2, Vh2 = Vs2
    w = mp.mpc('0.4')
    psi = mp.exp(-1j*w*Q.XMIN); phi = -1j*w*psi; w2 = w*w
    for j in range(0, Q._N):
        v0, vh, v1 = Vi2[j], Vh2[j], Vi2[j+1]
        k1p, k1q = phi, (v0-w2)*psi
        s2p = psi+Q.DX/2*k1p; s2q = phi+Q.DX/2*k1q
        k2p, k2q = s2q, (vh-w2)*s2p
        s3p = psi+Q.DX/2*k2p; s3q = phi+Q.DX/2*k2q
        k3p, k3q = s3q, (vh-w2)*s3p
        s4p = psi+Q.DX*k3p; s4q = phi+Q.DX*k3q
        k4p, k4q = s4q, (v1-w2)*s4p
        psi = psi + Q.DX/6*(k1p+2*k2p+2*k3p+k4p)
        phi = phi + Q.DX/6*(k1q+2*k2q+2*k3q+k4q)
    print(f"  dx={dxv}: Im(psi* psi') = {mp.nstr(2*mp.im(mp.conj(psi)*phi),10)}")
