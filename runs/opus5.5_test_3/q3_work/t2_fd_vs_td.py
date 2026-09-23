import numpy as np, time
from rwcore import *
from freqdomain import *
L=20.0; T=400.0; h=1/128
t0=time.time()
tau,h0,h1,rec = td_leapfrog(0.0, L, T, h, born=True)
print("TD time",time.time()-t0)
np.savez('td_L20_born.npz',tau=tau,h0=h0,h1=h1)
for w in [0.05,0.2,0.37,0.8,1.5,3.0]:
    t0=time.time()
    F=Freq(w); B=F.barrier(1e-3,L)
    lap0=np.trapezoid(np.exp(1j*w*tau)*h0,tau); lap1=np.trapezoid(np.exp(1j*w*tau)*h1,tau)
    print("w=%.2f  FD h0=%s TD=%s | FD h1=%s TD=%s | beta/eps=%s bBorn=%s  dt=%.2fs"%(w,np.round(F.h0,7),np.round(lap0,7),np.round(B['h1'],7),np.round(lap1,7),np.round(B['beta']/1e-3,6),np.round(B['b_born'],6),time.time()-t0))
