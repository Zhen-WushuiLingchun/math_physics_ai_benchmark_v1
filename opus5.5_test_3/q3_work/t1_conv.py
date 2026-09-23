import numpy as np, time
from rwcore import *
L=20.0; T=80.0
res={}
for h in [1/32,1/64,1/128,1/256]:
    t0=time.time(); tau,h0,du,rec = td_leapfrog(0.0, L, T, h, born=True)
    res[h]=(tau,h0,du)
    vals=[np.interp(t,tau,h0) for t in (10,20,30,50)]+[np.interp(t,tau,du) for t in (35,40,60)]
    print("LF h=%g t=%.2fs"%(h,time.time()-t0), " ".join("%.10f"%v for v in vals))
for h in [1/16,1/32,1/64]:
    t0=time.time(); tau,h0,du,rec = td_mol(0.0, L, T, h, born=True, rec_dt=1/64)
    vals=[np.interp(t,tau,h0) for t in (10,20,30,50)]+[np.interp(t,tau,du) for t in (35,40,60)]
    print("MOL h=%g t=%.2fs"%(h,time.time()-t0), " ".join("%.10f"%v for v in vals))
