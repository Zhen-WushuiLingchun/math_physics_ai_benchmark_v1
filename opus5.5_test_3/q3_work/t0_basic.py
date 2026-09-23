import numpy as np, time
from rwcore import *
print("C =", C_NORM)
print("W norms L1,L2,Linf:", W_norms())
print("E1 =", E1_energy())
print("rho(0)=", rho_of_x(np.array([0.0]))[0], " rho(10)=", rho_of_x(np.array([10.0]))[0], " rho(-1),rho(1)=", rho_of_x(np.array([-1.0,1.0])))
xs=np.linspace(-5,20,200001); v=V0(xs); i=np.argmax(v); print("V0 peak at x=",xs[i]," rho=",rho_of_x(xs[i:i+1])[0]," Vmax=",v[i])
print("check x_of_rho(rho_of_x(x)):", np.max(np.abs(x_of_rho(rho_of_x(np.array([-30.,-1,0,10,100,1e4,1e7])))-np.array([-30.,-1,0,10,100,1e4,1e7]))))
# convergence test of leapfrog, unperturbed and Born term, L=20
for h in [1/16,1/32,1/64]:
    t0=time.time(); tau,h0,du,rec = td_leapfrog(0.0, 20.0, 80.0, h, born=True);
    print(h, "time",time.time()-t0, " h0(10)=",np.interp(10.0,tau,h0), " h0(30)=",np.interp(30.0,tau,h0), " u1dot(40)=",np.interp(40.0,tau,du), " u1dot(35)=",np.interp(35.0,tau,du))
