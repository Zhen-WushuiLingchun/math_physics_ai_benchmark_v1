"""Independent reviewer calculation; not code or evidence supplied by a candidate."""
from pathlib import Path
import json
import sys
import time
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(__file__).parent

def bump(x):
    x = np.asarray(x)
    out = np.zeros_like(x, dtype=float)
    m = np.abs(x) < 1
    out[m] = np.exp(1 - 1/(1-x[m]**2))
    return out

def bumpprime(x):
    x = np.asarray(x)
    out = np.zeros_like(x, dtype=float)
    m = np.abs(x) < 1
    out[m] = -2*x[m]*bump(x[m])/(1-x[m]**2)**2
    return out

def potential(x):
    s = np.asarray(x)/2-1
    z = np.where(s <= 1, np.exp(np.minimum(s,1)), np.maximum(s-np.log(np.maximum(s,1)),.5))
    for _ in range(12):
        delta=(z+np.log(z)-s)/(1+1/z)
        z -= delta
    rho = 2*(1+z)
    return (z/(1+z))*(6/rho**2-6/rho**3)

nodes, weights = np.polynomial.legendre.leggauss(512)
energy = float(np.dot(weights,bumpprime(nodes)**2 + potential(nodes)*bump(nodes)**2/2))
C = energy**-0.5

def run(dx, dt_factor=.5, tmax=160, order=2, cases=None, interval=(-200,220)):
    if cases is None:
        cases=[(20.,.01),(20.,.005),(40.,.01)]
    dt = dt_factor*dx
    x = interval[0]+np.arange(round((interval[1]-interval[0])/dx)+1)*dx
    io = round((10-interval[0])/dx)
    steps = round(tmax/dt)
    labels=[{"kind":"background"}]
    for L, eps in cases:
        labels.extend([{"kind":"difference", "L":L, "epsilon":eps},
                       {"kind":"born", "L":L, "epsilon":eps}])
    V = np.broadcast_to(potential(x),(len(labels),x.size)).copy()
    B = np.zeros_like(V)
    for i, case in enumerate(labels):
        if i:
            WL = bump(x-case["L"])
            B[i] = WL*(case["epsilon"] if case["kind"]=="difference" else 1)
            if case["kind"]=="difference":
                V[i] += case["epsilon"]*WL
    def rhs(arr):
        ans = -V*arr-B*arr[0]
        if order==2:
            ans[:,1:-1] += (arr[:,2:]-2*arr[:,1:-1]+arr[:,:-2])/dx**2
            ans[:,[0,-1]]=0
        elif order==4:
            ans[:,2:-2] += (-arr[:,4:]+16*arr[:,3:-1]-30*arr[:,2:-2]+16*arr[:,1:-3]-arr[:,:-4])/(12*dx**2)
            ans[:,:2]=0
            ans[:,-2:]=0
        return ans
    state=np.zeros_like(V)
    vel=np.zeros_like(V)
    state[0]=C*bump(x)
    vel[0]=-C*bumpprime(x)
    acc=rhs(state)
    prev=state-dt*vel+dt**2*acc/2-dt**3*rhs(vel)/6+dt**4*rhs(acc)/24
    signals=np.zeros((steps,len(labels)))
    tic=time.perf_counter()
    for n in range(steps):
        nxt=2*state-prev+dt**2*rhs(state)
        signals[n]=(nxt[:,io]-prev[:,io])/(2*dt)
        prev,state=state,nxt
    times=np.arange(steps)*dt
    bg=signals[:,0]
    bg2=float(np.trapezoid(bg*bg,times))
    results=[]
    for i in range(1,len(labels),2):
        eps=labels[i]["epsilon"]
        L=labels[i]["L"]
        diff=signals[:,i]
        born=signals[:,i+1]
        diff2=float(np.trapezoid(diff*diff,times))
        cross=float(np.trapezoid(diff*bg,times))
        full2=bg2+2*cross+diff2
        # Stable angle formula, avoids subtracting nearly equal quantities.
        cos=(bg2+cross)/np.sqrt(bg2*full2)
        mismatch=(bg2*diff2-cross*cross)/(bg2*full2)/(1+abs(cos))
        rem=diff-eps*born
        results.append({**labels[i],"difference_norm_over_epsilon":np.sqrt(diff2)/eps,
            "E_finite_denominator_over_epsilon":np.sqrt(diff2/bg2)/eps,
            "remainder_norm_over_epsilon_squared":np.sqrt(np.trapezoid(rem*rem,times))/eps**2,
            "mismatch":mismatch,
            "precursor_max_tstar_minus_one":float(np.max(abs(diff[times<2*L-14]))),
            "born_norm":float(np.sqrt(np.trapezoid(born*born,times)))})
    return {"dx":dx,"dt":dt,"order":order,"tmax":tmax,"interval":interval,
            "C":C,"energy_C1":energy,"background_norm":np.sqrt(bg2),"background_norm_squared":bg2,
            "h0_50":float(np.interp(50,times,bg)) if times[-1]>=50 else None,
            "seconds":time.perf_counter()-tic,"cases":results}

if __name__=="__main__":
    if len(sys.argv)>1 and sys.argv[1]=="fourth":
        runs=[run(.025,.4,40,4,[(14,.01),(14,.02)],(-70,90)),
              run(.0125,.4,40,4,[(14,.01),(14,.02)],(-70,90))]
        tag="fourth"
    else:
        runs=[run(.05),run(.025)]
        tag="second"
    obj={"provenance":"Independent reviewer implementation, not candidate-produced evidence", "runs":runs}
    (OUT/("independent_wave_"+tag+".json")).write_text(json.dumps(obj,indent=2),encoding="utf-8")
    print(json.dumps(obj,indent=2))
