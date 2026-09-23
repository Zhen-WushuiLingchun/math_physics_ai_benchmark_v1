import numpy as np, json
import qnm
from freqdomain import source_integral, RHO_O
w0 = 0.3736716844180418-0.0889623156889357j
w1 = 0.34671099687916346-0.2739148752912348j
out={}
for name,wq in [('w0',w0),('w1',w1)]:
    J=qnm.jost(wq); hstep=1e-5
    dA=(qnm.jost(wq+hstep)['Ain']-qnm.jost(wq-hstep)['Ain'])/(2*hstep)
    dA2=(qnm.jost(wq+1j*hstep)['Ain']-qnm.jost(wq-1j*hstep)['Ain'])/(2j*hstep)
    I=source_integral(wq)
    # residue of hat h0 at wq: hat h0 = -i w * (-psi_+(x_o) I / (2 i w A_in)) => Res = -i w * (-psi_+ I/(2 i w A_in'))
    pp=J['pp'][0]
    res_h0 = -1j*wq*(-pp*I/(2j*wq*dA))
    # time-domain amplitude: h(t) contains -i*Res*e^{-i w t} (+ mirror)
    print(name, wq, " Aout=",J['Aout']," dAin=",dA," (check CR)",abs(dA-dA2)," I=",I," Res(h0)=",res_h0, " TD amp:",-1j*res_h0)
    out[name]=dict(w=[wq.real,wq.imag],Aout=[J['Aout'].real,J['Aout'].imag],dAin=[dA.real,dA.imag],res_h0=[res_h0.real,res_h0.imag])
    kap=[]
    for L in [13,15,20,25,30,40,50,60,80,100]:
        B=qnm.barrier_c(wq,0.0,float(L))
        k1=B['b_born']*J['Aout']/dA
        g=-wq.imag
        kap.append((L,k1))
        print("   L=%3d kappa=delta w^(1)/eps=%s |kappa|=%.4e  |kappa| e^{-2 gamma L}=%.5f  unpert_check=%.1e"%(L,np.round(k1,6),abs(k1),abs(k1)*np.exp(-2*g*L),B['unpert_check']))
    out[name]['kappa']=[[L,k.real,k.imag] for L,k in kap]
json.dump(out,open('first_order.json','w'),indent=1)
