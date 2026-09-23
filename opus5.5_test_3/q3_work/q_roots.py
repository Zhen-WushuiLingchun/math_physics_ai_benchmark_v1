import numpy as np, json, sys, time
import qnm
from scipy.special import lambertw
w0 = 0.3736716844180418-0.0889623156889357j
fo = json.load(open('first_order.json'))
Aout0 = complex(*fo['w0']['Aout']); dA0 = complex(*fo['w0']['dAin'])

def F(w, eps, L):
    J = qnm.jost(w)
    B = qnm.barrier_c(w, eps, L)
    beta = B['b'][0]/B['a'][0]
    return (J['Ain'] - beta*J['Aout'])/J['Aout']

def muller(f, x0, x1, x2, tol=1e-12, maxit=60):
    f0, f1, f2 = f(x0), f(x1), f(x2)
    for it in range(maxit):
        h1 = x1-x0; h2 = x2-x1
        d1 = (f1-f0)/h1; d2 = (f2-f1)/h2
        a = (d2-d1)/(h2+h1); b = a*h2+d2; c = f2
        disc = np.sqrt(b*b-4*a*c)
        den = b+disc if abs(b+disc) > abs(b-disc) else b-disc
        dx = -2*c/den
        x3 = x2+dx
        x0, x1, x2 = x1, x2, x3
        f0, f1, f2 = f1, f2, f(x3)
        if abs(dx) < tol*max(1, abs(x3)):
            return x3, it
    return x2, -1

def born_b(w, L):
    return qnm.barrier_c(w, 0.0, L)['b_born']

if __name__ == '__main__':
    rows = []
    for eps, L in [(1e-3, 20.), (1e-3, 30.), (1e-4, 40.), (1e-4, 50.), (1e-5, 50.), (1e-6, 60.), (1e-6, 70.), (1e-6, 80.), (1e-8, 100.)]:
        t0 = time.time()
        d1 = eps*born_b(w0, L)*Aout0/dA0
        z = -2j*L*d1
        Lam = abs(2*L*d1)
        d_lw = 1j*lambertw(z, 0)/(2*L)
        guess = w0 + (d_lw if Lam < 50 else d1 if Lam < 0.3 else d_lw)
        f = lambda w: F(w, eps, L)
        r, it = muller(f, guess, guess*(1+1e-4), guess*(1-1e-4)+1e-4j)
        rows.append(dict(eps=eps, L=L, Lambda=Lam, d1=[d1.real, d1.imag], dLW=[d_lw.real, d_lw.imag],
                         root=[r.real, r.imag], it=it))
        print("eps=%.0e L=%g Lambda=2L|d1|=%.3g | first-order w=%s | LambertW w=%s | exact root=%s it=%d |F|=%.1e (%.1fs)" % (
            eps, L, Lam, np.round(w0+d1, 6), np.round(w0+d_lw, 6), np.round(r, 6), it, abs(f(r)), time.time()-t0), flush=True)
    json.dump(rows, open('roots_near_w0.json', 'w'), indent=1)
