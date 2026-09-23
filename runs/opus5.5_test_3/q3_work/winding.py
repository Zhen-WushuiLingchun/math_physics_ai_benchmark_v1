import numpy as np
from multiprocessing import Pool
from q_roots import F, w0
def f(args):
    w, eps, L = args
    return F(w, eps, L)
if __name__ == '__main__':
    for eps, L, r in [(1e-3, 60., 0.035), (1e-6, 80., 0.012)]:
        th = np.linspace(0, 2*np.pi, 241)[:-1]
        ws = w0 + r*np.exp(1j*th)
        with Pool(20) as p:
            vals = np.array(p.map(f, [(w, eps, L) for w in ws]))
        ph = np.unwrap(np.angle(np.append(vals, vals[0])))
        wind = (ph[-1]-ph[0])/(2*np.pi)
        jumps = np.max(np.abs(np.diff(ph)))
        print("eps=%g L=%g radius=%.3f around w0: winding number = %.4f  (max phase step %.3f rad)" % (eps, L, r, wind, jumps))
