import numpy as np, time, sys
from multiprocessing import Pool
from freqdomain import Freq

def base(w):
    F = Freq(w)
    return (w, F.I, F.Ain, F.Aout, F.R, F.phiR, F.Aamp, F.h0, F.pm, F.pp)

def bar(args):
    w, epss, L = args
    F = Freq(w)
    B = F.barrier(epss, L)
    return (w, B['beta'], B['dh'], B['h1'], B['b_born'], F.h0)

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'base':
        ws = np.concatenate([np.arange(0.0025, 12.0, 0.0025), np.arange(12.0, 80.0, 0.02)])
        t0 = time.time()
        with Pool(22) as p:
            res = p.map(base, ws, chunksize=8)
        arr = {k: np.array([r[i] for r in res]) for i, k in enumerate(['w','I','Ain','Aout','R','phiR','Aamp','h0','pm','pp'])}
        np.savez('fd_base.npz', **arr)
        print('base done', time.time()-t0)
    elif mode == 'bar':
        Ls = [float(v) for v in sys.argv[2].split(',')]
        epss = [1e-1, 1e-2, 1e-3, 1e-4]
        for L in Ls:
            wmax = 12.0
            dw = min(0.0025, 0.05/L)
            ws = np.arange(dw, wmax, dw)
            t0 = time.time()
            with Pool(22) as p:
                res = p.map(bar, [(w, epss, L) for w in ws], chunksize=8)
            np.savez('fd_bar_L%g.npz' % L, w=ws, eps=np.array(epss),
                     beta=np.array([r[1] for r in res]), dh=np.array([r[2] for r in res]),
                     h1=np.array([r[3] for r in res]), bborn=np.array([r[4] for r in res]),
                     h0=np.array([r[5] for r in res]))
            print('L', L, 'done', time.time()-t0, flush=True)
