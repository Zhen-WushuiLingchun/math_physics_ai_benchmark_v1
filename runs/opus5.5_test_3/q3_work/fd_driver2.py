import numpy as np, time, sys
from multiprocessing import Pool
from freqdomain import Freq

LS = [13.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0]
EPSS = [1e-1, 1e-2, 1e-3, 1e-4]

def work(w):
    F = Freq(w)
    base = (F.I, F.Ain, F.Aout, F.R, F.phiR, F.Aamp, F.h0, F.pm, F.pp)
    bars = []
    for L in LS:
        B = F.barrier(EPSS, L)
        bars.append((B['beta'], B['dh'], B['h1'], B['b_born']))
    return base, bars

if __name__ == '__main__':
    ws = np.concatenate([np.arange(0.0008, 12.0, 0.0008)])
    t0 = time.time()
    with Pool(23) as p:
        res = p.map(work, ws, chunksize=16)
    out = dict(w=ws, L=np.array(LS), eps=np.array(EPSS))
    for i, k in enumerate(['I','Ain','Aout','R','phiR','Aamp','h0','pm','pp']):
        out[k] = np.array([r[0][i] for r in res])
    out['beta'] = np.array([[r[1][j][0] for j in range(len(LS))] for r in res])   # (nw, nL, neps)
    out['dh'] = np.array([[r[1][j][1] for j in range(len(LS))] for r in res])
    out['h1'] = np.array([[r[1][j][2] for j in range(len(LS))] for r in res])     # (nw, nL)
    out['bborn'] = np.array([[r[1][j][3] for j in range(len(LS))] for r in res])
    np.savez('fd_scan.npz', **out)
    print('done', time.time() - t0, flush=True)
