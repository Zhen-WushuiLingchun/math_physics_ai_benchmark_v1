"""Time-domain batch: waveform difference, onset time, E(T), M(T), Born comparison, two resolutions."""
import numpy as np, json, sys, time
from multiprocessing import Pool
from rwcore import td_leapfrog, td_mol

CASES = [(1e-1, 20.), (1e-2, 20.), (1e-3, 20.), (1e-3, 30.), (1e-4, 50.), (1e-6, 80.), (1e-2, 40.),
         (1e-2, 13.), (1e-3, 60.)]


def run(args):
    eps, L, h = args
    T = 2 * L + 200.0
    t0 = time.time()
    tau, h0, dh, rec = td_leapfrog(eps, L, T, h, born=False)
    _, _, h1, rec1 = td_leapfrog(eps, L, T, h, born=True)
    # exact onset in the discrete scheme: first level where the recorded u is nonzero
    u = rec[1]
    nz = np.nonzero(u)[0]
    onset = nz[0] * h if nz.size else np.inf
    return dict(eps=eps, L=L, h=h, tau=tau, h0=h0, dh=dh, h1=h1, onset=onset, secs=time.time() - t0)


if __name__ == '__main__':
    jobs = [(e, L, h) for (e, L) in CASES for h in (1 / 128, 1 / 256)]
    with Pool(min(len(jobs), 18)) as p:
        res = p.map(run, jobs, chunksize=1)
    out = {}
    for r in res:
        key = 'e%g_L%g_h%g' % (r['eps'], r['L'], 1 / r['h'])
        for k in ('tau', 'h0', 'dh', 'h1'):
            out[key + '_' + k] = r[k]
        out[key + '_onset'] = r['onset']
        print(key, 'onset', r['onset'], '2L-13 =', 2 * r['L'] - 13, 'secs', round(r['secs'], 1), flush=True)
    np.savez_compressed('td_batch.npz', **out)
