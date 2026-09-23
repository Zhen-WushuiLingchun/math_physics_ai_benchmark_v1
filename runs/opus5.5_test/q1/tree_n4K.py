import sys, pickle, itertools, time
from colK import *
from multiprocessing import Pool
n = 4
sig = orderings(n)
basis = [(1, 2), (2, 3), (3, 4), (1, 3), (2, 4)]
hels = [(h, 2) for h in itertools.product([1, -1], repeat=n) if 2 <= h.count(-1) <= 3] + \
       [(h, -2) for h in itertools.product([1, -1], repeat=n) if 1 <= h.count(-1) <= 2]
def work(seed):
    rng = random.Random(seed)
    pk = PointK(n, rng)
    rows, rhs = [], []
    for hel, hP in hels:
        M = eym_treeK(pk, hel, hP)
        row = []
        for s_ in sig:
            A = ym_col_treeK(pk, s_, hel, 1 if hP > 0 else -1)
            row += [pk.s(i, j)*A for (i, j) in basis]
        rows.append([str(fromK(v)) for v in row]); rhs.append(str(fromK(M)))
    return rows, rhs
if __name__ == '__main__':
    t0 = time.time()
    with Pool(12) as pool:
        res = pool.map(work, [2000+k for k in range(12)])
    pickle.dump(res, open('tree_n4_rows.pkl', 'wb'))
    print('rows computed', time.time()-t0)
