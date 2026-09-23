from colK import *
import pickle
rng = random.Random(99)
n = 3
sig = orderings(n)
basis = [(1, 2), (2, 3)]
hels = [(h, hP) for hP in (2, -2) for h in itertools.product([1, -1], repeat=n)]
rows, rhs = [], []
for k in range(8):
    pk = PointK(n, rng)
    for hel, hP in hels:
        M = eym_treeK(pk, hel, hP)
        if M == 0:
            # still include to test zero configs
            pass
        row = []
        for s_ in sig:
            A = ym_col_treeK(pk, s_, hel, 1 if hP > 0 else -1)
            for (i, j) in basis:
                row.append(pk.s(i, j)*A)
        rows.append(row); rhs.append(M)
rank, cons, part, null, piv = solve_affine(rows, rhs)
print('rank', rank, 'consistent', cons, 'pivots', piv)
labels = [(s_, b) for s_ in sig for b in basis]
print('particular solution (RREF):')
for l, v in zip(labels, part):
    if v != 0: print('  ', l, sp.factor(v))
print('nullspace basis:')
for v in null:
    print('  ', [(labels[i], sp.factor(v[i])) for i in range(len(v)) if v[i] != 0])
pickle.dump((sig, basis, part, null), open('tree_n3.pkl', 'wb'))
