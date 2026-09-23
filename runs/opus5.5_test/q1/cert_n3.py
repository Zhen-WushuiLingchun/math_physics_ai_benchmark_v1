import pickle, sympy as sp, itertools
from symx import *
from sympy.polys.matrices import DomainMatrix
rows, rhs = pickle.load(open('oneloop_n3_sys.pkl', 'rb'))
rows = [[KX.from_sympy(sp.sympify(v)) for v in r] for r in rows]; rhs = [KX.from_sympy(sp.sympify(v)) for v in rhs]
sect = ['++', '1-', '2-', '3-']
def ranks(idx):
    A = DomainMatrix([rows[i] for i in idx], (len(idx), len(rows[0])), KX)
    b = DomainMatrix([[rhs[i]] for i in idx], (len(idx), 1), KX)
    return A.rank(), A.hstack(b).rank()
# sector-wise
for s in range(4):
    idx = [i for i in range(len(rows)) if i % 4 == s]
    print('sector', sect[s], 'only:', ranks(idx))
idx = [i for i in range(len(rows)) if i % 4 != 0]
print('single-minus sectors only:', ranks(idx))
# minimal subsets
for size in range(1, 8):
    found = None
    for sub in itertools.combinations(range(len(rows)), size):
        r, rb = ranks(list(sub))
        if rb > r:
            found = sub; break
    if found:
        print('minimal inconsistent subset size', size, [(i//4, sect[i % 4]) for i in found]); break
