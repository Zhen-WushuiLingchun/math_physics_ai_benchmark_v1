import pickle, sys
from eymloop import *
mp.mp.dps = 50
rng = random.Random(777)
out=[]
for k in range(14):
    pt = Point(3, rng)
    gl = [mp_gluon(pt.lam[i], pt.lt[i], h, rng) for i, h in enumerate((-1,1,1))]
    pP, eP = mp_gluon(pt.lam[3], pt.lt[3], 1, rng)
    tot, dg = eym_one_loop_split(gl, eP, pP)
    out.append({'lam':pt.lam,'lt':pt.lt,'M':{(-1,1,1):tot}})
pickle.dump(out, open('data_n3_sm.pkl','wb'))
print('done')
