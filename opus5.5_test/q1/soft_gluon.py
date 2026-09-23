from eymloop import *
from eym_formulas import *
from ym_formulas import *
mp.mp.dps = 60
rng = random.Random(12)
# base 4-pt kinematics (1,2,3,P) + soft gluon 4 with lambda_4 fixed, lt_4 = delta*eta ; restore mom. cons. by adjusting lt_1, lt_3
lam, lt = massless_kin(4, rng)          # labels 0,1,2 = gluons 1,2,3 ; 3 = P
lam4 = (rnd(rng), rnd(rng)); eta = (rnd(rng), rnd(rng))
def point(delta):
    L = [lam[0], lam[1], lam[2], lam4, lam[3]]     # gluons 1,2,3,4 then P
    T = [lt[0], lt[1], lt[2], (eta[0]*delta, eta[1]*delta), lt[3]]
    # remove delta*|4>[eta| using |1>,|3>: |4> = c1|1> + c3|3>
    c1 = ang(lam4, lam[2])/ang(lam[0], lam[2]); c3 = ang(lam4, lam[0])/ang(lam[2], lam[0])
    T[0] = (T[0][0] - c1*eta[0]*delta, T[0][1] - c1*eta[1]*delta)
    T[2] = (T[2][0] - c3*eta[0]*delta, T[2][1] - c3*eta[1]*delta)
    assert all(v == 0 for v in vsum([bisp(L[i], T[i]) for i in range(5)], Fr(0)))
    return L, T
L0 = [mpv(l) for l in lam]; T0 = [mpv(t) for t in lt]
for m in [0, 1, 2]:          # minus on gluon 1,2,3 ; gluon 4 plus
    hel = [1, 1, 1, 1]; hel[m] = -1
    for dl in [Fr(1, 10**5), Fr(1, 10**7)]:
        L, T = point(dl)
        gl = [mp_gluon(L[i], T[i], hel[i], rng) for i in range(4)]
        pP, eP = mp_gluon(L[4], T[4], 1, rng)
        M4 = eym_one_loop_split(gl, eP, pP)[0]
        l1, l3, l4 = mpv(L[0]), mpv(L[2]), mpv(L[3])
        S = ang(l3, l1)/(ang(l3, l4)*ang(l4, l1))
        M3 = M1_sm_n3(L0, T0, m)
        print('minus on', m+1, 'delta', dl, ' M4/(S*M3) =', mp.nstr(M4/(S*M3)/mpcf(dl)**-1 if False else M4/(S*M3), 15))
