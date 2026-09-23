"""OPP integrand reduction in D dims (harmonic transverse basis) + rational-part assembly."""
from loop import *


def lin(coefs, vecs):
    out = (mp.mpc(0),)*4
    for c, v in zip(coefs, vecs):
        out = add(out, scal(c, v))
    return out


class CutGeom:
    def __init__(self, topo, S, rng):
        self.S = S
        k = len(S)
        q = [topo.offset(p) for p in S]
        self.q0 = q[0]
        v = [sub(q[i], q[0]) for i in range(1, k)]
        self.v = v
        m = len(v)
        if m > 0:
            G = mp.matrix(m, m)
            for i in range(m):
                for j in range(m):
                    G[i, j] = dot(v[i], v[j])
            self.Ginv = G**-1
            rhs = [-dot(v[i], v[i])/2 for i in range(m)]
            c = [sum(self.Ginv[i, j]*rhs[j] for j in range(m)) for i in range(m)]
        else:
            self.Ginv, c = None, []
        self.lpar = lin(c, v)
        self.lpar2 = dot(self.lpar, self.lpar)
        basis = []
        while len(basis) < 5 - k:
            r = tuple(mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(4))
            r = self.proj_perp(r)
            for b in basis:
                r = sub(r, scal(dot(r, b)/dot(b, b), b))
            if abs(dot(r, r)) > 1e-3:
                basis.append(r)
        self.basis = basis
        if k == 4:
            self.n = basis[0]
            self.n2 = dot(basis[0], basis[0])
        if k in (3, 2):
            u1, u2 = basis[-2], basis[-1]
            lam = mp.sqrt(-dot(u1, u1)/dot(u2, u2))
            self.na = add(u1, scal(lam, u2))
            self.nb = sub(u1, scal(lam, u2))
            self.nab = dot(self.na, self.nb)
        if k == 2:
            self.n0 = basis[0]
            self.n02 = dot(basis[0], basis[0])

    def proj_perp(self, r):
        m = len(self.v)
        if m == 0:
            return r
        vr = [dot(self.v[j], r) for j in range(m)]
        c = [sum(self.Ginv[i, j]*vr[j] for j in range(m)) for i in range(m)]
        return sub(r, lin(c, self.v))

    def rho(self, mu2):
        return mu2 - self.lpar2

    def point(self, lperp):
        return sub(add(self.lpar, lperp), self.q0)


def gram_degenerate(topo, S):
    q = [topo.offset(p) for p in S]
    v = [sub(q[i], q[0]) for i in range(1, len(S))]
    m = len(v)
    G = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            G[i, j] = dot(v[i], v[j])
    d = mp.det(G)
    scale = 1
    for i in range(m):
        scale *= abs(G[i, i]) + max(abs(G[i, j]) for j in range(m))
    return abs(d) < mp.mpf(10)**(-(mp.mp.dps - 15))*scale


class Reducer:
    def __init__(self, topo, mu2s, seed=1, Nt3=12, Nt2=8, Nz=5):
        self.topo, self.mu2s = topo, mu2s
        self.rng = random.Random(seed)
        self.Nt3, self.Nt2, self.Nz = Nt3, Nt2, Nz
        self.data = {}
        self.bub = {}

    def D(self, l, pr, mu2):
        q = add(l, self.topo.offset(pr))
        return dot(q, q) - mu2

    def residue(self, S, l, imu):
        d = self.data[S]
        g = d['geom']
        k = len(S)
        if k == 5:
            return d['e']
        lp = add(l, g.q0)
        if k == 4:
            return d['de'][imu] + d['do'][imu]*dot(lp, g.n)/g.n2
        if k == 3:
            c = d['c'][imu]
            rho = g.rho(self.mu2s[imu])
            ta = dot(lp, g.na)/g.nab
            tb = 2*dot(lp, g.nb)/rho
            L = (len(c)-1)//2
            tot = c[L]
            for j in range(1, L+1):
                tot += c[L+j]*ta**j + c[L-j]*tb**j
            return tot
        raise ValueError

    def subtracted(self, S, l, mu2, imu):
        F = self.topo.cut_trees(S, l, mu2)
        sS = set(S)
        for S2 in self.data:
            if len(S2) > len(S) and sS <= set(S2):
                den = mp.mpc(1)
                for p in S2:
                    if p not in sS:
                        den *= self.D(l, p, mu2)
                F -= self.residue(S2, l, imu)/den
        return F

    def run(self):
        topo = self.topo
        self.skipped = []
        for S in topo.cuts(5):
            if gram_degenerate(topo, S):
                self.skipped.append(S); continue
            g = CutGeom(topo, S, self.rng)
            mu2 = g.lpar2
            l = g.point((mp.mpc(0),)*4)
            self.data[S] = {'geom': g, 'e': topo.cut_trees(S, l, mu2)}
        for S in topo.cuts(4):
            if gram_degenerate(topo, S):
                self.skipped.append(S); continue
            g = CutGeom(topo, S, self.rng)
            de, do = [], []
            for imu, mu2 in enumerate(self.mu2s):
                z = mp.sqrt(g.rho(mu2)/g.n2)
                Fp = self.subtracted(S, g.point(scal(z, g.n)), mu2, imu)
                Fm = self.subtracted(S, g.point(scal(-z, g.n)), mu2, imu)
                de.append((Fp+Fm)/2)
                do.append((Fp-Fm)/(2*z))
            self.data[S] = {'geom': g, 'de': de, 'do': do}
        Nt = self.Nt3
        R3 = mp.mpc(0.83, 0.21)
        for S in topo.cuts(3):
            if gram_degenerate(topo, S):
                raise ValueError('degenerate triangle %s' % (S,))
            g = CutGeom(topo, S, self.rng)
            cs = []
            for imu, mu2 in enumerate(self.mu2s):
                rho = g.rho(mu2)
                ts = [R3*mp.expjpi(2*mp.mpf(j)/Nt) for j in range(Nt)]
                vals = []
                for t in ts:
                    lperp = add(scal(rho/(2*t*g.nab), g.na), scal(t, g.nb))
                    vals.append(self.subtracted(S, g.point(lperp), mu2, imu))
                L = Nt//2 - 1
                cs.append([sum(vals[j]*ts[j]**(-l_) for j in range(Nt))/Nt for l_ in range(-L, L+1)])
            self.data[S] = {'geom': g, 'c': cs}
        R2 = mp.mpc(0.77, -0.19)
        for S in topo.cuts(2):
            segs = topo.segments(S)
            if any((len(gl) + (1 if hp else 0)) == 1 for gl, hp in segs):
                continue      # massless external bubble: scaleless, zero
            g = CutGeom(topo, S, self.rng)
            mono, hi = [], []
            for imu, mu2 in enumerate(self.mu2s):
                rho = g.rho(mu2)
                zs = [mp.mpc(0.3*(j+1) - 0.8, 0.1*j) for j in range(self.Nz)]
                gz = []
                ts = [R2*mp.expjpi(2*mp.mpf(j)/self.Nt2) for j in range(self.Nt2)]
                for z in zs:
                    ab = (rho - z*z*g.n02)/(2*g.nab)
                    acc = mp.mpc(0)
                    for t in ts:
                        lperp = add(scal(z, g.n0), add(scal(ab/t, g.na), scal(t, g.nb)))
                        acc += self.subtracted(S, g.point(lperp), mu2, imu)
                    gz.append(acc/self.Nt2)
                V = mp.matrix([[z**p for p in range(self.Nz)] for z in zs])
                coef = mp.lu_solve(V, mp.matrix(gz))
                hi.append(coef[self.Nz-1])
                r2 = rho/g.n02
                mono.append(sum(coef[2*kk]*r2**kk/(2*kk+1) for kk in range((self.Nz+1)//2)))
            self.bub[S] = {'geom': g, 'mono': mono, 'zhi': hi}

    def polyfit(self, vals):
        n = len(self.mu2s)
        V = mp.matrix([[m**p for p in range(n)] for m in self.mu2s])
        return list(mp.lu_solve(V, mp.matrix(vals)))

    def Jmu(self, qs, r):
        """rational value of J_n[mu^{2r}] (eps->0), qs = propagator offsets"""
        import itertools, math
        nn = len(qs); k = r + 2 - nn
        if r == 0 or k < 0:
            return mp.mpc(0)
        m = {}
        for i in range(nn):
            for j in range(i+1, nn):
                d = sub(qs[i], qs[j]); m[(i, j)] = dot(d, d)
        # int dF_n Delta^k, Delta = -sum x_i x_j m_ij ; expand multinomially
        pairs = list(m.keys()); tot = mp.mpc(0)
        def rec(idx, left, expo, coef):
            nonlocal tot
            if idx == len(pairs) - 1:
                e = list(expo); i, j = pairs[idx]; e[i] += left; e[j] += left
                c = coef*m[pairs[idx]]**left/math.factorial(left)
                num = 1
                for a in e: num *= math.factorial(a)
                tot += c*mp.mpf(num)/math.factorial(sum(e) + nn - 1)
                return
            for t in range(left+1):
                e = list(expo); i, j = pairs[idx]; e[i] += t; e[j] += t
                rec(idx+1, left-t, e, coef*m[pairs[idx]]**t/math.factorial(t))
        rec(0, k, [0]*nn, mp.mpc(1))
        integ = tot*math.factorial(k)*(-1)**k
        return (-1)**(nn+k+1)*mp.mpf(math.factorial(r-1))/math.factorial(k)*integ

    def contributions_general(self):
        topo = self.topo; out = {}
        for S, d in list(self.data.items()) + list(self.bub.items()):
            k = len(S)
            if k == 5: continue
            if k == 4: vals = d['de']
            elif k == 3:
                L = (len(d['c'][0])-1)//2; vals = [c[L] for c in d['c']]
            else: vals = d['mono']
            p = self.polyfit(vals)
            qs = [topo.offset(pr) for pr in S]
            out[S] = sum(p[r]*self.Jmu(qs, r) for r in range(1, len(p)))
        return out

    def contributions(self):
        """rational contribution of each cut S (box mu^4, triangle mu^2/mu^4, bubble mu^2/mu^4)"""
        topo = self.topo
        out = {}
        for S, d in self.data.items():
            k = len(S)
            if k == 4:
                out[S] = self.polyfit(d['de'])[2]*mp.mpf(-1)/6
            if k == 3:
                L = (len(d['c'][0])-1)//2
                p = self.polyfit([c[L] for c in d['c']])
                out[S] = p[1]/2 + p[2]*sum(dot(K, K) for K in topo.seg_momenta(S))/24
        for S, d in self.bub.items():
            p = self.polyfit(d['mono'])
            K = topo.seg_momenta(S)[0]; K2 = dot(K, K)
            out[S] = p[1]*(-K2/6) + p[2]*(-K2*K2/60)
        return out

    def assemble(self):
        """returns (rational part, diagnostics)"""
        topo = self.topo
        tot = mp.mpc(0)
        diag = {'box_mu0': 0, 'tri_mu0': 0, 'bub_mu0': 0, 'high_mu': 0, 'tri_hiL': 0, 'bub_hiz': 0}
        parts = {'box': mp.mpc(0), 'tri': mp.mpc(0), 'bub': mp.mpc(0)}
        for S, d in self.data.items():
            k = len(S)
            if k == 4:
                p = self.polyfit(d['de'])
                parts['box'] += p[2]*mp.mpf(-1)/6
                diag['box_mu0'] = max(diag['box_mu0'], abs(p[0]))
                diag['high_mu'] = max([diag['high_mu']] + [abs(x) for x in p[3:]])
            if k == 3:
                L = (len(d['c'][0])-1)//2
                p = self.polyfit([c[L] for c in d['c']])
                Ks = topo.seg_momenta(S)
                parts['tri'] += p[1]/2 + p[2]*sum(dot(K, K) for K in Ks)/24
                diag['tri_mu0'] = max(diag['tri_mu0'], abs(p[0]))
                diag['high_mu'] = max([diag['high_mu']] + [abs(x) for x in p[3:]])
                diag['tri_hiL'] = max([diag['tri_hiL']] + [abs(c[0]) + abs(c[-1]) for c in d['c']])
        for S, d in self.bub.items():
            p = self.polyfit(d['mono'])
            K = topo.seg_momenta(S)[0]
            K2 = dot(K, K)
            parts['bub'] += p[1]*(-K2/6) + p[2]*(-K2*K2/60)
            diag['bub_mu0'] = max(diag['bub_mu0'], abs(p[0]))
            diag['high_mu'] = max([diag['high_mu']] + [abs(x) for x in p[3:]])
            diag['bub_hiz'] = max([diag['bub_hiz']] + [abs(x) for x in d['zhi']])
        tot = parts['box'] + parts['tri'] + parts['bub']
        return tot, parts, diag


MU2S = [mp.mpc(1.3, 0.2), mp.mpc(-0.7, 0.9), mp.mpc(2.1, -1.1), mp.mpc(0.45, -0.6), mp.mpc(-1.6, -0.35)]


def one_loop(gluons, grav=None, mu2s=MU2S, seed=1, weight=None):
    topo = Topology(gluons, grav, weight)
    red = Reducer(topo, mu2s, seed=seed)
    red.run()
    return red.assemble()


def eym_one_loop_split(gluons, eP, pP, mu2s=MU2S, seed=1, return_terms=False):
    """M^(1)(1..n;P) = sum_{g=1}^{n} int 1/2 eps_P.(l + K_{1..g}) I_YM(1..g, P, g+1..n)"""
    n = len(gluons)
    tot = mp.mpc(0)
    terms = []
    diags = []
    for g in range(1, n+1):
        objs = gluons[:g] + [(pP, eP)] + gluons[g:]
        K = (mp.mpc(0),)*4
        for (p, e) in gluons[:g]:
            K = add(K, p)
        t, parts, diag = one_loop(objs, mu2s=mu2s, seed=seed, weight=(mp.mpf(1)/2, eP, K))
        tot += t
        terms.append(t)
        diags.append(diag)
    if return_terms:
        return tot, terms, diags
    return tot, diags


def eym_split_by_Pcorner(gluons, eP, pP, mu2s=MU2S, seed=1):
    """returns (total, part from masters with a lone-P corner, rest)"""
    n = len(gluons)
    tot = mp.mpc(0); pc = mp.mpc(0)
    for g in range(1, n+1):
        objs = gluons[:g] + [(pP, eP)] + gluons[g:]
        K = (mp.mpc(0),)*4
        for (p, e) in gluons[:g]:
            K = add(K, p)
        topo = Topology(objs, None, (mp.mpf(1)/2, eP, K))
        red = Reducer(topo, mu2s, seed=seed); red.run()
        for S, v in red.contributions().items():
            tot += v
            segs = topo.segments(S)
            if any(gl == [g] for gl, hp in segs):     # object index g is P in this ordering
                pc += v
    return tot, pc, tot - pc
