import mpmath as mp
def pade_fit(pts, maxdeg=8, tol=1e-25):
    """pts: list of (r, g) real mp. find minimal-degree P/Q with g = P/Q. returns list of candidates."""
    out = []
    for tot in range(0, 2*maxdeg+1):
        for dp in range(0, tot+1):
            dq = tot - dp
            nunk = dp + dq + 2
            if nunk + 2 > len(pts):
                continue
            M = mp.matrix([[r**k for k in range(dp+1)] + [-g*r**k for k in range(dq+1)] for r, g in pts])
            # scale rows
            for i in range(M.rows):
                nr = max(abs(M[i, j]) for j in range(M.cols))
                for j in range(M.cols):
                    M[i, j] /= nr
            U, S, V = mp.svd_r(M)
            smin = S[len(S)-1]
            if smin/S[0] < tol:
                v = V[V.rows-1, :]
                P = [v[k] for k in range(dp+1)]; Q = [v[dp+1+k] for k in range(dq+1)]
                out.append((dp, dq, smin/S[0], P, Q))
        if out:
            return out
    return out
