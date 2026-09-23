"""Outgoing/horizon series + Wronskian root for the l=2 Regge-Wheeler QNM.

The pure-exponential start used in an earlier pass is not accurate enough:
truncation of V~6/x^2 moves the root by several hundredths, and the two
radii disagreed. This script starts from a power series and cross-checks
two (radius, step-size, horizon-order) pairs. It also integrates a few
waveforms past the second echo.
"""
from __future__ import annotations

import json
import math

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root

from q3_compute import V_scalar, calibrate_C, evolve


_OUT_FUNCS = None


def _build_outgoing_funcs(N):
    """a_n from the f-equation, solved once. a hand recurrence earlier
    matched a_1 only and was dropped after the residual stopped falling.
    """
    om, z = sp.symbols("omega z")
    aas = list(sp.symbols(f"a1:{N+1}"))
    f = 1 + sum(aas[k] * z ** (k + 1) for k in range(N))
    fp = sp.diff(f, z) * (-z ** 2)
    fpp = sp.diff(fp, z) * (-z ** 2)
    res = (1 - 4 * z + 4 * z ** 2) * fpp
    res += (2 * sp.I * om - 4 * sp.I * om * z + 2 * z ** 2 - 4 * z ** 3) * fp
    res += (-6 * z ** 2 + 18 * z ** 3 - 12 * z ** 4) * f
    poly = sp.Poly(sp.expand(res), z)
    subs = {}
    funcs = []
    for k in range(1, N + 1):
        sol = None
        for p in range(poly.degree() + 1):
            coef = sp.expand(poly.coeff_monomial(z ** p).subs(subs))
            if coef != 0 and aas[k - 1] in coef.free_symbols:
                sol = sp.simplify(sp.solve(coef, aas[k - 1])[0])
                break
        if sol is None:
            raise RuntimeError(f"a_{k} not determined")
        subs[aas[k - 1]] = sol
        funcs.append(sp.lambdify(om, sol, "numpy"))
    return funcs


def outgoing_coeffs(omega, N):
    """f = sum_{n=0}^N a_n rho^{-n}, a_0=1."""
    global _OUT_FUNCS
    if _OUT_FUNCS is None or len(_OUT_FUNCS) < N:
        _OUT_FUNCS = _build_outgoing_funcs(max(N, 12))
    a = np.zeros(N + 1, dtype=np.complex128)
    a[0] = 1.0 + 0.0j
    for n in range(1, N + 1):
        a[n] = complex(_OUT_FUNCS[n - 1](omega))
    return a


def f_ode_residual(omega, rho, a):
    """Residual of the cleared first-order form of the f equation, scaled by rho^0."""
    z = 1.0 / rho
    f = 0j
    fp = 0j  # d/drho
    fpp = 0j
    for n, an in enumerate(a):
        f += an * z ** n
        if n >= 1:
            fp += -n * an * z ** (n + 1)
            fpp += n * (n + 1) * an * z ** (n + 2)
    return (
        (1 - 4 * z + 4 * z ** 2) * fpp
        + (2j * omega - 4j * omega * z + 2 * z ** 2 - 4 * z ** 3) * fp
        + (-6 * z ** 2 + 18 * z ** 3 - 12 * z ** 4) * f
    )


def horizon_coeff_functions(N):
    """Symbolic b_1..b_N for g=1+sum b_k s^k, s=rho-2."""
    s, om = sp.symbols("s omega")
    bs = sp.symbols(f"b1:{N+1}")
    g = 1 + sum(bs[k] * s ** (k + 1) for k in range(N))
    gp = sp.diff(g, s)
    rho = s + 2
    lp = -sp.I * om - 2 * sp.I * om / s
    ps = lp * g + gp
    inner = (s / rho) * ps
    dxx = (s / rho) * (lp * inner + sp.diff(inner, s))
    V = (1 - 2 / rho) * (6 / rho ** 2 - 6 / rho ** 3)
    res = sp.together(dxx + (om ** 2 - V) * g)
    ser = sp.series(res, s, 0, N + 1).removeO()
    poly = sp.Poly(sp.expand(ser * s ** 0 + 0 * s), s) if False else sp.Poly(sp.expand(ser), s)
    subs = {}
    funcs = []
    for k in range(1, N + 1):
        coef = sp.expand(poly.coeff_monomial(s ** k).subs(subs))
        sol = sp.solve(sp.together(coef), bs[k - 1])
        if len(sol) != 1:
            raise RuntimeError(f"b_{k} solve failed: {coef}")
        expr = sp.simplify(sol[0])
        subs[bs[k - 1]] = expr
        funcs.append(sp.lambdify(om, expr, "numpy"))
    b1 = sp.factor(subs[bs[0]])
    return funcs, b1


def eval_g(omega, s, bfuncs):
    g = 1.0 + 0.0j
    gp = 0.0 + 0.0j
    for k, func in enumerate(bfuncs, start=1):
        bk = complex(func(omega))
        g += bk * s ** k
        gp += bk * k * s ** (k - 1)
    return g, gp


def outgoing_state(omega, rho, a):
    f = 0j
    fp = 0j
    for n, an in enumerate(a):
        f += an * rho ** (-n)
        if n:
            fp += -n * an * rho ** (-n - 1)
    pref = np.exp(1j * omega * rho) * (rho - 2) ** (2j * omega)
    dpref = pref * (1j * omega + 2j * omega / (rho - 2))
    dpsi_drho = dpref * f + pref * fp
    dpsi_dx = ((rho - 2) / rho) * dpsi_drho
    x = rho + 2.0 * math.log(rho / 2.0 - 1.0)
    return x, complex(pref * f), complex(dpsi_dx)


def horizon_state(omega, rho, bfuncs):
    s = rho - 2.0
    g, gp = eval_g(omega, s, bfuncs)
    pref = np.exp(-1j * omega * rho) * s ** (-2j * omega)
    dpref_drho = pref * (-1j * omega - 2j * omega / s)
    dpsi_drho = dpref_drho * g + pref * gp
    dpsi_dx = ((rho - 2) / rho) * dpsi_drho
    x = rho + 2.0 * math.log(s / 2.0)
    return x, complex(pref * g), complex(dpsi_dx)


def integrate_to(x0, u0, up0, x1, omega, eps, L, rtol, method):
    w2 = omega * omega

    def f(x, y):
        u = y[0] + 1j * y[1]
        acc = (V_scalar(float(x), eps, L) - w2) * u
        return (y[2], y[3], float(acc.real), float(acc.imag))

    sol = solve_ivp(
        f,
        (x0, x1),
        [u0.real, u0.imag, up0.real, up0.imag],
        method=method,
        rtol=rtol,
        atol=rtol * 1e-3,
    )
    if not sol.success:
        return None
    y = sol.y[:, -1]
    return y[0] + 1j * y[1], y[2] + 1j * y[3]


def rel_w(omega, a, bfuncs, rho_out, rho_hor, x_match, eps, L, rtol, method):
    xR, uR, upR = outgoing_state(omega, rho_out, a)
    xL, uL, upL = horizon_state(omega, rho_hor, bfuncs)
    if not (xL < x_match < xR):
        return None
    left = integrate_to(xL, uL, upL, x_match, omega, eps, L, rtol, method)
    right = integrate_to(xR, uR, upR, x_match, omega, eps, L, rtol, method)
    if left is None or right is None:
        return None
    ul, ulp = left
    ur, urp = right
    W = ul * urp - ulp * ur
    scale = abs(ul) * abs(urp) + abs(ulp) * abs(ur)
    return W / scale


def find_root(guess, aN, bfuncs, cfg, eps=0.0, L=30.0, x_match=3.0):
    a = outgoing_coeffs(guess, aN)  # rebuilt inside obj at the true omega

    def obj(v):
        om = complex(v[0], v[1])
        aa = outgoing_coeffs(om, aN)
        rel = rel_w(om, aa, bfuncs, cfg["rho_out"], cfg["rho_hor"], x_match, eps, L, cfg["rtol"], cfg["method"])
        if rel is None or not np.isfinite(rel.real):
            return [1.0, 1.0]
        return [rel.real, rel.imag]

    sol = root(obj, [guess.real, guess.imag], method="hybr")
    om = complex(sol.x[0], sol.x[1])
    aa = outgoing_coeffs(om, aN)
    rel = rel_w(om, aa, bfuncs, cfg["rho_out"], cfg["rho_hor"], x_match, eps, L, cfg["rtol"], cfg["method"])
    return {
        "success": bool(sol.success) and rel is not None and abs(rel) < 1e-6,
        "omega_r": om.real,
        "omega_i": om.imag,
        "rel": None if rel is None else abs(rel),
        "nfev": int(sol.nfev),
    }, aa


def main():
    print("horizon series...", flush=True)
    N_h = 10
    bfuncs, b1 = horizon_coeff_functions(N_h)
    print("b1", b1, flush=True)
    # outgoing residual check
    om_t = 0.4 - 0.1j
    for rho in (30.0, 60.0):
        aa = outgoing_coeffs(om_t, 12)
        print(f"outgoing residual rho={rho} {abs(f_ode_residual(om_t, rho, aa)):.3e} a1={aa[1]}", flush=True)
    # a1 closed form 3i/omega
    print("a1 vs 3i/omega", aa[1], 3j / om_t, flush=True)

    configs = {
        "A": {"rho_out": 30.0, "rho_hor": 2.2, "rtol": 1e-7, "method": "RK45", "N": 8},
        "B": {"rho_out": 55.0, "rho_hor": 2.08, "rtol": 1e-9, "method": "DOP853", "N": 12},
    }
    # use N_h funcs; extra coeffs beyond the build are absent. Rebuild per config if needed.
    # One series of order 10 is enough for both; truncation is checked by rho_hor.
    roots = {}
    omegaB = None
    for name, cfg in configs.items():
        rec, _ = find_root(0.37 - 0.09j, cfg["N"], bfuncs, cfg)
        roots[name] = rec
        print(f"root {name} {rec}", flush=True)
        if name == "B":
            omegaB = complex(rec["omega_r"], rec["omega_i"])
    results = {"b1": str(b1), "qnm0": roots, "series_checks": {}}
    if omegaB is None or not roots["B"]["success"]:
        json.dump(results, open("q3_qnm.json", "w"), indent=2, default=str)
        print("FAILED root", flush=True)
        return

    # shifts
    shifts = []
    cfg = configs["B"]
    for L, eps in ((14.0, 0.002), (16.0, 0.005), (20.0, 0.004), (24.0, 0.002), (28.0, 0.001), (16.0, 0.02)):
        growth = math.exp(2.0 * abs(omegaB.imag) * L)
        kappa = eps * growth
        deps = min(5e-4, 2e-3 / max(growth, 1.0))
        # derivative of rel at match L+2
        def rel(om, ep):
            aa = outgoing_coeffs(om, cfg["N"])
            return rel_w(om, aa, bfuncs, cfg["rho_out"], cfg["rho_hor"], L + 2.0, ep, L, cfg["rtol"], cfg["method"])

        d_om = 2e-5
        r0 = rel(omegaB, 0.0)
        rp = rel(omegaB + d_om, 0.0)
        rm = rel(omegaB - d_om, 0.0)
        re = rel(omegaB, deps)
        omega1 = -((re - r0) / deps) / ((rp - rm) / (2 * d_om))
        # continuation
        om = omegaB
        ok = True
        last = None
        for k in range(1, 7):
            ep = eps * k / 6
            rec, _ = find_root(om, cfg["N"], bfuncs, cfg, eps=ep, L=L, x_match=L + 2.0)
            if not rec["success"]:
                ok = False
                last = rec["rel"]
                break
            om = complex(rec["omega_r"], rec["omega_i"])
            last = rec["rel"]
        row = {
            "L": L,
            "eps": eps,
            "kappa": kappa,
            "deps": deps,
            "rel0": None if r0 is None else abs(r0),
            "pred_delta_r": float((eps * omega1).real),
            "pred_delta_i": float((eps * omega1).imag),
            "pred_abs": abs(eps * omega1),
            "continuation_ok": ok,
        }
        if ok:
            row["num_delta_r"] = om.real - omegaB.real
            row["num_delta_i"] = om.imag - omegaB.imag
            row["num_abs"] = abs(om - omegaB)
        shifts.append(row)
        print(
            f"L={L} eps={eps} kappa={kappa:.3e} pred={row['pred_abs']:.3e} num={row.get('num_abs')} ok={ok} rel0={row['rel0']}",
            flush=True,
        )
    results["shifts"] = shifts
    results["omega0"] = {"re": omegaB.real, "im": omegaB.imag}
    results["c_crit"] = 1.0 / (2.0 * abs(omegaB.imag))

    print("long waveforms...", flush=True)
    C, _ = calibrate_C(100000)
    long_rows = []
    for dx in (0.05, 0.025):
        times0, h0 = evolve(dx, 130.0, C, 0.0, 16.0, x_max=100.0)
        den = float(np.trapezoid(h0 ** 2, times0))
        tail = float(np.trapezoid(h0[times0 >= 100] ** 2, times0[times0 >= 100]))
        for L, eps in ((16.0, 0.05), (16.0, 0.1), (22.0, 0.05)):
            times, h = evolve(dx, 130.0, C, eps, L, x_max=max(100.0, L + 30.0))
            n = min(len(times), len(times0))
            tt, hh, href = times[:n], h[:n], h0[:n]
            marks = {}
            for Tcut in (15.0, 35.0, 55.0, 80.0, 105.0, 130.0):
                m = tt <= Tcut
                num = float(np.trapezoid((hh[m] - href[m]) ** 2, tt[m]))
                marks[str(Tcut)] = math.sqrt(num / den)
            long_rows.append(
                {
                    "dx": dx,
                    "L": L,
                    "eps": eps,
                    "den": den,
                    "tail_after_100": tail,
                    "E": marks,
                    "E_over_eps": marks["130.0"] / eps,
                }
            )
            print(f"  dx={dx} L={L} eps={eps} den={den:.6e} tail={tail:.3e} E={marks}", flush=True)
    results["long"] = long_rows
    results["C"] = C
    with open("q3_qnm.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print("wrote q3_qnm.json", flush=True)


if __name__ == "__main__":
    main()
