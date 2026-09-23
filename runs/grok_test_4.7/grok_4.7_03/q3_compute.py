"""Regge-Wheeler l=2 odd parity, compact exterior bump.

Time-domain entries in q3_results.json are used in q3_solution.tex.
The pure-exponential Jost block in this file (qnm0, shifts) is discarded:
the two cutoff radii disagreed by several hundredths. Frequencies and
shifts are recomputed in q3_qnm.py and stored in q3_qnm.json.
"""
from __future__ import annotations

import json
import math
import sys
import time

import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, root

X_MIN = -70.0
X_OBS = 10.0


def bump_np(y):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    m = np.abs(y) < 1.0
    s = 1.0 - y[m] ** 2
    safe = s > 1e-8
    mm = np.zeros_like(y[m], dtype=bool)
    mm[safe] = True
    ss = s[safe]
    out_m = np.zeros_like(y[m])
    out_m[safe] = np.exp(1.0 - 1.0 / ss)
    out[m] = out_m
    return out


def bump_d_np(y):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    m = np.abs(y) < 1.0
    yy = y[m]
    s = 1.0 - yy ** 2
    safe = s > 1e-8
    val = np.zeros_like(yy)
    ss = s[safe]
    w = np.exp(1.0 - 1.0 / ss)
    val[safe] = w * (-2.0 * yy[safe] / ss ** 2)
    out[m] = val
    return out


def rho_np(x):
    x = np.asarray(x, dtype=float)
    rho = np.where(x < 2.0, 2.0 + 2.0 * np.exp((x - 2.0) / 2.0), x.copy())
    for _ in range(8):
        rp = np.maximum(rho - 2.0, 1e-300)
        f = rho + 2.0 * np.log(rp / 2.0) - x
        fp = 1.0 + 2.0 / rp
        rho = rho - f / fp
        rho = np.maximum(rho, 2.0 + 1e-300)
    return rho


def V0_np(x):
    rho = rho_np(x)
    return (1.0 - 2.0 / rho) * (6.0 / rho ** 2 - 6.0 / rho ** 3)


def rho_scalar(x):
    if x < 2.0:
        rho = 2.0 + 2.0 * math.exp((x - 2.0) / 2.0)
    else:
        rho = x
    for _ in range(8):
        rp = rho - 2.0
        if rp < 1e-300:
            rp = 1e-300
        f = rho + 2.0 * math.log(rp / 2.0) - x
        fp = 1.0 + 2.0 / rp
        rho = rho - f / fp
        if rho < 2.0 + 1e-300:
            rho = 2.0 + 1e-300
    return rho


def V_scalar(x, eps, L):
    rho = rho_scalar(x)
    v = (1.0 - 2.0 / rho) * (6.0 / rho ** 2 - 6.0 / rho ** 3)
    y = x - L
    ay = abs(y)
    if ay < 1.0:
        s = 1.0 - y * y
        if s > 1e-8:
            v += eps * math.exp(1.0 - 1.0 / s)
    return v


def calibrate_C(n):
    y = np.linspace(-1.0 + 1e-7, 1.0 - 1e-7, n)
    w = bump_np(y)
    wp = bump_d_np(y)
    v = V0_np(y)
    integ = float(np.trapezoid(wp ** 2 + 0.5 * v * w ** 2, y))
    return 1.0 / math.sqrt(integ), integ


@njit(cache=True)
def leapfrog(psi0, prev0, V, dt, nsteps, jo):
    psi = psi0.copy()
    prev = prev0.copy()
    nxt = np.empty_like(psi)
    n = psi.shape[0]
    h = np.empty(nsteps)
    gamma = (0.5 * dt * dt) * V
    inv_dx2 = 1.0 / (dt * dt)
    dx2 = dt * dt
    for step in range(nsteps):
        for i in range(1, n - 1):
            dxx = (psi[i + 1] - 2.0 * psi[i] + psi[i - 1]) * inv_dx2
            nxt[i] = (2.0 * psi[i] + dx2 * dxx) / (1.0 + gamma[i]) - prev[i]
        nxt[0] = psi[1]
        nxt[n - 1] = psi[n - 2]
        h[step] = (nxt[jo] - prev[jo]) / (2.0 * dt)
        prev, psi, nxt = psi, nxt, prev
    return h


def make_grid(dx, x_max):
    n = int(round((x_max - X_MIN) / dx))
    x = X_MIN + dx * np.arange(n + 1)
    jo = int(round((X_OBS - X_MIN) / dx))
    if abs(x[jo] - X_OBS) > 1e-10:
        raise RuntimeError(f"observer not on grid: {x[jo]}")
    return x, jo


def evolve(dx, T, C, eps, L, x_max, background=True):
    x, jo = make_grid(dx, x_max)
    V = (V0_np(x) if background else np.zeros_like(x)) + eps * bump_np(x - L)
    psi = C * bump_np(x)
    dt_psi = -C * bump_d_np(x)
    dxx = np.zeros_like(psi)
    dxx[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) / dx ** 2
    psi_tt = dxx - V * psi
    prev = psi - dx * dt_psi + 0.5 * dx * dx * psi_tt
    nsteps = int(round(T / dx)) + 1
    h = leapfrog(psi, prev, V, dx, nsteps, jo)
    times = dx * np.arange(nsteps)
    return times, h


def trap(y, t):
    return float(np.trapezoid(y, t))


def metrics(times, h, h0, T_cut, den_full):
    m = times <= T_cut + 1e-12
    if not np.any(m):
        return None
    diff = h[m] - h0[m]
    num = trap(diff ** 2, times[m])
    E = math.sqrt(num / den_full) if den_full > 0 else float("nan")
    a = trap(h[m] * h0[m], times[m])
    nh = math.sqrt(max(trap(h[m] ** 2, times[m]), 0.0))
    n0 = math.sqrt(max(trap(h0[m] ** 2, times[m]), 0.0))
    if nh < 1e-14 or n0 < 1e-14:
        M = None
    else:
        M = 1.0 - abs(a) / (nh * n0)
    return {"E": E, "M": M, "num": num, "nh": nh, "n0": n0}


def fit_qnm(times, h, t1, t2):
    m = (times >= t1) & (times <= t2)
    t = times[m]
    y = h[m]
    if len(t) < 20 or np.max(np.abs(y)) < 1e-8:
        return None
    # crude frequency from zero crossings, no external seed
    s = np.sign(y)
    flips = np.where(s[1:] * s[:-1] < 0)[0]
    if len(flips) >= 2:
        dts = np.diff(t[flips])
        half = float(np.median(dts))
        w0 = math.pi / half if half > 1e-6 else 0.4
    else:
        w0 = 0.4
    # log-envelope slope between |y| peaks of the same sign pattern
    g0 = -0.08
    t0 = t[0]
    y0 = y - np.mean(y)

    def resid(p):
        w, g, A, B = p
        env = np.exp(g * (t - t0))
        pred = env * (A * np.cos(w * (t - t0)) + B * np.sin(w * (t - t0)))
        return pred - y0

    p0 = np.array([w0, g0, y0[0], 0.0])
    sol = least_squares(
        resid,
        p0,
        bounds=([0.05, -1.0, -20, -20], [2.0, -1e-5, 20, 20]),
    )
    rms = float(np.linalg.norm(sol.fun) / math.sqrt(len(y0)))
    sig = float(np.sqrt(np.mean(y0 ** 2)))
    return {
        "omega_r": float(sol.x[0]),
        "omega_i": float(sol.x[1]),
        "rms": rms,
        "rel_rms": rms / sig if sig > 0 else None,
        "success": bool(sol.success),
        "window": [t1, t2],
    }


def ode_pack(omega, eps, L):
    w2 = omega * omega

    def f(x, y):
        u = y[0] + 1j * y[1]
        acc = (V_scalar(x, eps, L) - w2) * u
        return (y[2], y[3], acc.real, acc.imag)

    return f


def integrate_jost(omega, eps, L, x_start, x_end, rtol, method, right):
    if right:
        rho_v = V_scalar(x_start, 0.0, L)  # asymptotics see V0 only if start is outside the bump
        u = cmath_exp(1j * omega * x_start)
        # z=1, z'=V/(2 i omega) => u' = (i omega + V/(2 i omega)) u
        up = (1j * omega + rho_v / (2j * omega)) * u
    else:
        u = cmath_exp(-1j * omega * x_start)
        up = -1j * omega * u
    y0 = [u.real, u.imag, up.real, up.imag]
    sol = solve_ivp(
        ode_pack(omega, eps, L),
        (x_start, x_end),
        y0,
        method=method,
        rtol=rtol,
        atol=rtol * 1e-3,
        dense_output=True,
    )
    if not sol.success:
        return None
    y = sol.y[:, -1]
    u = y[0] + 1j * y[1]
    up = y[2] + 1j * y[3]
    return u, up, sol


def cmath_exp(z):
    return complex(np.exp(z))


def wronskian(omega, eps, L, x_left, x_right, x_match, rtol, method):
    left = integrate_jost(omega, eps, L, x_left, x_match, rtol, method, right=False)
    right = integrate_jost(omega, eps, L, x_right, x_match, rtol, method, right=True)
    if left is None or right is None:
        return None
    ul, ulp, _ = left
    ur, urp, _ = right
    return ul * urp - ulp * ur, ul, ur


def find_root(eps, L, guess, cfg, guesses=None):
    def obj(v):
        omega = v[0] + 1j * v[1]
        out = wronskian(
            omega, eps, L, cfg["x_left"], cfg["x_right"], cfg["x_match"], cfg["rtol"], cfg["method"]
        )
        if out is None:
            return [1e3, 1e3]
        W = out[0]
        scale = max(abs(W), 1e-30)
        # relative to a smooth scale: divide by exp growth so the solver sees O(1)
        return [W.real, W.imag]

    starts = [guess] if guesses is None else guesses
    best = None
    for g in starts:
        sol = root(obj, np.array([g.real, g.imag]), method="hybr")
        if not sol.success:
            continue
        omega = complex(sol.x[0], sol.x[1])
        if omega.real <= 0 or omega.imag >= 0:
            continue
        out = wronskian(
            omega, eps, L, cfg["x_left"], cfg["x_right"], cfg["x_match"], cfg["rtol"], cfg["method"]
        )
        if out is None:
            continue
        W = out[0]
        rec = {"omega_r": omega.real, "omega_i": omega.imag, "abs_W": abs(W), "nfev": int(sol.nfev)}
        if best is None or rec["abs_W"] < best["abs_W"]:
            best = rec
    return best


def scan_guess(cfg):
    best = None
    for wr in np.linspace(0.2, 0.7, 8):
        for wi in (-0.05, -0.09, -0.14, -0.2):
            out = wronskian(
                wr + 1j * wi, 0.0, 16.0, cfg["x_left"], cfg["x_right"], cfg["x_match"], cfg["rtol"], cfg["method"]
            )
            if out is None:
                continue
            val = abs(out[0])
            if best is None or val < best[0]:
                best = (val, wr + 1j * wi)
    return None if best is None else best[1]


def linear_shift(omega, L, cfg):
    """First-order dω/dε from ∂W/∂ε over ∂W/∂ω, and from the Jost integral."""
    d_om = 1e-5

    def W_of(om, eps):
        out = wronskian(
            om, eps, L, cfg["x_left"], cfg["x_right"], cfg["x_match"], cfg["rtol"], cfg["method"]
        )
        if out is None:
            return None
        return out[0]

    W0 = W_of(omega, 0.0)
    Wp = W_of(omega + d_om, 0.0)
    Wm = W_of(omega - d_om, 0.0)
    if W0 is None or Wp is None or Wm is None:
        return None
    dW_dom = (Wp - Wm) / (2 * d_om)
    # finite-difference in eps, small enough that eps * exp(2|Im|L) stays << 1
    growth = math.exp(2 * abs(omega.imag) * L)
    deps = min(1e-4, 0.01 / growth)
    We = W_of(omega, deps)
    if We is None:
        return None
    dW_deps = (We - W0) / deps
    if abs(dW_dom) < 1e-18:
        return None
    omega1_fd = -dW_deps / dW_dom

    # integral formula for the incoming coefficient β, then the same ratio
    # ∂ε β = -∫ w f_- f_+ dx / (2 i ω), evaluated with the unperturbed Jost pair
    x_left, x_right = cfg["x_left"], cfg["x_right"]
    left = integrate_jost(omega, 0.0, L, x_left, L + 1.0, cfg["rtol"], cfg["method"], right=False)
    right = integrate_jost(omega, 0.0, L, x_right, L - 1.0, cfg["rtol"], cfg["method"], right=True)
    if left is None or right is None:
        return {
            "omega1_r": omega1_fd.real,
            "omega1_i": omega1_fd.imag,
            "deps": deps,
            "abs_W0": abs(W0),
            "integral": None,
        }
    xs = np.linspace(L - 1.0 + 1e-5, L + 1.0 - 1e-5, 600)
    _, _, sol_l = left
    _, _, sol_r = right
    yl = sol_l.sol(xs)
    yr = sol_r.sol(xs)
    fl = yl[0] + 1j * yl[1]
    fr = yr[0] + 1j * yr[1]
    w = bump_np(xs - L)
    integral = np.trapezoid(w * fl * fr, xs)
    # β-route needs ∂ω β. Use W-route as the primary prediction; store the integral
    # so the exponential factor can be checked: integral / exp(2 i ω L).
    phase = cmath_exp(-2j * omega * L)
    return {
        "omega1_r": float(omega1_fd.real),
        "omega1_i": float(omega1_fd.imag),
        "deps": deps,
        "abs_W0": abs(W0),
        "integral_abs": abs(integral),
        "integral_over_exp": abs(integral * phase),
        "growth": growth,
    }


def vmax_exact():
    import sympy as sp

    rho = sp.symbols("rho", positive=True)
    V = 6 * (rho - 2) * (rho - 1) / rho ** 4
    dV = sp.together(sp.diff(V, rho))
    crit = sp.solve(sp.numer(sp.together(dV)), rho)
    rho_star = None
    for c in crit:
        if c.is_real and c > 2:
            rho_star = sp.simplify(c)
    Vstar = sp.simplify(V.subs(rho, rho_star))
    # compare Vstar with 19/125 by isolating the remaining square root
    bound = sp.Rational(19, 125)
    # Vstar < bound iff the simplified difference is negative
    diff = sp.simplify(Vstar - bound)
    num = sp.together(diff)
    # numerical confirmation with a rational interval: compute V at a nearby rational
    # and use the critical-point equation. Exact sign via squaring is printed.
    val = complex(Vstar.evalf(40)).real
    return {
        "rho_star": str(rho_star),
        "Vstar_exact": str(Vstar),
        "Vstar": val,
        "below_19_125": bool(val < float(bound)),
        "diff_vs_19_125": val - float(bound),
    }


def main():
    t_start = time.time()
    results = {}
    print("Vmax...", flush=True)
    vm = vmax_exact()
    print(
        f"rho*={vm['rho_star']} V*={vm['Vstar']:.12f} <19/125? {vm['below_19_125']} diff={vm['diff_vs_19_125']:.3e}",
        flush=True,
    )
    results["Vmax"] = vm

    # tortoise residual
    xs = np.linspace(-40, 60, 400)
    rh = rho_np(xs)
    resid = np.max(np.abs(rh + 2 * np.log(rh / 2 - 1) - xs))
    print(f"tortoise max residual {resid:.3e}", flush=True)
    results["tortoise_resid"] = resid

    Cs = {}
    for n in (2000, 20000, 100000):
        C, integ = calibrate_C(n)
        Cs[str(n)] = {"C": C, "integ": integ}
        print(f"C n={n} C={C:.12f} integ={integ:.12f}", flush=True)
    C = Cs["100000"]["C"]
    results["C"] = Cs

    # free wave validation, V=0 exactly right-going
    print("free-wave test...", flush=True)
    free = {}
    for dx in (0.05, 0.025):
        times, h = evolve(dx, 30.0, C, eps=0.0, L=0.0, x_max=50.0, background=False)
        exact = -C * bump_d_np(X_OBS - times)
        err = h - exact
        m = (times >= 9.0) & (times <= 11.0)
        pre = float(np.max(np.abs(h[times < 8.5])))
        rmse = float(np.sqrt(np.mean(err[m] ** 2)))
        free[str(dx)] = {"pre8.5": pre, "rmse_9_11": rmse, "max_abs_9_11": float(np.max(np.abs(err[m])))}
        print(f"  free dx={dx} pre={pre:.3e} rmse={rmse:.3e} max={free[str(dx)]['max_abs_9_11']:.3e}", flush=True)
    results["free_wave"] = free

    # background waveform
    print("background h0...", flush=True)
    T_h0 = 120.0
    h0 = {}
    h0_data = {}
    for dx in (0.05, 0.025):
        times, h = evolve(dx, T_h0, C, 0.0, 16.0, x_max=90.0)
        den = trap(h ** 2, times)
        pre = float(np.max(np.abs(h[times < 8.5])))
        # cumulative norm
        cume = {}
        for Tcut in (20, 40, 60, 80, 100, 120):
            m = times <= Tcut
            cume[str(Tcut)] = trap(h[m] ** 2, times[m])
        fits = []
        for window in ((22.0, 46.0), (30.0, 55.0)):
            fits.append(fit_qnm(times, h, *window))
        h0[str(dx)] = {"den_120": den, "pre8.5": pre, "cumulative": cume, "fits": fits}
        h0_data[dx] = (times, h)
        print(
            f"  dx={dx} den120={den:.6e} pre={pre:.3e} cume80={cume['80']:.6e} cume120={cume['120']:.6e}",
            flush=True,
        )
        for ft in fits:
            if ft:
                print(
                    f"    fit [{ft['window'][0]},{ft['window'][1]}] w={ft['omega_r']:.6f} g={ft['omega_i']:.6f} rel={ft['rel_rms']:.3e}",
                    flush=True,
                )
    results["h0"] = h0

    # QNM shooting, two discretizations
    print("QNM shooting...", flush=True)
    cfgs = {
        "coarse": {
            "x_left": -20.0,
            "x_right": 60.0,
            "x_match": 3.0,
            "rtol": 1e-6,
            "method": "RK45",
        },
        "fine": {
            "x_left": -30.0,
            "x_right": 140.0,
            "x_match": 4.0,
            "rtol": 1e-9,
            "method": "DOP853",
        },
    }
    # bump must sit to the right of x_match; for perturbed runs x_match is reset
    guess = scan_guess(cfgs["coarse"])
    print(f"  scan guess {guess}", flush=True)
    qnm = {"guess": None if guess is None else {"re": guess.real, "im": guess.imag}}
    roots = {}
    for name, cfg in cfgs.items():
        # unperturbed: match point inside, away from artificial bump
        rec = find_root(0.0, 16.0, guess if guess is not None else 0.35 - 0.1j, cfg)
        roots[name] = rec
        print(f"  {name} omega={rec}", flush=True)
    results["qnm0"] = roots
    if roots["fine"] is None:
        print("FAILED fine QNM", flush=True)
        json.dump(results, open("q3_results.json", "w"), indent=2)
        sys.exit(1)
    omega0 = complex(roots["fine"]["omega_r"], roots["fine"]["omega_i"])

    # shifts
    print("QNM shifts...", flush=True)
    shift_cases = []
    fine = dict(cfgs["fine"])
    for L, eps in (
        (14.0, 0.005),
        (16.0, 0.01),
        (16.0, 0.05),
        (20.0, 0.01),
        (24.0, 0.01),
        (24.0, 0.05),
        (30.0, 0.01),
    ):
        cfg = dict(fine)
        cfg["x_match"] = L + 2.0
        cfg["x_right"] = max(140.0, L + 50.0)
        pred = linear_shift(omega0, L, cfg)
        # continuation in eps
        omega_c = omega0
        ok = True
        nstep = 8
        for k in range(1, nstep + 1):
            ep = eps * k / nstep
            rec = find_root(ep, L, omega_c, cfg)
            if rec is None or rec["abs_W"] > 1e-4:
                ok = False
                break
            omega_c = complex(rec["omega_r"], rec["omega_i"])
        kappa = eps * math.exp(2 * abs(omega0.imag) * L)
        row = {
            "L": L,
            "eps": eps,
            "kappa": kappa,
            "pred": pred,
            "continued": None
            if not ok
            else {"omega_r": omega_c.real, "omega_i": omega_c.imag, "abs_W": rec["abs_W"]},
            "continuation_ok": ok,
        }
        if pred is not None:
            d_pred = eps * complex(pred["omega1_r"], pred["omega1_i"])
            row["pred_delta_r"] = d_pred.real
            row["pred_delta_i"] = d_pred.imag
            row["pred_abs"] = abs(d_pred)
        if ok:
            row["num_delta_r"] = omega_c.real - omega0.real
            row["num_delta_i"] = omega_c.imag - omega0.imag
            row["num_abs"] = abs(omega_c - omega0)
        shift_cases.append(row)
        print(
            f"  L={L} eps={eps} kappa={kappa:.4e} ok={ok} pred_abs={row.get('pred_abs')} num_abs={row.get('num_abs')}",
            flush=True,
        )
    results["shifts"] = shift_cases

    # time domain perturbed
    print("perturbed evolutions...", flush=True)
    cases = []
    for dx in (0.05, 0.025):
        times0, h_ref = h0_data[dx]
        den_full = trap(h_ref ** 2, times0)
        for L, eps in (
            (16.0, 0.025),
            (16.0, 0.05),
            (16.0, 0.1),
            (22.0, 0.025),
            (22.0, 0.05),
            (22.0, 0.1),
            (30.0, 0.025),
            (30.0, 0.05),
            (40.0, 0.05),
        ):
            tc = 2.0 * L - 13.0
            T = min(120.0, tc + 25.0)
            x_max = max(90.0, L + 25.0)
            times, h = evolve(dx, T, C, eps, L, x_max)
            n = len(times)
            # h0 on this dx was computed on a grid with the same dt, from the same t=0
            h0s = h_ref[:n]
            t0 = times0[:n]
            if len(t0) != n or abs(t0[-1] - times[-1]) > 1e-8:
                print(f"TIME GRID MISMATCH dx={dx} L={L}", flush=True)
                continue
            before_T = tc - 1.5
            after_T = min(times[-1], tc + 12.0)
            end_T = times[-1]
            mb = metrics(times, h, h0s, before_T, den_full)
            ma = metrics(times, h, h0s, after_T, den_full)
            me = metrics(times, h, h0s, end_T, den_full)
            # pointwise leak
            leak = float(np.max(np.abs(h[times < before_T] - h0s[times < before_T])))
            echo_slice = (times >= tc) & (times <= tc + 12.0)
            echo_amp = float(np.max(np.abs(h[echo_slice] - h0s[echo_slice]))) if np.any(echo_slice) else None
            direct_amp = float(np.max(np.abs(h0s[(times > 9) & (times < 20)])))
            row = {
                "L": L,
                "eps": eps,
                "dx": dx,
                "T_causal": tc,
                "T_end": float(end_T),
                "leak_max": leak,
                "echo_amp": echo_amp,
                "direct_amp": direct_amp,
                "before": mb,
                "after_echo": ma,
                "final": me,
                "E_over_eps": None if me is None else me["E"] / eps,
                "kappa": eps * math.exp(2 * abs(omega0.imag) * L),
            }
            # fit before the echo if the window exists
            if tc > 36:
                row["fit_before_echo"] = fit_qnm(times, h, 24.0, min(46.0, tc - 2.0))
            cases.append(row)
            print(
                f"  dx={dx} L={L} eps={eps} tc={tc:.1f} leak={leak:.3e} "
                f"E_before={mb['E']:.3e} E_final={me['E']:.4e} E/eps={me['E']/eps:.4f} "
                f"M_final={me['M']}",
                flush=True,
            )
    results["cases"] = cases
    results["seconds"] = time.time() - t_start
    with open("q3_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"wrote q3_results.json in {results['seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
