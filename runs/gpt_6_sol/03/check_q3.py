"""Reproducible finite-time checks for question 03.

Run from this directory: python check_q3.py
Only files in this directory are written. The finite computational boundaries
are causally disconnected from the detector during the reported time window.
This script does not compute Schwarzschild QNM poles or the infinite-time norm.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from math import sqrt
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import lambertw


def bump(y: np.ndarray | float) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    mask = np.abs(y) < 1.0
    q = y[mask]
    out[mask] = np.exp(1.0 - 1.0 / (1.0 - q * q))
    return out


def bump_prime(y: np.ndarray | float) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    mask = np.abs(y) < 1.0
    q = y[mask]
    out[mask] = -2.0 * q * np.exp(1.0 - 1.0 / (1.0 - q * q)) / (1.0 - q * q) ** 2
    return out


def rho(x: np.ndarray | float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    z = lambertw(np.exp((x - 2.0) / 2.0)).real
    return 2.0 * (1.0 + z)


def v0(x: np.ndarray | float) -> np.ndarray:
    r = rho(x)
    return (1.0 - 2.0 / r) * (6.0 / r**2 - 6.0 / r**3)


def normalization() -> float:
    derivative_term = quad(lambda x: float(bump_prime(x) ** 2), -1, 1,
                           epsabs=2e-13, limit=200)[0]
    potential_term = quad(lambda x: float(v0(x) * bump(x) ** 2), -1, 1,
                          epsabs=2e-13, limit=200)[0]
    return 1.0 / sqrt(derivative_term + 0.5 * potential_term)


@dataclass
class Run:
    t: np.ndarray
    h: np.ndarray
    dx: float
    dt: float
    x_min: float
    x_max: float
    l: float
    eps: tuple[float, float]


def evolve(dx: float, t_max: float, l: float, eps: tuple[float, float],
           c: float) -> Run:
    x_min, x_max = -70.0, 90.0
    n_x = int(round((x_max - x_min) / dx)) + 1
    x = x_min + dx * np.arange(n_x)
    assert abs(x[-1] - x_max) < 1e-10
    i_obs = int(round((10.0 - x_min) / dx))
    assert abs(x[i_obs] - 10.0) < 1e-10
    n_t = int(np.ceil(t_max / (0.4 * dx)))
    dt = t_max / n_t
    assert dt / dx < 1.0

    w_init = bump(x)
    f = c * w_init
    g = -c * bump_prime(x)
    v = v0(x)[None, :] + np.asarray((0.0, *eps))[:, None] * bump(x - l)[None, :]
    u_prev = np.repeat(f[None, :], 3, axis=0)
    lap = np.zeros_like(u_prev)
    lap[:, 2:-2] = (-u_prev[:, 4:] + 16 * u_prev[:, 3:-1]
                    - 30 * u_prev[:, 2:-2] + 16 * u_prev[:, 1:-3]
                    - u_prev[:, :-4]) / (12 * dx**2)
    u_curr = u_prev + dt * g[None, :] + 0.5 * dt**2 * (lap - v * u_prev)
    u_prev[:, (0, -1)] = 0.0
    u_curr[:, (0, -1)] = 0.0

    h = np.zeros((3, n_t + 1))
    h[:, 0] = g[i_obs]
    for step in range(1, n_t):
        lap[:, 2:-2] = (-u_curr[:, 4:] + 16 * u_curr[:, 3:-1]
                        - 30 * u_curr[:, 2:-2] + 16 * u_curr[:, 1:-3]
                        - u_curr[:, :-4]) / (12 * dx**2)
        u_next = 2 * u_curr - u_prev + dt**2 * (lap - v * u_curr)
        u_next[:, (0, -1)] = 0.0
        h[:, step] = (u_next[:, i_obs] - u_prev[:, i_obs]) / (2 * dt)
        u_prev, u_curr = u_curr, u_next
    h[:, -1] = (u_curr[:, i_obs] - u_prev[:, i_obs]) / dt
    return Run(dt * np.arange(n_t + 1), h, dx, dt, x_min, x_max, l, eps)


def diagnostics(run: Run) -> dict:
    t, h = run.t, run.h
    h0 = h[0]
    base = np.trapezoid(h0**2, t)
    onset = 2 * run.l - 13
    before = t <= onset - 0.25
    out: dict = {
        "base_norm_0_T": sqrt(base),
        "onset": onset,
        "max_pre_onset_abs_difference": float(np.max(np.abs(h[1:, before] - h0[None, before]))),
        "max_boundary_return_risk_time": min(2 * abs(run.x_min) + 10 - 1,
                                               2 * run.x_max - 10 - 1),
    }
    for idx, e in enumerate(run.eps, start=1):
        d = h[idx] - h0
        d_norm = sqrt(np.trapezoid(d**2, t))
        f_norm = sqrt(np.trapezoid(h[idx] ** 2, t))
        overlap = np.trapezoid(h[idx] * h0, t)
        out[f"e{e}_difference_norm"] = d_norm
        out[f"e{e}_truncated_E_proxy"] = d_norm / sqrt(base)
        out[f"e{e}_mismatch_T"] = 1 - abs(overlap) / (f_norm * sqrt(base))
    return out


def fourier_bump(omega: complex) -> complex:
    real = quad(lambda y: float(bump(y) * np.exp(2j * omega * y).real),
                -1, 1, epsabs=1e-13, limit=200)[0]
    imag = quad(lambda y: float(bump(y) * np.exp(2j * omega * y).imag),
                -1, 1, epsabs=1e-13, limit=200)[0]
    return real + 1j * imag


def plot_runs(coarse: Run, fine: Run, target: str) -> None:
    os.environ["MPLCONFIGDIR"] = str(Path(__file__).resolve().parent / ".mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8.4, 4.6))
    plt.plot(fine.t, fine.h[0], color="#23395b", label="background $h_0$")
    for i, e in enumerate(fine.eps, 1):
        plt.plot(fine.t, (fine.h[i] - fine.h[0]) / e,
                 label=rf"$(h_{{\epsilon}}-h_0)/\epsilon$, $\epsilon={e:g}$")
    plt.axvline(2 * fine.l - 13, color="black", ls="--", lw=0.9,
                label="causal onset")
    plt.xlim(7, min(35, fine.t[-1]))
    plt.xlabel(r"$\tau$")
    plt.ylabel("signal / first-order response")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(target, dpi=180)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--l", type=float, default=14.0)
    parser.add_argument("--time", type=float, default=40.0)
    parser.add_argument("--coarse", type=float, default=0.025)
    parser.add_argument("--fine", type=float, default=0.0125)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    assert args.l > 12 and args.time < 100
    c = normalization()
    omega_ref = 0.37367168 - 0.08896232j  # literature value; not computed here
    fw = fourier_bump(omega_ref)
    print("C =", c)
    print("reference omega (literature only) =", omega_ref)
    print("W Fourier at 2 omega_ref =", fw)
    print("c_critical =", 1 / (-2 * omega_ref.imag))

    eps = (0.01, 0.02)
    coarse = evolve(args.coarse, args.time, args.l, eps, c)
    fine = evolve(args.fine, args.time, args.l, eps, c)
    for name, run in (("coarse", coarse), ("fine", fine)):
        print(name, "dx", run.dx, "dt", run.dt)
        for key, value in diagnostics(run).items():
            print(" ", key, value)

    t = coarse.t
    h_interp = np.vstack([np.interp(t, fine.t, fine.h[i]) for i in range(3)])
    base_err = sqrt(np.trapezoid((coarse.h[0] - h_interp[0]) ** 2, t))
    response_err = sqrt(np.trapezoid(
        ((coarse.h[1] - coarse.h[0]) - (h_interp[1] - h_interp[0])) ** 2, t))
    response = sqrt(np.trapezoid((h_interp[1] - h_interp[0]) ** 2, t))
    print("cross-grid background relative L2 error", base_err /
          sqrt(np.trapezoid(h_interp[0] ** 2, t)))
    print("cross-grid delta relative L2 error", response_err / response)

    assert diagnostics(fine)["max_pre_onset_abs_difference"] < 1e-5
    assert response_err / response < 0.1
    if args.plot:
        plot_runs(coarse, fine, "waveform_check.png")
        print("wrote waveform_check.png")


if __name__ == "__main__":
    main()
