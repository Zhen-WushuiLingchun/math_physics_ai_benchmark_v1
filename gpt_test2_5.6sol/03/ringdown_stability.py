"""Independent calculations for the perturbed Regge--Wheeler problem.

The time-domain solver uses a centered second-order scheme on a domain whose
boundaries are causally disconnected from the reported observation window.
The frequency-domain solver matches logarithmic derivatives of left/right
Jost solutions.  No QNM is treated as an L2 eigenfunction.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, root
from scipy.special import lambertw


XO = 10.0


def bump(y: np.ndarray | float) -> np.ndarray | float:
    values = np.asarray(y, dtype=float)
    out = np.zeros_like(values)
    mask = np.abs(values) < 1.0
    z = values[mask]
    out[mask] = np.exp(1.0 - 1.0 / (1.0 - z * z))
    return float(out) if out.ndim == 0 else out


def bump_prime(y: np.ndarray | float) -> np.ndarray | float:
    values = np.asarray(y, dtype=float)
    out = np.zeros_like(values)
    mask = np.abs(values) < 1.0
    z = values[mask]
    out[mask] = -2.0 * z * np.exp(1.0 - 1.0 / (1.0 - z * z)) / (1.0 - z * z) ** 2
    return float(out) if out.ndim == 0 else out


def rho_of_x(x: np.ndarray | float) -> np.ndarray | float:
    values = np.asarray(x, dtype=float)
    rho = 2.0 * (1.0 + lambertw(np.exp(values / 2.0 - 1.0)).real)
    return float(rho) if rho.ndim == 0 else rho


def rw_potential(x: np.ndarray | float) -> np.ndarray | float:
    rho = rho_of_x(x)
    return (1.0 - 2.0 / rho) * (6.0 / rho**2 - 6.0 / rho**3)


def rw_potential_complex(z: complex) -> complex:
    """Analytic continuation used only on the exterior-scaled rays."""
    rho = 2.0 * (1.0 + lambertw(np.exp(z / 2.0 - 1.0)))
    return complex((1.0 - 2.0 / rho) * (6.0 / rho**2 - 6.0 / rho**3))


def total_potential(x: np.ndarray | float, epsilon: float, length: float) -> np.ndarray | float:
    return rw_potential(x) + epsilon * bump(np.asarray(x) - length)


def normalization_constant(dx: float = 2.5e-4) -> tuple[float, float]:
    x = np.arange(-1.0, 1.0 + dx / 2.0, dx)
    w = bump(x)
    wp = bump_prime(x)
    raw_energy = np.trapezoid(wp * wp + 0.5 * rw_potential(x) * w * w, x)
    return 1.0 / math.sqrt(raw_energy), float(raw_energy)


def earliest_influence_time(length: float, xo: float = XO) -> float:
    """Support-to-support null travel time for initial [-1,1], bump [L-1,L+1]."""
    return 2.0 * (length - 1.0) - xo - 1.0


@dataclass
class Waveform:
    t: np.ndarray
    h: np.ndarray
    energy0: float
    dx: float
    dt: float
    xmin: float
    xmax: float
    epsilon: float
    length: float


def simulate_waveform(
    epsilon: float,
    length: float,
    *,
    dx: float = 0.04,
    tmax: float = 100.0,
    xmin: float = -140.0,
    xmax: float = 180.0,
    cfl: float = 0.45,
) -> Waveform:
    intervals = int(round((xmax - xmin) / dx))
    dx = (xmax - xmin) / intervals
    x = np.linspace(xmin, xmax, intervals + 1)
    observer_index = int(round((XO - xmin) / dx))
    if abs(x[observer_index] - XO) > 1e-10:
        raise ValueError("Grid must contain the observer x_o=10")

    c_norm, raw_energy = normalization_constant()
    u = c_norm * bump(x)
    pi = -c_norm * bump_prime(x)
    potential = total_potential(x, epsilon, length)

    steps = int(math.ceil(tmax / (cfl * dx)))
    dt = tmax / steps
    lap = np.zeros_like(u)
    lap[1:-1] = (u[2:] - 2.0 * u[1:-1] + u[:-2]) / dx**2
    u_next = u + dt * pi + 0.5 * dt**2 * (lap - potential * u)
    u_next[[0, -1]] = 0.0

    times = np.linspace(0.0, tmax, steps + 1)
    h = np.empty(steps + 1)
    h[0] = pi[observer_index]
    u_prev = u
    u_now = u_next

    for n in range(1, steps):
        lap[1:-1] = (u_now[2:] - 2.0 * u_now[1:-1] + u_now[:-2]) / dx**2
        u_new = 2.0 * u_now - u_prev + dt**2 * (lap - potential * u_now)
        u_new[[0, -1]] = 0.0
        h[n] = (u_new[observer_index] - u_prev[observer_index]) / (2.0 * dt)
        u_prev, u_now = u_now, u_new
    h[-1] = (u_now[observer_index] - u_prev[observer_index]) / dt

    return Waveform(times, h, c_norm * c_norm * raw_energy, dx, dt, xmin, xmax, epsilon, length)


def waveform_metrics(base: Waveform, perturbed: Waveform, cutoff: float) -> dict[str, float | None]:
    if base.t.shape != perturbed.t.shape or not np.allclose(base.t, perturbed.t):
        raise ValueError("Waveforms must share the same time grid")
    mask = base.t <= cutoff + 1e-12
    t = base.t[mask]
    f = base.h[mask]
    g = perturbed.h[mask]
    total_base = float(np.trapezoid(base.h * base.h, base.t))
    diff = float(np.trapezoid((g - f) ** 2, t))
    nf = float(np.trapezoid(f * f, t))
    ng = float(np.trapezoid(g * g, t))
    inner = float(np.trapezoid(g * f, t))
    mismatch = None if nf <= 0.0 or ng <= 0.0 else 1.0 - abs(inner) / math.sqrt(nf * ng)
    return {
        "T": float(cutoff),
        "absolute_error": math.sqrt(max(0.0, diff / total_base)),
        "mismatch": mismatch,
        "base_norm2_T": nf,
        "perturbed_norm2_T": ng,
    }


def wkb_fundamental_guess() -> complex:
    maximum = minimize_scalar(lambda z: -float(rw_potential(z)), bounds=(-5.0, 10.0), method="bounded")
    xp = float(maximum.x)
    step = 2.0e-3
    vp = float(rw_potential(xp))
    vpp = (float(rw_potential(xp + step)) - 2.0 * vp + float(rw_potential(xp - step))) / step**2
    omega_squared = vp - 0.5j * math.sqrt(-2.0 * vpp)
    guess = complex(np.sqrt(omega_squared))
    return guess if guess.imag < 0 else guess.conjugate()


def _jost_log_derivative(
    omega: complex,
    epsilon: float,
    length: float,
    *,
    side: str,
    xmatch: float,
    ecs_radius: float,
    ecs_length: float,
    ecs_angle: float,
    rtol: float,
    atol: float,
) -> complex:
    phase = np.exp(1j * ecs_angle)
    if side == "left":
        real_edge = -ecs_radius
        contour_sign = -phase
    elif side == "right":
        real_edge = ecs_radius
        contour_sign = phase
    else:
        raise ValueError("side must be left or right")

    # Continue rho from the real axis along the contour instead of repeatedly
    # selecting a Lambert-W branch at complex x.
    def rho_ode(_s: float, rho: np.ndarray) -> np.ndarray:
        return np.asarray([(1.0 - 2.0 / rho[0]) * contour_sign], dtype=complex)

    rho_path = solve_ivp(
        rho_ode,
        (0.0, ecs_length),
        np.asarray([complex(rho_of_x(real_edge))], dtype=complex),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    if not rho_path.success:
        raise RuntimeError(rho_path.message)

    def potential_from_rho(rho: complex) -> complex:
        return (1.0 - 2.0 / rho) * (6.0 / rho**2 - 6.0 / rho**3)

    far_potential = potential_from_rho(complex(rho_path.y[0, -1]))
    if side == "left":
        initial = -1j * omega
    else:
        # First WKB correction for the inverse-square far-zone tail.
        initial = 1j * omega - 1j * far_potential / (2.0 * omega)

    def contour_ode(s: float, y: np.ndarray) -> np.ndarray:
        rho = complex(rho_path.sol(s)[0])
        potential = potential_from_rho(rho)
        return np.asarray([(potential - omega * omega - y[0] * y[0]) * contour_sign], dtype=complex)

    contour = solve_ivp(
        contour_ode,
        (ecs_length, 0.0),
        np.asarray([initial], dtype=complex),
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not contour.success:
        raise RuntimeError(contour.message)

    def real_ode(x: float, y: np.ndarray) -> np.ndarray:
        potential = float(total_potential(x, epsilon, length))
        return np.asarray([potential - omega * omega - y[0] * y[0]], dtype=complex)

    solution = solve_ivp(
        real_ode,
        (real_edge, xmatch),
        np.asarray([contour.y[0, -1]], dtype=complex),
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return complex(solution.y[0, -1])


def resonance_mismatch(
    omega: complex,
    epsilon: float = 0.0,
    length: float = 20.0,
    *,
    xmatch: float = 0.0,
    ecs_radius: float = 28.0,
    ecs_length: float = 70.0,
    ecs_angle: float = 0.48,
    rtol: float = 2.0e-10,
    atol: float = 2.0e-11,
) -> complex:
    if ecs_radius <= length + 1.0 and epsilon != 0.0:
        raise ValueError("Exterior scaling must begin to the right of the added bump")
    left = _jost_log_derivative(
        omega, epsilon, length, side="left", xmatch=xmatch, ecs_radius=ecs_radius,
        ecs_length=ecs_length, ecs_angle=ecs_angle, rtol=rtol, atol=atol
    )
    right = _jost_log_derivative(
        omega, epsilon, length, side="right", xmatch=xmatch, ecs_radius=ecs_radius,
        ecs_length=ecs_length, ecs_angle=ecs_angle, rtol=rtol, atol=atol
    )
    return left - right


def find_resonance(
    guess: complex,
    epsilon: float = 0.0,
    length: float = 20.0,
    *,
    ecs_radius: float = 28.0,
) -> tuple[complex, float, bool]:
    def objective(parts: np.ndarray) -> np.ndarray:
        value = resonance_mismatch(complex(parts[0], parts[1]), epsilon, length, ecs_radius=ecs_radius)
        return np.asarray([value.real, value.imag])

    result = root(objective, np.asarray([guess.real, guess.imag]), method="hybr", tol=2.0e-9)
    omega = complex(result.x[0], result.x[1])
    residual = abs(resonance_mismatch(omega, epsilon, length, ecs_radius=ecs_radius))
    return omega, float(residual), bool(result.success and residual < 2.0e-6)


def resonance_linear_coefficient(
    omega0: complex,
    length: float,
    *,
    ecs_radius: float | None = None,
    h_omega: float = 2.0e-6,
    h_epsilon: float = 2.0e-7,
) -> complex:
    radius = max(28.0, length + 3.0) if ecs_radius is None else ecs_radius
    f_plus = resonance_mismatch(omega0 + h_omega, 0.0, length, ecs_radius=radius)
    f_minus = resonance_mismatch(omega0 - h_omega, 0.0, length, ecs_radius=radius)
    derivative_omega = (f_plus - f_minus) / (2.0 * h_omega)
    e_plus = resonance_mismatch(omega0, h_epsilon, length, ecs_radius=radius)
    e_minus = resonance_mismatch(omega0, -h_epsilon, length, ecs_radius=radius)
    derivative_epsilon = (e_plus - e_minus) / (2.0 * h_epsilon)
    return -derivative_epsilon / derivative_omega


def run_time_experiments(
    epsilons: list[float],
    c_value: float,
    resolutions: list[float],
    tmax: float,
) -> tuple[list[dict[str, object]], dict[tuple[float, float], tuple[Waveform, Waveform]]]:
    results: list[dict[str, object]] = []
    waveforms: dict[tuple[float, float], tuple[Waveform, Waveform]] = {}
    for dx in resolutions:
        base = simulate_waveform(0.0, 20.0, dx=dx, tmax=tmax)
        denominator = float(np.trapezoid(base.h * base.h, base.t))
        tail = float(np.trapezoid(base.h[base.t >= tmax - 20.0] ** 2, base.t[base.t >= tmax - 20.0]))
        for epsilon in epsilons:
            length = c_value * math.log(1.0 / epsilon)
            perturbed = simulate_waveform(epsilon, length, dx=dx, tmax=tmax)
            causal = earliest_influence_time(length)
            cutoffs = [max(0.0, causal - 1.0), min(tmax, causal + 12.0), tmax]
            metrics = [waveform_metrics(base, perturbed, cutoff) for cutoff in cutoffs]
            prec_mask = base.t <= max(0.0, causal - 0.5)
            precursor = float(np.max(np.abs((perturbed.h - base.h)[prec_mask]))) if np.any(prec_mask) else 0.0
            results.append(
                {
                    "epsilon": epsilon,
                    "c": c_value,
                    "L": length,
                    "earliest_exact_time": causal,
                    "dx": base.dx,
                    "dt": base.dt,
                    "domain": [base.xmin, base.xmax],
                    "initial_energy": base.energy0,
                    "base_norm2_to_tmax": denominator,
                    "last_20_time_units_fraction": tail / denominator,
                    "numerical_precursor_max": precursor,
                    "metrics": metrics,
                }
            )
            waveforms[(dx, epsilon)] = (base, perturbed)
    return results, waveforms


def make_plot(
    waveforms: dict[tuple[float, float], tuple[Waveform, Waveform]],
    finest_dx: float,
    epsilons: list[float],
    c_value: float,
    output: Path,
) -> None:
    fig, axes = plt.subplots(len(epsilons), 1, figsize=(8.0, 2.2 * len(epsilons)), sharex=True)
    if len(epsilons) == 1:
        axes = [axes]
    for axis, epsilon in zip(axes, epsilons):
        base, perturbed = waveforms[(finest_dx, epsilon)]
        difference = perturbed.h - base.h
        causal = earliest_influence_time(c_value * math.log(1.0 / epsilon))
        axis.plot(base.t, base.h, color="#17365D", lw=1.0, label=r"$h_0$")
        axis.plot(base.t, difference / epsilon, color="#B54A35", lw=1.0, label=r"$(h_\epsilon-h_0)/\epsilon$")
        axis.axvline(causal, color="black", ls="--", lw=0.8, label=r"$2L-13$")
        axis.set_ylabel(f"eps={epsilon:g}")
        axis.grid(alpha=0.2)
        axis.legend(loc="upper right", ncol=3, fontsize=8)
    axes[-1].set_xlabel(r"$\tau$")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmax", type=float, default=100.0)
    parser.add_argument("--coarse-dx", type=float, default=0.08)
    parser.add_argument("--fine-dx", type=float, default=0.04)
    parser.add_argument("--c", type=float, default=5.0)
    parser.add_argument("--output", type=Path, default=Path("results_03.json"))
    parser.add_argument("--plot", type=Path, default=Path("waveforms_03.png"))
    parser.add_argument("--skip-resonance", action="store_true")
    args = parser.parse_args()

    epsilons = [0.04, 0.02, 0.01]
    resolutions = [args.coarse_dx, args.fine_dx]
    time_results, waveforms = run_time_experiments(epsilons, args.c, resolutions, args.tmax)
    make_plot(waveforms, args.fine_dx, epsilons, args.c, args.plot)

    resonance_results: dict[str, object] = {"status": "NOT_RUN"}
    if not args.skip_resonance:
        guess = wkb_fundamental_guess()
        omega_28, residual_28, ok_28 = find_resonance(guess, ecs_radius=28.0)
        omega_32, residual_32, ok_32 = find_resonance(omega_28, ecs_radius=32.0)
        coefficients = []
        for length in [12.0, 16.0, 20.0]:
            coefficient = resonance_linear_coefficient(omega_28, length)
            check_epsilon = 1.0e-4
            radius = max(28.0, length + 3.0)
            predicted = omega_28 + check_epsilon * coefficient
            checked_root, checked_residual, checked_ok = find_resonance(
                predicted, check_epsilon, length, ecs_radius=radius
            )
            coefficients.append(
                {
                    "L": length,
                    "domega_depsilon": [coefficient.real, coefficient.imag],
                    "magnitude": abs(coefficient),
                    "scaled_magnitude": abs(coefficient) * math.exp(2.0 * omega_28.imag * length),
                    "check_epsilon": check_epsilon,
                    "exact_perturbed_root": [checked_root.real, checked_root.imag],
                    "exact_root_residual": checked_residual,
                    "exact_root_success": checked_ok,
                    "linear_prediction_error": abs(checked_root - predicted),
                }
            )
        resonance_results = {
            "status": "RUN",
            "wkb_guess": [guess.real, guess.imag],
            "root_ecs_radius_28": [omega_28.real, omega_28.imag],
            "residual_ecs_radius_28": residual_28,
            "success_ecs_radius_28": ok_28,
            "root_ecs_radius_32": [omega_32.real, omega_32.imag],
            "residual_ecs_radius_32": residual_32,
            "success_ecs_radius_32": ok_32,
            "ecs_radius_difference": abs(omega_28 - omega_32),
            "linear_coefficients": coefficients,
        }

    payload = {
        "normalization": {
            "C": normalization_constant()[0],
            "raw_energy_C1": normalization_constant()[1],
        },
        "time_domain": time_results,
        "resonance": resonance_results,
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
