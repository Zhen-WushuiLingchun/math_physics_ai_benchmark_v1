import math

import numpy as np

import ringdown_stability as r


def test_tortoise_inverse_and_potential():
    x = np.asarray([-20.0, 0.0, 10.0, 40.0])
    rho = r.rho_of_x(x)
    reconstructed = rho + 2.0 * np.log(rho / 2.0 - 1.0)
    assert np.max(np.abs(reconstructed - x)) < 2e-10
    assert np.all(r.rw_potential(x) >= 0.0)


def test_bump_support_smoothness_and_derivative():
    assert r.bump(-1.0) == 0.0
    assert r.bump(1.0) == 0.0
    assert r.bump(0.0) == 1.0
    points = np.asarray([-0.7, -0.2, 0.3, 0.8])
    step = 1e-6
    finite_difference = (r.bump(points + step) - r.bump(points - step)) / (2.0 * step)
    assert np.max(np.abs(finite_difference - r.bump_prime(points))) < 2e-8


def test_initial_energy_normalization_and_causal_time():
    constant, raw = r.normalization_constant()
    assert abs(constant * constant * raw - 1.0) < 2e-13
    assert abs(r.earliest_influence_time(20.0) - 27.0) < 1e-14


def test_time_solver_convergence_and_precursor_control():
    kwargs = dict(epsilon=0.02, length=14.0, tmax=24.0, xmin=-30.0, xmax=50.0)
    coarse = r.simulate_waveform(dx=0.08, **kwargs)
    fine = r.simulate_waveform(dx=0.04, **kwargs)
    base_coarse = r.simulate_waveform(0.0, 14.0, dx=0.08, tmax=24.0, xmin=-30.0, xmax=50.0)
    base_fine = r.simulate_waveform(0.0, 14.0, dx=0.04, tmax=24.0, xmin=-30.0, xmax=50.0)
    assert abs(coarse.energy0 - 1.0) < 2e-12
    assert abs(fine.energy0 - 1.0) < 2e-12
    causal = r.earliest_influence_time(14.0)
    coarse_metric = r.waveform_metrics(base_coarse, coarse, causal + 8.0)["absolute_error"]
    fine_metric = r.waveform_metrics(base_fine, fine, causal + 8.0)["absolute_error"]
    assert coarse_metric > 0.0 and fine_metric > 0.0
    assert abs(coarse_metric - fine_metric) < 0.012
    pre = base_fine.t <= causal - 0.5
    assert np.max(np.abs((fine.h - base_fine.h)[pre])) < 2e-7


def test_wkb_guess_is_damped_and_resonance_residual_decreases():
    guess = r.wkb_fundamental_guess()
    assert guess.real > 0.0 and guess.imag < 0.0
    before = abs(r.resonance_mismatch(guess, ecs_radius=26.0, rtol=2e-8, atol=2e-9))
    omega, residual, success = r.find_resonance(guess, ecs_radius=26.0)
    assert success
    assert omega.real > 0.0 and omega.imag < 0.0
    assert residual < 2e-6
    assert residual < before * 1e-3
    omega_second, residual_second, success_second = r.find_resonance(omega, ecs_radius=30.0)
    assert success_second and residual_second < 2e-6
    assert abs(omega_second - omega) < 2e-7
