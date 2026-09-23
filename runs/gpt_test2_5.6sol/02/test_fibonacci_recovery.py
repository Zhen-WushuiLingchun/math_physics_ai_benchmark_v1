from fractions import Fraction

import numpy as np

import fibonacci_recovery as f


def test_fibonacci_sector_dimensions_and_gluing():
    d32 = f.sector_data(3, 2)
    assert d32.m == (1, 2)
    assert d32.n == (1, 1)
    assert d32.q == (1, 2)
    assert d32.dimension == 3

    d24 = f.sector_data(2, 4)
    assert d24.m == (1, 1)
    assert d24.n == (2, 3)
    assert d24.q == (2, 3)
    assert d24.dimension == 5


def test_exact_entropy_and_replica_calibrations():
    d32 = f.sector_data(3, 2)
    assert f.mean_sector_shannon(d32) == Fraction(1, 2)
    assert f.mean_algebra_entropy(d32) == Fraction(1, 2)
    assert f.mean_algebra_purity(d32) == Fraction(2, 3)
    assert abs(
        f.mean_quantum_trace_entropy(d32)
        - (0.5 + Fraction(2, 3) * np.log(f.PHI))
    ) < 1e-14

    d24 = f.sector_data(2, 4)
    assert f.mean_algebra_entropy(d24) == Fraction(7, 12)
    assert f.mean_algebra_purity(d24) == Fraction(3, 5)


def test_exact_haar_stiefel_c2_formulas():
    d32 = f.sector_data(3, 2)
    assert f.mean_c2(d32, 2, "R") == Fraction(1, 8)
    assert f.mean_c2(d32, 2, "B") == Fraction(3, 8)
    assert f.mean_sector_label_c2(d32, 2) == Fraction(1, 8)

    d24 = f.sector_data(2, 4)
    assert f.mean_c2(d24, 2, "R") == Fraction(3, 8)
    assert f.mean_c2(d24, 2, "B") == Fraction(3, 40)
    assert f.mean_sector_label_c2(d24, 2) == Fraction(3, 40)


def test_haar_sampler_and_monte_carlo_holdout():
    data = f.sector_data(3, 2)
    rng = np.random.default_rng(92017)
    values = []
    for _ in range(1200):
        v = f.haar_isometry(data.dimension, 2, rng)
        assert np.linalg.norm(v.conj().T @ v - np.eye(2)) < 2e-12
        values.append(f.c2_sample(v, data, "R"))
    assert abs(np.mean(values) - float(Fraction(1, 8))) < 0.008


def test_decoder_sdp_and_monogamy_control():
    data = f.sector_data(3, 2)
    v = f.haar_isometry(data.dimension, 2, np.random.default_rng(123456))
    fr = f.optimal_recovery_fidelity(v, data, "R")
    fb = f.optimal_recovery_fidelity(v, data, "B")
    assert 0.25 - 1e-6 <= fr <= 1.0 + 1e-6
    assert 0.25 - 1e-6 <= fb <= 1.0 + 1e-6
    assert fr + fb <= 1.5 + 2e-5
    assert f.sector_label_delta(v, data, "R") <= 1.0 + 1e-12
