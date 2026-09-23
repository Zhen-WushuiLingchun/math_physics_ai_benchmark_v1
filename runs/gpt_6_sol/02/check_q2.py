"""Exact and Monte Carlo checks for the Fibonacci random-code problem.

All paths are relative to this directory.  Monte Carlo uses a fixed seed and
Haar/Stiefel matrices obtained by QR of a complex Gaussian matrix.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import sqrt

import numpy as np
import sympy as sp


def fib(n: int) -> int:
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def dimensions(nr: int, nb: int):
    assert nr >= 2 and nb >= 2
    m = (fib(nr - 1), fib(nr))
    n = (fib(nb - 1), fib(nb))
    d = sum(u * v for u, v in zip(m, n))
    assert d == fib(nr + nb - 1)
    return m, n, d


def entropy_mean(nr: int, nb: int):
    """Exact E[S_alg], E[S_qtr] (natural logarithms)."""
    m, n, d = dimensions(nr, nb)
    blocks = [u * v for u, v in zip(m, n)]
    alg = sp.harmonic(d)
    for u, v, h in zip(m, n, blocks):
        big, small = max(u, v), min(u, v)
        alg -= sp.Rational(h, d) * (
            sp.harmonic(big) + sp.Rational(small - 1, 2 * big)
        )
    phi = (1 + sp.sqrt(5)) / 2
    qtr = alg + sp.Rational(blocks[1], d) * sp.log(phi)
    return sp.simplify(alg), sp.simplify(qtr)


def pure_second_moment(nr: int, nb: int):
    """Exact E Tr(rho_R^2) and quantum-trace analogue."""
    m, n, d = dimensions(nr, nb)
    pieces = [sp.Rational(u * v * (u + v), d * (d + 1))
              for u, v in zip(m, n)]
    phi = (1 + sp.sqrt(5)) / 2
    return sp.simplify(sum(pieces)), sp.simplify(pieces[0] + pieces[1] / phi)


def mixed_moments(nr: int, nb: int, k: int):
    """Exact Stiefel moments of C_QB, C_QR, rho_QR^2 and center C_QA."""
    m, n, d = dimensions(nr, nb)
    assert 2 <= k <= d
    s_mnn = sum(u * v * v for u, v in zip(m, n))
    s_mmn = sum(u * u * v for u, v in zip(m, n))
    factor = sp.Rational(k * k - 1, k * k) / (d * d - 1)
    c_qb = sp.simplify(factor * (s_mnn - sp.Rational(s_mmn, d)))
    c_qr = sp.simplify(factor * (s_mmn - sp.Rational(s_mnn, d)))
    p_qr = sp.simplify(
        (s_mmn + sp.Rational(s_mnn, k)
         - sp.Rational(1, d) * (s_mnn + sp.Rational(s_mmn, k)))
        / (d * d - 1)
    )
    h = [u * v for u, v in zip(m, n)]
    center = sp.simplify(
        factor * (d - sp.Rational(sum(t * t for t in h), d))
    )
    return c_qb, c_qr, p_qr, center


def haar_isometry(d: int, k: int, rng: np.random.Generator):
    g = rng.normal(size=(d, k)) + 1j * rng.normal(size=(d, k))
    q, r = np.linalg.qr(g, mode="reduced")
    diag = np.diag(r)
    q = q * (diag / np.abs(diag))[None, :]
    return q


def block_matrices(v: np.ndarray, m, n):
    pos = 0
    blocks = []
    for u, w in zip(m, n):
        blocks.append(v[pos:pos + u * w, :].reshape(u, w, v.shape[1]))
        pos += u * w
    return blocks


def complement_centered_purity(blocks, k: int):
    """Direct matrix construction of C_QB; no independent-column shortcut."""
    total = 0.0
    for x in blocks:
        _, nb, _ = x.shape
        rho4 = np.einsum("rbi,rcj->ibjc", x, x.conj(), optimize=True) / k
        rho = rho4.reshape(k * nb, k * nb)
        rho_b = np.einsum("ibic->bc", rho4)
        c = rho - np.kron(np.eye(k) / k, rho_b)
        total += float(np.vdot(c, c).real)
    return total


def exact_decoder_value(blocks, k: int):
    """Optimum for the two selected small cases, from explicit branch decoders."""
    if all(x.shape[1] == 1 for x in blocks):
        # One Kraus operator per classical sector; polar recovery is optimal.
        total = 0.0
        for x in blocks:
            a = x[:, 0, :]
            vals = np.linalg.eigvalsh(a.conj().T @ a).clip(min=0)
            total += float(np.sqrt(vals).sum() ** 2)
        return total / (k * k)
    if all(x.shape[0] == 1 for x in blocks):
        # A classical-output measurement; prepare the top eigenvector per outcome.
        total = 0.0
        for x in blocks:
            a = x[0, :, :]
            vals = np.linalg.eigvalsh(a.conj().T @ a)
            total += float(vals[-1])
        return total / (k * k)
    raise ValueError("No closed-form decoder implemented for this cut")


def run_sample(nr: int, nb: int, k: int, samples: int, seed: int):
    m, n, d = dimensions(nr, nb)
    rng = np.random.default_rng(seed)
    fvals, cvals = [], []
    for _ in range(samples):
        blocks = block_matrices(haar_isometry(d, k, rng), m, n)
        fvals.append(exact_decoder_value(blocks, k))
        cvals.append(complement_centered_purity(blocks, k))
    return float(np.mean(fvals)), float(np.mean(cvals))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=6000)
    args = parser.parse_args()

    for nr, nb in ((2, 2), (3, 2), (2, 3)):
        print((nr, nb), "dims", dimensions(nr, nb),
              "E S_alg/qtr", entropy_mean(nr, nb),
              "E purity alg/qtr", pure_second_moment(nr, nb),
              "E C_QB^2/C_QR^2/P_QR/C_QA^2", mixed_moments(nr, nb, 2))

    assert entropy_mean(3, 2) == entropy_mean(2, 3)
    assert pure_second_moment(3, 2) == pure_second_moment(2, 3)
    assert mixed_moments(3, 2, 2)[0] == sp.Rational(1, 8)
    assert mixed_moments(2, 3, 2)[0] == sp.Rational(3, 8)

    # For both D=3 cuts p=Tr(V†P_1V) ~ Beta(2,1).
    f32 = sp.Rational(23, 30)
    f23 = sp.Rational(5, 12)
    print("exact E F_rec (3,2), (2,3):", f32, f23)
    for nr, nb, target in ((3, 2, float(f32)), (2, 3, float(f23))):
        f_mc, c_mc = run_sample(nr, nb, 2, args.samples, 20260923 + nr)
        c_exact = float(mixed_moments(nr, nb, 2)[0])
        print("Monte Carlo", (nr, nb), "samples", args.samples,
              "F", f_mc, "target", target,
              "C2", c_mc, "target", c_exact)
        assert abs(f_mc - target) < 0.025
        assert abs(c_mc - c_exact) < 0.025

    phi = (1 + sqrt(5)) / 2
    c1 = 1 / phi
    z1, z2, z3 = 1 + c1, 1 + c1**2, 1 + c1**3
    gamma = sqrt(z1 * z3) / z2
    print("asymptotic gamma:", gamma, "expected sqrt(2)/Z2:",
          sqrt(2) / z2)
    assert abs(gamma - sqrt(2) / z2) < 1e-12


if __name__ == "__main__":
    main()
