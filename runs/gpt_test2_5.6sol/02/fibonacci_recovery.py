"""Independent finite-dimensional checks for the Fibonacci recovery problem.

The script keeps the two topological-charge sectors as an ordinary direct sum.
Quantum dimensions are used only in the explicitly requested quantum-trace
entropy; they are never treated as Hilbert-space dimensions.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import cvxpy as cp
import numpy as np


PHI = (1.0 + math.sqrt(5.0)) / 2.0


def fibonacci(n: int) -> int:
    if n < 0:
        raise ValueError("Fibonacci index must be nonnegative")
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def harmonic(n: int) -> Fraction:
    return sum((Fraction(1, j) for j in range(1, n + 1)), Fraction(0))


@dataclass(frozen=True)
class SectorData:
    nr: int
    nb: int
    m: tuple[int, int]
    n: tuple[int, int]
    q: tuple[int, int]
    dimension: int


def sector_data(nr: int, nb: int) -> SectorData:
    if nr < 2 or nb < 2:
        raise ValueError("The problem requires nr, nb >= 2")
    m = (fibonacci(nr - 1), fibonacci(nr))
    n = (fibonacci(nb - 1), fibonacci(nb))
    q = (m[0] * n[0], m[1] * n[1])
    dimension = q[0] + q[1]
    if dimension != fibonacci(nr + nb - 1):
        raise AssertionError("Fibonacci gluing identity failed")
    return SectorData(nr, nb, m, n, q, dimension)


def page_entropy(m: int, n: int) -> Fraction:
    small, large = sorted((m, n))
    return harmonic(m * n) - harmonic(large) - Fraction(small - 1, 2 * large)


def mean_sector_shannon(data: SectorData) -> Fraction:
    d = data.dimension
    return harmonic(d) - sum(
        (Fraction(q, d) * harmonic(q) for q in data.q), Fraction(0)
    )


def mean_algebra_entropy(data: SectorData) -> Fraction:
    d = data.dimension
    within = sum(
        (
            Fraction(q, d) * page_entropy(m, n)
            for q, m, n in zip(data.q, data.m, data.n)
        ),
        Fraction(0),
    )
    return mean_sector_shannon(data) + within


def mean_quantum_trace_entropy(data: SectorData) -> float:
    return float(mean_algebra_entropy(data)) + data.q[1] / data.dimension * math.log(PHI)


def mean_algebra_purity(data: SectorData) -> Fraction:
    d = data.dimension
    numerator = sum(m * n * (m + n) for m, n in zip(data.m, data.n))
    return Fraction(numerator, d * (d + 1))


def mean_quantum_trace_purity(data: SectorData) -> float:
    d = data.dimension
    terms = [
        m * n * (m + n) / dim
        for m, n, dim in zip(data.m, data.n, (1.0, PHI))
    ]
    return sum(terms) / (d * (d + 1))


def mean_c2(data: SectorData, k: int, recover: str = "R") -> Fraction:
    """Exact Haar--Stiefel mean of Tr(C_QE^2).

    recover='R' means that B is the complementary system E in the decoupling
    test.  recover='B' swaps m and n.
    """
    if not (2 <= k <= data.dimension):
        raise ValueError("k must satisfy 2 <= k <= D")
    if recover not in {"R", "B"}:
        raise ValueError("recover must be R or B")
    m, n = (data.m, data.n) if recover == "R" else (data.n, data.m)
    d = data.dimension
    s1 = sum(mi * ni * ni for mi, ni in zip(m, n))
    s2 = sum(mi * mi * ni for mi, ni in zip(m, n))
    return Fraction(k * k - 1, k * k) * Fraction(d * s1 - s2, d * (d * d - 1))


def mean_sector_label_c2(data: SectorData, k: int) -> Fraction:
    """Exact mean Tr(C_QA^2) when the complement keeps only charge label A."""
    if not (2 <= k <= data.dimension):
        raise ValueError("k must satisfy 2 <= k <= D")
    d = data.dimension
    sum_q2 = sum(q * q for q in data.q)
    return Fraction(k * k - 1, k * k) * Fraction(d * d - sum_q2, d * (d * d - 1))


def trace_distance_bound(data: SectorData, k: int, recover: str = "R") -> float:
    complement_dims = data.n if recover == "R" else data.m
    rank_bound = k * sum(complement_dims)
    return 0.5 * math.sqrt(rank_bound * float(mean_c2(data, k, recover)))


def recovery_lower_bound(data: SectorData, k: int, recover: str = "R") -> float:
    eps = min(1.0, trace_distance_bound(data, k, recover))
    return (1.0 - eps) ** 2


def haar_isometry(d: int, k: int, rng: np.random.Generator) -> np.ndarray:
    z = (rng.normal(size=(d, k)) + 1j * rng.normal(size=(d, k))) / math.sqrt(2.0)
    q, r = np.linalg.qr(z, mode="reduced")
    phase = np.diag(r)
    phase = np.where(np.abs(phase) > 0, phase / np.abs(phase), 1.0)
    return q * phase.conj()


def split_isometry(v: np.ndarray, data: SectorData) -> list[np.ndarray]:
    blocks: list[np.ndarray] = []
    start = 0
    for m, n in zip(data.m, data.n):
        stop = start + m * n
        blocks.append(v[start:stop, :].reshape(m, n, v.shape[1]))
        start = stop
    if start != v.shape[0]:
        raise AssertionError("Sector dimensions do not exhaust the isometry")
    return blocks


def complementary_c_blocks(
    v: np.ndarray, data: SectorData, recover: str = "R"
) -> list[np.ndarray]:
    """Blocks of C_QE for the environment complementary to the recovered side."""
    k = v.shape[1]
    result: list[np.ndarray] = []
    for x in split_isometry(v, data):
        if recover == "R":
            # Trace R, retain B.
            output_dim = x.shape[1]
            reduced = lambda i, j: x[:, :, i].T @ x[:, :, j].conj()
        elif recover == "B":
            # Trace B, retain R.
            output_dim = x.shape[0]
            reduced = lambda i, j: x[:, :, i] @ x[:, :, j].conj().T
        else:
            raise ValueError("recover must be R or B")

        rho_qe = np.zeros((k * output_dim, k * output_dim), dtype=complex)
        for i in range(k):
            for j in range(k):
                rho_qe[
                    i * output_dim : (i + 1) * output_dim,
                    j * output_dim : (j + 1) * output_dim,
                ] = reduced(i, j) / k
        rho_e = sum(
            rho_qe[
                i * output_dim : (i + 1) * output_dim,
                i * output_dim : (i + 1) * output_dim,
            ]
            for i in range(k)
        )
        result.append(rho_qe - np.kron(np.eye(k) / k, rho_e))
    return result


def c2_sample(v: np.ndarray, data: SectorData, recover: str = "R") -> float:
    return float(
        sum(np.trace(block @ block).real for block in complementary_c_blocks(v, data, recover))
    )


def sector_label_delta(v: np.ndarray, data: SectorData, recover: str = "R") -> float:
    """Trace-distance lower bound obtained by measuring only the sector label."""
    k = v.shape[1]
    total = 0.0
    for block, out_dim in zip(
        complementary_c_blocks(v, data, recover),
        data.n if recover == "R" else data.m,
    ):
        reshaped = block.reshape(k, out_dim, k, out_dim)
        cq = np.einsum("irjr->ij", reshaped)
        total += np.linalg.norm(cq, ord="nuc")
    return 0.5 * float(total)


def channel_output_matrices(
    x: np.ndarray, recover: str
) -> tuple[int, list[list[np.ndarray]]]:
    """N(|i><j|) matrices for one charge sector."""
    k = x.shape[2]
    if recover == "R":
        out_dim = x.shape[0]
        matrix = lambda i, j: x[:, :, i] @ x[:, :, j].conj().T
    elif recover == "B":
        out_dim = x.shape[1]
        matrix = lambda i, j: x[:, :, i].T @ x[:, :, j].conj()
    else:
        raise ValueError("recover must be R or B")
    return out_dim, [[matrix(i, j) for j in range(k)] for i in range(k)]


def optimal_recovery_fidelity(
    v: np.ndarray,
    data: SectorData,
    recover: str = "R",
    solver: str = "CLARABEL",
) -> float:
    """Solve the direct-sum decoder SDP for the entanglement fidelity."""
    k = v.shape[1]
    variables: list[cp.Variable] = []
    constraints: list[cp.Constraint] = []
    objective_terms = []

    for x in split_isometry(v, data):
        out_dim, matrices = channel_output_matrices(x, recover)
        size = out_dim * k
        jdec = cp.Variable((size, size), hermitian=True)
        variables.append(jdec)
        constraints.append(jdec >> 0)
        for r in range(out_dim):
            for s in range(out_dim):
                constraints.append(
                    sum(jdec[r * k + ell, s * k + ell] for ell in range(k))
                    == (1.0 if r == s else 0.0)
                )

        coeff = np.zeros((size, size), dtype=complex)
        for i in range(k):
            for j in range(k):
                aij = matrices[i][j]
                for r in range(out_dim):
                    for s in range(out_dim):
                        coeff[r * k + i, s * k + j] = aij[r, s] / (k * k)
        # sum(coeff[p,q] * J[p,q]) = Tr(coeff^T J).
        objective_terms.append(cp.real(cp.trace(coeff.T @ jdec)))

    problem = cp.Problem(cp.Maximize(sum(objective_terms)), constraints)
    kwargs = {"verbose": False}
    if solver == "SCS":
        kwargs.update({"eps": 2e-7, "max_iters": 100_000})
    value = problem.solve(solver=solver, **kwargs)
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"Decoder SDP failed: {problem.status}")
    return float(np.real(value))


def run_experiment(
    nr: int,
    nb: int,
    k: int,
    samples: int,
    sdp_samples: int,
    seed: int,
) -> dict[str, object]:
    data = sector_data(nr, nb)
    rng = np.random.default_rng(seed)
    c2_r: list[float] = []
    c2_b: list[float] = []
    label_r: list[float] = []
    frec_r: list[float] = []
    frec_b: list[float] = []
    monogamy: list[float] = []
    for sample in range(samples):
        v = haar_isometry(data.dimension, k, rng)
        c2_r.append(c2_sample(v, data, "R"))
        c2_b.append(c2_sample(v, data, "B"))
        label_r.append(sector_label_delta(v, data, "R"))
        if sample < sdp_samples:
            fr = optimal_recovery_fidelity(v, data, "R")
            fb = optimal_recovery_fidelity(v, data, "B")
            frec_r.append(fr)
            frec_b.append(fb)
            monogamy.append(fr + fb)

    def summary(values: Iterable[float]) -> dict[str, float]:
        array = np.asarray(list(values), dtype=float)
        if not array.size:
            return {}
        return {
            "mean": float(array.mean()),
            "std": float(array.std(ddof=1)) if array.size > 1 else 0.0,
            "min": float(array.min()),
            "max": float(array.max()),
        }

    return {
        "sector_data": asdict(data),
        "k": k,
        "samples": samples,
        "sdp_samples": sdp_samples,
        "seed": seed,
        "exact": {
            "mean_S_alg_fraction": str(mean_algebra_entropy(data)),
            "mean_S_alg": float(mean_algebra_entropy(data)),
            "mean_S_qtr": mean_quantum_trace_entropy(data),
            "mean_purity_alg_fraction": str(mean_algebra_purity(data)),
            "mean_purity_qtr": mean_quantum_trace_purity(data),
            "mean_C2_R_fraction": str(mean_c2(data, k, "R")),
            "mean_C2_B_fraction": str(mean_c2(data, k, "B")),
            "mean_C2_sector_label_fraction": str(mean_sector_label_c2(data, k)),
            "delta_R_bound": trace_distance_bound(data, k, "R"),
            "Frec_R_lower_bound": recovery_lower_bound(data, k, "R"),
        },
        "monte_carlo": {
            "C2_R": summary(c2_r),
            "C2_B": summary(c2_b),
            "sector_label_delta_R": summary(label_r),
            "Frec_R_SDP": summary(frec_r),
            "Frec_B_SDP": summary(frec_b),
            "Frec_sum": summary(monogamy),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=4000)
    parser.add_argument("--sdp-samples", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260923)
    parser.add_argument("--output", type=Path, default=Path("results_02.json"))
    args = parser.parse_args()

    configs = [(3, 2, 2), (2, 4, 2), (4, 3, 2)]
    results = [
        run_experiment(nr, nb, k, args.samples, args.sdp_samples, args.seed + i)
        for i, (nr, nb, k) in enumerate(configs)
    ]
    args.output.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
