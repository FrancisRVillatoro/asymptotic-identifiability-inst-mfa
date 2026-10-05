#!/usr/bin/env python3
"""
Numerical asymptotic-identifiability toolkit.

Core input
----------
A coefficient list [A0, A1, ..., AM] representing

    A(epsilon) = A0 + epsilon A1 + ... + epsilon**M AM,

where A is a weighted sensitivity matrix.

Core output
-----------
A recursive Schur-complement deflation that determines how many singular
values scale as epsilon**m at each order m.  The corresponding Fisher
information eigenvalues scale as epsilon**(2m).

The coefficient matrices may be supplied analytically or estimated from
evaluations of A(epsilon) with `fit_matrix_series`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np


Array = np.ndarray


@dataclass
class DeflationLevel:
    order: int
    rank_lifted: int
    singular_values: Array
    leading_directions: Array


@dataclass
class DeflationResult:
    levels: list[DeflationLevel]
    unresolved_directions: Array

    @property
    def singular_orders(self) -> list[Optional[int]]:
        orders: list[Optional[int]] = []
        for level in self.levels:
            orders.extend([level.order] * level.rank_lifted)
        orders.extend([None] * self.unresolved_directions.shape[1])
        return orders

    @property
    def fim_exponents(self) -> list[Optional[int]]:
        return [
            None if order is None else 2 * order
            for order in self.singular_orders
        ]


def fit_matrix_series(
    matrix_function: Callable[[float], Array],
    *,
    degree: int = 5,
    epsilon_scale: float = 2.0e-3,
    scaled_nodes: Optional[Array] = None,
) -> list[Array]:
    """
    Estimate a matrix power series from evaluations at positive epsilon.

    The fitted variable is z = epsilon / epsilon_scale.  Therefore the
    returned coefficients describe A(epsilon_scale*z) as a series in z.
    Orders in z and epsilon are identical, which is all the deflation
    algorithm requires.

    This is an audit utility, not a computer-assisted proof.  For rigorous
    work, supply analytically derived or interval-certified coefficients.
    """
    if degree < 0:
        raise ValueError("degree must be non-negative")
    if epsilon_scale <= 0:
        raise ValueError("epsilon_scale must be positive")

    if scaled_nodes is None:
        scaled_nodes = np.linspace(0.2, 1.0, max(degree + 5, 14))
    scaled_nodes = np.asarray(scaled_nodes, dtype=float)

    matrices = [
        np.asarray(matrix_function(epsilon_scale * z), dtype=float)
        for z in scaled_nodes
    ]
    shape = matrices[0].shape
    if any(matrix.shape != shape for matrix in matrices):
        raise ValueError("matrix_function returned inconsistent shapes")

    vandermonde = np.vander(
        scaled_nodes, N=degree + 1, increasing=True
    )
    samples = np.stack(
        [matrix.reshape(-1) for matrix in matrices], axis=0
    )
    coefficients, *_ = np.linalg.lstsq(
        vandermonde, samples, rcond=None
    )

    return [
        coefficients[j].reshape(shape)
        for j in range(degree + 1)
    ]


def _series_multiply(
    left: list[Array],
    right: list[Array],
    degree: int,
) -> list[Array]:
    rows = left[0].shape[0]
    cols = right[0].shape[1]
    result: list[Array] = []

    for n in range(degree + 1):
        coefficient = np.zeros((rows, cols), dtype=float)
        for k in range(n + 1):
            if k < len(left) and n - k < len(right):
                coefficient += left[k] @ right[n - k]
        result.append(coefficient)

    return result


def _series_inverse(
    series: list[Array],
    degree: int,
) -> list[Array]:
    inverse0 = np.linalg.inv(series[0])
    inverse = [inverse0]

    for n in range(1, degree + 1):
        accumulation = np.zeros_like(series[0], dtype=float)
        for k in range(1, n + 1):
            if k < len(series):
                accumulation += series[k] @ inverse[n - k]
        inverse.append(-inverse0 @ accumulation)

    return inverse


def deflate_matrix_series(
    coefficients: list[Array],
    *,
    relative_rank_tolerance: float = 1.0e-5,
    absolute_rank_tolerance: float = 1.0e-10,
) -> DeflationResult:
    """
    Recursively determine asymptotic singular-value orders.

    At each level the nonzero constant block is removed by a matrix-series
    Schur complement.  Division by epsilon produces the next deflated
    problem.

    If r_m directions are lifted at level m, then:
        singular values ~ epsilon**m,
        FIM eigenvalues ~ epsilon**(2m).

    Directions unresolved after all supplied coefficients are classified
    as structural/unresolved at the available expansion order.
    """
    if not coefficients:
        raise ValueError("at least one coefficient matrix is required")

    current = [np.asarray(c, dtype=float) for c in coefficients]
    rows, parameter_count = current[0].shape
    if any(c.shape != (rows, parameter_count) for c in current):
        raise ValueError("all coefficient matrices must have equal shape")

    parameter_basis = np.eye(parameter_count)
    levels: list[DeflationLevel] = []
    order = 0

    while current and current[0].shape[1] > 0:
        constant = current[0]
        u, singular_values, vt = np.linalg.svd(
            constant, full_matrices=True
        )

        reference = (
            singular_values[0]
            if singular_values.size and singular_values[0] > 0
            else 1.0
        )
        tolerance = max(
            absolute_rank_tolerance,
            relative_rank_tolerance * reference,
        )
        rank = int(np.sum(singular_values > tolerance))

        v = vt.T
        lifted = (
            parameter_basis @ v[:, :rank]
            if rank
            else np.zeros((parameter_count, 0))
        )
        levels.append(
            DeflationLevel(
                order=order,
                rank_lifted=rank,
                singular_values=singular_values.copy(),
                leading_directions=lifted,
            )
        )

        column_count = constant.shape[1]
        if rank == column_count:
            return DeflationResult(
                levels=levels,
                unresolved_directions=np.zeros((parameter_count, 0)),
            )

        null_basis = parameter_basis @ v[:, rank:]
        degree = len(current) - 1
        if degree == 0:
            return DeflationResult(levels, null_basis)

        transformed = [u.T @ coefficient @ v for coefficient in current]

        b11 = [coefficient[:rank, :rank] for coefficient in transformed]
        b12 = [coefficient[:rank, rank:] for coefficient in transformed]
        b21 = [coefficient[rank:, :rank] for coefficient in transformed]
        b22 = [coefficient[rank:, rank:] for coefficient in transformed]

        if rank == 0:
            schur = b22
        else:
            inverse11 = _series_inverse(b11, degree)
            product = _series_multiply(
                _series_multiply(b21, inverse11, degree),
                b12,
                degree,
            )
            schur = [
                b22[j] - product[j]
                for j in range(degree + 1)
            ]

        # The constant coefficient is zero after exact elimination.
        # Removing it is equivalent to division by epsilon.
        current = schur[1:]
        parameter_basis = null_basis
        order += 1

    return DeflationResult(
        levels=levels,
        unresolved_directions=np.zeros((parameter_count, 0)),
    )


def exact_fim_exponents(
    matrix_function: Callable[[float], Array],
    epsilon_values: Array,
    *,
    fit_count: int = 8,
    structural_threshold: float = 1.0e-22,
) -> list[Optional[float]]:
    """
    Fit log-log slopes of exact FIM eigenvalues.

    Eigenvalues are ordered from largest to smallest, so their order agrees
    with the deflation output: O(1) directions first, weakest directions
    last.  A direction whose eigenvalue remains below structural_threshold
    is reported as None.
    """
    epsilon_values = np.asarray(epsilon_values, dtype=float)
    spectra = []

    for epsilon in epsilon_values:
        matrix = np.asarray(matrix_function(float(epsilon)), dtype=float)
        singular_values = np.linalg.svd(matrix, compute_uv=False)
        spectra.append(np.sort(singular_values**2)[::-1])

    spectra = np.asarray(spectra)
    exponents: list[Optional[float]] = []

    for column in range(spectra.shape[1]):
        values = spectra[:, column]
        if np.max(values) < structural_threshold:
            exponents.append(None)
            continue

        x = np.log(epsilon_values[-fit_count:])
        y = np.log(values[-fit_count:])
        exponents.append(float(np.polyfit(x, y, 1)[0]))

    return exponents


def format_orders(values: list[Optional[float]], digits: int = 3) -> str:
    formatted = []
    for value in values:
        if value is None:
            formatted.append("structural")
        elif isinstance(value, (int, np.integer)):
            formatted.append(str(int(value)))
        else:
            formatted.append(f"{float(value):.{digits}f}")
    return ",".join(formatted)
