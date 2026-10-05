#!/usr/bin/env python3
"""Cross-validation of the asymptotic-identifiability deflation algorithm."""

from __future__ import annotations

from pathlib import Path
import csv
import importlib.util
import json

import numpy as np
import matplotlib.pyplot as plt

from asymptotic_identifiability import (
    deflate_matrix_series,
    exact_fim_exponents,
    fit_matrix_series,
    format_orders,
)


ROOT = Path(__file__).resolve().parent


def import_local(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


guppy = import_local("guppy", "phaseI_guppy_fim.py")
minimal = import_local("minimal", "phaseII_minimal_emu.py")
condensation = import_local(
    "condensation", "phaseII_condensation_emu.py"
)
cascade = import_local("cascade", "phaseII_coupled_fast_emu.py")
exchange = import_local(
    "exchange", "phaseII_reversible_exchange_emu.py"
)


def guppy_matrix(epsilon: float) -> np.ndarray:
    theta = np.array([
        np.log(guppy.K1),
        guppy.ALPHA1,
        guppy.BETA,
        np.log(guppy.K1 / epsilon),
        guppy.ALPHA2,
    ])
    blocks = []
    for eta in guppy.ETAS:
        _, sensitivity = guppy.state_and_sensitivities(
            theta, eta / guppy.K1
        )
        blocks.append(sensitivity)
    return np.vstack(blocks)


def minimal_matrix(
    epsilon: float,
    tracers: tuple[int, ...],
) -> np.ndarray:
    theta = np.array([
        np.log(minimal.KX),
        minimal.RHO,
        np.log(minimal.KX / epsilon),
        minimal.ETA,
        minimal.ALPHA,
    ])
    return minimal.stacked_jacobian(theta, tracers)


def condensation_matrix(epsilon: float) -> np.ndarray:
    theta = np.array([
        np.log(condensation.K_A),
        np.log(condensation.K_B),
        np.log(condensation.K_A / epsilon),
        condensation.RHO,
    ])
    return np.vstack([
        condensation.jacobian(
            theta, time, ("A", "B"), True
        )
        for time in condensation.TIMES
    ])


def cascade_matrix(
    epsilon: float,
    observation: str,
) -> np.ndarray:
    kx = cascade.K_A / epsilon
    ky = cascade.CHI * kx
    theta = np.array([
        np.log(cascade.K_A),
        np.log(cascade.K_B),
        np.log(kx),
        np.log(ky),
        cascade.RHO,
        cascade.ETA,
        cascade.ALPHA,
    ])
    full = cascade.stacked_jacobian(theta, observation)
    return full[:, 2:4]


def exchange_matrix(
    epsilon: float,
    experiments: tuple[str, ...],
    observation: str,
) -> np.ndarray:
    theta = np.log([
        exchange.K_A / epsilon,
        exchange.PHI,
        exchange.PSI,
    ])
    return np.vstack([
        exchange.jacobian(
            theta, time, experiments, observation
        )
        for time in exchange.TIMES
    ])


CASES = [
    (
        "Guppy reversible KFP",
        guppy_matrix,
        [0, 0, 0, 1, 2],
    ),
    (
        "Minimal EMU, one tracer",
        lambda eps: minimal_matrix(eps, (1,)),
        [0, 0, 0, 1, None],
    ),
    (
        "Minimal EMU, two tracers",
        lambda eps: minimal_matrix(eps, (1, 2)),
        [0, 0, 0, 0, 1],
    ),
    (
        "Condensation EMU, A+B",
        condensation_matrix,
        [0, 0, 0, 1],
    ),
    (
        "Two-stage cascade, Y only",
        lambda eps: cascade_matrix(eps, "Y"),
        [1, 2],
    ),
    (
        "Two-stage cascade, X+Y",
        lambda eps: cascade_matrix(eps, "XY"),
        [1, 1],
    ),
    (
        "Reversible exchange, one input/Z",
        lambda eps: exchange_matrix(eps, ("X",), "Z"),
        [0, 1, 2],
    ),
    (
        "Reversible exchange, two inputs/Z",
        lambda eps: exchange_matrix(eps, ("X", "Y"), "Z"),
        [0, 0, 1],
    ),
    (
        "Reversible exchange, one input/XYZ",
        lambda eps: exchange_matrix(eps, ("X",), "XYZ"),
        [0, 1, 1],
    ),
]


def main() -> None:
    scaled_nodes = np.linspace(0.2, 1.0, 20)
    epsilon_values = np.logspace(-2, -4, 16)
    rows = []
    all_passed = True

    for case_name, matrix_function, expected_orders in CASES:
        coefficients = fit_matrix_series(
            matrix_function,
            degree=5,
            epsilon_scale=2.0e-3,
            scaled_nodes=scaled_nodes,
        )
        result = deflate_matrix_series(
            coefficients,
            relative_rank_tolerance=1.0e-5,
            absolute_rank_tolerance=1.0e-10,
        )
        predicted_orders = result.singular_orders
        observed_fim = exact_fim_exponents(
            matrix_function,
            epsilon_values,
            fit_count=8,
        )

        passed = predicted_orders == expected_orders
        all_passed = all_passed and passed

        rows.append({
            "case": case_name,
            "expected_singular_orders": format_orders(expected_orders),
            "predicted_singular_orders": format_orders(predicted_orders),
            "predicted_fim_exponents": format_orders(
                result.fim_exponents
            ),
            "observed_fim_exponents": format_orders(observed_fim),
            "status": "PASS" if passed else "FAIL",
        })

    # Direction audit for the Guppy model.
    guppy_coefficients = fit_matrix_series(
        guppy_matrix,
        degree=5,
        epsilon_scale=2.0e-3,
        scaled_nodes=scaled_nodes,
    )
    guppy_result = deflate_matrix_series(
        guppy_coefficients,
        relative_rank_tolerance=1.0e-5,
        absolute_rank_tolerance=1.0e-10,
    )

    q_fast = np.array([0.0, 0.0, 0.0, 1.0, 0.0])
    c = 1.0 - guppy.ALPHA2
    r = 1.0 - guppy.BETA * c
    x_star = (
        guppy.ALPHA1 + guppy.BETA * guppy.ALPHA2
    ) / r
    y_star = c * x_star + guppy.ALPHA2
    q_degenerate = np.array([c / r, -y_star, 1.0, 0.0, 0.0])
    q_degenerate /= np.linalg.norm(q_degenerate)

    order1_direction = guppy_result.levels[1].leading_directions[:, 0]
    order1_direction /= np.linalg.norm(order1_direction)
    order2_direction = guppy_result.levels[2].leading_directions[:, 0]
    order2_direction /= np.linalg.norm(order2_direction)

    fast_alignment = abs(float(order1_direction @ q_fast))
    degenerate_alignment = abs(
        float(order2_direction @ q_degenerate)
    )

    rows.append({
        "case": "Guppy direction audit",
        "expected_singular_orders": "q_fast@1;q_d@2",
        "predicted_singular_orders": (
            f"alignment={fast_alignment:.9f};"
            f"{degenerate_alignment:.9f}"
        ),
        "predicted_fim_exponents": "2;4",
        "observed_fim_exponents": "",
        "status": (
            "PASS"
            if fast_alignment > 0.999 and degenerate_alignment > 0.999
            else "FAIL"
        ),
    })

    with open(
        ROOT / "asymptotic_identifiability_validation.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Plot observed versus predicted FIM exponents for nonstructural modes.
    predicted = []
    observed = []
    labels = []

    for case_name, matrix_function, expected_orders in CASES:
        exact = exact_fim_exponents(
            matrix_function, epsilon_values, fit_count=8
        )
        expected_fim = [
            None if order is None else 2 * order
            for order in expected_orders
        ]
        for index, (theory, empirical) in enumerate(
            zip(expected_fim, exact), start=1
        ):
            if theory is not None and empirical is not None:
                predicted.append(theory)
                observed.append(empirical)
                labels.append(f"{case_name}:{index}")

    plt.figure(figsize=(6.8, 5.2))
    plt.scatter(predicted, observed)
    minimum = min(predicted + observed) - 0.2
    maximum = max(predicted + observed) + 0.2
    plt.plot([minimum, maximum], [minimum, maximum])
    plt.xlabel("Predicted FIM exponent")
    plt.ylabel("Measured log-log exponent")
    plt.tight_layout()
    plt.savefig(
        ROOT / "asymptotic_identifiability_validation.png",
        dpi=200,
    )
    plt.close()

    print("ALL_PASS =", all_passed)
    print("Guppy fast-direction alignment =", fast_alignment)
    print("Guppy degenerate-direction alignment =", degenerate_alignment)
    for row in rows:
        print(
            row["status"],
            "|",
            row["case"],
            "| predicted:",
            row["predicted_singular_orders"],
            "| observed FIM:",
            row["observed_fim_exponents"],
        )


if __name__ == "__main__":
    main()
