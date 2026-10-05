#!/usr/bin/env python3
"""
Phase II: reversible fast-exchange EMU block.

Slow labeled precursor:
    da/dt = kA(1-a).

Fast reversible block:
    dX/dt = kF[rho*a*I_X -(1+phi)X + psi*Y]
    dY/dt = kF[rho*a*I_Y + phi*X -(1+psi)Y]

I_X,I_Y select one of two complementary tracer-entry experiments:
    E_X: tracer enters X;
    E_Y: tracer enters Y.

Slow downstream reporter:
    dZ/dt = kZ(eta*Y-Z).

Conditional fast parameters:
    (log kF, log phi, log psi)

Designs:
    X_Z       : E_X, observe Z only
    XY_Z      : E_X+E_Y, observe Z only
    X_XYZ     : E_X, observe X,Y,Z
    XY_XYZ    : E_X+E_Y, observe X,Y,Z

Exact local sensitivities are evaluated by complex-step differentiation
of an augmented matrix exponential.  The reported Fisher matrix is J.T@J,
equivalent to independent observations with equal fixed variance.
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm, svdvals

K_A = 0.04
K_Z = 0.06
PHI = 0.80
PSI = 0.50
RHO = 0.80
ETA = 0.70
TAU = np.linspace(0.10, 5.0, 15)
TIMES = TAU / K_A
Q0 = np.array([0.0, 0.0, 0.0, 0.0, 1.0])

def augmented_matrix(theta, experiment):
    lkF, lphi, lpsi = theta
    kF, phi, psi = np.exp(lkF), np.exp(lphi), np.exp(lpsi)

    # state = (a, X, Y, Z, 1)
    M = np.zeros((5,5), dtype=np.result_type(theta))
    M[0,0], M[0,4] = -K_A, K_A

    M[1,1] = -kF*(1.0+phi)
    M[1,2] =  kF*psi
    M[2,1] =  kF*phi
    M[2,2] = -kF*(1.0+psi)

    if experiment == "X":
        M[1,0] = kF*RHO
    elif experiment == "Y":
        M[2,0] = kF*RHO
    else:
        raise ValueError(experiment)

    M[3,2], M[3,3] = K_Z*ETA, -K_Z
    return M

def response(theta, t, experiments, observation):
    values = []
    for experiment in experiments:
        q = expm(augmented_matrix(theta, experiment)*t) @ Q0
        if observation == "Z":
            values.append(q[3])
        elif observation == "XYZ":
            values.extend((q[1],q[2],q[3]))
        else:
            raise ValueError(observation)
    return np.asarray(values)

def jacobian(theta, t, experiments, observation):
    h = 1e-28
    y = response(theta, t, experiments, observation)
    J = np.zeros((len(y), len(theta)))
    for j in range(len(theta)):
        z = theta.astype(complex)
        z[j] += 1j*h
        J[:,j] = np.imag(
            response(z, t, experiments, observation)
        )/h
    return J

def spectrum(epsilon, experiments, observation):
    kF = K_A/epsilon
    theta = np.log([kF, PHI, PSI])
    J = np.vstack([
        jacobian(theta, t, experiments, observation)
        for t in TIMES
    ])
    return np.sort(svdvals(J)**2)

def main():
    eps_values = np.logspace(-1, -4, 25)
    designs = [
        ("X_Z", ("X",), "Z"),
        ("XY_Z", ("X","Y"), "Z"),
        ("X_XYZ", ("X",), "XYZ"),
        ("XY_XYZ", ("X","Y"), "XYZ"),
    ]

    rows = []
    spectra = {}
    for name, experiments, observation in designs:
        arr = np.asarray([
            spectrum(eps, experiments, observation)
            for eps in eps_values
        ])
        spectra[name] = arr
        for eps, eig in zip(eps_values, arr):
            rows.append([eps, name, *eig])

    with open("phaseII_reversible_exchange_spectrum.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "epsilon", "design",
            "lambda1", "lambda2", "lambda3"
        ])
        writer.writerows(rows)

    plt.figure(figsize=(7.4,5.3))
    plot_items = [
        ("X_Z", 0, "o", "one input, Z: smallest"),
        ("X_Z", 1, "o", "one input, Z: second"),
        ("XY_Z", 0, "s", "two inputs, Z: smallest"),
        ("X_XYZ", 0, "^", "one input, X+Y+Z: smallest"),
    ]
    for name, column, marker, label in plot_items:
        plt.loglog(
            eps_values, spectra[name][:,column],
            marker=marker, label=label
        )

    plt.loglog(
        eps_values, 1e-2*eps_values**4,
        linestyle=":", label=r"reference $\varepsilon^4$"
    )
    plt.loglog(
        eps_values, 1e-1*eps_values**2,
        linestyle="--", label=r"reference $\varepsilon^2$"
    )
    plt.xlabel(r"$\varepsilon=k_A/k_F$")
    plt.ylabel("conditional FIM eigenvalue")
    plt.legend()
    plt.tight_layout()
    plt.savefig("phaseII_reversible_exchange_scaling.png", dpi=200)

    for name, arr in spectra.items():
        print(name, "final spectrum:", arr[-1])
        slopes = [
            np.polyfit(
                np.log(eps_values[-10:]),
                np.log(arr[-10:,j]), 1
            )[0]
            for j in range(3)
        ]
        print("  slopes:", slopes)

if __name__ == "__main__":
    main()
