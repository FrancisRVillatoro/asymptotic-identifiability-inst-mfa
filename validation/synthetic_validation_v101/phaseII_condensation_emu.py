#!/usr/bin/env python3
"""
Phase II: minimal condensation EMU audit.

Network:
    A(a) + B(b) -> X(ab)

The labeled fractions a(t), b(t) satisfy independent step-response ODEs.
The formation MID of X is the Bernoulli convolution
    M+0 = (1-a)(1-b)
    M+1 = a(1-b) + (1-a)b
    M+2 = ab
and X has turnover kX.

Parameters:
    theta = (log kA, log kB, log kX, rho)

Designs:
    A       : label A only
    B       : label B only
    AB      : label A and B simultaneously
    A+B     : combine two separate single-precursor experiments

Exact local sensitivities are evaluated by complex-step differentiation.
The reported Fisher matrix is J.T @ J, corresponding to equal fixed
observation variance.
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import svdvals

K_A = 0.04
K_B = 0.07
RHO = 0.80
TAU = np.linspace(0.10, 5.0, 15)
TIMES = TAU / K_A

def lowpass_exp(t, kx, lam):
    """kX convolution of exp(-lam*t), zero initial output."""
    if abs(kx-lam) < 1e-10:
        return kx*t*np.exp(-kx*t)
    return kx*(np.exp(-lam*t)-np.exp(-kx*t))/(kx-lam)

def response(theta, t, design, observe_m2=True):
    kA, kB, kX = np.exp(theta[0]), np.exp(theta[1]), np.exp(theta[2])
    rho = theta[3]
    h0 = 1.0 - np.exp(-kX*t)
    LA = lowpass_exp(t, kX, kA)
    LB = lowpass_exp(t, kX, kB)
    LAB = lowpass_exp(t, kX, kA+kB)

    values = []
    for experiment in design:
        if experiment == "A":
            m1 = rho*(h0-LA)
            m2 = 0.0*m1
        elif experiment == "B":
            m1 = rho*(h0-LB)
            m2 = 0.0*m1
        elif experiment == "AB":
            m1 = rho*(LA+LB-2.0*LAB)
            m2 = rho*(h0-LA-LB+LAB)
        else:
            raise ValueError(f"Unknown experiment: {experiment}")
        values.append(m1)
        if observe_m2:
            values.append(m2)
    return np.asarray(values)

def jacobian(theta, t, design, observe_m2=True):
    h = 1e-28
    y = response(theta, t, design, observe_m2)
    J = np.zeros((len(y), len(theta)))
    for j in range(len(theta)):
        z = theta.astype(complex)
        z[j] += 1j*h
        J[:, j] = np.imag(response(z, t, design, observe_m2))/h
    return J

def spectrum(epsilon, design, observe_m2=True, kB=K_B):
    kX = K_A/epsilon
    theta = np.array([np.log(K_A), np.log(kB), np.log(kX), RHO])
    J = np.vstack([
        jacobian(theta, t, design, observe_m2)
        for t in TIMES
    ])
    s = svdvals(J)
    return np.sort(s*s)

def reduced_response(theta, t, design):
    kA, kB, rho = np.exp(theta[0]), np.exp(theta[1]), theta[2]
    a = 1.0-np.exp(-kA*t)
    b = 1.0-np.exp(-kB*t)
    values = []
    for experiment in design:
        if experiment == "A":
            values.extend([rho*a, 0.0*a])
        elif experiment == "B":
            values.extend([rho*b, 0.0*b])
        elif experiment == "AB":
            values.extend([rho*(a+b-2*a*b), rho*a*b])
    return np.asarray(values)

def reduced_jacobian(theta, t, design):
    h = 1e-28
    y = reduced_response(theta, t, design)
    J = np.zeros((len(y), 3))
    for j in range(3):
        z = theta.astype(complex)
        z[j] += 1j*h
        J[:,j] = np.imag(reduced_response(z, t, design))/h
    return J

def reduced_spectrum(kB, design):
    theta = np.array([np.log(K_A), np.log(kB), RHO])
    J = np.vstack([reduced_jacobian(theta, t, design) for t in TIMES])
    s = svdvals(J)
    return np.sort(s*s)

def main():
    eps_values = np.logspace(-1, -4, 25)
    configurations = [
        ("A_only", ("A",), True),
        ("A_plus_B", ("A","B"), True),
        ("AB_M1_M2", ("AB",), True),
        ("AB_M1_only", ("AB",), False),
        ("A_B_AB", ("A","B","AB"), True),
    ]

    rows = []
    spectra = {}
    for name, design, observe_m2 in configurations:
        values = []
        for eps in eps_values:
            eig = spectrum(eps, design, observe_m2)
            values.append(eig)
            rows.append([eps, name, *eig])
        spectra[name] = np.asarray(values)

    with open("phaseII_condensation_spectrum.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "epsilon", "design",
            "lambda1", "lambda2", "lambda3", "lambda4"
        ])
        writer.writerows(rows)

    plt.figure(figsize=(7.2,5.2))
    for name, marker in [
        ("A_plus_B","o"),
        ("AB_M1_M2","s"),
        ("AB_M1_only","^"),
    ]:
        plt.loglog(
            eps_values, spectra[name][:,0],
            marker=marker, label=name
        )
    plt.loglog(
        eps_values, 1.0*eps_values**2,
        linestyle="--", label=r"$\propto\varepsilon^2$"
    )
    plt.xlabel(r"$\varepsilon=k_A/k_X$")
    plt.ylabel("smallest FIM eigenvalue")
    plt.legend()
    plt.tight_layout()
    plt.savefig("phaseII_condensation_scaling.png", dpi=200)

    deltas = np.logspace(-5, -0.3, 35)
    dual = []
    complementary = []
    for delta in deltas:
        kB = K_A*(1.0+delta)
        dual.append(reduced_spectrum(kB, ("AB",))[0])
        complementary.append(reduced_spectrum(kB, ("A","B"))[0])

    plt.figure(figsize=(7.2,5.2))
    plt.loglog(deltas, dual, marker="o", label="dual tracer AB")
    plt.loglog(
        deltas, complementary, marker="s",
        label="separate A and B experiments"
    )
    plt.loglog(
        deltas, 0.15*deltas**2,
        linestyle="--", label=r"$\propto\delta^2$"
    )
    plt.xlabel(r"$\delta=(k_B-k_A)/k_A$")
    plt.ylabel("smallest reduced-FIM eigenvalue")
    plt.legend()
    plt.tight_layout()
    plt.savefig("phaseII_condensation_symmetry.png", dpi=200)

if __name__ == "__main__":
    main()
