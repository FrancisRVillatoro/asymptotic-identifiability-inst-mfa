#!/usr/bin/env python3
"""
Minimal two-carbon EMU model for Phase II.

Two input tracer mixtures have identical mass-isotopomer distributions:
  T1 = (M00, M10, M01, M11) = (1-e-u, e, 0, u)
  T2 = (M00, M10, M01, M11) = (1-e-u, 0, e, u)

A two-carbon pool X is supplied by tracer and unlabeled dilution.
A fast one-carbon pool Y receives carbon 1 or carbon 2 of X through two routes.

Parameters:
  theta = (log kX, rho, log kY, eta, alpha)

Observations at each time:
  X M+1, X M+2, Y labeled fraction.

The script computes exact local sensitivities by complex-step differentiation
and compares one-tracer and two-tracer Fisher spectra as epsilon=kX/kY -> 0.
"""

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

KX = 0.04
RHO = 0.80
ETA = 0.70
ALPHA = 0.30
E = 0.60
U = 0.20
TAU = np.linspace(0.10, 5.0, 12)
TIMES = TAU / KX

def phi(t, kx, ky):
    if abs(ky-kx) < 1e-10:
        return 1.0 - (1.0 + kx*t)*np.exp(-kx*t)
    return 1.0 - (
        ky*np.exp(-kx*t) - kx*np.exp(-ky*t)
    )/(ky-kx)

def response(theta, t, tracers=(1,2)):
    lkx, rho, lky, eta, alpha = theta
    kx, ky = np.exp(lkx), np.exp(lky)

    gx = 1.0 - np.exp(-kx*t)
    m1 = rho*E*gx
    m2 = rho*U*gx
    ph = phi(t, kx, ky)

    output = []
    for tracer in tracers:
        branch = alpha if tracer == 1 else (1.0-alpha)
        amplitude = rho*eta*(U + branch*E)
        output.extend((m1, m2, amplitude*ph))
    return np.asarray(output)

def jacobian_complex_step(theta, t, tracers=(1,2)):
    h = 1e-28
    y = response(theta, t, tracers)
    J = np.zeros((len(y), len(theta)))
    for j in range(len(theta)):
        perturbed = theta.astype(complex)
        perturbed[j] += 1j*h
        J[:,j] = np.imag(response(perturbed, t, tracers))/h
    return y, J

def stacked_jacobian(theta, tracers):
    return np.vstack([
        jacobian_complex_step(theta, t, tracers)[1]
        for t in TIMES
    ])

def spectrum(epsilon, tracers):
    ky = KX/epsilon
    theta = np.array([
        np.log(KX), RHO, np.log(ky), ETA, ALPHA
    ])
    J = stacked_jacobian(theta, tracers)
    singular_values = np.linalg.svd(J, compute_uv=False)
    return np.sort(singular_values**2)

def main():
    eps_values = np.logspace(-1, -4, 25)
    rows = []
    spectra = {}

    for name, tracers in (
        ("tracer_1", (1,)),
        ("tracers_1_2", (1,2)),
    ):
        values = []
        for eps in eps_values:
            eig = spectrum(eps, tracers)
            values.append(eig)
            rows.append([eps, name, *eig])
        spectra[name] = np.asarray(values)

    with open("phaseII_minimal_emu_spectrum.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "epsilon", "design",
            "lambda1", "lambda2", "lambda3", "lambda4", "lambda5"
        ])
        writer.writerows(rows)

    plt.figure(figsize=(7.2, 5.2))
    one = spectra["tracer_1"]
    two = spectra["tracers_1_2"]

    # The single-tracer smallest eigenvalue is structurally zero, so plot
    # its smallest nonzero eigenvalue.
    plt.loglog(
        eps_values, one[:,1], marker="o",
        label="one tracer: smallest nonzero"
    )
    plt.loglog(
        eps_values, two[:,0], marker="s",
        label="two tracers: smallest"
    )
    plt.loglog(
        eps_values, 0.1*eps_values**2, linestyle="--",
        label=r"reference $\varepsilon^2$"
    )
    plt.xlabel(r"$\varepsilon=k_X/k_Y$")
    plt.ylabel("FIM eigenvalue")
    plt.legend()
    plt.tight_layout()
    plt.savefig("phaseII_minimal_emu_scaling.png", dpi=200)

    for name, values in spectra.items():
        print(name)
        print("last spectrum:", values[-1])
        for j in range(5):
            positive = values[-10:,j] > 1e-25
            if positive.sum() >= 4:
                slope = np.polyfit(
                    np.log(eps_values[-10:][positive]),
                    np.log(values[-10:,j][positive]), 1
                )[0]
                print(f"  lambda{j+1} slope = {slope:.6f}")
            else:
                print(f"  lambda{j+1}: structurally zero")

if __name__ == "__main__":
    main()
