#!/usr/bin/env python3
"""
Phase II: coupled fast EMU block.

Slow precursor labeling:
    A(a), B(b), q=ab.

Fast isotopomer pool X(ab):
    x10, x01, x11.

Fast downstream one-carbon pool Y:
    y receives carbon a or b from X.

Parameters:
    theta = (log kA, log kB, log kX, log kY, rho, eta, alpha).

Experiments:
    A-only, B-only, and A+B labeling.

Observations:
    XY: X M+1=x10+x01, X M+2=x11, and Y labeled fraction.
    Y : Y labeled fraction only.

The script reports:
  1. the conditional 2x2 Fisher spectrum for (log kX, log kY)
     with the remaining parameters fixed;
  2. the full 7x7 Fisher spectrum for the XY design.

All local sensitivities are evaluated by complex-step differentiation
of an exact augmented matrix exponential.
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm, svdvals

K_A = 0.04
K_B = 0.07
RHO = 0.80
ETA = 0.70
ALPHA = 0.30
CHI = 1.70                  # kY/kX
TAU = np.linspace(0.10, 5.0, 15)
TIMES = TAU / K_A
EXPERIMENTS = ((1,0), (0,1), (1,1))
Q0 = np.zeros(8)
Q0[-1] = 1.0

def augmented_matrix(theta, experiment):
    lkA, lkB, lkX, lkY, rho, eta, alpha = theta
    kA, kB, kX, kY = np.exp(lkA), np.exp(lkB), np.exp(lkX), np.exp(lkY)
    sA, sB = experiment

    # state = (a, b, q=ab, x10, x01, x11, y, 1)
    M = np.zeros((8,8), dtype=np.result_type(theta))
    M[0,0], M[0,7] = -kA, kA*sA
    M[1,1], M[1,7] = -kB, kB*sB
    M[2,0] = kB*sB
    M[2,1] = kA*sA
    M[2,2] = -(kA+kB)

    M[3,0] = kX*rho
    M[3,2] = -kX*rho
    M[3,3] = -kX

    M[4,1] = kX*rho
    M[4,2] = -kX*rho
    M[4,4] = -kX

    M[5,2] = kX*rho
    M[5,5] = -kX

    M[6,3] = kY*eta*alpha
    M[6,4] = kY*eta*(1-alpha)
    M[6,5] = kY*eta
    M[6,6] = -kY
    return M

def observed_response(theta, t, observation):
    values = []
    for experiment in EXPERIMENTS:
        q = expm(augmented_matrix(theta, experiment)*t) @ Q0
        if observation == "XY":
            values.extend((q[3]+q[4], q[5], q[6]))
        elif observation == "Y":
            values.append(q[6])
        else:
            raise ValueError(observation)
    return np.asarray(values)

def jacobian(theta, t, observation):
    h = 1e-28
    y = observed_response(theta, t, observation)
    J = np.zeros((len(y), len(theta)))
    for j in range(len(theta)):
        z = theta.astype(complex)
        z[j] += 1j*h
        J[:,j] = np.imag(observed_response(z, t, observation))/h
    return J

def stacked_jacobian(theta, observation):
    return np.vstack([jacobian(theta, t, observation) for t in TIMES])

def spectra(epsilon, observation):
    kX = K_A/epsilon
    kY = CHI*kX
    theta = np.array([
        np.log(K_A), np.log(K_B), np.log(kX), np.log(kY),
        RHO, ETA, ALPHA
    ])
    J = stacked_jacobian(theta, observation)
    conditional = np.sort(svdvals(J[:,2:4])**2)
    full = np.sort(svdvals(J)**2)
    return conditional, full

def main():
    eps_values = np.logspace(-1, -4, 25)
    rows = []
    cond = {}
    full = {}

    for observation in ("Y", "XY"):
        cvals, fvals = [], []
        for eps in eps_values:
            c, f = spectra(eps, observation)
            cvals.append(c)
            fvals.append(f)
            rows.append([eps, observation, *c, *f])
        cond[observation] = np.asarray(cvals)
        full[observation] = np.asarray(fvals)

    with open("phaseII_coupled_fast_spectrum.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "epsilon", "observation",
            "conditional_lambda1", "conditional_lambda2",
            "full_lambda1", "full_lambda2", "full_lambda3",
            "full_lambda4", "full_lambda5", "full_lambda6",
            "full_lambda7"
        ])
        writer.writerows(rows)

    plt.figure(figsize=(7.3,5.3))
    plt.loglog(
        eps_values, cond["Y"][:,0], marker="o",
        label="Y only: smallest"
    )
    plt.loglog(
        eps_values, cond["Y"][:,1], marker="o", linestyle="--",
        label="Y only: second"
    )
    plt.loglog(
        eps_values, cond["XY"][:,0], marker="s",
        label="X+Y: smallest"
    )
    plt.loglog(
        eps_values, cond["XY"][:,1], marker="s", linestyle="--",
        label="X+Y: second"
    )
    plt.loglog(
        eps_values, 1e-3*eps_values**4, linestyle=":",
        label=r"reference $\varepsilon^4$"
    )
    plt.loglog(
        eps_values, 1e-1*eps_values**2, linestyle="-.",
        label=r"reference $\varepsilon^2$"
    )
    plt.xlabel(r"$\varepsilon=k_A/k_X$")
    plt.ylabel("conditional FIM eigenvalue")
    plt.legend()
    plt.tight_layout()
    plt.savefig("phaseII_coupled_fast_scaling.png", dpi=200)

    for observation in ("Y", "XY"):
        print(observation)
        print("  final conditional spectrum:", cond[observation][-1])
        for j in range(2):
            slope = np.polyfit(
                np.log(eps_values[-10:]),
                np.log(cond[observation][-10:,j]), 1
            )[0]
            print(f"  conditional lambda{j+1} slope = {slope:.6f}")

    print("XY full spectrum at epsilon=1e-4:")
    print(full["XY"][-1])

if __name__ == "__main__":
    main()
