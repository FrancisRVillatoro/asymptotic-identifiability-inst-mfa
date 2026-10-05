#!/usr/bin/env python3
"""
Phase I audit for the reversible two-metabolite KFP model of
Guppy, Mitchell & Taylor (Bull. Math. Biol. 87, 7, 2025).

Model:
  x1dot = k1(-x1 + beta*x2 + alpha1)
  x2dot = k2((1-alpha2)*x1 - x2 + alpha2)
  x1(0)=x2(0)=1

Parameter coordinates used for the Fisher information matrix:
  theta = (log k1, alpha1, beta12, log k2, alpha2)

Audit design:
  - baseline k1=1/25, alpha1=1/4, beta12=3/20, alpha2=3/10
  - epsilon=k1/k2 varied
  - 10 samples equally spaced in slow time eta=k1*t over [0.25,5]
  - both x1 and x2 observed
  - 3 replicates
  - 5% relative Gaussian noise; weights frozen at the true trajectory

The sensitivities are computed using the Frechet derivative of the matrix
exponential, avoiding finite-difference errors in the stiff limit.
"""

import numpy as np
from scipy.linalg import expm, expm_frechet, eigh
import csv
import matplotlib.pyplot as plt

K1 = 1/25
ALPHA1 = 1/4
BETA = 3/20
ALPHA2 = 3/10
ETAS = np.linspace(0.25, 5.0, 10)
REL_NOISE = 0.05
NREP = 3
QINIT = np.array([1.0, 1.0, 1.0])

def augmented_matrix(theta):
    lk1, a1, beta, lk2, a2 = theta
    k1 = np.exp(lk1)
    k2 = np.exp(lk2)
    A = np.zeros((3,3))
    A[0,0] = -k1
    A[0,1] = k1*beta
    A[0,2] = k1*a1
    A[1,0] = k2*(1-a2)
    A[1,1] = -k2
    A[1,2] = k2*a2
    return A

def augmented_derivatives(theta):
    lk1, a1, beta, lk2, a2 = theta
    k1 = np.exp(lk1)
    k2 = np.exp(lk2)
    D = []

    X = np.zeros((3,3))
    X[0,0], X[0,1], X[0,2] = -k1, k1*beta, k1*a1
    D.append(X)

    X = np.zeros((3,3)); X[0,2] = k1
    D.append(X)

    X = np.zeros((3,3)); X[0,1] = k1
    D.append(X)

    X = np.zeros((3,3))
    X[1,0], X[1,1], X[1,2] = k2*(1-a2), -k2, k2*a2
    D.append(X)

    X = np.zeros((3,3)); X[1,0], X[1,2] = -k2, k2
    D.append(X)

    return D

def state_and_sensitivities(theta, t):
    A = augmented_matrix(theta)
    E = expm(A*t)
    q = E @ QINIT
    S = np.empty((2,5))
    for j, D in enumerate(augmented_derivatives(theta)):
        L = expm_frechet(A*t, D*t, compute_expm=False)
        S[:,j] = (L @ QINIT)[:2]
    return q[:2], S

def fisher(theta):
    k1 = np.exp(theta[0])
    F = np.zeros((5,5))
    for eta in ETAS:
        t = eta/k1
        y, S = state_and_sensitivities(theta, t)
        W = np.diag(NREP/(REL_NOISE*y)**2)
        F += S.T @ W @ S
    return F

def reduced_null_direction():
    c = 1-ALPHA2
    r = 1-BETA*c
    B = ALPHA1+BETA*ALPHA2
    xs = B/r
    ys = c*xs+ALPHA2
    q = np.array([c/r, -ys, 1.0, 0.0, 0.0])
    return q/np.linalg.norm(q)

def slow_manifold(eps):
    c = 1-ALPHA2
    r = 1-BETA*c
    B = ALPHA1+BETA*ALPHA2
    xs = B/r
    ys = c*xs+ALPHA2
    disc = (1-eps)**2 + 4*eps*BETA*c
    ms = (-(1-eps)+np.sqrt(disc))/(2*eps*BETA)
    mismatch = 1-(ys+ms*(1-xs))
    return xs, ys, ms, mismatch

def main():
    epsvals = np.logspace(-1, -3, 17)
    rows = []
    eigvecs = {}
    for eps in epsvals:
        k2 = K1/eps
        theta = np.array([np.log(K1), ALPHA1, BETA, np.log(k2), ALPHA2])
        F = fisher(theta)
        w, V = eigh(F)
        rows.append([eps, *w])
        if np.isclose(eps, 0.04, rtol=0.15):
            eigvecs[eps] = V

    with open("phaseI_guppy_fim_spectrum.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["epsilon","lambda1","lambda2","lambda3","lambda4","lambda5"])
        writer.writerows(rows)

    arr = np.array(rows)
    plt.figure(figsize=(7,5))
    for j in range(1,6):
        plt.loglog(arr[:,0], arr[:,j], marker="o", label=f"lambda{j}")
    # reference slopes, normalized visually
    ref = arr[:,0]
    plt.loglog(ref, 1.0e2*ref**4, linestyle="--", label="slope 4")
    plt.loglog(ref, 2.5e2*ref**2, linestyle="--", label="slope 2")
    plt.xlabel(r"$\varepsilon=k_1/k_2$")
    plt.ylabel("FIM eigenvalue")
    plt.legend()
    plt.tight_layout()
    plt.savefig("phaseI_guppy_fim_scaling.png", dpi=200)

    eps = 0.04
    theta = np.array([np.log(K1), ALPHA1, BETA, np.log(K1/eps), ALPHA2])
    F = fisher(theta)
    w, V = eigh(F)
    print("FIM at epsilon=0.04:")
    print(F)
    print("eigenvalues:", w)
    print("reduced null direction:", reduced_null_direction())
    print("slow manifold quantities:", slow_manifold(eps))

if __name__ == "__main__":
    main()
