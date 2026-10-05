#!/usr/bin/env python3
"""Exact symbolic derivation of the OAA-Fum weak direction."""
import sympy as sp

v1   = sp.Rational(10)
v6f  = sp.Rational(25,2)
v6b  = sp.Rational(15,2)
v7   = sp.Rational(5)
v5   = sp.Rational(5)
CO   = sp.Rational(1,10)
CF   = sp.Rational(1,5)

tauO = sp.simplify(CO/(v1+v6b))
tauF = sp.simplify(CF/v6f)
alpha = sp.simplify(v6f/(v1+v6b))
beta  = sp.simplify(v6b/v6f)
Delta = sp.simplify(1-alpha*beta)

D = sp.symbols("D")
den = sp.expand((1+tauO*D)*(1+tauF*D)-alpha*beta)

q = sp.Matrix([14,-5])
M1 = sp.Matrix([[tauO,tauF]])

print("alpha =", alpha)
print("beta =", beta)
print("Delta =", Delta)
print("tauO =", tauO)
print("tauF =", tauF)
print("tauO/tauF =", sp.simplify(tauO/tauF))
print("denominator =", den)
print("M1*q =", sp.simplify(M1*q))
print("delta log(tauO*tauF) along q =", q[0]+q[1])
print(
    "second-order coefficient =",
    sp.simplify(-9*tauO*tauF/Delta)
)
