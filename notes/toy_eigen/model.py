"""Shared model definition for the eigenvalue toy (units: cm, 1/cm)."""
# fuel on [0, a], reflector on [a, b], vacuum at x = 0 and x = b
# capture = st - ss - sf;  k ~ 1.019, dominance ratio ~ 0.57 (S_N, see results.md)
FUEL = dict(st=1.0, ss=0.6, sf=0.22, nu=2.4)
REFL = dict(st=1.5, ss=1.4, sf=0.0, nu=0.0)
A0 = 3.0
B0 = 4.5
