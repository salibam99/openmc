"""Shared definitions for the 3D eigenvalue toy (units: cm, 1/cm)."""
# same one-group materials as the 1D toy (../model.py)
FUEL = dict(st=1.0, ss=0.6, sf=0.22, nu=2.4)
REFL = dict(st=1.5, ss=1.4, sf=0.0, nu=0.0)

CASES = {
    # 1: bare fuel sphere (k ~ 1.000)
    1: dict(R=4.2),
    # 2: fuel sphere in a concentric reflector shell (k ~ 1.023)
    2: dict(a=3.5, b=5.5),
    # 3: off-centre fuel sphere (centre x0, 0, 0) in a reflector cylinder
    #    x^2 + y^2 < c^2, zlo < z < zhi (planes symmetric about z = 0)
    3: dict(a=3.5, x0=1.0, c=5.5, zlo=-6.0, zhi=6.0),
}
