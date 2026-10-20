"""
S_N reference for the FINITE-L expectation of variant (a) (and of LR).

For fixed L, variant (a) estimates (derivation.md, section 4)

    Q_L(theta) = <psi_L+, H_theta psi> = <1, H_0^L H_theta psi>,
    E[(a), L]  = d ln Q_L / d theta,

i.e. only ONE generation's kernel is perturbed, and the source psi and the L-step
importance are frozen as functions of PHYSICAL position. This is independent of
the reparameterization, so S_N can compute it. As L -> inf it tends to d ln k.

To keep psi and the importance fixed in physical space, use a fixed uniform mesh
on [0, B0 + hmax], with void (sigma = 0) beyond b, and move a or b by whole cells.
"""
import json
from pathlib import Path

import numpy as np
import scipy.linalg as sla

from model import FUEL, REFL, A0, B0
from sn import response_matrix

LMAX = 10


def operator(dx, M, a, b, rho_f=1.0, rho_r=1.0, nquad=64):
    """Dense next-generation operator H: emission density -> fission-site density."""
    mu, wt = np.polynomial.legendre.leggauss(nquad)
    pos = mu > 0
    mu, wt = mu[pos].copy(), wt[pos].copy()
    xc = (np.arange(M) + 0.5) * dx
    fuel = xc < a
    refl = (xc > a) & (xc < b)
    st = np.where(fuel, rho_f * FUEL["st"], np.where(refl, rho_r * REFL["st"], 0.0))
    ss = np.where(fuel, rho_f * FUEL["ss"], np.where(refl, rho_r * REFL["ss"], 0.0))
    nsf = np.where(fuel, rho_f * FUEL["nu"] * FUEL["sf"], 0.0)
    h = np.full(M, dx)
    T = response_matrix(h, st, mu, wt)
    A = np.eye(M) - T * ss[None, :]
    return nsf[:, None] * sla.solve(A, T, check_finite=False)


def frozen(dx, nquad, steps=(2, 4)):
    hmax = max(steps) * dx
    M = int(round((B0 + hmax) / dx))
    assert abs(A0 / dx - round(A0 / dx)) < 1e-9 and abs(B0 / dx - round(B0 / dx)) < 1e-9
    H0 = operator(dx, M, A0, B0, nquad=nquad)
    # fundamental fission source
    psi = np.ones(M)
    for _ in range(400):
        p2 = H0 @ psi
        k = p2.sum() / psi.sum()
        psi = p2 / k
    # importance vectors v_L = (H0^T)^L 1   (cell volumes are uniform)
    v = [np.ones(M)]
    for _ in range(LMAX):
        v.append(H0.T @ v[-1])
    v = np.array(v)
    Q0 = v @ (H0 @ psi)
    res = {"k": k}
    for name in ("a", "b", "rho_f", "rho_r"):
        est = []
        for m in steps:
            hh = m * dx if name in ("a", "b") else 0.01 * m / steps[0]
            kw_p = dict(a=A0, b=B0)
            kw_m = dict(a=A0, b=B0)
            if name in ("a", "b"):
                kw_p[name] += hh
                kw_m[name] -= hh
            else:
                kw_p[name] = 1 + hh
                kw_m[name] = 1 - hh
            Qp = v @ (operator(dx, M, nquad=nquad, **kw_p) @ psi)
            Qm = v @ (operator(dx, M, nquad=nquad, **kw_m) @ psi)
            est.append((Qp - Qm) / (2 * hh * Q0))
        # O(h^2) with step ratio 2 -> Richardson
        res["dlnQ_d" + name] = ((4 * est[0] - est[1]) / 3).tolist()
        res["dlnQ_d" + name + "_raw"] = [e.tolist() for e in est]
    return res


if __name__ == "__main__":
    out = {}
    for dx, nq in ((0.01, 32), (0.005, 32), (0.005, 64), (0.0025, 64), (0.005, 128)):
        r = frozen(dx, nq)
        out[f"dx{dx}_S{nq}"] = r
        print(f"dx={dx} S{nq}: k={r['k']:.7f}")
        for name in ("a", "b", "rho_f", "rho_r"):
            print(f"  d lnQ_L/d{name:5s} L=0..10: " + " ".join(f"{x:+.5f}" for x in r["dlnQ_d" + name]))
    # best estimate: angular Richardson (O(1/N^2)) on S64/S128 at dx = 0.005,
    # plus the mesh correction dx 0.005 -> 0.0025 at S64
    best = {}
    for name in ("a", "b", "rho_f", "rho_r"):
        key = "dlnQ_d" + name
        lo, hi = np.array(out["dx0.005_S64"][key]), np.array(out["dx0.005_S128"][key])
        mesh = np.array(out["dx0.0025_S64"][key]) - lo
        best[key] = ((4 * hi - lo) / 3 + mesh).tolist()
        print(f"best {name:5s}: " + " ".join(f"{x:+.5f}" for x in best[key]))
    out["best"] = best
    json.dump(out, open(Path(__file__).parent / "sn_frozen_results.json", "w"), indent=1)
