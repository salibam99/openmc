"""
Finite-L (frozen-importance) spherical S_N reference for cases 1-2, as
../sn_frozen.py does for the slab:

    Q_L(theta) = <1, H_0^L H_theta psi>,   target of variant (a) and of LR at window L
    d ln Q_L/d theta -> d ln k/d theta as L -> inf.

Fixed uniform radial mesh with void shells beyond the outer radius, so psi and the
importance stay fixed in physical space. Radii move by whole shells (2 and 4,
Richardson). Also prints the radial dominance ratio |k2|/k1. It is valid for
spherically symmetric perturbations only; translations excite the l = 1 modes,
which this 1D solver does not contain.
"""
import json
from pathlib import Path

import numpy as np
import scipy.linalg as sla

from geom3d import FUEL, REFL, CASES
from sn_sphere import response_matrix_sph

LMAX = 12


def operator(dr, M, a, b, rho_f=1.0, rho_r=1.0, nquad=64):
    rc = (np.arange(M) + 0.5) * dr
    fuel = rc < a
    refl = (rc > a) & (rc < b)
    st = np.where(fuel, rho_f * FUEL["st"], np.where(refl, rho_r * REFL["st"], 0.0))
    ss = np.where(fuel, rho_f * FUEL["ss"], np.where(refl, rho_r * REFL["ss"], 0.0))
    nsf = np.where(fuel, rho_f * FUEL["nu"] * FUEL["sf"], 0.0)
    r = np.arange(M + 1) * dr
    mu, wt = np.polynomial.legendre.leggauss(nquad)
    T = response_matrix_sph(r, st, mu, wt)
    A = np.eye(M) - T * ss[None, :]
    return nsf[:, None] * sla.solve(A, T, check_finite=False), r


def frozen(case, dr, nquad, steps=(2, 4)):
    if case == 1:
        a, b = CASES[1]["R"], CASES[1]["R"]
        params = ["R", "rho_f"]
    else:
        a, b = CASES[2]["a"], CASES[2]["b"]
        params = ["a", "b", "rho_f", "rho_r"]
    M = int(round((b + max(steps) * dr) / dr))
    H0, r = operator(dr, M, a, b, nquad=nquad)
    vol = (r[1:] ** 3 - r[:-1] ** 3)
    ev = np.linalg.eigvals(H0)
    ev = ev[np.argsort(-np.abs(ev))]
    psi = np.ones(M)
    for _ in range(600):
        p2 = H0 @ psi
        k = (vol @ p2) / (vol @ psi)
        psi = p2 / k
    w = [vol.copy()]
    for _ in range(LMAX):
        w.append(H0.T @ w[-1])
    w = np.array(w)
    Q0 = w @ (H0 @ psi)
    res = {"k": float(k), "dr_ratio": float(abs(ev[1]) / abs(ev[0]))}
    for p in params:
        est = []
        for m in steps:
            kw_p, kw_m = dict(a=a, b=b), dict(a=a, b=b)
            if p in ("rho_f", "rho_r"):
                h = 0.01 * m / steps[0]
                kw_p[p], kw_m[p] = 1 + h, 1 - h
            else:
                h = m * dr
                keys = ("a", "b") if p == "R" else (p,)     # bare sphere: a = b = R
                for kk in keys:
                    kw_p[kk] += h
                    kw_m[kk] -= h
            Qp = w @ (operator(dr, M, nquad=nquad, **kw_p)[0] @ psi)
            Qm = w @ (operator(dr, M, nquad=nquad, **kw_m)[0] @ psi)
            est.append((Qp - Qm) / (2 * h * Q0))
        res["dlnQ_d" + p] = ((4 * est[0] - est[1]) / 3).tolist()
    return res


if __name__ == "__main__":
    out = {}
    for case in (1, 2):
        for dr, nq in ((0.01, 64), (0.005, 64), (0.005, 128)):
            r = frozen(case, dr, nq)
            out[f"case{case}_dr{dr}_S{nq}"] = r
            print(f"case {case} dr={dr} S{nq}: k={r['k']:.6f}  radial DR={r['dr_ratio']:.4f}")
            for key, v in r.items():
                if key.startswith("dlnQ"):
                    print(f"   {key:12s} L=0..{LMAX}: " + " ".join(f"{x:+.5f}" for x in v), flush=True)
    json.dump(out, open(Path(__file__).parent / "sn_sphere_frozen_results.json", "w"), indent=1)
