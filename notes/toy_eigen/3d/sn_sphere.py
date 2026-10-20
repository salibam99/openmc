"""
Deterministic reference for 3D cases 1-2: 1D spherical S_N, one group, isotropic
scattering, vacuum outside.

Discretization (Lewis & Miller): diamond difference in r and in mu, Gauss-Legendre
directions sorted by mu, angular-redistribution coefficients
alpha_{n+1/2} = alpha_{n-1/2} - w_n mu_n, starting direction mu = -1, mirror
condition at r = 0. No negative-flux fixup, so the discrete operator stays linear
and k_h(theta) smooth. Shells are uniform within each region and stretch with the
radii (fixed number per region).

As in ../sn.py, the scattering problem is solved exactly with a dense response
matrix, and the fission-source matrix on fuel shells is power-iterated.

Independent check (bare sphere only): for a homogeneous sphere, r*phi(r) solves
the slab equation on [-R, R] with odd symmetry. So k_sphere(R) is the largest
odd-mode eigenvalue of the slab [-R, R] (../sn.py response matrix).
"""
import json
import sys
from pathlib import Path

import numba as nb
import numpy as np
import scipy.linalg as sla

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sn import response_matrix as slab_response  # noqa: E402

from geom3d import FUEL, REFL, CASES  # noqa: E402


@nb.njit(parallel=True, cache=True)
def response_matrix_sph(r, sig, mu, wt):
    """T[:, j]: scalar flux from unit isotropic emission density in shell j."""
    M = sig.size
    N = mu.size
    A = 4.0 * np.pi * r ** 2
    V = 4.0 * np.pi / 3.0 * (r[1:] ** 3 - r[:-1] ** 3)
    alpha = np.zeros(N + 1)
    for n in range(N):
        alpha[n + 1] = alpha[n] - wt[n] * mu[n]
    T = np.zeros((M, M))
    for j in nb.prange(M):
        col = np.zeros(M)
        half = np.zeros(M)          # psi_{n-1/2, i}
        center = np.zeros(N)        # psi_{n, 1/2} for mu_n < 0
        # starting direction mu = -1: slab-like DD along the diameter, inward
        pin = 0.0
        for i in range(M - 1, -1, -1):
            h = r[i + 1] - r[i]
            s = 0.5 if i == j else 0.0
            pc = (s + 2.0 * pin / h) / (sig[i] + 2.0 / h)
            half[i] = pc
            pin = 2.0 * pc - pin
        for n in range(N):
            m = mu[n]
            ap = alpha[n + 1]
            am = alpha[n]
            if m < 0.0:
                pin = 0.0                                   # vacuum at r = R
                for i in range(M - 1, -1, -1):
                    s = 0.5 if i == j else 0.0
                    dA = A[i + 1] - A[i]
                    num = s * V[i] - m * (A[i + 1] + A[i]) * pin + dA / wt[n] * (ap + am) * half[i]
                    den = -2.0 * m * A[i] + 2.0 * dA * ap / wt[n] + sig[i] * V[i]
                    pc = num / den
                    pin = 2.0 * pc - pin
                    half[i] = 2.0 * pc - half[i]
                    col[i] += wt[n] * pc
                center[n] = pin
            else:
                pin = center[N - 1 - n]                     # mirror at r = 0
                for i in range(M):
                    s = 0.5 if i == j else 0.0
                    dA = A[i + 1] - A[i]
                    num = s * V[i] + m * (A[i + 1] + A[i]) * pin + dA / wt[n] * (ap + am) * half[i]
                    den = 2.0 * m * A[i + 1] + 2.0 * dA * ap / wt[n] + sig[i] * V[i]
                    pc = num / den
                    pin = 2.0 * pc - pin
                    half[i] = 2.0 * pc - half[i]
                    col[i] += wt[n] * pc
        T[:, j] = col
    return T


def solve(radii, mats, rho=None, ncell=400, nquad=64):
    """k for concentric shells radii[0] < radii[1] < ...; mats[i] in {FUEL, REFL}."""
    rho = rho if rho is not None else [1.0] * len(radii)
    edges = [0.0]
    st, ss, nsf = [], [], []
    for i, R in enumerate(radii):
        e = np.linspace(edges[-1], R, ncell + 1)[1:]
        edges.extend(e)
        st += [rho[i] * mats[i]["st"]] * ncell
        ss += [rho[i] * mats[i]["ss"]] * ncell
        nsf += [rho[i] * mats[i]["nu"] * mats[i]["sf"]] * ncell
    r = np.array(edges)
    st, ss, nsf = np.array(st), np.array(ss), np.array(nsf)
    mu, wt = np.polynomial.legendre.leggauss(nquad)
    T = response_matrix_sph(r, st, mu, wt)
    Amat = np.eye(st.size) - T * ss[None, :]
    fis = np.nonzero(nsf > 0)[0]
    X = sla.solve(Amat, T[:, fis], check_finite=False)
    G = nsf[fis][:, None] * X[fis, :]
    vol = (r[1:] ** 3 - r[:-1] ** 3)[fis]
    s = np.ones(fis.size)
    k = 1.0
    for it in range(4000):
        s2 = G @ s
        k2 = np.dot(s2, vol) / np.dot(s, vol)
        s = s2 / k2
        if abs(k2 - k) < 1e-15 * k2 and it > 20:
            k = k2
            break
        k = k2
    return k


def slab_odd_k(R, st, ss, nsf, ncell=2000, nquad=64):
    """Largest odd-mode eigenvalue of the homogeneous slab [-R, R]."""
    mu, wt = np.polynomial.legendre.leggauss(nquad)
    pos = mu > 0
    mu, wt = mu[pos].copy(), wt[pos].copy()
    M = 2 * ncell
    h = np.full(M, 2 * R / M)
    sig = np.full(M, st)
    T = slab_response(h, sig, mu, wt)
    X = sla.solve(np.eye(M) - T * ss, T, check_finite=False)
    G = nsf * X
    # odd subspace on the right half: s(-x) = -s(x)
    right = np.arange(ncell, M)
    left = M - 1 - right              # mirror image of each right-half cell
    Godd = G[np.ix_(right, right)] - G[np.ix_(right, left)]
    return np.max(np.linalg.eigvals(Godd).real)


def case_k(case, ncell, nquad, **pert):
    c = dict(CASES[case])
    c.update(pert)
    if case == 1:
        return solve([c["R"]], [FUEL], [c.get("rho_f", 1.0)], ncell, nquad)
    return solve([c["a"], c["b"]], [FUEL, REFL], [c.get("rho_f", 1.0), c.get("rho_r", 1.0)],
                 ncell, nquad)


def derivs(case, ncell, nquad, rel=(2e-3, 1e-3)):
    params = ["R", "rho_f"] if case == 1 else ["a", "b", "rho_f", "rho_r"]
    out = {"k": case_k(case, ncell, nquad)}
    for p in params:
        base = CASES[case].get(p, 1.0)
        est = []
        for rr in rel:
            h = rr * base
            est.append((case_k(case, ncell, nquad, **{p: base + h})
                        - case_k(case, ncell, nquad, **{p: base - h})) / (2 * h))
        out["dk_d" + p] = (4 * est[1] - est[0]) / 3
    return out


if __name__ == "__main__":
    import time
    t0 = time.time()
    if len(sys.argv) > 1 and sys.argv[1] == "check":
        # spherical S_N vs slab odd-mode equivalence (bare sphere)
        R = CASES[1]["R"]
        for nq in (16, 32, 64):
            print(f"S{nq}: sphere k = {solve([R], [FUEL], ncell=400, nquad=nq):.8f}   "
                  f"slab-odd k = {slab_odd_k(R, FUEL['st'], FUEL['ss'], FUEL['nu'] * FUEL['sf'], nquad=nq):.8f}",
                  flush=True)
        sys.exit()
    res = {}
    for case in (1, 2):
        for ncell, nq in ((200, 32), (400, 32), (800, 32), (400, 64), (400, 128), (800, 128), (800, 256)):
            d = derivs(case, ncell, nq)
            res[f"case{case}_n{ncell}_S{nq}"] = d
            print(f"case {case} cells/region={ncell} S{nq}: " +
                  "  ".join(f"{k_}={v:+.8f}" for k_, v in d.items()) + f"  [{time.time()-t0:.0f}s]",
                  flush=True)
        # Richardson in angle (O(1/N^2), checked in the log) at 800 cells
        lo, hi = res[f"case{case}_n800_S128"], res[f"case{case}_n800_S256"]
        ref = {k_: (4 * hi[k_] - lo[k_]) / 3 for k_ in hi}
        for k_ in list(ref):
            if k_.startswith("dk_d"):
                ref["dlnk_d" + k_[4:]] = ref[k_] / ref["k"]
        res[f"case{case}_reference"] = ref
        print(f"case {case} reference:", json.dumps(ref), flush=True)
    json.dump(res, open(Path(__file__).parent / "sn_sphere_results.json", "w"), indent=1)
