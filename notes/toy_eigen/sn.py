"""
Deterministic reference: 1D slab, one group, isotropic scattering, vacuum both
sides. Gauss-Legendre S_N, diamond difference (no negative-flux fixup, so the
discrete operator stays linear and k_h(theta) smooth), uniform mesh with a fixed
number of cells per region (cells stretch with a and b).

The scattering problem is solved exactly with a dense matrix: T[:, j] is the
scalar-flux response to a unit isotropic emission density in cell j. Then the
fission-source matrix on fuel cells, G = nuSf * [(I - T Ss)^-1 T]_ff, is
diagonalized. That gives k and the dominance ratio |k2|/k1 to round-off.
"""
import numpy as np
import numba as nb
import scipy.linalg as sla

from model import FUEL, REFL, A0, B0


@nb.njit(parallel=True, cache=True)
def response_matrix(h, sig, mu, wt):
    M = h.size
    T = np.zeros((M, M))
    for j in nb.prange(M):
        col = np.zeros(M)
        for k in range(mu.size):
            m = mu[k]
            c = 2.0 * m / h
            # mu > 0, left to right, vacuum inflow; source q_j/2 in cell j
            psi = 0.0
            for i in range(M):
                s = 0.5 if i == j else 0.0
                pc = (s + c[i] * psi) / (sig[i] + c[i])
                psi = 2.0 * pc - psi
                col[i] += wt[k] * pc
            # mu < 0, right to left (symmetric quadrature)
            psi = 0.0
            for i in range(M - 1, -1, -1):
                s = 0.5 if i == j else 0.0
                pc = (s + c[i] * psi) / (sig[i] + c[i])
                psi = 2.0 * pc - psi
                col[i] += wt[k] * pc
        T[:, j] = col
    return T


def solve(a=A0, b=B0, rho_f=1.0, rho_r=1.0, nf=1000, nr=1000, nquad=32,
          want_modes=False):
    """k (and optionally dominance ratio, fission source) for the slab."""
    mu, wt = np.polynomial.legendre.leggauss(nquad)
    pos = mu > 0
    mu, wt = mu[pos].copy(), wt[pos].copy()       # each direction pair, weights sum to 1 per half
    h = np.concatenate([np.full(nf, a / nf), np.full(nr, (b - a) / nr)])
    st = np.concatenate([np.full(nf, rho_f * FUEL["st"]), np.full(nr, rho_r * REFL["st"])])
    ss = np.concatenate([np.full(nf, rho_f * FUEL["ss"]), np.full(nr, rho_r * REFL["ss"])])
    nsf = rho_f * FUEL["nu"] * FUEL["sf"]
    T = response_matrix(h, st, mu, wt)
    A = np.eye(h.size) - T * ss[None, :]
    X = sla.solve(A, T[:, :nf], check_finite=False)     # flux response to fuel emission
    G = nsf * X[:nf, :]                                  # fission source -> fission source
    if not want_modes:
        # power iteration started from the flat source; converges to round-off
        s = np.ones(nf)
        k = 1.0
        for it in range(5000):
            s2 = G @ s
            k2 = np.dot(s2, h[:nf]) / np.dot(s, h[:nf])
            s = s2 / k2
            if abs(k2 - k) < 1e-15 * k2 and it > 20:
                k = k2
                break
            k = k2
        return k
    ev, vec = np.linalg.eig(G)
    order = np.argsort(-np.abs(ev))
    ev = ev[order]
    k = ev[0].real
    src = np.abs(vec[:, order[0]].real)
    src /= np.dot(src, h[:nf])
    return dict(k=k, dr=abs(ev[1]) / k, k2=ev[1], src=src, x=np.cumsum(h[:nf]) - h[:nf] / 2)


def derivs(nf, nr, nquad, rel=(2e-3, 1e-3)):
    """Central differences in a, b, rho_f, rho_r with Richardson in the step."""
    out = {}
    k0 = solve(nf=nf, nr=nr, nquad=nquad)
    out["k"] = k0
    for name, base in (("a", A0), ("b", B0), ("rho_f", 1.0), ("rho_r", 1.0)):
        est = []
        for r in rel:
            hst = r * base
            kp = solve(nf=nf, nr=nr, nquad=nquad, **{name: base + hst})
            km = solve(nf=nf, nr=nr, nquad=nquad, **{name: base - hst})
            est.append((kp - km) / (2 * hst))
        # error ~ h^2 -> Richardson with ratio 2
        out["dk_d" + name] = (4 * est[1] - est[0]) / 3
        out["dk_d" + name + "_raw"] = est
    return out


if __name__ == "__main__":
    import json
    import sys
    import time

    if len(sys.argv) > 1 and sys.argv[1] == "modes":
        r = solve(nf=400, nr=400, nquad=32, want_modes=True)
        print(f"k = {r['k']:.8f}  k2 = {r['k2']:.6f}  DR = {r['dr']:.4f}")
        sys.exit()

    results = {}
    t0 = time.time()
    # mesh study at S32, quadrature study on fine meshes (mesh error ~1e-8)
    configs = [(32, 250), (32, 500), (32, 1000), (32, 2000),
               (16, 2000), (64, 2000), (128, 2000), (128, 1000), (256, 1000)]
    for nq, ncell in configs:
        d = derivs(ncell, ncell, nq)
        results[f"S{nq}_n{ncell}"] = d
        print(f"S{nq:3d} cells/region={ncell:5d}  k={d['k']:.10f}  "
              f"dk/da={d['dk_da']:+.8f}  dk/db={d['dk_db']:+.8f}  "
              f"dk/drf={d['dk_drho_f']:+.8f}  dk/drr={d['dk_drho_r']:+.8f}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    # angular error is O(1/N^2) (see log): Richardson on S128/S256
    lo, hi = results["S128_n1000"], results["S256_n1000"]
    ref = {key: (4 * hi[key] - lo[key]) / 3
           for key in ("k", "dk_da", "dk_db", "dk_drho_f", "dk_drho_r")}
    ref["richardson_correction_dk_da"] = ref["dk_da"] - hi["dk_da"]
    for key in ("a", "b", "rho_f", "rho_r"):
        ref["dlnk_d" + key] = ref["dk_d" + key] / ref["k"]
    ref["scaling_lhs"] = A0 * ref["dk_da"] + B0 * ref["dk_db"]
    ref["scaling_rhs"] = ref["dk_drho_f"] + ref["dk_drho_r"]
    results["reference"] = ref
    m = solve(nf=1000, nr=1000, nquad=64, want_modes=True)
    results["modes"] = dict(k=m["k"], dr=m["dr"], k2_real=m["k2"].real, k2_imag=m["k2"].imag)
    np.savez("sn_source.npz", x=m["x"], src=m["src"])
    print("reference:", json.dumps(ref, indent=1))
    print(f"DR = {m['dr']:.5f}")
    with open("sn_results.json", "w") as f:
        json.dump(results, f, indent=1)
