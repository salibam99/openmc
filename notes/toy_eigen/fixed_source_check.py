"""
Fixed-source sanity check of the eigenvalue toy's transport kernel
(mc.transport_one), on the problem of notes/reparam_slab_check_openmc_style.py:
beam at x = 0, mu = 1; A on [0, a), B on [a, L), vacuum; analog, isotropic
scattering, no fission. theta = a, and also theta = L (outer boundary, not in the
reference script) and sigma_t^B at fixed albedo (LR via the density slot).

One generation, active slot only, so variants (a), (b), (c) must agree exactly.
Reference numbers: the reference script (run separately, output pasted in
results.md) and CRN central differences with this kernel.
"""
import numpy as np
import numba as nb

import mc

L_SLAB, ST_A, C_A, ST_B, C_B, A_INT = 2.0, 2.0, 0.5, 0.6, 0.9, 0.9


@nb.njit(parallel=True, cache=True)
def run_chunk(seed, chunk, N, a, Lslab, st, ss, sf, nu, deriv, out):
    # out[i] = (T, R, D_a(a), D_b(a), D_c(a), D_a(L), D_b(L), D_c(L), D_rho_B)
    Da = np.zeros((N, 2, 1))
    ga = np.zeros((N, 2, 1))
    Db = np.zeros((N, 2, 1))
    Dc = np.zeros((N, 2, 1))
    gc = np.zeros((N, 2))
    Dr = np.zeros((N, 2, 1))
    for i in nb.prange(N):
        s = mc.stream(seed, chunk, i)
        res, x, nsite, s = mc.transport_one(i, 0.0, 1.0, 0, 0, s, st, ss, sf, nu, a, Lslab,
                                            Da, ga, Db, Dc, gc, Dr, deriv)
        out[i, 0] = 1.0 if res == mc.LEAK_RIGHT else 0.0
        out[i, 1] = 1.0 if res == mc.LEAK_LEFT else 0.0
        out[i, 2] = Da[i, 0, 0]
        out[i, 3] = Db[i, 0, 0]
        out[i, 4] = Dc[i, 0, 0]
        out[i, 5] = Da[i, 1, 0]
        out[i, 6] = Db[i, 1, 0]
        out[i, 7] = Dc[i, 1, 0]
        out[i, 8] = Dr[i, 1, 0]


def simulate(n_hist, seed, a=A_INT, Lslab=L_SLAB, st_b=ST_B, deriv=True, nchunk=1_000_000):
    st = np.array([ST_A, st_b])
    ss = np.array([C_A * ST_A, C_B * st_b])
    z = np.zeros(2)
    # per-history scores: T, R and their derivative scores
    s1 = np.zeros(14)
    s2 = np.zeros(14)
    out = np.zeros((nchunk, 9))
    for c in range(n_hist // nchunk):
        run_chunk(seed, c, nchunk, a, Lslab, st, ss, z, z, deriv, out)
        T, R = out[:, 0], out[:, 1]
        sc = np.column_stack([T, R,
                              T * out[:, 2], T * out[:, 3], T * out[:, 4],
                              R * out[:, 2],
                              T * out[:, 5], T * out[:, 6], T * out[:, 7],
                              R * out[:, 5],
                              T * out[:, 8] / st_b, R * out[:, 8] / st_b,
                              T * out[:, 2] - T * out[:, 4], T * out[:, 5] - T * out[:, 7]])
        s1 += sc.sum(0)
        s2 += (sc ** 2).sum(0)
    n = (n_hist // nchunk) * nchunk
    mean = s1 / n
    se = np.sqrt((s2 / n - mean ** 2) / n)
    return mean, se


def fd(param, n_hist, seed, hs=np.linspace(-0.06, 0.06, 7)):
    """CRN central difference via quadratic fit, as in the reference script."""
    base = dict(a=A_INT, Lslab=L_SLAB, st_b=ST_B)
    T, R = [], []
    for h in hs:
        kw = dict(base)
        kw[param] += h
        m, _ = simulate(n_hist, seed, deriv=False, **kw)
        T.append(m[0])
        R.append(m[1])
    return np.polyfit(hs, T, 2)[1], np.polyfit(hs, R, 2)[1]


if __name__ == "__main__":
    import json
    N = 40_000_000
    m, se = simulate(N, 12345)
    names = ["T", "R", "dT/da (a)", "dT/da (b)", "dT/da (c)", "dR/da",
             "dT/dL (a)", "dT/dL (b)", "dT/dL (c)", "dR/dL",
             "dT/dsig_B", "dR/dsig_B", "(a)-(c) dT/da", "(a)-(c) dT/dL"]
    res = {}
    for nm, v, e in zip(names, m, se):
        print(f"{nm:16s} {v:+.5f} +- {e:.5f}")
        res[nm] = (v, e)
    # independent FD seeds; replicate the FD 4x to get its own error bar
    for p in ("a", "Lslab", "st_b"):
        reps = np.array([fd(p, 10_000_000, 1000 + r) for r in range(4)])
        mt, mr = reps.mean(0)
        et, er = reps.std(0, ddof=1) / 2
        print(f"FD d/d{p}:  dT = {mt:+.5f} +- {et:.5f}   dR = {mr:+.5f} +- {er:.5f}")
        res[f"FD dT/d{p}"] = (mt, et)
        res[f"FD dR/d{p}"] = (mr, er)
    with open("fixed_source_results.json", "w") as f:
        json.dump(res, f, indent=1)
