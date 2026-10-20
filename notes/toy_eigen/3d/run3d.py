"""
Driver for the 3D toy.

    python run3d.py sens CASE R N n_inactive n_active [first_seed]
        -> runs/caseC/sens_seedXXX.npz   (all params, variants (a)/(c), L = 0..8)
    python run3d.py fd CASE R N n_inactive n_active h PARAM[,PARAM...] [first_seed]
        -> runs/caseC/fd_PARAM_seedXXXX.npz   (CRN pair k(+h), k(-h), k-only runs)
Existing files are skipped, so runs can be extended.
"""
import sys
import time
from pathlib import Path

import numpy as np

import mc3d
from geom3d import CASES

OUT = Path(__file__).parent / "runs"


def pert_for(case, param, h):
    """Keyword perturbation for mc3d.run (translations are offsets from the base)."""
    if param in ("tx", "ty", "tz"):
        return {param: h}
    if param in ("rho_f", "rho_r"):
        return {param: 1.0 + h}
    return {param: CASES[case][param] + h}


def sens(case, R, N, n_in, n_act, first=1):
    d = OUT / f"case{case}"
    d.mkdir(parents=True, exist_ok=True)
    for seed in range(first, first + R):
        f = d / f"sens_seed{seed:03d}.npz"
        if f.exists():
            continue
        t = time.time()
        r = mc3d.run(case, seed, N, n_in, n_act)
        dt = time.time() - t
        np.savez(f, N=N, seconds=dt, **r)
        est = r["S1"].sum(0) / r["den"].sum()
        P = len(r["names"])
        print(f"case {case} seed {seed:3d}  k = {r['kgen'][n_in:].mean():.5f}  "
              + "  ".join(f"{r['names'][p]}: a {est[p, 8]:+.4f} c {est[P + p, 8]:+.4f}" for p in range(P))
              + f"  [{dt:.0f}s]", flush=True)


def fd(case, R, N, n_in, n_act, h, params, first=1001):
    d = OUT / f"case{case}"
    d.mkdir(parents=True, exist_ok=True)
    for seed in range(first, first + R):
        for p in params:
            f = d / f"fd_{p}_seed{seed:04d}.npz"
            if f.exists():
                continue
            t = time.time()
            kp = mc3d.run(case, seed, N, n_in, n_act, deriv=False, **pert_for(case, p, +h))["kgen"][n_in:].mean()
            km = mc3d.run(case, seed, N, n_in, n_act, deriv=False, **pert_for(case, p, -h))["kgen"][n_in:].mean()
            dt = time.time() - t
            np.savez(f, kp=kp, km=km, h=h, N=N, n_inactive=n_in, n_active=n_act, seconds=dt)
            print(f"case {case} seed {seed} {p}: dk/d{p} = {(kp - km) / (2 * h):+.5f}  [{dt:.0f}s]", flush=True)


if __name__ == "__main__":
    mode, case = sys.argv[1], int(sys.argv[2])
    R, N, n_in, n_act = (int(v) for v in sys.argv[3:7])
    if mode == "sens":
        sens(case, R, N, n_in, n_act, *(int(v) for v in sys.argv[7:8]))
    else:
        fd(case, R, N, n_in, n_act, float(sys.argv[7]), sys.argv[8].split(","),
           *(int(v) for v in sys.argv[9:10]))
