"""
Driver: independent-seed eigenvalue runs (sensitivities, all variants, L = 0..10)
and CRN central-difference runs for k.

    python run_eigen.py sens  R N n_inactive n_active   -> runs/sens_seedXXX.npz
    python run_eigen.py fd    R N n_inactive n_active h -> runs/fd_seedXXX.npz
"""
import sys
import time
from pathlib import Path

import numpy as np

import mc
from model import A0, B0, FUEL, REFL

OUT = Path(__file__).parent / "runs"


def sens(R, N, n_in, n_act, first_seed=1):
    OUT.mkdir(exist_ok=True)
    for seed in range(first_seed, first_seed + R):
        f = OUT / f"sens_seed{seed:03d}.npz"
        if f.exists():
            continue
        t = time.time()
        r = mc.run(seed, N, n_in, n_act, A0, B0, FUEL, REFL, deriv=True)
        dt = time.time() - t
        np.savez(f, kgen=r["kgen"], S1=r["S1"], S2=r["S2"], den=r["den"],
                 N=N, n_inactive=n_in, n_active=n_act, seconds=dt)
        est = r["S1"].sum(0) / r["den"].sum()
        print(f"seed {seed:3d}  k = {r['kgen'].mean():.5f}  "
              f"(a) dlnk/da[L=8] = {est[0, 8]:+.5f}  (c) = {est[4, 8]:+.5f}  "
              f"dlnk/db[L=8] = {est[1, 8]:+.5f}  [{dt:.0f}s]", flush=True)


def fd(R, N, n_in, n_act, h, first_seed=1001):
    OUT.mkdir(exist_ok=True)
    for seed in range(first_seed, first_seed + R):
        f = OUT / f"fd_seed{seed:04d}.npz"
        if f.exists():
            continue
        t = time.time()
        k = {}
        for name, (a, b) in {"a+": (A0 + h, B0), "a-": (A0 - h, B0),
                             "b+": (A0, B0 + h), "b-": (A0, B0 - h)}.items():
            k[name] = mc.run(seed, N, n_in, n_act, a, b, FUEL, REFL, deriv=False)["kgen"].mean()
        np.savez(f, h=h, N=N, n_inactive=n_in, n_active=n_act, **{kk.replace("+", "p").replace("-", "m"): v for kk, v in k.items()})
        print(f"seed {seed}  dk/da = {(k['a+'] - k['a-']) / (2 * h):+.5f}  "
              f"dk/db = {(k['b+'] - k['b-']) / (2 * h):+.5f}  [{time.time() - t:.0f}s]", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    R, N, n_in, n_act = (int(v) for v in sys.argv[2:6])
    if mode == "sens":
        sens(R, N, n_in, n_act)
    else:
        fd(R, N, n_in, n_act, float(sys.argv[6]))
