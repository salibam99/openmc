"""
Side probe (not part of the 1D toy): does a PASSIVE g (V_n = 0 on every surface)
decay or grow across generations in 3D? In 1D every crossing zeroes g. In 3D a
piece maps g -> (I - u omega n^T / (n.omega)) g, with eigenvalues {1, 1, 1-u}.
Tangential components are not damped and grazing hits amplify them, so with
several non-parallel surfaces the product of maps could grow. That would decide
whether a never-reset g (variant c) or long passive windows (variant a) are
usable in 3D.

Same physics as the 1D toy (same cross sections), analog, combing. In generation
0 every source particle gets a random unit g. Then g is propagated passively and
inherited at fission (variant-c style, never reset). We report RMS |g| over
fission sites per generation. Geometries: concentric spheres (R1 = 3, R2 = 4.5)
and nested cubes (half-widths 3 and 4.5).
"""
import math
import sys

import numba as nb
import numpy as np

import mc
from model import FUEL, REFL

R1, R2 = 3.0, 4.5


@nb.njit(inline="always")
def sphere_hit(px, py, pz, ox, oy, oz, R, outward):
    # distance along omega to |x| = R (smallest positive root); inf if none
    b = px * ox + py * oy + pz * oz
    c = px * px + py * py + pz * pz - R * R
    disc = b * b - c
    if disc < 0.0:
        return np.inf
    sq = math.sqrt(disc)
    t = -b + sq if outward else -b - sq
    return t if t > 1e-12 else np.inf


@nb.njit(inline="always")
def box_exit(p, o, A):
    # inside box |x|_inf < A: distance to exit and axis of the exit face
    tbest = np.inf
    ax = 0
    for k in range(3):
        if o[k] > 0.0:
            t = (A - p[k]) / o[k]
        elif o[k] < 0.0:
            t = (-A - p[k]) / o[k]
        else:
            t = np.inf
        if t < tbest:
            tbest = t
            ax = k
    return tbest, ax


@nb.njit(inline="always")
def box_enter(p, o, A):
    # outside box: distance to entry (slab method) and entry axis; inf if missed
    tmin = -np.inf
    tmax = np.inf
    ax = 0
    for k in range(3):
        if o[k] != 0.0:
            t1 = (-A - p[k]) / o[k]
            t2 = (A - p[k]) / o[k]
            lo = min(t1, t2)
            hi = max(t1, t2)
            if lo > tmin:
                tmin = lo
                ax = k
            if hi < tmax:
                tmax = hi
        elif abs(p[k]) >= A:
            return np.inf, 0
    if tmax < tmin or tmin <= 1e-12:
        return np.inf, 0
    return tmin, ax


@nb.njit(parallel=True, cache=True)
def generation(n, seed, geom, sx, sg, fx, fg, fnu, st, ss, sf, nu):
    N = sx.shape[0]
    for i in nb.prange(N):
        s = mc.stream(seed, n, i)
        p = sx[i].copy()
        g = sg[i].copy()
        s, u1 = mc.rnd(s)
        s, u2 = mc.rnd(s)
        o = np.empty(3)
        cz = 2 * u1 - 1
        sz = math.sqrt(1 - cz * cz)
        o[0] = sz * math.cos(2 * math.pi * u2)
        o[1] = sz * math.sin(2 * math.pi * u2)
        o[2] = cz
        if geom == 0:
            r = 0 if p[0] ** 2 + p[1] ** 2 + p[2] ** 2 < R1 * R1 else 1
        else:
            r = 0 if max(abs(p[0]), abs(p[1]), abs(p[2])) < R1 else 1
        fnu[i] = 0
        while True:
            nrm = np.zeros(3)
            leak_if_cross = False
            if geom == 0:
                if r == 0:
                    l = sphere_hit(p[0], p[1], p[2], o[0], o[1], o[2], R1, True)
                    rnext = 1
                else:
                    li = sphere_hit(p[0], p[1], p[2], o[0], o[1], o[2], R1, False)
                    lo = sphere_hit(p[0], p[1], p[2], o[0], o[1], o[2], R2, True)
                    if li < lo:
                        l = li
                        rnext = 0
                    else:
                        l = lo
                        rnext = -1
                        leak_if_cross = True
                q = p + l * o
                nrm[:] = q / math.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2)
            else:
                if r == 0:
                    l, ax = box_exit(p, o, R1)
                    rnext = 1
                else:
                    li, axi = box_enter(p, o, R1)
                    lo, axo = box_exit(p, o, R2)
                    if li < lo:
                        l, ax = li, axi
                        rnext = 0
                    else:
                        l, ax = lo, axo
                        rnext = -1
                        leak_if_cross = True
                q = p + l * o
                nrm[ax] = 1.0 if q[ax] > 0 else -1.0
            sig = st[r]
            s, xi = mc.rnd(s)
            d = -math.log(xi) / sig
            ndo = nrm[0] * o[0] + nrm[1] * o[1] + nrm[2] * o[2]
            ndg = nrm[0] * g[0] + nrm[1] * g[1] + nrm[2] * g[2]
            dl = -ndg / ndo                       # passive: V_n = 0
            if d < l:
                g += o * (d / l) * dl
                p = p + d * o
                s, xi = mc.rnd(s)
                if xi * sig < ss[r]:
                    s, u1 = mc.rnd(s)
                    s, u2 = mc.rnd(s)
                    cz = 2 * u1 - 1
                    sz = math.sqrt(1 - cz * cz)
                    o[0] = sz * math.cos(2 * math.pi * u2)
                    o[1] = sz * math.sin(2 * math.pi * u2)
                    o[2] = cz
                elif xi * sig < ss[r] + sf[r]:
                    s, xi = mc.rnd(s)
                    k = int(nu[r])
                    fnu[i] = k + (1 if xi < nu[r] - k else 0)
                    fx[i] = p
                    fg[i] = g
                    break
                else:
                    break
            else:
                g += o * dl
                p = q
                if leak_if_cross:
                    break
                r = rnext


@nb.njit(cache=True)
def comb(n, seed, fx, fg, fnu, sx, sg):
    N = sx.shape[0]
    M = 0
    for i in range(N):
        M += fnu[i]
    s = mc.stream(seed, n, 0xFFFFFFFF)
    s, xi = mc.rnd(s)
    i = 0
    cum = fnu[0]
    for k in range(N):
        pos = (k + xi) * M / N
        while cum <= pos:
            i += 1
            cum += fnu[i]
        sx[k] = fx[i]
        sg[k] = fg[i]
    return M


def probe(geom, N=200_000, n_in=60, n_g=25, seed=5):
    st = np.array([FUEL["st"], REFL["st"]])
    ss = np.array([FUEL["ss"], REFL["ss"]])
    sf = np.array([FUEL["sf"], REFL["sf"]])
    nu = np.array([FUEL["nu"], REFL["nu"]])
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(N, 3))
    sx = R1 * 0.5 * v / np.linalg.norm(v, axis=1)[:, None] * rng.random((N, 1)) ** (1 / 3)
    sg = np.zeros((N, 3))
    fx, fg, fnu = np.zeros((N, 3)), np.zeros((N, 3)), np.zeros(N, dtype=np.int64)
    rms, p99, mx, ks = [], [], [], []
    for n in range(n_in + n_g):
        if n == n_in:      # start: random unit g on every source particle
            w = rng.normal(size=(N, 3))
            sg[:] = w / np.linalg.norm(w, axis=1)[:, None]
        generation(n, seed, geom, sx, sg, fx, fg, fnu, st, ss, sf, nu)
        if n >= n_in:
            gm = np.linalg.norm(fg[fnu > 0], axis=1)
            rms.append(np.sqrt(np.mean(gm ** 2)))
            p99.append(np.quantile(gm, 0.99))
            mx.append(gm.max())
        M = comb(n, seed, fx, fg, fnu, sx, sg)
        ks.append(M / N)
    return np.array(rms), np.array(p99), np.array(mx), np.mean(ks[n_in:])


if __name__ == "__main__":
    for geom, name in ((0, "concentric spheres"), (1, "nested cubes")):
        rms, p99, mx, k = probe(geom)
        print(f"{name}: k = {k:.4f}")
        print("  gen  RMS|g|     99%|g|     max|g|")
        for j in range(len(rms)):
            print(f"  {j:3d}  {rms[j]:.3e}  {p99[j]:.3e}  {mx[j]:.3e}")
        # per-generation decay factor of RMS from a log-linear fit over gens 2..
        j = np.arange(2, len(rms))
        sl = np.polyfit(j, np.log(rms[2:]), 1)[0]
        print(f"  RMS decay factor per generation (fit gens 2-{len(rms)-1}): {np.exp(sl):.3f}")
        sys.stdout.flush()
