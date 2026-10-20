"""
3D Monte Carlo eigenvalue toy: one group, isotropic scattering, analog, vacuum
outside. Distance reparameterization + IFP, variants (a) per-origin slots and
(c) single inherited g (../derivation.md), carried on the same random walk.

Surfaces (id): 0 fuel sphere (centre (x0, y0, z0), radius a), 1 outer sphere
(centre 0, radius b), 2 cylinder x^2 + y^2 = c^2, 3 plane z = zlo, 4 plane
z = zhi. Each F is oriented with n = grad F/|grad F| (outward for spheres and
cylinder, +z for planes). Cells: 0 fuel, 1 reflector, -1 outside.

Parameter kinds and V_n = -(dF/dtheta)/|grad F| on the parameter's surface:
    radius / plane offset (kind 0): V_n = 1
    sphere-0 translation along e (kind 1): V_n = n.e
Per-piece rules (unchanged from 1D, now with 3-vectors):
    dl = (V_n - n.g) / (n.omega)
    crossing / leak:         D += -st dl,              g += omega dl
    collision at u = d/l:    D += -st u dl + dl/l,     g += omega u dl
(a): slot n active (V_n), slots o < n passive (V_n = 0); (c): one g, increments
into slot n. Densities by LR: D -= st * distance, +1 per collision.
"""
import math
import sys
from pathlib import Path

import numba as nb
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mc import stream, rnd  # noqa: E402

from geom3d import CASES, FUEL, REFL  # noqa: E402

LMAX = 12
W = LMAX + 1
NB = 180                     # log10 histogram bins on [HLO, HHI] (0.1 decade)
HLO, HHI = -10.0, 8.0
KTOP = 6                     # top fission histories kept per generation, variant, parameter
NTF = 9                      # fields of a top record, see score_generation
INF = 1e300


def setup(case, **pert):
    """Geometry vector gp and parameter tables for a case (optionally perturbed)."""
    c = dict(CASES[case])
    c.update(pert)
    # gp = [a, x0, y0, z0, b, c, zlo, zhi]
    if case == 1:
        gp = [c["R"], c.get("tx", 0.0), c.get("ty", 0.0), c.get("tz", 0.0), 0, 0, 0, 0]
        names = ["R", "tx"]
        kind, surf = [0, 1], [0, 0]
        pdir = [[0, 0, 0], [1, 0, 0]]
    elif case == 2:
        gp = [c["a"], 0, 0, 0, c["b"], 0, 0, 0]
        names = ["a", "b"]
        kind, surf = [0, 0], [0, 1]
        pdir = [[0, 0, 0], [0, 0, 0]]
    else:
        gp = [c["a"], c["x0"] + c.get("tx", 0.0), c.get("ty", 0.0), c.get("tz", 0.0),
              0, c["c"], c["zlo"], c["zhi"]]
        names = ["a", "tx", "ty", "tz", "c", "zlo", "zhi"]
        kind, surf = [0, 1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 2, 3, 4]
        pdir = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    rho = [c.get("rho_f", 1.0), c.get("rho_r", 1.0)]
    st = np.array([rho[0] * FUEL["st"], rho[1] * REFL["st"]])
    ss = np.array([rho[0] * FUEL["ss"], rho[1] * REFL["ss"]])
    sf = np.array([rho[0] * FUEL["sf"], rho[1] * REFL["sf"]])
    nu = np.array([FUEL["nu"], REFL["nu"]])
    return dict(gp=np.array(gp, float), names=names, kind=np.array(kind, np.int64),
                surf=np.array(surf, np.int64), pdir=np.array(pdir, float),
                st=st, ss=ss, sf=sf, nu=nu)


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------
@nb.njit(inline="always")
def _sph(px, py, pz, ox, oy, oz, cx, cy, cz, R, inside):
    dx, dy, dz = px - cx, py - cy, pz - cz
    b = dx * ox + dy * oy + dz * oz
    cc = dx * dx + dy * dy + dz * dz - R * R
    disc = b * b - cc
    if inside:
        return -b + math.sqrt(max(disc, 0.0))
    if disc <= 0.0 or b >= 0.0:
        return INF
    t = -b - math.sqrt(disc)
    return t if t > 1e-10 else INF


@nb.njit(inline="always")
def next_surface(case, gp, px, py, pz, ox, oy, oz, cell):
    """Distance to the bounding surface, its id and the cell beyond (-1: leak)."""
    if cell == 0:
        l = _sph(px, py, pz, ox, oy, oz, gp[1], gp[2], gp[3], gp[0], True)
        return l, 0, (-1 if case == 1 else 1)
    # reflector
    l = _sph(px, py, pz, ox, oy, oz, gp[1], gp[2], gp[3], gp[0], False)
    sid = 0
    nxt = 0
    if case == 2:
        t = _sph(px, py, pz, ox, oy, oz, 0.0, 0.0, 0.0, gp[4], True)
        if t < l:
            l, sid, nxt = t, 1, -1
        return l, sid, nxt
    a2 = ox * ox + oy * oy
    if a2 > 0.0:
        b2 = px * ox + py * oy
        c2 = px * px + py * py - gp[5] * gp[5]
        t = (-b2 + math.sqrt(max(b2 * b2 - a2 * c2, 0.0))) / a2
        if t < l:
            l, sid, nxt = t, 2, -1
    if oz > 0.0:
        t = (gp[7] - pz) / oz
        if t < l:
            l, sid, nxt = t, 4, -1
    elif oz < 0.0:
        t = (gp[6] - pz) / oz
        if t < l:
            l, sid, nxt = t, 3, -1
    return l, sid, nxt


@nb.njit(inline="always")
def normal(sid, qx, qy, qz, gp):
    if sid == 0:
        return (qx - gp[1]) / gp[0], (qy - gp[2]) / gp[0], (qz - gp[3]) / gp[0]
    if sid == 1:
        return qx / gp[4], qy / gp[4], qz / gp[4]
    if sid == 2:
        rr = math.sqrt(qx * qx + qy * qy)
        return qx / rr, qy / rr, 0.0
    return 0.0, 0.0, 1.0


@nb.njit(inline="always")
def surf_dist(sid, px, py, pz, gp):
    """Exact distance from (px, py, pz) to surface sid."""
    if sid == 0:
        return abs(math.sqrt((px - gp[1]) ** 2 + (py - gp[2]) ** 2 + (pz - gp[3]) ** 2) - gp[0])
    if sid == 1:
        return abs(math.sqrt(px * px + py * py + pz * pz) - gp[4])
    if sid == 2:
        return abs(math.sqrt(px * px + py * py) - gp[5])
    if sid == 3:
        return abs(pz - gp[6])
    return abs(pz - gp[7])


@nb.njit(inline="always")
def hbin(v):
    if v <= 0.0:
        return 0
    b = int((math.log10(v) - HLO) / (HHI - HLO) * NB)
    return min(max(b, 0), NB - 1)


@nb.njit(inline="always")
def iso(s):
    s, u1 = rnd(s)
    s, u2 = rnd(s)
    cz = 2.0 * u1 - 1.0
    sz = math.sqrt(max(1.0 - cz * cz, 0.0))
    ph = 2.0 * math.pi * u2
    return s, sz * math.cos(ph), sz * math.sin(ph), cz


# --------------------------------------------------------------------------
# one generation
# --------------------------------------------------------------------------
@nb.njit(parallel=True, cache=True)
def transport_generation(n, seed, case, gp, kind, surf, pdir, st, ss, sf, nu,
                         sx, sDa, sga, sDc, sgc, sDr,
                         fx, fnu, Da, ga, Dc, gc, Dr, Jd, Jx,
                         hl, hs, hj, deriv, diag):
    N = sx.shape[0]
    P = kind.size
    Wl = Da.shape[2]
    act = n % Wl
    for i in nb.prange(N):
        tid = nb.get_thread_id()
        s = stream(seed, n, i)
        Da[i] = sDa[i]
        ga[i] = sga[i]
        Dc[i] = sDc[i]
        gc[i] = sgc[i]
        Dr[i] = sDr[i]
        Jd[i] = 0.0
        Jx[i] = 0.0
        px, py, pz = sx[i, 0], sx[i, 1], sx[i, 2]
        s, ox, oy, oz = iso(s)
        dx_, dy_, dz_ = px - gp[1], py - gp[2], pz - gp[3]
        cell = 0 if dx_ * dx_ + dy_ * dy_ + dz_ * dz_ < gp[0] * gp[0] else 1
        fnu[i] = 0
        while True:
            l, sid, nxt = next_surface(case, gp, px, py, pz, ox, oy, oz, cell)
            sig = st[cell]
            s, xi = rnd(s)
            d = -math.log(xi) / sig
            collided = d < l
            qx, qy, qz = px + l * ox, py + l * oy, pz + l * oz
            if deriv:
                nx, ny, nz = normal(sid, qx, qy, qz, gp)
                ndo = nx * ox + ny * oy + nz * oz
                u = d / l if collided else 1.0
                if collided and diag:
                    hl[tid, hbin(l)] += 1
                    hs[tid, hbin(l * abs(ndo))] += 1
                for p in range(P):
                    V = 0.0
                    if sid == surf[p]:
                        if kind[p] == 0:
                            V = 1.0
                        else:
                            V = nx * pdir[p, 0] + ny * pdir[p, 1] + nz * pdir[p, 2]
                    # (a) per-origin slots
                    for sl in range(Wl):
                        Vs = V if sl == act else 0.0
                        num = Vs - (nx * ga[i, p, sl, 0] + ny * ga[i, p, sl, 1] + nz * ga[i, p, sl, 2])
                        if num == 0.0:
                            continue
                        dl = num / ndo
                        if collided:
                            Da[i, p, sl] += -sig * u * dl + dl / l
                            ga[i, p, sl, 0] += ox * u * dl
                            ga[i, p, sl, 1] += oy * u * dl
                            ga[i, p, sl, 2] += oz * u * dl
                            if sl == act:
                                Jd[i, 0, p] += dl / l
                                if abs(dl / l) > abs(Jx[i, 0, p, 0]):
                                    Jx[i, 0, p, 0] = dl / l
                                    Jx[i, 0, p, 1] = l
                                    Jx[i, 0, p, 2] = l * abs(ndo)
                                    Jx[i, 0, p, 3] = surf_dist(sid, px, py, pz, gp)
                                    Jx[i, 0, p, 4] = sid
                                if diag:
                                    hj[tid, 0, p, hbin(abs(dl / l))] += 1
                        else:
                            Da[i, p, sl] += -sig * dl
                            ga[i, p, sl, 0] += ox * dl
                            ga[i, p, sl, 1] += oy * dl
                            ga[i, p, sl, 2] += oz * dl
                    # (c) single inherited g
                    num = V - (nx * gc[i, p, 0] + ny * gc[i, p, 1] + nz * gc[i, p, 2])
                    if num != 0.0:
                        dl = num / ndo
                        if collided:
                            Dc[i, p, act] += -sig * u * dl + dl / l
                            gc[i, p, 0] += ox * u * dl
                            gc[i, p, 1] += oy * u * dl
                            gc[i, p, 2] += oz * u * dl
                            Jd[i, 1, p] += dl / l
                            if abs(dl / l) > abs(Jx[i, 1, p, 0]):
                                Jx[i, 1, p, 0] = dl / l
                                Jx[i, 1, p, 1] = l
                                Jx[i, 1, p, 2] = l * abs(ndo)
                                Jx[i, 1, p, 3] = surf_dist(sid, px, py, pz, gp)
                                Jx[i, 1, p, 4] = sid
                            if diag:
                                hj[tid, 1, p, hbin(abs(dl / l))] += 1
                        else:
                            Dc[i, p, act] += -sig * dl
                            gc[i, p, 0] += ox * dl
                            gc[i, p, 1] += oy * dl
                            gc[i, p, 2] += oz * dl
                # densities (LR)
                if collided:
                    Dr[i, cell, act] += 1.0 - sig * d
                else:
                    Dr[i, cell, act] -= sig * l
            if collided:
                px, py, pz = px + d * ox, py + d * oy, pz + d * oz
                s, xi = rnd(s)
                if xi * sig < ss[cell]:
                    s, ox, oy, oz = iso(s)
                elif xi * sig < ss[cell] + sf[cell]:
                    s, xi = rnd(s)
                    k = int(nu[cell])
                    fnu[i] = k + (1 if xi < nu[cell] - k else 0)
                    fx[i, 0], fx[i, 1], fx[i, 2] = px, py, pz
                    break
                else:
                    break
            else:
                px, py, pz = qx, qy, qz
                if nxt < 0:
                    break
                cell = nxt


@nb.njit(parallel=True, cache=True)
def score_generation(n, fnu, Da, Dc, Dr, Jd, Jx, ga, gc, S1, S2, Smax, jstat, gstat, ghist, top):
    """IFP scores of generation n (w = 1 per fission event) and diagnostics.

    S1/S2/Smax[ch, L]: sum, sum of squares and max |.| of D_{n-L} over fission
    events. Channels: (a) params, (c) params, densities.
    jstat[v, p]: [count, n(|J| > |D_own - J|), sum J^2, sum D_own^2, sum (D_own-J)^2],
    where D_own is the own-generation slot and J its dl/l part.
    gstat[v, p, age]: [sum |g|^2, max |g|, count] at fission sites; (c) age 0 only.
    top[v, p, k]: the KTOP fission histories with largest |D_own| in this generation:
    [D_own, J (all dl/l terms), largest single dl/l term, its piece length l,
    normal distance s = l |n.omega|, exact distance from the piece start to that
    surface, surface id, generation, particle index].
    """
    N = fnu.size
    P = Da.shape[1]
    Wl = Da.shape[2]
    nL = S1.shape[1]
    nch = S1.shape[0]
    act = n % Wl
    for job in nb.prange(nch * nL):
        ch = job // nL
        L = job % nL
        sl = (n - L) % Wl
        acc = 0.0
        acc2 = 0.0
        mx = 0.0
        for i in range(N):
            if fnu[i] == 0:
                continue
            if ch < P:
                v = Da[i, ch, sl]
            elif ch < 2 * P:
                v = Dc[i, ch - P, sl]
            else:
                v = Dr[i, ch - 2 * P, sl]
            acc += v
            acc2 += v * v
            if abs(v) > mx:
                mx = abs(v)
        S1[ch, L] += acc
        S2[ch, L] += acc2
        Smax[ch, L] = max(Smax[ch, L], mx)
    for job in nb.prange(2 * P):
        v = job // P
        p = job % P
        for i in range(N):
            if fnu[i] == 0:
                continue
            Do = Da[i, p, act] if v == 0 else Dc[i, p, act]
            J = Jd[i, v, p]
            jstat[v, p, 0] += 1
            if abs(J) > abs(Do - J):
                jstat[v, p, 1] += 1
            jstat[v, p, 2] += J * J
            jstat[v, p, 3] += Do * Do
            jstat[v, p, 4] += (Do - J) ** 2
            # top-K by |D_own| (insertion into a short sorted list)
            K = top.shape[2]
            if abs(Do) > abs(top[v, p, K - 1, 0]):
                k = K - 1
                while k > 0 and abs(Do) > abs(top[v, p, k - 1, 0]):
                    top[v, p, k] = top[v, p, k - 1]
                    k -= 1
                top[v, p, k, 0] = Do
                top[v, p, k, 1] = J
                for f in range(5):
                    top[v, p, k, 2 + f] = Jx[i, v, p, f]
                top[v, p, k, 7] = n
                top[v, p, k, 8] = i
            nage = Wl if v == 0 else 1
            for age in range(nage):
                if v == 0:
                    sl = (n - age) % Wl
                    g2 = ga[i, p, sl, 0] ** 2 + ga[i, p, sl, 1] ** 2 + ga[i, p, sl, 2] ** 2
                else:
                    g2 = gc[i, p, 0] ** 2 + gc[i, p, 1] ** 2 + gc[i, p, 2] ** 2
                gm = math.sqrt(g2)
                gstat[v, p, age, 0] += g2
                gstat[v, p, age, 1] = max(gstat[v, p, age, 1], gm)
                gstat[v, p, age, 2] += 1
                if gm > 0.0:
                    ghist[v, p, age, hbin(gm)] += 1
    nf = 0
    for i in range(N):
        if fnu[i] > 0:
            nf += 1
    return nf


@nb.njit(cache=True)
def comb(n, seed, fx, fnu, Da, ga, Dc, gc, Dr, sx, sDa, sga, sDc, sgc, sDr):
    """Comb to N equal-weight sites; daughters copy slots, open slot (n+1) % W."""
    N = sx.shape[0]
    M = 0
    for i in range(N):
        M += fnu[i]
    s = stream(seed, n, 0xFFFFFFFF)
    s, xi = rnd(s)
    nxt = (n + 1) % Da.shape[2]
    i = 0
    cum = fnu[0]
    for k in range(N):
        pos = (k + xi) * M / N
        while cum <= pos:
            i += 1
            cum += fnu[i]
        sx[k] = fx[i]
        sDa[k] = Da[i]
        sga[k] = ga[i]
        sDc[k] = Dc[i]
        sgc[k] = gc[i]
        sDr[k] = Dr[i]
        sDa[k, :, nxt] = 0.0
        sga[k, :, nxt] = 0.0
        sDc[k, :, nxt] = 0.0
        sDr[k, :, nxt] = 0.0
    return M


def run(case, seed, N, n_in, n_act, deriv=True, diag=True, **pert):
    """Power iteration. Returns per-generation k and, if deriv, IFP sums/diagnostics."""
    c = setup(case, **pert)
    gp, P = c["gp"], len(c["names"])
    Wd = W if deriv else 1
    Pd = P if deriv else 1
    nT = nb.config.NUMBA_NUM_THREADS
    src = [np.zeros((N, 3)), np.zeros((N, Pd, Wd)), np.zeros((N, Pd, Wd, 3)),
           np.zeros((N, Pd, Wd)), np.zeros((N, Pd, 3)), np.zeros((N, 2, Wd))]
    rec = [np.zeros((N, 3)), np.zeros(N, dtype=np.int64), np.zeros((N, Pd, Wd)),
           np.zeros((N, Pd, Wd, 3)), np.zeros((N, Pd, Wd)), np.zeros((N, Pd, 3)),
           np.zeros((N, 2, Wd)), np.zeros((N, 2, Pd)), np.zeros((N, 2, Pd, 5))]
    hl, hs = np.zeros((nT, NB)), np.zeros((nT, NB))
    hj = np.zeros((nT, 2, Pd, NB))
    # initial source: uniform in the (possibly perturbed) fuel sphere; CRN-consistent
    rng = np.random.default_rng(seed + 7919)
    v = rng.normal(size=(N, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    src[0][:] = gp[1:4] + gp[0] * v * rng.random((N, 1)) ** (1 / 3)
    ngen = n_in + n_act
    nch = 2 * P + 2
    kgen = np.zeros(ngen)
    S1 = np.zeros((n_act, nch, W))
    S2 = np.zeros((n_act, nch, W))
    Smax = np.zeros((n_act, nch, W))
    den = np.zeros(n_act)
    jstat = np.zeros((2, P, 5))
    gstat = np.zeros((ngen, 2, P, W, 3))
    ghist = np.zeros((2, P, W, NB))
    ghist_g = np.zeros((2, P, W, NB))
    top = np.zeros((n_act, 2, P, KTOP, NTF))
    args = (case, gp, c["kind"][:Pd], c["surf"][:Pd], c["pdir"][:Pd], c["st"], c["ss"], c["sf"], c["nu"])
    for n in range(ngen):
        transport_generation(n, seed, *args, *src, *rec, hl, hs, hj, deriv,
                             diag and deriv and n >= n_in)
        if deriv:
            g = n - n_in
            if g >= 0:
                den[g] = score_generation(n, rec[1], rec[2], rec[4], rec[6], rec[7], rec[8], rec[3], rec[5],
                                          S1[g], S2[g], Smax[g], jstat, gstat[n], ghist, top[g])
            else:   # inactive: g statistics only (evolution from g = 0)
                score_generation(n, rec[1], rec[2], rec[4], rec[6], rec[7], rec[8], rec[3], rec[5],
                                 np.zeros((nch, W)), np.zeros((nch, W)), np.zeros((nch, W)),
                                 np.zeros((2, P, 5)), gstat[n], ghist_g, np.zeros((2, P, KTOP, NTF)))
        M = comb(n, seed, rec[0], rec[1], rec[2], rec[3], rec[4], rec[5], rec[6], *src)
        kgen[n] = M / N
    out = dict(kgen=kgen, n_inactive=n_in, names=np.array(c["names"]))
    if deriv:
        out.update(S1=S1, S2=S2, Smax=Smax, den=den, jstat=jstat, gstat=gstat, ghist=ghist,
                   hl=hl.sum(0), hs=hs.sum(0), hj=hj.sum(0), top=top.astype(np.float32))
    return out
