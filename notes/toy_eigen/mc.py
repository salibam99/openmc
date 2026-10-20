"""
Monte Carlo eigenvalue toy: 1D slab, one group, isotropic scattering, analog,
fuel on [0, a], reflector on [a, b], vacuum on both sides.

Geometric parameters theta1 = a (interface), theta2 = b (outer boundary), with
distance reparameterization. Density parameters rho_f, rho_r by likelihood ratio.
k sensitivities by IFP with per-generation slots, see derivation.md.

All three variants (a) per-origin slots, (b) g reset at birth, (c) single
inherited g are carried on the same random walk, as is every window
L = 0..LMAX. They differ only in derivative bookkeeping, never in transport.

Per-piece rules (1D: n = +x, n.omega = mu, n.g = g):
    dl = (V - g) / mu
    crossing / leak:            D += -st dl,              g += mu dl
    collision at u = d/l:       D += -st u dl + dl / l,   g += mu u dl
Density (LR, rho = 1):          D += -st * (distance), +1 per collision.
"""
import math

import numba as nb
import numpy as np

LMAX = 10
W = LMAX + 1                 # slots o = n-LMAX .. n, ring-indexed by o % W

# outcomes of one history
LEAK_LEFT, LEAK_RIGHT, FISSION, CAPTURE = 0, 2, 3, 4

# score channels: variant x parameter, then the two LR density parameters
CHANNELS = ["a:theta_a", "a:theta_b", "b:theta_a", "b:theta_b",
            "c:theta_a", "c:theta_b", "rho_f", "rho_r"]
NCH = len(CHANNELS)


# --------------------------------------------------------------------------
# counter-based RNG (splitmix64), one stream per (seed, generation, particle)
# --------------------------------------------------------------------------
@nb.njit(inline="always")
def _mix64(z):
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return z ^ (z >> np.uint64(31))


@nb.njit(inline="always")
def stream(seed, gen, pid):
    z = _mix64(np.uint64(seed) + np.uint64(0x9E3779B97F4A7C15))
    z = _mix64(z ^ (np.uint64(gen) * np.uint64(0xD1B54A32D192ED03)))
    return _mix64(z ^ (np.uint64(pid) * np.uint64(0x8CB92BA72F3D8DD7)))


@nb.njit(inline="always")
def rnd(s):
    """Return (new state, uniform in (0, 1))."""
    s = s + np.uint64(0x9E3779B97F4A7C15)
    z = _mix64(s)
    return s, ((z >> np.uint64(11)) + 0.5) * (1.0 / 9007199254740992.0)


# --------------------------------------------------------------------------
# one history
# --------------------------------------------------------------------------
@nb.njit(inline="always")
def transport_one(i, x, mu, r, act, s, st, ss, sf, nu, a, b,
                  Da, ga, Db, Dc, gc, Dr, deriv):
    """Track particle i from (x, mu) in region r until it dies.

    Derivative state lives in row i of the slot arrays and is updated in place;
    act is the ring index of the active (own-generation) slot. Returns
    (outcome, x, number of sites, rng state). On FISSION the slots hold the
    state after the fission collision's piece update.
    """
    gb0 = 0.0     # variant (b): g of the active slot only, zero at birth
    gb1 = 0.0
    while True:
        # bounding surface of this piece: 0 -> x=0, 1 -> x=a, 2 -> x=b
        if mu > 0.0:
            if r == 0:
                xsurf = a
                sid = 1
            else:
                xsurf = b
                sid = 2
        else:
            if r == 0:
                xsurf = 0.0
                sid = 0
            else:
                xsurf = a
                sid = 1
        l = (xsurf - x) / mu
        sig = st[r]
        s, xi = rnd(s)
        d = -math.log(xi) / sig
        collided = d < l

        if deriv:
            uu = d / l if collided else 1.0
            for p in range(2):
                # V_n = -(dF/dtheta)/|gradF| = 1 on the surface that is theta_p
                V = 1.0 if sid == p + 1 else 0.0
                # (a) per-origin slots: active slot sees V, passive slots V = 0
                for sl in range(Da.shape[2]):
                    Vs = V if sl == act else 0.0
                    g = ga[i, p, sl]
                    if Vs == 0.0 and g == 0.0:
                        continue
                    dl = (Vs - g) / mu
                    if collided:
                        Da[i, p, sl] += -sig * uu * dl + dl / l
                        ga[i, p, sl] = g + mu * uu * dl
                    else:
                        Da[i, p, sl] += -sig * dl
                        ga[i, p, sl] = g + mu * dl
                # (b) g reset at birth, increments to the active slot
                g = gb0 if p == 0 else gb1
                dl = (V - g) / mu
                if collided:
                    Db[i, p, act] += -sig * uu * dl + dl / l
                    g = g + mu * uu * dl
                else:
                    Db[i, p, act] += -sig * dl
                    g = g + mu * dl
                if p == 0:
                    gb0 = g
                else:
                    gb1 = g
                # (c) single inherited g, increments to the active slot
                g = gc[i, p]
                dl = (V - g) / mu
                if collided:
                    Dc[i, p, act] += -sig * uu * dl + dl / l
                    gc[i, p] = g + mu * uu * dl
                else:
                    Dc[i, p, act] += -sig * dl
                    gc[i, p] = g + mu * dl
            # density of region r (LR): attenuation, +1 per collision
            if collided:
                Dr[i, r, act] += 1.0 - sig * d
            else:
                Dr[i, r, act] -= sig * l

        if collided:
            x = x + mu * d
            s, xi = rnd(s)
            if xi * sig < ss[r]:
                s, xi = rnd(s)
                mu = 2.0 * xi - 1.0
            elif xi * sig < ss[r] + sf[r]:
                s, xi = rnd(s)
                n_int = int(nu[r])
                nsite = n_int + (1 if xi < nu[r] - n_int else 0)
                return FISSION, x, nsite, s
            else:
                return CAPTURE, x, 0, s
        else:
            x = xsurf
            if sid == 0:
                return LEAK_LEFT, x, 0, s
            if sid == 2:
                return LEAK_RIGHT, x, 0, s
            r = 1 - r


# --------------------------------------------------------------------------
# one generation
# --------------------------------------------------------------------------
@nb.njit(parallel=True, cache=True)
def transport_generation(n, seed, src_x, sDa, sga, sDb, sDc, sgc, sDr,
                         fx, fnu, Da, ga, Db, Dc, gc, Dr,
                         st, ss, sf, nu, a, b, deriv):
    N = src_x.size
    act = n % Da.shape[2]
    for i in nb.prange(N):
        s = stream(seed, n, i)
        Da[i] = sDa[i]
        ga[i] = sga[i]
        Db[i] = sDb[i]
        Dc[i] = sDc[i]
        gc[i] = sgc[i]
        Dr[i] = sDr[i]
        s, xi = rnd(s)
        mu = 2.0 * xi - 1.0
        x = src_x[i]
        r = 0 if x < a else 1
        out, x, nsite, s = transport_one(i, x, mu, r, act, s, st, ss, sf, nu, a, b,
                                         Da, ga, Db, Dc, gc, Dr, deriv)
        fx[i] = x
        fnu[i] = nsite if out == FISSION else 0


@nb.njit(parallel=True, cache=True)
def score_generation(n, fnu, Da, Db, Dc, Dr, S1, S2):
    """IFP scores of generation n: S1[ch, L] += w D_{n-L}, S2 += (w D)^2, w = 1.

    Only fission events score. Returns the number of fission events (= den).
    """
    N = fnu.size
    Wl = Da.shape[2]
    nL = S1.shape[1]
    for job in nb.prange(NCH * nL):
        ch = job // nL
        L = job % nL
        sl = (n - L) % Wl
        acc = 0.0
        acc2 = 0.0
        for i in range(N):
            if fnu[i] == 0:
                continue
            if ch < 2:
                v = Da[i, ch, sl]
            elif ch < 4:
                v = Db[i, ch - 2, sl]
            elif ch < 6:
                v = Dc[i, ch - 4, sl]
            else:
                v = Dr[i, ch - 6, sl]
            acc += v
            acc2 += v * v
        S1[ch, L] += acc
        S2[ch, L] += acc2
    nf = 0
    for i in range(N):
        if fnu[i] > 0:
            nf += 1
    return nf


@nb.njit(cache=True)
def comb(n, seed, fx, fnu, Da, ga, Db, Dc, gc, Dr,
         src_x, sDa, sga, sDb, sDc, sgc, sDr):
    """Comb the fission bank to N equal-weight sites and build generation n+1.

    Daughters copy the parent's slots at the fission point, then open their own
    active slot (n+1) % W = {0, 0}, which drops slot n - LMAX.
    """
    N = src_x.size
    M = 0
    for i in range(N):
        M += fnu[i]
    s = stream(seed, n, 0xFFFFFFFF)
    s, xi = rnd(s)
    nxt = (n + 1) % Da.shape[2]
    i = 0
    cum = fnu[0]          # sites of parents 0..i
    for k in range(N):
        pos = (k + xi) * M / N
        while cum <= pos:
            i += 1
            cum += fnu[i]
        src_x[k] = fx[i]
        sDa[k] = Da[i]
        sga[k] = ga[i]
        sDb[k] = Db[i]
        sDc[k] = Dc[i]
        sgc[k] = gc[i]
        sDr[k] = Dr[i]
        sDa[k, :, nxt] = 0.0
        sga[k, :, nxt] = 0.0
        sDb[k, :, nxt] = 0.0
        sDc[k, :, nxt] = 0.0
        sDr[k, :, nxt] = 0.0
    return M


def run(seed, N, n_inactive, n_active, a, b, fuel, refl, deriv=True):
    """Power iteration. Returns per-generation k and, if deriv, IFP sums.

    S1[g, ch, L], S2[g, ch, L], den[g] for active generation g; L = 0..LMAX.
    """
    st = np.array([fuel["st"], refl["st"]])
    ss = np.array([fuel["ss"], refl["ss"]])
    sf = np.array([fuel["sf"], refl["sf"]])
    nu = np.array([fuel["nu"], refl["nu"]])
    shp = (N, 2, W) if deriv else (N, 2, 1)
    src = [np.zeros(N)] + [np.zeros(shp) for _ in range(4)] + [np.zeros((N, 2)), np.zeros(shp)]
    rec = [np.zeros(N), np.zeros(N, dtype=np.int64)] + [np.zeros(shp) for _ in range(4)] \
        + [np.zeros((N, 2)), np.zeros(shp)]
    # initial source: uniform in the fuel, scaled with a (keeps CRN meaningful)
    rng = np.random.default_rng(seed + 7919)
    src[0][:] = a * rng.random(N)
    ngen = n_inactive + n_active
    kgen = np.zeros(ngen)
    S1 = np.zeros((n_active, NCH, LMAX + 1))
    S2 = np.zeros((n_active, NCH, LMAX + 1))
    den = np.zeros(n_active)
    for n in range(ngen):
        transport_generation(n, seed, *src, *rec, st, ss, sf, nu, a, b, deriv)
        if deriv and n >= n_inactive:
            g = n - n_inactive
            den[g] = score_generation(n, rec[1], rec[2], rec[4], rec[5], rec[7],
                                      S1[g], S2[g])
        M = comb(n, seed, *rec, *src)
        kgen[n] = M / N
    return dict(kgen=kgen[n_inactive:], S1=S1, S2=S2, den=den)
