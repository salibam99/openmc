"""
OPENMC-STYLE (per-piece) VARIANT. Each advance = one piece, ending at a collision
or at the nearest boundary, exactly like OpenMC's event loop. Track-length scoring
uses: l_v * D (current), face terms with the CURRENT piece start motion, and a
Leibniz endpoint term +dl added to the bin containing the piece end whenever the
piece ends at a surface (crossing or leak). Without that endpoint term the bin
containing a moving interface is badly wrong (-1.44 vs -0.61 in this test).

Toy check of Deng et al. (SIGGRAPH '26) distance-space reparameterization,
applied as a *post-hoc* derivative on an ordinary analog random walk
(i.e. what you could do with OpenMC track / collision_track data, or inside
OpenMC's existing per-particle derivative accumulator).

Geometry: slab, material A on [0, a), material B on [a, L), vacuum outside.
Beam at x=0, mu=1. One group, isotropic scattering, analog absorption.
theta = a (interface position).  Also sigma_t of B (albedo fixed) via LR.

Responses: transmission T (leak through x=L), track-length flux in bins.
Reference: finite differences with common random numbers + polynomial fit.
"""
import math
import numpy as np
import numba as nb


@nb.njit(cache=True)
def score_flux(x, mu, s0, s1, Dw, g, nbins, w, flux, dflux, dend):
    # track-length part: overlap of [X0,X1] with bins, weighted by d log S
    X0 = x + mu * s0
    X1 = x + mu * s1
    lo = min(X0, X1)
    hi = max(X0, X1)
    amu = abs(mu)
    j0 = max(0, int(lo / w))
    j1 = min(nbins - 1, int(hi / w))
    for j in range(j0, j1 + 1):
        ov = min(hi, (j + 1) * w) - max(lo, j * w)
        if ov > 0.0:
            Lb = ov / amu
            flux[j] += Lb
            dflux[j] += Dw * Lb
    # start-point motion term: bin faces crossed in (s0, s1]
    if g != 0.0:
        for jf in range(nbins + 1):
            xf = jf * w
            sf = (xf - x) / mu
            if sf > s0 and sf < s1 - 1e-12:
                if jf < nbins:
                    dflux[jf] += g / amu      # lo face of bin jf
                if jf >= 1:
                    dflux[jf - 1] -= g / amu  # hi face of bin jf-1
    if dend != 0.0:
        Xe = x + mu * s1 - mu * 1e-12
        je = min(nbins - 1, max(0, int(Xe / w)))
        dflux[je] += dend


@nb.njit(cache=True)
def simulate(N, a, L, st_a, c_a, st_b, c_b, nbins, seed):
    np.random.seed(seed)
    w = L / nbins
    T = 0.0
    dT = 0.0      # dT/da
    dTs = 0.0     # dT/dsigma_t^B  (albedo fixed)
    flux = np.zeros(nbins)
    dflux = np.zeros(nbins)
    ps = np.zeros(2)
    dps = np.zeros(2)
    ms = np.zeros(2, dtype=np.int64)
    for n in range(N):
        x = 0.0
        mu = 1.0
        g = 0.0    # dx/da of current start point
        D = 0.0    # d/da log(path density), reparameterized coordinates
        Ds = 0.0   # d/dsigma_B log(path density)
        while True:
            # material-change boundaries ahead along the ray
            inA = x < a
            if mu > 0.0:
                if inA:
                    nbd = 2; ps[0] = a; dps[0] = 1.0; ms[0] = 0; ps[1] = L; dps[1] = 0.0; ms[1] = 1
                else:
                    nbd = 1; ps[0] = L; dps[0] = 0.0; ms[0] = 1
            else:
                if inA:
                    nbd = 1; ps[0] = 0.0; dps[0] = 0.0; ms[0] = 0
                else:
                    nbd = 2; ps[0] = a; dps[0] = 1.0; ms[0] = 1; ps[1] = 0.0; dps[1] = 0.0; ms[1] = 0
            tau = -math.log(np.random.random())
            prev = 0.0; dprev = 0.0
            A = 0.0    # -sum_{m<J} sigma_m dl_m
            As = 0.0   # -(length in B so far)
            collided = False
            m = 0
            for k in range(nbd):
                p = ps[k]; dp = dps[k]; m = ms[k]
                sig = st_a if m == 0 else st_b
                b = (p - x) / mu
                db = (dp - g) / mu          # implicit-function derivative of hit distance
                ell = b - prev
                dell = db - dprev
                if tau < sig * ell:
                    s = tau / sig
                    tend = prev + s
                    score_flux(x, mu, prev, tend, D + A, g + mu*dprev, nbins, w, flux, dflux, 0.0)
                    u = s / ell
                    D += A - sig * u * dell + dell / ell   # attenuation + Jacobian
                    Ds += As
                    if m == 1:
                        Ds += -s + 1.0 / sig
                    x = x + mu * tend
                    g = g + mu * (dprev + u * dell)
                    collided = True
                    break
                else:
                    score_flux(x, mu, prev, b, D + A, g + mu*dprev, nbins, w, flux, dflux, dell)
                    tau -= sig * ell
                    A += -sig * dell
                    if m == 1:
                        As += -ell
                    prev = b; dprev = db
            if not collided:
                if mu > 0.0:
                    T += 1.0
                    dT += D + A
                    dTs += Ds + As
                break
            c = c_a if m == 0 else c_b
            if np.random.random() > c:
                break
            mu = 2.0 * np.random.random() - 1.0
    return T / N, dT / N, dTs / N, flux / (N * w), dflux / (N * w)


if __name__ == "__main__":
    L, st_a, c_a, st_b, c_b = 2.0, 2.0, 0.5, 0.6, 0.9
    nbins = 10
    a0, sb0 = 0.9, st_b   # keep interface mid-bin: bin flux has a kink in a when interface sits on a bin face
    N = 4_000_000

    T, dT, dTs, phi, dphi = simulate(N, a0, L, st_a, c_a, st_b, c_b, nbins, 12345)

    # ---- FD references (common random numbers, quadratic fit) ----
    hs = np.linspace(-0.06, 0.06, 7)
    Ta, phia = [], []
    for h in hs:
        r = simulate(N, a0 + h, L, st_a, c_a, st_b, c_b, nbins, 999)
        Ta.append(r[0]); phia.append(r[3])
    Ta = np.array(Ta); phia = np.array(phia)
    dT_fd = np.polyfit(hs, Ta, 2)[1]
    dphi_fd = np.array([np.polyfit(hs, phia[:, j], 2)[1] for j in range(nbins)])

    Ts = []
    hs2 = np.linspace(-0.06, 0.06, 7)
    for h in hs2:
        Ts.append(simulate(N, a0, L, st_a, c_a, sb0 + h, c_b, nbins, 777)[0])
    dTs_fd = np.polyfit(hs2, np.array(Ts), 2)[1]

    print(f"T = {T:.5f}")
    print(f"dT/da      reparam = {dT:+.5f}   FD = {dT_fd:+.5f}")
    print(f"dT/dsig_B  LR      = {dTs:+.5f}   FD = {dTs_fd:+.5f}")
    print("bin  phi      dphi/da(track-len reparam)  FD")
    for j in range(nbins):
        print(f"{j:3d}  {phi[j]:.4f}   {dphi[j]:+.4f}   {dphi_fd[j]:+.4f}")
