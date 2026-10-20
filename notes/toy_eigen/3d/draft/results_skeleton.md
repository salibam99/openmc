# Step 0b: 3D eigenvalue toy, variant (a) vs (c)

@@SUMMARY@@

---

## 1. Setup

One energy group, isotropic scattering, analog transport, vacuum outside
(`geom3d.py`, `mc3d.py`). Fuel: Σt 1.0, Σs 0.6, Σf 0.22, ν 2.4 [1/cm].
Reflector: Σt 1.5, Σs 1.4. Units are cm.

| case | geometry | parameters (MC) | k |
|---|---|---|---|
| 1 | bare fuel sphere, R = 4.2 | R, tx (translation of the sphere) | 1.00018 (S_N 1.00021) |
| 2 | fuel sphere a = 3.5 in reflector shell b = 5.5 | a, b | 1.02306 (S_N 1.02309) |
| 3 | fuel sphere a = 3.5 at x0 = 1.0 in reflector cylinder c = 5.5, z ∈ [−6, 6] | a, tx, ty, tz, c, zlo, zhi | 1.02467 |

All cases also give ρ_f and ρ_r (density scalings, by LR) on the same walk.

- **Estimators.** (a) and (c) as in `../derivation.md`, plus LR, all carried on the
  same random walk.
  - Per-piece rule: dl = (V_n − n·g)/(n·ω). A crossing adds D += −Σ dl and g += ω dl.
    A collision at fraction u of the piece adds D += −Σ u dl + dl/l and g += ω u dl.
  - (a): per-origin slots {D_o, g_o}. The own slot is active (V_n of the hit
    surface); slots o < n are passive (V_n = 0).
  - (c): one inherited g whose increments go into the current slot.
  - IFP with OpenMC indexing (L = 1 uses the parent's slot). Denominator: w per
    fission event. L = 1, 2, 3, 5, 8, 10, 12 (all L ≤ 12 are stored).
- **Runs.**
  - Derivative runs: 100 independent seeds per case, N = 1e5 per generation,
    100 inactive + 200 active (H = 2e7 active histories per seed), all
    parameters and both variants in one run.
  - FD runs: CRN central differences, 100 seeds per parameter (seeds 1001–1100),
    same N and generations, h = 0.1 cm, k only.
  - 20 threads.
- **Cost per seed.**

  | | derivative run (all channels) | FD pair (k only, per parameter) |
  |---|---|---|
  | case 1 | 18.5 s | @@FDT1@@ s |
  | case 2 | 18.2 s | @@FDT2@@ s |
  | case 3 | 55.5 s | @@FDT3@@ s |

### Which reference each number is compared against
- **S_N exact** (`sn_sphere.py`, cases 1 and 2): d ln k/dθ from the spherical S_N
  eigenvalue (DD in r and μ, mesh-converged to 1e-8). On the bare sphere it was
  validated against the exact slab-odd-mode equivalence to 2e-7. This is the
  target of (c), and of every variant as L → ∞.
- **S_N finite-L** (`sn_sphere_frozen.py`, cases 1 and 2, radial parameters and
  densities): d ln Q_L/dθ, the frozen-importance functional. This is the exact
  target of (a) and LR at each finite L. For case 1, (c) R is exactly LR ρ_f/R
  pathwise (global dilation), so its target is also S_N finite-L ρ_f/R.
- **CRN FD** (MC, all cases, parameters R; a, b; a, tx, c, zhi): central difference
  of k with common random numbers, 100 seeds, h = 0.1 cm. Its O(h²) truncation,
  estimated in S_N, is +0.03% (R), +0.04% (a), +0.30% (b), far below its
  statistical error. It is the only reference for case 3.
- **Symmetry** (exact zeros): case 1 tx; case 3 ty, tz, and k_zlo + k_zhi.
- **Identities:**
  - case 3 scaling identity, a k_a + x0 k_tx + c k_c + zlo k_zlo + zhi k_zhi =
    k_ρf + k_ρr (exact for the true derivatives, against the LR densities on the
    same walk);
  - cases 1 and 2: R k_R = k_ρf, and a k_a + b k_b = k_ρf + k_ρr.

The tables are in `runs/analysis_log.txt`, which `analyze3d.py` regenerates; this
file quotes them.

@@BODY@@
