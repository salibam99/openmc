# Step 0 results: offline k-eigenvalue toy (1D slab), distance reparameterization + IFP

Definitions are in derivation.md: estimator, variants (a)/(b)/(c), window L, and
§5.1, which was worked out after the first runs. Units: cm and 1/cm. All MC error
bars are ±1σ from independent seeds. Batch statistics are shown only as a
diagnostic in §6.

## 0. Summary

- **(a) per-origin slots is unbiased at every L.** It matches an independent S_N
  calculation of its exact finite-L target (one generation's kernel perturbed,
  source and L-generation importance frozen in physical space) within 1.5σ for
  all L = 0..10. It converges to d ln k/da, but slowly. Its target is 3.1% low at L = 3, 0.29% low
  at L = 8 and 0.10% low at L = 10. The variance is heavy-tailed, in effect
  log-infinite. One history in 4.8e9 moved a seed's L = 5 estimate by 10 typical
  σ (§6).
- **(b) g reset at birth is biased.** For θ_a it converges to 0.1297, against
  0.08415 exact (+54%, hundreds of σ), for every L. Dropping the ∇ψ†·g term is
  not a small effect.
- **(c) single inherited g matches its finite-L target within 1.5σ at all L and
  converges much faster.** It is within 0.1% of exact for L ≥ 2, has about 3.5×
  lower per-history variance than (a) at L = 8, and no heavy tails. But in this 1D problem (c)
  is degenerate: in stationarity it is pathwise identical to a fixed
  combination of density LR scores (§7). Its good behaviour here therefore says
  nothing yet about 3D.
- θ_b (outer boundary) cannot tell the variants apart in 1D. All three are
  identical to round-off, and pathwise equal to LR(ρ_r)/(b−a). It is correct at
  every finite L (within 1.5σ).
- Scaling identity: S_N satisfies it to 3e-11. MC LR (right-hand side) converges
  to S_N. (a) violates it at finite L by exactly the amount the finite-L S_N
  predicts (within 1.7σ). (c) satisfies it pathwise, so the test cannot check (c).
- Nothing disagreed by more than 3σ with its correct reference. The largest
  deviation was −2.0σ (ρ_f, L = 0), highly correlated across L. Before S128 was
  added to the finite-L reference it was −2.6σ.
- A 3D side probe (spheres and nested cubes): passive g does not blow up. It
  decays only by a factor 0.86–0.92 per generation (in 1D it is zeroed at every
  crossing), with heavy-tailed |g|. This is the main open risk for C++.

## 1. Model and run parameters

| region | interval [cm] | Σ_t | Σ_s | Σ_f | ν̄ (integer sampling) |
|---|---|---|---|---|---|
| fuel | [0, a], a = 3.0 | 1.0 | 0.6 | 0.22 | 2.4 |
| reflector | [a, b], b = 4.5 | 1.5 | 1.4 | 0 | – |

Vacuum at x = 0 and x = b. One group, isotropic scattering, analog. S_N:
k = 1.018887, dominance ratio |k₂|/k = 0.569.

MC: N = 1e5 per generation, combed to N equal weights, 60 inactive + 400 active
generations. 120 independent seeds for the sensitivities (4.8e9 active histories)
and 60 seeds for CRN central differences. All variants, L = 0..10 and all
parameters come from the same random walks (derivation.md §4).

Cost on 20 threads: 15 s per seed with all derivative channels, about 3.2e-7 s
per history wall-clock. A k-only run is about 1.0e-7 s per history.

## 2. References

### 2.1 S_N (sn.py)

Gauss–Legendre S_N, diamond difference without fixup, uniform cells per region
(they stretch with a and b), and a dense exact inner solve. Central differences
with steps 2e-3 and 1e-3 relative, Richardson in h.

| S_N | cells/region | k | dk/da | dk/db | dk/dρ_f | dk/dρ_r |
|---|---|---|---|---|---|---|
| 32 | 250 | 1.0188571599 | 0.08575420 | 0.01483299 | 0.30176159 | 0.02224949 |
| 32 | 500 | 1.0188581943 | 0.08575435 | 0.01483288 | 0.30176166 | 0.02224931 |
| 32 | 1000 | 1.0188584529 | 0.08575438 | 0.01483285 | 0.30176168 | 0.02224927 |
| 32 | 2000 | 1.0188585176 | 0.08575439 | 0.01483284 | 0.30176169 | 0.02224926 |
| 16 | 2000 | 1.0187681251 | 0.08579657 | 0.01484277 | 0.30191803 | 0.02226416 |
| 64 | 2000 | 1.0188801465 | 0.08574451 | 0.01483041 | 0.30172475 | 0.02224561 |
| 128 | 2000 | 1.0188854511 | 0.08574211 | 0.01482981 | 0.30171574 | 0.02224471 |
| 128 | 1000 | 1.0188853865 | 0.08574210 | 0.01482981 | 0.30171573 | 0.02224472 |
| 256 | 1000 | 1.0188867008 | 0.08574150 | 0.01482966 | 0.30171350 | 0.02224450 |
| **Richardson S128/S256** | | **1.0188871** | **0.0857413** | **0.0148296** | **0.3017128** | **0.0222444** |

The mesh is converged to about 1e-8. Angular error is O(1/N²): the S16→32→64→128
differences shrink by about 4× each step. The extrapolation correction is 2e-7 on
dk/da, negligible against the MC errors (≥ 3e-5).

Reference: **d ln k/da = 0.084152, d ln k/db = 0.014555, d ln k/dρ_f = 0.29612,
d ln k/dρ_r = 0.021832.** Scaling identity in S_N:
a·dk/da + b·dk/db = 0.32395718 and Σρ dk/dρ = 0.32395718 (difference −3e-11).

### 2.2 Finite-L S_N reference (sn_frozen.py), added during the investigation

Variant (a) at window L estimates d ln Q_L/dθ, with Q_L(θ) = ⟨1, H₀ᴸ H_θ ψ⟩. Only
one generation is perturbed; ψ and the L-step importance are frozen in physical
space (derivation.md §4). This does not depend on the reparameterization, so S_N
computes it directly. It uses a fixed uniform mesh with void beyond b, a or b
moved by whole cells (h = 2, 4 cells, Richardson), dense H.

Convergence: dx = 0.01, 0.005, 0.0025 at S32/S64, and S128 at dx = 0.005. Changes
are ≤ 2e-5 at S64 and ≤ 1e-6 between meshes. "Best" = S64/S128 Richardson + mesh
correction. For L = 10 it reaches d ln Q/da = 0.08407 (exact 0.08415), d ln Q/db
= 0.01455, d ln Q/dρ_f = 0.29612. So it tends to the reference as L grows.

The finite-L prediction for (c) on θ_a follows from §7:
[LR_f + LR_r − b·θ_b](L)/a.

## 3. Fixed-source sanity check (fixed_source_check.py)

This is the problem of notes/reparam_slab_check_openmc_style.py, run through the
eigenvalue toy's own per-piece kernel (mc.transport_one): beam at x = 0, A on
[0, 0.9) with Σ_t = 2, c = 0.5; B on [0.9, 2) with Σ_t = 0.6, c = 0.9. 4e7
histories, iid, so within-run errors are valid here. FD is CRN quadratic fit over
h ∈ [−0.06, 0.06], 4 independent replicas of 1e7 each.

| quantity | toy kernel | FD (this kernel) | reference script (4e6, its FD) |
|---|---|---|---|
| T | 0.15700 ± 0.00006 | – | 0.15675 |
| dT/da | −0.22488 ± 0.00010 | −0.22494 ± 0.00013 | −0.22430 (FD −0.22500) |
| dT/dσ_B (LR) | −0.09113 ± 0.00009 | −0.09118 ± 0.00010 | −0.09121 (FD −0.09188) |
| dT/dL (outer boundary) | −0.04971 ± 0.00005 | −0.04975 ± 0.00008 | not in script |
| dR/da | −0.01490 ± 0.00009 | −0.01488 ± 0.00015 | – |
| dR/dL | +0.00419 ± 0.00001 | +0.00419 ± 0.00003 | – |
| dR/dσ_B | +0.00769 ± 0.00003 | +0.00771 ± 0.00005 | – |

(a), (b) and (c) give bit-identical fixed-source scores (difference 0), as they
must with a single generation. The reference script's own values carry about
√10 × larger errors (4e6 histories), so its dT/da = −0.22430 is within 2σ. The
script was run unmodified (a copy in the scratchpad); its track-length flux bins
were not reproduced, because the eigenvalue toy has no track-length scoring.

## 4. Main results: d ln k/dθ vs L

![d ln k/da](fig_dlnk_da.png)

### θ₁ = a (interface)

z is measured against the exact reference and against the finite-L target of the
same variant (§2.2).

| L | (a) MC | z exact | (a) S_N finite-L | z | (b) MC | z exact | (c) MC | z exact | (c) S_N finite-L | z |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | +0.13124 ± 0.00011 | +429 | +0.13130 | −0.5 | +0.13124 ± 0.00011 | +429 | +0.08131 ± 0.00001 | −235 | +0.08133 | −1.5 |
| 1 | +0.08588 ± 0.00006 | +30 | +0.08580 | +1.3 | +0.12858 ± 0.00012 | +374 | +0.08353 ± 0.00001 | −43 | +0.08355 | −1.3 |
| 2 | +0.08147 ± 0.00011 | −25 | +0.08141 | +0.5 | +0.12873 ± 0.00020 | +226 | +0.08407 ± 0.00002 | −4.5 | +0.08409 | −1.1 |
| 3 | +0.08146 ± 0.00011 | −26 | +0.08155 | −0.9 | +0.12903 ± 0.00015 | +309 | +0.08419 ± 0.00002 | +1.7 | +0.08421 | −1.0 |
| 5 | +0.08276 ± 0.00018 | −7.9 | +0.08300 | −1.3 | +0.12936 ± 0.00013 | +353 | +0.08419 ± 0.00003 | +1.5 | +0.08420 | −0.2 |
| 8 | +0.08395 ± 0.00006 | −3.6 | +0.08391 | +0.6 | +0.12968 ± 0.00015 | +301 | +0.08419 ± 0.00003 | +1.1 | +0.08416 | +0.7 |
| 10 | +0.08404 ± 0.00006 | −1.8 | +0.08407 | −0.6 | +0.12966 ± 0.00016 | +281 | +0.08418 ± 0.00003 | +1.0 | +0.08416 | +0.8 |
| exact | 0.084152 | | | | | | | | | |

L = 0 means no adjoint weighting (derivation.md §4). (a) and (b) coincide there,
but (c) does not, because its own slot already holds the inherited-g increments.

### θ₂ = b, and density LR (variants identical)

![d ln k/db and LR](fig_dlnk_db_rho.png)

| L | d ln k/db | z exact | finite-L | z | d ln k/dρ_f (LR) | z exact | finite-L | z | d ln k/dρ_r (LR) | z exact | finite-L | z |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | +0.01427 ± 0.00001 | −36 | +0.01427 | −0.1 | +0.28672 ± 0.00003 | −338 | +0.28678 | −2.0 | +0.02140 ± 0.00001 | −36 | +0.02140 | −0.1 |
| 1 | +0.01421 ± 0.00001 | −38 | +0.01422 | −0.2 | +0.29322 ± 0.00004 | −77 | +0.29328 | −1.6 | +0.02132 ± 0.00001 | −38 | +0.02132 | −0.2 |
| 2 | +0.01430 ± 0.00001 | −23 | +0.01431 | −0.4 | +0.29512 ± 0.00005 | −22 | +0.29519 | −1.5 | +0.02145 ± 0.00002 | −23 | +0.02146 | −0.4 |
| 3 | +0.01438 ± 0.00001 | −15 | +0.01439 | −0.8 | +0.29572 ± 0.00005 | −7.6 | +0.29581 | −1.8 | +0.02158 ± 0.00002 | −15 | +0.02159 | −0.8 |
| 5 | +0.01449 ± 0.00001 | −4.6 | +0.01450 | −0.7 | +0.29603 ± 0.00006 | −1.4 | +0.29609 | −0.8 | +0.02173 ± 0.00002 | −4.6 | +0.02174 | −0.7 |
| 8 | +0.01451 ± 0.00002 | −2.0 | +0.01454 | −1.5 | +0.29610 ± 0.00007 | −0.3 | +0.29612 | −0.3 | +0.02177 ± 0.00003 | −2.0 | +0.02182 | −1.5 |
| 10 | +0.01452 ± 0.00002 | −1.5 | +0.01455 | −1.4 | +0.29611 ± 0.00008 | −0.1 | +0.29612 | −0.1 | +0.02178 ± 0.00003 | −1.5 | +0.02183 | −1.4 |
| exact | 0.014555 | | | | 0.29612 | | | | 0.021832 | | | |

For θ_b, max |(a)−(b)| = max |(a)−(c)| = 3.5e-18 over all seeds and L. The θ_b and
ρ_r columns have identical z-scores because the two estimators are pathwise
proportional (§7). Estimates at different L come from the same runs and are
strongly correlated, so the −1.4/−1.5σ at L = 8, 10 is one fluctuation, not two.

### dk/dθ = k · d ln k/dθ (per seed, then averaged)

| | S_N | (a) L=8 | (a) L=10 | (b) L=10 | (c) L=8 | (c) L=10 |
|---|---|---|---|---|---|---|
| dk/da | 0.085741 | 0.08553 ± 0.00006 | 0.08562 ± 0.00006 | 0.13211 ± 0.00016 | 0.08578 ± 0.00003 | 0.08577 ± 0.00003 |
| dk/db | 0.014830 | 0.01479 ± 0.00002 | 0.01479 ± 0.00002 | = | = | = |
| dk/dρ_f | 0.301713 | 0.30169 ± 0.00007 | 0.30170 ± 0.00008 | | | |
| dk/dρ_r | 0.022244 | 0.02218 ± 0.00003 | 0.02219 ± 0.00003 | | | |

k_MC = 1.018887 ± 0.000017 against S_N 1.018887 (z = −0.01).

## 5. Scaling identity a ∂k/∂a + b ∂k/∂b = Σ_m ρ_m ∂k/∂ρ_m

| L | (a) LHS | RHS (LR) | (a) LHS − RHS | finite-L S_N prediction of the gap | z |
|---|---|---|---|---|---|
| 0 | 0.46658 ± 0.00033 | 0.31394 ± 0.00003 | +0.15264 ± 0.00032 | +0.15274 | −0.3 |
| 1 | 0.32767 ± 0.00018 | 0.32048 ± 0.00004 | +0.00718 ± 0.00017 | +0.00690 | +1.7 |
| 2 | 0.31460 ± 0.00032 | 0.32255 ± 0.00005 | −0.00795 ± 0.00031 | −0.00818 | +0.7 |
| 3 | 0.31494 ± 0.00032 | 0.32329 ± 0.00006 | −0.00835 ± 0.00032 | −0.00812 | −0.7 |
| 5 | 0.31938 ± 0.00053 | 0.32377 ± 0.00007 | −0.00438 ± 0.00054 | −0.00368 | −1.3 |
| 8 | 0.32314 ± 0.00017 | 0.32388 ± 0.00008 | −0.00073 ± 0.00014 | −0.00076 | +0.2 |
| 10 | 0.32345 ± 0.00020 | 0.32389 ± 0.00009 | −0.00045 ± 0.00016 | −0.00025 | −1.2 |
| S_N exact | 0.32396 | 0.32396 | −3e-11 | | |

- (a) violates the identity at finite L by 5–43σ for L = 1–8. The investigation
  shows this is correct behaviour, not a wrong term. The geometric and LR
  estimators are first-order perturbation formulas with the same approximate
  importance ψ†_L. The two operator derivatives differ by a commutator
  [H, x·∇], which vanishes only against the exact ψ†. The finite-L S_N
  reproduces the gap at every L.
- (c): LHS − RHS is below 3e-13 for every seed and every L (the z-scores of ±1
  are round-off noise over round-off). This is a pathwise identity (§7), so the
  test does not discriminate (c).
- The RHS (LR) converges to S_N: 0.32389 ± 0.00009 at L = 10, against 0.32396.

## 6. Variance per history, tails, error bars

Per-history variance σ²_hist = Var_seeds(run estimate) × H, with H = 4e7 active
histories per seed. This is the number for later FOM comparisons: σ of a run with
H histories is √(σ²_hist/H). Units 1/cm² (dimensionless for ρ).

![variance](fig_variance.png)

| channel | L=0 | L=1 | L=2 | L=3 | L=5 | L=8 | L=10 | eff/naive (L=8) |
|---|---|---|---|---|---|---|---|---|
| (a) θ_a | 57.7 | 16 | 55.5 | 53.8 | **149** | 16.2 | 19.4 | 0.44 |
| (b) θ_a | 57.7 | 67.9 | 187 | 101 | 78.7 | 110 | 126 | 5.75 |
| (c) θ_a | 0.70 | 1.03 | 1.52 | 2.09 | 3.24 | 4.59 | 5.16 | 6.26 |
| θ_b (all) | 0.31 | 0.38 | 0.58 | 0.65 | 1.08 | 1.91 | 2.40 | 6.55 |
| ρ_f (LR) | 3.71 | 6.87 | 9.81 | 13.5 | 18.3 | 23.4 | 27.3 | 6.17 |
| ρ_r (LR) | 0.69 | 0.85 | 1.31 | 1.46 | 2.43 | 4.29 | 5.41 | 6.55 |

Relative per-history variance (σ_hist/value)² at L = 8: (a) θ_a 2.3e3,
(c) θ_a 6.5e2, θ_b 9.0e3, ρ_f 2.7e2, (b) θ_a 1.6e4.

"eff/naive" divides by the per-history variance an iid-history formula would give
for the same ratio estimator.
- About 6 for (c), θ_b and LR. Within a generation, all descendants of one
  ancestor carry the same D_{n−L}, so the scores are clustered.
- 0.44 for (a). See the next item.

Heavy tails in (a) and (b).
- The excess kurtosis of the per-generation estimates is +221 for (a) and +220
  for (b), against about 0 for every other channel.
- Example, seed 40: one history in generation 395 had a collision extremely close
  to x = a. The term dl/l = (V − g)/s with g = 0 in the active slot gave D of order
  1e5. The same history also wrote about −D·g_o into its passive slots.
- (a) scores those slots at different generations, with different progeny
  counts, so the ± spikes do not cancel. They appear as alternating ±0.0105,
  ±0.021 shifts in generations 395–399 of that run, for L = 0..5. Some of the
  partners fell after the last generation.
- This one seed lifts the (a), L = 5 per-history variance from 15.6 (without it)
  to 149.
- The same mechanism explains (a)'s lag-1 autocorrelation of −0.41 between
  generation estimates, and why batch error bars overestimate (a)'s error by 3×.
- The variance is log-divergent, as warned in the notes. In practice, (a)'s
  variance estimates are reliable only to within a factor of about 3, and its
  error bars should be treated as indicative.
- (c) has no such spikes in 1D. Its g equals the stretching field, so
  V − g(P) = s/a near the surface, and dl/l stays bounded.

Error-bar methodology. Batch-SE (from the 400 generations of one run) divided by
the seed SD:
- about 1.0 (0.9–1.1) for (c), θ_b and LR. Generation estimates are
  uncorrelated at lag 1, so for these channels batch errors would have been fine;
- 0.5–0.9 for (b) at small L, so batches underestimate;
- up to 3.2 for (a) (anticorrelation), so batches overestimate.

The variance estimates are stable for all channels except (a): first half vs
second half of the seeds gives ratios 0.9–1.2, and half-length runs give 0.9–1.2.
Per-seed estimates pass Shapiro–Wilk (p ≥ 0.12), and the largest |z| among 120
seeds is 3.1.

Figure of merit vs CRN central differences. FOM = 1/(σ²T), per parameter. The
derivative run produces all 4 parameters and all L at once; CRN-FD needs 2 runs
per parameter.

| | (c) L=8 | (a) L=8 | CRN-FD (h = 0.1 cm) |
|---|---|---|---|
| dk/da FOM [1/s] | 5.7e5 (2.1× FD) | 1.6e5 (0.6× FD) | 2.7e5, plus an O(h²) bias of +3.4e-5 (S_N) |
| dk/db FOM [1/s] | 1.4e6 (5.7× FD) | same | 2.4e5, plus bias +4.8e-5 |

The CRN gain over independent runs is only about 2× in eigenvalue mode, because
combing decorrelates the perturbed and unperturbed populations. The timings are
toy timings, so take the ratios as rough.

## 7. 1D degeneracies: read before generalizing

- Fission sites are only in the fuel, and reaching one requires crossing x = a.
  That crossing zeroes the θ_b displacement. So all passive g_b are zero and
  (a) = (b) = (c) for θ_b.
- More strongly, inside the reflector g_b = (x − a)/(b − a) exactly, for every
  history (it holds piece by piece). So D_b = D_ρr/(b − a) pathwise, which is the
  1D optical-thickness identity (b − a)·∂k/∂b = ρ_r ∂k/∂ρ_r.
- Variant (c), θ_a: the never-reset g converges to the piecewise-linear stretching
  field, x/a in the fuel and (b − x)/(b − a) in the reflector. Every crossing sets
  it exactly and collisions contract the error by (1 − u). Then
  D_a^(c) = D_ρf/a − D_ρr/(b − a) per slot. Checked: (c) at L = 0 is 0.08131, and
  0.28672/3 − 0.02140/1.5 = 0.08131.
- So in 1D, (c) is LR in disguise. It inherits LR's variance and LR's fast
  L-convergence. Dilating the whole system reproduces Σρ ∂/∂ρ exactly
  (derivation.md §5.1). None of this holds in 3D, where g is not a function of
  position.
- What 1D does test for (c): its expectation is the frozen functional in
  coordinates moving with the stationary displacement field ḡ (§5.1 of the
  derivation), and its limit is exact. What it cannot test: its variance and
  L-convergence when ḡ is not a smooth field.

## 8. 3D side probe: does passive g decay? (g_growth_3d.py)

Same cross sections. Concentric spheres (R = 3, 4.5) and nested cubes (half-widths
3, 4.5), analog with combing. A random unit g is placed on every source particle,
then propagated passively (V_n = 0) and inherited across fissions (never reset).
RMS |g| over fission sites:

| generation | 0 | 2 | 5 | 10 | 15 | 20 | 24 | decay factor / generation |
|---|---|---|---|---|---|---|---|---|
| spheres | 0.92 | 0.71 | 0.45 | 0.23 | 0.13 | 0.053 | 0.031 | 0.865 |
| cubes | 1.16 | 1.05 | 0.84 | 0.54 | 0.38 | 0.23 | 0.18 | 0.924 |

The maximum |g| among 2e5 sites stays in the range 3–77. Passive g decays: no
positive Lyapunov exponent is visible in these geometries. But the decay is slow,
because a crossing only resets the normal component and tangential components
survive. In 3D, a passive slot therefore still carries O(1) memory after
L = 5–10 generations, and (c)'s g is a random, heavy-tailed field rather than
the smooth map it is in 1D. This probe propagated g only, not D or k
sensitivities.

## 9. Assessment

Which variants are unbiased:
- (a) per-origin slots: unbiased for its finite-L target at every L, verified by
  an independent deterministic calculation. That target converges to d ln k/dθ
  with an IFP truncation error that is large at small L: −3.1% at L = 3, −0.29%
  at L = 8, −0.10% at L = 10 (S_N finite-L), decaying roughly like DR^L with
  DR = 0.57. The
  truncation is in the importance gradient, which carries about 50% of
  d ln k/da (compare (b)).
- (b): biased by +54% on θ_a for every L. Rejected.
- (c): unbiased for its finite-L target. Its limit is exact (derivation.md §3 and
  §5.1). In this toy it converges faster than (a): 0.1% for L ≥ 2 versus L ≈ 10,
  because in 1D it reduces to LR. Its derivation and its exact limit are general;
  its variance and convergence advantages are not yet shown to be.

Convergence in L, compared with |exact| for θ_a and in units of the S_N value:
(a) needs L ≈ 10 for about 0.1%; (c) is within 0.1% from L = 2; LR(ρ_f) needs
L ≈ 5; θ_b and ρ_r need L ≈ 8. All are monotone or flat from L = 3 on in this
problem.

Variance per history, θ_a at L = 8: (c) about 4.6, (a) about 16 but heavy-tailed
and unstable. For the well-behaved channels ((c), θ_b, LR) the variance grows
roughly linearly with L, 7–8× from L = 0 to L = 10. That is the usual cost of
IFP.

Before C++, it would be worth deciding:
1. Which variant to implement. (a) is the clean "exact frozen functional" and its
   memory is bounded (slots die after L), but it has heavy tails and needs L + 1
   g-vectors per parameter per particle. (c) needs one g-vector per parameter
   plus L + 1 scalars, has lower variance here and is exact in the limit. Its 3D
   behaviour is the open question: g decays at only about 0.9 per generation and
   is heavy-tailed.
2. A 3D toy (concentric spheres, θ = radius), against 1D-spherical S_N or
   high-statistics CRN FD, is the obvious next discriminating test of (a) vs (c).
   It is not done here, as requested.
3. The log-divergent dl/l variance of (a) is caused by g at a piece start being
   inconsistent with V_n on its bounding surface. (c) shows that a g field that
   matches V at the surface removes the 1/s spikes. This may be worth exploring
   as a variance-reduction idea (a smooth extension of the surface velocity used
   as the initial g, with the corresponding correction), but it is not derived
   here.

## 10. Files and reproduction

| file | content |
|---|---|
| derivation.md | estimator, slot rules, variants, finite-L characterization |
| model.py | cross sections and geometry |
| mc.py | numba MC: per-piece rules, slots for (a)/(b)/(c), LR, IFP scoring, combing, counter-based RNG |
| sn.py | S_N reference, convergence study → sn_results.json, sn_log.txt |
| sn_frozen.py | finite-L (frozen-importance) S_N reference → sn_frozen_results.json |
| run_eigen.py | seed driver: `sens R N n_in n_act`, `fd R N n_in n_act h` → runs/*.npz (66 MB) |
| fixed_source_check.py | fixed-source check → fixed_source_log.txt / .json |
| analyze.py | tables, figures, analysis_results.json, tables.md, analysis_log.txt |
| g_growth_3d.py | 3D passive-g probe → g_growth_3d_log.txt |

Environment: a venv with numba 0.68, numpy 2.3.5, scipy 1.18, matplotlib 3.11. The
anaconda base has numpy 2.5, which this numba does not support. Commands, run
from notes/toy_eigen:

    python sn.py; python sn_frozen.py; python fixed_source_check.py
    python run_eigen.py sens 120 100000 60 400; python run_eigen.py fd 60 100000 60 400 0.1
    python analyze.py; python g_growth_3d.py
