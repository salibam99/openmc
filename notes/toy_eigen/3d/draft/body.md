## 2. d ln k/dθ vs L

±1σ is the standard seed SE (100 seeds). The robust error bars are in §5. The
reference column says which reference each row is compared against.

### Case 1, bare sphere (reference: S_N exact 0.094291; finite-L S_N for (a))

| | L=1 | L=2 | L=3 | L=5 | L=8 | L=10 | L=12 |
|---|---|---|---|---|---|---|---|
| (a) R | 0.15176(5) | 0.12165(4) | 0.10921(3) | 0.09942(4) | 0.09548(4) | 0.09472(4) | 0.09445(4) |
| finite-L S_N (target of (a)) | 0.15191 | 0.12173 | 0.10925 | 0.09947 | 0.09550 | 0.09476 | 0.09448 |
| (a) vs exact | +61% | +29% | +16% | +5.4% | +1.26% | +0.45% | +0.17% |
| (c) R | 0.09106(1) | 0.09251(1) | 0.09325(2) | 0.09389(2) | 0.09419(2) | 0.09422(3) | 0.09422(3) |
| (c) vs exact | −3.4% | −1.9% | −1.1% | −0.43% | −0.11% | −0.08% | −0.08% |
| (a) tx (symmetry 0) | −2(3)e-5 | −3(3)e-5 | +1(2)e-5 | −1(2)e-5 | −3(2)e-5 | −1(2)e-5 | +2(2)e-5 |

- (a) vs its finite-L target: z = −2.9, −2.3, −1.4, −1.2, −0.3, −0.9, −0.7. At
  L = 1–2 the gap (−1.5e-4, −0.8e-4) is the size of the S_N mesh effect at small L:
  going from dr = 0.01 to 0.005 moves the target by −1.2e-4 at L = 1 and −5e-5 at
  L = 2. At L ≥ 5 the mesh effect is ≤ 1e-5 and |z| ≤ 1.2.
- (c) R is pathwise LR ρ_f/R (global dilation), so its target is finite-L S_N
  ρ_f/R: z = +0.3, +0.3, −0.3, −1.2, −0.6, −1.3, −1.8. LR ρ_f gives the same z.
  (c) tx ≡ 0 (rigid translation).
- vs CRN FD at L = 12: (a) @@C1FDZA@@σ, (c) @@C1FDZC@@σ.

### Case 2, sphere in a reflector shell (reference: S_N exact; finite-L S_N for (a))

| | L=1 | L=2 | L=3 | L=5 | L=8 | L=10 | L=12 |
|---|---|---|---|---|---|---|---|
| (a) a | 0.13379(17) | 0.10776(18) | 0.09738(20) | 0.08951(43) | 0.08758(54) | 0.08762(75) | 0.08751(69) |
| finite-L S_N a | 0.13370 | 0.10760 | 0.09725 | 0.08971 | 0.08708 | 0.08669 | 0.08656 |
| (c) a | 0.08493(11) | 0.08568(14) | 0.08602(13) | 0.08668(28) | 0.08706(44) | 0.08712(50) | 0.08712(40) |
| (c) a vs exact 0.086492 | −1.8% | −0.94% | −0.54% | +0.22% | +0.66% | +0.73% | +0.73% |
| (a) b | 0.00982(6) | 0.00936(9) | 0.00919(9) | 0.00876(22) | 0.00858(32) | 0.00854(36) | 0.00861(23) |
| finite-L S_N b | 0.00983 | 0.00945 | 0.00925 | 0.00909 | 0.00902 | 0.00901 | 0.00901 |
| (c) b | 0.00904(7) | 0.00900(9) | 0.00902(9) | 0.00876(18) | 0.00860(28) | 0.00855(31) | 0.00856(25) |
| (c) b vs exact 0.009010 | +0.3% | −0.1% | +0.1% | −2.8% | −4.6% | −5.1% | −5.0% |

- (a) and LR are within 2σ of their finite-L targets at every L (max |z| = 1.9).
- (c) is within 1.8σ of the exact value for L ≥ 5. Its deviation at L ≤ 3 is its
  truncation bias (a: −1.8% at L = 1, −0.5% at L = 3).
- At L ≥ 5, a and b of both variants are dominated by **one event in seed 54**:
  a re-entrant chord in inactive generation 95 (§4.3). It enters every window
  with L ≥ 5 and shifts (a) a and (c) a up by +6.4e-4 and +3.3e-4, and both b's
  down by about 2e-4. Without that seed, (c) a at L = 12 is 0.08679 and (a) a is
  0.08686. The shift explains the apparent drift of (c) b at L ≥ 5.
- The LR densities (no heavy tail) at L = 12: ρ_f 0.30636(13) vs finite-L
  0.30653, ρ_r 0.04561(6) vs 0.04573 (z −1.3, −1.9).
- Scaling identity a k_a + b k_b = k_ρf + k_ρr (S_N 0.360409):
  - (c) holds to 1e-10 at every L (degenerate, see §3);
  - (a) misses by +4.2e-2 at L = 3, +1.7e-3 at L = 8 and +1.7e-3 ± 1.4e-3 at L = 12.

### Case 3, off-centre sphere in a reflector cylinder (reference: CRN FD; symmetry; scaling identity)

| | L=1 | L=2 | L=3 | L=5 | L=8 | L=10 | L=12 | CRN FD (d ln k) |
|---|---|---|---|---|---|---|---|---|
| (a) a | 0.13529(15) | 0.10908(16) | 0.09878(14) | 0.09131(16) | 0.08881(28) | 0.08849(25) | 0.08787(16) | @@C3FDa@@ |
| (c) a | 0.08658(7) | 0.08714(8) | 0.08759(11) | 0.08804(12) | 0.08821(19) | 0.08821(20) | 0.08813(15) | |
| (a) tx | −0.00395(9) | −0.00407(8) | −0.00390(9) | −0.00410(10) | −0.00424(11) | −0.00418(9) | −0.00423(11) | @@C3FDtx@@ |
| (c) tx | −0.00469(4) | −0.00458(4) | −0.00451(4) | −0.00446(4) | −0.00441(6) | −0.00443(6) | −0.00443(6) | |
| (a) c | 0.00762(4) | 0.00732(4) | 0.00709(5) | 0.00697(5) | 0.00689(6) | 0.00690(7) | 0.00689(7) | @@C3FDc@@ |
| (c) c | 0.00702(4) | 0.00703(4) | 0.00696(5) | 0.00692(6) | 0.00688(6) | 0.00689(6) | 0.00692(6) | |
| (a) zhi | 0.000284(13) | 0.000283(14) | 0.000289(15) | 0.000287(16) | 0.000277(20) | 0.000288(22) | 0.000277(24) | @@C3FDzhi@@ |
| (c) zhi | 0.000217(14) | 0.000252(15) | 0.000273(17) | 0.000273(19) | 0.000280(21) | 0.000288(22) | 0.000275(24) | |

@@C3FDTEXT@@

- **Symmetry.** ty, tz and k_zlo + k_zhi are consistent with 0 at every L for both
  variants: max |z| = 2.2 ((a) ty at L = 5) and 2.1 ((a) k_zlo + k_zhi at L = 5),
  out of 56 checks, all others ≤ 1.7.
- **Scaling identity**, LHS − RHS:

  | L | 1 | 2 | 3 | 5 | 8 | 10 | 12 |
  |---|---|---|---|---|---|---|---|
  | (a) | +0.180 | +0.081 | +0.042 | +0.0124 | +0.0023(6) | +0.0012(5) | −0.0009(5) |

  (c) holds to 1e-10. The (a) defect decays like the finite-L bias of (a) and is
  consistent with 0 (−1.9σ) at L = 12. This is the only identity test of (a) in
  case 3, and (a) passes it.

### dk/dθ
k changes the numbers by at most 2.5%, so only L = 8 and 12 are tabulated (in
`runs/analysis_log.txt`). Examples at L = 12:

| | (a) | (c) | reference |
|---|---|---|---|
| case 1 dk/dR | 0.09446(4) | 0.09424(3) | S_N 0.094311 |
| case 2 dk/da | 0.08953(71) | 0.08912(41) | S_N 0.088488 |
| case 2 dk/db | 0.00881(24) | 0.00876(26) | S_N 0.009218 |
| case 3 dk/da | 0.09004(16) | 0.09030(15) | @@C3FDdka@@ |
| case 3 dk/dtx | −0.00434(12) | −0.00454(7) | @@C3FDdktx@@ |

## 3. Degeneracies of (c) (predicted and confirmed)

Along a global dilation, and for a rigid translation of the whole system, the
(c) g field is exactly the displacement field (g = x/R, or g = e). (c) then
reduces pathwise to LR (or to 0).

- **Case 1:** (c) R = LR ρ_f/R exactly. (c) tx ≡ 0. The (c) g is bounded by 1
  (max |g| = 1).
- **Case 3:** the scaling identity holds for (c) to 1e-10. It therefore tests only
  (a).
- **Not degenerate:** case 2 (a, b) and the individual case-3 parameters. They are
  the real tests of (c).

## 4. Variance and tails

### 4.1 What the "x^-1.8" was
It was the slope of the tail probability (CCDF) P(|dl/l| > x) per collision
piece, not of the density. It is a local slope; the asymptotic behaviour depends
on the geometry:

- **Collision near a moving surface** (piece from inside to the surface,
  s = l|n·ω| the normal distance). P(s < ε) ∝ ε, and dl/l ≈ V_n/s, so
  P(|dl/l| > x) ∝ x^-2 ln(x/x0). The local slope is −2 + 1/ln(x/x0) ≈ −1.8 at
  x/x0 ~ 100. Variance diverges only logarithmically; CLT holds with the
  effective variance growing like log² of the sample size. **Case 1 is in this
  class** (only this mechanism is present).
- **Re-entrant chord of a convex body (new, dominant in cases 2 and 3).** A
  particle enters the fuel sphere with μ = |n·ω| and leaves through the same
  sphere after a chord l = 2aμ. Composing the piece rules for the entry crossing
  and the exit collision piece:

      dl_exit = (V_1 + V_2 − (n_1 + n_2)·g) / μ,      dl_exit/l = A / (2 a μ²),

  with A = V_1 + V_2 − (n_1 + n_2)·g and g the value before entry. The collision
  lands in the chord with probability ≈ 2aΣμ. Entries have density ∝ μ dμ
  (cosine law). Hence

      P(|dl/l| > x) ∝ ∫_0^{μ*} μ · μ dμ ∝ μ*³,   μ* = (|A|/(2a x))^{1/2}
                    ∝ E|A|^{3/2} · x^{−3/2}.

  So **α = 3/2 exactly, with no log**: the variance is power-law infinite, while
  the mean is finite. This holds for any parameter, as long as A ≠ 0. A ≠ 0 holds
  for V on the sphere, and for any nonzero g (from any other moving surface,
  including the planes and the cylinder). It holds for both variants and for any
  smooth convex body with curvature along the flight (spheres, cylinders: fuel
  pins, TRISO), but not for planar facets.

### 4.2 Measurements (top records: largest |D_own| per generation; `analysis_log.txt`)

**Mechanism check.**
- In cases 2 and 3, the records whose collision piece starts on the surface it ends
  on (re-entrant chords) all satisfy l²/(2s) = 3.500 cm = a. That is the chord
  geometry exactly.
- They are 52–90% of the top-100 |D_own| in every channel, and 33–78% of the top
  5000.
- In case 1 there are none.

**Hill estimates of α** (± bootstrap over seeds):

| | chords k = 100 / 300 / 1000 | rest k = 100 / 300 / 1000 | all, k = 1000 | CCDF of dl/l: pure / log-corrected |
|---|---|---|---|---|
| case 1 (a) R | (none) | 1.90(21) / 1.92(11) / 1.84(6) | 1.84(6) | 1.90(2) / 2.00(5) |
| case 1 (a) tx | (none) | 1.91(21) / 1.81(10) / 1.89(7) | 1.89(6) | 1.90(2) / 2.01(6) |
| case 2 (a) a | 1.46(16) / 1.59(10) / 1.51(5) | 2.11(27) / 1.75(9) / 1.87(6) | 1.75(5) | 1.77(1) / 1.86(2) |
| case 2 (a) b | 1.55(18) / 1.62(9) / 1.52(4) | 1.98(21) / 1.95(12) / 1.91(5) | 1.60(5) | 1.65(2) / 1.75(3) |
| case 2 (c) a | 1.66(20) / 1.53(9) / 1.55(5) | 2.09(19) / 1.89(11) / 1.78(8) | 1.63(5) | 1.67(2) / 1.76(2) |
| case 2 (c) b | 1.66(21) / 1.53(9) / 1.54(5) | 2.09(20) / 1.89(12) / 1.78(8) | 1.63(6) | 1.66(2) / 1.76(4) |
| case 3 (a), 7 params | 1.48–1.63 (k = 1000, ±0.04–0.06) | 1.82–1.94 | 1.59–1.70 | 1.65–1.78 / 1.75–1.87 |
| case 3 (c), 7 params | 1.51–1.68 (k = 1000) | 1.83–1.89 | 1.58–1.73 | 1.66–1.69 / 1.75–1.79 |

- **Case 1:** α = 2 with the log correction. The log-corrected CCDF fit gives
  2.00 ± 0.05, and the Hill values sit at the 1.80–1.87 expected for
  x^-2 ln(x/x0) at u/x0 = 1e2–1e3.
- **Cases 2 and 3:**
  - the chord class gives α = 1.5 within ±0.05 at k = 1000 in every channel, with
    no drift from k = 300 to k = 1000;
  - the rest class gives about 1.85, i.e. the x^-2 ln x class;
  - the pooled estimate and the CCDF fits are mixtures. They drift down toward
    1.5 as the threshold rises (see `fig_case2_tail.png`).
  - The (c) case 1 R channel is bounded (Hill > 13) because of the degeneracy.
- **|g| at fission sites:** tail index 1.8–2.1 for both variants in cases 2 and 3,
  i.e. α ≈ 2. g is not where the infinite variance comes from, but it feeds A in
  the chord term.

### 4.3 How single events become 9σ seeds (worst-seed traces)
- **Inheritance.** An event in generation o enters slot o and is inherited by all
  descendants. It is scored once per descendant fission at lag L, i.e. on the
  diagonal n − L = o of the per-generation tallies S1[n, L].
- **Case 3 seed 25 (z = +9.1 for zlo, −8.8 for (c) tz).** One chord event in
  generation 167 (D_own = 1.4e4, l = 0.028 cm, s = 1.1e-4 cm). Its diagonal
  carries 95–100% of the seed's excess. Descendant multiplicity at L = 0..12:
  1 1 1 1 2 2 2 3 6 7 7 6 4.
- **Case 3 seed 68 ((a) a, tx, z = ±4.5).** A lineage grows to ×25–29 by L = 12.
- **Case 2 seed 54 (z = +9.3 (a) a, +8.2 (c) a, −7.9/−8.1 for b).** One chord
  event in inactive generation 95. The variance over 100 seeds drops to 12% of
  its value ((a) a) or 32% ((c) a) when this seed is removed.
- **(a) vs (c) on the same event:**
  - (c) puts the event on **one** diagonal (share 1.02 of the excess).
  - (a) spreads it over **8** diagonals with alternating signs. The chord term
    A_o = −(n_1 + n_2)·g_o hits every passive slot o with a nonzero g_o at once.
  - More generally, in (a) the slot scored at lag L collects passive terms from
    L + 1 generations of the lineage; in (c) it collects terms from one
    generation. That is why (a)'s variance grows faster with L.
- **Fission-score shift is not the measure.** The largest single fission-event
  score shifts its seed by 0.4–2 seed SDs in most channels, at most 5.2 (case 3
  (a) tz). The 9σ comes from multiplicity
  and, for (a), from the passive slots.

### 4.4 Per-history variance σ²_hist = Var_seed × H

With α = 3/2 these are sample-dependent numbers, not constants. Where they are
dominated by one seed this is marked †.

| | L=1 | L=3 | L=5 | L=8 | L=12 | L=12 without the worst seed / with |
|---|---|---|---|---|---|---|
| case 1 (a) R | 5.04 | 1.81 | 2.94 | 3.1 | 3.44 | 0.95 |
| case 1 (c) R | 0.23 | 0.46 | 0.74 | 1.18 | 1.62 | 0.93 |
| case 1 LR ρ_f | 4.09 | 8.04 | 13.0 | 20.8 | 28.5 | 0.93 |
| case 2 (a) a | 57.9 | 79.4 | 362† | 577† | 959† | 0.12 |
| case 2 (c) a | 22.6 | 35.4 | 156† | 383† | 323† | 0.32 |
| case 2 (a) b | 7.8 | 15.3 | 96.7† | 203† | 107† | 0.37 |
| case 2 (c) b | 8.8 | 14.9 | 62.2† | 153† | 128† | 0.33 |
| case 3 (a) a | 44.1 | 41.4 | 48.4 | 156 | 51.6 | 0.80 |
| case 3 (c) a | 9.2 | 24.7 | 28.6 | 69 | 42.7 | 0.58 |
| case 3 (a) tx | 17.5 | 16.3 | 20.7 | 23.4 | 25.2 | 0.81 |
| case 3 (c) tx | 2.8 | 3.1 | 4.0 | 6.4 | 8.1 | 0.87 |
| case 3 (a) ty | 20.2 | 15.0 | 14.3 | 15.5 | 19.9 | 0.86 |
| case 3 (c) ty | 2.5 | 2.1 | 2.9 | 5.9 | 5.1 | 0.94 |
| case 3 (a) c | 3.3 | 4.3 | 5.2 | 8.0 | 8.6 | 0.93 |
| case 3 (c) c | 3.7 | 5.3 | 6.5 | 8.4 | 8.0 | 0.92 |
| case 3 (a) zhi | 0.33 | 0.45 | 0.54 | 0.82 | 1.11 | 0.89 |
| case 3 (c) zhi | 0.37 | 0.57 | 0.69 | 0.90 | 1.19 | 0.91 |
| case 3 (a)/(c) zlo | 1.2/1.2 | 1.7/1.8 | 3.2/3.3 | 20†/21† | 11†/11† | 0.16/0.15 |

Where no single seed dominates, (c) has 2–4× smaller per-history variance at
L = 12 (up to 8× at small L) for translations and for the case 1 radius, about 1.2× for case 3 a, and the same
variance for parameters the fuel does not touch (c, zhi, zlo). It is never
larger beyond noise.

## 5. Error bars: standard vs robust

**What α = 3/2 implies.**
- The sum of n histories is in the domain of attraction of a stable law with index
  3/2, not of a Gaussian. The error falls as n^(1/α − 1) = n^(−1/3), not
  n^(−1/2), once the rare chord events dominate the core.
- FOM = 1/(error² T) then falls like T^(−1/3) instead of staying constant.
- **Seed-based SEs are unreliable.** They usually underestimate, because most
  samples hold no large event, and occasionally one seed carries most of the
  variance. Batch means and the bootstrap only re-weight the same sample
  variance. The bootstrap of the mean is inconsistent for α < 2 (Athreya 1987).

**What is in the table.** Two robust estimates are shown next to the standard SE,
batch-means SE and bootstrap SE:
- **SE core** = IQR/1.349/√R: the error the run would have without rare events
  (optimistic).
- **SE stable**: m-out-of-R subsampling (m = 25 seeds, without replacement, IQR
  of the subsample means), rescaled to R with the stable rate (m/R)^(1/3). This
  is consistent in the stable domain. When the variance is finite it is
  conservative by (R/m)^(1/6) = 1.26.
- **β**: the measured exponent of IQR(mean of m seeds) ~ m^-β, with
  finite-population correction.

| channel (L = 12) | SE std | batch | bootstrap | core | stable | β |
|---|---|---|---|---|---|---|
| case 1 (a) R | 4.15e-5 | 3.74e-5 | 4.03e-5 | 4.65e-5 | 5.09e-5 | 0.53 |
| case 1 (c) R | 2.84e-5 | 2.07e-5 | 2.88e-5 | 2.76e-5 | 3.57e-5 | 0.50 |
| case 2 (a) a | 69.3e-5 | 66.1e-5 | 69.4e-5 | 12.3e-5 | 109e-5 | −0.01 |
| case 2 (c) a | 40.2e-5 | 36.5e-5 | 39.8e-5 | 10.8e-5 | 57.0e-5 | 0.06 |
| case 2 (a) b | 23.1e-5 | 22.0e-5 | 23.1e-5 | 5.95e-5 | 32.0e-5 | 0.03 |
| case 2 (c) b | 25.4e-5 | 24.3e-5 | 25.0e-5 | 5.68e-5 | 35.6e-5 | 0.00 |
| case 2 LR ρ_f | 12.8e-5 | 10.5e-5 | 12.8e-5 | 12.4e-5 | 15.7e-5 | 0.49 |
| case 3 (a) a | 16.1e-5 | 18.3e-5 | 16.4e-5 | 14.3e-5 | 20.3e-5 | 0.47 |
| case 3 (c) a | 14.6e-5 | 15.7e-5 | 14.6e-5 | 10.3e-5 | 15.9e-5 | 0.43 |
| case 3 (a) tx | 11.2e-5 | 14.0e-5 | 10.9e-5 | 10.1e-5 | 14.4e-5 | 0.46 |
| case 3 (c) tx | 6.36e-5 | 5.16e-5 | 6.38e-5 | 5.01e-5 | 6.82e-5 | 0.49 |
| case 3 (a) zlo | 7.38e-5 | 7.79e-5 | 7.21e-5 | 2.72e-5 | 7.52e-5 | 0.30 |
| case 3 (c) zlo | 7.42e-5 | 7.92e-5 | 7.39e-5 | 2.23e-5 | 10.0e-5 | 0.20 |

- **Case 1:** std ≈ core ≈ batch ≈ bootstrap and β ≈ 0.5. The variance is
  effectively finite at this sample size (log-divergence only) and the seed error
  bars are fine.
- **Case 2:** std/core = 3.7–5.7 and β ≈ 0. The error is set by one event, so the
  standard SE is a sample of one large-deviation event, not an error bar.
  Batch means and bootstrap reproduce the standard SE (batch/std = 0.95) and add
  nothing.
- **Case 3:** a, tx, ty, c and zhi are still core-dominated (std/core ≈ 1.0–1.4,
  β ≈ 0.45–0.5). zlo and (c) tz are single-event dominated (std/core 2.4–3.3,
  β 0.2–0.35).
- **Calibration against the exact references** (z of the L = 12 mean, using SE
  std / core / stable):
  - case 2 (a) a +1.4 / +7.8 / +0.9;
  - case 2 (c) a +1.6 / +5.8 / +1.1;
  - case 2 (a) b −1.7 / −6.8 / −1.3;
  - case 2 (c) b −1.8 / −8.0 / −1.3;
  - case 1 channels within 2.4σ for all three.

  The core SE alone is badly anti-conservative in case 2. The stable SE is the
  usable one there. Under α = 3/2 the standard SE happens to be adequate here only
  because the dominant event is inside the sample.

## 6. Equal-CPU efficiency vs CRN central differences (L = 12, dk/dθ)

FOM = 1/(SE² × total CPU). The derivative run is charged in full to every
parameter, although one run gives all of them; FD needs one pair per parameter.
The FD SE uses std / batch / core / core, since the k-only FD has no heavy tail.
For the derivative runs the four values use SE std / batch / core / stable.

@@EFF@@

## 7. The dl/l term and the l distribution

**dl/l in the own-generation slot, fission histories.**

| | (a) P(\|J\| > \|D_own − J\|) | (a) E[J²]/E[D²] | (c) P(…) | (c) E[J²]/E[D²] |
|---|---|---|---|---|
| case 1 R / tx | 0.71 / 0.66 | 1.41 / 1.28 | 0.70 / 0.53 | 4.68 / 2.64 |
| case 2 a / b | 0.70 / 0.09 | 1.04 / 1.14 | 0.69 / 0.68 | 1.07 / 1.15 |
| case 3 a, tx, ty, tz | 0.67–0.69 | 1.01–1.06 | 0.66–0.69 | 1.09–1.20 |
| case 3 c / zlo / zhi | 0.08 / 0.03 / 0.03 | 1.23 / 1.04 / 1.28 | 0.65–0.67 | 1.06–1.26 |

- dl/l dominates the own-generation D in about 2/3 of fission histories wherever
  the slot has a nonzero g or V_n.
- (a) is lower only where the own generation rarely touches the moving surface
  (b, c, z planes). The passive slots carry the rest.
- E[J²] ≥ E[D²]: the remaining terms (−Σ dl) are anti-correlated with J and
  partly cancel it. In the (c) case 1 degenerate channel, LR's terms cancel J
  almost entirely (4.7).

**Collision pieces: quantiles of l and of s = l|n·ω| [cm].**

| quantile | 1e-6 | 1e-5 | 1e-4 | 1e-3 | 0.01 | 0.1 | 0.5 | 0.9 |
|---|---|---|---|---|---|---|---|---|
| case 1 l | 0.0056 | 0.018 | 0.056 | 0.18 | 0.56 | 1.4 | 3.5 | 5.6 |
| case 1 s | 0.0022 | 0.0071 | 0.022 | 0.089 | 0.28 | 1.1 | 2.8 | 5.6 |
| case 2 l | 0.0028 | 0.0089 | 0.028 | 0.089 | 0.28 | 1.1 | 2.8 | 5.6 |
| case 2 s | 0.0007 | 0.0022 | 0.0089 | 0.035 | 0.14 | 0.71 | 2.2 | 4.5 |
| case 3 l | 0.0028 | 0.0089 | 0.028 | 0.089 | 0.35 | 1.1 | 3.5 | 5.6 |
| case 3 s | 0.0007 | 0.0028 | 0.0089 | 0.035 | 0.14 | 0.71 | 2.2 | 4.5 |

- The histograms are 0.1-decade bins, so the values are bin centres.
- The low-quantile tails of cases 2 and 3 are 2–3× shorter than in case 1. These
  are the re-entrant chords.
- The top-100 |D_own| events have:
  - median l = 0.02–0.1 cm and median s = 2e-4 to 1.5e-3 cm in cases 2 and 3;
  - median l = 0.006 cm and s = 6e-4 cm in case 1 (collisions within 1e-3 cm of
    the surface).

## 8. g statistics across generations

| | (a) RMS g_o at age 0 / 1 / 4 / 8 / 12 | (a) decay/gen | (a) max | (c) stationary RMS | (c) max | (c) generations to 1% of stationary RMS |
|---|---|---|---|---|---|---|
| case 1 R | 0.59 / 0.51 / 0.32 / 0.19 / 0.12 | 0.877 | 11.8 | 0.68 | 1 (= x/R) | 17 |
| case 1 tx | 0.37 / 0.32 / 0.20 / 0.12 / 0.07 | 0.878 | 7.6 | 1 (= e) | 1 | 28 |
| case 2 a | 1.01 / 1.05 / 0.72 / 0.47 / 0.37 | 0.909 | 127 | 1.15 | 93 | 7 |
| case 2 b | 0.29 / 0.24 / 0.18 / 0.13 / 0.08 | 0.910 | 48 | 0.57 | 59 | 12 |
| case 3 a | 1.02 / 1.07 / 0.68 / 0.40 / 0.23 | 0.872 | 161 | 1.08 | 112 | 0 |
| case 3 tx | 0.72 / 0.73 / 0.44 / 0.25 / 0.15 | 0.870 | 132 | 1.05 | 62 | 16 |
| case 3 c | 0.27 / 0.22 / 0.15 / 0.08 / 0.05 | 0.862 | 53 | 0.50 | 71 | — (pooled RMS noise > 1%) |
| case 3 zlo / zhi | 0.10 / 0.09 / 0.05 / 0.03 / 0.01 | 0.85 / 0.89 | 18 / 11 | 0.18 / 0.17 | 23 / 18 | 79 / 29 |

- (a) passive slots decay geometrically, at 0.85–0.91 per generation.
- The RMS above is pooled over seeds and dominated by rare inherited large g. A
  robust measure, the median over seeds of each seed's RMS per generation,
  changes by ≤ 1% over the 200 active generations in every channel of both
  variants: (c) case 2 a 0.928 → 0.931, b 0.379 → 0.382; case 3 a 0.917 → 0.918.
  So **(c)'s g is stationary and does not grow**. The pooled RMS (for example
  case 2 a 1.13 → 1.22) moves with individual outliers.
- RMS g is O(1) cm/cm. Max |g| is about 100. The |g| tail index is ≈ 2 for both
  variants.

## 9. Variance reduction for the dl/l term (proposals, not implemented), by priority

**Key lemma (why nothing piece-local can work).** A reparameterization of a single
flight piece must map [0, l] onto [0, l + dl]. Its Jacobian therefore averages
dl/l = (V_n − n·g)/s over the piece, so the 1/s singularity (and the 1/μ² chord
singularity, which is the same identity applied twice) remains whatever happens
inside the piece. The singularity is removable only if the displacement used at
the surface is consistent with V_n there, n·g = V_n, for every particle. That
requires a deterministic, smooth displacement field rather than a per-particle
g, which is what (c)'s g approximates and what it becomes exactly in the
degenerate cases.

### P1 (top priority): smooth warp field, i.e. a domain-mapping reparameterization

- **Construction.** Choose a deterministic C¹ field 𝒱(x) with compact support
  such that:
  - 𝒱·n = V_n on each moving surface;
  - 𝒱·n = 0 on the fixed surfaces;
  - 𝒱 is affine in a band of width δ around each moving surface.

  Examples: a constant e on and around a translated sphere or pin; a dilation
  about the centre or axis for a radius; e_z for a plane offset. Each is
  multiplied by a smooth cutoff in the distance to the surface.
- **Particle state.** The walk keeps only D slots, as LR does; there is no g.
- **Unbiasedness.**
  - Let Φ_θ(y) = y + θ𝒱(y), a diffeomorphism that maps every reference surface
    onto its perturbed position. The material regions are then θ-independent in
    reference coordinates y.
  - The pulled-back transport-fission operator Φ_θ* H_θ (Φ_θ*)^-1 is similar to
    H_θ, so it has the same eigenvalue k(θ).
  - In y-coordinates the point-to-point kernel of a free flight y → y′ is

        Σ(y′) exp(−τ_θ) / (4π |Φ_θ(y′) − Φ_θ(y)|²) · det ∇Φ_θ(y′),

    with τ_θ the optical length along the physical segment. At θ = 0 this is
    the analog kernel. The walk is sampled at θ = 0, so the derivative estimator
    is the ordinary score function (LR) of this kernel. It is unbiased for the
    frozen-L functional of the pulled-back operator, and converges to d ln k/dθ
    as L → ∞, exactly as LR does.
  - That finite-L target differs from (a)'s. How fast it converges in L is
    untested.
- **Per-segment score terms**, with Δ𝒱 = 𝒱(y′) − 𝒱(y) and ω the flight direction:
  - attenuation: −Σ ω·Δ𝒱 in a homogeneous segment;
  - geometric spreading: −2 ω·Δ𝒱/|y′ − y|;
  - collision-site Jacobian: ∇·𝒱(y′);
  - at an interface crossing at y_c: (Σ_1 − Σ_2) n·[𝒱(y_c) − 𝒱_lin(y_c)]/(n·ω),
    where 𝒱_lin is the linear interpolation of 𝒱 between the endpoints.
  The first three are bounded by |∇𝒱|. The crossing term vanishes if 𝒱 is affine
  along the segment. It is O(l²|∇²𝒱|)/μ otherwise.
- **Effect on α.**
  - **No dl/l term at all**: the near-surface collision (x^-2 ln x) and the
    re-entrant chord (x^-3/2) both disappear. For a chord, l = 2aμ, so even a
    non-affine 𝒱 gives a crossing term O(μ²)/μ = O(μ), which is bounded.
  - What remains is the grazing crossing term for long segments (l = O(1)) that
    leave the affine band. It has α = 2 (density μ dμ, term ∝ 1/μ) and an
    amplitude ∝ |∇²𝒱|, controlled by δ.
  - The variance should then behave like LR's: σ²_hist ≈ 20–30 at L = 12 here,
    growing linearly with L.
- **Cost and risk.**
  - One evaluation of 𝒱 per event and parameter, plus ∇·𝒱 at collisions, plus
    𝒱 at each interface crossing.
  - The hard part is constructing 𝒱 in general CSG near surface intersections,
    where 𝒱 must be tangent to the fixed surfaces. It is easy for pin-type
    perturbations (rigid translation, radius) and for planes perpendicular to the
    other surfaces.

### P2 (second): forced collision on short re-entrant chords, with weighted fission sites

- **Method.** On entering a convex region with a short chord (Σl < ε, or μ < μ_c),
  follow both branches instead of sampling collide or not:
  - collided, with weight w_c = 1 − e^{−Σl} and the collision point from the
    truncated exponential on [0, l];
  - uncollided, with weight w_u = e^{−Σl}.
- **Unbiasedness.**
  - This is an exact expectation over one Bernoulli variable, with both futures
    sampled: d/dθ E = Σ_b [w_b′ E_b + w_b E_b′].
  - Each branch carries D + d ln w_b + d ln p_b. For the collided branch,
    d ln w_c + d ln p(u | c) = dl/l − Σ u dl, the same as the analog increment,
    but now multiplied by w_c ≈ Σl = 2aΣμ.
- **Effect on α.** The weighted chord contribution is w_c·A/(2aμ²) ≈ ΣA/μ, and
  every entry contributes, not only the collided ones. So P(|wD| > x) ∝ x^-2:
  **α 3/2 → 2**, with a log-divergent variance, i.e. the case 1 class. It does
  nothing for near-surface collisions.
- **Caveat.** It only works if the collided branch's fission sites keep the
  fractional weight w_c, in the IFP numerator and denominator and in the fission
  bank. With OpenMC-style unit-weight sites (number of sites by stochastic
  rounding of w ν Σf/Σt/k), a branch of weight w_c yields a site with probability
  ∝ w_c that carries the full 1/μ² D, and the tail returns unchanged. The same
  objection kills plain particle splitting at chord entry.
- **Cost.** Cheap: P(μ < μ_c) ∝ μ_c², so this is a rare extra branch. It needs
  weighted fission banking, which is the invasive part.

### P3 (low): baseline (control variate) on the score increments

- **Method.**
  - Each score-function increment ΔD (a dl/l term, a crossing term) is later
    multiplied by its future F, the lag-L descendant score of the lineage.
  - Tally additionally −b(x)·ΔD at the event, where b(x) is a deterministic or
    previous-run estimate of E[F | state before the event], e.g. an IFP
    importance map.
- **Unbiased** because E[ΔD | state before the event] = 0, since a score function
  has zero conditional mean. Then E[b(x) ΔD] = E[b · E[ΔD | x]] = 0.
- **Effect.** The event's contribution becomes ΔD·(F − b). It reduces the
  amplitude where F is predictable from x. It does not change α, because the
  singular ΔD is still multiplied by an O(1) random F − b, which is mostly −b
  since lineages die out.
- **Cost.** Expected gain is modest. It needs an importance map.

### Rejected
- **Rao-Blackwellizing D over the collision position:** D multiplies the whole
  future, so the expectation cannot be taken without that future.
- **Clipping or trimming large D:** biased. The bias is of the order of the
  clipped tail mass, and it is precisely the part the estimator needs.
- **Piece-local reparameterization fixes:** ruled out by the lemma.
- **"(c) with g reset to 𝒱 at birth":** needs passive defect slots, i.e. (a)'s
  bookkeeping, and keeps the chord term for g − 𝒱 ≠ 0.
- **Median-of-means as the tally estimator:** gives better-behaved confidence
  intervals, but with only moments < 3/2 finite the rate is still n^(−1/3). It
  is a reporting tool, not a fix.

## 10. Pilot episode
In the pilot, the 4 seeds of case 2 were about 3.5σ high on ρ_r (13σ apparent
with 4-seed error bars). This was not reproduced in 52 other seeds. Seed 1 was
bit-reproducible, and there was no seed-to-seed correlation. The long run (seeds
1–100, an independent set) gives ρ_r at L = 12 of 0.04561(6) against 0.04573
(−1.9σ), so the pilot excursion was a fluctuation of the LR channel. It is
unrelated to the heavy tails of §4: LR has σ² ≈ 8 and no dl/l term.
See `pilot/pilot_report.md`.

## 11. Files and reproduction
- **Environment.**

      python -m venv notes/toy_eigen/.venv
      notes/toy_eigen/.venv/bin/pip install -r notes/toy_eigen/requirements.txt

  Python 3.13.9, numba 0.68.0, numpy 2.3.5, scipy 1.18.1, matplotlib 3.11.2.
- **Run** (from `notes/toy_eigen/3d/`):
  - `../.venv/bin/python sn_sphere.py` writes `sn_sphere_results.json` (S_N exact).
  - `../.venv/bin/python sn_sphere_frozen.py` writes `sn_sphere_frozen_results.json`
    (finite-L S_N).
  - `./run_long.sh` runs `run3d.py sens CASE 100 100000 100 200 1` for cases 1–3,
    then `run3d.py fd CASE 100 100000 100 200 0.1 PARAMS 1001`. Output goes to
    `runs/caseC/`. About 3.5 h on 20 threads.
  - `../.venv/bin/python analyze3d.py runs --figs > runs/analysis_log.txt` writes
    the tables, `runs/analysis.json` and `runs/fig_case{1,2,3}_{L,tail}.png`.
- **Raw data per seed** (`runs/caseC/sens_seedNNN.npz`):
  - per-generation S1, S2 and Smax for every channel and L;
  - k per generation;
  - jstat, gstat, ghist, and the dl/l, l and s histograms;
  - top records: the 6 largest |D_own| per generation, variant and parameter, as
    [D_own, J, largest single dl/l term, its piece length l, s = l|n·ω|,
    distance from the piece start to that surface, surface id, generation,
    particle index].
- **Kernel and model:** `mc3d.py` (kernel), `geom3d.py` (model), `run3d.py`
  (drivers).
- **Pilot:** `pilot/`, with W = 9. It is bit-identical to the current kernel for
  L ≤ 8.
