# Step 0b pilot report (3D toy)

The full printouts are in pilot_analysis.txt (seeds 1–4, all cases) and
check2_analysis.txt (case 2, seeds 101–140). Pilot runs used N = 1e5,
100 inactive + 100 active generations, and CRN FD at h = 0.1 cm.

## References
- Spherical S_N (sn_sphere.py): diamond difference in r and μ, starting
  direction, dense inner solve. On the bare sphere it matches the exact
  slab-odd-mode equivalence to 2e-7 after O(1/N²) extrapolation. The mesh is
  converged to 1e-8.
  - Case 1: d ln k/dR = 0.094291.
  - Case 2: d ln k/da = 0.086492, d ln k/db = 0.009010,
    d ln k/dρ_f = 0.306553, d ln k/dρ_r = 0.045723.
  - The scaling identity holds in S_N to 1e-9.
- Finite-L S_N (sn_sphere_frozen.py) for radial parameters: the exact targets
  of (a) and LR at each L. Radial dominance ratio: 0.628 (case 1), 0.577 (case 2).
- S_N FD truncation at h = 0.1: +0.03% (R), +0.04% (a), +0.30% (b).

## Correctness
- Case 1 (4 seeds): (a) R and LR ρ_f match finite-L S_N at L = 0..8 (|z| ≤ 1.3).
  (c) R matches the exact value at L ≥ 5 (z = −0.5 at L = 8). tx is 0 by
  symmetry for (a) (|z| ≤ 1.1); for (c) it is identically 0.
- Case 2 (40 seeds): every channel is within 1.7σ of its finite-L S_N target,
  or of the exact value for (c) at L = 8. The CRN FD agrees.
- Case 3 (4 seeds): a, tx, c, zhi all agree with CRN FD within 1.5σ (the FD is
  noisy). ty, tz and k_zlo + k_zhi are consistent with 0.
- Pilot excursion. The 4 pilot seeds of case 2 were jointly about 3.5σ high on
  ρ_r (13σ apparent with 4-seed error bars). It was not reproduced by 52 other
  seeds. Seed 1 is bit-reproducible, and there is no correlation between
  consecutive seeds. Conclusion: a statistical fluctuation, and 4-seed error
  bars are useless here.

## Predicted degeneracies (confirmed)
(c) is pathwise LR along the global dilation, and for a rigid translation of
the whole system. Hence case 1 (c) R = LR/R exactly, case 1 (c) tx ≡ 0, and the
case 3 scaling identity holds for (c) to 1e-10. The case 3 identity therefore
tests only (a).

## Early (a) vs (c) signals (case 2, 40 seeds, θ = a)
- L-convergence: (a) +12.8% at L = 3, +4.1% at L = 5, +0.6% at L = 8.
  (c) −0.6% at L = 3, −0.1% at L = 5, +0.1% at L = 8.
- σ²_hist at L = 8: (a) 70, (c) 29. For b they are equal (about 11).
- Seed tails: (a) kurtosis +8.9 with max |z| 4.6; (c) kurtosis +0.4 with max |z| 2.9.
- g: (a) passive decays about 0.86 per generation. (c) has a stationary RMS of
  about 1 (translations) and a max of about 100. Heavy-tailed, but no growth.
- dl/l: it dominates the own-generation D in about 2/3 of fission histories for
  every parameter whose surface the fuel touches, and E[J²]/E[D²] ≈ 1.1–1.9.
  The CCDF of |dl/l| falls about as x^−1.8; (c)'s amplitude is about 4× lower than (a)'s.

## Timings (20 threads, N = 1e5, 200 generations)
| | derivative run, all channels | FD pair, k-only |
|---|---|---|
| case 1 | 14 s | 4 s |
| case 2 | 9 s | 5 s |
| case 3 | 25 s | 5 s |
