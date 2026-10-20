"""Fill the markers of results.md from the final analysis (runs/analysis_log.txt)."""
from pathlib import Path

here = Path(__file__).resolve().parent.parent
body = (here / "draft/body.md").read_text()
res = (here / "results.md").read_text()

def sub(text, old, new, count=1):
    assert text.count(old) == count, (old, text.count(old))
    return text.replace(old, new)

# --- §2 FD comparisons -------------------------------------------------------
body = sub(body, "- vs CRN FD at L = 12: (a) @@C1FDZA@@σ, (c) @@C1FDZC@@σ.",
           "- vs CRN FD (+0.09443 ± 0.00013) at L = 12: (a) +0.1σ, (c) −1.6σ.")
body = sub(body, "@@C3FDa@@", "0.08829(13)")
body = sub(body, "@@C3FDtx@@", "−0.00426(13)")
body = sub(body, "@@C3FDc@@", "0.00697(14)")
body = sub(body, "@@C3FDzhi@@", "0.00028(13)")
body = sub(body, "@@C3FDTEXT@@", """\
- **vs CRN FD at L = 12** (z with SE std; core and stable give the same picture,
  |Δz| ≤ 0.3):
  - (a): a −2.1, tx +0.1, c −0.5, zhi 0.0;
  - (c): a −0.8, tx −1.2, c −0.4, zhi 0.0.

  No channel is off by more than 2.1σ. The open question from the pilot, (a) tx
  −0.00423(11) vs (c) tx −0.00443(6), is not resolved by FD: the difference is
  1.6σ combined, and FD (±0.00013) sits between them. Its error is twice that of
  (c).
- The FD error is about 1.3e-4 for every parameter (CRN noise of the k
  difference divided by 2h), so FD does not resolve zhi (0.00028 ± 0.00013) at
  all, while both derivative variants give it to ±9%.
- Truncation bias against FD is visible at small L: (c) a is −1.9% at L = 1 and
  −0.3% at L = 5; (a) a is +53% at L = 1 and +3.4% at L = 5.""")
body = sub(body, "@@C3FDdka@@", "FD 0.09047(13)")
body = sub(body, "@@C3FDdktx@@", "FD −0.00436(13)")

# --- §5 robust table and text from the final analysis --------------------------
old_tab = body[body.index("| channel (L = 12) | SE std"):body.index("- **Case 1:** std ≈ core")]
new_tab = """\
| channel (L = 12) | SE std | batch | bootstrap | core | stable | β |
|---|---|---|---|---|---|---|
| case 1 (a) R | 4.15e-5 | 3.74e-5 | 4.03e-5 | 4.65e-5 | 5.09e-5 | 0.53 |
| case 1 (c) R | 2.84e-5 | 2.07e-5 | 2.88e-5 | 2.76e-5 | 3.57e-5 | 0.50 |
| case 2 (a) a | 69.3e-5 | 66.1e-5 | 70.5e-5 | 12.3e-5 | 102e-5 | −0.02 |
| case 2 (c) a | 40.2e-5 | 36.5e-5 | 38.6e-5 | 10.8e-5 | 54.9e-5 | 0.05 |
| case 2 (a) b | 23.1e-5 | 22.0e-5 | 23.0e-5 | 5.95e-5 | 33.2e-5 | 0.02 |
| case 2 (c) b | 25.4e-5 | 24.3e-5 | 25.5e-5 | 5.68e-5 | 36.8e-5 | −0.01 |
| case 2 LR ρ_f | 12.8e-5 | 10.5e-5 | 12.5e-5 | 12.4e-5 | 15.6e-5 | 0.50 |
| case 3 (a) a | 16.1e-5 | 18.3e-5 | 15.9e-5 | 14.3e-5 | 20.3e-5 | 0.47 |
| case 3 (c) a | 14.6e-5 | 15.7e-5 | 14.6e-5 | 10.5e-5 | 18.8e-5 | 0.39 |
| case 3 (a) tx | 11.2e-5 | 14.0e-5 | 11.3e-5 | 10.1e-5 | 14.1e-5 | 0.45 |
| case 3 (c) tx | 6.36e-5 | 5.16e-5 | 6.32e-5 | 5.31e-5 | 7.76e-5 | 0.46 |
| case 3 (a) zlo | 7.38e-5 | 7.79e-5 | 7.49e-5 | 2.72e-5 | 9.81e-5 | 0.25 |
| case 3 (c) zlo | 7.42e-5 | 7.92e-5 | 7.47e-5 | 2.25e-5 | 11.0e-5 | 0.17 |

"""
body = body.replace(old_tab, new_tab)
body = sub(body, """\
- **Case 3:** a, tx, ty, c and zhi are still core-dominated (std/core ≈ 1.0–1.4,
  β ≈ 0.45–0.5). zlo and (c) tz are single-event dominated (std/core 2.4–3.3,
  β 0.2–0.35).""", """\
- **Case 3:** a, tx, ty, c and zhi are still core-dominated (std/core ≈ 1.0–1.4,
  β ≈ 0.39–0.52). zlo and (c) tz are single-event dominated (std/core 2.2–3.3,
  β 0.17–0.32).""")
body = sub(body, "  - case 2 (a) b −1.7 / −6.8 / −1.3;", "  - case 2 (a) b −1.7 / −6.8 / −1.2;")
body = sub(body, "  - case 2 (c) b −1.8 / −8.0 / −1.3;", "  - case 2 (c) b −1.8 / −8.0 / −1.2;")

# --- §6 efficiency -------------------------------------------------------------
body = sub(body, "@@EFF@@", """\
Ratios of FOM, each cell **std / batch / core / stable** (`runs/analysis_log.txt`
has the absolute FOMs). A ratio > 1 means the first estimator is more efficient.

| case: parameter | (a)/FD | (c)/FD | (c)/(a) |
|---|---|---|---|
| 1: R | 3.4 / 6.7 / 2.6 / 2.1 | 7.0 / 21 / 7.4 / 4.7 | 2.1 / 3.1 / 2.8 / 2.3 |
| 2: a | 0.016 / 0.015 / 0.40 / 0.0071 | 0.048 / 0.050 / 0.51 / 0.020 | 3.0 / 3.3 / 1.3 / 2.8 |
| 2: b | 0.097 / 0.16 / 1.3 / 0.046 | 0.081 / 0.13 / 1.4 / 0.037 | 0.83 / 0.83 / 1.1 / 0.81 |
| 3: a | 0.081 / 0.052 / 0.088 / 0.045 | 0.097 / 0.071 / 0.18 / 0.053 | 1.2 / 1.4 / 2.0 / 1.2 |
| 3: tx | 0.18 / 0.073 / 0.19 / 0.093 | 0.55 / 0.54 / 0.67 / 0.29 | 3.1 / 7.4 / 3.6 / 3.1 |
| 3: c | 0.60 / 0.24 / 0.56 / 0.37 | 0.64 / 0.25 / 0.62 / 0.40 | 1.07 / 1.07 / 1.1 / 1.09 |
| 3: zhi | 4.1 / 13 / 6.0 / 3.2 | 3.8 / 15 / 5.7 / 3.1 | 0.93 / 1.12 / 0.94 / 0.99 |

Cost per seed: derivative run 18.5 / 18.2 / 55.5 s (cases 1 / 2 / 3); FD pair
6.7 / 7.5 / 7.2 s per parameter.

**How to read it.**
- **Per parameter vs amortized.** The table charges the whole derivative run to
  each parameter. One run gives all of them, so the amortized ratio vs FD is the
  table value × P, with P = 2 / 2 / 7 geometric parameters (plus the densities).
  Amortized, with the stable SE, (c)/FD is 9.3 (case 1 R), 0.04 and 0.07
  (case 2 a, b), and 0.37 / 2.0 / 2.8 / 22 (case 3 a, tx, c, zhi).
- **The derivative cost is an upper bound.** It includes both variants, LR and all
  channels. A production run of one variant would be cheaper (not measured).
- **FD depends on h.** h = 0.1 cm keeps the truncation ≤ 0.3%. FD's SE is about
  the same for every parameter, so it is relatively worst for small derivatives
  (zhi), where the derivative estimators win by 3–15×.
- **Batch means** use only 10 groups, so the batch column is noisy (it moves the
  FD FOM by up to 3×). Base conclusions on std and stable.
- **Heavy tail.** Where α = 3/2 dominates (case 2, case 3 zlo), the derivative FOM
  falls like T^(−1/3) with run length while FD's stays constant. The case 2 gaps
  vs FD widen with longer runs; the table is at H = 2e9 active histories per
  parameter.

**Does the conclusion change between standard and robust error bars?**
- **(c) vs (a): no.** (c)/(a) is ≥ 1.07 in every channel and under every estimate,
  except case 2 b (0.81–0.83 with std and stable, 1.1 with core) and case 3 zhi
  (0.93–1.12, parity). The two exceptions are surfaces the fuel does not touch,
  where both variants carry the same g-transport term. The ordering is the same
  with std and with stable; the margins change by ≤ 10% (case 1, case 3) and
  ≤ 7% (case 2).
- **vs FD: it depends on the case.**
  - Case 1 (no chords): no change. The derivative wins, (c) by 4.7–21×.
  - Case 3 (core-dominated): no change in direction; the stable SE lowers the
    ratios by up to 2×. Per parameter FD wins for a, tx and c, (c) wins for zhi.
    Amortized over 7 parameters, (c) wins for tx, c and zhi and loses for a.
  - Case 2 (one chord event dominates): **yes.** With the core SE alone, (c) is at
    parity with FD (0.5–1.4). With std or stable, FD is 20–50× more efficient per
    parameter (12–27× amortized). The core SE is 6–8σ anti-conservative here
    (§5), so the stable numbers are the honest ones: in re-entrant geometry the
    current estimators lose to CRN FD, and lose more the longer the run.""")

# --- summary ------------------------------------------------------------------
summary = """\
## Recommendation: variant (c), and a variance-reduction fix for the chord term before production

**Choose (c).** Evidence, all at matched histories on the same walk:

1. **Both are unbiased.**
   - (a) matches its finite-L S_N target at every L in cases 1 and 2 (|z| ≤ 1.9
     beyond the mesh effect at L ≤ 2).
   - (c) matches exact S_N and CRN FD once L ≥ 5.
   - Case 3: no channel is off FD by more than 2.1σ. All symmetry checks pass
     (max |z| 2.2 of 56), and (a) passes the scaling identity at L = 12.
2. **(c) converges much faster in L.**
   - Case 1 R: (c) is −1.1% at L = 3 and −0.08% at L = 12. (a) is +16% at L = 3
     and still +0.17% at L = 12.
   - Cases 2–3: (c) is within 0.5% at L = 3; (a) is +12% (a) at L = 3.
3. **(c) has equal or lower variance.**
   - σ²_hist is 2–4× smaller for translations and the case 1 radius, about 1.2×
     for the case 3 radius, and equal for surfaces the fuel does not touch.
   - Equal-CPU (c)/(a) is 1.1–3.6 except case 2 b (0.81) and zhi (parity). This
     is the same under the standard and the robust error bars.
4. **(c) is structurally simpler and better behaved.**
   - A rare event enters one tally slot in (c). In (a) it enters about 8 slots,
     with alternating signs, through the passive g_o.
   - (c) keeps one g per parameter instead of W slots of g.
   - (c)'s g is stationary: ≤ 1% change of the robust RMS over 200 generations.
5. **(c)'s known degeneracy is harmless.** It reduces exactly to LR for a global
   dilation and to 0 for a rigid translation of the whole system (§3), so the
   scaling identity cannot test it. The non-degenerate channels (case 2, case 3
   individual parameters) are where it was tested, and it passes there.

**Both variants share a heavy tail that must be fixed before production (§4, §9).**
- Re-entrant chords of the convex fuel sphere give dl/l = A/(2aμ²), so
  P(|dl/l| > x) ∝ x^(−3/2). This is the tail probability (CCDF), not the density.
  α = 3/2 exactly: finite mean, infinite variance.
  - Measured: chord-class Hill α = 1.5 ± 0.05 in every channel of both variants
    in cases 2 and 3.
  - Case 1 has no chords and is in the α = 2, x^-2 ln x class: Hill 1.84–1.92,
    log-corrected CCDF 2.00 ± 0.05.
  - The pilot's "x^-1.8" was a local CCDF slope of a mixture of the two classes.
- Consequences:
  - The error falls as n^(−1/3) once rare events dominate.
  - Seed SEs, batch means and the bootstrap are unreliable there; one case 2 seed
    holds 68–88% of the variance.
  - The derivative FOM decays like T^(−1/3).
  - In case 2, CRN central differences are 20–50× more efficient per parameter by
    the robust estimate.
- Robust error bars: the subsampling "stable" SE (§5) gives |z| ≤ 1.3 against the
  exact references, where the IQR "core" SE gives 6–8σ.
- **Variance-reduction priorities** (not implemented):
  - **P1:** a smooth warp field (domain mapping). It removes the dl/l term
    entirely and should give α = 2 or better.
  - **P2:** forced collision on short chords. α 3/2 → 2, but only with
    fractional-weight fission sites.
  - **P3:** a baseline (control variate). It is unbiased but leaves α unchanged.

---"""
res = sub(res, "@@SUMMARY@@", summary)
for k, v in (("@@FDT1@@", "6.7"), ("@@FDT2@@", "7.5"), ("@@FDT3@@", "7.2")):
    res = sub(res, k, v)
res = sub(res, "@@BODY@@", body.rstrip() + "\n")
assert "@@" not in res, [l for l in res.splitlines() if "@@" in l]
(here / "results.md").write_text(res)
print("ok", len(res.splitlines()), "lines")
