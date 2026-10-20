import sys
from pathlib import Path
p = Path(sys.argv[1]); s = p.read_text()
R = [
("| finite-L S_N a | 0.13370 | 0.10760 | 0.09725 | 0.08971 | 0.08752 | 0.08683 | 0.08656 |",
 "| finite-L S_N a | 0.13370 | 0.10760 | 0.09725 | 0.08971 | 0.08708 | 0.08669 | 0.08656 |"),
("| finite-L S_N b | 0.00983 | 0.00945 | 0.00925 | 0.00909 | 0.00903 | 0.00902 | 0.00901 |",
 "| finite-L S_N b | 0.00983 | 0.00945 | 0.00925 | 0.00909 | 0.00902 | 0.00901 | 0.00901 |"),
("- All rows are within 2σ of their targets at every L (max |z| = 1.9).",
 "- (a) and LR are within 2σ of their finite-L targets at every L (max |z| = 1.9).\n"
 "- (c) is within 1.8σ of the exact value for L ≥ 5. Its deviation at L ≤ 3 is its\n"
 "  truncation bias (a: −1.8% at L = 1, −0.5% at L = 3)."),
("| (a) zhi | 0.000284(14) | 0.000283(14) | 0.000289(15) | 0.000287(17) | 0.000277(20) | 0.000287(22) | 0.000276(24) |",
 "| (a) zhi | 0.000284(13) | 0.000283(14) | 0.000289(15) | 0.000287(16) | 0.000277(20) | 0.000288(22) | 0.000277(24) |"),
("| (c) zhi | 0.000217(13) | 0.000252(15) | 0.000273(17) | 0.000273(18) | 0.000280(21) | 0.000288(22) | 0.000275(24) |",
 "| (c) zhi | 0.000217(14) | 0.000252(15) | 0.000273(17) | 0.000273(19) | 0.000280(21) | 0.000288(22) | 0.000275(24) |"),
("- **Symmetry.** ty, tz and k_zlo + k_zhi are consistent with 0 at every L for both\n"
 "  variants (max |z| 2.2, (a) ty at L = 5; all others ≤ 1.7).",
 "- **Symmetry.** ty, tz and k_zlo + k_zhi are consistent with 0 at every L for both\n"
 "  variants: max |z| = 2.2 ((a) ty at L = 5) and 2.1 ((a) k_zlo + k_zhi at L = 5),\n"
 "  out of 56 checks, all others ≤ 1.7."),
("- **Case 3 seed 68 ((a) a, tx, z = ±4.5).** A lineage grows to ×25 by L = 12.",
 "- **Case 3 seed 68 ((a) a, tx, z = ±4.5).** A lineage grows to ×25–29 by L = 12."),
("- **Fission-score shift is not the measure.** The largest single fission-event\n"
 "  score shifts its seed by at most 1–3.6 seed SDs.",
 "- **Fission-score shift is not the measure.** The largest single fission-event\n"
 "  score shifts its seed by 0.4–2 seed SDs in most channels, at most 5.2 (case 3\n"
 "  (a) tz)."),
("Where no single seed dominates, (c) has 2–7× smaller per-history variance for\n"
 "translations and for the case 1 radius,",
 "Where no single seed dominates, (c) has 2–4× smaller per-history variance at\n"
 "L = 12 (up to 8× at small L) for translations and for the case 1 radius,"),
("| case 3 c | 0.27 / 0.22 / 0.15 / 0.08 / 0.05 | 0.862 | 53 | 0.50 | 71 | (< 1% from the start) |",
 "| case 3 c | 0.27 / 0.22 / 0.15 / 0.08 / 0.05 | 0.862 | 53 | 0.50 | 71 | — (pooled RMS noise > 1%) |"),
("- (a) passive slots decay geometrically, at 0.85–0.91 per generation. That is\n"
 "  slower than the radial dominance ratio would suggest (0.58–0.63), because g is\n"
 "  transported, not projected.\n"
 "- (c)'s single g reaches a stationary distribution within the 100 inactive\n"
 "  generations in every channel. It does not grow over the 200 active\n"
 "  generations.",
 "- (a) passive slots decay geometrically, at 0.85–0.91 per generation.\n"
 "- The RMS above is pooled over seeds and dominated by rare inherited large g. A\n"
 "  robust measure, the median over seeds of each seed's RMS per generation,\n"
 "  changes by ≤ 1% over the 200 active generations in every channel of both\n"
 "  variants: (c) case 2 a 0.928 → 0.931, b 0.379 → 0.382; case 3 a 0.917 → 0.918.\n"
 "  So **(c)'s g is stationary and does not grow**. The pooled RMS (for example\n"
 "  case 2 a 1.13 → 1.22) moves with individual outliers."),
]
for a, b in R:
    assert s.count(a) == 1, a[:70]
    s = s.replace(a, b)
i0 = s.index("### P3 (low)"); i1 = s.index("### Rejected")
s = s[:i0] + """### P3 (low): baseline (control variate) on the score increments

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

""" + s[i1:]
p.write_text(s); print("ok")
