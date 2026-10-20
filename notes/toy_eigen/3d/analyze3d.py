"""
Analysis of the 3D toy runs.

    python analyze3d.py RUNDIR [--figs]

RUNDIR contains case1/, case2/, case3/ with sens_seed*.npz and fd_*_seed*.npz.
Prints markdown tables and writes RUNDIR/analysis.json (and figures with --figs).

Error bars:
- standard: SD over independent seeds / sqrt(R);
- batch: 10 groups of R/10 seeds, SD of the group means / sqrt(10);
- bootstrap: resampling seeds (SE, and the 95% percentile interval as half-width / 1.96);
- core: IQR of the seed values / 1.349 / sqrt(R), the error without rare events;
- stable: subsampling without replacement, m = 25 seeds, IQR of the subsample means
  / 1.349, rescaled to R by (m/R)^(1 - 1/alpha) with alpha = 3/2 (consistent for
  the mean in the domain of attraction of a stable law; conservative by
  (R/m)^(1/6) = 1.26 when the variance is finite).
The batch and bootstrap SEs are re-weightings of the same sample variance and
are not consistent when alpha < 2.
Tails: |dl/l| per collision piece (histograms), and the per-history own-generation D
of fission histories (top records), Hill estimator with a log-corrected model,
split into re-entrant chord pieces (piece starts on the surface it ends on) and
the rest.
"""
import json
import sys
from pathlib import Path

import numpy as np

import mc3d
from geom3d import CASES

HERE = Path(__file__).parent
LS = [1, 2, 3, 5, 8, 10, 12]
LREF = 12                       # window used for efficiency and tail summaries
NGROUP = 10
RNG = np.random.default_rng(12345)
EDGES = np.linspace(mc3d.HLO, mc3d.HHI, mc3d.NB + 1)


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def load_sn():
    out = {}
    f = HERE / "sn_sphere_results.json"
    if f.exists():
        r = json.load(open(f))
        out = {c: r[f"case{c}_reference"] for c in (1, 2)}
    f = HERE / "sn_sphere_frozen_results.json"
    if f.exists():
        r = json.load(open(f))
        for c in (1, 2):
            out.setdefault(c, {})["frozen"] = r[f"case{c}_dr0.005_S128"]
    return out


def load_sens(d):
    fs = sorted(d.glob("sens_seed*.npz"))
    if not fs:
        return None
    keys_sum = ("jstat", "ghist", "hl", "hs")
    out = {k: [] for k in ("est", "est_q", "est_h", "k", "secs", "hj", "top", "smax", "den")}
    out["files"] = fs
    for k in keys_sum:
        out[k] = 0
    out["gstat"] = 0
    for f in fs:
        z = np.load(f)
        n_in = int(z["n_inactive"])
        S1, den = z["S1"], z["den"]
        na = S1.shape[0]
        out["est"].append(S1.sum(0) / den.sum())
        out["est_h"].append(S1[: na // 2].sum(0) / den[: na // 2].sum())
        out["est_q"].append(S1[: na // 4].sum(0) / den[: na // 4].sum())
        out["k"].append(z["kgen"][n_in:].mean())
        out["secs"].append(float(z["seconds"]))
        out["hj"].append(z["hj"])
        out["top"].append(z["top"])
        out["smax"].append(z["Smax"])
        out["den"].append(den.sum())
        for k in keys_sum:
            out[k] = out[k] + z[k]
        out["gstat"] = out["gstat"] + z["gstat"]
        out["names"] = [str(x) for x in z["names"]]
        out["N"], out["n_act"], out["n_in"] = int(z["N"]), na, n_in
        out["W"] = S1.shape[2]
    for k in ("est", "est_q", "est_h", "k", "secs", "hj", "top", "smax", "den"):
        out[k] = np.array(out[k])
    out["R"] = len(fs)
    out["H"] = out["N"] * out["n_act"]
    out["gstat"] = out["gstat"] / out["R"]
    return out


def load_fd(d):
    res = {}
    for f in sorted(d.glob("fd_*_seed*.npz")):
        p = f.name.split("_")[1]
        z = np.load(f)
        res.setdefault(p, []).append(((z["kp"] - z["km"]) / (2 * z["h"]), float(z["seconds"]), float(z["h"])))
    out = {}
    for p, v in res.items():
        x = np.array([a for a, _, _ in v])
        out[p] = dict(x=x, secs=np.mean([s for _, s, _ in v]), h=v[0][2], R=len(x))
    return out


# ---------------------------------------------------------------------------
# error estimates
# ---------------------------------------------------------------------------
def errors(x):
    """x: per-seed values (R,). Returns mean, standard SE, batch SE, bootstrap SE, bootstrap pct SE."""
    R = x.size
    se = x.std(ddof=1) / np.sqrt(R)
    g = R // NGROUP
    if g >= 2:
        gm = x[: g * NGROUP].reshape(NGROUP, g).mean(1)
        seb = gm.std(ddof=1) / np.sqrt(NGROUP)
    else:
        seb = np.nan
    bs = x[RNG.integers(0, R, size=(2000, R))].mean(1)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return x.mean(), se, seb, bs.std(ddof=1), (hi - lo) / (2 * 1.96)


def iqr(x):
    q1, q3 = np.percentile(x, [25, 75])
    return q3 - q1


def robust_errors(x, m=25, alpha=1.5, nsub=4000):
    """Core SE (IQR of the seeds), stable-rate subsampling SE, and the subsampling
    exponent beta of IQR(mean of m seeds) ~ m^-beta for m = 1..25 (finite-population
    corrected; 1/2 for finite variance, 1 - 1/alpha in the stable regime, ~0 when a
    single seed dominates)."""
    R = x.size
    core = iqr(x) / 1.349 / np.sqrt(R)
    ms = np.array([1, 2, 4, 8, 16, m])
    q = []
    for mm in ms:
        idx = np.argsort(RNG.random((nsub, R)), axis=1)[:, :mm]
        q.append(iqr(x[idx].mean(1)) / np.sqrt(1 - (mm - 1) / (R - 1)))
    q = np.array(q)
    beta = -np.polyfit(np.log(ms), np.log(np.maximum(q, 1e-300)), 1)[0] if q.min() > 0 else np.nan
    stable = q[-1] * (m / R) ** (1 - 1 / alpha) / 1.349
    return core, stable, beta


# ---------------------------------------------------------------------------
# tails
# ---------------------------------------------------------------------------
def hill(vals, ks):
    """Hill estimates of the tail index for the k largest |vals|."""
    a = np.sort(np.abs(vals[vals != 0]))[::-1]
    out = []
    for k in ks:
        if k + 1 >= a.size:
            out.append((np.nan, np.nan))
            continue
        out.append((k / np.sum(np.log(a[:k] / a[k])), a[k]))
    return out


def hill_model_fit(alphas, us):
    """Fit 1/alpha_hat(u) = 1/alpha + 1/(alpha^2 ln(u/x0)), the Hill expectation for
    P(|X| > x) ~ x^-alpha ln(x/x0). Grid search over x0, linear in the rest."""
    best = None
    inv = 1 / np.asarray(alphas)
    for lx0 in np.linspace(-3, np.log(min(us)) - 0.05, 400):
        t = 1 / (np.log(us) - lx0)
        A = np.column_stack([np.ones_like(t), t])
        c, res, *_ = np.linalg.lstsq(A, inv, rcond=None)
        r = np.sum((A @ c - inv) ** 2)
        if best is None or r < best[0]:
            best = (r, c, lx0)
    _, c, lx0 = best
    return 1 / c[0], np.exp(lx0)


def ccdf_fit(h, x_lo, min_count=20):
    """Fits of the CCDF of a histogram above x_lo (bins with >= min_count counts above).
    Returns (alpha_pure, alpha_log, x0_log, x_hi, (x, P)):
    alpha_pure from log P = c - alpha log x;
    alpha_log, x0_log from P = C x^-alpha ln(x/x0) (x0 < x_lo, grid search)."""
    tail = np.cumsum(h[::-1])[::-1]               # counts above each lower edge
    x = 10 ** EDGES[:-1]
    P = tail / max(h.sum(), 1)
    sel = (x >= x_lo) & (tail >= min_count)
    if sel.sum() < 4:
        return np.nan, np.nan, np.nan, np.nan, (x, P)
    lx, lp = np.log(x[sel]), np.log(P[sel])
    a_pure = -np.polyfit(lx, lp, 1)[0]
    best = None
    for lx0 in np.linspace(np.log(x_lo) - 8, np.log(x_lo) - 0.2, 300):
        y = lp - np.log(lx - lx0)
        c = np.polyfit(lx, y, 1)
        res = np.sum((np.polyval(c, lx) - y) ** 2)
        if best is None or res < best[0]:
            best = (res, -c[0], np.exp(lx0))
    return a_pure, best[1], best[2], x[sel].max(), (x, P)


# ---------------------------------------------------------------------------
# report per case
# ---------------------------------------------------------------------------
def case_report(case, d, sn, figs):
    r = load_sens(d)
    if r is None:
        return {}
    fd = load_fd(d)
    names, P, R, H, W = r["names"], len(r["names"]), r["R"], r["H"], r["W"]
    Ls = [L for L in LS if L < W]
    est, k = r["est"], r["k"]
    kbar = k.mean()
    dk = est * k[:, None, None]
    m, se = est.mean(0), est.std(0, ddof=1) / np.sqrt(R)
    rep = dict(R=R, H=H, secs=float(r["secs"].mean()), k=(kbar, k.std(ddof=1) / np.sqrt(R)))
    print(f"\n## Case {case}: {CASES[case]}")
    print(f"seeds R = {R}, N = {r['N']}, {r['n_in']} inactive + {r['n_act']} active "
          f"(H = {H:.2g} active histories per seed), {r['secs'].mean():.1f} s/seed for all "
          f"channels; k = {kbar:.6f} ± {rep['k'][1]:.6f}"
          + (f" (S_N {sn[case]['k']:.6f})" if case in sn and "k" in sn[case] else ""))

    chans = [(f"(a) {nm}", p, nm) for p, nm in enumerate(names)] + \
            [(f"(c) {nm}", P + p, nm) for p, nm in enumerate(names)] + \
            [("LR rho_f", 2 * P, "rho_f"), ("LR rho_r", 2 * P + 1, "rho_r")]
    if case == 1:
        chans = chans[:-1]

    # references per channel: list of (label, value, se, kind)
    def refs_for(lab, nm):
        out = []
        if case in sn and ("dlnk_d" + nm) in sn[case]:
            out.append(("S_N exact", sn[case]["dlnk_d" + nm], 0.0, "exact"))
        fz = sn.get(case, {}).get("frozen", {})
        key = "dlnQ_d" + nm
        if key in fz and (lab.startswith("(a)") or lab.startswith("LR")):
            out.append(("S_N finite-L", np.array(fz[key]), 0.0, "finite"))
        # (c) along the global dilation is LR (case 1: R only)
        if case == 1 and lab == "(c) R" and "dlnQ_drho_f" in fz:
            out.append(("S_N finite-L (= LR/R)", np.array(fz["dlnQ_drho_f"]) / CASES[1]["R"], 0.0, "finite"))
        if nm in fd:
            x = fd[nm]["x"]
            out.append(("CRN FD", x.mean() / kbar, x.std(ddof=1) / np.sqrt(x.size) / kbar, "exact"))
        if nm in ("ty", "tz") or (case == 1 and nm == "tx"):
            out.append(("symmetry", 0.0, 0.0, "exact"))
        return out

    print("\n### d ln k/dθ [1/cm] (densities: dimensionless) vs L, ±1σ standard (seeds)")
    print("z: E = vs exact S_N, F = vs CRN FD, f = vs finite-L S_N target, s = vs symmetry 0\n")
    print("| channel | " + " | ".join(f"L={L}" for L in Ls) + " | references |")
    print("|---|" + "---|" * (len(Ls) + 1))
    tab = {}
    for lab, ch, nm in chans:
        refs = refs_for(lab, nm)
        cells = []
        for L in Ls:
            zs = []
            for rn, rv, rs, kind in refs:
                v = rv[L] if kind == "finite" else rv
                tag = {"S_N exact": "E", "CRN FD": "F", "symmetry": "s"}.get(rn, "f")
                zs.append(f"{tag}{(m[ch, L] - v) / np.hypot(se[ch, L], rs):+.1f}")
            cells.append(f"{m[ch, L]:+.5f} ± {se[ch, L]:.5f}" + (" (" + ",".join(zs) + ")" if zs else ""))
        rs_ = "; ".join(f"{rn} {(rv if kind != 'finite' else rv[LREF]):+.5f}"
                        + (f" ± {rs:.5f}" if rs else "") + (" (L=12)" if kind == "finite" else "")
                        for rn, rv, rs, kind in refs)
        print(f"| {lab} | " + " | ".join(cells) + f" | {rs_} |")
        tab[lab] = dict(mean=m[ch].tolist(), se=se[ch].tolist(), dk_mean=dk[:, ch].mean(0).tolist(),
                        dk_se=(dk[:, ch].std(0, ddof=1) / np.sqrt(R)).tolist(),
                        refs=[(rn, (rv.tolist() if kind == "finite" else rv), rs, kind) for rn, rv, rs, kind in refs])
    rep["table"] = tab

    # dk/dtheta at L = 8, 12 (with references in dk units)
    print("\n### dk/dθ at L = 8 and 12 (standard SE)\n")
    print("| channel | L=8 | L=12 | reference |")
    print("|---|---|---|---|")
    for lab, ch, nm in chans:
        x8, x12 = dk[:, ch, 8], dk[:, ch, 12] if W > 12 else dk[:, ch, -1]
        ref_s = []
        if case in sn and ("dk_d" + nm) in sn[case]:
            ref_s.append(f"S_N {sn[case]['dk_d' + nm]:+.5f}")
        if nm in fd:
            x = fd[nm]["x"]
            ref_s.append(f"FD {x.mean():+.5f} ± {x.std(ddof=1) / np.sqrt(x.size):.5f}")
        print(f"| {lab} | {x8.mean():+.5f} ± {x8.std(ddof=1) / np.sqrt(R):.5f} | "
              f"{x12.mean():+.5f} ± {x12.std(ddof=1) / np.sqrt(R):.5f} | {'; '.join(ref_s)} |")

    # ---- error estimates: standard vs robust ----
    print(f"\n### error estimates at L = {LREF}: standard vs robust (all ×1e5)\n")
    print("core = IQR/1.349/√R (no rare events); stable = m-out-of-R subsampling (m = 25) rescaled "
          "with (m/R)^(1/3) for α = 3/2; β = subsampling exponent of IQR(mean of m seeds) ~ m^-β "
          "(1/2: finite variance; 1/3: stable α = 3/2; ≈ 0: a single seed dominates).\n")
    print("| channel | mean | SE std | SE batch (10 groups) | SE bootstrap | bootstrap 95% half-width/1.96 "
          "| SE core | SE stable | β | std/core | stable/std |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    errs = {}
    for lab, ch, nm in chans:
        x = est[:, ch, LREF]
        mu, s1, sb, sbs, spc = errors(x)
        if s1 > 0:
            core, stab, beta = robust_errors(x)
        else:
            core = stab = beta = np.nan
        errs[lab] = dict(mean=mu, std=s1, batch=sb, boot=sbs, boot_pct=spc, core=core, stable=stab, beta=beta)
        print(f"| {lab} | {mu:+.5f} | {s1 * 1e5:.2f} | {sb * 1e5:.2f} | {sbs * 1e5:.2f} | {spc * 1e5:.2f} | "
              f"{core * 1e5:.2f} | {stab * 1e5:.2f} | {beta:.2f} | {s1 / core:.2f} | {stab / s1:.2f} |")
    rep["errors_L12"] = errs

    # z-scores of the L = 12 means against the references, with each error estimate
    zr = []
    for lab, ch, nm in chans:
        for rn, rv, rs, kind in refs_for(lab, nm):
            if rn == "symmetry" or (kind == "exact" and rn == "S_N exact" and lab.startswith("(a)")):
                continue      # (a) targets the finite-L value; symmetry zeros are in the main table
            v = rv[LREF] if kind == "finite" else rv
            e = errs[lab]
            if not e["std"] > 0:
                continue
            zr.append((lab, rn, *[(e["mean"] - v) / np.hypot(e[k_], rs) for k_ in ("std", "core", "stable")]))
    if zr:
        print(f"\nz of the L = {LREF} mean vs its reference with SE std / core / stable:")
        print("; ".join(f"{lab} vs {rn}: {a:+.1f} / {b:+.1f} / {c:+.1f}" for lab, rn, a, b, c in zr))
    rep["z_L12"] = zr

    # ---- per-history variance, growth with sample size, seed outliers ----
    var = est.var(0, ddof=1) * H
    var_h = r["est_h"].var(0, ddof=1) * H / 2
    var_q = r["est_q"].var(0, ddof=1) * H / 4
    print("\n### per-history variance σ²_hist = Var_seed × H, growth with run length, seed outliers\n")
    print("| channel | " + " | ".join(f"L={L}" for L in Ls) +
          f" | L={LREF}: quarter / half / full run | max \\|z\\| over seeds | var without worst seed / var |")
    print("|---|" + "---|" * (len(Ls) + 3))
    tails = {}
    for lab, ch, nm in chans:
        x = est[:, ch, LREF]
        sd = x.std(ddof=1)
        z = (x - x.mean()) / sd if sd > 0 else np.zeros(R)
        iw = np.argmax(np.abs(z))
        loo = np.delete(x, iw).var(ddof=1) / x.var(ddof=1) if sd > 0 else np.nan
        tails[lab] = dict(var=var[ch].tolist(), var_q=var_q[ch, LREF], var_h=var_h[ch, LREF],
                          zmax=float(np.abs(z).max()), loo=loo)
        print(f"| {lab} | " + " | ".join(f"{var[ch, L]:.3g}" for L in Ls) +
              f" | {var_q[ch, LREF]:.3g} / {var_h[ch, LREF]:.3g} / {var[ch, LREF]:.3g} | "
              f"{np.abs(z).max():.2f} | {loo:.2f} |")
    rep["variance"] = tails

    # largest single fission-event score, as a shift of its run estimate in run-sigmas
    print(f"\nlargest single fission-event score at L = {LREF}, as a shift of its seed's estimate "
          "(units of the seed-to-seed SD):")
    print(", ".join(f"{lab} {(np.abs(r['smax'][:, :, ch, LREF]).max(1) / r['den']).max() / est[:, ch, LREF].std(ddof=1):.2f}"
                    for lab, ch, nm in chans if est[:, ch, LREF].std(ddof=1) > 0))
    print("(one fission event only; an event's D is inherited by all its descendants, see the trace below)")

    # ---- trace of the worst seed: which origin slots carry its excess ----
    # An event in generation o enters slot o (and, for (a), every passive slot with
    # nonzero g); it is scored at lag L in generation o + L, i.e. on a diagonal of
    # S1[n, L]. Excess = S1 - median over generations.
    print(f"\n### worst seed per channel (|z| > 4 at L = {LREF}): origin of its excess\n")
    print("| channel | seed | z | share of the excess at L=12 from the largest diagonal | diagonals with > 10% "
          "| origin generation o (active index) | E(o, L)/E(o, 0) for L = 0..12 ((c): descendants carrying the event's D; (a): also its passive follow-up terms) |")
    print("|---|---|---|---|---|---|---|")
    trace = {}
    for lab, ch, nm in chans:
        x = est[:, ch, LREF]
        sd = x.std(ddof=1)
        if not sd > 0:
            continue
        z = (x - x.mean()) / sd
        w = int(np.argmax(np.abs(z)))
        if abs(z[w]) < 4:
            continue
        zz = np.load(r["files"][w])
        S1 = zz["S1"][:, ch, :]
        E = S1 - np.median(S1, 0)
        na = S1.shape[0]
        diag = {}                    # origin o -> excess at L = 12
        for n in range(na):
            diag[n - LREF] = E[n, LREF]
        tot = sum(diag.values())
        o_best = max(diag, key=lambda o: abs(diag[o]))
        share = diag[o_best] / tot
        nbig = sum(abs(v) > 0.1 * abs(tot) for v in diag.values())
        mult = ""
        if o_best >= 0 and abs(E[o_best, 0]) > 0:
            mult = " ".join(f"{E[o_best + L, L] / E[o_best, 0]:.1f}" for L in range(W) if o_best + L < na)
        trace[lab] = dict(seed=r["files"][w].name, z=float(z[w]), share=float(share), ndiag=int(nbig), o=int(o_best))
        print(f"| {lab} | {r['files'][w].stem[-3:]} | {z[w]:+.1f} | {share:.2f} | {nbig} | {o_best} | {mult or '(origin before the active cycles)'} |")
    rep["trace"] = trace

    # ---- equal-CPU efficiency, standard and robust ----
    if fd:
        T = r["secs"].mean()
        print(f"\n### equal-CPU efficiency at L = {LREF}: FOM = 1/(SE² · total CPU), dk/dθ, per parameter. The "
              "derivative run is charged in full to every parameter (one run yields all of them). "
              "Four error estimates: std / batch / core / stable (FD: std / batch / core / core, "
              "the k-only FD has no heavy tail).\n")
        print("| parameter | FOM (a) | FOM (c) | FOM FD | (a)/FD | (c)/FD | (c)/(a) |")
        print("|---|---|---|---|---|---|---|")
        eff = {}

        def four(x, secs, heavy=True):
            _, s1, sb, _, _ = errors(x)
            core, stab, _ = robust_errors(x)
            ses = (s1, sb, core, stab if heavy else core)
            return np.array([1 / (s_ ** 2 * x.size * secs) for s_ in ses])

        for p, nm in enumerate(names):
            if nm not in fd:
                continue
            fa, fc = four(dk[:, p, LREF], T), four(dk[:, P + p, LREF], T)
            ff = four(fd[nm]["x"], fd[nm]["secs"], heavy=False)
            eff[nm] = dict(a=fa.tolist(), c=fc.tolist(), fd=ff.tolist())
            f3 = lambda v: " / ".join(f"{u:.3g}" for u in v)
            print(f"| {nm} | {f3(fa)} | {f3(fc)} | {f3(ff)} | {f3(fa / ff)} | {f3(fc / ff)} | {f3(fc / fa)} |")
        rep["fom"] = eff
        print(f"(T derivative run = {T:.1f} s/seed; FD pair = " +
              ", ".join(f"{nm} {fd[nm]['secs']:.1f} s" for nm in fd) + f"; FD h = {fd[next(iter(fd))]['h']} cm)")

    # ---- dl/l dominance ----
    js = r["jstat"]
    print("\n### dl/l term J in the own-generation slot, fission histories\n")
    print("| param | (a) P(\\|J\\| > \\|D_own − J\\|) | (a) E[J²]/E[D²] | (c) P(\\|J\\| > \\|D_own − J\\|) | (c) E[J²]/E[D²] |")
    print("|---|---|---|---|---|")
    jrep = {}
    for p, nm in enumerate(names):
        row = []
        for v in (0, 1):
            c_, dom, j2, d2, o2 = js[v, p]
            row.append((dom / c_, j2 / d2 if d2 > 0 else np.nan))
        jrep[nm] = row
        print(f"| {nm} | {row[0][0]:.3f} | {row[0][1]:.2f} | {row[1][0]:.3f} | {row[1][1]:.2f} |")
    rep["dl_over_l"] = jrep

    # ---- l and s at collisions ----
    ctr = 10 ** (0.5 * (EDGES[1:] + EDGES[:-1]))

    def quant(h, qs):
        c = np.cumsum(h) / h.sum()
        return [ctr[np.searchsorted(c, q)] for q in qs]

    qs = [1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 0.5, 0.9]
    print("\n### collision pieces: quantiles of the piece length l and the normal distance s = l |n·ω| [cm]\n")
    print("| quantile | " + " | ".join(f"{q:g}" for q in qs) + " |")
    print("|---|" + "---|" * len(qs))
    print("| l | " + " | ".join(f"{x:.2g}" for x in quant(r["hl"], qs)) + " |")
    print("| s | " + " | ".join(f"{x:.2g}" for x in quant(r["hs"], qs)) + " |")
    rep["l_quant"], rep["s_quant"] = quant(r["hl"], qs), quant(r["hs"], qs)

    # ---- tail index ----
    print("\n### tail index α of P(|X| > x)\n")
    print("Per collision piece, X = dl/l (histograms of all collision pieces): log-log CCDF fit above "
          "x = 10, (i) pure power law, (ii) P ∝ x^-α ln(x/x0) with x0 free. "
          "Per fission history, X = D_own (top records): Hill estimator from the k largest, "
          "u = the k-th largest value. ± = bootstrap over seeds (200 resamples).\n")
    print("| variant: param | CCDF pure α | CCDF log-corr. α (x0) | fit range x | Hill k=100 | Hill k=300 | Hill k=1000 (u) | Hill k=3000 (u) |")
    print("|---|---|---|---|---|---|---|---|")
    ks = [100, 300, 1000, 3000]
    trep = {}
    for v, vn in ((0, "a"), (1, "c")):
        for p, nm in enumerate(names):
            h = r["hj"][:, v, p].sum(0)
            vals = r["top"][:, :, v, p, :, 0].ravel()
            if h.sum() == 0 or np.abs(vals).max() < 1.0:
                print(f"| ({vn}) {nm} | degenerate (no dl/l tail) | | | | | | |")
                continue
            ap, al, x0, xh, _ = ccdf_fit(h, 10.0)
            bs = []
            for _ in range(200):
                idx = RNG.integers(0, R, R)
                bs.append(ccdf_fit(r["hj"][idx, v, p].sum(0), 10.0)[:2])
            bs = np.array(bs)
            hh = hill(vals, ks)
            hbs = np.array([[a for a, _ in hill(r["top"][RNG.integers(0, R, R), :, v, p, :, 0].ravel(), ks)]
                            for _ in range(200)])
            hs_ = np.nanstd(hbs, 0)
            trep[f"{vn}:{nm}"] = dict(ccdf_pure=(ap, np.nanstd(bs[:, 0])), ccdf_log=(al, np.nanstd(bs[:, 1]), x0),
                                      x_hi=xh, hill=[(a, s_, u) for (a, u), s_ in zip(hh, hs_)])
            print(f"| ({vn}) {nm} | {ap:.2f} ± {np.nanstd(bs[:, 0]):.2f} | {al:.2f} ± {np.nanstd(bs[:, 1]):.2f} ({x0:.2g}) | "
                  f"10–{xh:.0e} | {hh[0][0]:.2f} ± {hs_[0]:.2f} | {hh[1][0]:.2f} ± {hs_[1]:.2f} | "
                  f"{hh[2][0]:.2f} ± {hs_[2]:.2f} ({hh[2][1]:.3g}) | {hh[3][0]:.2f} ± {hs_[3]:.2f} ({hh[3][1]:.3g}) |")
    rep["tail"] = trep
    print("\nFor reference: if P ∝ x^-2 ln(x/x0), the Hill estimator at threshold u has expectation "
          "2 / (1 + 1/(2 ln(u/x0))), e.g. 1.80 at u/x0 = 100 and 1.87 at u/x0 = 1e3.")

    # what drives the top tail: dl/l share and geometry of the top records
    print("\n### the 100 largest |D_own| per variant/param: what produced them\n")
    print("| variant: param | median \\|D_own\\| | share with \\|largest dl/l term\\| > 0.5 \\|D_own\\| | median s [cm] | median distance to surface [cm] | median l [cm] | surfaces (count) |")
    print("|---|---|---|---|---|---|---|")
    drivers = {}
    for v, vn in ((0, "a"), (1, "c")):
        for p, nm in enumerate(names):
            t = r["top"][:, :, v, p].reshape(-1, mc3d.NTF)
            t = t[t[:, 0] != 0]
            if t.size == 0:
                continue
            t = t[np.argsort(-np.abs(t[:, 0]))][:100]
            share = np.mean(np.abs(t[:, 2]) > 0.5 * np.abs(t[:, 0]))
            sids, cnt = np.unique(t[:, 6].astype(int), return_counts=True)
            drivers[f"{vn}:{nm}"] = dict(med=float(np.median(np.abs(t[:, 0]))), share=share,
                                         s=float(np.median(t[:, 4])), dist=float(np.median(t[:, 5])),
                                         l=float(np.median(t[:, 3])))
            print(f"| ({vn}) {nm} | {np.median(np.abs(t[:, 0])):.3g} | {share:.2f} | {np.median(t[:, 4]):.2g} | "
                  f"{np.median(t[:, 5]):.2g} | {np.median(t[:, 3]):.2g} | "
                  + ", ".join(f"{s_}:{c_}" for s_, c_ in zip(sids, cnt)) + " |")
    rep["tail_drivers"] = drivers

    # chord split: the piece starts on the surface it ends on (re-entrant chord of the
    # convex fuel sphere: l = 2 a mu, s = 2 a mu^2, so l^2 / (2 s) = a). Theory:
    # dl/l ~ 1/(a mu^2), P(collision in chord) ~ mu, entry density ~ mu dmu  ->  alpha = 3/2.
    print("\n### re-entrant chords vs the rest: Hill α per class (top records, ± bootstrap over seeds)\n")
    print("| variant: param | chord share of top 100 / 1000 / 5000 | chords: median l²/(2s) [cm] "
          "| Hill chords k=100 / 300 / 1000 | Hill rest k=100 / 300 / 1000 |")
    print("|---|---|---|---|---|")
    ks3 = [100, 300, 1000]
    crep = {}

    def split_hill(top_vp):
        t = top_vp.reshape(-1, mc3d.NTF)
        t = t[t[:, 0] != 0]
        ch_ = (t[:, 5] < 1e-6) & (t[:, 6] == 0)
        return [a for a, _ in hill(t[ch_, 0], ks3)], [a for a, _ in hill(t[~ch_, 0], ks3)]

    for v, vn in ((0, "a"), (1, "c")):
        for p, nm in enumerate(names):
            tv = r["top"][:, :, v, p]
            t = tv.reshape(-1, mc3d.NTF)
            t = t[t[:, 0] != 0]
            if t.size == 0 or np.abs(t[:, 0]).max() < 1.0:
                continue
            t = t[np.argsort(-np.abs(t[:, 0]))]
            chord = (t[:, 5] < 1e-6) & (t[:, 6] == 0)
            shares = [chord[:K].mean() for K in (100, 1000, 5000)]
            aeff = np.median(t[:1000][chord[:1000], 3] ** 2 / (2 * t[:1000][chord[:1000], 4])) if chord[:1000].any() else np.nan
            hc, hn = split_hill(tv)
            bs = np.array([np.concatenate(split_hill(r["top"][RNG.integers(0, R, R), :, v, p])) for _ in range(200)])
            sd_ = np.nanstd(bs, 0)
            crep[f"{vn}:{nm}"] = dict(shares=shares, aeff=aeff, chord=list(zip(hc, sd_[:3])), rest=list(zip(hn, sd_[3:])))
            fmt = lambda h, s: " / ".join(f"{a:.2f} ± {e:.2f}" for a, e in zip(h, s))
            print(f"| ({vn}) {nm} | " + " / ".join(f"{s_:.2f}" for s_ in shares) + f" | {aeff:.3f} | "
                  f"{fmt(hc, sd_[:3])} | {fmt(hn, sd_[3:])} |")
    rep["chords"] = crep

    # ---- g statistics ----
    gs = r["gstat"]
    n_in = r["n_in"]
    act = slice(n_in, gs.shape[0])
    print("\n### g at fission sites: (a) passive g_o by age n − o; (c) single g\n")
    print("| param | (a) RMS age 0 / 1 / 4 / 8 / 12 | (a) decay/gen | (a) max | (c) stationary RMS | (c) max | (c) RMS gen 0 → 1% of stationary |")
    print("|---|---|---|---|---|---|---|")
    grep = {}
    for p, nm in enumerate(names):
        rms_a = [np.sqrt(gs[act, 0, p, age, 0].sum() / gs[act, 0, p, age, 2].sum()) for age in range(W)]
        mx_a = max(gs[act, 0, p, age, 1].max() for age in range(W))
        rms_c = np.sqrt(gs[:, 1, p, 0, 0] / np.maximum(gs[:, 1, p, 0, 2], 1))
        ages = np.arange(1, W)
        ra = np.array(rms_a[1:])
        fac = np.exp(np.polyfit(ages, np.log(np.maximum(ra, 1e-300)), 1)[0]) if (ra > 0).all() else 0.0
        stat = rms_c[act].mean()
        dev = np.abs(rms_c[:n_in] - stat) / stat
        g1 = int(np.argmax(dev < 0.01)) if (dev < 0.01).any() else -1
        grep[nm] = dict(rms_a=rms_a, decay_a=fac, max_a=mx_a, rms_c=stat, max_c=gs[act, 1, p, 0, 1].max(),
                        gen_1pct=g1, rms_c_by_gen=rms_c.tolist())
        print(f"| {nm} | " + " / ".join(f"{rms_a[a]:.3g}" for a in (0, 1, 4, 8, 12) if a < W) +
              f" | {fac:.3f} | {mx_a:.3g} | {stat:.3g} | {gs[act, 1, p, 0, 1].max():.3g} | {g1} |")
    # tail index of |g| from the histograms (active generations), (c) and (a) age 1
    gt = []
    for p, nm in enumerate(names):
        for v, vn, age in ((0, "a", 1), (1, "c", 0)):
            h = r["ghist"][v, p, age]
            if h.sum() == 0:
                continue
            q99 = quant(h, [0.99])[0]
            gt.append(f"({vn}{', age 1' if v == 0 else ''}) {nm} {ccdf_fit(h, q99, 10)[0]:.2f}")
    print("\ntail index of |g| at fission sites (pure power-law CCDF fit above the 99th percentile; "
          "nan = too few counts or bounded g): " + ", ".join(gt))
    rep["g"] = grep

    # ---- identities ----
    if case == 3:
        c = CASES[3]
        coef = {"a": c["a"], "tx": c["x0"], "c": c["c"], "zlo": c["zlo"], "zhi": c["zhi"]}
        print("\n### scaling identity: a k_a + x0 k_tx + c k_c + zlo k_zlo + zhi k_zhi vs k_ρf + k_ρr\n")
        print("| L | variant | LHS | RHS (LR) | LHS − RHS (z) |")
        print("|---|---|---|---|---|")
        sc = {}
        for vn, off in (("a", 0), ("c", P)):
            for L in Ls:
                lhs = sum(coef[nm] * dk[:, off + p, L] for p, nm in enumerate(names) if nm in coef)
                rhs = dk[:, 2 * P, L] + dk[:, 2 * P + 1, L]
                dd = lhs - rhs
                s_ = dd.std(ddof=1) / np.sqrt(R)
                sc[f"{vn} L={L}"] = (lhs.mean(), rhs.mean(), dd.mean(), s_)
                print(f"| {L} | ({vn}) | {lhs.mean():.5f} ± {lhs.std(ddof=1) / np.sqrt(R):.5f} | "
                      f"{rhs.mean():.5f} ± {rhs.std(ddof=1) / np.sqrt(R):.5f} | {dd.mean():+.2e} ± {s_:.1e} ({dd.mean() / s_:+.1f}) |")
        rep["scaling"] = sc
        izlo, izhi = names.index("zlo"), names.index("zhi")
        print("\nsymmetry z → −z: k_zhi + k_zlo = 0")
        for vn, off in (("a", 0), ("c", P)):
            s_ = dk[:, off + izlo, :] + dk[:, off + izhi, :]
            print(f"({vn}): " + " ".join(f"L={L}: {s_[:, L].mean():+.1e} ± {s_[:, L].std(ddof=1) / np.sqrt(R):.1e}" for L in Ls))
    if case in (1, 2) and case in sn and "dk_drho_f" in sn[case]:
        coef = {"R": CASES[1]["R"]} if case == 1 else {"a": CASES[2]["a"], "b": CASES[2]["b"]}
        nr = 1 if case == 1 else 2
        print(f"\n### scaling identity (S_N: LHS {sum(v * sn[case]['dk_d' + k_] for k_, v in coef.items()):.6f}, "
              f"RHS {sum(sn[case][k_] for k_ in ('dk_drho_f', 'dk_drho_r')[:nr]):.6f})")
        for vn, off in (("a", 0), ("c", P)):
            for L in (3, 8, 12):
                if L >= W:
                    continue
                lhs = sum(v * dk[:, off + names.index(k_), L] for k_, v in coef.items())
                rhs = sum(dk[:, 2 * P + j, L] for j in range(nr))
                dd = lhs - rhs
                print(f"  ({vn}) L={L}: LHS {lhs.mean():.5f}, RHS {rhs.mean():.5f}, diff {dd.mean():+.2e} ± {dd.std(ddof=1) / np.sqrt(R):.1e}")

    if figs:
        figures(case, r, m, se, chans, refs_for, fd, kbar, d.parent)
    return rep


# ---------------------------------------------------------------------------
# figures (palette: dataviz skill reference instance, light mode)
# ---------------------------------------------------------------------------
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
COL = {"a": "#2a78d6", "c": "#1baf7a"}


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def figures(case, r, m, se, chans, refs_for, fd, kbar, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = r["names"]
    main = {1: ["R"], 2: ["a", "b"], 3: ["a", "tx", "c", "zhi"]}[case]
    W = r["W"]
    Lall = np.arange(1, W)
    fig, axes = plt.subplots(1, len(main), figsize=(4.2 * len(main), 3.6), facecolor=SURF, squeeze=False)
    for ax, nm in zip(axes[0], main):
        style(ax)
        p = names.index(nm)
        P = len(names)
        for vn, ch, off in (("a", p, -0.1), ("c", P + p, 0.1)):
            ax.errorbar(Lall + off, m[ch, 1:], yerr=se[ch, 1:], color=COL[vn], lw=1.5, marker="o", ms=4,
                        capsize=2, mec=SURF, mew=0.8, label=f"({vn})")
        for rn, rv, rs, kind in refs_for(f"(a) {nm}", nm):
            if rn == "S_N exact":
                ax.axhline(rv, color=INK2, lw=1.2, ls="--", label="S_N exact")
            elif rn == "S_N finite-L":
                ax.plot(Lall, rv[1:W], color=COL["a"], lw=1.0, ls=":", label="(a) finite-L S_N")
            elif rn == "CRN FD":
                ax.axhspan(rv - rs, rv + rs, color=INK2, alpha=0.15, lw=0, label="CRN FD ±1σ")
        ax.set_title(f"case {case}: θ = {nm}", color=INK, fontsize=10, loc="left")
        ax.set_xlabel("window L (generations)", color=INK, fontsize=9)
        ax.set_xticks(Lall)
        # zoom on the converged region
        lo = min(m[[p, len(names) + p], 3:].min(), m[[p, len(names) + p], 3:].min()) - 3 * se[[p, len(names) + p], 3:].max()
        hi = max(m[[p, len(names) + p], 3:].max(), m[[p, len(names) + p], 3:].max()) + 3 * se[[p, len(names) + p], 3:].max()
        ax.set_ylim(lo, hi)
    axes[0][0].set_ylabel("d ln k / dθ  [1/cm]", color=INK, fontsize=9)
    axes[0][-1].legend(frameon=False, fontsize=7, labelcolor=INK2, loc="best")
    fig.tight_layout()
    fig.savefig(outdir / f"fig_case{case}_L.png", dpi=150, facecolor=SURF)
    plt.close(fig)

    # tail figure: CCDF of |dl/l| per collision piece and Hill plot for the first main parameter
    nm = main[0]
    p = names.index(nm)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), facecolor=SURF)
    for ax in axes:
        style(ax)
    x = 10 ** EDGES[:-1]
    for v, vn in ((0, "a"), (1, "c")):
        h = r["hj"][:, v, p].sum(0)
        if h.sum() == 0:
            continue
        tail = np.cumsum(h[::-1])[::-1] / r["hl"].sum()
        sel = tail > 0
        axes[0].loglog(x[sel], tail[sel], color=COL[vn], lw=1.5, label=f"({vn})")
        t = r["top"][:, :, v, p].reshape(-1, mc3d.NTF)
        t = t[t[:, 0] != 0]
        chord = (t[:, 5] < 1e-6) & (t[:, 6] == 0)
        for sel_, ls_, tag in ((np.ones(len(t), bool), "-", "all"), (chord, "--", "chords"), (~chord, ":", "rest")):
            if tag != "all" and not chord.any():
                continue
            vals = np.sort(np.abs(t[sel_, 0]))[::-1]
            if vals.size < 20:
                continue
            kk = np.unique(np.logspace(1, np.log10(min(5000, vals.size - 2)), 40).astype(int))
            al = [k_ / np.sum(np.log(vals[:k_] / vals[k_])) for k_ in kk]
            axes[1].semilogx(kk, al, color=COL[vn], lw=1.5 if tag == "all" else 1.0, ls=ls_,
                             marker="o" if tag == "all" else None, ms=3, label=f"({vn}) {tag}")
    h0 = r["hj"][:, 0, p].sum(0)
    p10 = (np.cumsum(h0[::-1])[::-1] / r["hl"].sum())[np.searchsorted(x, 10.0)]
    xx = np.logspace(1, 4, 20)
    axes[0].loglog(xx, p10 * (xx / 10) ** -2.0 * np.log(xx) / np.log(10.0), color=INK2, lw=1.0, ls="--",
                   label="∝ x⁻² ln x (through (a) at x = 10)")
    if case > 1:
        axes[0].loglog(xx, p10 * 0.1 * (xx / 10) ** -1.5, color=INK2, lw=1.0, ls=":", label="∝ x⁻³ᐟ² (slope only)")
        axes[1].axhline(1.5, color=INK2, lw=1.0, ls=":")
    axes[0].set_xlabel("x", color=INK, fontsize=9)
    axes[0].set_ylabel("P(|dl/l| > x) per collision piece", color=INK, fontsize=9)
    axes[0].set_title(f"case {case}, θ = {nm}: tail of the dl/l term", color=INK, fontsize=10, loc="left")
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK2)
    axes[1].axhline(2.0, color=INK2, lw=1.0, ls="--")
    axes[1].set_xlabel("k (number of largest |D_own| used)", color=INK, fontsize=9)
    axes[1].set_ylabel("Hill estimate of α", color=INK, fontsize=9)
    axes[1].set_title("per-history own-generation D", color=INK, fontsize=10, loc="left")
    lo_, hi_ = axes[1].get_ylim()
    axes[1].set_ylim(lo_, max(hi_, 2.0) + 0.35 * (max(hi_, 2.0) - lo_))
    axes[1].legend(frameon=False, fontsize=7, labelcolor=INK2, loc="upper right", ncol=3)
    fig.tight_layout()
    fig.savefig(outdir / f"fig_case{case}_tail.png", dpi=150, facecolor=SURF)
    plt.close(fig)


def main():
    rundir = Path(sys.argv[1])
    figs = "--figs" in sys.argv
    sn = load_sn()
    allrep = {}
    for case in (1, 2, 3):
        d = rundir / f"case{case}"
        if d.exists() and any(d.glob("sens_seed*.npz")):
            allrep[case] = case_report(case, d, sn, figs)
    json.dump(allrep, open(rundir / "analysis.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
