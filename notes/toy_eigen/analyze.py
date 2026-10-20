"""
Analysis of runs/: tables (markdown, printed) and figures for results.md.
Error bars are from independent seeds only. Batch (per-generation) error bars
inside one run are computed only to show how far off they are.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import mc
import sn
from model import A0, B0

HERE = Path(__file__).parent
RUNS = HERE / "runs"
LS = [0, 1, 2, 3, 5, 8, 10]

# reference palette (dataviz skill, light mode), categorical slots 1-3
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
COL = {"a": "#2a78d6", "b": "#eb6834", "c": "#1baf7a"}
LABEL = {"a": "(a) per-origin slots", "b": "(b) g reset at birth", "c": "(c) single inherited g"}


def load_sens():
    files = sorted(RUNS.glob("sens_seed*.npz"))
    est, k, s1, s2, den, gen_est, secs = [], [], [], [], [], [], []
    for f in files:
        d = np.load(f)
        S1, S2, dn = d["S1"], d["S2"], d["den"]
        est.append(S1.sum(0) / dn.sum())
        k.append(d["kgen"].mean())
        s1.append(S1.sum(0))
        s2.append(S2.sum(0))
        den.append(dn.sum())
        gen_est.append(S1 / dn[:, None, None])
        secs.append(float(d["seconds"]))
        N, n_act = int(d["N"]), int(d["n_active"])
    return dict(est=np.array(est), k=np.array(k), s1=np.array(s1), s2=np.array(s2),
                den=np.array(den), gen=np.array(gen_est), N=N, n_act=n_act,
                H=N * n_act, secs=np.array(secs), R=len(files))


def mean_se(x):
    return x.mean(0), x.std(0, ddof=1) / np.sqrt(x.shape[0])


def main():
    ref = json.load(open(HERE / "sn_results.json"))["reference"]
    dr = json.load(open(HERE / "sn_results.json"))["modes"]["dr"]
    r = load_sens()
    R, H = r["R"], r["H"]
    est, k = r["est"], r["k"]
    out = {"R": R, "H_per_seed": H, "N": r["N"], "n_active": r["n_act"],
           "seconds_per_seed": float(r["secs"].mean())}

    kmean, kse = k.mean(), k.std(ddof=1) / np.sqrt(R)
    print(f"seeds R = {R}, active histories per seed H = {H:.3g}, "
          f"{r['secs'].mean():.0f} s/seed")
    print(f"k_MC = {kmean:.6f} +- {kse:.6f}   k_SN = {ref['k']:.6f}   "
          f"z = {(kmean - ref['k']) / kse:+.2f}   DR_SN = {dr:.4f}\n")
    out["k"] = (kmean, kse)

    m, se = mean_se(est)
    dk = est * k[:, None, None]                  # dk/dtheta per seed
    mk, sek = mean_se(dk)
    refs = {"theta_a": ("dlnk_da", "dk_da"), "theta_b": ("dlnk_db", "dk_db"),
            "rho_f": ("dlnk_drho_f", "dk_drho_f"), "rho_r": ("dlnk_drho_r", "dk_drho_r")}

    # ---------------- finite-L S_N predictions (sn_frozen.py) ----------------
    fz = json.load(open(HERE / "sn_frozen_results.json"))["best"]
    fin = {"a:theta_a": np.array(fz["dlnQ_da"]), "theta_b": np.array(fz["dlnQ_db"]),
           "rho_f": np.array(fz["dlnQ_drho_f"]), "rho_r": np.array(fz["dlnQ_drho_r"])}
    # (c) obeys the dilation identity pathwise (results.md), and (c)_b = (a)_b here:
    fin["c:theta_a"] = (fin["rho_f"] + fin["rho_r"] - B0 * fin["theta_b"]) / A0
    for nm in ("a:theta_b", "b:theta_b", "c:theta_b"):
        fin[nm] = fin["theta_b"]
    out["finite_L_sn"] = {kk: v.tolist() for kk, v in fin.items()}

    # ---------------- main tables ----------------
    tables = {}
    for par in ("theta_a", "theta_b", "rho_f", "rho_r"):
        rl, rk = refs[par]
        chans = [c for c, nm in enumerate(mc.CHANNELS) if nm.endswith(par)]
        if par == "theta_b":
            chans = chans[:1]          # (a), (b), (c) identical to the last digit, see below
        lines = []
        hdr = ["L"]
        for c in chans:
            v = mc.CHANNELS[c].split(":")[0] if ":" in mc.CHANNELS[c] else "LR"
            if par == "theta_b":
                v = "a = b = c"
            hdr += [f"{v}: MC", "z vs exact"]
            if mc.CHANNELS[c] in fin or par in fin:
                hdr += ["S_N finite-L", "z vs finite-L"]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("|" + "---|" * len(hdr))
        for L in LS:
            row = [f"{L}"]
            for c in chans:
                nm = mc.CHANNELS[c]
                z = (m[c, L] - ref[rl]) / se[c, L]
                row += [f"{m[c, L]:+.5f} ± {se[c, L]:.5f}", f"{z:+.1f}"]
                fk = nm if nm in fin else (par if par in fin else None)
                if fk is not None:
                    zf = (m[c, L] - fin[fk][L]) / se[c, L]
                    row += [f"{fin[fk][L]:+.5f}", f"{zf:+.1f}"]
            lines.append("| " + " | ".join(row) + " |")
        lines.append(f"| exact (S_N) | {ref[rl]:+.5f} |" + " |" * (len(hdr) - 2))
        tables[par] = "\n".join(lines)
        print(f"### d ln k / d{par}  [S_N {ref[rl]:+.6f}]")
        print(tables[par], "\n")
        for c in chans:
            for L in (8, 10):
                out[f"dk/d{par} {mc.CHANNELS[c]} L={L}"] = (mk[c, L], sek[c, L])
        print("dk/d" + par + f" [S_N {ref[rk]:+.6f}]: " + "; ".join(
            f"{mc.CHANNELS[c]} L={L} {mk[c, L]:+.5f}±{sek[c, L]:.5f}"
            for c in chans for L in (8, 10)) + "\n")
    # largest |difference| between variants for theta_b (identical in 1D, see derivation)
    print("max |(a)-(b)|, |(a)-(c)| for theta_b over seeds and L: "
          f"{np.abs(est[:, 1] - est[:, 3]).max():.2e}, {np.abs(est[:, 1] - est[:, 5]).max():.2e}\n")

    out["dlnk_mean"] = m.tolist()
    out["dlnk_se"] = se.tolist()
    out["dk_mean"] = mk.tolist()
    out["dk_se"] = sek.tolist()

    # ---------------- per-history variance ----------------
    var_eff = est.var(0, ddof=1) * H                           # from seeds
    # naive iid per-history variance of the linearised ratio score
    Rr = est
    Ey = r["den"][:, None, None] / H
    naive = (r["s2"] / H - 2 * Rr * r["s1"] / H + Rr ** 2 * Ey) / Ey ** 2
    naive = naive.mean(0)
    print("### per-history variance sigma^2_hist = Var_seed(estimate) x H  [cm^-2]")
    print("(naive = iid-history formula; ratio eff/naive = generation-correlation inflation)")
    print("| channel | " + " | ".join(f"L={L}" for L in LS) + " | eff/naive (L=8) |")
    print("|" + "---|" * (len(LS) + 2))
    vt = []
    for c, nm in enumerate(mc.CHANNELS):
        row = f"| {nm} | " + " | ".join(f"{var_eff[c, L]:.3g}" for L in LS) + \
            f" | {var_eff[c, 8] / naive[c, 8]:.2f} |"
        vt.append(row)
        print(row)
    print()
    out["var_hist_eff"] = var_eff.tolist()
    out["var_hist_naive"] = naive.tolist()
    # relative variance per history (FOM ingredient): sigma^2/ref^2
    print("relative per-history variance (sigma_hist/ref)^2 at L=8: " + ", ".join(
        f"{nm} {var_eff[c, 8] / ref[refs[nm.split(':')[-1]][0]] ** 2:.3g}"
        for c, nm in enumerate(mc.CHANNELS)) + "\n")

    # ---------------- batch vs seed error bars ----------------
    g = r["gen"]                              # [seed, gen, ch, L]
    batch_se = g.std(1, ddof=1) / np.sqrt(g.shape[1])        # naive per-run SE from gens
    true_sd = est.std(0, ddof=1)                             # per-run sd across seeds
    ratio = batch_se.mean(0) / true_sd
    print("### batch-SE / seed-SD (per run; <1 means batch error bars are too small)")
    for c, nm in enumerate(mc.CHANNELS):
        print(f"{nm:10s} " + " ".join(f"L={L}:{ratio[c, L]:.2f}" for L in LS))
    out["batch_over_seed"] = ratio.tolist()
    # lag-1 autocorrelation of per-generation estimates (averaged over seeds)
    gc_ = g - g.mean(1, keepdims=True)
    ac1 = (gc_[:, 1:] * gc_[:, :-1]).sum(1) / (gc_ ** 2).sum(1)
    out["gen_autocorr_lag1"] = ac1.mean(0).tolist()
    print("lag-1 autocorr of generation estimates, L=8: " + ", ".join(
        f"{nm} {ac1.mean(0)[c, 8]:+.2f}" for c, nm in enumerate(mc.CHANNELS)) + "\n")

    # ---------------- tails ----------------
    print("### tails: per-seed estimates (Shapiro p, excess kurtosis, max |z| among seeds)")
    tails = {}
    for c, nm in enumerate(mc.CHANNELS):
        x = est[:, c, 8]
        zmax = np.max(np.abs((x - x.mean()) / x.std(ddof=1)))
        p = stats.shapiro(x).pvalue
        kurt = stats.kurtosis(x)
        # per-generation level (n_seeds * n_gens samples): heavier tails show here
        xg = g[:, :, c, 8].ravel()
        kurt_g = stats.kurtosis(xg)
        tails[nm] = dict(shapiro_p=p, kurt_seed=kurt, zmax_seed=zmax, kurt_gen=kurt_g)
        print(f"{nm:10s} shapiro p = {p:.2f}  kurt(seed) = {kurt:+.2f}  max|z| = {zmax:.2f}  "
              f"kurt(generation) = {kurt_g:+.1f}")
    out["tails_L8"] = tails
    # stability of the variance estimate: first half vs second half of the seeds
    h = R // 2
    v1 = est[:h].var(0, ddof=1)
    v2 = est[h:].var(0, ddof=1)
    print("variance ratio (seeds 1st half / 2nd half), L=8: " + ", ".join(
        f"{nm} {v1[c, 8] / v2[c, 8]:.2f}" for c, nm in enumerate(mc.CHANNELS)))
    # variance vs run length: per-run estimate over the first half of the generations
    nh = g.shape[1] // 2
    est_half = r["s1"] * 0
    S1h, denh = [], []
    for f in sorted(RUNS.glob("sens_seed*.npz")):
        d = np.load(f)
        S1h.append(d["S1"][:nh].sum(0) / d["den"][:nh].sum())
    est_half = np.array(S1h)
    vh = est_half.var(0, ddof=1) * H / 2
    print("per-history variance from half-length runs / full-length, L=8: " + ", ".join(
        f"{nm} {vh[c, 8] / var_eff[c, 8]:.2f}" for c, nm in enumerate(mc.CHANNELS)) + "\n")
    out["var_halfrun_over_full_L8"] = (vh[:, 8] / var_eff[:, 8]).tolist()

    # ---------------- scaling identity ----------------
    print("### scaling identity  a dk/da + b dk/db - (dk/drho_f + dk/drho_r)")
    print(f"S_N: lhs = {ref['scaling_lhs']:.8f}, rhs = {ref['scaling_rhs']:.8f}")
    scal = {}
    rows = ["| L | variant | LHS | RHS | LHS - RHS (z) |", "|---|---|---|---|---|"]
    for var, (ca, cb) in {"a": (0, 1), "c": (4, 5)}.items():
        for L in LS:
            lhs = A0 * dk[:, ca, L] + B0 * dk[:, cb, L]
            rhs = dk[:, 6, L] + dk[:, 7, L]
            dm, dse = mean_se(lhs - rhs)
            lm, lse = mean_se(lhs)
            rm, rse = mean_se(rhs)
            scal[f"{var} L={L}"] = (lm, lse, rm, rse, dm, dse)
            rows.append(f"| {L} | ({var}) | {lm:.5f} ± {lse:.5f} | {rm:.5f} ± {rse:.5f} | "
                        f"{dm:+.5f} ± {dse:.5f} ({dm / dse:+.1f}) |")
    rows.append(f"| S_N | | {ref['scaling_lhs']:.5f} | {ref['scaling_rhs']:.5f} | "
                f"{ref['scaling_lhs'] - ref['scaling_rhs']:+.1e} |")
    print("\n".join(rows) + "\n")
    tables["scaling"] = "\n".join(rows)
    out["scaling"] = scal

    # ---------------- MC central differences (CRN) ----------------
    fdf = sorted(RUNS.glob("fd_seed*.npz"))
    if fdf:
        hh = float(np.load(fdf[0])["h"])
        ka = np.array([(np.load(f)["ap"] - np.load(f)["am"]) / (2 * hh) for f in fdf])
        kb = np.array([(np.load(f)["bp"] - np.load(f)["bm"]) / (2 * hh) for f in fdf])
        # S_N central difference with the same h (isolates the O(h^2) truncation)
        nq, nc = 128, 1000
        sna = (sn.solve(a=A0 + hh, nf=nc, nr=nc, nquad=nq) - sn.solve(a=A0 - hh, nf=nc, nr=nc, nquad=nq)) / (2 * hh)
        snb = (sn.solve(b=B0 + hh, nf=nc, nr=nc, nquad=nq) - sn.solve(b=B0 - hh, nf=nc, nr=nc, nquad=nq)) / (2 * hh)
        d128 = json.load(open(HERE / "sn_results.json"))["S128_n1000"]
        trunc_a, trunc_b = sna - d128["dk_da"], snb - d128["dk_db"]
        # spread of a single (uncorrelated) k estimate, to show the CRN gain
        sd_k = k.std(ddof=1)
        fdres = dict(h=hh, R=len(fdf),
                     dk_da=(ka.mean(), ka.std(ddof=1) / np.sqrt(len(ka))),
                     dk_db=(kb.mean(), kb.std(ddof=1) / np.sqrt(len(kb))),
                     sn_fd_trunc_a=trunc_a, sn_fd_trunc_b=trunc_b,
                     crn_gain_a=(np.sqrt(2) * sd_k / (2 * hh)) / ka.std(ddof=1),
                     crn_gain_b=(np.sqrt(2) * sd_k / (2 * hh)) / kb.std(ddof=1))
        out["mc_fd"] = fdres
        print(f"### MC CRN central differences, h = {hh} cm, {len(fdf)} seeds")
        for nm, x, rk, tr in (("a", ka, "dk_da", trunc_a), ("b", kb, "dk_db", trunc_b)):
            mm, ss = x.mean(), x.std(ddof=1) / np.sqrt(len(x))
            print(f"dk/d{nm}: MC-FD {mm:+.5f} ± {ss:.5f}; S_N exact {ref[rk]:+.5f}, "
                  f"S_N FD(h) {ref[rk] + tr:+.5f} (truncation {tr:+.1e}); "
                  f"z vs S_N FD(h) = {(mm - ref[rk] - tr) / ss:+.1f}")
        print(f"CRN sd reduction factor vs independent runs: a {fdres['crn_gain_a']:.1f}, "
              f"b {fdres['crn_gain_b']:.1f}\n")

    json.dump(out, open(HERE / "analysis_results.json", "w"), indent=1, default=float)
    figures(m, se, var_eff, ref, refs, fin)
    Path(HERE / "tables.md").write_text(
        "\n\n".join(f"### {k_}\n\n{v}" for k_, v in tables.items()) + "\n")


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def figures(m, se, var_eff, ref, refs, fin):
    Lall = np.arange(mc.LMAX + 1)
    off = {"a": -0.12, "b": 0.0, "c": 0.12}

    # Figure 1: d ln k / da vs L, full range and zoom
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURF)
    for ax, variants, title in ((axes[0], "abc", "all variants"),
                                (axes[1], "ac", "zoom: (a) and (c)")):
        style(ax)
        ax.axhline(ref["dlnk_da"], color=INK2, lw=1.2, ls="--", zorder=1)
        for v in variants:
            c = {"a": 0, "b": 2, "c": 4}[v]
            if f"{v}:theta_a" in fin:
                ax.plot(Lall, fin[f"{v}:theta_a"], color=COL[v], lw=1.0, ls=":", zorder=2,
                        label=f"{v}: S_N finite-L prediction" if ax is axes[1] else None)
            ax.errorbar(Lall + off[v], m[c], yerr=se[c], color=COL[v], lw=1.5, marker="o",
                        ms=5, capsize=2.5, mec=SURF, mew=1.0, label=LABEL[v])
            if ax is axes[1] or v == "b":     # (a), (c) coincide at this scale
                ax.annotate(f"({v})", (Lall[-1] + off[v], m[c, -1]), xytext=(8, 0),
                            textcoords="offset points", color=INK2, fontsize=9, va="center")
        ax.annotate("S_N", (0, ref["dlnk_da"]), xytext=(0, 4), textcoords="offset points",
                    color=INK2, fontsize=9)
        ax.set_xlabel("window L (generations)")
        ax.set_ylabel("d ln k / da  [1/cm]")
        ax.set_title(title, color=INK, fontsize=10, loc="left")
        ax.set_xticks(Lall)
        ax.set_xlim(-0.5, mc.LMAX + 1.2)
    lo = min(m[[0, 4], 1:].min() - 3 * se[[0, 4], 1:].max(), ref["dlnk_da"] - 0.002)
    hi = max(m[[0, 4], 1:].max() + 3 * se[[0, 4], 1:].max(), ref["dlnk_da"] + 0.002)
    axes[1].set_ylim(lo, hi)
    axes[0].legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper right")
    axes[1].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower right")
    fig.suptitle("Interface a: IFP-weighted d ln k/da vs window L (±1σ over seeds)",
                 color=INK, fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(HERE / "fig_dlnk_da.png", dpi=150, facecolor=SURF)
    plt.close(fig)

    # Figure 2: d ln k / db and the density LR derivatives (variants coincide for b)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), facecolor=SURF)
    for ax, (c, key, ttl) in zip(axes, ((1, "dlnk_db", "outer boundary b (a, b, c identical)"),
                                        (6, "dlnk_drho_f", "fuel density ρ_f (LR)"),
                                        (7, "dlnk_drho_r", "reflector density ρ_r (LR)"))):
        style(ax)
        ax.axhline(ref[key], color=INK2, lw=1.2, ls="--")
        fk = {1: "theta_b", 6: "rho_f", 7: "rho_r"}[c]
        ax.plot(Lall, fin[fk], color=COL["a"], lw=1.0, ls=":", label="S_N finite-L")
        ax.errorbar(Lall, m[c], yerr=se[c], color=COL["a"], lw=1.5, marker="o", ms=5,
                    capsize=2.5, mec=SURF, mew=1.0)
        ax.annotate("S_N", (0, ref[key]), xytext=(0, 4), textcoords="offset points",
                    color=INK2, fontsize=9)
        ax.set_title(ttl, color=INK, fontsize=10, loc="left")
        ax.set_xlabel("window L (generations)")
        ax.set_xticks(Lall)
    axes[0].set_ylabel("d ln k / dθ  [1/cm or 1]")
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower right")
    fig.tight_layout()
    fig.savefig(HERE / "fig_dlnk_db_rho.png", dpi=150, facecolor=SURF)
    plt.close(fig)

    # Figure 3: per-history variance vs L for theta_a
    fig, ax = plt.subplots(figsize=(6, 4), facecolor=SURF)
    style(ax)
    for v in "abc":
        c = {"a": 0, "b": 2, "c": 4}[v]
        ax.plot(Lall, var_eff[c], color=COL[v], lw=1.5, marker="o", ms=5, mec=SURF,
                label=LABEL[v])
        ax.annotate(f"({v})", (Lall[-1], var_eff[c, -1]), xytext=(8, 0),
                    textcoords="offset points", color=INK2, fontsize=9, va="center")
    ax.plot(Lall, var_eff[1], color=INK2, lw=1.2, ls=":", marker="s", ms=4, label="θ = b (all variants)")
    ax.set_yscale("log")
    ax.set_xlabel("window L (generations)")
    ax.set_ylabel("per-history variance  [1/cm²]")
    ax.set_xticks(Lall)
    ax.set_xlim(-0.5, mc.LMAX + 1.2)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
    ax.set_title("Var(seed estimate) × histories per seed", color=INK, fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(HERE / "fig_variance.png", dpi=150, facecolor=SURF)
    plt.close(fig)


if __name__ == "__main__":
    main()
