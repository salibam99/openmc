import sys, numpy as np
from pathlib import Path
d = Path(sys.argv[1]); chs = [int(c) for c in sys.argv[2].split(",")]
fs = sorted(d.glob("sens_seed*.npz"))
Z = [np.load(f) for f in fs]
names = list(Z[0]["names"]); P = len(names)
lab = lambda ch: (["(a)"]*P+["(c)"]*P+["LR"]*2)[ch]+" "+(names*2+["rho_f","rho_r"])[ch]
for ch in chs:
    est = np.array([z["S1"][:, ch, 12].sum()/z["den"].sum() for z in Z])
    zs = (est-est.mean())/est.std(ddof=1); w = int(np.argmax(abs(zs)))
    z = Z[w]; S1 = z["S1"][:, ch, :]; den = z["den"]
    med = np.median(S1, 0)
    E = S1 - med
    nstar = int(np.argmax(abs(E[:, 12])))
    o = nstar - 12
    print(f"\n{lab(ch)}: worst seed {fs[w].name} z={zs[w]:+.2f}; est {est[w]:+.5f} vs mean {est.mean():+.5f} (sd {est.std(ddof=1):.5f})")
    print(f"  max |excess| at L=12 in active gen {nstar}; origin gen o = {o}")
    # excess of seed estimate explained by this origin
    tot_ex = (S1[:,12].sum()-med[12]*len(den))/den.sum()
    if o >= 0:
        prof = [E[o+L, L] for L in range(13) if o+L < len(den)]
        print("  E(o, L) = S1[o+L, L] - median, L=0..12:", " ".join(f"{x:+.3g}" for x in prof))
        print(f"  E(o,12)/den = {E[o+12,12]/den[o+12]*len(den)/len(den):+.3g} per-gen;"
              f" contribution to seed estimate {E[o+12,12]/den.sum():+.5f}; total excess {tot_ex:+.5f}")
        v = 0 if ch < P else 1; p = ch % P
        t = z["top"][o, v, p]  # gen o (active index) top records
        print("  top D_own records in gen o: [D_own, J, Jmax, l, s, dist, sid, gen, i]")
        for r in t[:3]: print("   ", " ".join(f"{x:.3g}" for x in r))
    print("  per-gen |E| at L=12, top 5 gens:", sorted([(round(float(abs(E[n,12])),1), n) for n in range(len(den))])[-5:])
