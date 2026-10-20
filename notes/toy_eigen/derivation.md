# k-eigenvalue sensitivity to surface parameters: estimator and variants

Context: notes/geometric_sensitivity.md (per-piece rules, "Eigenvalue" section).
This file fixes the definitions used by the toy in this directory. Units: lengths
in cm, cross sections in 1/cm, so d ln k/dθ and dk/dθ are in 1/cm.

## 1. Target

Let H(y ← x; θ) be the next-generation operator: the expected density of fission
sites at y produced by one neutron born (isotropically) at x. The fundamental mode
is H ψ = k ψ and the adjoint is H† ψ† = k ψ†. Differentiating and projecting on ψ†:

    dk/dθ ⟨ψ†, ψ⟩ = ⟨ψ†, (∂H/∂θ) ψ⟩
    d ln k/dθ     = ⟨ψ†, (∂H/∂θ) ψ⟩ / ⟨ψ†, H ψ⟩            (1)
    dk/dθ         = k · d ln k/dθ.

Only the kernel of the generation is differentiated; ψ and ψ† are held fixed. No
mean subtraction appears anywhere in (1).

## 2. Pathwise derivative of one generation

In the distance reparameterization (Deng et al.) a history is a sequence of
θ-independent reference variables (piece fractions u, directions, reaction
choices). Positions depend on θ and the start point. For any fixed weight function
φ over daughter sites y,

    ∂/∂θ (H†φ)(x) = E_x[ Σ_sites ( φ(y) D + ∇φ(y)·g ) ]              (2)

D = ∂ ln p/∂θ is the derivative of the log path density in reference coordinates,
accumulated with the per-piece rules. g = ∂y/∂θ is the motion of the site.

If the start point is displaced by g0 and θ is held fixed, the same rules with
V_n = 0 and initial g = g0 give

    d/dε (H†φ)(x + ε g0) = E_x[ Σ_sites ( φ(y) D^pass + ∇φ(y)·g^pass ) ].   (3)

Per-piece rules. The piece starts at P, has direction ω and length l to the
nearest surface F(x; θ) = 0, with n = ∇F/|∇F| and V_n = -(∂F/∂θ)/|∇F|:

    dl = (V_n - n·g) / (n·ω)
    ends on a surface (crossing or leak):  D += -Σ_t dl,            g += ω dl
    ends in a collision at u = d/l:        D += -Σ_t u dl + dl/l,   g += ω u dl

Reaction choice (scatter, capture, fission), ν sampling and the outgoing direction
do not depend on θ, so they add no term. The collision update of the piece that
ends in a fission happens *before* the sites are banked. The sites therefore
inherit D and g at the fission point.

1D slab specialisation. Surfaces are x = 0, x = a, x = b, all with n = +x̂, so
n·ω = μ and n·g = g (a scalar, dx/dθ). F = x - θ gives V_n = 1 on the surface
that is θ, and 0 otherwise:

| piece bounded by | V for θ1 = a | V for θ2 = b |
|------------------|--------------|--------------|
| x = 0            | 0            | 0            |
| x = a            | 1            | 0            |
| x = b            | 0            | 1            |

Crossing a surface sets g to V exactly (g + μ·(V - g)/μ = V). A passive g is
therefore zeroed at every crossing and decays as (1-u) at collisions. Its memory
over generations decays geometrically.

## 3. Importance recursion and the window L

ψ† = H†ψ†/k. Differentiate along a displacement g of the point where ψ† is
evaluated, using (3) with φ = ψ†:

    ∇ψ†(y)·g = (1/k) E_y[ Σ_sites ( ψ†(y') D^pass + ∇ψ†(y')·g'^pass ) ]   (4)

So the ∇ψ†·g term in (2) unrolls into passive contributions in all later
generations. Each one is weighted by the importance of the sites produced at that
point.

IFP replaces ψ† at a site born in generation o+1 by the expected number of fission
events among its descendants in generation o+L. Call this ψ†_L. This is OpenMC's
convention: with `ifp_n_generation = L`, a particle of generation n uses the value
of its ancestor from generation n-L (`lifetimes[0]`). With L = 1, the importance
of a daughter is its own probability of causing a fission. Proportionality
constants (k^L, population control) cancel in the ratio.

## 4. Estimator (variant a, per-origin slots)

### State
A particle of generation n carries slots {D_o, g_o} (one per parameter j) for
o = n-L_max, ..., n:
- slot n is ACTIVE: pieces use the V_n of the bounding surface;
- slots o < n are PASSIVE: V_n = 0, so dl_o = -g_o/μ. The same D/g update rules
  apply.

All slots of a particle see the same piece and the same sampled u. They differ
only through their own g_o and V_n.

### Birth: what goes into a daughter's slots
A source particle of the initial generation has all slots = 0. At a fission in
generation n (after the collision update), each of the ν daughters (generation
n+1) gets
- D_o, g_o copied from the parent, for o = n-L_max+1, ..., n. Slot n becomes
  passive in the daughter and keeps the parent's active D_n and g_n at the
  fission point;
- slot n+1 = {0, 0}, the daughter's own active slot;
- slot n-L_max is dropped (it has left the window).

Slots travel with the site through the bank and the combing step. Combing selects
sites independently of θ and preserves expected weight.

### What a slot accumulates
Slot o, from its creation (birth of a generation-o particle) to generation o+L:
- gen o: active derivative of that particle's path up to and including the
  fission collision. This is the D term of (2).
- gens o+1, ..., o+L: passive increments from displacing the generation-o site
  (g_o), carried down the descendant chain. This is the ∇ψ†_L·g term of (2),
  unrolled via (3).

### Scoring at a fission event in generation n
For each fission event of a particle with a full window (weight w = weight at the
collision, as OpenMC's `p.wgt_last()`), after that collision's piece update:

    num_L += w · D_{n-L}        (the oldest slot for window L)
    den   += w                  (plain IFP denominator, SCORE_IFP_DENOM)

    d ln k/dθ ≈ Σ_active gens num_L / Σ_active gens den,   dk/dθ = k · (that).

For a fixed L, E[num_L]/E[den] = ⟨ψ†_L, ∂H ψ⟩/⟨ψ†_L, H ψ⟩. This is the exact
derivative of the generation-o kernel weighted by a fixed function ψ†_L. The
passive slot differentiates ψ†_L along the site motion. The finite-L error is
first order in ψ†_L - ψ† and decays as DR^L (DR = dominance ratio).

The content of slot o at generation n does not depend on how much larger the
window is. So one run with L_max = 10 gives every L ≤ 10 exactly: score
D_{n-L} for each L from the same slots. The estimates for different L are
correlated, which is fine. Error bars come from independent seeds.

L = 0 is also reported as a baseline: score the own active slot D_n at its own
fission. For (a) and (b) this is ⟨1, ∂H ψ⟩/⟨1, H ψ⟩ (no adjoint weighting), and
the two coincide. (c) does not coincide with them at L = 0: its own slot also holds
the increments driven by the inherited g (see §5).

## 5. Variants

All three use the same transport, the same per-piece rules and the same scoring
(D of slot n-L at fission events in generation n). They differ only in what is
carried and which g drives each D increment.

(a) Per-origin slots, as in §4: g_o per slot; the active slot uses V_n,
passive slots use V_n = 0.

(b) g reset to 0 at birth. One g per particle, set to 0 at birth. It is driven
with V_n (active only), and its increments go into the active slot D_n. Daughters
inherit D_o for o in the window but no g. Passive slots are frozen after their
generation. This drops the ∇ψ†·g term in (2).
Prediction: biased for all L, with no convergence to the reference. The size of
the bias depends on how strongly ψ† varies where sites are displaced.

(c) A single g per particle, inherited and never reset. g is the total position
derivative dx/dθ, summed over all origins and never truncated. It is driven with
V_n of the bounding surface. Every D increment in generation n, including the part
caused by the inherited g, goes into the current slot D_n. Daughters inherit g
and the windowed D slots.

Linearity of dl in (V_n, g) gives g^(c) = Σ_{o ≤ n} g_o^(a), with no window. The
increment of D_n^(c) in generation n is Σ_o (increment of slot o in gen n)^(a).
So (c) re-attributes each passive increment to the generation in which it happens
instead of its origin generation:
- (a) weights the passive increment of origin o in generation m with the fission
  count in generation o+L, and drops it for m > o+L.
- (c) weights it with the fission count in generation m+L, and never drops it.

Prediction: by (4), (c) also converges to (1) as L → ∞, because each term of the
unrolled recursion is weighted by an L-generation importance. Its finite-L bias
differs from (a). §5.1 identifies the fixed functional it estimates. This
is an expectation from the derivation, to be tested, not a result.

### 5.1 Addendum, worked out after the first runs: what (c) estimates at finite L

Let ḡ(x) = E[g | site born at x] be the stationary mean of the inherited total g
of variant (c). Given x, a history's future is independent of its past, and the
passive increments are linear in g. So only ḡ matters. In stationarity a
daughter's g is the parent's g at the fission site, so E[g at birth | y] = ḡ(y).
Hence the term ∇φ(y)·(g_end - ḡ(y)) that (c) never credits to slot m averages to
zero, and

    E[(c), L] = ∂/∂θ ⟨ψ†_L, Φ_θ⁻¹ H_θ Φ_θ ψ⟩ / ⟨ψ†_L, H ψ⟩,   Φ_θ(x) = x + θ ḡ(x).

So (c) is the same frozen-importance functional as (a), but written in coordinates
that move with the stationary displacement field ḡ instead of physical
coordinates. Conjugation does not change the spectrum, and
⟨ψ†, [H, ḡ·∇] ψ⟩ = 0, so the L → ∞ limit is d ln k/dθ again. The finite-L biases
of (a) and (c) differ. In the toy, (c)'s is much smaller.

Dilation direction (θ = (a, b), dθ = (a, b) ε): V on each surface is x_surface.
With g = x, a piece gives dl = l. The increments become -Σ_t l (crossing) and
1 - Σ_t d (collision), which are exactly the density LR increments. g = x is
preserved, and any error e = g - x is zeroed at crossings and multiplied by (1-u)
at collisions. So in stationarity a·D_a^(c) + b·D_b^(c) = Σ_m D_ρm per slot,
history by history. Two consequences:
- the scaling identity cannot test (c), since it holds pathwise;
- (c)'s finite-L value for θ_a follows from the LR and θ_b finite-L values,
  using (c)_b = (a)_b, which holds in this 1D problem because x = a must be
  crossed to reach a fission site.

Density / cross-section parameters ρ_m (all macroscopic cross sections of region
m scaled by ρ_m, evaluated at ρ_m = 1) use likelihood ratio. Per piece in m,
D -= Σ_t,m · (distance travelled in the piece). Per collision in m, D += 1. There
is no g and no passive term, so (a), (b) and (c) are identical for these
parameters. They use the same windowed slots and the same IFP ratio.

## 6. Scaling identity (extra check)

k depends only on optical thicknesses, so k(λa, λb; ρ) = k(a, b; λρ).
Differentiating at λ = 1:

    a ∂k/∂a + b ∂k/∂b = Σ_m ρ_m ∂k/∂ρ_m              (the same with ln k).

In discrete S_N with a fixed number of cells per region this holds exactly, since
diamond difference depends only on Σh/μ. In MC both sides use the same IFP
machinery; the right side uses LR.

## 7. Variance remarks (expected)

- Collision pieces: dl/l = (V - g)/(±s), where s is the distance from the piece
  start to the bounding surface along the normal. The density of collision sites
  near a surface is finite, so E[(dl/l)^2] diverges logarithmically. The variance
  per history is (log-)infinite. The mean is finite. Seed-to-seed scatter should
  be checked for outliers and for stability of the variance estimate as the
  number of histories grows.
- Crossing pieces: -Σ_t dl ∝ 1/μ, with crossing probability ∝ |μ| at small μ.
  This is also log-divergent.
- IFP: num_L carries a sum over L generations of passive increments times a
  progeny count. The variance grows with L, roughly linearly in the toy, because
  passive g dies fast.

## 8. Simulation conventions in the toy

- Analog transport: capture and fission end the history. ν is an integer
  floor(ν̄) or floor(ν̄)+1. Isotropic scattering and birth directions.
- Each generation is combed to exactly N sites of weight 1. The analog k of a
  generation is (number of sites banked)/N.
- At most one fission per history (analog), so scoring is done from the
  per-history fission record. That record holds the slots after the collision
  update, which are the same slots the daughters copy.
- Inactive generations ≥ (source convergence) + L_max, so every scored particle
  has a full window.
- Each particle has its own counter-based random stream, keyed by (seed,
  generation, particle index). This makes CRN central differences meaningful and
  the parallel run deterministic.
