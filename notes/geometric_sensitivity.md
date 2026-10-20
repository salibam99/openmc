# Geometric (and cross-section) sensitivities in OpenMC — working notes

## Goal
Compute derivatives of tallies (leakage, reaction rates, mesh flux, later k-eff)
with respect to geometric surface parameters (plane offset, sphere/cylinder radius,
translation), in ONE run, without finite differences. Primary use case: geometric
uncertainty budgets of benchmarks (ICSBEP/SINBAD), sandwich rule S^T C_geom S,
and correlations between benchmarks. Secondary: gradient-based local optimization
of parametric CSG designs (tens of parameters). Freeform/mesh shape optimization is
out of scope.

Method: pathwise distance-space reparameterization (Deng et al., refs/), computed in
forward mode by hand (no autodiff), accumulated per particle like the existing
TallyDerivative machinery.

## References (refs/, local only — consult specific sections, not whole papers)
- Deng et al. 2026, "Differentiable Neutron Transport": §4.1–4.3 reparameterization,
  §5.2 track-length gradient. Their Eqs. 19–20 are loosely written; OUR CONVENTIONS
  BELOW OVERRIDE THE PAPER'S NOTATION.
- Ebiwonjumi, Forget, Peterson 2026, fixed-source sensitivity in OpenMC: Eqs. 10,
  19–21 for how per-event sensitivity terms are tallied in OpenMC.
- Chen & Wang 2025 (arXiv 2507.12757): HMF-058/066 Be-reflected spheres modelled in
  OpenMC. Candidate realistic test case (nested spherical shells).

## Reference implementations (ground truth — notes/)
- reparam_slab_check.py: collision-frame version, validated vs finite differences.
- reparam_slab_check_openmc_style.py: per-piece version matching OpenMC's event
  loop, validated vs the above. THIS is the one to mirror in C++.

## Method — per-piece rules (OpenMC framing)
Per particle and per geometric parameter j: scalar D_j = d log(weight)/d theta_j,
and 3-vector g_j = d(current position)/d theta_j. Direction omega is held fixed.
Source particles: D=0, g=0 (unless the source moves with theta).

Each advance is one "piece": from start P (collision site or crossing point) along
omega, length l = d_boundary (distance to nearest surface, any surface).

Surface the piece is bounded by: F(x; theta) = 0, n = gradF/|gradF|,
V_n = -(dF/dtheta)/|gradF|  (zero if the surface does not depend on theta_j).

    dl = (V_n - n.g) / (n.omega)         # = -(dF/dtheta + gradF.g)/(gradF.omega)

Ending at a surface (crossing or leak), with Sigma_t of the piece's material:

    D += -Sigma_t * dl
    g += omega * dl

Ending in a collision at distance d < l, u = d/l:

    D += -Sigma_t * u * dl + dl / l       # attenuation + Jacobian
    g += omega * u * dl

Collision outcome, nuclide choice and outgoing angle/energy: theta-independent, no term.
Fission daughters (eigenvalue phase): per-origin slots, see "Eigenvalue" below.

Scoring:
- Collision and leakage scores: score * D (D after this piece's update).
- Track-length scores for a piece, per tally bin v:
    1. l_v * D                      (D at the start of the piece)
    2. for every bin face crossed strictly inside the piece, with outward normal
       n_out of bin v:  -(n_out.g)/|n_out.omega|      (g at piece start)
    3. Leibniz endpoint term: if the piece ends at a surface (crossing or leak),
       add +dl to the bin containing the end point.
  Omitting term 3 is silently wrong in bins containing a moving interface
  (tested: -1.44 instead of -0.61).

Cross-section parameter (likelihood ratio, Sigma_t of material m, albedo fixed):
per piece in m: D -= distance travelled; per collision in m: D += 1/Sigma_t.
Validated for leakage only. Track-length scoring for xs parameters needs the
within-piece survival term (not just l_v*D); check how existing TallyDerivative
handles this before implementing.

## Eigenvalue (k) — decisions (supersede earlier "inherit D, g" / mean-subtraction text)
Target, with H the next-generation fission-source operator, H psi = k psi:

    d ln k/dtheta = <psi+, (dH/dtheta) psi> / <psi+, H psi>,   dk/dtheta = k * d ln k/dtheta

Units: theta in cm, so both are per cm. Output both.

Per-origin slots. A particle of generation n carries slots {D_o, g_o} for
o = n-L, ..., n (L inherited + 1 own):
- slot n (own generation) is ACTIVE: uses V_n of the surface bounding the piece;
- slots o < n are PASSIVE: V_n = 0, so dl_o = -(n.g_o)/(n.omega); same D/g rules.
At a fission in generation n (after that collision's piece update): each daughter
(generation n+1) copies slots n-L+1..n, drops slot n-L, opens slot n+1 = {0, 0}.
Slots travel with the site through the bank and combing.

Scoring (IFP, OpenMC convention: ifp_n_generation = L, ancestor index 0 is from
generation n-L): at each fission event in generation n of a particle with a full
window, num += w * D_{n-L}, den += w (the plain IFP denominator). Estimator
d ln k/dtheta = sum(num)/sum(den). NO mean subtraction, no cumulative-D
normalization. Passive slots give the del(psi+).g term (importance of a displaced
daughter); dropping them (g reset at birth) is biased.

xs/density parameters (LR) use the same windowed D slots, no g (no passive term).
Validation status: offline toy in notes/toy_eigen/ (variants (a) per-origin,
(b) g reset at birth, (c) single inherited g).

## Known caveats
- Variance: dl/l is ~1/l for collisions in short pieces → log-divergent variance;
  grazing angles add 1/|n.omega|. Finely subdivided cells are the worst case.
  Batch error bars can be misleading; validate with >=1e8 histories near interfaces.
- Never place a moving interface on a tally-bin face in tests (bin flux has a kink
  in theta there; symmetric FD is meaningless).
- Invalid when topology changes (surfaces touching/crossing, layers -> 0 thickness).
- Validation runs: analog or survival biasing only. No weight windows (position-
  dependent splitting/roulette introduces bias). Source rejection sampling against
  a moving surface needs a source term (not handled).

## Phases (do not skip ahead)
1. Offline Python: post-process track output (MG mode, analog) for the slab and a
   concentric-sphere leakage problem; match the reference scripts and central
   differences.
2. In-code, fixed source, on a feature branch from `develop`:
   - Python: extend openmc.TallyDerivative with a surface-parameter variable
     (surface id, parameter type, direction for translations); make XML read/write
     branch on variable (material is currently always written).
   - C++: per-particle g next to flux_derivs; hooks in advance / cross-surface /
     collide (history AND event-based); Surface method for dF/dtheta (planes,
     spheres, cylinders first); face + endpoint terms for mesh track-length.
   - Regression tests in MG mode vs central differences.
   - Kill test: equal-CPU figure of merit vs FD with correlated sampling.
3. Eigenvalue: per-origin {D_o, g_o} slots carried through the fission bank over
   a window of L generations, IFP-weighted ratio, no mean subtraction (see
   "Eigenvalue" above). Step 0: offline Python toy (notes/toy_eigen/) vs S_N.
   Then reuse IFP machinery (PR #3133, settings.ifp_n_generation). Validate on
   Godiva / reflected sphere vs central differences and published MCNP6
   interface sensitivities.
4. Many parameters per cell / parameter axis on tallies; variance stress tests.

## Working rules for the agent
- Follow AGENTS.md for build, test and style conventions.
- Start every new area read-only; report findings with file:line references and
  flag anything uncertain. Several claims in these notes about OpenMC internals
  (function names, eigenvalue handling of existing derivatives, what track output
  records) are UNVERIFIED — verify before relying on them.
- Before coding each update rule, derive it from the dl formula in a code comment.
- One change per session; each ends with a passing test and numbers with error
  bars vs FD or the reference scripts. "Tests pass" alone is not a result.
- Never modify AGENTS.md, .gitignore, or commit anything from notes/ or refs/.
- Ask before refactoring shared structures (particle state, source bank, tallies).

## Current task
(Edit this section each session.)
Phase 3, step 0: offline eigenvalue toy (Python/numba, no C++), notes/toy_eigen/:
1D two-region slab (fuel [0,a], reflector [a,b], vacuum), theta1 = a, theta2 = b.
derivation.md (estimator, variants a/b/c), MC vs S_N central differences for
L = 1..10, fixed-source sanity check, scaling identity
a dk/da + b dk/db = sum_m rho_m dk/drho_m. Results in results.md.

### Pending from an earlier session
Phase 1: read-only analysis. Answer, with file:line references, written to
notes/derivative_analysis.md:
1. TallyDerivative path: Python XML -> C++ parsing -> per-particle storage ->
   every update/apply site.
2. One derivative per tally? How could a parameter axis fit?
3. Do existing derivatives account for the perturbed fission source in eigenvalue mode?
4. In history- and event-based loops: where are d_boundary, collision distance,
   surface crossing and leakage handled? Where would the hooks go?
5. Where are surface gradients implemented per surface type?
6. IFP bank: how is ancestor data stored, normalized, and sent under MPI?
7. Track output: which events add a state? Are crossings recorded? Fission
   parent-daughter linkage?
