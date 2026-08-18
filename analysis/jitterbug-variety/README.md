# Jitterbug configuration-variety derivation record

Bead `inviscid-qvf.1`, 2026-08-11.

**This is a derivation record, not part of the build.** Maven never sees it. It exists because
these scripts are the evidence base for four retractions of claims that had been recorded as
established, and a retraction whose evidence has evaporated is just folklore.

Superseded as the *living* harness by the Java port under
`src/{main,test}/java/com/chiralbehaviors/inviscid/jitterbug/`. Kept because it is an
independent second implementation: when the Java and the Python disagree, one of them is wrong,
and that is worth more than either alone.

## What was overturned

Four claims recorded in T2/T3 and in the `inviscid-qvf` epic as "computed and verified":

1. **The free-dynamics conservation law was backwards.** `a` is not a cyclic coordinate
   (`dL/da = ½·M_eff'(a)·ȧ² ≠ 0`), so `p = M_eff·ȧ` is *not* conserved. `E` is. Measured at
   rtol 1e-12: `p` spread √3, `E` spread 1.0000000, `ȧ` spread √3 (not 3).
2. **The 1-DOF premise was never measured.** Shared vertices take 48 body DOF to 12, not to 1.
   The linkage has **6 internal DOF**, all extending to finite motions.
3. **`a=60` is not a bifurcation** — rank is constant 36 through it. But genuine branch points
   do lie on the symmetric path, at **a=90 and a=270** (local dimension 7 vs 6). The instinct
   was right; the location and the reasoning were not.
4. **The ħ derivation is dead** four independent ways, taking the one-formable-speed and
   depth-as-refractive-index results with it into "unsupported".

Plus two findings nobody was looking for: the model **self-intersects on (60°, 120°)**, making
a=60 the *contact boundary* and the admissible symmetric path the interval [−60, 60] with no
2π circuit; and there is **no potential energy anywhere in the model**, so its six DOF are six
zero-frequency modes and no dispersion relation can exist yet (bead `inviscid-qvf.2`).

> **Amended 2026-08-12 (USER DECISION 16).** The clause above — "the admissible symmetric path
> is the interval [−60, 60] with no 2π circuit" — is a *verdict* on the self-intersection
> measurement, and it has been WITHDRAWN. Interference is permitted; the configuration space is
> the whole variety. The measurement itself (struts cross throughout (60,120) ∪ (240,300)) is
> unaffected and stands. Sweeps in `jb_h`…`jb_n` therefore exclude no region.

## Files

Primary derivation (`jb_*`), in dependency order:

| file | what it establishes |
|---|---|
| `jb_a_family.py` | symmetric family; forced chirality; a_ico = 22.238756093°; volumes 20/18.512296/4/1 |
| `jb_b_variety.py` | linkage constraints, Jacobian, rank/nullity — the 6-DOF result |
| `jb_c_branches.py` | rank scan, finite-motion continuation, 1344 tetrahedron seatings |
| `jb_d_tet.py` | explicit tetrahedron construction, first reachability pass |
| `jb_e_tighten.py` | reachability with exact Newton projection at every step |
| `jb_f_components.py` | random-restart connectivity for the targets the direct walk missed |
| `verify_critic.py` | independent re-derivation of the conservation error and the collision boundary |
| `verify_pathcollide.py` | interpenetration count along the demonstrated VE→tetrahedron path |

Potential-energy candidate survey (`jb_h` … `jb_n`, 2026-08-12, bead `inviscid-6dp`).
These test whether any of the shape functionals on the table can serve as the `V` the model
does not have. Written to *break* `V = −k·Vol_hull`, which was proposed because it agrees with
Fuller (440.05, 466.14) — and agreement with a prior expectation is when confirmation bias is
cheapest.

| file | what it measures |
|---|---|
| `jb_h_hull_control.py` | M0 control: reproduces C5's hull landmarks under Qhull; reproduces the non-simplicial **overcount trap** on purpose (brute-force hull gives V(0)=32, not 20); verifies cosphericity independently; measures the square faces' out-of-plane term |
| `jb_i_hull_smoothness.py` | M1: is `Vol_hull` C2 at the VE? Central second difference vs step size, the analytic `|a|` decomposition, one-sided derivatives, and the corner slope to 1e-9 |
| `jb_j_internal_frame.py` | the shared machinery: a **6-D chart** on the linkage variety (exact rigid parameterisation, explicit linear gauge removing the 6 global modes, Newton projection). Everything in M2–M5 differentiates through this |
| `jb_k_hull_hessian.py` | M2/M3: gradient and Hessian **inertia** of `−Vol_hull` on all six internal directions at the VE, the hull maximum, the icosahedron, and the octahedron; Sylvester-invariance and zero-mode audits |
| `jb_l_vertex_potentials.py` | M4/M5: Thomson `Σ1/r_ij` raw and radius-normalised, and strut repulsion `Σ1/d_ij` over the `jb_g` clearance machinery, with the divergence stated rather than regularised |
| `jb_m_kink_crosscheck.py` | cross-checks the two non-smoothness claims against independent machinery: a reference segment-segment distance, and the hull's **facet combinatorics** (14 planes at a=0, 20 everywhere else) |
| `jb_n_global_search.py` | random-restart projected ascent/descent on the variety — is each claimed ground state global, or only a minimum of the 1-D slice? |

Kernel-family survey (`jb_o` … `jb_q`, 2026-08-12, bead `inviscid-qvf.2`). The prior survey killed
two candidates and left one structural hypothesis: that the non-smooth ones (`Vol_hull`, strut
clearance) fail *because* they are **witness-selection** functionals, so the surviving family is
smooth all-pairs kernels. These files sweep that family and attack the hypothesis.

| file | what it measures |
|---|---|
| `jb_o_kernel_family.py` | K0/K1/K2: `V_f = Σ f(‖r_i−r_j‖²)` over the 12 shared vertices for `f` = `1/r^p` (p=1,2,3,6,12), Gaussians, and quadratic spread — raw **and** centroid-normalised. Argmin over the whole circle, smoothness by h-scaling against two opposite controls, 6-D gradient, Hessian **inertia**. Also establishes the two structural facts that decide how the table reads: 24 of the 66 vertex pairs are **frozen struts**, and **a ↔ 180−a is an exact isometry** |
| `jb_p_witness_falsify.py` | K3: five attacks on the witness-selection hypothesis, including the one that lands |
| `jb_q_strut_kernels.py` | K4: the same kernels on the 24 strut **midpoints**, plus strut-axis and midpoint+axis variants — does the choice of primitive move the ground state? |

Headlines. **No smooth all-pairs kernel measured moves the raw ground state off the vector
equilibrium** (10 kernels × vertices and struts, all at the VE, all smooth, all inertia (6,0,0)).
The witness-selection hypothesis is **REFUTED as stated** — hull *surface area* is a witness
functional and is smooth at the VE, and `λ_max` vs `λ_sum` of the same eigen-decomposition split
on smoothness. The operative property is narrower: selecting a **proper subset** of a tied orbit,
which is just non-smoothness of max/min over a tied family. Details and verdicts in T2
`inviscid/V-kernel-family-survey.md`.

Frequency spectrum (`jb_r` … `jb_t`, 2026-08-14, bead `inviscid-qvf.2`). Everything before this was
**inertia only** — the ±/0 eigenvalue count, which is metric-independent by Sylvester's law and says
nothing about how fast anything moves. USER DECISION 17 fixed the ground state (the VE) and the
variant (RAW, on real distances) and left the **kernel** free, so the kernel is exactly what the
spectrum measures.

| file | what it measures |
|---|---|
| `jb_r_mass_metric.py` | the internal mass metric `M = Pᵀ diag(m) P` restricted to the six internal directions, built independently of `src/test/java/.../InternalMassMetric.java`. **Validation anchor:** contracted with the symmetric-path generator it must reproduce `M_eff(a) = (2/3)sin²a + 1/3` (point masses) and `+ 1/12` (uniform laminae). Plus total-mass, positive-definiteness, chart-covariance, and gauge-leak checks |
| `jb_s_frequency_spectrum.py` | S0 control (re-measured, not inherited), then `H v = ω² M v` at the VE for nine raw kernels × both mass models; zero-mode count with a tolerance sweep; the symmetry-block decomposition that reduces the spectrum to three numbers; ratio tables; the exponent at which the point-mass mode ordering flips |
| `jb_t_modes_primitive_offpath.py` | S3 mode geometry as block **traces** (invariant under the eigensolver's arbitrary basis inside a degenerate block); S4 raw vertex vs raw strut-midpoint spectra; S5 off-equilibrium, using a second genuinely different chart (`OriginFrame`, rotation pivoted at the origin) to *measure* the chart-dependence rather than quote one chart's answer |

Headlines. Six real frequencies exist, with **zero** zero-modes at every tolerance from 1e-4 to
1e-12. The 2+3+1 degeneracy makes `H` and `M` both scalar on each symmetry block, so
`ω_b = √(h_b/m_b)` and the whole spectrum is three numbers. **The frequency ratios are NOT
kernel-invariant** (ω_S/ω_D spans 1.475×, ω_T/ω_D 1.184× across the nine kernels) — the kernel is a
real physical parameter, not an overall stiffness scale. The mass model rescales each block by a
fixed kernel-independent factor (√(8/5), √2, 2) and **changes the mode ordering** for five of nine
kernels. Raw vertex and raw strut-midpoint kernels share a ground state but **not** a spectrum. And
every inverse power — plus the narrow Gaussian — has a **second stable equilibrium inside the
interference band (60,120)**, legal only under DECISION 16. Details and verdicts in T2
`inviscid/first-frequency-spectrum.md`.

Riemannian Hessian (`jb_u`, 2026-08-15, bead `inviscid-qvf.9`). Off a critical point the chart
Hessian is not a tensor, so no spectrum away from an equilibrium was quotable. This builds the
Riemannian Hessian of the mass metric, `(Hess V)_ab = ∂_a∂_b V − Γ^c_ab ∂_c V`, over a genuine local
chart `q → Newton-project(x0 + Bq)` that stays on the variety.

| file | what it measures |
|---|---|
| `jb_u_riemannian_hessian.py` | U0 analytic control (polar coordinates, where Γ is textbook and a linear potential must have Riemannian Hessian **exactly zero** while its chart Hessian is O(1)); U1 control at the VE, reproducing the recorded ω² to 1.153e-7 relative and the recorded chart Hessian to 1.835e-7; U2a the gauge diagnostic; U2b the deliverable at the icosahedron across nine kernels × two primitives × two mass models; U2c the measurement that the metric-form fix is a no-op in the measuring chart, for both projections tested; U3 nonlinear-reparameterisation invariance **with the Γ:=0 arm printed beside it** so the test's teeth are measured; U4 connection checks, including U4d off the symmetric path and U4e a negative result; U5 step-size sweep; a **GATE** block that sets the process exit code |

**These are not frequencies away from the VE.** At the icosahedron `|dV| = 3.06`, so the
configuration is not an equilibrium and nothing oscillates about it. The generalised eigenvalues of
(Hess V, g) there are chart-invariant **local curvature scales** of the pair (V, g) — which is
exactly why making them chart-invariant was worth doing, and exactly why they may not be read as a
vibrational spectrum, a dispersion relation or a wave speed. At a critical point the two readings
coincide, so the **VE** numbers in the record are frequencies and are unaffected. They are still
written ω² below, for continuity with the record.

Headlines, each with the row it comes from — the script prints these same numbers and the same
definitions. At the icosahedron the two charts' naive Hessians disagree by **51.2% of the spectrum
span** (0.852682 absolute); their Riemannian ω² agree to **3.138e-6 relative, worst over all 36
kernel × primitive × mass-model combinations** (2.270e-7 for the headline point-mass row) — at
`h = 1e-3`, which is the **largest** step size at which that worst-of-36 meets the 1e-5 criterion:
the same statistic is 2.8e-5 at `h = 3e-3` and 4.1e-7 at `h = 3e-4`, pure O(h²) truncation, swept in
U2b(iv). The criterion is a threshold on that h, not a property of the method. A
reparameterisation with `Dφ = I` — the rows that isolate the Christoffel term, because they leave the
metric untouched — moves ω² by **3.738e-6**, a teeth ratio of **6.73e+05** against its own Γ:=0 arm.
The *linear* remix the earlier covariance check used gives a teeth ratio of **1.03** and therefore
proves nothing. The `Dφ ≠ I` nonlinear row sits at 3.236e-4, but so does the pure-linear control
(3.554e-4 with Γ, 3.670e-4 without): that magnitude is linear-remix conditioning present in both
arms, not a Christoffel residual, and it is not the headline.

**And a finding the brief did not anticipate:** a chart here is a *section* of the 12-D constrained
variety (6 internal + 6 rigid), and two sections can be tilted differently relative to the rigid
orbit. They are — at the icosahedron **three** of the origin-pivot slice's six principal angles with
the centroid-pivot slice exceed 1°, and **the largest of them is 31.9°** (both computed; the count
and the maximum are two different numbers and an earlier version of this line quoted one as the
other). So part of the recorded "45% chart
disagreement" is two **different 6-D subspaces** being compared, not the (dV)·D²ψ anomaly. The
invariant object is the Riemannian Hessian of the **momentum-free (mechanical-connection)** metric
`Wh = W − WZ(ZᵀWZ)⁻¹ZᵀW`; the section metric is *not* invariant (21–31% residual, entirely in the
triplet). The centroid-pivot gauge happens to be exactly momentum-free (5e-17) **on the symmetric
path** and is not off it (1e-4), so every spectrum in the record is safe and anything transverse is
not.

**The fix did not move the deliverable** (U2c). In the centroid chart at every on-path point the
section and momentum-free forms give the *same* metric to identity level — at the icosahedron, point
mass, `max |g_sec − g_horiz| = 6.9e-18`, `max |dg_sec − dg_horiz| = 5.2e-11` and
`max |Γ_sec − Γ_horiz| = 1.6e-9` against `|Γ| = 8.0e-2` — because every term of
`d(WZ(ZᵀWZ)⁻¹ZᵀW)` contracted as `Dᵀ(·)D` carries a factor `ZᵀWD`, which vanishes there. So the
whole "section 2.1e-1" column lives in the *origin* chart: the momentum-free form repairs the
**confirming** chart so it agrees with the **measuring** one. The quoted ratios are **measured to be
independent of the metric-form choice for the two projections tested** — the W-orthogonal one and
U4e's Euclidean rival — to 2.6e-8 relative in Γ and 1.2e-10 in the spectra. That is a measurement
over two forms, *not* a theorem over every equivariant projection: the two vanish against the section
form by **different** factors (`ZᵀWD` = 1.8e-10 and `ZᵀD` = 8.8e-9 respectively), and the second
vanishes only because `jb_j`'s linear gauge happens to be Euclidean. An earlier version of this line
said "provably independent", resting on a claim that every such projection carries the same vanishing
factor, which is false. Off the symmetric path the same statistic is 8.9e-3 — five and a half decades
larger — so none of this transfers to `inviscid-qvf.4`.

**A correction to a recorded result.** The two mass models disagree on the *sign* of the ratio drift
from the VE to the icosahedron: point T/D 1.022237 → 1.044917 (rises), lamina T/D 1.142895 →
1.141515 (falls). Equivalently `(lamina T/D)/(point T/D)` is exactly `√2/√(8/5) = 1.118034` at the VE
and 1.092446 at the icosahedron. That ratio-of-ratios is only 2.3% off, but **it is not the size of
the error** — most of it cancels. Per block, jb_s S2c's factors `√(8/5), √2, 2` overstate the
measured lamina/point ratio at the icosahedron by **+13.7% (D) / +16.4% (T) / +48.7% (S)**.

And the cause is sharper than "the block values are VE-specific". `ω²_b = h_b/m_b`, and S2c's factor
is `√(m_pt/m_lam)`, which assumes the block stiffness `h_b` does not depend on the mass model.
Recomputing `√(m_pt/m_lam)` *live* from the icosahedron's own mass metric still misses by
**+13.7% / +18.3% / +15.2%**, so no re-derived values rescue the formula. What actually happens is
that **off a critical point the mass model enters twice** — through `m_b` *and* through `Γ`, since
`Hess = H_naive − Γ·dV` and `Γ` is built from `g`. At a critical point `dV = 0`, the Riemannian
Hessian *is* the chart Hessian of V, which carries no mass model at all: measured `h_lam/h_pt =
1.000000` in every block at the VE, and 0.714–0.773 at the icosahedron. So the recorded factors hold
**at critical points only**; converting between mass models anywhere with `dV ≠ 0` requires
recomputing the Riemannian Hessian in that model.

**What the file cannot decide** (U4e): chart agreement does *not* select the mechanical connection.
A Euclidean-orthogonal projection along span(Z) — same kernel, same equivariance, not the kinetic
energy — reaches the same chart agreement to five significant figures (2.5747560e-7 vs 2.5747506e-7)
while giving spectra 1.99e-6 apart. The
mechanical connection is chosen by *derivation* (Riemannian submersion at zero momentum / the Eckart
frame), not by any measurement here. Trap for anyone re-testing: under the point model
`W = (1/48)·I₇₂` exactly, so the two projections coincide identically — use lamina.

Details, with all four declarations on every number, in T2 `inviscid/qvf.9-riemannian-hessian.md`.

Transverse curvature (`jb_v`, 2026-08-15, bead `inviscid-qvf.4`). `jb_u` made an off-equilibrium
spectrum quotable but produced none along the path. This sweeps the fundamental domain and asks the
question `qvf.4` was opened for: is the symmetric 1-DOF sector a transverse **valley floor** or a
**ridge**, and where does it change.

> **The file is called `transverse_curvature`, not `transverse_stability`, and that is a scope
> correction rather than a naming preference.** Off a critical point an eigenvalue of `(Hess V, g)`
> is a curvature invariant, not a frequency (`qvf.9` corollary i), and everywhere on this path
> except its critical points the system is *moving*. A dynamical stability verdict for a moving
> trajectory is the normal variational equation along it — a Floquet-type problem needing an
> equation of motion integrated, which this project has never done. What is delivered is the sign
> of the curvature of `V` in the five directions transverse to the path. The bead's acceptance
> criterion also asks for `a ∈ [−60, 60]`; that range predates USER DECISION 16 and quotes an
> admissibility verdict since withdrawn, so the sweep runs the fundamental domain `[0, 90)`.

| file | what it measures |
|---|---|
| `jb_v_transverse_curvature.py` | V0 the chiral tetrahedral group, *derived* from the forced chirality rather than tabulated, its 72×72 ambient action, its push-down to the chart, and the three character projectors, with the multiplication law, the isometry property, the projector algebra, the ranks, the properness of the 12 rotations and the group-invariance of `V` itself all checked; V1 the singlet **is** the path tangent at every swept angle, with the same statistic on a transverse direction as its teeth; V1b the gradient has **no transverse component anywhere on the path**; V2 why `jb_u`'s mass-metric labelling cannot be used here; V3 the sweep and the 36-combination **sign map**; V3c the sign changes located by bisection; V3d their error bar — reproducibility across `h`, the slope that makes a root well posed, and the absolute scalar deviation converted to an angle; V4 the ratios, and the square-root cross-check onto the record's `ω` ratios; V5 the critical points, the only rows carrying a stability reading in the ordinary sense; V6 the `a = 60` pole approached at three step sizes and V6b all twenty inverse powers swept *inside* the band; V7 the `a = 90` chart failure, the same coverage inside *its* band, the `a → 180−a` isometry, and a second pole; V8 metric form and chart invariance re-measured for this sweep; V9 the `h` sweep with a per-threshold range; a **GATE** of 56 rows that sets the exit code |

**The block labelling had to be rebuilt, and that is not a stylistic choice.** `jb_u` labels blocks
by the eigenspaces of the mass metric, which is kernel-independent and immune to the D/T ordering
flips that defeat sort-position labelling. It is strictly better than sorting and it is *still not
enough here*: in closed form — no stencil, no truncation — the momentum-free mass metric at `a = 60`
has spectrum **`[1/32 ×5, 5/96]`**. The doublet mass is `1/32` at every `a` and the triplet mass
falls through it *exactly at the octahedron*, so around `a = 60` the mass metric cannot separate the
two blocks at all and `blocks_by_irrep` correctly **refuses** (measured: it raises for `a` in
`[59.9, 60.5]`). The replacement is a **symmetry-character** labelling: the 12 rotations are
constructed, their ambient action is built by *matching* rather than from a labelling convention
(residual `2.2e-16`, and the same action verified to fix the configuration at **all 180 swept
angles** plus the eleven guard-band angles the sweep excludes, worst `2.2e-16`), and the irrep
projectors come from the character table. Ranks `(1, 2, 3)` hold at every angle including
`a = 60`, and in both charts — the origin-chart ranks are now *asserted* rather than printed, which
they were not when this paragraph was first written. Everything in the labelling is built from the **closed-form** chart
Jacobian and its own metric, so it carries no truncation: the representation law, the isometry and
the projector algebra all sit at roundoff (`1.0e-15`, `9.3e-16`, `6.7e-16`). Mixing the closed-form
`D` with the stencil's `g` put all three at `3e-9` — the `O(h²)` gap between the two metrics — which
is how the inconsistency was found.

**The singlet is the path tangent at every swept angle**, verified rather than assumed: the record
had `overlap² = 1.000000000` at the VE and nowhere else. `jb_k`'s `aligned_frame` pins chart
direction 0 to the symmetric-path tangent, so the test is whether the A projector fixes `e₀`.
Worst `|P_A e₀ − e₀|_g / |e₀|_g` over 180 angles × 2 mass models is **7.8e-09**; the same statistic
on a transverse chart direction is **exactly 1**. So the five transverse directions *are* the
doublet and the triplet, and "transverse curvature" means precisely those five.

**The path is a critical manifold in the transverse directions, and that is what makes a sign worth
reporting.** `dV` is an invariant covector at a point the group fixes, so it lies in the A-isotypic
component — which V0/V1 measure to be one-dimensional and equal to the path tangent. So `dV` can have
no doublet or triplet part at any angle. Measured over 180 angles × 36 combinations, with the split
taken invariantly (`Qᵀ dV` in the `g`-orthonormal block bases, not chart directions): worst
transverse/path ratio **4.14e-04**, on the worst row in the sweep, falling `4.59e-03 → 4.14e-04 →
4.60e-05` at `h = 1e-3 / 3e-4 / 1e-4` — pure `O(h²)` on a quantity that is exactly zero. The 42 rows
(of 6480) whose *path* gradient is itself below `1e-3` are the critical points, where there is no
direction to be transverse to; they are reported by absolute size instead of being dropped.

> **Two presentation corrections in that section, both earned by review.** First, `|dV|` and the
> path/transverse columns are norms of *different objects in different metrics*: `|dV|` is the
> Euclidean norm of the covector in chart coordinates (the record's convention, and what the
> icosahedron anchor `3.057691` is), while the path and transverse columns are `g`-norms of the
> gradient *vector* `g⁻¹dV`. At the icosahedron the "part" is 4.7× the "whole" for exactly that
> reason. Only the ratio of the last two is used, and both of its terms are `g`-norms, so no
> measurement moves — but two conventions were sitting in one row under one name. Second, the
> non-vacuity teeth ("the path part does not vanish") was a **max over the whole sweep**, attained
> on a `1/r^12` strut row adjacent to the pole where the gradient is `2.4e+24`. An independent
> mutation audit annihilated both gradient components on every well-behaved row and the gate stayed
> green: the teeth was *measured vacuous*, not merely weak. A second row now asserts the same teeth
> at one fixed well-behaved configuration — the icosahedron, path gradient `1.441e+01`.
> **Consequence for what the doublet and triplet blocks mean.** They are the second variation of `V`
> **normal to an invariant submanifold**, not generic off-critical curvature — `V` is stationary
> under a transverse displacement at fixed `a`, at every `a`. That is a stronger object than "the
> curvature happens to be negative here". It still does **not** make them frequencies: the *path*
> component of the gradient is emphatically not zero, so the system accelerates along the path and
> the transverse motion is driven by a time-dependent coefficient. The transverse Hessian is the
> potential term of the normal variational equation — a necessary ingredient of the dynamical
> verdict, and not the equation.

**Headline, and it is a negative result about invariance.** The transverse sign structure is **not**
invariant across the four declarations.

- At the **ground state** (`a = 0`) the sector is a transverse **valley floor** for all 36
  kernel × primitive × mass-model combinations — the recorded inertia `(6,0,0)`, recovered through a
  labelling the record did not use. *The saddle scenario the bead was opened to rule out does not
  hold at the ground state.*
- Every **inverse power** (`1/r^p`, p = 1,2,3,6,12), both primitives, both mass models — 20
  combinations — keeps **both** transverse blocks positive at every swept angle.
- …and **inside the two guard bands the main sweep excludes**, which together are 3.3% of the
  fundamental domain: all twenty, at two step sizes, at six angles across `|a − 60| < 1` (V6b) and
  ten across `a > 89` (V7 a-ter) — in each case the grid points the main sweep would have visited
  plus several closer to the pole than the grid ever gets. `a = 60` exactly is not among them and
  cannot be: the vertices merge there, every inverse power is infinite, and that is the pole itself
  rather than a gap in the coverage. That evidence did not exist when this
  section first said "valley floor throughout the fundamental domain": the coverage inside the first
  band was one combination of twenty and inside the second was one kernel of five. **The claim was
  true and the evidence for it was not in the file.** With the bands entered, the residue is stated
  exactly rather than implied: most rows in the `a = 60` band do *not* meet this file's own scalarity
  criterion at any step size the chart supports (`1/r^12` on strut midpoints a tenth of a degree from
  a pole is not a converged measurement), and at `h = 3e-4` within a hundredth of a degree of the
  pole **eight of the twenty come back negative** — every one with a per-block deviation above
  `1e+01`, and every one positive again at `h = 1e-4`. So the bands carry **signs only**, read at the
  finest step, exactly as V6 already said for the single combination it covered.
- Every **Gaussian** (s = 0.5, 1.0, 1.5, 2.5), both primitives, both mass models — 16 combinations —
  **changes sign**. The sector becomes a transverse **ridge**, and the doublet and the triplet do
  not turn over together. The located turnovers span `a = 23.2455 ± 0.0001` (gauss 2.5 / strut /
  lamina, doublet) to `a = 78.0748 ± 0.0001` (gauss 1.0 / vertex / point, triplet). There is no
  kernel-independent angle to record. The mass model alone moves a turnover by 11° (gauss 2.5 /
  vertex: doublet at `42.5641` under point masses, `31.7507` under laminae) and can change the sign
  *pattern* at one configuration (gauss 1.0 / vertex at `a = 75.262042`: `(−,+)` point, `(−,−)`
  lamina).
- These two claims gate each other. "No negative transverse curvature was found" alone would be
  satisfied by a build incapable of producing one; the Gaussian arm is what shows the sweep can see
  a negative.

> **Four decimals, and why the record previously carried six.** The bisection now runs to `1e-7`
> degrees, which is *not* the precision of the answer: V3d measures each root moving with the step
> size by up to `1.05e-04` degrees, which is the real error bar, and gates it — together with the
> slope `|dλ/da| ≥ 3.8e-01` per degree that makes a root a well-posed thing to reproduce, and the
> *absolute* per-block deviation converted into an angular uncertainty (`4.6e-06` degrees). That last
> row exists because the file's ordinary scalarity falsifier is **structurally blind here**: it
> divides by `|λ|` and a turnover is exactly where `λ → 0`. Previously the bisection stopped at
> `1e-4` while the table printed six decimals, nothing in the gate read the roots at all, and the two
> headline angles that reached this file and T2 were wrong in their fifth and sixth digits
> (`78.074799 → 78.074820`, `23.245514 → 23.245493`). Every turnover is now quoted as four decimals
> **with an explicit `± 0.0001`** — the fourth decimal is the uncertain digit, which is where a
> measurement's last quoted digit belongs. **Treat a difference below `1e-03` degrees between two
> turnovers as no difference.**

> **The guard band is applied by kernel family, and two turnovers came back from it.** The band's
> justification is an *inverse-power* pole — the twelve shared vertices merge in pairs at `a = 60`
> and `1/r^p` diverges. Gaussians do not diverge there; V6's own table has them at `dev ≈ 2e-08`
> straight through the band at every step size while the Thomson rows blow up to `4.13e-01`. Applying
> the guard by *angle* therefore suppressed measurements the file's own machinery resolves cleanly,
> and it did: two turnovers were reported as `STRADDLES a=60 / not refined`. Bisected through the
> band they are **gauss 1.0 / strut / point (triplet) at `a = 59.6881`** and **gauss 1.5 / vertex /
> point (triplet) at `a = 60.3360`**, at absolute deviations of `7.6e-06` and `2.5e-06`. One of them
> was never near the pole at all — it sits at 59.69, unlocated only because 59.5 had been struck from
> the grid.

So the kernel — already known to set the spectrum without moving the ground state (DECISION 17a) —
also sets the **sign of the topography** around the path. Choosing it is choosing whether the
symmetric sector is a valley or a ridge away from the ground state.

**What this does and does not license about the existing 1-DOF results.** `M_eff`'s 3:1 and 9:1, and
`ȧ` peaking at the vector equilibrium, are statements about the symmetric sector's own kinematics.
That sector is the fixed-point set of the chiral tetrahedral symmetry — which is the same fact V1
measures — so a trajectory started exactly on it stays on it exactly, whatever the transverse
curvature does. **Those results are correct as statements about that invariant sub-model,
unconditionally, and nothing here touches them.** What was at issue is whether the sub-model is
*representative*. At the ground state, yes, for every kernel. Away from it, it depends on the
kernel. Turning that into "the 1-DOF motion is unrepresentative" requires the dynamical calculation
above, which is not this bead and does not exist yet.

> **And the one piece of empirical evidence on that question is not a computation:
> `inviscid-qvf.11`.** The owner built this bead's 1-DOF sector in wood and wire — jitterbug units
> wired at the triangle vertices, each triangle constrained by a dowel through its centre, which
> mechanically removes five of the six internal degrees of freedom. In phase, the array **locks
> one-sidedly at the icosahedral phase**: it cannot expand, it can contract. That bead's fork is
> (A) a tension-only wire going taut, an artefact of the rig, against (B) a genuine boundary or cusp
> of the array's configuration variety. **What this file can say is narrow and is measured**: at
> `a = 22.2387561` the *single unit's* chart has local dimension 6, its momentum-free mass metric's
> smallest eigenvalue is `1/32` with **multiplicity 2** — the ordinary doublet — against
> **multiplicity 5** at `a = 60` where the doublet and triplet masses coincide exactly, and the
> smallest transverse `λ` over all 36 combinations is `+5.62e+00`. (Multiplicity, not the ratio of
> extreme eigenvalues: at `a = 60` that ratio is `5/3`, a perfectly ordinary number, because the
> metric there is degenerate but not singular.) The instrument
> is not blind — it detects both places where this variety genuinely misbehaves, the metric
> degeneracy at 60 and the dimension jump to 7 at 90 — and it sees nothing at the icosahedral phase.
> **That is evidence for fork (A) and it is not proof of it**: an array is a different variety from
> one unit, a boundary can be created by the coupling without existing in any single unit, and
> qvf.11's own precedent applies — the same instinct said `a = 60` must be singular and R3 measured
> the rank constant 36 straight through it. qvf.11's P0 (change the wire lengths, see whether the
> lock angle moves) settles this more cheaply than any computation. **Do not let a smooth single
> unit close that bead.**

**Two boundaries are measured rather than omitted, and the guards are now bounded from above.** At
`a = 60` the twelve shared vertices merge in pairs and every inverse power diverges for both
primitives; the main sweep guards `|a − 60| < 1` and V6 measures the band at three step sizes.
A gate row caps what the two guards may hide: together they remove `3.33%` of the fundamental domain
against a stated policy bound of `5%`. That row exists because both guards were measured to be
*unconstrained in that direction* — `POLE_GUARD = 3` and `BRANCH_GUARD = 5` each exited 0 with zero
red rows **and improved the reported worst scalarity**, so the gate as it stood rewarded hiding more
of the domain. There is no derivation for a particular width (a pole is a point and any positive
width is a fit), so what is asserted is that the fit stays small; it is labelled a policy bound, not
a derived one. That section also contains the run's most useful failure:
at `a = 59.99` with `h = 1e-3` the triplet comes back **negative** with a per-block scalar deviation
of **2.4** — larger than the quantity it qualifies — while every smaller step says positive. The
falsifier caught a spurious sign, and the claim it refuted ("a relative deviation of 1e-2 cannot
move a sign") was deleted rather than softened. Two gate rows now assert the fine-`h` sign *and*
that a coarse-`h` sign error must arrive with a deviation that announces it. At `a = 90` the
chart's Newton solve fails (local dimension 7; measured) and nothing is claimed there.

**And a pole the record does not have.** Approaching `a = 90`, the minimum **strut-midpoint**
separation falls linearly to zero as `2(90 − a)` in radians (ratio to that prediction: 0.999949 at
89°, 0.99999999 at 89.99°; **now asserted**, `|ratio − 1| = 5.1e-09` at 89.99, with the minimum
vertex separation over the guard band as its teeth — a mutation replacing the pole law `2(90−a)` by
`3(90−a)` previously exited 0 with every row green), while the minimum **vertex** separation stays
at 0.8165. So every
inverse power on the raw strut-midpoint primitive diverges at `a = 90` as well as at `a = 60`, and
the vertex primitive does not. The record has the `a = 60` pole for both primitives; this one
belongs to one of them — exactly the clause DECISION 17a left live, that the two primitives share
the ground state and not the landscape. It is a by-product here and deserves its own bead.

**Two other findings that were not the deliverable.** `a = 90` is a **critical point of `V`
restricted to the path, forced by symmetry** — and that qualifier is load-bearing rather than
cautious. It is the fixed point of the exact `a → 180−a` isometry, so `V` *restricted to the path*
is even about it and the path component of `dV` vanishes there. The isometry licenses only the
path-tangential derivative; V1b covers the five *chart*-transverse directions and is measured up to
`a = 89`; but the local dimension at `a = 90` is **seven**, so there is a direction the 6-D chart
never sees, and invariance of `dV` under the involution induced on a 7-D tangent cone forces only
the `(−1)`-eigenspace to vanish, not the whole covector. "`a = 90` is a critical point of `V`" with
no restriction is *not* what these rows show and is not claimed. The transverse component is
already zero at every angle (V1b), so it is the path component alone that has to collapse, and it
does — linearly, which is what an even function's derivative does at its symmetry point (Thomson /
raw vertex / point, path gradient in the `g`-norm: 2.944311 → 1.479068 → 0.296257 → 0.029628 at
`a` = 89 / 89.5 / 89.9 / 89.99, with the ratio to `(90−a)` converging 2.944311 → 2.958136 →
2.962570 → 2.962753 rather than sitting constant, which is the `O((90−a)²)` correction). **That
collapse law is now asserted in two independent ways** — its linearity (`|r(89.99) − r(89.9)| =
1.8e-04`) and its *value* against a local anchor `A90_RATE = 2.9628` (`1.6e-05` relative). It
previously fed nothing: a mutation halving the ratio's denominator, i.e. reporting a collapse twice
as slow as the measured one, exited 0 with every gate row green. The linearity row alone does not
catch that — halving the denominator halves every ratio uniformly and leaves the convergence
untouched — so the value row is the one that does. The two kernel
families disagree about the character of that point — Thomson approaches it with the path direction
negative and both transverse blocks positive, gauss 1.0 with all six negative — but the chart is 6-D
and the point is 7-D, so *no inertia count at the branch point is offered*. That is `qvf.3` /
`inviscid-yli`. And the same isometry, applied to the curvature blocks, is a real check that passes
(`1.3e-07` relative worst, against `4.1e-01` on a deliberately mismatched pair).

**What was re-measured rather than inherited, and one hole the mutation matrix found in it.** The
seed asked for one thing to be confirmed for this sweep: the base points are on the path, where `jb_u` U2c measures the section and momentum-free
forms to coincide — but `Γ` needs `dg` in all six directions and the stencil producing it steps
*transverse* to the path, so the coincidence is a question about the stencil, not the base point,
and `jb_u` measured it at two angles. Over all 180 angles × 2 mass models, `|Γ_section − Γ_horizontal|
/ |Γ| = 7.2e-09`, and the closed-form `|ZᵀWD|` on the path is `7.6e-16`. Chart invariance was also
re-run on *this file's* statistic — a chart-invariant construction with a chart-dependent labelling
would still give chart-dependent numbers: two charts agree to `5.3e-07` relative under the
momentum-free form while the section form fails the same criterion at `2.6e-01`.
> The hole: the first version of that comparison built a section geometry and compared it against
> *the sweep's own*, which is correct only while the sweep's own is the momentum-free one. A
> mutation switching the default made the section compare a thing to itself, report `0.0` as
> confirmation, and exit 0 with every gate row green — and because the two forms genuinely coincide
> on the path, no other row noticed. Both arms are now built explicitly and a row asserts, at the
> source, that the ambient form in use *differs* from the section form (`1.0e-01` relative; the
> projector removes an `O(1)` piece, so the threshold is derived). **When two things are measured to
> agree, "we used the right one" becomes unfalsifiable from the results and has to be asserted where
> the choice is made.**

The `h` sweep is reported as a **range per threshold**, not a margin at one step. `h = 3e-4` is
`jb_u` U5's measured minimum for the same chart and stencil, which is why it is used; at that step
the per-block scalarity over all 36 combinations is `1.6e-04` and it falls by `1.2e+03` across a 30×
range of `h`, against the `9.0e+02` that pure `O(h²)` predicts. The scalarity criterion of this file
does **not** hold at `h = 1e-3` — `jb_u`'s own quoted step — on the worst combination, which is why
the step size was not simply inherited.

> **The window is two step sizes wide and it fails on *both* sides.** V9 sweeps four thresholds, not
> three, and the fourth is why: an independent whole-gate `h` sweep found the gate **red at
> `h = 1e-4`, finer than the quoted step**, on V7's isometry row (`1.224e-06` against `1e-06`,
> roundoff-limited rather than truncation-limited) while V9's three original thresholds all passed
> there and V9 therefore reported `1e-4` as acceptable. With the isometry swept, V9 reproduces the
> whole-gate window independently: green at **`h = 6e-4` and `3e-4`**, red at `3e-3`, `1e-3` and
> `1e-4`. Two criteria do the excluding and they pull from opposite ends — scalarity (truncation)
> rules out the coarse steps, isometry (roundoff) rules out the finest. The count is 2 against a
> criterion of 2: **zero margin**. Two further disclosures the section now carries: the
> singlet-tangent leak is `1.855e-09` at *every* step size — flat, so the threshold annotation that
> called it `O(h²)` was refuted by the column it pointed at and has been corrected — and V9's
> scalarity column is **one combination**, the worst at `H_MAIN`, which does not stay the worst as
> `h` moves (the full-sweep worst at `1e-4` is `6.654e-04` against this probe's `1.259e-05`, a factor
> of 53).

Runtime ~220 s, exit 0, byte-identical across runs, 56 gate rows. Details, with all four declarations
on every number, in T2 `inviscid/qvf.4-transverse-curvature.md`; the independent gate audit that
found several of the holes closed above is in T2 `inviscid/qvf.4-check-suite-validation.md`.

Prior headline: `Vol_hull` is **not differentiable** at the vector equilibrium. It has a cone point —
`V/V_tet = 20 + 4√3·|a_rad| − …` along the path, and a strictly positive first-order rise in
*all six* internal directions. So `V = −k·Vol_hull` has no Hessian at the VE and cannot yield a
frequency there. Details and verdicts in T2 `inviscid/V-candidate-ground-state-measurements.md`.

Independent verification by two reviewing agents, written from scratch against the same claims —
`rv_*.py` (finite-difference Jacobian check, tolerance sweep, local-dimension estimate that
located the real branch points, labelling-free endpoint classification) and `crit_*.py`
(dynamics, collision geometry, path collisions). `DumpJb.java` dumps the real
`PhiCoordinates.Octahedrons[4]` face data from compiled Java to check model fidelity.

## Known limitations, recorded so they are not rediscovered as news

- Two checks **cannot fail** and must not be read as verification: `a_ico`'s "exactly 30 short
  edges" slices to 30 before counting, and "struts equal at every a" is a tautology of a
  rigid-body parameterisation.
- Rank constancy near a=60 is dense sampling with σ bounded away from zero, not a proof.
- Reachability is one constructed path per target, not a component decomposition of the variety.
- `continue_along`'s arclength is a binary survived/died signal, not a magnitude — it caps by
  construction.
- Sign convention: these scripts use `σ = +(sx·sy·sz)`; the Java uses `σ = −(sx·sy·sz)`. The two
  families are related by `a → −a` exactly (set-Hausdorff 4e-16). Immaterial — do not "fix" it.
- A **pole manufactures a false stationary point** in any dV/da sign-change scan: an inverse power
  runs +∞ just below a=60 and −∞ just above, and bisection refines that textbook sign change to
  a=60, where V is not finite. `jb_t.path_critical` verifies every root (V finite, |dV/da| actually
  small) and prints the rejections rather than dropping them.
- Off a critical point the Hessian is **not** a tensor, so "the spectrum at the icosahedron" is
  chart-dependent — measured at 45% of the spectrum span between two charts that agree to 6e-8 at
  the VE. **The chart-dependence is resolved by `jb_u`** (Riemannian Hessian of the momentum-free
  mass metric; worst chart disagreement 3.138e-6 relative over 36 combinations). What is *not*
  resolved, and is not claimed to be: the resulting eigenvalues off a critical point are curvature
  scales, not frequencies, and no spectrum along the path was produced — so "how the curvature scales
  vary along the motion" is **enabled but not delivered**. Four follow-on limits replace the old
  bullet: (i) the *section* mass metric is not chart-invariant off a critical point and must not be
  used there — use `HorizontalForm`, and note that Java `InternalMassMetric` and `jb_r.metric()` are
  section metrics; (ii) the centroid-pivot gauge is momentum-free only ON the symmetric path, so any
  transverse work (`inviscid-qvf.4`) must project the rigid orbit out explicitly; (iii) the `d_c W`
  term in `dg` is invisible on the path and real off it (`jb_u` U4d gates it there); (iv) no
  measurement discriminates among equivariant projections along span(Z).
- A LINEAR basis remix **cannot test a Christoffel term** — a linear map has zero second derivative.
  `jb_r` R5 and `jb_s` S1(b) are covariance checks, not connection checks. Measured teeth ratio of
  the linear test: 1.03. Of the `Dφ = I` nonlinear test: 6.73e+05.
- Per-block mass factors `√(8/5), √2, 2` between the point and lamina models are **VE-specific**
  (`jb_u` U2b ii-ter). Do not carry them off a critical point.
- The mass metric's eigenspaces are **not** a usable block labelling near `a = 60`: doublet and
  triplet masses are equal *exactly* at the octahedron (`jb_v` V2, closed form). `jb_u`'s
  `blocks_by_irrep` refuses there, correctly. Use `jb_v`'s character projectors.
- `jb_v` sweeps `[0, 90)` and guards **two** neighbourhoods, `|a − 60| < 1` and `a > 89`. Inside
  them a finite-difference second derivative stops resolving, and inside the first one a coarse
  step produces a **spurious sign** (measured). The guards are priced from the run's own numbers,
  not assumed; nothing inside them is quoted as a ratio.
- `jb_v` re-gates only what it adds. The correctness of the Riemannian Hessian itself rests on
  `jb_u`'s 17-row gate and is not re-derived.
- `jb_v` gives **no dynamical stability verdict**. No equation of motion is integrated, no Floquet
  analysis exists, and a negative transverse *curvature* at a point the system passes through at
  speed is not an unstable *mode*. That calculation needs its own bead.
- `numpy` 1.26 on this BLAS emits spurious `RuntimeWarning: divide by zero / overflow / invalid
  encountered in matmul`, reproducibly for `np.zeros((72,48)) @ np.zeros(48)`. Verified spurious;
  `jb_r`…`jb_t` silence it with `np.seterr`. Do not read it as a numerical event.

Requires `numpy` and `scipy`. Run any file directly: `python3 -W ignore jb_a_family.py`.
