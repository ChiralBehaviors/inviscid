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
- `numpy` 1.26 on this BLAS emits spurious `RuntimeWarning: divide by zero / overflow / invalid
  encountered in matmul`, reproducibly for `np.zeros((72,48)) @ np.zeros(48)`. Verified spurious;
  `jb_r`…`jb_t` silence it with `np.seterr`. Do not read it as a numerical event.

Requires `numpy` and `scipy`. Run any file directly: `python3 -W ignore jb_a_family.py`.
