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
  the VE. The canonical fix (the Riemannian Hessian of the mass metric) is not implemented.
- `numpy` 1.26 on this BLAS emits spurious `RuntimeWarning: divide by zero / overflow / invalid
  encountered in matmul`, reproducibly for `np.zeros((72,48)) @ np.zeros(48)`. Verified spurious;
  `jb_r`…`jb_t` silence it with `np.seterr`. Do not read it as a numerical event.

Requires `numpy` and `scipy`. Run any file directly: `python3 -W ignore jb_a_family.py`.
