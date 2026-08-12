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

Headline: `Vol_hull` is **not differentiable** at the vector equilibrium. It has a cone point —
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

Requires `numpy` and `scipy`. Run any file directly: `python3 -W ignore jb_a_family.py`.
