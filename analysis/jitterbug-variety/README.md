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
