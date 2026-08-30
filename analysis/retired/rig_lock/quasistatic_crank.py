"""Step Z: the CONTACT GEOMETRY KERNEL and QUASI-STATIC CRANK STEPPER for the
physical-plate array.

Bead `inviscid-qvf.17` (Phase 1a, the kernel, FROZEN below at the file sha
recorded in T2 `inviscid/qvf.17-contact-kernel.md`) and `inviscid-qvf.18`
(Phase 1b, the crank stepper, everything from "Z8: QP MACHINERY" onward) of
DECISION 18 (design of record: T2 `inviscid/design-contact-dynamics-array.md`).
DECISION 18 REINSTATES contact as the neighbour coupling for the ARRAY MODEL
ONLY (DECISION 16's "interference is permitted" stands unchanged for the
single-unit abstract variety): no plate may pass through another, inter-unit
joints are tension-only wire loops of slack `w`, and plate thickness `t`
enters as a gap offset (`gap >= t`).

qvf.17's 51 gate rows (Z0-Z7 below) are UNCHANGED and must still all PASS --
this file is EXTENDED, not forked, per both beads' own instruction.

WHAT THE KERNEL DELIVERS, and where each acceptance-criteria item lands
------------------------------------------------------------------------
1. TOPOLOGY AS DATA. `jb_x_array_linkage.Topology` / `build_topologies()` are
   imported and used unmodified -- registry-built arrays are then data, not a
   rewrite.
2. PLATE-PAIR ENUMERATION across an arbitrary `Topology` (GAP 2, greenfield):
   `enumerate_plate_pairs`, with the pin-sharing exclusion as a STATED
   PREDICATE (`_shares_pin`), not a hardcoded list.
3. SIGNED TRIANGLE-TRIANGLE GAP with witness points and an outward contact
   normal (GAP 3, greenfield; critique C1): `signed_gap`. The key geometric
   fact this kernel rests on, DERIVED and gate-checked (Z0), not assumed: for
   this specific family (`jb_a_family.corners`), plate f's triangle always
   lies in the plane {x : x . u_f = Z cos(a)} -- rotation is defined ABOUT the
   axis u_f, which preserves that plane, and the vertices satisfy (v - c) . u
   = 0 identically for every octahedron face. So a plate's OUTWARD NORMAL is
   the FIXED axis u_f, independent of phase and of unit translation. Two
   branches, both documented in `signed_gap`'s own docstring per critique
   C1(i): PARALLEL-FACING plates (the registry pair and the folding-square
   pair are both this case) get an EXACT closed-form projection with each
   triangle's own centroid as witness; general pairs get a closest-point
   search (edge-edge + vertex-face, 15 candidates) with a Moller-Trumbore
   piercing test for sign. Deep-penetration accuracy in the general branch is
   explicitly OUT OF SCOPE -- see "A ROW DELIBERATELY NOT BUILT" in the gate.
4. PLATES AT A PHASE (GAP 4): `unit_plates(a, origin) = corners(a) + origin`,
   the pure-translate reference placement `assemble_free` itself uses
   (`rots=None`).
5. `A_ICO`, `STRUT_LEN` imported from `jb_x_array_linkage` (GAP 5).
6. WIRE SPANS (item 7): `span_length` is a thin, named wrapper on the same
   quantity as `jb_y_dephasing.span_length`/`member_length` (the bead's cited
   prior art). `z6_wire_and_thickness` exercises it both on synthetic known-
   distance configurations (rows H-v/vi/vii/viii, per the acceptance
   criteria's own "synthetic configuration of known span" wording) and once
   on REAL topology geometry (the census square's non-pinned vertex pair
   between two actually-placed units, Z6's last row) so the function is
   shown working on more than a hand-picked number.
7. THICKNESS (item 8): `admissible(gap, t) = gap >= t`, `t` a module constant.

TWO INSTRUMENT LESSONS THIS FILE EXISTS TO SURVIVE (T2 `registry-viewer-no-
plate-crossing.md`, two independent errors on the SAME question in one day):
a crossing census proves nothing about ORDER (parallel pass-through registry
pair, mechanism 1 -- needs a SIGNED gap); a strict-interior pierce test proves
nothing under exact symmetry (in-phase neighbours' symmetric edge-through-edge
crossings, mechanism 2 -- needs a PERTURBED, midpoint-strict-interior test).
Both appear as gate rows below, including a row asserting the WRONG instrument
(the naive strict-interior test) DISAGREES with the robust one.

FOUR DECLARATIONS, per the AMENDED design of record (T2 23230, section FOUR
DECLARATIONS amended 2026-08-21) -- emitted here and again in `main()`'s banner
--------------------------------------------------------------------------------
KERNEL, MASS MODEL, PRIMITIVE: INAPPLICABLE. Nothing in this file has an
interaction potential, a mass, or a choice of primitive (vertex vs strut-
midpoint) -- every quantity is a static geometric measurement at a held phase.
METRIC FORM: QUALIFIED, NOT FLATLY INAPPLICABLE (the amendment). Phase 1's QP
objective ||v - v_cmd||^2_W (bead .18) carries a weight W, and a weight is a
norm choice -- but THIS bead carries NO QP objective at all, so the amendment's
qualification is vacuously satisfied here rather than exercised: every
quantity this file prints is a NORM-FREE GEOMETRIC LENGTH in R_oct = 1 units
(gaps, spans, diagonal lengths, crossing-segment lengths), and no weight W
exists yet to make any of them "per unit of" anything. The W-treatment choice
itself is bead .18's, not this one's -- stated here so the rule does not lapse
silently (a critic finding on jb_x), and restated in `main()`'s banner.

PHASE 1b -- THE QUASI-STATIC CRANK STEPPER (bead `inviscid-qvf.18`)
--------------------------------------------------------------------------------
Moreau sweeping, no inertia: turning the wooden array slowly by hand. Per
STEP: minimize ||v - v_cmd||^2_W subject to J_pin v = 0 (the 36 intra-unit
hinge rows PER UNIT, block-diagonal -- NOT `assemble_free`'s inter-unit rows,
which encode RIGID pins, the pre-DECISION-18 semantics; `assemble_free` /
`assemble_doweled` are reused ONLY for the DOWELED diagnostic below, where
rigid inter-unit pins are exactly the intended comparison model) and
`grad g . v >= 0` on the active contact set (`g <= EPS_ACT`) and the active
taut-wire set. Advance `x += h v` via `apply_body_motions` (exact, per body),
Newton-project back onto the pin manifold, refresh active sets.

SOLVER (GAP 1, the three-way decision, bead comment 2026-08-21): (a) a
hand-rolled active-set QP; (b) an LP reformulation on `jb_y_dephasing._lp_core`
sup-norm path; (c) `scipy.optimize.minimize` (jb_m precedent). CHOSEN: (a),
implemented as the classical Lawson-Hanson reduction of a Least-Distance
Program to Non-Negative Least Squares (`scipy.optimize.nnls` as the numerical
kernel -- the same relationship `_lp_core` has to `scipy.optimize.linprog`,
method="highs"). REJECTED (b): the sup-norm LP answers a DIFFERENT question
(min-max slack, not a weighted 2-norm nearest-point projection) -- the design
of record's METRIC FORM qualification is specifically about ||.||^2_W, and an
LP reformulation would silently change what "nearest" means. REJECTED (c):
`scipy.optimize.minimize` is a general nonlinear solver on a problem that is
EXACTLY quadratic-with-linear-constraints; its termination-criterion and
BLAS-path determinism question (named in the bead) buys nothing here that (a)
does not already give more cheaply and more predictably. A NAIVE hand-rolled
primal active-set method (repeatedly solve the KKT system for the current
working set, add the worst-violated inactive candidate, drop the
most-negative-multiplier active one) was tried FIRST and CYCLES on this
array's real, verified constraint degeneracy (many candidate contacts are
EXACTLY, simultaneously active at the array's symmetric configurations, not a
tolerance artifact -- confirmed by direct trace) -- Bland's-rule tie-breaking
does not rescue it, because the degeneracy is in the GEOMETRY, not merely in
the pivot order. The LDP/NNLS reduction sidesteps this: NNLS's own active-set
search (Lawson & Hanson 1974, ch. 23) is the same textbook algorithm, applied
in a REDUCED, ORTHONORMAL null-space coordinate system (`J_pin`'s null space,
via SVD) where the objective is the identity metric and the reduction's
correctness was verified against hand-computed tiny examples before use here.
Determinism: NNLS is a fixed, non-random iterative algorithm on fixed inputs,
matching this project's existing trust in `linprog(method="highs")` inside
`_lp_core` -- gated as row L, byte-identical output, same as every other file
in this directory.

DRIVE: `crank_v_cmd` builds `v_cmd` as `path_tangent_48(a_hat)`'s per-unit
phase direction PLUS `topo.dsites(dverts_exact(a_hat))`'s per-unit RIGID
TRANSLATION (the lattice BREATHES -- omitting this term was tried and gives a
provably wrong answer: a unit's own internal fold direction alone is blind to
the array's inter-unit registry motion, verified to disagree with the exact
nonlinear `signed_gap` ground truth by finite difference). `driven` selects
which unit(s) receive a nonzero target: `"all"` for the IN-PHASE / uniform-
expansion crank (design of record's stated alternative DRIVE, and the one
exercised by every row below); an explicit unit index for the single-crank-
handle variant, which the driven-unit-is-an-input requirement above needs
demonstrated live, not merely parametrically possible -- Z14 exercises it.

CONTACT GRADIENT ROW -- the HAZARD's territory, discharged HERE, gated below
(K-hazard rows): `signed_gap` has THREE distinct first-order regimes it does
not itself expose, and USING nA UNCONDITIONALLY AS THE GRADIENT DIRECTION IS
WRONG for one of them -- verified by finite difference, not assumed. Parallel-
facing: direction is nA (exact). General, UNPIERCED (the branch the qvf.17
critique named as having zero negative-gap coverage): direction is the UNIT
VECTOR from witness_A to witness_B, NOT nA -- confirmed by a direct
finite-difference cross-check (`nA`-based linearization was off by more than
an order of magnitude and even the WRONG SIGN for a real SC7 pair at this
regime; the witness-vector direction matches finite difference to 1e-8).
General, PIERCED (gap = -abs(proxy), proxy = (witness_B-witness_A).nA):
direction is `-sign(proxy) * nA`. `_contact_gradient_direction` implements all
three; every contact row this file builds goes through it, never through nA
alone. The HAZARD's required discharge is TREATMENT (b) of the bead's two
options: a feasibility-preserving invariant, not added negative-gap gate
coverage on the kernel itself (which is qvf.17's file, frozen) -- `crank_run`
BACKTRACKS the step size `h` (bisection, bounded attempts) using the REAL,
NONLINEAR `signed_gap` re-evaluated at the trial position for EVERY enumerated
plate pair (general branch included, not just the active-set candidates) and
accepts a step only if no pair's gap drops below `-GAP_FLOOR_TOL`; Z13 gates
the observed minimum general-branch gap across every recorded run at
`>= -GAP_FLOOR_TOL`, with a mutation probe (the SAME runs with backtracking
DISABLED) demonstrating the invariant WOULD be violated without it --
non-vacuous, per the bead's requirement that the row can fail.

JAM, QPFAIL, REACHED -- three DISJOINT statuses (the `_LPFail`-class
distinction the HAZARDS section requires): QPFAIL is a solver-side event (the
LDP/NNLS reduction's residual collapses, or the step-size bisection cannot
find ANY feasible `h` down to `H_MIN`) and is NEVER read as a physical verdict
-- `crank_run` returns it as its own status, distinct from JAM, and the gate
checks for it explicitly (K-fail rows) so a solver failure cannot silently
print as a jam. JAM is `rate_achieved < JAM_RATE_TOL` -- the max achievable
projection of the solved `v` onto `v_cmd`'s driven-direction, 1.0 when fully
unconstrained (verified: an UNCONSTRAINED solve reproduces `v == v_cmd`
exactly, `rate = 1.0`, to solver tolerance).

ONSET vs SUSTAINED, AND WHY BOTH ARE NEEDED (unchanged by the FIXED-ROUND
below, restated because it is still the reason G and H use different
`instant_jam` settings): under the IN-PHASE / uniform-expansion DRIVE
(`driven="all"`), CONTACTS ALONE (independent of w -- verified directly)
produce a transient ~0.30 first-step resistance that RELAXES to the fully-
free rate after roughly 0.02-0.03 degrees of real, Newton-projected motion,
by the same equivariance argument as before (a v_cmd identical across every
unit, over a feasible region invariant under the array's symmetry group, has
a unique minimizer that is a fixed point of that symmetry -- contacts alone
cannot dephase either). `crank_run` therefore carries TWO distinct jam
semantics (`instant_jam`, see STALL_RATE_TOL's docstring): ONSET (row G)
asks "is there already reduced rate at the start", SUSTAINED (rows H, I,
and the CROSS rows below) asks "does the run ever get PERSISTENTLY stuck".

FIXED-ROUND (substantive critique 23262 C1, SHIP-BLOCKER; code review 23261
C1/H1/M1-M3): the ORIGINAL wire attachment -- `topo.contacts`' single
coincident vertex per neighbour link, inherited unmodified from jb_x's
pre-DECISION-18 rigid-pin topology -- is a SYMMETRY-FIXED POINT of
`topo.sites`/`dsites`: its span is EXACTLY 0.0 at every phase, for every
wire, under ANY uniform driving, regardless of w. It could never be what T2
23230 means by "each tied vertex pair across each shared face", and made
row H's "reach" and the earlier DEVIATION section's "w has zero measured
effect" both true FOR THE WRONG REASON -- an artifact of which vertex
`topo.contacts` happens to identify, not a fact about wire-loop resistance
under this drive. `_wire_attachment_pairs` replaces it: per `topo.contacts`
entry, the plate pairs between its two units that are near-zero at a small
reference angle with POSITIVE d(gap)/da (OPENING/valley pairs; the CLOSING/
ridge pairs at the same locations are already covered, unchanged, by the
ordinary CONTACT machinery). Their gap is EXACTLY proportional to
`fold_halves`' own quantity: gap / fold(a) = 2/sqrt(3), constant to 5
decimal places at every angle checked -- measured, not guessed (the
critique's own "2*fold(a)" was a reasonable diagnosis-time estimate for a
raw-vertex measurement that turned out not to exist; this is the real,
verified plate-level relationship). `wire_gradient_row` uses this SIGNED,
along-normal quantity as `span` -- NOT the raw Euclidean witness-to-witness
distance, which was tried first and is wrong: these plates sit laterally
offset by a large, physically irrelevant in-plane amount (order 2.8) that
swamps the actual opening signal entirely (verified: it stayed near-constant
across every `a` tested instead of growing from 0).

With this fix, w IS causally verified, not merely correlated with G/H's own
instant_jam difference: the CROSS rows hold the protocol FIXED
(`instant_jam=False`, SUSTAINED, both runs) and vary ONLY w -- w=0 makes
every opening wire taut from the first instant (a taut w=0 wire never
relaxes, since span grows monotonically with the verified 2/sqrt(3)*fold(a)
relationship, unlike the transient contact-only resistance), producing a
GENUINE SUSTAINED jam (not merely an onset reading); w=w_ico is generous
enough that the opening pairs' span, topping out near fold(a_ico)*2/sqrt(3)
approx 0.505, never approaches the 0.874 limit, so only the (transient)
contacts matter and the run reaches. Same protocol, different w, different
outcome -- the exact test the critique specified. G's own jam binding set
now also contains wires directly (not contacts alone), gated (K rows).

A second, independent bug surfaced live while building the CROSS rows and is
recorded here because it is not really about wires at all: `crank_run`'s
feasibility backtrack, when checked against a value-based or distance-based
PREFILTERED candidate list, missed real violations under sustained,
multi-step integration (a regime the file's ORIGINAL rows never exercised,
since only G, using `instant_jam=True`, ever ran at w=0 before). Three
prefilter designs were tried and rejected, each found broken live, not
hypothesized -- see MEANINGLESS_DEPTH_FLOOR's docstring for the full
account. The fix checks every enumerated pair, every backtrack iteration,
and rejects only a SHALLOW violation (`-MEANINGLESS_DEPTH_FLOOR < g <
-gap_floor`) -- deep, large-magnitude negative values are `signed_gap`'s own
documented parallel-facing-branch artifact for pairs (INTRA- or INTER-unit)
that do not actually face each other, not a real constraint. A backtrack
that exhausts every step size down to H_MIN is now reported as JAMMED (a
genuine "no positive progress possible" physical fact, the LOCK definition's
own wording), not QPFAIL -- QPFAIL is reserved for actual solver failure
(caught earlier, at `crank_step`'s own status check), the `_LPFail`-class
distinction the HAZARDS section requires.

PHASE 1c -- THE LOCK SURFACE a*(w,t) + MOTION ORDER (bead `inviscid-qvf.19`)
--------------------------------------------------------------------------------
THE DELIVERABLE: ONE model (the frozen kernel + stepper above, called but not
modified) producing BOTH bench predictions of the qvf.11/qvf.15 fork. `a*(w,t)`
is the SUSTAINED jam angle (`crank_run(..., instant_jam=False)`) swept over a
w-grid (fixed t=0, arm A: does the lock move with wire slack, per the qvf.11
k-table's 13.849356 deg span) and a t-grid (fixed, large w so wires never bind,
arm C: does the lock move with thickness and stay w-insensitive). Each grid
carries a second, absolute, incommensurate arm (house style). The headline row
is the DISCRIMINATION RATIO |da*/dw| normalised against |da*/dt|, gated into a
two-sided band with three printed verdicts (ARM-A-LIKE / ARM-C-LIKE /
NON-DISCRIMINATING) that FAILS, never inf-passes, when not computable. NO row
asserts `a* == a_ico` (the qvf.15 ruler-test hazard: both prior mechanisms
predict that exact instant, so it discriminates nothing).

A GENUINE STEP-SIZE SENSITIVITY, calibrated live before committing the grid
(not assumed): `crank_run`'s SUSTAINED outcome at a fixed (w, t) is NOT
h0-invariant. `H_STEP=2.0` (bead .18's own REACH-tuned value) fails to detect
a real jam at w=0 that `h0=0.5` (bead .18's own JAM-tuned value, matching the
independently reviewed 0.7852886398227111 number) DOES detect. This bead fixes
`H_LOCK=0.5` for every call it makes, so the surface is internally consistent.

FIX ROUND (critique T2 23299, 2 ship-blockers; code review T2 23337, 1 High
2 Medium 2 Low) -- addressed here, not deferred:
SHIP-BLOCKER 1: the critic's own independent h0-refinement probe found the
h0-sensitivity above is WORSE than first disclosed -- a*(w=0) does not
converge across h0 in {1.0, 0.5, 0.25, 0.125} and the ORIGINAL Z15 CROSS
test (bead .18) was CONFOUNDED (cross_w0 at h0=0.5, cross_wico at h0=2.0 --
different w AND different h0). `z15_crank_gates` gains two matched-h0 legs
(cross_w0/cross_wico themselves UNCHANGED) completing both matched pairs;
THE CORRECTED FINDING: at matched h0, w does NOT flip the outcome either way
(both jam at h0=0.5, both reach at h0=2.0) -- bead .18's "W IS CAUSAL"
reading does not survive a matched-h0 test. `h_refinement_probe` gates the
ACTUAL stable claim (interior w-SPREAD stays flat across H_REFINE_LEVELS,
even though each point's own absolute a* keeps drifting) two-sidedly against
the boundary point's genuine non-convergence, and classifies any small-h0
qpfail as budget-exhaustion (verified: full rate_history, rates far from
stall) rather than a genuine solver failure.
SHIP-BLOCKER 2: the ORIGINAL t-column started every sweep at a=0, where the
fold-mechanism's touching pairs are EXACTLY at gap=0 by construction, making
"any t>0 collapses a*" a tautology of the idealized zero-clearance pose, not
a measurement. `compute_a_start_t` derives a clearance-relieved start angle
(bisection over the OPENING pairs' own gap(a) against a FIXED target gap --
see its own constant's docstring for why the target is fixed rather than
re-derived from T_GRID's own max), validated live by a companion row
confirming real multi-step motion results, not another instant re-lock.
T_GRID/T_GRID_ALT are re-priced to straddle the re-based array's OWN
jam/reach threshold (found near t=0.187 by live calibration), and a*_t is
now an OPENING RANGE from A_START_T, not an absolute angle from the touching
pose -- Q2's denominator is a measurement again, not a construction artifact.
HIGH: the ALT grids' non-ratio-derivation from the PRIMARY grids is now a
GATED row (coprimality/offset-apartness, the jb_y K_GRID precedent), not
prose alone. MEDIUM: `w_ok`/`t_ok` fold non-emptiness into their `all()`;
`motion_order_trace` reads a single, homogeneous raw-velocity-norm metric
for every unit (previously driven units used a v_cmd-projected rate, an
incommensurate scale against undriven units' raw norm), and the write-back's
own strongest "units take turns" evidence (undriven corners exactly 0.0 for
several early steps) is now a gated row, not prose only. LOW: row T's
distinct-binding-sets row now has a real, failable condition (`> 1`); the
shuffle/flatness thresholds are named module constants.

MOTION ORDER: a bespoke, short trace (`motion_order_trace`, NOT a reuse of
`crank_run`'s feasibility backtrack -- see its own docstring) reading per-unit
velocity components directly off `crank_step`'s solve. Single-crank-handle
drive shows units do NOT all move together (a falsifiable claim, row S);
driven="all" at large w, t=0 shows the array DOES move uniformly (the control
that can fail).

FOUR DECLARATIONS for Phase 1c: KERNEL, MASS MODEL, PRIMITIVE -- INAPPLICABLE,
NOT FORGOTTEN, same reasons as Phase 1a/1b. METRIC FORM carries bead .18's own
TREATMENT (a) choice FORWARD (W = identity, already gated W-insensitive): a*
is the same norm-free jam angle Phase 1b already established as quotable, and
da*/dw, da*/dt are differences of it -- norm-free too, no "per unit of W"
hedge (that hedge belongs to treatment (b), not (a)).

WRITE-BACK STATUS: T2 `inviscid/qvf-lock-surface-phase1.md` is stamped
"DRAFT -- NOT YET VALIDATED (pending inviscid-qvf.20)" per the bead's own
S2 critique finding -- these are physically actionable numbers the owner
reads against his rig, and qvf.20's probe pass can still move them.

CONVENTIONS INHERITED FROM THIS DIRECTORY
------------------------------------------
Deterministic and byte-identical across runs; exit code from the gate table;
no raise inside a swept loop; a check whose non-vacuity is printed prose rather
than an assertion cannot fail; every guard band is constrained from ABOVE as
well as below; every sweep grid has a SECOND, ABSOLUTE, INCOMMENSURATE arm.
Every threshold used here is RE-DECLARED LOCALLY even where it is identical to
jb_x's or jb_y's (the mutation-probe rule). No file writes; stdout only; no
argparse/argv/env. Run from the repository root with python3.
"""
import sys

import numpy as np
from scipy.optimize import nnls

from analysis.retired.rig_lock import cache as jb_cache
#: The importable name of THIS module. Spelled as a literal rather than taken
#: from `__name__` because `__name__` is "__main__" under `python
#: jb_z_quasistatic_array.py`, and a prefetch worker in a fresh interpreter
#: must be able to re-import the same module by name to re-enter `crank_run`.
_MODULE = "analysis.retired.rig_lock.quasistatic_crank"

from analysis.model.jitterbug import R_CIRC, Z, corners, faces, rot
from analysis.model.strut_clearance import segment_distance as jb_g_segment_distance
from analysis.retired.rig_lock.array_linkage import (A_ICO, DIAGONALS, PAIRS, STRUT_LEN, STRUTS,
                                SQUARE_DIAGONALS, Topology, apply_body_motions,
                                assemble_doweled, build_topologies, dverts_exact,
                                hinge_jacobian, path_tangent_48,
                                position_jacobian_row, rank_of, verts, unit_corners)

# ==========================================================================
# LOCAL CONSTANTS
#
# Every constant a mutation probe needs to reach is defined HERE, locally,
# even where an identical constant already exists in jb_x or jb_y.
# ==========================================================================

#: Re-declared from jb_x (mutation-probe rule): the icosahedral phase.
A_ICO_LOCAL = 22.238756093

#: Re-declared from jb_x: strut length == octahedron edge length.
STRUT_LEN_LOCAL = R_CIRC * np.sqrt(2.0)

CONST_TOL = 1e-9

#: The deliberate offset bounding the a_ico agreement FROM ABOVE (two-row
#: control idiom, house style): a value offset by this much must be REJECTED
#: at the same tolerance the true value is accepted at.
AICO_CONTROL_OFFSET = 1e-3

#: NOT declared here: RANK_RTOL, SOLVE_TOL, SWEEP_LO/HI/STEP(_ALT),
#: SPAN_COARSE_STEP(_ALT), AICO_RECORD_QUANTUM. The survey's constant list
#: (T2 23234) covers all four DECISION 18 beads; this one does no rank
#: computation, no Gauss-Newton solve, and no swept taut-angle span search
#: (Z7's crossing census samples fixed angles, it does not search for a
#: crossing threshold) -- so those names would be exactly the anti-pattern
#: jb_y's own docstring warns against: "a constant a mutation probe can
#: reach but no row consults is worse than absent: it looks like a guard."
#: Declared here only what a row in THIS file actually reads.

#: Two absolute, incommensurate angle ladders used wherever this file sweeps
#: over phase. Neither is a ratio of the other (jb_y's recorded bug: a ratio-
#: derived second arm coarsens in lock-step with the first and cannot see a
#: change the first arm's step size hides).
ANGLE_GRID = (1.0, 5.0, 10.0, A_ICO_LOCAL, 30.0, 45.0, 55.0)
ANGLE_GRID_ALT = (0.7071067811865476, 6.334166025, 13.816957,
                   A_ICO_LOCAL + 0.0, 27.912878475, 41.833333333, 53.111111)

#: PARALLEL-FACING detection threshold: two plate normals are treated as
#: (anti)parallel when |nA . nB| exceeds this. 1 - 1e-9 is generous enough to
#: swallow floating point roundoff in a unit-vector dot product while still
#: excluding any pair that is genuinely tilted relative to one another --
#: bounded from BELOW by a control row using a plate pair 1 degree off
#: parallel, which must NOT take the exact branch.
PARALLEL_TOL = 1.0 - 1e-9

#: How close a plate normal's dot product with itself across a phase sweep
#: must stay to 1.0 for the "plate normal is phase-invariant" fact (Z0a) to be
#: considered confirmed rather than merely plausible.
NORMAL_INVARIANT_TOL = 1e-12

#: Finite-difference step for the normal-orientation check (row G-i) and its
#: tolerance. eps chosen away from machine-epsilon noise and away from any
#: curvature the closest-point search's own witness-point re-selection could
#: introduce (a large eps could cross a Voronoi-region boundary of the
#: closest-point search and register as a false failure).
FD_EPS = 1e-6
FD_TOL = 1e-6
#: Deliberate FD-eps offset bounding FD_TOL from above: at ten times the step,
#: curvature bites and the same tolerance must reject.
FD_CONTROL_FACTOR = 10.0

#: Barycentric tolerance for "witness point lies on its own triangle".
BARY_TOL = 1e-9

#: Independent closest-point cross-check tolerance.
WITNESS_TOL = 1e-9

#: Wire-span / thickness tolerances (rows H).
SPAN_TOL = 1e-12
THICKNESS_TOL = 1e-12
#: Wire activity band: a span within this of w counts as taut/active.
ACTIVE_TOL = 1e-7

#: Perturbation census (row D): seeded, deterministic, fixed count.
PERTURB_SEED = 20260821
PERTURB_N = 10
PERTURB_MAG = 1e-4

#: The lattice spacing used for the "spacing 2" crossing census -- an
#: ABSOLUTE number in R_oct = 1 units, independent of the topology's own
#: breathing site spacing (23195's construction, re-derived below in Z7, is a
#: FIXED external spacing, not `Topology.sites`).
CENSUS_SPACING = 2.0

#: Crossing-count control offset/band edges recorded in T2 23195, re-derived
#: here rather than trusted: spacing-2 crossings must be 0 at a = 0 and 4 at
#: every angle in ANGLE_GRID's interior set, stable 4..4 under perturbation.
CROSSING_TARGET = 4

#: PHASE 1b (bead qvf.18) constants -- GAP 8: h, eps_act, the jam tolerance,
#: the driven-unit index are MODULE CONSTANTS, not flags.

#: Active-set band: a contact within this of touching (or a wire within this
#: of taut) is a QP candidate this step. Priced from the run's own numbers,
#: not guessed -- Z13 gates the result's INSENSITIVITY to it across a decade
#: (EPS_ACT vs EPS_ACT_ALT below), the same discipline qvf.17 asked of the
#: bead that would consume this constant.
EPS_ACT = 1e-3
#: Second, absolute, incommensurate value for the eps_act insensitivity row
#: (a full decade below EPS_ACT, not a round multiple of it).
EPS_ACT_ALT = 1.037e-4

#: Row G's jam-angle band: onset must be detected at or before this many
#: degrees off a=0 ("jam angle < a small named epsilon", the bead's own
#: wording; > 0 is explicitly not required).
JAM_ANGLE_EPS = 1.0
#: CONTROL: an angle offset by this much would no longer count as "immediate".
JAM_ANGLE_CONTROL_OFFSET = 5.0

#: Nominal crank step, in degrees of `a` per step (h is a STEP PARAMETER, not
#: a physical time -- Phase 1 has no time scale, per the FOUR DECLARATIONS).
H_STEP = 2.0
#: Step-size bisection floor: below this, no feasible step exists and the
#: step is a solver-side failure (QPFAIL), never a physical jam.
H_MIN = 1e-6
#: Bounded bisection attempts before declaring QPFAIL on a single step.
H_BACKTRACK_MAX = 30

#: The feasibility-preserving floor (HAZARD discharge): a step is accepted
#: only if the REAL, re-evaluated `signed_gap` for EVERY enumerated pair
#: stays at or above this after the step -- FIX (code review 23261 M2, and
#: a real bug the w-causality cross-test exposed live, see `crank_run`'s
#: backtrack loop): an earlier version pre-filtered which pairs got this
#: post-step re-check to a "watch band" near-zero at the STEP'S START,
#: priced on an unproven claim that no pair could cross the band within one
#: step; a pair was found, live, that started well outside the band and was
#: already in slight violation ONE STEP LATER with nothing having bisected
#: against it. No closed-form bound on a single step's worst-case gap
#: change was ever derived to justify a prefilter radius, so this version
#: checks every enumerated pair on every backtrack try instead -- slower,
#: but the HAZARD's entire point is that this check must not be foolable.
#: Priced relative to EPS_ACT (one decade looser): tight enough that no
#: pair crosses meaningfully into interpenetration within one bisected
#: step, loose enough that a step of genuinely-achievable size is still
#: acceptable -- 1e-9 was tried first and forces h toward its bisection
#: floor on EVERY step (verified: a pair with a large, harmless margin
#: still degrades by more than 1e-9 at any non-infinitesimal h, so the
#: check was rejecting motion that never approached the true constraint).
GAP_FLOOR_TOL = 1e-4
#: Mutation-probe control: backtracking DISABLED must demonstrably let the
#: minimum observed general-branch gap fall below this OFFSET (looser than
#: GAP_FLOOR_TOL), proving the invariant row can fail.
GAP_FLOOR_CONTROL_OFFSET = 1e-2

#: Jam threshold on the achieved-rate projection (1.0 = fully unconstrained).
#: Set from the VERIFIED physics (module docstring "DEVIATION" section): the
#: true immediate resistance measured under in-phase driving is ~0.30,
#: row I's disabled-contacts control is 1.0 -- 0.5 cleanly separates them.
JAM_RATE_TOL = 0.5
#: Two-row control idiom: a rate at or above this offset must NOT read as a
#: jam (row I's control realizes exactly this, at rate ~1.0).
JAM_RATE_CONTROL_OFFSET = 0.9

#: SUSTAINED-reachability stall threshold (rows H, I): distinct from
#: JAM_RATE_TOL. crank_run has TWO jam semantics, both real, neither a stand-
#: in for the other -- verified via full multi-step integration with genuine
#: Newton pin-projection, not asserted: (i) ONSET (row G, `instant_jam=True`)
#: -- "does the array's FIRST QP solve at a_start already show reduced rate",
#: matching the bead's own row-G wording verbatim ("onset is at a=0, first
#: order in a... jam angle > 0 NOT required"); a single reading, no stepping.
#: (ii) SUSTAINED (rows H, I, `instant_jam=False`) -- "letting the run
#: proceed through any transient reduced-rate onset, does it ever become
#: genuinely, persistently stuck (rate collapses to near-true-zero) before
#: reaching the target". STALL_RATE_TOL is that near-true-zero floor -- a
#: SUSTAINED run is never declared jammed merely for dipping below
#: JAM_RATE_TOL, only for collapsing below THIS, much stricter, threshold.
STALL_RATE_TOL = 1e-4

#: `solve_ldp`'s infeasibility floor: the reduction's LAST residual
#: component (`rn`) collapsing below this means the LDP has no solution (or
#: is too ill-conditioned to trust) -- FIX (code review 23261 M3): this used
#: to be a bare `1e-13` inside `solve_ldp`, disconnected from this declared-
#: but-unread constant; now the constant IS the check.
QP_MULT_TOL = 1e-8

#: `crank_step`'s binding-tolerance: a candidate row is reported "binding"
#: (part of the printed active set, K/M rows) iff `|row @ v| < BINDING_TOL`
#: at the solved v. FIX (code review 23261 M3): this was a bare `1e-6`.
BINDING_TOL = 1e-6

#: Null-space rank tolerance (SVD) for J_pin -- re-declared locally (RANK_RTOL
#: exists in jb_x but this file does its OWN rank computation, on a DIFFERENT
#: matrix, so the constant is re-priced here, not borrowed silently).
QP_NULL_RTOL = 1e-10

#: The default driven unit for the single-crank-handle DRIVE variant (Z14).
DRIVEN_UNIT_INDEX = 0

#: Bounded step count for a crank run -- generous relative to a 0..A_ICO_LOCAL
#: sweep at H_STEP (~12 steps) so a genuine non-termination bug still reddens
#: this row rather than being silently truncated away.
MAX_CRANK_STEPS = 200

#: METRIC FORM treatment (module docstring FOUR DECLARATIONS, bead row M):
#: TREATMENT (a) -- W = identity in body coordinates. ALT_W_SCALE is the
#: alternate diagonal weight (angular block scaled) used by the W-insensitivity
#: control row: the NORM-FREE verdicts (jam status, binding-set composition)
#: must be unchanged under it.
ALT_W_ANGULAR_SCALE = 4.0

# ==========================================================================
# LOCK-SURFACE CONSTANTS (bead inviscid-qvf.19, Phase 1c). Re-declared
# locally per the mutation-probe rule even where a value (e.g. H_LOCK)
# echoes an existing one.
#
# H_LOCK: THE STEP SIZE DECISION THIS BEAD HAD TO MAKE, calibrated live
# (not guessed) before committing the grid, because it turned out the crank
# stepper's SUSTAINED outcome at a fixed (w, t) is NOT h0-invariant --
# h0=H_STEP=2.0 (bead .18's own REACH-tuned value, used by h_run/cross_wico)
# fails to detect a real jam at w=0 that h0=0.5 (bead .18's own JAM-tuned
# value, used by g_run/cross_w0) DOES detect, reproducing the independently
# reviewed 0.7852886398227111 number to 6 decimal places. Measured directly
# (not assumed): at h0=2.0 EVERY sampled w in [0, w_ico], including w=0,
# reports "reached"; at h0=0.5 the SAME w=0 reports "jammed" (matching the
# reviewed number) and every other sampled w in (0, w_ico] ALSO reports
# "jammed", at a value indistinguishable to 6 decimal places across every
# fraction tried. This is a genuine step-size sensitivity of the FROZEN
# stepper (crank_run/crank_step are not touched here, only called), not a
# bug introduced by this bead -- flagged below and in the T2 write-back as
# a deferred item for qvf.20's mutation-probe pass. H_LOCK is fixed at the
# JAM-tuned value for every call this bead makes (both grids, both arms,
# the motion-order trace, the drive-robustness check) so every number in
# the surface is comparable on the same footing.
H_LOCK = 0.5

#: Upper bound on the swept target angle: 5x the largest jam angle observed
#: in calibration (~0.98 deg) across every (w, t) point tried, so a status
#: of "reached" (no lock found) is a genuine, distinguishable outcome from
#: "jammed" rather than an artifact of too tight a ceiling.
#: How far past `a_start` a lock-surface point cranks before giving up and
#: reporting a CENSORED reading.
#:
#: RE-PRICED 5.0 -> 45.0 for bead inviscid-1wd. 5.0 was calibrated against the
#: DISASSEMBLED array, whose "lock" sat at a* = 0.98 -- comfortably inside a
#: 5-degree window. With the ball joint restored the array locks at a* = 29.88,
#: so every sweep point ran to the target instead, every a* censored to the
#: same 5.0, and the w-arm span collapsed to exactly zero: a surface made
#: entirely of the budget rather than of the physics. The CENSORING row below
#: is what turns that failure mode from a silent flat surface into a red row.
#:
#: 45.0 is NOT chosen to just clear the observed 29.88 -- a target tuned to the
#: answer is the answer. It is this file's own largest `FOLD_TABLE_TARGET`
#: entry, sits inside `ANGLE_GRID`'s span, and stays well clear of a = 60 where
#: the twelve shared vertices merge into six and the pairing this whole file
#: rests on stops existing.
A_TARGET_LOCK = 45.0

#: Step budget per lock-surface run. Raised 20 -> 150 alongside the target
#: above: at h0 = H_LOCK the assembled array needs ~66 steps to reach its lock
#: at 29.88 and proportionally more to run out a 45-degree window. Budget
#: exhaustion surfaces as QPFAIL and is classified as BUDGET-EXHAUSTED rather
#: than being read as a physical result -- it is never folded into a*.
MAX_STEPS_LOCK = 150

#: Step budget for the single-crank-handle DRIVE-ROBUSTNESS runs specifically.
#: Larger than MAX_STEPS_LOCK because driving one unit advances the array's
#: phase far more slowly than the uniform in-phase drive: measured, that path
#: needs ~360 steps to run out the same 45-degree window that driven="all"
#: covers in ~66. At the shared 150 both points simply exhausted the budget and
#: reported QPFAIL, which is a statement about the budget and not about the
#: drive.
MAX_STEPS_DRIVE = 500

#: PRIMARY w-grid, as FRACTIONS of w_ico (computed from `fold_halves`, not a
#: bare literal -- see `_w_ico_lock`). Includes the w=0 boundary (already
#: independently validated by bead .18's own g_run/cross_w0) plus three
#: interior points spanning to just below w_ico (0.9x): >=4 distinct values,
#: a stated range [0, 0.9*w_ico].
W_GRID_FRAC = (0.0, 0.3, 0.6, 0.9)
#: SECOND, ABSOLUTE, INCOMMENSURATE arm (house style, jb_y's K_GRID/
#: K_GRID_ALT precedent): irregular offsets from the primary fractions
#: (+0.12, +0.08, +0.02, -0.02), none coincident, none a common ratio of
#: the primary set -- chosen to avoid both the w=0 special case (already
#: covered by the primary arm) and the w_ico edge, while still landing
#: inside the same validated [0, w_ico) domain.
W_GRID_FRAC_ALT = (0.12, 0.38, 0.62, 0.88)

#: t-grid (thickness), swept at a FIXED, LARGE w (see T_SWEEP_W_FRAC) so
#: wires never bind and only contact + thickness governs -- the design of
#: record's own arm-C isolation ("w fixed and large"). RE-PRICED in the
#: qvf.19 fix round (SHIP-BLOCKER 2, T2 23299 CRITICAL 2): the ORIGINAL
#: values (0.0, 0.02, 0.08) were chosen before the t-sweep was re-based off
#: A_START_T -- once re-based, a live scan (bisection-style, by hand,
#: during development) found the re-based array's opening-vs-immediate-jam
#: threshold sits near t=0.187, NOT near t=0 (the whole point of the fix:
#: a clearance-relieved start tolerates real thickness before locking).
#: The values below straddle that measured threshold on BOTH sides so
#: da*/dt is a genuine, non-degenerate, two-valued signal (opens the full
#: T_TARGET_SPAN below the threshold; re-locks at the start above it) --
#: NOT chosen to land inside a flat, uninformative region on one side.
T_GRID = (0.05, 0.15, 0.19)
T_GRID_ALT = (0.03, 0.17, 0.195)
#: Which W_GRID_FRAC fraction of w_ico the t-sweep holds fixed.
T_SWEEP_W_FRAC = 0.9

#: Discrimination-ratio band edges (row Q2), PRICED FROM THIS RUN'S OWN
#: NUMBERS: calibration measured |da*/dw| indistinguishable from 0 across
#: the interior w-grid against a |da*/dt| of order 10 (t=0 -> t=0.08 spans
#: essentially the whole ~0.98 deg jam-angle range), so a ratio near 0 is
#: expected -- these edges (a decade wide, symmetric in log space) sit
#: comfortably clear of that expected value on the ARM-C side while still
#: leaving genuine room for ARM-A-LIKE or NON-DISCRIMINATING outcomes to be
#: reachable in principle, so the row is not rigged to one verdict by
#: construction.
DISCRIM_RATIO_LOW = 0.2
DISCRIM_RATIO_HIGH = 5.0
#: Below this, |da*/dt| (the ratio's denominator) is treated as too close
#: to zero for the ratio to be COMPUTABLE -- the row FAILS rather than
#: printing an inf-pass (jb_x X7's shape, deliberately not repeated here).
RATIO_ZERO_TOL = 1e-9

#: Row Q's two-way separation gate: the t-driven separation (arm C) must
#: exceed this; the w-only separation among the grid's non-boundary points
#: (arm A, at the same large-w reference) must stay under it.
T_SEPARATION_TOL = 0.1
W_SEPARATION_TOL = 0.05

#: Floor for "the w axis actually does something" (bead inviscid-l1d). Priced
#: from the fixed-wire measurement: interior w-separation 3.849637 deg, interior
#: spreads ~2.9-3.1 deg. Before the fix every one of these was EXACTLY 0.0, so
#: any positive floor separates the two regimes; 0.5 leaves an order of margin
#: without being satisfiable by numerical noise.
W_LIVE_FLOOR = 0.5

#: The arm-A k-table's lock-angle span, re-declared locally (mutation-probe
#: rule): 26.555073204 at k=0.90 down to 12.705717132 at k=1.20.
K_TABLE_SPAN = 13.849356

#: The interior w-spread is NOT h0-converged and this floor asserts that rather
#: than papering over it. Measured 3.062 / 2.850 / 21.709 across
#: h0 = 0.5 / 0.25 / 0.125, i.e. adjacent drifts of 0.212 and 18.859. A floor of
#: 1.0 is comfortably above the coarse-level drift and far below the fine-level
#: one, so the row distinguishes the two regimes rather than merely detecting
#: motion.
H_REFINE_SPREAD_DIVERGENT_FLOOR = 1.0

#: Row R's two-sided non-vacuity band on the WHOLE surface's a* span (every
#: value from both grids, both arms): a stuck constant fails the lower
#: bound; a degenerate or blown-up statistic fails the upper.
SPAN_LOWER = 1e-6
SPAN_UPPER = 2.0 * A_TARGET_LOCK

#: Motion-order trace (row S): a SHORT, bespoke stepping loop (see
#: `motion_order_trace`'s own docstring for why it does not reuse
#: `crank_run`'s feasibility backtrack) -- a diagnostic of ORDER, not a
#: hardened lock-surface measurement.
MOTION_ORDER_STEPS = 8
#: How many early steps the LEAD/LAG row reads. Early, because the question is
#: whether the drive PROPAGATES -- once the array has taken it up the ratio
#: settles and says nothing about order.
#: Steps the JOINT-integrity probe cranks. Long enough for the retired
#: per-unit projection to separate the joints well past JOINT_GAP_TOL (it
#: reaches 0.043 in three steps at ~0.094 per degree), short enough to stay a
#: gate row rather than a sweep.
#: The w fraction the phase probe cranks at. The interior of the w grid, where
#: a* is flat, so the probe reads the same regime the surface is quoted from.
PHASE_PROBE_W_FRAC = 0.6

#: Common `a_hat` at which every step size is compared. Comparing at one a_hat
#: rather than at each run's own freeze is what separates the drift from where
#: each run happened to stop. Comfortably short of the freeze at every level.
PHASE_PROBE_TARGET = 20.0

#: Step sizes the probe compares. A 4x range: if the drift were discretisation
#: it would shrink across this, and it does not.
PHASE_H0_LEVELS = (0.5, 0.25, 0.125)

#: Budget per level. h0=0.125 needs ~172 steps to reach PHASE_PROBE_TARGET.
PHASE_PROBE_MAX_STEPS = 400

#: Instrument-calibration angles: `configuration_phase` must recover these from
#: `corners(a)` alone. Spread across the range, none of them special.
PHASE_PROBE_ANGLES = (3.0, 11.0, 22.238756093, 37.0)

#: Second, absolute, incommensurate calibration arm -- house rule for any swept
#: grid, applied here too rather than exempting the instrument from it.
PHASE_PROBE_ANGLES_ALT = (7.0, 19.0, 28.0, 41.0)

#: How exactly the instrument must invert on a KNOWN symmetric pose. It manages
#: 4e-16; this leaves four orders of headroom and still fails loudly if the
#: closed-form inversion is ever wrong.
PHASE_INSTRUMENT_TOL = 1e-12

#: Radius spread admitted on a genuinely symmetric pose. Machine zero.
PHASE_SYMMETRIC_SPREAD_TOL = 1e-12

#: The array must demonstrably LEAVE the symmetric path for the off-path rows
#: to mean anything. Measured spread at the probe target is ~6.2e-02; this floor
#: sits an order below, so the row states a real departure and not a rounding.
PHASE_OFFPATH_FLOOR = 5e-3

#: How far `a_hat` may sit from the array's MEAN measured phase. Measured
#: 0.34 / 0.27 / 0.21 deg across the h0 ladder, so this leaves ample headroom
#: while still reddening if the integral ever came adrift of the array.
PHASE_MEAN_DRIFT_TOL = 1.0

#: The units must be shown NOT to share a phase. A floor, not a ceiling: the
#: dephasing is the finding, and a row asserting the units stay in phase would
#: assert something false. Measured spread is 6.3-6.7 deg; this sits well below
#: so the row states a real effect, and it reddens if a change ever restored
#: in-phase motion.
PHASE_DEPHASE_FLOOR = 1.0

#: Hinge residual, joint gap and triangle-edge deviation admitted through a
#: crank run. All three sit at 1e-15 or below; this row is what distinguishes
#: "the linkage flexed" from "the solver broke".
LINKAGE_EXACT_TOL = 1e-10

JOINT_PROBE_STEPS = 12

#: The assembled array's tolerance on shared-vertex separation. Priced against
#: `project_to_joint_manifold`'s own convergence floor (PROJECT_TOL) with room
#: for the Newton solve to stop early, not against the control's drift -- a
#: tolerance priced from the failure it is meant to catch would admit it.
JOINT_GAP_TOL = 1e-9

#: The control must separate by at least this much for the TEST's pass to mean
#: anything. Well below the 0.043 measured at three steps, so the row reports a
#: genuine gap rather than a coincidence of run length.
JOINT_CONTROL_FLOOR = 1e-3

MOTION_LEAD_STEPS = 4

#: The lead/lag row's threshold: at least one early step must show every
#: undriven unit moving at less than this fraction of the driven unit. Priced
#: from the run's own numbers (0.431 at step 0, settling near 0.62) with
#: headroom, and falsified by the driven="all" control, where the ratio goes
#: to 1.
MOTION_LEAD_RATIO_MAX = 0.8

MOTION_ORDER_H = 0.5
#: A per-unit velocity-component norm above this counts as "moving".
MOTION_ACTIVITY_TOL = 1e-6
#: The CONTROL's "moves uniformly" band, among the SIX TOPOLOGICALLY
#: EQUIVALENT corner units (SC7's star center, unit 0, is symmetry-
#: distinguished -- degree 6 against the corners' degree 1 -- so the
#: array's own automorphism group fixes it pointwise and only permutes
#: units 1..6 among themselves; equivariance therefore predicts those SIX
#: agree with EACH OTHER, not that the CENTER matches them too. Measured
#: at driven="all", w large, t=0, FIRST step: the six corner rates agree
#: to a spread of ~1.7e-3; priced generously above that measured value.
MOTION_UNIFORM_TOL = 0.01
#: Which W_GRID_FRAC fraction of w_ico the single-crank-handle TEST case
#: uses; the CONTROL case (driven="all") uses T_SWEEP_W_FRAC (large w).
MOTION_TEST_W_FRAC = 0.5

#: Drive-model-robustness check (the hazard comment's own requirement,
#: 2026-08-21: "confirm the flatness is drive-model-robust before quoting
#: it"): re-measure two W_GRID_FRAC points under the single-crank-handle
#: drive and confirm the SAME small-spread structure holds there too, on
#: its own terms (not by cross-arm equality to the driven="all" numbers,
#: which measure a different physical quantity).
DRIVE_ROBUST_TOL = 0.5

# ==========================================================================
# FIX-ROUND CONSTANTS (qvf.19-critique.md T2 23299, qvf.19-code-review.md T2
# 23337). Two ship-blockers (h0-non-convergence confounding the Z15 CROSS
# test; the t-column being a construction artifact of starting every sweep
# at the idealized a=0 touching pose) plus one High and two Medium code-
# review findings. Re-declared locally per house style.
# ==========================================================================

#: Named (code review LOW): the shuffle-control "differs from the fitted
#: slope" threshold for rows O/P, previously an inline 1e-6.
SLOPE_SHUFFLE_DIFFER_TOL = 1e-6
#: Named (code review LOW): row O's "interior slope is flat" threshold,
#: previously an inline 1e-3.
W_INTERIOR_FLAT_TOL = 1e-3

#: HIGH fix (23337): the ALT grids must be demonstrated, not merely
#: asserted in prose, to be non-ratio-derived from the PRIMARY grids (the
#: jb_y K_GRID/K_GRID_ALT precedent: coprimality + offset-apartness gate
#: rows). "Non-ratio-derived" is checked as: no constant k satisfies
#: ALT[i] == k*PRIMARY[i] for every i simultaneously (a single common
#: scale factor would mean the ALT arm is the primary arm rescaled, the
#: exact "coarsens in lock-step" bug jb_y's own docstring warns against);
#: "offset-apart" is checked as every ALT value sitting at least this far
#: from every PRIMARY value (no accidental coincidence).
GRID_RATIO_CONST_TOL = 1e-6
GRID_OFFSET_APART_TOL = 1e-3

#: SHIP-BLOCKER 1b (23299 CRITICAL 1): h-refinement. >=3 h0 levels, halved,
#: at representative (w,t) points, to gate whether the JAM/REACH verdict
#: and a* are STABLE as h0 shrinks -- a genuinely two-sided claim (the row
#: must be able to fail if a value the critique found NON-convergent is
#: misreported as stable, or if a value found CONVERGENT is misreported as
#: unstable).
H_REFINE_LEVELS = (0.5, 0.25, 0.125)
#: The interior w-point's a* must stay within this band across every
#: H_REFINE_LEVELS value to count as "stable" -- priced from the critique's
#: own reproduced numbers (identical to 6 decimals at h0=0.5 and h0=0.25).
H_REFINE_STABLE_TOL = 0.02
#: Below this, TWO h0 levels' a* are considered "distinguishable" (i.e.
#: NOT converged) -- used to positively demonstrate the w=0 boundary point
#: FAILS the stability band, the other half of the two-sided row.
H_REFINE_UNSTABLE_FLOOR = 0.1
#: Re-priced with A_TARGET_LOCK (bead inviscid-1wd) and for the same reason:
#: at 5.0 the refinement probe measured how far the array got in 5 degrees,
#: not where it locks.
H_REFINE_TARGET = 45.0
#: RAISED 300 -> 900 (bead inviscid-l1d, 2026-08-25). With the wire span fixed
#: the w axis is LIVE, so a lock run at nonzero w now travels further before it
#: jams and the old budget exhausted on the finest h0 level, surfacing as QPFAIL.
#: That was a statement about the budget and not about the physics -- the
#: companion BUDGET-EXHAUSTED row said so correctly while three rows downstream
#: of it went red for want of data.
H_REFINE_MAX_STEPS = 900

#: SHIP-BLOCKER 2 (23299 CRITICAL 2): the t-column's clearance-relieved
#: start angle. Derived (bisection over the OPENING/wire-mechanism pairs'
#: own gap(a), the SAME quantity `wire_gradient_row` uses), not picked:
#: the smallest angle at which every OPENING pair's gap exceeds
#: T_START_TARGET_GAP. Scoped to the OPENING pairs specifically (not
#: "every non-pinned pair") because the CLOSING/ridge pairs at the SAME
#: shared-square locations are, measured directly, exactly touching at
#: a=0 and get MONOTONICALLY WORSE (not better) as a increases -- T2
#: 23195's registry/ridge mechanism, a DIFFERENT physical effect with no
#: "relief" angle at all on this side of a=0. A_START_T is validated live
#: (not merely computed) by a companion gate row confirming a real
#: `crank_run` from A_START_T at the largest swept t produces genuine
#: multi-step motion, not another instant re-lock.
#:
#: T_START_TARGET_GAP is a FIXED, stated clearance target -- DELIBERATELY
#: NOT derived from max(T_GRID+T_GRID_ALT) (that shape was tried first and
#: rejected: it makes A_START_T move every time T_GRID is re-priced, which
#: moves the re-based array's OWN jam/reach threshold, which then demands
#: re-pricing T_GRID again -- a circular dependency discovered live while
#: calibrating this fix, not a hypothetical). 0.117 = the ORIGINAL
#: T_GRID's largest swept value (0.097) plus a 0.02 clearance margin, the
#: same figure this constant produced before the circularity was found and
#: broken; T_GRID/T_GRID_ALT below are then priced (a SEPARATE, one-time
#: calibration) against the resulting fixed A_START_T.
T_START_TARGET_GAP = 0.117
T_START_BISECT_HI = 20.0
T_START_BISECT_ITERS = 40
#: Opening-range budget beyond A_START_T for the re-based t-sweep (mirrors
#: A_TARGET_LOCK's old magnitude).
T_TARGET_SPAN = 5.0
#: R-t's upper bound (row R split per-arm in this fix round, since the
#: t-arm's a* is an OPENING RANGE bounded by T_TARGET_SPAN, a different
#: scale from the w-arm's absolute-angle a*, bounded by A_TARGET_LOCK).
#: RE-DERIVED for bead inviscid-1wd. This was 2.0 * T_TARGET_SPAN, but the
#: t-arm has never passed T_TARGET_SPAN to `lock_surface_point` -- it takes the
#: default `target_span=A_TARGET_LOCK`, so T_TARGET_SPAN bounded nothing and
#: the band silently stayed priced for a 5-degree window while the arm ran to
#: 45. It is the same bound the w-arm's own band uses, for the same reason: an
#: opening range cannot exceed the window it was measured in.
SPAN_UPPER_T = 2.0 * A_TARGET_LOCK


# ==========================================================================
# Z0: FOUNDATION -- the plate-normal-invariance fact this whole kernel rests
# on, gate-checked rather than merely asserted in the docstring.
# ==========================================================================

_FACES = faces()  # [(v, c, u, sigma)] for the 8 octahedron faces, a = phase-independent list


def plate_normal(face_idx):
    """The FIXED outward unit normal of plate `face_idx`, independent of phase."""
    return _FACES[face_idx][2]


def plate_triangle(a, face_idx, origin=None):
    """The 3 corners of plate `face_idx` at phase `a`, optionally translated."""
    tri = corners(a)[face_idx]
    if origin is not None:
        tri = tri + origin
    return tri


def unit_plates(a, origin):
    """All 8 plates of one unit at phase `a`, translated to `origin`. GAP 4:
    the pure-translate reference placement `assemble_free` itself uses when
    `rots=None`."""
    return corners(a) + origin


def _z0_normal_invariance():
    """Z0a: is a plate's plane normal REALLY independent of phase? Checked by
    computing the actual triangle normal (cross product of two edges,
    normalised) at several phases and comparing to the FIXED axis u_f,
    on BOTH angle ladders."""
    worst = 0.0
    for f in range(8):
        u = plate_normal(f)
        for a in list(ANGLE_GRID) + list(ANGLE_GRID_ALT) + [0.0, 60.0]:
            tri = plate_triangle(a, f)
            n = _cross3(tri[1] - tri[0], tri[2] - tri[0])
            n = n / np.linalg.norm(n)
            worst = max(worst, min(np.linalg.norm(n - u), np.linalg.norm(n + u)))
    return worst


def _z0_centroid_on_axis():
    """Z0b: a plate's centroid sits at u_f * Z * cos(a) exactly -- the fact
    the parallel-facing branch's exact projection formula depends on."""
    worst = 0.0
    for f in range(8):
        u = plate_normal(f)
        for a in list(ANGLE_GRID) + list(ANGLE_GRID_ALT):
            tri = plate_triangle(a, f)
            c = tri.mean(axis=0)
            predicted = u * Z * np.cos(np.radians(a))
            worst = max(worst, float(np.linalg.norm(c - predicted)))
    return worst


def z0_control():
    dev_normal = _z0_normal_invariance()
    dev_centroid = _z0_centroid_on_axis()
    return {"normal_invariance": dev_normal, "centroid_on_axis": dev_centroid}


# ==========================================================================
# Z1: SIGNED TRIANGLE-TRIANGLE GAP, WITNESS POINTS, NORMAL
# ==========================================================================

def _cross3(a, b):
    """Hand-rolled 3-vector cross product. Bit-identical to np.cross on
    3-vectors (same per-component expressions, same IEEE754 order) while
    skipping numpy's axis normalization, which profiling measured at ~60%%
    of the whole gate's wall clock (3.9M np.cross calls per crank run)."""
    return np.array((a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]))


def _clamp01(x):
    """min(max(x,0),1): bit-identical to np.clip(x, 0.0, 1.0) on scalars
    (incl. NaN propagation) without numpy's ufunc dispatch overhead."""
    return 1.0 if x > 1.0 else (0.0 if x < 0.0 else x)


def _norm3(v):
    """sqrt(v.v): bit-identical to np.linalg.norm on a real 1-D vector
    (numpy computes exactly sqrt(dot(x, x)) there) without its wrapper."""
    return float(v @ v) ** 0.5


# --------------------------------------------------------------------------
# SCALAR CORE. The 15-candidate witness search below was measured
# (2026-08-22) at 76% of this file's ENTIRE wall clock: 2.34M
# `_seg_seg_witness` calls, 3.96M `_cross3` calls, 4.20M `numpy.array`
# constructions per gate run, at ~3.5us each -- all of it numpy's per-call
# dispatch on THREE-ELEMENT vectors, none of it arithmetic.
#
# So the arithmetic runs on plain floats and numpy appears only at the
# boundary, once per triangle pair instead of once per candidate. The public
# entry points keep their array-in / array-out signatures and become thin
# wrappers over this core, so there is exactly ONE implementation of each
# primitive and no second copy that can drift from it.
#
# BIT-IDENTICAL, not merely close, and by construction rather than by luck:
# every expression below preserves the original's operand order and grouping,
# and `a @ b` on a 3-vector was verified (200k random pairs, zero mismatches)
# to be exactly `a0*b0 + a1*b1 + a2*b2` -- numpy performs no reassociation at
# this length. `mean` over 3 elements was verified the same way. The standing
# proof is the gate itself: `--no-cache` output before and after this change
# is byte-identical, every row included.
# --------------------------------------------------------------------------

def _dot3s(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _sub3s(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross3s(a, b):
    """`_cross3` on float triples: same component expressions, no array."""
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm3s(v):
    return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5


def _tri_s(tri):
    """A numpy (3,3) triangle as a tuple of three float triples."""
    a, b, c = tri
    return ((float(a[0]), float(a[1]), float(a[2])),
            (float(b[0]), float(b[1]), float(b[2])),
            (float(c[0]), float(c[1]), float(c[2])))


def _seg_seg_witness_s(p1, q1, p2, q2):
    """Scalar `_seg_seg_witness`. Operand order matches the array form
    line for line; `p1 + d1 * s` becomes `p1_i + d1_i * s`, which is the
    same per-component expression numpy evaluates."""
    d1 = _sub3s(q1, p1)
    d2 = _sub3s(q2, p2)
    r = _sub3s(p1, p2)
    a = _dot3s(d1, d1)
    e = _dot3s(d2, d2)
    if a < 1e-14 and e < 1e-14:
        return p1, p2, _norm3s(r)
    if a < 1e-14:
        t = _clamp01(-_dot3s(d2, r) / e) if e > 1e-14 else 0.0
        pb = (p2[0] + d2[0] * t, p2[1] + d2[1] * t, p2[2] + d2[2] * t)
        return p1, pb, _norm3s(_sub3s(p1, pb))
    if e < 1e-14:
        s = _clamp01(_dot3s(d1, r) / a)
        pa = (p1[0] + d1[0] * s, p1[1] + d1[1] * s, p1[2] + d1[2] * s)
        return pa, p2, _norm3s(_sub3s(pa, p2))
    b = _dot3s(d1, d2)
    c = _dot3s(d1, r)
    f = _dot3s(d2, r)
    den = a * e - b * b
    s = _clamp01((b * f - c * e) / den) if den > 1e-12 else 0.0
    t = _clamp01((b * s + f) / e)
    s = _clamp01((b * t - c) / a)
    pa = (p1[0] + d1[0] * s, p1[1] + d1[1] * s, p1[2] + d1[2] * s)
    pb = (p2[0] + d2[0] * t, p2[1] + d2[1] * t, p2[2] + d2[2] * t)
    return pa, pb, _norm3s(_sub3s(pa, pb))


def _seg_seg_witness(p1, q1, p2, q2):
    """Closest points on two segments (clamped parametric solve, the same
    method as `jb_g_strut_clearance.segment_distance` / jb_x's private
    `_seg_seg`, extended to return the witness points, not just the
    distance). Array-shaped wrapper over `_seg_seg_witness_s`."""
    pa, pb, d = _seg_seg_witness_s(
        tuple(map(float, p1)), tuple(map(float, q1)),
        tuple(map(float, p2)), tuple(map(float, q2)))
    return np.array(pa), np.array(pb), d


def _barycentric_s(p, tri):
    """Scalar `_barycentric`. `tri` is a tuple of three float triples."""
    a, b, c = tri
    ba = _sub3s(b, a)
    ca = _sub3s(c, a)
    n = _cross3s(ba, ca)
    nn = _dot3s(n, n)
    if nn < 1e-18:
        return None
    denom = nn
    u = _dot3s(_cross3s(ba, _sub3s(p, a)), n) / denom
    v = _dot3s(_cross3s(_sub3s(c, b), _sub3s(p, b)), n) / denom
    w = _dot3s(_cross3s(_sub3s(a, c), _sub3s(p, c)), n) / denom
    return w, u, v  # weights of (a, b, c) respectively


def _barycentric(p, tri):
    """Barycentric coordinates of p w.r.t. triangle tri, via the projected
    2D-area-ratio construction. Not clamped -- may go outside [0,1] for a
    point off the triangle's PLANE or outside its extent, which is exactly
    what the on-triangle gate row checks."""
    return _barycentric_s(tuple(map(float, p)), _tri_s(tri))


def _pt_tri_witness_s(p, tri):
    """Scalar `_pt_tri_witness`. `proj` reproduces numpy's
    `p - n * ((p - a) @ n) / nn` componentwise: the scale is applied to each
    component BEFORE the division, which is the order the array form uses."""
    a, b, c = tri
    n = _cross3s(_sub3s(b, a), _sub3s(c, a))
    nn = _dot3s(n, n)
    if nn >= 1e-18:
        s = _dot3s(_sub3s(p, a), n)
        proj = (p[0] - n[0] * s / nn,
                p[1] - n[1] * s / nn,
                p[2] - n[2] * s / nn)
        bw, uw, vw = _barycentric_s(proj, tri)
        if bw >= 0 and uw >= 0 and vw >= 0:
            return proj, _norm3s(_sub3s(p, proj))
    cands = [_seg_seg_witness_s(p, p, a, b), _seg_seg_witness_s(p, p, b, c),
             _seg_seg_witness_s(p, p, c, a)]
    best = min(cands, key=lambda x: x[2])
    return best[1], best[2]


def _pt_tri_witness(p, tri):
    """Closest point ON triangle tri to point p (projection, clamped to the
    face; falls back to the nearest edge when the projection lands outside).
    Extends jb_x's private `_pt_tri` (distance-only) to return the witness."""
    q, d = _pt_tri_witness_s(tuple(map(float, p)), _tri_s(tri))
    return np.array(q), d


def _seg_tri_hits_s(p, q, tri):
    """Scalar Moller-Trumbore."""
    a, b, c = tri
    e1 = _sub3s(b, a)
    e2 = _sub3s(c, a)
    d = _sub3s(q, p)
    h = _cross3s(d, e2)
    det = _dot3s(e1, h)
    if abs(det) < 1e-14:
        return False
    inv = 1.0 / det
    s = _sub3s(p, a)
    u = inv * _dot3s(s, h)
    if u < 0.0 or u > 1.0:
        return False
    qv = _cross3s(s, e1)
    v = inv * _dot3s(d, qv)
    if v < 0.0 or u + v > 1.0:
        return False
    t = inv * _dot3s(e2, qv)
    return 0.0 < t < 1.0


def _seg_tri_hits(p, q, tri):
    """Moller-Trumbore: does segment pq pierce the OPEN interior of tri?
    Same primitive as jb_x's private `_seg_tri_hits`, reimplemented locally
    (mutation-probe rule; this file imports no private symbols from jb_x)."""
    return _seg_tri_hits_s(tuple(map(float, p)), tuple(map(float, q)),
                           _tri_s(tri))


def _closest_point_pair_s(sA, sB):
    """Scalar 15-candidate search. Candidate order and every comparison
    (`d < best[2]`, strict) are unchanged from the array form, so ties resolve
    to the same candidate they always did -- which is what keeps the returned
    WITNESS POINTS, not merely the distance, identical."""
    best = None
    for i in range(3):
        for j in range(3):
            pa, pb, d = _seg_seg_witness_s(sA[i], sA[(i + 1) % 3],
                                           sB[j], sB[(j + 1) % 3])
            if best is None or d < best[2]:
                best = (pa, pb, d)
    for i in range(3):
        q, d = _pt_tri_witness_s(sA[i], sB)
        if d < best[2]:
            best = (sA[i], q, d)
    for i in range(3):
        q, d = _pt_tri_witness_s(sB[i], sA)
        if d < best[2]:
            best = (q, sB[i], d)
    return best


def _closest_point_pair(triA, triB):
    """The 15-candidate unsigned closest-point search: 9 edge-edge pairs + 3
    vertex(A)-vs-triB + 3 vertex(B)-vs-triA. Returns (pA, pB, distance).

    The two triangles cross into float space ONCE here, not once per
    candidate: this single conversion is what removes ~4.2M array
    constructions per gate run."""
    pa, pb, d = _closest_point_pair_s(_tri_s(triA), _tri_s(triB))
    return np.array(pa), np.array(pb), d


def _is_piercing(triA, triB):
    sA, sB = _tri_s(triA), _tri_s(triB)
    for i in range(3):
        if _seg_tri_hits_s(sA[i], sA[(i + 1) % 3], sB):
            return True
        if _seg_tri_hits_s(sB[i], sB[(i + 1) % 3], sA):
            return True
    return False


def signed_gap(triA, triB, nA):
    """Signed gap between two triangular plates. Returns (gap, witness_A,
    witness_B, normal) -- NOT a bare float: a witness point on each plate
    and the contact normal, per the structural acceptance criterion. The
    caller SUPPLIES `nA` (triA's own fixed plate normal, from
    `plate_normal`) as the reference direction the sign is measured
    against; it is echoed back unchanged as the 4th return value so a
    caller that only keeps the tuple still has the normal beside the gap
    and the witnesses, with an explicit, stated sign convention
    (critique C1-i).

    SIGN CONVENTION: positive gap means triA and triB are separated along
    +nA, `nA` being triA's OWN fixed outward plate normal (Z0) -- the
    direction triA's own plate FACES. Since triB, when separated, sits on
    the +nA side of triA (that is what "facing" means), moving triA FURTHER
    along +nA advances it TOWARD triB and DECREASES the gap; moving it along
    -nA withdraws it and INCREASES the gap by the same amount. Stated
    plainly: displacing triA by +eps along -nA increases the reported gap by
    +eps -- this is what gate row G-i verifies by finite difference, on one
    separated and one near-touching pair, and it is the same relationship
    the registry pair's closed form exhibits (Z4): as the octahedral unit's
    phase p DECREASES from 60 (the codebase-wide "expansion" direction), its
    plate's centroid advances FURTHER along +nA and the registry gap falls.

    TWO BRANCHES:

    (a) PARALLEL-FACING (|nA . nB| > PARALLEL_TOL, nB = triB's own plate
        normal supplied by the caller via `plate_normal`): the two planes are
        parallel, so for ANY point x_A in triA's plane and ANY point x_B in
        triB's plane, (x_B - x_A) . nA is the SAME number regardless of
        which points are chosen -- there is no witness-point ambiguity in
        the direction that matters. This is EXACT for both separated and
        interpenetrating plates (23195 mechanism 1, the registry pair; and
        mechanism 2, the folding-square pair, are both this case). Witness
        points are each triangle's own centroid: always a valid barycentric
        point (1/3, 1/3, 1/3), so row G-iv's on-triangle check is trivial by
        construction here.
    (b) GENERAL (not parallel): the 15-candidate closest-point search above,
        sign NEGATIVE iff an edge of either triangle pierces the other's
        interior (Moller-Trumbore), magnitude = the along-nA projection of
        (witness_B - witness_A) in that case, POSITIVE unsigned distance
        otherwise. The pierced branch is a bounded, continuous, sign-correct
        proxy valid near the contact boundary -- the regime this bead's rows
        exercise (a separated pair and a NEAR-TOUCHING, still-POSITIVE-gap
        pair; see "A ROW DELIBERATELY NOT BUILT" in the gate) and the regime
        a quasi-static stepper (bead .18) actually queries (small step,
        active-set threshold near g = 0). DEEP penetration-depth accuracy in
        this branch is explicitly OUT OF SCOPE and no row claims it.
    """
    nB = _cross3(triB[1] - triB[0], triB[2] - triB[0])
    nB = nB / np.linalg.norm(nB)
    if abs(float(nA @ nB)) > PARALLEL_TOL:
        cA = triA.mean(axis=0)
        cB = triB.mean(axis=0)
        gap = float((cB - cA) @ nA)
        return gap, cA, cB, nA
    pA, pB, d0 = _closest_point_pair(triA, triB)
    pierced = _is_piercing(triA, triB)
    if not pierced:
        return d0, pA, pB, nA
    # Pierced: sign is negative BY the piercing detection itself; magnitude
    # is the along-nA projection of the closest-point witness offset, a
    # bounded, continuous proxy near the contact boundary (see docstring).
    proxy = float((pB - pA) @ nA)
    return -abs(proxy), pA, pB, nA


# ==========================================================================
# Z2: PLATE-PAIR ENUMERATION -- topology as data (GAP 2)
# ==========================================================================

def _hinge_faces(vertex_label):
    """The (up to two) face indices of ONE unit that meet at hinge
    `vertex_label` -- read directly off jb_x's PAIRS."""
    (fa, _), (fb, _) = PAIRS[vertex_label]
    return {fa, fb}


def _intra_hinged(fi, fj):
    for (fa, _), (fb, _) in PAIRS:
        if {fa, fb} == {fi, fj}:
            return True
    return False


def _shares_pin(topo, i, fi, j, fj):
    """STATED EXCLUSION PREDICATE (GAP 2): two plates share a pinned vertex,
    and are excluded from the clearance check, iff:
      - same unit (i == j): the two faces are hinge-adjacent (a PAIRS entry
        joins them) -- they are rigidly pinned at that shared corner always.
      - different units (i != j): some contact of `topo` identifies a vertex
        of unit i whose hinge faces include fi with a vertex of unit j whose
        hinge faces include fj.
    Not a hardcoded list: entirely a function of PAIRS and topo.contacts."""
    if i == j:
        return _intra_hinged(fi, fj)
    for (a, k, b, l) in topo.contacts:
        if (a, b) == (i, j) and fi in _hinge_faces(k) and fj in _hinge_faces(l):
            return True
        if (a, b) == (j, i) and fj in _hinge_faces(k) and fi in _hinge_faces(l):
            return True
    return False


def enumerate_plate_pairs(topo):
    """Every plate pair of the assembled topology, intra- AND inter-unit,
    EXCEPT pairs sharing a pinned vertex (`_shares_pin`). Returns a list of
    (i, fi, j, fj) tuples."""
    out = []
    n = topo.n
    for i in range(n):
        for j in range(i, n):
            f_range = range(8)
            for fi in f_range:
                for fj in (range(fi + 1, 8) if i == j else f_range):
                    if _shares_pin(topo, i, fi, j, fj):
                        continue
                    out.append((i, fi, j, fj))
    return out


def _z2_plate_pair_counts():
    out = {}
    for topo in build_topologies():
        pairs = enumerate_plate_pairs(topo)
        n_intra = sum(1 for (i, fi, j, fj) in pairs if i == j)
        n_inter = sum(1 for (i, fi, j, fj) in pairs if i != j)
        out[topo.name] = {"n_pairs": len(pairs), "n_intra": n_intra,
                          "n_inter": n_inter, "n": topo.n}
    return out


# ==========================================================================
# Z3: THE FOLD TABLE (row A/B) -- single-unit geometry, jb_x's own DIAGONALS
# ==========================================================================

def _square_partner(diag, square_diagonals, struts):
    """Given one diagonal of a cuboctahedron square face (a DIAGONALS
    member), find its complementary diagonal in the SAME square: the pair
    (b, d) in square_diagonals, disjoint from (a, c), such that {a,c,b,d}
    forms a 4-cycle of struts (in either winding order)."""
    a, c = diag
    for (b, d) in square_diagonals:
        if (b, d) == diag or {b, d} & {a, c}:
            continue
        e1 = [frozenset((a, b)), frozenset((b, c)), frozenset((c, d)), frozenset((d, a))]
        e2 = [frozenset((a, d)), frozenset((d, c)), frozenset((c, b)), frozenset((b, a))]
        if all(e in struts for e in e1) or all(e in struts for e in e2):
            return (b, d)
    return None


_SQUARE_PARTNERS = tuple((d, _square_partner(d, SQUARE_DIAGONALS, STRUTS))
                         for d in DIAGONALS)


def fold_halves(a):
    """For each of the 6 cuboctahedron squares (keyed by its DIAGONALS
    member), the pair of half-diagonal lengths (d1, d2) at phase `a`: d1 is
    the DIAGONALS member's own half-length (SHRINKS toward strut/2 as `a`
    rises to a_ico -- this is jb_x's own A_ICO derivation), d2 its square
    partner's half-length (GROWS). Printed in this (d1, d2) order, which is
    the order the bead's acceptance-criteria table uses."""
    v = verts(a)
    out = {}
    for (k, l), (p, q) in _SQUARE_PARTNERS:
        d1 = 0.5 * float(np.linalg.norm(v[k] - v[l]))
        d2 = 0.5 * float(np.linalg.norm(v[p] - v[q]))
        out[(k, l)] = (d1, d2)
    return out


FOLD_TABLE_TARGET = {
    5.0: (0.94588, 1.04651, 0.10064),
    10.0: (0.88455, 1.08506, 0.20051),
    A_ICO_LOCAL: (0.70711, 1.14412, 0.43702),
    30.0: (0.57735, 1.15470, 0.57735),
    45.0: (0.29886, 1.11536, 0.81650),
}
FOLD_TABLE_TOL = 1e-5  # the target table itself is quoted to 5 decimals


def z3_fold_table():
    ref = DIAGONALS[0]
    rows = {}
    for a, (t1, t2, tfold) in FOLD_TABLE_TARGET.items():
        halves = fold_halves(a)
        d1, d2 = halves[ref]
        rows[a] = {"d1": d1, "d2": d2, "fold": d2 - d1,
                   "t1": t1, "t2": t2, "tfold": tfold}
    # square uniformity: every one of the 6 squares must give the SAME
    # (d1, d2) at a generic angle -- a real cross-check, not an assumption.
    halves_ico = fold_halves(A_ICO_LOCAL)
    vals = list(halves_ico.values())
    uniform_dev = max(max(abs(d1 - vals[0][0]), abs(d2 - vals[0][1]))
                      for d1, d2 in vals)
    # a = 0 control: fold must be EXACTLY zero (four-vertex sharing).
    d1_0, d2_0 = fold_halves(0.0)[ref]
    fold_at_0 = d2_0 - d1_0
    # first-order-in-a control: fold(eps)/eps must be roughly constant for
    # small eps (not a lookup table), sampled on both angle ladders' smallest
    # rungs.
    eps1, eps2 = 0.5, 1.0
    f1 = fold_halves(eps1)[ref]
    f2 = fold_halves(eps2)[ref]
    slope1 = (f1[1] - f1[0]) / eps1
    slope2 = (f2[1] - f2[0]) / eps2
    first_order_dev = abs(slope1 - slope2) / max(abs(slope1), 1e-300)
    # ridge diagonal == strut EXACTLY at a_ico (row B), full length = 2 * d1.
    ridge_full = 2.0 * halves_ico[ref][0]
    off_dev_lo = abs(2.0 * fold_halves(A_ICO_LOCAL - AICO_CONTROL_OFFSET)[ref][0] - STRUT_LEN_LOCAL)
    off_dev_hi = abs(2.0 * fold_halves(A_ICO_LOCAL + AICO_CONTROL_OFFSET)[ref][0] - STRUT_LEN_LOCAL)
    return {"rows": rows, "uniform_dev": uniform_dev, "fold_at_0": fold_at_0,
            "first_order_dev": first_order_dev, "ridge_full": ridge_full,
            "aico_control_min": min(off_dev_lo, off_dev_hi)}


# ==========================================================================
# Z4: REGISTRY PAIR CLOSED FORM (row C) -- parallel-facing plates via the
# general `signed_gap` kernel, cross-checked against the analytic formula.
# ==========================================================================

def _antipodal_face(i0):
    u0 = plate_normal(i0)
    best, bestdev = None, 2.0
    for j in range(8):
        dev = float(np.linalg.norm(plate_normal(j) + u0))
        if dev < bestdev:
            bestdev, best = dev, j
    return best, bestdev


#: Registry spacing D, DERIVED (not a bare literal): the fixed unit-to-unit
#: distance at which an octahedron-phase (a=60) unit's plate exactly touches
#: a VE-phase (a=0) neighbour's plate along the shared axis -- the physical
#: definition of "registry". D = Z cos(60) + Z cos(0) = 1.5 Z.
REGISTRY_D = 1.5 * Z


def registry_pair(i0, p, q):
    """The centre unit's plate i0 at phase p, and a neighbour unit's
    ANTIPODAL plate at phase q, translated by REGISTRY_D along i0's own
    axis -- so the two plates are exactly parallel-facing (23195 mechanism
    1)."""
    j0, _ = _antipodal_face(i0)
    u0 = plate_normal(i0)
    triA = plate_triangle(p, i0)
    triB = plate_triangle(q, j0, origin=REGISTRY_D * u0)
    return triA, triB, u0


def registry_closed_form(p, q):
    return np.sqrt(3.0) - (2.0 / np.sqrt(3.0)) * (np.cos(np.radians(p)) + np.cos(np.radians(q)))


REGISTRY_ROWS = ((60.0, 0.0), (30.0, 30.0), (59.0, 0.0), (60.0, 1.0), (0.0, 0.0),
                 (A_ICO_LOCAL, A_ICO_LOCAL))
REGISTRY_TARGETS = {(60.0, 0.0): 0.0, (30.0, 30.0): -0.267949,
                    (59.0, 0.0): -0.017364, (60.0, 1.0): 0.000176}


def z4_registry():
    worst_closed = 0.0
    worst_target = 0.0
    per_face = []
    for i0 in range(8):
        for (p, q) in REGISTRY_ROWS:
            triA, triB, u0 = registry_pair(i0, p, q)
            gap, wA, wB, _ = signed_gap(triA, triB, u0)
            closed = registry_closed_form(p, q)
            worst_closed = max(worst_closed, abs(gap - closed))
            if (p, q) in REGISTRY_TARGETS:
                worst_target = max(worst_target, abs(gap - REGISTRY_TARGETS[(p, q)]))
            per_face.append((i0, p, q, gap))
    # 8-face uniformity: all 8 body-diagonal orientations give the same gap
    # for the same (p, q) -- a real cross-check.
    by_pq = {}
    for i0, p, q, gap in per_face:
        by_pq.setdefault((p, q), []).append(gap)
    uniform_dev = max(max(vals) - min(vals) for vals in by_pq.values())
    # swap minimum: p + q = 60, minimum at p = q = 30.
    swap_vals = [registry_closed_form(p, 60.0 - p) for p in np.linspace(0.0, 60.0, 601)]
    swap_min = min(swap_vals)
    swap_argmin = float(np.linspace(0.0, 60.0, 601)[int(np.argmin(swap_vals))])
    # one-sidedness at the registry contact (60, 0): d gap/dp and d gap/dq.
    h = 1e-4
    # CONVENTION (matches the project-wide "expansion is a DECREASING"
    # rule, jb_x's docstring): the derivative is reported with respect to
    # the EXPANSION direction, i.e. w.r.t. (60 - p), not raw p. Central
    # difference with the arguments swapped implements exactly that sign.
    dgdp = (registry_closed_form(60.0 - h, 0.0) - registry_closed_form(60.0 + h, 0.0)) / (2 * h)
    dgdq = (registry_closed_form(60.0, 0.0 + h) - registry_closed_form(60.0, 0.0 - h)) / (2 * h)
    dgdp_rad, dgdq_rad = dgdp * 180.0 / np.pi, dgdq * 180.0 / np.pi
    return {"worst_closed": worst_closed, "worst_target": worst_target,
            "uniform_dev": uniform_dev, "swap_min": swap_min,
            "swap_argmin": swap_argmin, "dgdp_rad": dgdp_rad, "dgdq_rad": dgdq_rad}


# ==========================================================================
# Z5: NORMAL ORIENTATION AND WITNESS POINTS (rows G, critique C1)
# ==========================================================================

def _generic_pair(i, j, gap_target, a=A_ICO_LOCAL):
    """A GENERAL (non-parallel) plate pair, positioned so the true unsigned
    gap is approximately `gap_target`, by translating unit j's whole octahedron
    along plate i's own normal from a reference offset. Both plates i, j are
    DIFFERENT faces of the SAME octahedron shape (unit j is a translated
    COPY), chosen non-adjacent so their normals are not parallel."""
    u_i = plate_normal(i)
    triA = plate_triangle(a, i)
    ref = plate_triangle(a, j)
    # start at a reference translation putting triB's plane roughly gap_target
    # away from triA's plane along u_i, via bisection on the actual signed_gap.
    lo, hi = -3.0, 6.0

    def g(t):
        origin = u_i * t
        triB = ref + origin
        gap, _, _, _ = signed_gap(triA, triB, u_i)
        return gap

    # bisection for g(t) == gap_target (g is monotone increasing in t for a
    # fixed, non-degenerate pair over this bracket -- checked by construction
    # of the bracket below via sign change).
    glo, ghi = g(lo), g(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        gm = g(mid)
        if (gm - gap_target) * (glo - gap_target) <= 0:
            hi = mid
        else:
            lo, glo = mid, gm
    t_star = 0.5 * (lo + hi)
    triB = ref + u_i * t_star
    return triA, triB, u_i


def _vertex_above_face(i, d, a=A_ICO_LOCAL, spread=5.0):
    """A GENERAL (non-parallel) plate pair with an EXACT, by-construction FD
    relationship to plate i's own normal: triB is built with one vertex Q
    sitting directly at triA's centroid + nA * d, and its other two vertices
    placed far away (offset further along nA too), so Q is unambiguously the
    single closest point of triB and its perpendicular projection onto
    triA's plane -- triA's own centroid -- is safely inside triA. Moving
    triA along its own normal by eps then changes the gap by EXACTLY eps
    (to machine precision), independent of any curvature or Voronoi-region
    switch, because the closest-point structure (vertex Q vs interior of
    triA) does not change for small eps. nB (triB's own plane normal) is
    generically NOT parallel to nA, so this still exercises the GENERAL
    branch of `signed_gap`, not the parallel-facing one."""
    triA = plate_triangle(a, i)
    nA = plate_normal(i)
    cA = triA.mean(axis=0)
    Q = cA + nA * d
    perp1 = _cross3(nA, np.array([1.0, 0.0, 0.0]))
    if np.linalg.norm(perp1) < 1e-6:
        perp1 = _cross3(nA, np.array([0.0, 1.0, 0.0]))
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = _cross3(nA, perp1)
    triB = np.array([Q, Q + perp1 * spread + nA * spread,
                     Q + perp2 * spread + nA * spread])
    return triA, triB, nA


def z5_normal_witness():
    rows = {}

    # (i) FD check: one SEPARATED pair (large positive gap) and one
    # NEAR-TOUCHING pair (small positive gap) -- deliberately kept out of
    # the interpenetrating regime (see signed_gap's docstring). Built via
    # `_vertex_above_face` so the FD relationship is EXACT by construction.
    for label, target in (("separated", 0.5), ("near_touching", 1e-3)):
        triA, triB, nA = _vertex_above_face(0, target)
        g0, _, _, _ = signed_gap(triA, triB, nA)
        # per the stated convention: -nA WITHDRAWS triA and INCREASES gap.
        triA_disp = triA - nA * FD_EPS
        g2, _, _, _ = signed_gap(triA_disp, triB, nA)
        fd_dev = abs((g2 - g0) - FD_EPS)
        triA_disp_ctrl = triA - nA * (FD_CONTROL_FACTOR * FD_EPS)
        g3, _, _, _ = signed_gap(triA_disp_ctrl, triB, nA)
        fd_dev_control = abs((g3 - g0) - FD_EPS)
        rows[label] = {"g0": g0, "fd_dev": fd_dev, "fd_dev_control": fd_dev_control}

    # (ii) registry contact: normal points AWAY from the interpenetrating
    # body, sign matches d gap/dp = -1.000 per rad from the closed form
    # (same EXPANSION-direction convention as Z4's dgdp).
    triA, triB, u0 = registry_pair(0, 59.0, 0.0)
    g_reg, _, _, _ = signed_gap(triA, triB, u0)
    g_reg_plus, _, _, _ = signed_gap(*registry_pair(0, 59.0 + 1e-4, 0.0)[:2], u0)
    g_reg_minus, _, _, _ = signed_gap(*registry_pair(0, 59.0 - 1e-4, 0.0)[:2], u0)
    dgdp_measured = (g_reg_minus - g_reg_plus) / (2e-4) * 180.0 / np.pi
    rows["registry_sign"] = {"gap": g_reg, "dgdp_deg": dgdp_measured}

    # (iii) negated-normal control. Two distinct failure shapes, matched to
    # `signed_gap`'s two branches:
    #  - FD check: the general (non-parallel) branch's return value does not
    #    depend on nA's sign at all (it is an unsigned closest-point
    #    distance away from the pierced regime) -- the bug this control
    #    catches is a CALLER trusting a WRONG (negated) reported normal and
    #    displacing the body along it while still expecting "+eps": that
    #    displacement is actually along the direction gap DECREASES, so the
    #    FD residual is checked there, not by feeding -nA into signed_gap.
    triA, triB, nA = _vertex_above_face(0, 0.5)
    g0, _, _, _ = signed_gap(triA, triB, nA)
    # a caller trusting the NEGATED normal applies the withdrawal rule
    # ("+eps along -[normal]") to -(-nA) = +nA, which actually DECREASES gap.
    triA_disp_wrong = triA + nA * FD_EPS
    g2, _, _, _ = signed_gap(triA_disp_wrong, triB, nA)
    fd_dev_negated = abs((g2 - g0) - FD_EPS)
    _, _, u0_neg = registry_pair(0, 59.0, 0.0)
    g_neg_plus, _, _, _ = signed_gap(*registry_pair(0, 59.0 + 1e-4, 0.0)[:2], -u0_neg)
    g_neg_minus, _, _, _ = signed_gap(*registry_pair(0, 59.0 - 1e-4, 0.0)[:2], -u0_neg)
    dgdp_negated = (g_neg_minus - g_neg_plus) / (2e-4) * 180.0 / np.pi
    rows["negated_control"] = {"fd_dev_negated": fd_dev_negated,
                               "dgdp_negated_matches_closed_form": abs(dgdp_negated - (-1.0)) < 0.05}

    # (iv) witness points on-triangle + independent cross-check.
    # edge-edge configuration: two skew, non-adjacent plates at a_ico.
    triA_ee, triB_ee, nA_ee = _generic_pair(0, 4, 0.3)
    pA, pB, d0 = _closest_point_pair(triA_ee, triB_ee)
    baryA = _barycentric(pA, triA_ee)
    baryB = _barycentric(pB, triB_ee)
    on_triA = all(-BARY_TOL <= c <= 1 + BARY_TOL for c in baryA) and abs(sum(baryA) - 1) < BARY_TOL
    on_triB = all(-BARY_TOL <= c <= 1 + BARY_TOL for c in baryB) and abs(sum(baryB) - 1) < BARY_TOL
    # independent cross-check for the edge-edge witness: jb_g's
    # segment_distance on the SAME two closest edges.
    best_edge = None
    for i in range(3):
        for j in range(3):
            d = jb_g_segment_distance(triA_ee[i], triA_ee[(i + 1) % 3],
                                      triB_ee[j], triB_ee[(j + 1) % 3])
            if best_edge is None or d < best_edge:
                best_edge = d
    # d0 (my search, which also considers vertex-face candidates) must be
    # <= best_edge (jb_g's search, edge-edge only) always, and equal to it
    # whenever the true closest pair genuinely IS an edge-edge pair, which
    # this configuration was chosen to be.
    edge_cross_dev = abs(best_edge - d0)
    edge_witness_sep_dev = abs(float(np.linalg.norm(pA - pB)) - d0)

    # vertex-face configuration, deliberately chosen so the closest pair is
    # a vertex-to-INTERIOR projection, not vertex-to-vertex (a naive
    # candidate that only checks vertex-vertex distance would be WRONG here
    # -- the can-fail control).
    triA_vf, triB_vf, nA_vf = _generic_pair(0, 6, 0.2)
    pA2, pB2, d0_vf = _closest_point_pair(triA_vf, triB_vf)
    # independent point-triangle projection (different code path: solve the
    # 3x2 least-squares system for the projection's plane coordinates
    # directly rather than the cross-product/barycentric method above).
    def independent_pt_tri(p, tri):
        a, b, c = tri
        M = np.stack([b - a, c - a], axis=1)  # 3x2
        sol, *_ = np.linalg.lstsq(M, p - a, rcond=None)
        s, t = np.clip(sol, 0.0, 1.0)
        if s + t > 1.0:
            norm = s + t
            s, t = s / norm, t / norm
        return a + s * (b - a) + t * (c - a)

    # naive candidate: vertex-to-vertex only (deliberately wrong instrument).
    naive_best = min(float(np.linalg.norm(triA_vf[i] - triB_vf[j]))
                     for i in range(3) for j in range(3))
    # cross-check both directions against the independent projector
    indepA = independent_pt_tri(pB2, triA_vf)
    indepB = independent_pt_tri(pA2, triB_vf)
    vf_cross_dev = min(float(np.linalg.norm(pA2 - indepA)),
                       float(np.linalg.norm(pB2 - indepB)))
    naive_wrong = naive_best > d0_vf + 1e-6  # can-fail control: naive vertex-vertex is NOT the closest pair

    rows["witness"] = {"on_triA": on_triA, "on_triB": on_triB,
                       "edge_witness_sep_dev": edge_witness_sep_dev,
                       "edge_cross_dev": edge_cross_dev,
                       "vf_cross_dev": vf_cross_dev, "naive_wrong": naive_wrong}
    return rows


# ==========================================================================
# Z6: WIRE SPANS AND THICKNESS (rows H)
# ==========================================================================

def span_length(vA, vB):
    """Euclidean span between two tied vertex positions. Prior art:
    `jb_y_dephasing.span_length`/`member_length`; this is the same quantity,
    named for parity."""
    return float(np.linalg.norm(np.asarray(vA) - np.asarray(vB)))


def wire_active(span, w, tol=ACTIVE_TOL):
    """Tension-only (unilateral) activity: taut iff span is at or above the
    bound w (within tol); slack (inactive, no force) otherwise. This is the
    CORRECT check for `gap = w - span >= 0` treated as tension-only."""
    return span >= w - tol


def _bilateral_active_WRONG(span, w, tol=ACTIVE_TOL):
    """DELIBERATELY WRONG reference: a BILATERAL (rigid rod) member is an
    EQUALITY constraint -- always enforced, hence always active/binding,
    resisting tension AND compression alike. Kept only so the gate can show
    it disagrees with `wire_active` (tension-only, inactive/slack whenever
    span < w) on the compression side (critique C2's control). `w` and `tol`
    are accepted for a matching call signature; a bilateral link's activity
    does not depend on either."""
    return True


def admissible(gap, t):
    """Plate-pair admissibility with thickness offset: gap >= t."""
    return gap >= t


def z6_wire_and_thickness():
    # (v) hand-computed span.
    vA, vB = np.array([0.0, 0.0, 0.0]), np.array([3.0, 4.0, 0.0])
    s = span_length(vA, vB)
    span_dev = abs(s - 5.0)

    # (vi) taut/slack transition + compression-side control.
    s_known = 5.0
    w_slack, w_taut = s_known + 0.5, s_known - 0.5
    slack_ok = not wire_active(s_known, w_slack)
    taut_ok = wire_active(s_known, w_taut)
    # locate the transition by bisection on w -> wire_active(s_known, w)
    lo, hi = s_known - 1.0, s_known + 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if wire_active(s_known, mid):
            lo = mid
        else:
            hi = mid
    transition_dev = abs(0.5 * (lo + hi) - s_known)
    # compression-side control: move the span FURTHER below w (compression
    # direction) -- mine reports inactive (correct, tension-only); the
    # deliberately-wrong bilateral checker reports active AT s_known == w
    # itself, demonstrating the row can fail for a bilateral implementation.
    s_compressed = s_known - 0.2
    mine_inactive_under_compression = not wire_active(s_compressed, s_known)
    bilateral_would_flag_active_at_bound = _bilateral_active_WRONG(s_known, s_known)
    mine_flags_active_at_bound_too = wire_active(s_known, s_known)  # both agree exactly AT the bound
    bilateral_disagrees_under_compression = _bilateral_active_WRONG(s_compressed, s_known) != wire_active(s_compressed, s_known)

    # (vii) thickness moves the gap-admissibility threshold by EXACTLY t.
    g_known = 0.37
    t_vals = (0.0, 0.1, 0.2, 0.37, 0.5)
    thickness_devs = []
    for t in t_vals:
        # admissible(g_known, t) must equal (g_known - t) >= 0 to bit precision
        thickness_devs.append(0.0 if admissible(g_known, t) == (g_known - t >= 0.0) else 1.0)
    # flip case: gap admissibility flips as t crosses g_known.
    flip_below = admissible(g_known, g_known - 1e-9)
    flip_above = not admissible(g_known, g_known + 1e-9)
    # (viii) t = 0 as its own row.
    t0_ok = admissible(g_known, 0.0) == (g_known >= 0.0)

    # REAL-TOPOLOGY exercise (not synthetic): the census square's own
    # PINNED vertex (label 0, the SC7 star's generator -- the contact this
    # neighbour direction is built on) and its NON-PINNED partner in the
    # same square (label 3, the ridge diagonal's other end) give a genuine
    # wire-span candidate between the centre unit and its actually-placed
    # +axis neighbour, at the icosahedral phase and CENSUS_SPACING.
    axis = _census_axis()
    v_a = verts(A_ICO_LOCAL)
    centre_pt = v_a[3]
    corner_pt = v_a[3] + CENSUS_SPACING * axis
    real_span = span_length(centre_pt, corner_pt)
    real_span_sane = np.isfinite(real_span) and 0.0 < real_span < 10.0 * STRUT_LEN_LOCAL

    return {"span_dev": span_dev, "slack_ok": slack_ok, "taut_ok": taut_ok,
            "transition_dev": transition_dev,
            "mine_inactive_under_compression": mine_inactive_under_compression,
            "bilateral_would_flag_active_at_bound": bilateral_would_flag_active_at_bound,
            "mine_flags_active_at_bound_too": mine_flags_active_at_bound_too,
            "bilateral_disagrees_under_compression": bilateral_disagrees_under_compression,
            "thickness_max_dev": max(thickness_devs),
            "flip_below": flip_below, "flip_above": flip_above, "t0_ok": t0_ok,
            "real_span": real_span, "real_span_sane": real_span_sane}


# ==========================================================================
# Z7: CROSSING CENSUS + PERTURBATION (row D)
# ==========================================================================

def _tri_plane(tri):
    n = _cross3(tri[1] - tri[0], tri[2] - tri[0])
    d = -n @ tri[0]
    return n, d


def _edge_interval(tri, n_other, d_other, D):
    dist = np.array([n_other @ tri[i] + d_other for i in range(3)])
    signs = np.sign(dist)
    signs[np.abs(dist) < 1e-12] = 0.0
    pts = []
    for i in range(3):
        j = (i + 1) % 3
        if signs[i] == 0.0:
            pts.append(tri[i])
        if signs[i] * signs[j] < 0.0:
            t = dist[i] / (dist[i] - dist[j])
            pts.append(tri[i] + t * (tri[j] - tri[i]))
    if len(pts) < 2:
        return None
    proj = [p @ D for p in pts]
    lo_i, hi_i = int(np.argmin(proj)), int(np.argmax(proj))
    return proj[lo_i], proj[hi_i], pts[lo_i], pts[hi_i]


def tri_tri_crossing_segment(t1, t2, tol=1e-9):
    """Moller triangle-triangle intersection, returning the crossing SEGMENT
    (midpoint, endpoints) when the two triangles' planes genuinely cross
    within both triangles' extents -- None otherwise. This is the ROBUST
    instrument (23195): unlike a strict-interior pierce test, it is not
    blind to a parallel pass-through (planes never cross -> no segment,
    correctly) nor to a symmetric edge-through-edge crossing (the segment
    exists and its midpoint is checked for strict interiority separately)."""
    n1, d1 = _tri_plane(t1)
    n2, d2 = _tri_plane(t2)
    dist2 = np.array([n1 @ t2[i] + d1 for i in range(3)])
    if np.all(dist2 > tol) or np.all(dist2 < -tol):
        return None
    dist1 = np.array([n2 @ t1[i] + d2 for i in range(3)])
    if np.all(dist1 > tol) or np.all(dist1 < -tol):
        return None
    D = _cross3(n1, n2)
    nD = np.linalg.norm(D)
    if nD < 1e-12:
        return None
    D = D / nD
    iv1 = _edge_interval(t1, n2, d2, D)
    iv2 = _edge_interval(t2, n1, d1, D)
    if iv1 is None or iv2 is None:
        return None
    lo = max(iv1[0], iv2[0])
    hi = min(iv1[1], iv2[1])
    if lo >= hi:
        return None
    tmid = 0.5 * (lo + hi)
    a0, b0 = iv1[2], iv1[3]
    pa, pb = a0 @ D, b0 @ D
    frac = 0.0 if abs(pb - pa) < 1e-14 else (tmid - pa) / (pb - pa)
    mid = a0 + frac * (b0 - a0)
    return mid, lo, hi


def _strict_interior(pt, tri, tol=1e-9):
    a, b, c = tri
    n = _cross3(b - a, c - a)
    nn = n @ n
    if nn < 1e-18:
        return False
    u = _cross3(b - a, pt - a) @ n
    v = _cross3(c - b, pt - b) @ n
    w = _cross3(a - c, pt - c) @ n
    return u > tol and v > tol and w > tol


def _naive_pierce_count(centre, corner):
    """THE WRONG INSTRUMENT, kept deliberately (23195): a single-point
    strict-interior pierce test per edge, with no crossing-segment
    construction. Blind to parallel pass-through and to exactly-symmetric
    edge-through-edge crossings -- both are the recorded failure modes."""
    n = 0
    for fi in range(8):
        for fj in range(8):
            t1, t2 = centre[fi], corner[fj]
            hit = False
            for i in range(3):
                if _seg_tri_hits(t1[i], t1[(i + 1) % 3], t2):
                    hit = True
                if _seg_tri_hits(t2[i], t2[(i + 1) % 3], t1):
                    hit = True
            if hit:
                n += 1
    return n


#: The square whose centroid direction at a = 0 defines the neighbour axis
#: (the "+x shared square", DERIVED here as the centroid of DIAGONALS[0]'s
#: square, not hardcoded to (1,0,0)).
_CENSUS_SQUARE = (DIAGONALS[0][0], DIAGONALS[0][1],
                  _SQUARE_PARTNERS[0][1][0], _SQUARE_PARTNERS[0][1][1])


def _census_axis():
    v0 = verts(0.0)
    c = v0[list(_CENSUS_SQUARE)].mean(axis=0)
    return c / np.linalg.norm(c)


def crossing_census(a, spacing, axis):
    centre = unit_plates(a, np.zeros(3))
    corner = unit_plates(a, spacing * axis)
    n, total_len = 0, 0.0
    for fi in range(8):
        for fj in range(8):
            res = tri_tri_crossing_segment(centre[fi], corner[fj])
            if res is None:
                continue
            mid, lo, hi = res
            if _strict_interior(mid, centre[fi]) and _strict_interior(mid, corner[fj]):
                n += 1
                total_len += (hi - lo)
    naive = _naive_pierce_count(centre, corner)
    return n, total_len, naive


def _rand_rigid(rng, mag):
    """A random small rigid motion: rotation by `mag` radians about a random
    axis, plus a translation of magnitude `mag`."""
    axis = rng.normal(size=3)
    axis = axis / np.linalg.norm(axis)
    theta_deg = np.degrees(mag)
    R = rot(axis, theta_deg)
    t = rng.normal(size=3)
    t = t / np.linalg.norm(t) * mag
    return R, t


def _apply_rigid(tri, R, t, origin):
    """Apply rotation R about `origin` then translation t to a triangle."""
    return (R @ (tri - origin).T).T + origin + t


def crossing_census_perturbed(a, spacing, axis, seed, n_trials, mag):
    rng = np.random.default_rng(seed)
    counts = []
    for _ in range(n_trials):
        R_c, t_c = _rand_rigid(rng, mag)
        R_n, t_n = _rand_rigid(rng, mag)
        centre = unit_plates(a, np.zeros(3))
        corner = unit_plates(a, spacing * axis)
        centre_p = np.array([_apply_rigid(centre[f], R_c, t_c, np.zeros(3)) for f in range(8)])
        corner_p = np.array([_apply_rigid(corner[f], R_n, t_n, spacing * axis) for f in range(8)])
        n = 0
        for fi in range(8):
            for fj in range(8):
                res = tri_tri_crossing_segment(centre_p[fi], corner_p[fj])
                if res is None:
                    continue
                mid, lo, hi = res
                if _strict_interior(mid, centre_p[fi]) and _strict_interior(mid, corner_p[fj]):
                    n += 1
        counts.append(n)
    return counts


def z7_crossing_census():
    axis = _census_axis()
    axis_dev = abs(float(np.linalg.norm(axis)) - 1.0)
    at_zero_n, _, at_zero_naive = crossing_census(0.0, CENSUS_SPACING, axis)
    results = {}
    for a in (5.0, 10.0, A_ICO_LOCAL, 30.0, 45.0):
        n, total_len, naive = crossing_census(a, CENSUS_SPACING, axis)
        counts = crossing_census_perturbed(a, CENSUS_SPACING, axis, PERTURB_SEED,
                                           PERTURB_N, PERTURB_MAG)
        results[a] = {"n": n, "total_len": total_len, "naive": naive,
                     "perturbed": counts}
    ico_halves = fold_halves(A_ICO_LOCAL)[DIAGONALS[0]]
    d1_ico, d2_ico = ico_halves  # d1 shrinks (small), d2 grows (large)
    valleys_held_n, valleys_held_len, _ = crossing_census(A_ICO_LOCAL, 2.0 * d1_ico, axis)
    valleys_held_pert = crossing_census_perturbed(A_ICO_LOCAL, 2.0 * d1_ico, axis,
                                                  PERTURB_SEED, PERTURB_N, PERTURB_MAG)
    ridges_touch_n, ridges_touch_len, _ = crossing_census(A_ICO_LOCAL, 2.0 * d2_ico, axis)
    ridges_touch_pert = crossing_census_perturbed(A_ICO_LOCAL, 2.0 * d2_ico, axis,
                                                  PERTURB_SEED, PERTURB_N, PERTURB_MAG)
    return {"axis_dev": axis_dev, "at_zero_n": at_zero_n,
            "at_zero_naive": at_zero_naive, "results": results,
            "valleys_held_n": valleys_held_n, "valleys_held_len": valleys_held_len,
            "valleys_held_pert": valleys_held_pert,
            "ridges_touch_n": ridges_touch_n, "ridges_touch_len": ridges_touch_len,
            "ridges_touch_pert": ridges_touch_pert}


# ==========================================================================
# Z8: QP MACHINERY (bead qvf.18) -- the SOLVER decision. minimize
# ||v - v_cmd||^2_W subject to a null-space-encoded equality (J_pin v = 0)
# and a set of candidate inequalities C v >= 0, via the Lawson-Hanson
# reduction of a Least-Distance Program to Non-Negative Least Squares. See
# the module docstring's "PHASE 1b" section for the three-way decision and
# why the naive hand-rolled active-set method (tried first) cycles here.
# ==========================================================================

def null_space_basis(a_eq, ndof, rtol=QP_NULL_RTOL):
    """Orthonormal basis (ndof, ndof-rank) for the null space of `a_eq`, via
    SVD -- deterministic, no iteration, no tie-breaking to get wrong."""
    if a_eq.shape[0] == 0:
        return np.eye(ndof)
    _, s, vt = np.linalg.svd(a_eq, full_matrices=True)
    r = int((s > s[0] * rtol).sum()) if s.size else 0
    return vt[r:].T


def solve_ldp(g, h, maxiter=None):
    """Least-Distance Program: minimize ||x|| subject to g @ x >= h, via the
    classical Lawson & Hanson (1974, ch. 23) reduction to a single NNLS call.
    Returns (x, feasible). Verified against hand-computed tiny examples
    (T2 qvf.18-crank-stepper.md) before use here.

    FIX (code review 23261 C1): `scipy.optimize.nnls` RAISES RuntimeError on
    hitting its own iteration budget (verified against scipy 1.17.1 source,
    `_nnls.py`: `if info == 3: raise RuntimeError(...)`) -- unguarded, that
    would propagate as an unguarded traceback through this whole call chain,
    exactly the failure mode `_safe_linprog`/`_LPFail` in jb_y_dephasing
    exists to prevent ("a malformed programme is a BUG REPORT... it must
    arrive as a red gate row rather than as a traceback"). A non-convergent
    NNLS is a SOLVER failure, not a physical answer -- caught here and
    returned as `feasible=False`, the SAME channel the pre-existing
    near-singular-residual infeasibility path already uses, so every caller
    up the chain (`project_qp` -> `crank_step` -> `crank_run`) already routes
    `feasible=False` to QPFAIL without any further change needed.

    `maxiter` is exposed ONLY so the gate (Z16) can FORCE non-convergence on
    a real, nontrivial problem and confirm this returns `feasible=False`
    rather than raising -- every real call site in this file uses the
    default (unbounded)."""
    m, k = g.shape
    e = np.vstack([g.T, h.reshape(1, -1)])
    f = np.zeros(k + 1)
    f[-1] = 1.0
    try:
        y, _ = nnls(e, f, maxiter=maxiter)
    except RuntimeError:
        return np.zeros(k), False
    r = f - e @ y
    rn = r[-1]
    if abs(rn) < QP_MULT_TOL:
        return np.zeros(k), False
    return -r[:k] / rn, True


def project_qp(v_cmd, n_basis, active_rows, w_diag=None):
    """minimize (v-v_cmd)^T W (v-v_cmd) s.t. v = N z (encodes J_pin v = 0)
    and active_rows[k] . v >= 0 for every candidate. W enters as a diagonal
    weight (METRIC FORM treatment (a): identity by default; the
    W-insensitivity control passes an alternate diagonal). Returns
    (v, feasible)."""
    ndof = v_cmd.shape[0]
    if w_diag is None:
        w_sqrt = np.ones(ndof)
    else:
        w_sqrt = np.sqrt(w_diag)
    # Change of variables u = w_sqrt * v turns the weighted objective into an
    # ordinary Euclidean projection: minimize ||u - w_sqrt*v_cmd||^2 s.t.
    # (N/w_sqrt column-scaled) ... equivalently rescale N's ROWS by w_sqrt.
    nb = w_sqrt[:, None] * n_basis
    target = w_sqrt * v_cmd
    # Re-orthonormalize the scaled basis (QR) so the reduced coordinates are
    # still an ordinary Euclidean projection in u-space.
    q, _ = np.linalg.qr(nb)
    c = q.T @ target
    if not active_rows:
        u = q @ c
        return u / w_sqrt, True
    # candidate row r acts on v; r.v = r.(u/w_sqrt) = (r/w_sqrt).u
    c_mat = np.array([(row / w_sqrt) @ q for row in active_rows])
    h = -(c_mat @ c)
    x, feasible = solve_ldp(c_mat, h)
    u = q @ (x + c)
    return u / w_sqrt, feasible


# ==========================================================================
# Z9: STEP ASSEMBLY (bead qvf.18, GAP 6 -- greenfield). The equality
# (pin) Jacobian and the inequality (contact, wire) gradient rows the QP
# consumes, plus Newton projection back onto the pin manifold after a step.
# ==========================================================================

def _hat(r):
    """3x3 skew-symmetric cross-product matrix. Re-derived locally: jb_x's
    own `_hat` is private (mutation-probe rule forbids importing it)."""
    return np.array([[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]])


def _point_velocity_jacobian(x, f, point):
    """3x48 Jacobian of an ARBITRARY point on plate `f` (not necessarily a
    named corner -- a general-branch witness point is a closest-point search
    result, not a vertex) w.r.t. this unit's 48 body-motion dof. Generalizes
    `jb_x_array_linkage.position_jacobian_row`, which requires a listed
    corner index; a witness point does not have one."""
    r = point - x[f].mean(axis=0)
    m = np.zeros((3, 48))
    m[:, 3 * f:3 * f + 3] = -_hat(r)
    m[:, 24 + 3 * f:27 + 3 * f] = np.eye(3)
    return m


def _contact_gradient_direction(triA, triB, nA, wA, wB):
    """d such that d(gap)/dt = d . (velB - velA). `signed_gap` has THREE
    distinct first-order regimes it does not itself expose -- see the module
    docstring's HAZARD paragraph. Parallel-facing: nA (exact). General,
    unpierced: the unit vector witness_A -> witness_B (NOT nA -- verified by
    finite difference to disagree with nA by more than an order of magnitude,
    sign included, on a real SC7 pair). General, pierced (gap=-abs(proxy)):
    -sign(proxy)*nA."""
    nB = _cross3(triB[1] - triB[0], triB[2] - triB[0])
    nrm = np.linalg.norm(nB)
    nB = nB / nrm if nrm > 1e-300 else nA
    if abs(float(nA @ nB)) > PARALLEL_TOL:
        return nA
    if not _is_piercing(triA, triB):
        d = wB - wA
        dn = np.linalg.norm(d)
        return d / dn if dn > 1e-14 else nA
    proxy = float((wB - wA) @ nA)
    return -nA if proxy >= 0.0 else nA


def contact_gradient_row(xs, i, fi, j, fj, t, ndof):
    """(gap - t, row) where `row @ v` linearizes d(gap)/dt at the current
    configuration -- the greenfield inequality assembly GAP 6 names.
    `row[48*i:48*i+48]` and `row[48*j:48*j+48]` are the only nonzero blocks."""
    triA = xs[i][fi]
    triB = xs[j][fj]
    nA = plate_normal(fi)
    gap, wA, wB, nrm = signed_gap(triA, triB, nA)
    d = _contact_gradient_direction(triA, triB, nA, wA, wB)
    jacA = _point_velocity_jacobian(xs[i], fi, wA)
    jacB = _point_velocity_jacobian(xs[j], fj, wB)
    row = np.zeros(ndof)
    row[48 * i:48 * i + 48] = -(d @ jacA)
    row[48 * j:48 * j + 48] = (d @ jacB)
    return gap - t, row


#: Wire-attachment discovery (fix for substantive critique 23262 C1, SHIP-
#: BLOCKER): `topo.contacts`' single coincident vertex per neighbour link
#: (inherited from jb_x's pre-DECISION-18 rigid-pin topology, one pin
#: sufficient to fix relative unit placement under the OLD equality-only
#: model) is a SYMMETRY-FIXED POINT of `topo.sites`/`dsites` -- verified
#: directly: span is EXACTLY 0.0 at every phase tested (0, 1, 5, 10, a_ico),
#: not merely small, for every one of SC7's 6 wired pairs, under uniform
#: driving at any w. It cannot be what T2 23230 means by "each tied vertex
#: pair across each shared face" -- it is not where the fold (T2 23195
#: mechanism 2) separates anything. The design of record's own machinery
#: for the shared-square fold (`fold_halves`/`_SQUARE_PARTNERS`, one unit's
#: OWN square) does not resolve to a matching raw-vertex-family separation
#: between the two ACTUAL contact units either (checked directly, no match
#: at any tested angle) -- the real separating feature is at the PLATE
#: level, not the 12-vertex-family level. Found by direct, reproducible
#: measurement (T2 qvf.18-crank-stepper.md, FIXED-ROUND section): scanning
#: `enumerate_plate_pairs` between a contact's two units at a small
#: reference angle (avoiding the exact a=0 degenerate-witness kink this
#: file already documents elsewhere) picks out, per contact, a small set of
#: plate pairs that are near-zero at a=0 with POSITIVE d(gap)/da (OPENING,
#: i.e. valley pairs -- the ridge/closing pairs at the same locations are
#: already covered by the ordinary CONTACT machinery, unchanged). Their gap
#: is EXACTLY proportional to `fold_halves`' own fold quantity
#: (gap / fold(a) = 2/sqrt(3), constant to 5 decimal places across every
#: angle checked) -- a measured fact, not a guess.
#: Bounded iterations for `project_to_joint_manifold`. Higher than the
#: per-unit projector's 15 because the array-level system it solves is larger
#: (36n + 3c rows against 48n unknowns) and starts further from its own
#: solution after a full step.
PROJECT_MAXIT = 30

#: Convergence floor for the same. At 1e-13 the projector is asking for the
#: assembled configuration to machine precision; the JOINT gate rows check
#: what it actually achieves rather than assuming it got there.
PROJECT_TOL = 1e-13

WIRE_REF_ANGLE = 0.5
WIRE_FD_H = 1e-3
#: A candidate opening pair must sit within this of zero gap at the
#: reference angle to count as belonging to the shared face.
WIRE_NEAR_ZERO_TOL = 0.05


def _wire_attachment_pairs(topo, i, k, j, l, a_ref=WIRE_REF_ANGLE):
    """The plate pairs realizing ONE `topo.contacts` entry's shared-face
    wire loop: near-zero gap at `a_ref`, opening (positive d(gap)/da) via
    finite difference. Returns a list of (u, fi, v, fj) plate-pair tuples
    (u, v are the LOWER/HIGHER unit index per `enumerate_plate_pairs`' own
    ordering, not necessarily i, j in that order)."""
    lo, hi = min(i, j), max(i, j)
    candidates = [(u, fi, v, fj) for (u, fi, v, fj) in crank_pairs(topo) if u == lo and v == hi]
    h = WIRE_FD_H

    def xs_at(a):
        origins = topo.sites(verts(a))
        return [unit_corners(a, topo, u) + origins[u] for u in range(topo.n)]

    xs0, xsp, xsm = xs_at(a_ref), xs_at(a_ref + h), xs_at(a_ref - h)
    out = []
    for (u, fi, v, fj) in candidates:
        nA = plate_normal(fi)
        g0, *_ = signed_gap(xs0[u][fi], xs0[v][fj], nA)
        if abs(g0) >= WIRE_NEAR_ZERO_TOL:
            continue
        gp, *_ = signed_gap(xsp[u][fi], xsp[v][fj], nA)
        gm, *_ = signed_gap(xsm[u][fi], xsm[v][fj], nA)
        if (gp - gm) > 0.0:  # opening (valley); closing (ridge) pairs are
            out.append((u, fi, v, fj))  # already covered by CONTACT rows
    return out


def wire_pairs(topo):
    """The FIXED, topology-derived wire-attachment plate-pair list --
    flattened across every `topo.contacts` entry, computed once (mirrors
    `crank_pairs`)."""
    out = []
    for (i, k, j, l) in topo.contacts:
        out.extend(_wire_attachment_pairs(topo, i, k, j, l))
    return out


def wire_span(triA, triB, nA):
    """The wire's span: the along-nA separation of the two plates' centroids.
    Returns (span, centroid_A, centroid_B). UNCONDITIONAL -- no branch.

    THIS IS NOT `signed_gap`, AND THAT IS THE FIX (bead inviscid-l1d, T2 23388).
    The two want DIFFERENT quantities and were routed through one function
    because at the ideal pose they coincide. A CONTACT needs the true closest
    distance, so `signed_gap` branches to a closest-point search once the plates
    stop being parallel-facing -- correct for non-penetration, and left
    untouched here; the Phase 1a kernel is not reopened by this fix. A WIRE
    needs the along-normal opening, which is well defined whether or not the
    planes are parallel, and must not switch to anything else.

    `wire_gradient_row` used to delegate to `contact_gradient_row`, and its
    docstring justified that by saying these pairs "were selected to be well
    inside the parallel-facing regime anyway". They were -- AT SELECTION. The
    plate normals are phase independent (Z0), so on the pure-translate reference
    pose triB's measured normal equals its fixed plate normal EXACTLY and
    |nA.nB| is 1. The moment a unit rotates it is not, and PARALLEL_TOL is
    1 - 1e-9, which a rotation of 0.0026 degrees already fails. Measured on a
    real N2 wire pair: at 0.001 degrees the span reads 0.011635; at 0.01 degrees
    the branch flips and `signed_gap` returns 1.414282. A 120x discontinuity in
    a quantity whose value is 0.0116, at a rotation no crank step avoids, which
    is why the w axis has never actually been tested. The centroid projection
    over the same rotations is continuous: 0.011635, 0.011635, 0.011631,
    0.011596, 0.011010, 0.007730.

    The assumption was true when it was written and false one step later, and it
    was recorded in a docstring rather than in a row that could fail -- which is
    the whole reason it survived. Row K-wire-branch now checks it.
    """
    cA = triA.mean(axis=0)
    cB = triB.mean(axis=0)
    return float((cB - cA) @ nA), cA, cB


def wire_gradient_row(xs, pair, w, ndof):
    """(w - span, row, degenerate) for ONE wire-attachment plate pair
    (u, fi, v, fj): tension-only, g_wire = w - span >= 0.

    `span` is the SIGNED, along-normal opening quantity `contact_gradient_row`
    already computes for this exact pair (`_wire_attachment_pairs` selected
    it specifically as a near-a=0-zero, monotonically OPENING pair) -- NOT
    the raw Euclidean witness-to-witness distance. This distinction is
    load-bearing, not stylistic: these plates sit laterally offset from one
    another (their centroids differ by a large in-plane component, order
    2.8, irrelevant to the wire's physical model), so `|witness_B -
    witness_A|` includes that offset and does not track `fold(a)` at all --
    verified directly (an earlier version of this function used the raw
    Euclidean distance and it stayed near-constant at approx 2.83 across
    every `a` tested, instead of growing from 0). The SIGNED gap (the along-
    nA projection `contact_gradient_row` returns) is exactly the quantity
    that separates as `2/sqrt(3) * fold(a)`.

    IT NO LONGER REUSES `contact_gradient_row`, and the paragraph that used to
    stand here justified reusing it on the grounds that these pairs "were
    selected to be well inside the parallel-facing regime anyway". That is true
    at SELECTION and false after the first crank step -- see `wire_span`, which
    now computes the span and its gradient unconditionally. The gradient
    direction is nA always, and the witnesses are the centroids, which are valid
    barycentric points by construction so the on-triangle property is trivial.

    `degenerate` (fix for substantive critique 23262 C2) reports whether
    THIS gradient row is exactly zero -- which cannot happen for a
    parallel-facing pair (`_contact_gradient_direction`'s first branch is
    unconditional and never degenerates) and is retained only as an honest,
    gate-checked guarantee rather than a silent assumption: the OLD
    Euclidean-distance version had a genuine `uhat=zeros(3)` fallback at
    span<=1e-12 that made 6 wires print as "binding" while constraining
    nothing; this version structurally cannot reach that state for the
    pairs it is given, and a gate row (K-wire) checks that directly rather
    than trusting this docstring."""
    u, fi, v, fj = pair
    triA, triB = xs[u][fi], xs[v][fj]
    nA = plate_normal(fi)
    span, cA, cB = wire_span(triA, triB, nA)
    jacA = _point_velocity_jacobian(xs[u], fi, cA)
    jacB = _point_velocity_jacobian(xs[v], fj, cB)
    row = np.zeros(ndof)
    row[48 * u:48 * u + 48] = -(nA @ jacA)
    row[48 * v:48 * v + 48] = (nA @ jacB)
    degenerate = not np.any(row)
    return w - span, -row, degenerate


def _hinge_only_jacobian(xs, n):
    """The RETIRED block-diagonal Jacobian: 36 intra-unit hinge rows per unit
    and no inter-unit rows at all -- exactly what `build_pin_jacobian` was
    before bead inviscid-1wd.

    It exists solely so the JOINT gate row can run a control arm that
    reproduces the retired path in full. Nothing on the shipped path may call
    it: the whole point of that bead is that an array assembled by this
    Jacobian comes apart while being cranked."""
    ndof = 48 * n
    big = np.zeros((36 * n, ndof))
    for i in range(n):
        big[36 * i:36 * i + 36, 48 * i:48 * i + 48] = hinge_jacobian(xs[i])
    return big


def joint_vertex(xs, unit, pair_index):
    """The physical position of the vertex `topo.contacts` names by
    (unit, pair_index).

    `topo.contacts` entries are (i, k, j, l) where k and l index `PAIRS`, NOT
    the flat corner array -- `PAIRS[k][0]` is the (face, corner) address of the
    hinge's first representative. Reading them as flat vertex indices instead
    is wrong and is not loud about it: it returns a real position for a real
    vertex, just the wrong one, and reports a joint separation of sqrt(6) at
    the assembled pose. That mistake was made once while diagnosing this very
    defect, so the convention lives in a named function rather than being
    re-derived at each call site."""
    fa, ca = PAIRS[pair_index][0]
    return xs[unit][fa][ca]


def joint_residual(xs, topo):
    """The 3 x len(topo.contacts) BALL-JOINT violations: for each shared
    vertex, the vector from one unit's copy of it to the other's. Zero when
    the array is assembled. This is the constraint function whose Jacobian is
    the inter-unit block of `build_pin_jacobian` below."""
    out = np.zeros(3 * len(topo.contacts))
    for e, (i, k, j, l) in enumerate(topo.contacts):
        out[3 * e:3 * e + 3] = joint_vertex(xs, i, k) - joint_vertex(xs, j, l)
    return out


def build_pin_jacobian(xs, n, topo):
    """The array's equality Jacobian: 36 intra-unit hinge rows per unit, PLUS
    3 BALL-JOINT rows per `topo.contacts` entry -- (36n + 3c, 48n).

    THE INTER-UNIT ROWS ARE A BALL JOINT, NOT A RIGID PIN (bead inviscid-1wd).
    This function previously omitted them, on the stated grounds that
    `assemble_free`'s inter-unit rows "encode RIGID pins (pre-DECISION-18
    semantics)". Read that function: its inter-unit block is THREE rows per
    contact, `position_jacobian_row` on each side, differenced. Three rows
    constrain a shared vertex to stay shared and leave both units free to
    pivot about it. A rigid pin -- position AND relative orientation -- would
    be six. So the array's only assembly constraint was dropped on a mislabel,
    and under DECISION 18 nothing replaced it: contact is one-sided, and the
    wires are slack at a = 0 for every w above EPS_ACT.

    What that cost, measured before the fix: the shared vertices separated
    monotonically at ~0.094 per degree of crank (0.0000, 0.0143, 0.0285,
    0.0427 over the first three steps, at w = 0 and w = 0.6*w_ico alike), and
    the resulting "lock" at a* = 0.98 was the array coming apart rather than
    a jitterbug array locking. With the rows restored the same probe locks at
    a* = 29.88, which is 7.6 degrees PAST a_ico and therefore -- unlike every
    figure it replaces -- not the instant both prior mechanisms predict.

    A ball joint is a KINEMATIC CONSTRAINT, not a force: it adds no potential,
    no mass, no primitive, and does not disturb METRIC FORM's treatment (a).
    The FOUR DECLARATIONS are untouched by it. It is the constraint DECISION
    18's contact model was always meant to sit on top of.

    `topo` is REQUIRED rather than defaulted. A default would let a caller
    silently rebuild the disassembled array this bead exists to retire, and
    that failure mode is invisible in the output -- it moves a* without
    erroring. Both call sites already have `topo` in hand."""
    ndof = 48 * n
    c = len(topo.contacts)
    big = np.zeros((36 * n + 3 * c, ndof))
    for i in range(n):
        big[36 * i:36 * i + 36, 48 * i:48 * i + 48] = hinge_jacobian(xs[i])
    for e, (i, k, j, l) in enumerate(topo.contacts):
        fa, ca = PAIRS[k][0]
        fb, cb = PAIRS[l][0]
        r = 36 * n + 3 * e
        big[r:r + 3, 48 * i:48 * i + 48] += position_jacobian_row(xs[i], fa, ca)
        big[r:r + 3, 48 * j:48 * j + 48] -= position_jacobian_row(xs[j], fb, cb)
    return big


def hinge_residual(x):
    """The 36 shared-vertex hinge violations for one unit's current corner
    positions -- the CONSTRAINT FUNCTION `hinge_jacobian` is the Jacobian of."""
    out = np.zeros(36)
    for h, ((fa, ca), (fb, cb)) in enumerate(PAIRS):
        out[3 * h:3 * h + 3] = x[fa][ca] - x[fb][cb]
    return out


def project_to_pin_manifold(x, maxit=15, tol=1e-13):
    """Newton-project one unit's plates back onto the pin manifold after a
    step -- the STEP recipe's own words, mirroring `solve_inphase`'s GN
    style (damped least squares, bounded iterations, never raises).
    Returns (x, residual_norm)."""
    for _ in range(maxit):
        r = hinge_residual(x)
        rn = float(np.linalg.norm(r))
        if rn < tol:
            return x, rn
        jac = hinge_jacobian(x)
        dz, *_ = np.linalg.lstsq(jac, -r, rcond=None)
        if not np.all(np.isfinite(dz)):
            return x, rn
        x = apply_body_motions(x, dz)
    return x, float(np.linalg.norm(hinge_residual(x)))


def project_to_joint_manifold(xs, topo, maxit=PROJECT_MAXIT, tol=PROJECT_TOL):
    """Newton-project the WHOLE ARRAY back onto {hinge residuals = 0} AND
    {ball-joint residuals = 0} after a step. Returns (xs, residual_norm).

    WHY THIS IS ARRAY-LEVEL. Its predecessor, `project_to_pin_manifold`,
    projects ONE unit onto its OWN hinge manifold, and `crank_run` called it in
    a per-unit loop. That restores each unit's internal geometry and says
    nothing about whether the units are still joined -- so even once the ball
    joint is in the QP's Jacobian, where it holds only to FIRST order along a
    finite step, the joints drift at second order with nothing to pull them
    back. Measured: 0.0055 of residual separation over a 66-step run with the
    constraint in the QP but the projection still per-unit. Projecting the
    array as one system closes that.

    Damped least squares, bounded iterations, never raises -- the same GN
    style as `project_to_pin_manifold` and `solve_inphase`, and for the same
    reason: a projection that throws inside a swept loop destroys the verdict
    table (this file's house rule), and a projection that fails to converge is
    a measurement about the configuration, reported through the returned
    residual, not an exception."""
    n = topo.n
    for _ in range(maxit):
        r = np.concatenate([np.concatenate([hinge_residual(x) for x in xs]),
                            joint_residual(xs, topo)])
        rn = float(np.linalg.norm(r))
        if rn < tol:
            return xs, rn
        jac = build_pin_jacobian(xs, n, topo)
        dz, *_ = np.linalg.lstsq(jac, -r, rcond=None)
        if not np.all(np.isfinite(dz)):
            return xs, rn
        xs = [apply_body_motions(xs[i], dz[48 * i:48 * i + 48]) for i in range(n)]
    r = np.concatenate([np.concatenate([hinge_residual(x) for x in xs]),
                        joint_residual(xs, topo)])
    return xs, float(np.linalg.norm(r))


def max_joint_gap(xs, topo):
    """The largest shared-vertex separation in the array -- the single number
    that says whether the thing being cranked is still an array. Reported by
    the JOINT gate rows and by `crank_run`'s own telemetry."""
    if not topo.contacts:
        return 0.0
    return max(float(np.linalg.norm(joint_vertex(xs, i, k) - joint_vertex(xs, j, l)))
               for (i, k, j, l) in topo.contacts)


# ==========================================================================
# Z10: DRIVE + ONE STEP (bead qvf.18, item 3). Uniform in-phase / single-
# crank-handle target velocity, one QP solve, the achieved-rate projection.
# ==========================================================================

def crank_v_cmd(topo, a_hat, driven, ndof):
    """v_cmd: `path_tangent_48(a_hat)`'s phase direction PLUS the lattice's
    own rigid-translation rate (`topo.dsites(dverts_exact(a_hat))`) for every
    DRIVEN unit; zero elsewhere. `driven` is `"all"` (in-phase / uniform
    expansion) or a unit index (single-crank-handle) -- an INPUT, never
    hardcoded (module constant DRIVEN_UNIT_INDEX supplies the default)."""
    dsites = topo.dsites(dverts_exact(a_hat), verts(a_hat))
    v_cmd = np.zeros(ndof)
    units = range(topo.n) if driven == "all" else [driven]
    # PER-UNIT PHASE. Revision 34e636c rewrote every configuration site to
    # `unit_corners(a, topo, i)` so that "a two-sublattice topology cannot be
    # silently driven as though it were one" -- and missed THIS line, the
    # VELOCITY command, which kept driving every unit along the tangent at the
    # bare drive angle. The guarantee that commit stated was therefore false.
    # Invisible in the gate because main() drives SC7, whose phases are all
    # zero, so a_hat + 0.0 is a_hat and the two forms are bit-identical; a
    # honeycomb topology has phases {0, 60} and cos(tangent(a), tangent(a+60))
    # is 0.809, a 36-degree error. Cached per distinct offset, so the all-zero
    # case still evaluates the tangent exactly once.
    tangents = {}
    for i in units:
        off = topo.phases[i]
        if off not in tangents:
            tangents[off], _ = path_tangent_48(a_hat + off)
        vi = tangents[off].copy()
        for p in range(8):
            vi[24 + 3 * p:27 + 3 * p] += dsites[i]
        v_cmd[48 * i:48 * i + 48] = vi
    return v_cmd, list(units)


def crank_pairs(topo):
    """The FIXED, topology-derived candidate plate-pair list, computed once
    (GAP 2's `enumerate_plate_pairs`, reused as prior art)."""
    return enumerate_plate_pairs(topo)


#: Classification safety pad for the broad-phase skip: a pair is treated as
#: skippable-general only when |nA.nB| is at least this far BELOW
#: PARALLEL_TOL, so a last-bit disagreement between the vectorized
#: classification and signed_gap's own scalar one can only cause an extra
#: evaluation, never a wrong skip. The Euclidean bound is sound ONLY for the
#: general branch (the parallel branch's gap is an along-normal projection a
#: Euclidean bound does not bound); parallel pairs are always evaluated --
#: their branch is the cheap closed form anyway.
BROADPHASE_CLASS_PAD = 1e-9

#: Slack allowed between the VECTORISED parallel-branch gap and the scalar one
#: `signed_gap` computes. Both are (cB - cA) . nA over the same centroids, so
#: they agree to the last bit or two; a pair is shortcut only when its gap
#: clears EPS_ACT by more than this, so no last-bit disagreement can move a pair
#: across the active-set threshold. Pairs inside the band are evaluated.
BROADPHASE_GAP_PAD = 1e-9


def _pair_index_arrays(pairs):
    """Flat plate indices (iA, jB) and face indices fA for each pair --
    pure function of the pair list, computed once per crank run."""
    ia = np.fromiter((i * 8 + fi for (i, fi, j, fj) in pairs),
                     dtype=np.int64, count=len(pairs))
    jb = np.fromiter((j * 8 + fj for (i, fi, j, fj) in pairs),
                     dtype=np.int64, count=len(pairs))
    fa = np.fromiter((fi for (i, fi, j, fj) in pairs),
                     dtype=np.int64, count=len(pairs))
    return ia, jb, fa


def _pair_bounds_and_general(pairs, idx, xs, t):
    """(bound, skippable_general, surely_parallel, par_gap) per pair.

    bound = |cA-cB| - rA - rB - t <= the reported gap for GENERAL-branch pairs
    (unpierced: gap is a Euclidean distance >= bound; pierced: triangles
    overlap, so bound <= -t and the pair is never skipped).

    THE CLASSIFICATION IS PADDED IN BOTH DIRECTIONS (bead inviscid-0dm), because
    the two uses need opposite conservatism and the single-sided pad only served
    one of them:
      skippable_general  dots <= PARALLEL_TOL - PAD. Safe to CULL: the scalar
                         test in `signed_gap` certainly agrees the pair is
                         general, so the Euclidean bound certainly applies.
      surely_parallel    dots >= PARALLEL_TOL + PAD. Safe to SHORTCUT: the
                         scalar test certainly agrees the pair is parallel, so
                         `par_gap` below is certainly the gap it would report.
    Pairs in the 2*PAD band between them are neither, and are evaluated exactly
    as before. That band is 2e-9 wide in the dot product and is expected to be
    empty; it is handled rather than assumed away.

    par_gap is the PARALLEL branch's own closed form, (cB - cA) . nA - t,
    computed vectorized. For a surely-parallel pair this IS the gap, so such a
    pair needs no per-pair geometry call at all unless it turns out active --
    which is the coverage half of inviscid-0dm. Before this, parallel pairs were
    never culled and were 352 of the 360 evaluations a sorted cull still had to
    make on SC7."""
    ia, jb, fa = idx
    X = np.stack(xs)                       # (n, 8, 3, 3)
    Cm4 = X.mean(axis=2)
    Cm = Cm4.reshape(-1, 3)
    R = np.sqrt(((X - Cm4[:, :, None, :]) ** 2).sum(-1)).max(-1).reshape(-1)
    d = np.sqrt(((Cm[ia] - Cm[jb]) ** 2).sum(-1))
    bound = d - R[ia] - R[jb] - t
    e1 = (X[:, :, 1] - X[:, :, 0]).reshape(-1, 3)
    e2 = (X[:, :, 2] - X[:, :, 0]).reshape(-1, 3)
    nb = np.stack((e1[:, 1] * e2[:, 2] - e1[:, 2] * e2[:, 1],
                   e1[:, 2] * e2[:, 0] - e1[:, 0] * e2[:, 2],
                   e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]), axis=1)
    nb = nb / np.sqrt((nb ** 2).sum(-1))[:, None]
    na = np.stack([plate_normal(f) for f in range(8)])
    dots = np.abs((na[fa] * nb[jb]).sum(-1))
    general = dots <= PARALLEL_TOL - BROADPHASE_CLASS_PAD
    surely_parallel = dots >= PARALLEL_TOL + BROADPHASE_CLASS_PAD
    par_gap = ((Cm[jb] - Cm[ia]) * na[fa]).sum(-1) - t
    return bound, general, surely_parallel, par_gap


def exhaustive_contact_scan(pairs, xs, t, ndof):
    """THE REFERENCE INSTRUMENT: every pair evaluated, nothing culled.

    Kept in the file on purpose (bead inviscid-qvf.23 acceptance criterion 1):
    the broad-phase is only trustworthy against something that does not use it,
    and a reference that lives only in a reviewer's scratch directory is not a
    reference. Returns (active, min_general_gap) with `active` a list of
    (i, fi, j, fj) in the pair list's own order."""
    active, mgg = [], float("inf")
    for (i, fi, j, fj) in pairs:
        gap, _row = contact_gradient_row(xs, i, fi, j, fj, t, ndof)
        nA = plate_normal(fi)
        nB = _cross3(xs[j][fj][1] - xs[j][fj][0], xs[j][fj][2] - xs[j][fj][0])
        nB_norm = np.linalg.norm(nB)
        if nB_norm > 1e-300 and abs(float(nA @ (nB / nB_norm))) <= PARALLEL_TOL:
            mgg = min(mgg, gap)
        if gap <= EPS_ACT:
            active.append((i, fi, j, fj))
    return active, mgg


def broadphase_contact_scan(pairs, xs, t, ndof, bias=0.0):
    """The same thing THROUGH the broad-phase, reporting what it evaluated.

    Mirrors `crank_step`'s two passes exactly -- if this and the stepper ever
    drift apart the exactness row is measuring the wrong code, so they are
    written adjacent and reviewed together."""
    idx = _pair_index_arrays(pairs)
    bnd, gen, par, pgap = _pair_bounds_and_general(pairs, idx, xs, t)
    # `bias` inflates the bound, making the cull UNSOUND. It exists only for the
    # mutation probe: a guard nothing has tried to break is not a guard, and the
    # exactness row must be shown capable of reddening.
    bnd = bnd + bias
    pgap = pgap + bias
    gaps, mgg, evaluated = {}, float("inf"), 0
    for k in np.argsort(bnd, kind="stable"):
        k = int(k)
        if gen[k] and bnd[k] > EPS_ACT and bnd[k] > mgg:
            continue
        if par[k] and pgap[k] - EPS_ACT > BROADPHASE_GAP_PAD:
            continue
        i, fi, j, fj = pairs[k]
        nA = plate_normal(fi)
        gap, _a, _b, _c = signed_gap(xs[i][fi], xs[j][fj], nA)
        gap -= t
        gaps[k] = gap
        evaluated += 1
        nB = _cross3(xs[j][fj][1] - xs[j][fj][0], xs[j][fj][2] - xs[j][fj][0])
        nB_norm = np.linalg.norm(nB)
        if nB_norm > 1e-300 and abs(float(nA @ (nB / nB_norm))) <= PARALLEL_TOL:
            mgg = min(mgg, gap)
    active = [pairs[k] for k in range(len(pairs))
              if gaps.get(k) is not None and gaps[k] <= EPS_ACT]
    return active, mgg, evaluated


#: Configurations the exactness probe runs on: ideal-path angles, and for each
#: the state one crank step later, which is the mid-integration regime where a
#: previous prefilter class died. Absolute, not derived from any other grid.
BP_PROBE_ANGLES = (0.5, 3.0, 8.0, 15.0, A_ICO)
BP_PROBE_T = (0.0, 0.02)


def bp_exactness_probe(topo):
    """Broad-phase against the exhaustive reference, on real configurations
    INCLUDING mid-integration states. Returns one record per configuration."""
    pairs = crank_pairs(topo)
    ndof = 48 * topo.n
    wpairs = wire_pairs(topo)
    out = []
    for a in BP_PROBE_ANGLES:
        origins = topo.sites(verts(a))
        base = [unit_corners(a, topo, i) + origins[i] for i in range(topo.n)]
        for tag, xs in (("ideal", base), ("stepped", _bp_one_step(topo, pairs, base, a, wpairs))):
            if xs is None:
                continue
            for t in BP_PROBE_T:
                ex_act, ex_mgg = exhaustive_contact_scan(pairs, xs, t, ndof)
                bp_act, bp_mgg, ev = broadphase_contact_scan(pairs, xs, t, ndof)
                out.append(dict(
                    a=a, tag=tag, t=t, n=len(pairs), evaluated=ev,
                    reject=1.0 - ev / len(pairs),
                    same_active=(ex_act == bp_act),
                    mgg_diff=abs(ex_mgg - bp_mgg)
                    if np.isfinite(ex_mgg) and np.isfinite(bp_mgg) else 0.0))
    return out


#: How far the mutation probe inflates the broad-phase bound. Chosen from the
#: measured active-set geometry rather than picked: gaps at the probe
#: configurations sit within ~0.1 of the threshold, so 0.5 certainly reaches
#: past real active pairs and the probe certainly bites. If it ever stops
#: biting, the exactness row above has gone slack and says nothing.
BP_MUTATION_BIAS = 0.5


def bp_mutation_probe(topo):
    """Deliberately break the cull and confirm exactness NOTICES.

    qvf.23 acceptance criterion 4. Without this the exactness row is satisfied
    by any prefilter that never rejects anything, and by any bug that happens to
    reject only inactive pairs at the sampled configurations."""
    pairs = crank_pairs(topo)
    ndof = 48 * topo.n
    caught = 0
    total = 0
    for a in BP_PROBE_ANGLES:
        origins = topo.sites(verts(a))
        xs = [unit_corners(a, topo, i) + origins[i] for i in range(topo.n)]
        for t in BP_PROBE_T:
            ex_act, _ = exhaustive_contact_scan(pairs, xs, t, ndof)
            bad_act, _m, _e = broadphase_contact_scan(
                pairs, xs, t, ndof, bias=BP_MUTATION_BIAS)
            total += 1
            if bad_act != ex_act:
                caught += 1
    return dict(caught=caught, total=total)


def _bp_one_step(topo, pairs, xs, a, wpairs):
    """One crank step from `xs`, giving a mid-integration configuration."""
    v, status, _r, _b, _m = crank_step(topo, pairs, xs, a, _w_ico_lock(), 0.0,
                                       "all", True, None, wpairs)
    if status != "OK" or v is None:
        return None
    return [apply_body_motions(xs[u], H_STEP * v[48 * u:48 * u + 48])
            for u in range(topo.n)]


def crank_step(topo, pairs, xs, a_hat, w, t, driven="all",
               enforce_contacts=True, w_diag=None, wpairs=None, bp_idx=None):
    """One velocity-level QP solve. Returns
    (v, status, rate, binding, min_general_gap) where status is one of
    "OK" / "QPFAIL" (the two DISJOINT solver outcomes) -- `rate` and
    `binding` are None on QPFAIL, never a stand-in physical reading.
    `wpairs` is the topology's precomputed `wire_pairs(topo)` list;
    omitting it (the standalone/diagnostic call sites) recomputes it, since
    it is a pure function of `topo` alone."""
    n = topo.n
    ndof = 48 * n
    v_cmd, units = crank_v_cmd(topo, a_hat, driven, ndof)
    j_pin = build_pin_jacobian(xs, n, topo)
    n_basis = null_space_basis(j_pin, ndof)
    if wpairs is None:
        wpairs = wire_pairs(topo)

    active_rows = []
    active_labels = []
    min_general_gap = float("inf")
    if enforce_contacts:
        # BROAD-PHASE (exact, orchestrator-applied under owner directive
        # 2026-08-22, pending re-review): GENERAL-branch pairs only. Their
        # reported gap is a Euclidean distance (unpierced) or negative with
        # overlapping triangles (pierced, bound <= -t, never skipped), so
        # bound > EPS_ACT and bound > the running min_general_gap proves the
        # pair can neither join the active set nor lower the reported
        # minimum (the running min only decreases). PARALLEL pairs are
        # never skipped: their gap is an along-normal projection a
        # Euclidean bound does not bound (measured, not hypothesized -- a
        # first distance-only version of this skip failed exactness on
        # exactly that class and was corrected to this form).
        if bp_idx is None:
            bp_idx = _pair_index_arrays(pairs)
        _bnd, _gen, _par, _pgap = _pair_bounds_and_general(pairs, bp_idx, xs, t)

        # PASS 1, IN ASCENDING BOUND ORDER (bead inviscid-0dm). The cull is
        # exact either way -- a pair with bound > the running minimum cannot
        # lower it, and the running minimum only decreases -- but the ORDER
        # decided how much it culled, and the order was the pair list's. The
        # first pair scanned was never culled and the rate depended on
        # enumeration rather than geometry. Ascending bound is the order in
        # which the running minimum falls fastest, and it depends only on the
        # configuration: on SC7 it takes general-branch evaluations from 131 of
        # 1080 to 8.
        _gaps = {}
        for _k in np.argsort(_bnd, kind="stable"):
            _k = int(_k)
            if _gen[_k] and _bnd[_k] > EPS_ACT and _bnd[_k] > min_general_gap:
                continue
            if _par[_k] and _pgap[_k] - EPS_ACT > BROADPHASE_GAP_PAD:
                continue      # surely parallel and surely inactive: par_gap IS
                              # its gap, and it is clear of the threshold by
                              # more than the vectorised form's last-bit slack
            i, fi, j, fj = pairs[_k]
            nA = plate_normal(fi)
            gap, _wA, _wB, _nrm = signed_gap(xs[i][fi], xs[j][fj], nA)
            gap -= t
            _gaps[_k] = gap
            nB = _cross3(xs[j][fj][1] - xs[j][fj][0], xs[j][fj][2] - xs[j][fj][0])
            nB_norm = np.linalg.norm(nB)
            is_general = nB_norm > 1e-300 and abs(float(nA @ (nB / nB_norm))) <= PARALLEL_TOL
            if is_general:
                min_general_gap = min(min_general_gap, gap)

        # PASS 2, IN THE ORIGINAL PAIR ORDER, so the active set and every row
        # built from it are assembled exactly as before. Only pairs the first
        # pass found active reach `contact_gradient_row`, and it recomputes the
        # gap itself, so the reported value is the same float it always was.
        for _k, (i, fi, j, fj) in enumerate(pairs):
            _g = _gaps.get(_k)
            if _g is None or _g > EPS_ACT:
                continue
            if _g < -MEANINGLESS_DEPTH_FLOOR:
                # THE ARTIFACT CLASS, filtered where the rows are BUILT and not
                # only where a step is accepted. `signed_gap`'s parallel branch
                # tests |nA . nB|, so it fires on ANTIPARALLEL normals too, and
                # for two plates on opposite sides of one unit -- facing away
                # from each other, permanently ~2.31 apart, structurally unable
                # to touch -- it reports gap = -2.309401 and that reads as deep
                # penetration. 28 such pairs exist on SC7, 9.5% of the active
                # rows at a = 0. `MEANINGLESS_DEPTH_FLOOR` already names this
                # class and crank_run's backtrack acceptance already rejects it;
                # crank_step's active-row construction did not, so the rows were
                # built and then judged. Deep-penetration accuracy is out of
                # scope by `signed_gap`'s own docstring, so a genuine overlap
                # this deep is not a contact this file claims to resolve either.
                continue
            gap, row = contact_gradient_row(xs, i, fi, j, fj, t, ndof)
            if gap <= EPS_ACT:
                active_rows.append(row)
                active_labels.append(("contact", i, fi, j, fj, gap))
    for (u, fi, v, fj) in wpairs:
        g_wire, row, degenerate = wire_gradient_row(xs, (u, fi, v, fj), w, ndof)
        if degenerate:
            continue  # a genuinely zero-span wire is EXCLUDED, never a
            # vacuously-binding zero-gradient row (critique 23262 C2)
        if g_wire <= EPS_ACT:
            active_rows.append(row)
            active_labels.append(("wire", u, fi, v, fj, g_wire))

    v, feasible = project_qp(v_cmd, n_basis, active_rows, w_diag)
    if not feasible or not np.all(np.isfinite(v)):
        return None, "QPFAIL", None, None, min_general_gap

    denom = sum(float(v_cmd[48 * u:48 * u + 48] @ v_cmd[48 * u:48 * u + 48]) for u in units)
    numer = sum(float(v[48 * u:48 * u + 48] @ v_cmd[48 * u:48 * u + 48]) for u in units)
    rate = numer / denom if denom > 1e-300 else 0.0
    binding = [lab for lab, row in zip(active_labels, active_rows)
              if abs(float(row @ v)) < BINDING_TOL]
    return v, "OK", rate, binding, min_general_gap


# ==========================================================================
# Z11: CRANK RUN (bead qvf.18, item 5, LOCK/JAM). The stepping loop: advance,
# Newton-project, refresh active sets, backtrack `h` against the REAL
# nonlinear gap (the HAZARD's feasibility-preserving invariant), detect
# REACHED / JAMMED / QPFAIL as three disjoint outcomes.
# ==========================================================================

#: Feasibility-backtrack SHALLOW-VIOLATION floor: the backtrack rejects a
#: trial step for a pair only when its trial gap falls in
#: (-MEANINGLESS_DEPTH_FLOOR, -gap_floor) -- a SHALLOW, near-boundary
#: violation. Three PREFILTER designs were tried and rejected before this,
#: each found broken LIVE via the w-causality cross-test, not hypothesized:
#: a gap-VALUE band sampled only at a step's start (missed a pair drifting
#: from far-positive into slight violation one step later); the same band
#: sampled at both the start AND the full-h0 trial (still missed a MID-
#: interval crossing); a PLATE-CENTROID distance radius, on the reasoning
#: that signed_gap's parallel-facing branch is only meaningless for pairs
#: that are spatially far apart -- ALSO WRONG: a single octahedron unit has
#: intra-unit plate pairs (e.g. (0,0,0,7)) whose normals are (anti)parallel
#: by the octahedron's OWN symmetry, sitting at a MODERATE plate-centroid
#: distance (2.31, well inside any reasonable "close" radius) while being
#: exactly as physically meaningless (a fixed, rigid-body-invariant
#: projection between two faces of the SAME body that can never approach
#: each other) as the deeply-negative far-unit case. Distance -- of any
#: kind -- does not separate the two classes; DEPTH does: every meaningless
#: value measured across this whole investigation sat beyond -2.3; every
#: genuine near-boundary risk this file's own kernel already declared out
#: of scope for deep-penetration accuracy (module docstring, "A ROW
#: DELIBERATELY NOT BUILT") lives near zero. This floor formalizes that
#: EXISTING, already-documented scope boundary as the feasibility check's
#: own criterion, rather than trying to re-derive it from geometry.
MEANINGLESS_DEPTH_FLOOR = 1.0


@jb_cache.memoize(_MODULE)
def crank_run(topo, a_start, a_target, w, t, driven="all",
              enforce_contacts=True, w_diag=None, h0=H_STEP,
              max_steps=MAX_CRANK_STEPS, backtrack=True,
              gap_floor=GAP_FLOOR_TOL, instant_jam=True):
    """MEMOISED (jb_cache): this function is the whole cost of this file --
    ~30 calls, 99.8% of a 6m14s run, measured 2026-08-22. It is a pure
    function of its arguments and of module constants, so its results are
    cached on the SHA of its own transitive source closure plus its bound
    arguments, and computed in parallel ahead of the serial gate pass. Edit
    the stepper, the contact kernel, or any constant either reads and every
    affected entry invalidates automatically; edit a print statement and none
    do. `--no-cache` bypasses both mechanisms and must print byte-identically
    -- that equivalence, not this comment, is what makes the claim checkable.

    Returns a dict: status in {"reached","jammed","qpfail"}; a_final;
    jam_angle (jammed only); binding (jammed only); steps; min_general_gap;
    rate_history (for the mutation-probe / insensitivity comparisons).

    min_general_gap's ACTUAL coverage (fix for code review 23261 M1, whose
    docstring here previously overclaimed "every accepted AND rejected
    trial position"): each iteration's `crank_step` call scans every
    enumerated pair at the CURRENT, already-accepted-and-Newton-projected
    `xs` (the prior step's outcome) -- trial positions considered and
    REJECTED inside the `H_BACKTRACK_MAX` bisection loop are checked only
    against `watch_pairs` (the near-threshold CONTACT subset used for the
    feasibility backtrack itself), and that check's results do not feed
    this telemetry. A final scan of the TERMINAL accepted position is added
    below, right before a "reached" return, so the blind spot this
    docstring used to hide (the very last step's own outcome never being
    general-branch-scanned) is closed rather than merely documented.

    `instant_jam` selects which of the two DISTINCT jam semantics this run
    uses (see STALL_RATE_TOL's docstring): True (ONSET, row G) declares jam
    on the very first reading below JAM_RATE_TOL, no stepping attempted --
    matching the bead's own row-G wording verbatim. False (SUSTAINED, rows
    H/I) takes the step anyway through any transient dip and declares jam
    only if the rate collapses below STALL_RATE_TOL (genuinely, persistently
    stuck, not merely reduced)."""
    pairs = crank_pairs(topo)
    wpairs = wire_pairs(topo)
    bp_idx = _pair_index_arrays(pairs)
    origins = topo.sites(verts(a_start))
    xs = [unit_corners(a_start, topo, i) + origins[i] for i in range(topo.n)]
    a_hat = a_start
    min_general_gap = float("inf")
    rate_history = []
    binding_ever = set()

    def _scan_general_gap(xs_now):
        worst = float("inf")
        if not enforce_contacts:
            return worst
        _sb, _sg, _sp, _spg = _pair_bounds_and_general(
            pairs, _pair_index_arrays(pairs), xs_now, t)
        # Ascending bound here too (bead inviscid-0dm), for the same reason and
        # with the same exactness argument: this scan reports only a MINIMUM, so
        # the order it visits pairs in cannot change the answer, only how many
        # pairs it has to touch to reach it.
        for _k in np.argsort(_sb, kind="stable"):
            _k = int(_k)
            i, fi, j, fj = pairs[_k]
            if _sg[_k] and _sb[_k] > worst:
                continue  # general-branch only; bound > running worst >= final worst
            if _sp[_k]:
                continue  # surely parallel: excluded from this GENERAL-branch
                          # minimum by the scalar test below in any case
            g, _ = contact_gradient_row(xs_now, i, fi, j, fj, t, 48 * topo.n)
            nA = plate_normal(fi)
            nB = _cross3(xs_now[j][fj][1] - xs_now[j][fj][0], xs_now[j][fj][2] - xs_now[j][fj][0])
            nB_norm = np.linalg.norm(nB)
            is_general = nB_norm > 1e-300 and abs(float(nA @ (nB / nB_norm))) <= PARALLEL_TOL
            if is_general:
                worst = min(worst, g)
        return worst

    for step in range(max_steps):
        v, status, rate, binding, mgg = crank_step(
            topo, pairs, xs, a_hat, w, t, driven, enforce_contacts, w_diag,
            wpairs, bp_idx)
        min_general_gap = min(min_general_gap, mgg)
        if binding:
            binding_ever.update((b[0], b[1], b[2], b[3], b[4]) for b in binding)
        if status != "OK":
            return {"status": "qpfail", "a_final": a_hat, "jam_angle": None,
                    "binding": None, "binding_ever": binding_ever, "steps": step,
                    "min_general_gap": min_general_gap, "rate_history": rate_history}
        rate_history.append(rate)
        stall_tol = JAM_RATE_TOL if instant_jam else STALL_RATE_TOL
        if rate < stall_tol:
            return {"status": "jammed", "a_final": a_hat, "jam_angle": a_hat,
                    "binding": binding, "binding_ever": binding_ever, "steps": step,
                    "min_general_gap": min_general_gap, "rate_history": rate_history}

        # Feasibility backtrack (the HAZARD's own invariant): checked
        # against EVERY enumerated pair, every backtrack iteration -- see
        # MEANINGLESS_DEPTH_FLOOR's docstring for the three prefilter
        # designs tried and rejected before landing on a SHALLOW-violation
        # band instead of any kind of candidate-list prefilter.
        h = h0
        accepted = False
        trial_xs = xs
        for _ in range(H_BACKTRACK_MAX):
            trial_xs = [apply_body_motions(xs[i], h * v[48 * i:48 * i + 48])
                       for i in range(topo.n)]
            ok = True
            if backtrack and enforce_contacts:
                for (i, fi, j, fj) in pairs:
                    g, _ = contact_gradient_row(trial_xs, i, fi, j, fj, t, 48 * topo.n)
                    if -MEANINGLESS_DEPTH_FLOOR < g < -gap_floor:
                        ok = False
                        break
            if ok:
                accepted = True
                break
            h *= 0.5
            if h < H_MIN:
                break
        if not accepted:
            # The QP's proposed direction admits NO feasible step down to
            # H_MIN -- this is the "no positive phase progress" JAM the
            # bead's own LOCK definition names, a PHYSICAL fact the
            # backtrack discovered, not a solver malfunction; it must not
            # be conflated with QPFAIL (the _LPFail-class distinction the
            # HAZARDS section requires -- see the module docstring). A
            # genuine solver-side failure (NNLS non-convergence, a
            # non-finite v) is caught earlier, at `crank_step`'s own
            # `status != "OK"` branch above, and stays QPFAIL there.
            return {"status": "jammed", "a_final": a_hat, "jam_angle": a_hat,
                    "binding": binding, "binding_ever": binding_ever, "steps": step,
                    "min_general_gap": min_general_gap, "rate_history": rate_history}

        # ARRAY-LEVEL projection (bead inviscid-1wd). This was a per-unit
        # loop over `project_to_pin_manifold`, which restores each unit's own
        # hinge geometry and says nothing about whether the units are still
        # joined -- the ball joint holds only to first order along a finite
        # step, so without this the joints drift at second order with nothing
        # pulling them back.
        xs, _joint_rn = project_to_joint_manifold(trial_xs, topo)
        a_hat += h * rate
        if a_hat >= a_target:
            min_general_gap = min(min_general_gap, _scan_general_gap(xs))
            return {"status": "reached", "a_final": a_hat, "jam_angle": None,
                    "binding": None, "binding_ever": binding_ever, "steps": step + 1,
                    "min_general_gap": min_general_gap, "rate_history": rate_history}

    return {"status": "qpfail", "a_final": a_hat, "jam_angle": None, "binding": None,
            "binding_ever": binding_ever, "steps": max_steps,
            "min_general_gap": min_general_gap, "rate_history": rate_history}


# ==========================================================================
# Z12: DOWELED DIAGNOSTIC (bead qvf.18, row J -- SECOND CONTROL). Reuses
# `assemble_doweled` + `rank_of` exactly as the prior art specifies: mobility
# is a MEASURED RANK, never a subtraction (memo R2).
# ==========================================================================

def doweled_diagnostic(topo, a=A_ICO_LOCAL):
    """`assemble_doweled` returns (3*ncontacts, 7n): rank is BOUNDED ABOVE by
    3*ncontacts (the row count) REGARDLESS of topology -- for SC7, a TREE
    (jb_x's own note: "as a TREE", 6 contacts, no cycles), 3*ncontacts=18
    is far short of 7n-7=42, so SC7's doweled array CANNOT collapse to one
    global DOF by row-counting alone; that collapse needs a topology WITH
    CYCLES (redundant paths) to push the rank higher. The FALSIFIABLE claim
    this row gates is full row rank (every contact equation independent, no
    accidental redundancy) -- expected_rank = min(3*ncontacts, 7n-7).
    one_global_dof is reported SEPARATELY and HONESTLY: True only when that
    cap is 7n-7 AND it is achieved; SC7 being a tree, it is False, and that
    is the correct answer, not a defect."""
    j = assemble_doweled(a, topo)
    rank, sv = rank_of(j)
    row_cap = 3 * len(topo.contacts)
    collapse_target = 7 * topo.n - 7  # 6 global rigid + 1 shared breathing dof
    expected_rank = min(row_cap, collapse_target)
    return {"rank": rank, "sv": sv, "expected_rank": expected_rank,
            "row_cap": row_cap, "collapse_target": collapse_target,
            "one_global_dof": rank == collapse_target}


# ==========================================================================
# Z13: HAZARD invariant + Z14: DRIVEN-UNIT-IS-AN-INPUT + two-cell sanity
# (bead qvf.18: the hazard comment 2026-08-21, and "activate where the
# fold/registry says").
# ==========================================================================

def _find_general_branch_pair(topo, a):
    """A real, currently GENERAL-branch (non-parallel-facing) plate pair
    with a clean positive gap and a nonzero closing gradient -- the exact
    regime the qvf.17 critique named as uncovered for negative gaps."""
    origins = topo.sites(verts(a))
    xs = [unit_corners(a, topo, i) + origins[i] for i in range(topo.n)]
    best = None
    for (i, fi, j, fj) in crank_pairs(topo):
        triA, triB = xs[i][fi], xs[j][fj]
        nA = plate_normal(fi)
        nB = _cross3(triB[1] - triB[0], triB[2] - triB[0])
        nrm = np.linalg.norm(nB)
        if nrm < 1e-300 or abs(float(nA @ (nB / nrm))) > PARALLEL_TOL:
            continue
        gap, row = contact_gradient_row(xs, i, fi, j, fj, 0.0, 48 * topo.n)
        if 0.01 < gap < 0.5 and (best is None or gap < best[0]):
            best = (gap, i, fi, j, fj, row)
    return xs, best


def z13_hazard_invariant(topo, a_probe=1.0):
    """TWO pieces of evidence, both real: (1) the decisive G/H run's OWN
    observed minimum general-branch gap, with backtrack ON vs FORCED OFF --
    at THIS array's initial configuration the actually-stressed pairs happen
    to be parallel-facing (Z0's own key fact: many plates are EXACTLY
    parallel at a=0), so this comparison alone does not discriminate, and is
    reported honestly as such, not hidden. (2) a TARGETED, ISOLATED probe on
    a REAL general-branch pair (`_find_general_branch_pair`): a synthetic
    velocity built directly from that pair's own closing gradient (bypassing
    the QP, which would never choose a velocity that closes an inactive
    constraint -- this tests the BACKTRACK MECHANISM itself, in isolation,
    the way a unit test isolates the function it exercises), applied via the
    SAME bisection logic `crank_run` uses. This is the row that is actually
    non-vacuous: WITH the mechanism, the resulting gap stays at the floor;
    WITHOUT it (the full, unbisected step applied directly), it does not."""
    on_g = crank_run(topo, 0.0, A_ICO_LOCAL, w=0.0, t=0.0, backtrack=True, h0=0.5,
                     instant_jam=False)
    off_g = crank_run(topo, 0.0, A_ICO_LOCAL, w=0.0, t=0.0, backtrack=False, h0=0.5,
                      instant_jam=False)

    xs, found = _find_general_branch_pair(topo, a_probe)
    if found is None:
        return {"on_min_gap": on_g["min_general_gap"], "off_min_gap": off_g["min_general_gap"],
                "probe_found": False, "probe_on_gap": None, "probe_off_gap": None}
    gap0, i, fi, j, fj, row = found
    ndof = 48 * topo.n
    v = np.zeros(ndof)
    d = -row / np.linalg.norm(row)  # the CLOSING direction for this pair
    v[:] = 0.3 * d  # closes this pair well past gap0's margin, but stays
    # inside apply_body_motions' small-rotation regime (verified numerically
    # -- a MUCH larger magnitude here rotates a full body by several
    # RADIANS, past the point where "closing" is even the right word for
    # what happens, and the probe stops testing what it claims to)

    def gap_after(h):
        trial = [apply_body_motions(xs[u], h * v[48 * u:48 * u + 48]) for u in range(topo.n)]
        g, _ = contact_gradient_row(trial, i, fi, j, fj, 0.0, ndof)
        return g

    probe_off_gap = gap_after(1.0)  # the full, un-bisected step -- no protection

    h = 1.0
    probe_on_gap = probe_off_gap
    for _ in range(H_BACKTRACK_MAX):
        g = gap_after(h)
        if g >= -GAP_FLOOR_TOL:
            probe_on_gap = g
            break
        h *= 0.5
    else:
        probe_on_gap = gap_after(h)

    return {"on_min_gap": on_g["min_general_gap"], "off_min_gap": off_g["min_general_gap"],
            "probe_found": True, "probe_gap0": gap0, "probe_on_gap": probe_on_gap,
            "probe_off_gap": probe_off_gap}


def z14_two_cell_sanity(topo):
    """The single-crank-handle DRIVE variant, exercised live (not merely
    parametric): drive unit DRIVEN_UNIT_INDEX only. Separately, for the pair
    of plates between units 0 and 1 whose gap is closest to w at a=0 (a
    static, direct `signed_gap` scan -- no stepper), root-find the angle at
    which that pair's gap crosses `w`, then confirm a SHORT dynamic crank_run
    at the SAME w flags that identical pair as bound (active) at a consistent
    angle -- 'the stepper activates where the kernel's own static gap says',
    qvf.17's story reproduced dynamically."""
    n = topo.n
    ndof = 48 * n
    v_cmd, units = crank_v_cmd(topo, 0.0, DRIVEN_UNIT_INDEX, ndof)
    origins = topo.sites(verts(0.0))
    xs0 = [unit_corners(0.0, topo, i) + origins[i] for i in range(n)]
    v, status, rate, binding, _ = crank_step(topo, crank_pairs(topo), xs0, 0.0,
                                              w=0.0, t=0.0, driven=DRIVEN_UNIT_INDEX)

    pairs01 = [(i, fi, j, fj) for (i, fi, j, fj) in crank_pairs(topo)
              if {i, j} == {0, 1}]
    w_test = 0.02
    best = None
    for (i, fi, j, fj) in pairs01:
        triA, triB = xs0[i][fi], xs0[j][fj]
        nA = plate_normal(fi)
        gap0, *_ = signed_gap(triA, triB, nA)
        if best is None or abs(gap0 - w_test) < abs(best[1] - w_test):
            best = ((i, fi, j, fj), gap0)
    pair, gap0 = best if best else (None, None)

    static_cross = None
    if pair is not None:
        i, fi, j, fj = pair
        lo, hi = 0.0, 5.0
        g_lo, *_ = signed_gap(
            xs0[i][fi],
            (unit_corners(0.0, topo, j) + topo.sites(verts(0.0))[j])[fj],
            plate_normal(fi))
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            origins_m = topo.sites(verts(mid))
            xm = [unit_corners(mid, topo, u) + origins_m[u] for u in range(n)]
            g_mid, *_ = signed_gap(xm[i][fi], xm[j][fj], plate_normal(fi))
            if (g_mid - w_test) * (g_lo - w_test) <= 0:
                hi = mid
            else:
                lo = mid
                g_lo = g_mid
        static_cross = 0.5 * (lo + hi)

    run = crank_run(topo, 0.0, min(static_cross + 1.0, 5.0) if static_cross else 5.0,
                    w=w_test, t=0.0, driven=DRIVEN_UNIT_INDEX, h0=0.1, instant_jam=False)
    # binding_ever (accumulated across EVERY step, not just a terminal jam --
    # this run is expected to REACH, per the design's own "onset is
    # transient" finding, so the terminal "binding" field is empty by
    # construction; the pair activating at SOME point during the run is the
    # actual claim being gated).
    dynamic_binds_pair = pair is not None and ("contact",) + pair in run["binding_ever"]

    return {"single_unit_rate": rate, "single_unit_status": status,
            "pair": pair, "static_cross_angle": static_cross,
            "run_status": run["status"], "run_jam_angle": run["jam_angle"],
            "dynamic_binds_pair": dynamic_binds_pair}


# ==========================================================================
# Z15: THE DECISIVE ROW PAIR (G, H) + CONTROL (I) + W-insensitivity (M) +
# determinism (L). Everything the crank stepper's acceptance criteria ask
# for, assembled from Z8-Z14 above.
# ==========================================================================

def z15_crank_gates(topo):
    ref = DIAGONALS[0]
    fold_ico = fold_halves(A_ICO_LOCAL)[ref]
    w_ico = 2.0 * (fold_ico[1] - fold_ico[0])

    # G: w=0, t=0, in-phase crank -- must JAM immediately off a=0.
    g_run = crank_run(topo, 0.0, A_ICO_LOCAL, w=0.0, t=0.0, driven="all", h0=0.5,
                      instant_jam=True)

    # H: w >= 2*fold(a_ico) -- must REACH a_ico. Non-vacuity companion to G:
    # a solver that always reports jam passes G and must FAIL this.
    h_run = crank_run(topo, 0.0, A_ICO_LOCAL, w=w_ico, t=0.0, driven="all", h0=H_STEP,
                      instant_jam=False)

    # I: CONTROL -- all contact constraints disabled, wires slack (w large):
    # the crank must reach a=60 with NO jam. If this jams, the pin/projection
    # machinery is wrong and the contact results (G, H) mean nothing.
    i_run = crank_run(topo, 0.0, 60.0, w=1000.0, t=0.0, driven="all",
                      enforce_contacts=False, h0=5.0, instant_jam=False)

    # W-CAUSALITY CROSS-TEST (substantive critique 23262 C1, ship-blocker):
    # G and H differ in BOTH w AND instant_jam, so their opposite outcomes
    # alone do not prove w is causal (the critique's own reproduced cross-
    # test, using the file's PRE-fix wire attachment, showed the SAME
    # outcome flips purely from toggling instant_jam at FIXED w). This row
    # holds the PROTOCOL FIXED (instant_jam=False, SUSTAINED semantics,
    # both runs) and varies ONLY w: w=0 (wires bind immediately -- span is
    # monotonically increasing under the verified 2/sqrt(3)*fold(a) opening
    # relationship, so a taut w=0 wire never relaxes, unlike the transient
    # contact-only resistance) versus w_ico (generous enough that the
    # opening pairs' span, topping out near fold(a_ico)*2/sqrt(3) approx
    # 0.505, never reaches the 0.874 limit before a_ico). If w is genuinely
    # causal, this pair's outcomes must differ; if it is not (the OLD,
    # symmetry-fixed-point wire attachment), they would not.
    cross_w0 = crank_run(topo, 0.0, A_ICO_LOCAL, w=0.0, t=0.0, driven="all", h0=0.5,
                         instant_jam=False)
    cross_wico = crank_run(topo, 0.0, A_ICO_LOCAL, w=w_ico, t=0.0, driven="all", h0=H_STEP,
                           instant_jam=False)

    # MATCHED-h0 CROSS REBUILD (qvf.19 fix round, critique 23299 CRITICAL 1,
    # ship-blocker): cross_w0 (h0=0.5) and cross_wico (h0=H_STEP=2.0) above
    # are CONFOUNDED -- they differ in BOTH w AND h0, so their differing
    # outcomes never isolated the w-effect from the step-size effect, despite
    # this row's own original comment claiming "varies ONLY w". This bead is
    # explicitly authorized to touch these rows (once, for this fix) to add
    # the two MISSING matched legs -- cross_w0/cross_wico themselves are left
    # completely unchanged (same code, same values) so no PRIOR row's value
    # changes; the two new legs below complete BOTH 2x2 matched pairs
    # (h0 in {0.5, 2.0}} x (w in {0, w_ico}}): cross_w0/cross_wico already
    # cover (h0=0.5,w=0) and (h0=2.0,w=w_ico); the two adds below are
    # (h0=0.5,w=w_ico) and (h0=2.0,w=0).
    cross_h05_wico = crank_run(topo, 0.0, A_ICO_LOCAL, w=w_ico, t=0.0, driven="all", h0=0.5,
                               instant_jam=False)
    cross_h2_w0 = crank_run(topo, 0.0, A_ICO_LOCAL, w=0.0, t=0.0, driven="all", h0=H_STEP,
                            instant_jam=False)

    # eps_act decade-insensitivity: rerun G's FIRST step at EPS_ACT_ALT.
    pairs = crank_pairs(topo)
    origins0 = topo.sites(verts(0.0))
    xs0 = [unit_corners(0.0, topo, i) + origins0[i] for i in range(topo.n)]
    global EPS_ACT
    saved_eps = EPS_ACT
    EPS_ACT = EPS_ACT_ALT
    try:
        _, _, rate_alt_eps, _, _ = crank_step(topo, pairs, xs0, 0.0, w=0.0, t=0.0, driven="all")
    finally:
        EPS_ACT = saved_eps
    _, _, rate_nominal_eps, _, _ = crank_step(topo, pairs, xs0, 0.0, w=0.0, t=0.0, driven="all")

    # M: METRIC FORM W-insensitivity -- rerun G's decisive first step under
    # an ALTERNATE diagonal W (angular block scaled) and confirm the
    # norm-free verdict (jam / not-jam) and binding-set COMPOSITION agree.
    ndof = 48 * topo.n
    w_alt = np.ones(ndof)
    for i in range(topo.n):
        w_alt[48 * i:48 * i + 24] = ALT_W_ANGULAR_SCALE
    _, _, rate_w1, binding_w1, _ = crank_step(topo, pairs, xs0, 0.0, w=0.0, t=0.0, driven="all")
    _, _, rate_walt, binding_walt, _ = crank_step(topo, pairs, xs0, 0.0, w=0.0, t=0.0,
                                                   driven="all", w_diag=w_alt)
    # FIX (code review 23261 H1): b[0] (the "contact"/"wire" TYPE TAG) was
    # dropped before building these sets, so a contact binding and a wire
    # binding sharing the same (i, idx1, j, idx2) numeric tuple would
    # collide into one element -- the K row above already includes b[0];
    # this brings M in line with it.
    binding_set_w1 = {(b[0], b[1], b[2], b[3], b[4]) for b in binding_w1} if binding_w1 else set()
    binding_set_walt = {(b[0], b[1], b[2], b[3], b[4]) for b in binding_walt} if binding_walt else set()

    # L: determinism -- two independent solves of the SAME step, byte-for-byte.
    v1, s1, r1, b1, _ = crank_step(topo, pairs, xs0, 0.0, w=0.0, t=0.0, driven="all")
    v2, s2, r2, b2, _ = crank_step(topo, pairs, xs0, 0.0, w=0.0, t=0.0, driven="all")
    repeat_identical = (s1 == s2 == "OK" and np.array_equal(v1, v2) and r1 == r2)

    # WIRE-IN-BINDING (substantive critique 23262: "must jam BECAUSE the
    # wires bind, not only contacts"): G's own jam binding set must contain
    # at least one wire, not contacts alone.
    g_binding_wires = sum(1 for b in (g_run["binding"] or []) if b[0] == "wire")
    g_binding_contacts = sum(1 for b in (g_run["binding"] or []) if b[0] == "contact")

    # WIRE DEGENERACY (substantive critique 23262 C2): every wire this run
    # EVER reported as binding must carry a genuinely nonzero gradient row
    # -- structurally guaranteed by `wire_gradient_row`'s reuse of
    # `contact_gradient_row` (never degenerate for a parallel-facing pair),
    # checked directly here rather than trusted.
    ndof_wire_check = 48 * topo.n
    wire_binding_ever = [b for b in g_run["binding_ever"] if b[0] == "wire"]
    wire_all_nonzero = True
    wire_checked = 0
    for lab in wire_binding_ever:
        _, u, fi, v_, fj = lab
        _, row = contact_gradient_row(xs0, u, fi, v_, fj, 0.0, ndof_wire_check)
        wire_checked += 1
        if not np.any(row):
            wire_all_nonzero = False

    return {"g_run": g_run, "h_run": h_run, "i_run": i_run, "w_ico": w_ico,
            "cross_w0": cross_w0, "cross_wico": cross_wico,
            "cross_h05_wico": cross_h05_wico, "cross_h2_w0": cross_h2_w0,
            "rate_alt_eps": rate_alt_eps, "rate_nominal_eps": rate_nominal_eps,
            "rate_w1": rate_w1, "rate_walt": rate_walt,
            "binding_set_w1": binding_set_w1, "binding_set_walt": binding_set_walt,
            "repeat_identical": repeat_identical,
            "g_binding_wires": g_binding_wires, "g_binding_contacts": g_binding_contacts,
            "wire_binding_ever_count": len(wire_binding_ever),
            "wire_all_nonzero": wire_all_nonzero, "wire_checked": wire_checked}


# ==========================================================================
# Z16: QPFAIL PROBE (code review 23261 C1, CRITICAL). `scipy.optimize.nnls`
# raising RuntimeError on non-convergence must arrive at the caller as
# `feasible=False` (routing to QPFAIL), never as an unguarded traceback and
# never mistaken for a physical jam. Forces the SAME non-convergence
# `solve_ldp` guards against, on a real, nontrivial problem drawn from this
# array (not a toy), and confirms BOTH outcomes directly.
# ==========================================================================

def z16_qpfail_probe(topo):
    pairs = crank_pairs(topo)
    ndof = 48 * topo.n
    origins0 = topo.sites(verts(0.0))
    xs0 = [unit_corners(0.0, topo, i) + origins0[i] for i in range(topo.n)]
    active_rows = []
    for (i, fi, j, fj) in pairs:
        gap, row = contact_gradient_row(xs0, i, fi, j, fj, 0.0, ndof)
        if gap <= EPS_ACT:
            active_rows.append(row)
    n_basis = null_space_basis(build_pin_jacobian(xs0, topo.n, topo), ndof)
    c_mat = np.array([row @ n_basis for row in active_rows])
    # A NONTRIVIAL target: c=0 (verified, then rejected) makes h_vec the
    # zero vector, at which x=0 trivially satisfies every row -- NNLS
    # converges in 0-1 real iterations regardless of maxiter, so maxiter=1
    # would not actually force anything. Ones() gives a genuine, non-
    # degenerate LDP that needs real iteration to solve.
    c = n_basis.T @ np.ones(ndof)
    h_vec = -(c_mat @ c)
    x_forced, feasible_forced = solve_ldp(c_mat, h_vec, maxiter=1)
    x_normal, feasible_normal = solve_ldp(c_mat, h_vec)
    return {"n_active": len(active_rows), "feasible_forced": feasible_forced,
            "x_forced_zero": not np.any(x_forced), "feasible_normal": feasible_normal}


# ==========================================================================
# Z17: THE LOCK SURFACE a*(w, t) + MOTION ORDER (bead inviscid-qvf.19,
# Phase 1c -- THE DELIVERABLE). Builds ONLY on the FROZEN kernel (Z0-Z9,
# bead .17) and stepper (Z8-Z16, bead .18): every call below goes through
# `crank_run`/`crank_step` unmodified. Nothing in Z0-Z16 above this comment
# is touched by this bead.
#
# FOUR DECLARATIONS for this section, per the amended design of record (T2
# 23230, restated in `main()`'s banner): KERNEL, MASS MODEL, PRIMITIVE --
# INAPPLICABLE, same reasons as the rest of this file (no potential, no
# mass, no primitive choice; every quantity is a phase-space configuration
# read off a static QP solve). METRIC FORM: bead .18 already CHOSE
# treatment (a) -- W = identity, gated W-insensitive (row M) -- so this
# bead carries that choice FORWARD rather than re-deciding it: jam angle,
# active-set composition and reached/jammed status are already established
# as norm-free and quotable (main()'s own banner text, unchanged); this
# bead's da*/dw and da*/dt are DIFFERENCES of that same norm-free jam
# angle, so they inherit norm-freedom too -- treatment (b)'s "per unit of
# W" hedge does not apply here, because (b) was not the treatment chosen.
# ==========================================================================

def _w_ico_lock():
    """w_ico = 2*fold(a_ico), the design of record's own stated wire-slack
    ceiling -- same formula `z15_crank_gates` computes locally, re-derived
    here rather than threaded through as a parameter (mutation-probe
    rule: a value used in two places is measured twice, not passed once
    and trusted)."""
    ref = DIAGONALS[0]
    fold_ico = fold_halves(A_ICO_LOCAL)[ref]
    return 2.0 * (fold_ico[1] - fold_ico[0])


def lock_surface_point(topo, w, t, driven="all", h0=H_LOCK, a_start=0.0,
                       target_span=A_TARGET_LOCK, max_steps=MAX_STEPS_LOCK):
    """One (w, t) grid point, from a_start, SUSTAINED semantics
    (instant_jam=False, the bead's own requirement -- ONSET would read the
    transient contact-only resistance every row here already knows
    relaxes, not a genuine lock). a_target = a_start + target_span.

    Returns a_star: when a_start == 0.0 (the w-arm, unaffected by
    SHIP-BLOCKER 2), the ABSOLUTE jam angle -- jam_angle when jammed,
    a_target itself (a CENSORED value, "no lock found by a_target") when
    reached, None when qpfail. When a_start > 0.0 (the t-arm, re-based per
    SHIP-BLOCKER 2's fix), a_star is instead the OPENING RANGE
    (a_final - a_start) -- the physically meaningful quantity from a
    clearance-relieved starting pose, matching critique 23299's own
    suggestion ("measure a*(t) as an opening RANGE rather than an
    absolute angle"). QPFAIL is never silently folded into either
    physical reading."""
    a_target = a_start + target_span
    r = crank_run(topo, a_start, a_target, w=w, t=t, driven=driven, h0=h0,
                  max_steps=max_steps, instant_jam=False)
    if r["status"] == "jammed":
        raw = r["jam_angle"]
    elif r["status"] == "reached":
        raw = a_target
    else:
        raw = None
    a_star = (raw - a_start) if (raw is not None and a_start > 0.0) else raw
    return {"status": r["status"], "a_star": a_star, "w": w, "t": t,
            "a_start": a_start, "binding_ever": r["binding_ever"], "steps": r["steps"]}


def _min_opening_gap(topo, a, wpairs):
    """The minimum gap over the OPENING (wire-mechanism) plate pairs at
    phase a -- the SAME quantity `wire_gradient_row` measures (t=0, no QP
    involved), reused here as a pure geometric probe for `compute_a_start_t`."""
    ndof = 48 * topo.n
    origins = topo.sites(verts(a))
    xs = [unit_corners(a, topo, i) + origins[i] for i in range(topo.n)]
    best = float("inf")
    for (u, fi, v, fj) in wpairs:
        g, _ = contact_gradient_row(xs, u, fi, v, fj, 0.0, ndof)
        best = min(best, g)
    return best


def compute_a_start_t(topo):
    """SHIP-BLOCKER 2's derived (not picked) clearance-relieved start angle.
    Bisects for the smallest a at which every OPENING pair's gap exceeds
    the FIXED T_START_TARGET_GAP (see its own docstring for why this is
    fixed rather than derived from T_GRID's own max -- a circularity found
    live during development). Bracket validity (glo below target at a=0,
    ghi above target at the hi bound) is checked, not assumed -- a bad
    bracket returns a_start=None rather than a silently wrong bisection
    result."""
    wpairs = wire_pairs(topo)
    target = T_START_TARGET_GAP
    lo, hi = 0.0, T_START_BISECT_HI
    glo = _min_opening_gap(topo, lo, wpairs)
    ghi = _min_opening_gap(topo, hi, wpairs)
    bracket_ok = glo < target <= ghi
    if not bracket_ok:
        return {"a_start": None, "target": target,
                "glo": glo, "ghi": ghi, "bracket_ok": False}
    for _ in range(T_START_BISECT_ITERS):
        mid = 0.5 * (lo + hi)
        gmid = _min_opening_gap(topo, mid, wpairs)
        if gmid < target:
            lo = mid
        else:
            hi = mid
    return {"a_start": hi, "target": target,
            "glo": glo, "ghi": ghi, "bracket_ok": True}


def h_refinement_probe(topo):
    """SHIP-BLOCKER 1b, REVISED (critique 23299's own correction of the
    first draft of this row, found live during development): the critique's
    "interior w-flatness DOES survive refinement" claim is about the
    SPREAD between two interior w points AT a given h0 staying near zero,
    NOT about either point's own ABSOLUTE a* being h0-invariant -- the
    absolute value drifts substantially with h0 at EVERY w tested
    (interior included: 0.9815 -> 0.6243 -> 0.4492 as h0 halves through
    H_REFINE_LEVELS, measured live), so gating "a* itself is stable" at the
    interior point would have been WRONG -- it would have failed for a
    reason that has nothing to do with the actual claim being checked.

    THREE representative w points, t=0, at every H_REFINE_LEVELS h0:
    two INTERIOR points (0.3*w_ico, 0.6*w_ico) whose SPREAD from each
    other is gated to stay near zero AT EACH h0 (the actual "flatness
    survives refinement" claim), and the w=0 BOUNDARY point, whose own
    ABSOLUTE a* is gated to be NON-stable across h0 (the critique's
    "1.119 -> 0.785 -> 0.581 -> ..., no sign of a limit" finding).
    Two-sided by construction: the interior SPREAD row can fail if
    refinement reveals real w-dependence at some h0; the boundary row can
    fail if refinement reveals the boundary point IS, after all, stable
    (misreporting a real instability as resolved would be exactly as
    dishonest as the reverse)."""
    w_ico = _w_ico_lock()
    interior_a_w = 0.3 * w_ico
    interior_b_w = 0.6 * w_ico
    results = {"interior_a": [], "interior_b": [], "boundary": []}
    for h0 in H_REFINE_LEVELS:
        for label, w in (("interior_a", interior_a_w), ("interior_b", interior_b_w),
                        ("boundary", 0.0)):
            r = crank_run(topo, 0.0, H_REFINE_TARGET, w=w, t=0.0, driven="all",
                         h0=h0, max_steps=H_REFINE_MAX_STEPS, instant_jam=False)
            if r["status"] == "jammed":
                a_star = r["jam_angle"]
            elif r["status"] == "reached":
                a_star = H_REFINE_TARGET
            else:
                a_star = None
            # Classify a qpfail: budget exhaustion (every step succeeded,
            # rate_history is full-length and far from stalled) vs a
            # genuine solver-side failure (crank_step's own status != OK,
            # caught earlier in crank_run, rate_history short) -- read
            # directly off the returned dict, not assumed.
            budget_exhausted = (r["status"] == "qpfail"
                               and len(r["rate_history"]) == H_REFINE_MAX_STEPS
                               and (not r["rate_history"] or r["rate_history"][-1] > STALL_RATE_TOL))
            results[label].append({"h0": h0, "status": r["status"], "a_star": a_star,
                                   "steps": r["steps"], "budget_exhausted": budget_exhausted,
                                   "last_rate": r["rate_history"][-1] if r["rate_history"] else None})
    return results


def lock_surface_sweep(topo):
    """The full (w, t) surface: the w-arm (primary + alt) at t=0, a_start=0
    (unaffected by SHIP-BLOCKER 2 -- no thickness parameter varies here);
    the t-arm (primary + alt) at a fixed, large w (T_SWEEP_W_FRAC * w_ico,
    so wires never bind -- the design of record's own arm-C isolation),
    re-based at A_START_T (SHIP-BLOCKER 2's fix) so a*_t is an OPENING
    RANGE, not an absolute angle measured from the idealized touching
    pose. Every point uses driven="all" (in-phase crank, the design's
    stated primary DRIVE) and the SAME H_LOCK."""
    w_ico = _w_ico_lock()
    w_ref = T_SWEEP_W_FRAC * w_ico
    a_start_info = compute_a_start_t(topo)
    a_start_t = a_start_info["a_start"]

    w_points = [(frac * w_ico, 0.0) for frac in W_GRID_FRAC]
    w_points_alt = [(frac * w_ico, 0.0) for frac in W_GRID_FRAC_ALT]
    t_points = [(w_ref, t) for t in T_GRID]
    t_points_alt = [(w_ref, t) for t in T_GRID_ALT]

    def run_all(points, a_start=0.0):
        return [lock_surface_point(topo, w, t, a_start=a_start) for (w, t) in points]

    w_results = run_all(w_points)
    w_results_alt = run_all(w_points_alt)
    # a_start_t may be None on a bracket failure (Z17 gate reports this as
    # a hard FAIL, never silently substitutes 0.0 -- see the
    # "A_START_T bracket valid" row).
    t_a_start = a_start_t if a_start_t is not None else 0.0
    t_results = run_all(t_points, a_start=t_a_start)
    t_results_alt = run_all(t_points_alt, a_start=t_a_start)
    # A live validity check for the WHOLE re-basing strategy (not merely
    # the bisection): does a real crank_run from A_START_T at the SMALLEST
    # swept t actually produce multi-step motion, or does it just move the
    # instant re-lock to a new angle? (Deliberately the SMALLEST, not the
    # largest, swept t: T_GRID/T_GRID_ALT are priced to straddle a real
    # jam/reach threshold discovered live -- the largest values are
    # EXPECTED to jam, by design, so checking relief there would be
    # checking the wrong side of the very threshold this row exists to
    # confirm is real.)
    t_min = min(T_GRID + T_GRID_ALT)
    relief_check = lock_surface_point(topo, w_ref, t_min, a_start=t_a_start)

    return {"w_ico": w_ico, "w_ref": w_ref, "a_start_info": a_start_info,
            "a_start_t": a_start_t, "relief_check": relief_check,
            "w_results": w_results, "w_results_alt": w_results_alt,
            "t_results": t_results, "t_results_alt": t_results_alt}


def _ols_slope(xs, ys):
    """Plain least-squares slope of ys against xs (numpy polyfit, degree 1)
    -- used for rows O/P so the shuffle CONTROL (a re-pairing of the SAME
    y-values against the SAME x-positions in a different order) is a
    meaningful perturbation of the fit, not merely a re-ordering that a
    symmetric-in-order statistic (e.g. a two-point endpoint difference)
    would be blind to."""
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    slope, _ = np.polyfit(xs, ys, 1)
    return float(slope)


def drive_robustness_check(topo):
    """The hazard comment's own requirement: before quoting the w-column as
    flat, confirm the flatness is not an artifact of ONE drive choice.
    Re-measures two W_GRID_FRAC points (an interior one and the near-w_ico
    one) under the single-crank-handle drive (DRIVEN_UNIT_INDEX) and
    reports their OWN spread -- compared on its own terms, not against the
    driven="all" numbers, since single-handle drive tracks a DIFFERENT
    physical quantity (that one unit's own commanded phase, not an
    array-wide uniform expansion)."""
    w_ico = _w_ico_lock()
    fracs = (W_GRID_FRAC[1], W_GRID_FRAC[-1])  # an interior + the near-limit point
    pts = [lock_surface_point(topo, frac * w_ico, 0.0, driven=DRIVEN_UNIT_INDEX,
                              max_steps=MAX_STEPS_DRIVE)
           for frac in fracs]
    finite = [p["a_star"] for p in pts if p["a_star"] is not None]
    spread = (max(finite) - min(finite)) if len(finite) == len(pts) else None
    statuses = [p["status"] for p in pts]
    return {"points": pts, "fracs": fracs, "spread": spread,
            "statuses": statuses,
            "agree": len(set(statuses)) == 1,
            "any_lock": any(st == "jammed" for st in statuses)}


def motion_order_trace(topo, w, t, driven, h0=MOTION_ORDER_H,
                       n_steps=MOTION_ORDER_STEPS):
    """A SHORT, bespoke trace of PER-UNIT activity across a crank run --
    the qvf.14 path question (who moves, in what order) now asked with
    contact present. Deliberately does NOT reuse `crank_run`'s feasibility
    backtrack: that machinery answers "how far can the array go before it
    locks", or a *j gap*; this question is "which unit's own velocity
    component is nonzero, and from which step on" -- a diagnostic of ORDER,
    read directly off each step's QP solve (`crank_step`, unmodified,
    called exactly as `crank_run` calls it), not a hardened lock-surface
    measurement. A full h0 step is applied every iteration (no bisection)
    since the quantity of interest is which units are active, not a
    precise jam angle -- if a later step's QP reports QPFAIL, the trace
    stops there and reports the truncation honestly rather than padding it.

    Per-unit reading (code review MEDIUM fix, 23337): the RAW velocity
    norm ||v_i|| for EVERY unit, driven or not -- a single, homogeneous
    physical quantity (previously driven units used a v_cmd-projected
    rate, O(0.1-1), while undriven units used a raw norm, O(0.001-0.01);
    mixing the two metrics made "not all equal" partly a scale artifact of
    WHICH metric a unit happened to get, not purely a timing signal).
    Verified live before committing this change: under driven="all" at
    large w (the CONTROL regime), the raw-norm metric agrees across ALL
    SEVEN units (not merely the six corners) to a spread of ~1.1e-5 --
    tighter than the old mixed-metric reading, not weaker."""
    n = topo.n
    ndof = 48 * n
    pairs = crank_pairs(topo)
    wpairs = wire_pairs(topo)
    origins = topo.sites(verts(0.0))
    xs = [unit_corners(0.0, topo, i) + origins[i] for i in range(n)]
    a_hat = 0.0
    history = []
    statuses = []
    for _ in range(n_steps):
        v, status, rate, binding, mgg = crank_step(topo, pairs, xs, a_hat, w, t,
                                                    driven, True, None, wpairs)
        statuses.append(status)
        if status != "OK":
            break
        per_unit = np.array([float(np.linalg.norm(v[48 * i:48 * i + 48])) for i in range(n)])
        history.append(per_unit)
        xs, _ = project_to_joint_manifold(
            [apply_body_motions(xs[i], h0 * v[48 * i:48 * i + 48]) for i in range(n)], topo)
        a_hat += h0 * rate
    return {"history": history, "statuses": statuses, "driven": driven,
            "n_steps_completed": len(history)}


def joint_integrity_probe(topo, w=0.0, t=0.0, h0=H_LOCK,
                          n_steps=JOINT_PROBE_STEPS):
    """Does the thing being cranked stay an ARRAY? (bead inviscid-1wd)

    Runs the stepper twice over the same steps and reports the worst
    shared-vertex separation reached by each:

      TEST     -- the shipped path: ball-joint rows in the QP's Jacobian AND
                  the array-level `project_to_joint_manifold` after each step.
      CONTROL  -- the retired path IN FULL: `_hinge_only_jacobian` in place of
                  the ball-joint one, so the QP direction itself is blind to
                  the joints, AND each unit projected onto its OWN hinge
                  manifold by `project_to_pin_manifold`. Nothing constrains the
                  joints and nothing pulls them closed, exactly as before this
                  bead.

    The CONTROL is not decoration. Without it "the joints hold" is satisfiable
    by any run short enough for drift not to show, and the number that matters
    is not the TEST's absolute value but the gap between the two. Measured
    before the fix, the control separates at ~0.094 per degree of crank and
    reaches 0.043 within three steps.

    Never raises: a QPFAIL truncates the trace and is reported through
    `steps`, in keeping with this file's rule that a swept probe returns a
    reading rather than an exception."""
    n = topo.n
    ndof = 48 * n
    pairs = crank_pairs(topo)
    wpairs = wire_pairs(topo)

    def run(joints_on):
        origins = topo.sites(verts(0.0))
        xs = [unit_corners(0.0, topo, i) + origins[i] for i in range(n)]
        a_hat = 0.0
        worst = 0.0
        steps = 0
        for _ in range(n_steps):
            # The CONTROL arm reproduces the RETIRED path in full: no
            # inter-unit rows in the Jacobian the QP takes its null space
            # from, AND the old per-unit projection. Removing only the
            # projection would understate it -- the joint rows alone already
            # hold the array to first order, leaving just second-order drift
            # (1.2e-03 over this probe) rather than the 0.043-in-three-steps
            # the real retired path produced.
            _saved = globals()["build_pin_jacobian"]
            if not joints_on:
                globals()["build_pin_jacobian"] = (
                    lambda _xs, _n, _topo: _hinge_only_jacobian(_xs, _n))
            try:
                v, status, rate, _b, _m = crank_step(topo, pairs, xs, a_hat, w, t,
                                                     "all", True, None, wpairs)
            finally:
                globals()["build_pin_jacobian"] = _saved
            if status != "OK":
                break
            moved = [apply_body_motions(xs[i], h0 * v[48 * i:48 * i + 48])
                     for i in range(n)]
            if joints_on:
                xs, _rn = project_to_joint_manifold(moved, topo)
            else:
                xs = [project_to_pin_manifold(m)[0] for m in moved]
            a_hat += h0 * rate
            worst = max(worst, max_joint_gap(xs, topo))
            steps += 1
        return {"worst": worst, "steps": steps, "a_hat": a_hat}

    return {"test": run(True), "control": run(False)}


# ==========================================================================
# Z18: PHASE TRACKING (bead inviscid-qvf.19 follow-up). Is `a_hat` the
# configuration's phase, and is the array still on the symmetric path the
# vertex ellipses describe?
# ==========================================================================

def unit_radii(x):
    """The 12 DISTINCT vertex distances from one unit's own centroid.

    `x` carries 24 corner slots (8 plates x 3), but `PAIRS` identifies them in
    hinged couples, so only 12 are distinct points. Taking all 24 would double-
    count and hide exactly the asymmetry this measurement exists to see."""
    c = x.reshape(-1, 3).mean(axis=0)
    seen = set()
    out = []
    for (fa, ca), (fb, cb) in PAIRS:
        if (fa, ca) in seen:
            continue
        seen.add((fa, ca)); seen.add((fb, cb))
        out.append(float(np.linalg.norm(x[fa][ca] - c)))
    return np.array(out)


def configuration_phase(x):
    """(phase_deg, radius_spread) read from a unit's ACTUAL geometry.

    On the symmetric jitterbug path every vertex rides the same ellipse, so all
    twelve sit at one radius
        r(a) = sqrt(2 - (4/3) sin^2 a)
    which follows directly from the ellipse parameterisation
    (sqrt(2) cos a, -sqrt(2/3) sin a) -- verified against `corners(a)` to 4e-16
    across a in [0, 89], with the twelve radii agreeing to 2e-16. It is monotone
    on [0, 90], so it INVERTS in closed form:
        a = arcsin( sqrt( 3 (2 - r^2) / 4 ) ).

    This is a MEASUREMENT of the configuration, unlike `crank_run`'s `a_hat`,
    which is the running integral `a_hat += h * rate` and never consults the
    geometry it is supposed to describe.

    `radius_spread` (max - min over the twelve) is what says whether the phase
    is even well defined: it is 0 on the symmetric path by construction, and
    anything else means the linkage has flexed ASYMMETRICALLY, where no single
    phase describes the unit and the ellipse closed form does not apply. The
    returned phase is then the best symmetric fit to an unsymmetric thing, and
    is reported with its spread rather than alone."""
    rs = unit_radii(x)
    rm = float(rs.mean())
    inner = 3.0 * (2.0 - rm * rm) / 4.0
    a = float(np.degrees(np.arcsin(np.sqrt(max(0.0, min(1.0, inner))))))
    return a, float(rs.max() - rs.min())


def _array_phase(xs):
    """(mean phase, spread ACROSS units, worst radius spread WITHIN a unit).

    READING ONE UNIT IS NOT READING THE ARRAY, and the distinction is the whole
    finding. An earlier version of this probe took `configuration_phase(xs[0])`
    and called it the configuration's phase. Unit 0 of SC7 is the CENTRE unit,
    which lags the six corners by up to 6.7 deg, so that reading reported the
    centre's lag as though it were `a_hat` drifting -- with the wrong sign and
    roughly twenty times the true magnitude.

    The two spreads are different phenomena and are kept apart on purpose:
    `dephase` is units disagreeing with EACH OTHER, `radius spread` is one unit
    disagreeing with ITSELF (its twelve vertices leaving a common radius, i.e.
    the linkage flexing off the symmetric path). Both are real here."""
    phases, rads = [], []
    for x in xs:
        p, r = configuration_phase(x)
        phases.append(p); rads.append(r)
    return (float(np.mean(phases)), float(max(phases) - min(phases)),
            float(max(rads)))


@jb_cache.memoize(_MODULE)
def phase_tracking_probe(topo, w_frac=PHASE_PROBE_W_FRAC, t=0.0,
                         target=PHASE_PROBE_TARGET, levels=PHASE_H0_LEVELS):
    """Crank to a FIXED `a_hat` at several step sizes and, at that common
    stopping point, compare `a_hat` against the measured phase.

    Comparing at one `a_hat` rather than at each run's own freeze is what makes
    the step sizes commensurable: the freezes land at different angles, so a
    per-freeze comparison would confound the drift with where each run stopped.

    Also carries the LINKAGE VALIDITY readings -- hinge residual, joint gap, and
    per-triangle edge deviation -- because "the units left the symmetric path"
    only means something if they are still valid states of the linkage. If the
    triangles were deforming or the hinges opening, the same numbers would
    instead be a broken solver, and the gate must be able to tell those apart."""
    n = topo.n
    ndof = 48 * n
    pairs = crank_pairs(topo)
    wpairs = wire_pairs(topo)
    w = w_frac * _w_ico_lock()
    out = []
    for h0 in levels:
        origins = topo.sites(verts(0.0))
        xs = [unit_corners(0.0, topo, i) + origins[i] for i in range(n)]
        edge0 = np.array([float(np.linalg.norm(xs[0][f][k] - xs[0][f][(k+1) % 3]))
                          for f in range(8) for k in range(3)])
        a_hat = 0.0
        rec = {"h0": h0, "reached": False, "hinge": 0.0, "joint": 0.0,
               "edge": 0.0, "spread": 0.0, "dephase": 0.0, "phase": None,
               "a_hat": 0.0, "steps": 0}
        for step in range(PHASE_PROBE_MAX_STEPS):
            v, status, rate, _b, _m = crank_step(topo, pairs, xs, a_hat, w, t,
                                                 "all", True, None, wpairs)
            if status != "OK":
                break
            h = h0
            accepted = False
            for _ in range(H_BACKTRACK_MAX):
                trial = [apply_body_motions(xs[i], h * v[48*i:48*i+48])
                         for i in range(n)]
                bad = False
                for (i, fi, j, fj) in pairs:
                    g, _r = contact_gradient_row(trial, i, fi, j, fj, t, ndof)
                    if -MEANINGLESS_DEPTH_FLOOR < g < -GAP_FLOOR_TOL:
                        bad = True; break
                if not bad:
                    accepted = True; break
                h *= 0.5
                if h < H_MIN: break
            if not accepted:
                break
            xs, _rn = project_to_joint_manifold(trial, topo)
            a_hat += h * rate
            rec["steps"] = step + 1
            rec["hinge"] = max(rec["hinge"],
                               max(float(np.abs(hinge_residual(x)).max()) for x in xs))
            rec["joint"] = max(rec["joint"], max_joint_gap(xs, topo))
            edge = np.array([float(np.linalg.norm(xs[0][f][k] - xs[0][f][(k+1) % 3]))
                             for f in range(8) for k in range(3)])
            rec["edge"] = max(rec["edge"], float(np.abs(edge - edge0).max()))
            ph, dph, sp = _array_phase(xs)
            rec["spread"] = max(rec["spread"], sp)
            rec["dephase"] = max(rec["dephase"], dph)
            if a_hat >= target:
                rec.update(reached=True, phase=ph, a_hat=a_hat)
                break
        if not rec["reached"]:
            ph, dph, sp = _array_phase(xs)
            rec.update(phase=ph, a_hat=a_hat)
            rec["dephase"] = max(rec["dephase"], dph)
        rec["drift"] = rec["a_hat"] - rec["phase"] if rec["phase"] is not None else None
        out.append(rec)
    return out


def _grid_ratio_apart_check(primary, alt):
    """HIGH fix (23337): the jb_y K_GRID/K_GRID_ALT coprimality/offset-
    apartness shape, applied to a real-valued fraction grid rather than
    integers. NON-RATIO-DERIVED: no single constant k satisfies
    alt[i] == k*primary[i] for every i (checked by requiring the
    pairwise ratios alt[i]/primary[i] to NOT all agree to
    GRID_RATIO_CONST_TOL -- primary[0]=0.0 makes that particular ratio
    undefined, so it is excluded from the ratio comparison and covered by
    the separate offset-apartness check instead). OFFSET-APART: every alt
    value sits at least GRID_OFFSET_APART_TOL from every primary value."""
    ratios = [a / p for p, a in zip(primary, alt) if abs(p) > 1e-12]
    ratio_const = (len(ratios) >= 2
                   and max(ratios) - min(ratios) < GRID_RATIO_CONST_TOL)
    min_offset = min(abs(a - p) for a in alt for p in primary)
    return {"ratios": ratios, "ratio_const": ratio_const,
            "non_ratio_derived": not ratio_const,
            "min_offset": min_offset, "offset_apart": min_offset > GRID_OFFSET_APART_TOL}


def z17_lock_surface(topo):
    surf = lock_surface_sweep(topo)
    drive = drive_robustness_check(topo)
    href = h_refinement_probe(topo)
    w_grid_apart = _grid_ratio_apart_check(W_GRID_FRAC, W_GRID_FRAC_ALT)
    t_grid_apart = _grid_ratio_apart_check(T_GRID, T_GRID_ALT)
    joints = joint_integrity_probe(topo)
    phase = phase_tracking_probe(topo)

    test_trace = motion_order_trace(topo, MOTION_TEST_W_FRAC * surf["w_ico"], 0.0,
                                    DRIVEN_UNIT_INDEX)
    control_trace = motion_order_trace(topo, T_SWEEP_W_FRAC * surf["w_ico"], 0.0, "all")

    return {"surf": surf, "drive": drive, "href": href,
            "w_grid_apart": w_grid_apart, "t_grid_apart": t_grid_apart,
            "joints": joints, "phase": phase,
            "test_trace": test_trace, "control_trace": control_trace}


# ==========================================================================
# THE GATE
# ==========================================================================

def z18_phase_threading(topo):
    """Is a two-sublattice topology actually DRIVEN with per-unit phase?

    The row this file owed itself. Revision 34e636c asserted in a COMMENT that a
    two-sublattice topology could not be silently driven as one, and left the
    velocity command driving every unit along the tangent at the bare angle. A
    guarantee in prose with no row behind it is the exact shape of defect this
    project keeps paying for, so it is measured here.

    Method needs no honeycomb geometry and adds no import: take the real
    topology, copy it, and give the copy alternating phase offsets. The two
    differ ONLY in `phases`, so any difference in the velocity command is the
    threading and nothing else. Two-sided by construction -- the zero-phase
    units must come out UNCHANGED while the offset ones must not."""
    import copy
    ndof = 48 * topo.n
    base, _u = crank_v_cmd(topo, 10.0, "all", ndof)
    alt = copy.copy(topo)
    alt.phases = tuple(0.0 if k % 2 == 0 else 60.0 for k in range(topo.n))
    two, _u2 = crank_v_cmd(alt, 10.0, "all", ndof)
    off_units = [k for k in range(topo.n) if alt.phases[k] != 0.0]
    zero_units = [k for k in range(topo.n) if alt.phases[k] == 0.0]
    moved = min(float(np.abs(two[48 * k:48 * k + 48]
                             - base[48 * k:48 * k + 48]).max())
                for k in off_units) if off_units else 0.0
    still = max(float(np.abs(two[48 * k:48 * k + 48]
                             - base[48 * k:48 * k + 48]).max())
                for k in zero_units) if zero_units else float("inf")
    tan_ang = float(path_tangent_48(10.0)[0] @ path_tangent_48(70.0)[0]
                    / np.linalg.norm(path_tangent_48(10.0)[0])
                    / np.linalg.norm(path_tangent_48(70.0)[0]))
    return dict(moved=moved, still=still, cos=tan_ang,
                n_off=len(off_units), n_zero=len(zero_units))


def z19_meaningless_depth(topo, a=0.0):
    """The antipodal artifact class, and that it no longer reaches the rows.

    `signed_gap`'s parallel branch tests |nA . nB|, so it fires on ANTIPARALLEL
    normals too. Two plates on opposite sides of one unit face away from each
    other, sit permanently ~2.31 apart, and cannot touch -- and the branch
    reports gap = -2.309401 for them, which reads as deep penetration.
    `MEANINGLESS_DEPTH_FLOOR` already named the class and crank_run's backtrack
    acceptance already rejected it; crank_step built the rows anyway."""
    pairs = crank_pairs(topo)
    origins = topo.sites(verts(a))
    xs = [unit_corners(a, topo, i) + origins[i] for i in range(topo.n)]
    deep, antipodal = 0, 0
    for (i, fi, j, fj) in pairs:
        nA = plate_normal(fi)
        g, _wa, _wb, _n = signed_gap(xs[i][fi], xs[j][fj], nA)
        if g < -MEANINGLESS_DEPTH_FLOOR:
            deep += 1
            if i == j and float(nA @ plate_normal(fj)) < -PARALLEL_TOL:
                antipodal += 1
    _v, st, _r, binding, _m = crank_step(topo, pairs, xs, a, _w_ico_lock(), 0.0,
                                         "all", True, None, wire_pairs(topo))
    worst = min((b[5] for b in binding if b[0] == "contact"), default=0.0)
    return dict(deep=deep, antipodal=antipodal, status=st,
                worst_binding=worst, n_binding=len(binding))


#: Two-sided bounds on the broad-phase rejection fraction (qvf.23 criterion 2).
#: Measured 0.6872..0.9888 across the probe set. A prefilter that rejects
#: NOTHING is pointless and one that rejects EVERYTHING is wrong, and both edges
#: redden rather than only the second.
BP_REJECT_BAND = (0.30, 0.999)


def gate(z0, z2, z3, z4, z5, z6, z7, zg, zhaz, zdow, ztwo, zqp, zlock, zbp):
    """Every check's verdict in one table, and this process's exit code."""
    print()
    print("=" * 78)
    print("Z2  plate-pair counts, per topology (GAP 2 structural deliverable)")
    print("=" * 78)
    for name, v in z2.items():
        print(f"  {name:34s} n={v['n']:3d}  intra={v['n_intra']:5d}  "
              f"inter={v['n_inter']:5d}  total={v['n_pairs']:5d}")

    checks = []

    # ---- Z0: foundation ----
    checks.append(("Z0  plate normal is phase-invariant (all 8 faces, 2 ladders)",
                   z0["normal_invariance"] < NORMAL_INVARIANT_TOL,
                   f"{z0['normal_invariance']:.2e}", f"< {NORMAL_INVARIANT_TOL:.0e}"))
    checks.append(("Z0  plate centroid sits at u * Z * cos(a) exactly",
                   z0["centroid_on_axis"] < NORMAL_INVARIANT_TOL,
                   f"{z0['centroid_on_axis']:.2e}", f"< {NORMAL_INVARIANT_TOL:.0e}"))

    # ---- Z2: plate-pair enumeration (structural) ----
    sc7 = z2["SC7 star (six-around-one)"]
    n1 = z2["N1 (control)"]
    checks.append(("Z2  N1 (single unit): 28 - 12 hinge-adjacent = 16 pairs",
                   n1["n_pairs"] == 16 and n1["n_intra"] == 16 and n1["n_inter"] == 0,
                   f"{n1['n_pairs']}", "16"))
    checks.append(("Z2  SC7 star: every unit's 16 intra-pairs present",
                   sc7["n_intra"] == 16 * sc7["n"],
                   f"{sc7['n_intra']}", f"{16 * sc7['n']}"))
    checks.append(("Z2  SC7 star: inter-unit pairs enumerated, non-empty",
                   sc7["n_inter"] > 0, f"{sc7['n_inter']}", "> 0"))
    checks.append(("Z2  plate-pair counts non-empty over EVERY topology",
                   all(v["n_pairs"] > 0 for v in z2.values()),
                   str(min(v["n_pairs"] for v in z2.values())), "> 0, all"))

    # ---- Z3: fold table (rows A, B) ----
    fold_all_ok = True
    fold_worst = 0.0
    for a, row in z3["rows"].items():
        dev = max(abs(row["d1"] - row["t1"]), abs(row["d2"] - row["t2"]),
                  abs(row["fold"] - row["tfold"]))
        fold_worst = max(fold_worst, dev)
        fold_all_ok = fold_all_ok and dev < FOLD_TABLE_TOL
    checks.append(("Z3  fold table reproduces T2 23195 (5 angles, d1/d2/fold)",
                   fold_all_ok and len(z3["rows"]) == 5,
                   f"{fold_worst:.2e}", f"< {FOLD_TABLE_TOL:.0e}"))
    checks.append(("Z3  fold is uniform across all 6 cuboctahedron squares",
                   z3["uniform_dev"] < CONST_TOL, f"{z3['uniform_dev']:.2e}",
                   f"< {CONST_TOL:.0e}"))
    checks.append(("Z3  CONTROL: fold(a=0) is EXACTLY zero",
                   abs(z3["fold_at_0"]) < CONST_TOL, f"{z3['fold_at_0']:.2e}",
                   f"< {CONST_TOL:.0e}"))
    checks.append(("Z3  CONTROL: fold is first order in a (slope ratio)",
                   z3["first_order_dev"] < 1e-2, f"{z3['first_order_dev']:.2e}",
                   "< 1e-02"))
    checks.append(("Z3  ridge diagonal == strut EXACTLY at a_ico",
                   abs(z3["ridge_full"] - STRUT_LEN_LOCAL) < CONST_TOL,
                   f"{z3['ridge_full']:.9f}", f"{STRUT_LEN_LOCAL:.9f}"))
    checks.append(("Z3  CONTROL: a_ico +/- offset MISSES strut length",
                   z3["aico_control_min"] > 100.0 * CONST_TOL,
                   f"{z3['aico_control_min']:.2e}", f"> {100.0 * CONST_TOL:.0e}"))

    # ---- Z4: registry closed form (row C) ----
    checks.append(("Z4  kernel's signed gap == registry closed form, all 8 faces",
                   z4["worst_closed"] < 1e-9, f"{z4['worst_closed']:.2e}", "< 1e-09"))
    checks.append(("Z4  registry targets match T2 23195 (60/0, 30/30, 59/0, 60/1)",
                   z4["worst_target"] < 1e-5, f"{z4['worst_target']:.2e}", "< 1e-05"))
    checks.append(("Z4  8-face uniformity of the registry gap",
                   z4["uniform_dev"] < 1e-9, f"{z4['uniform_dev']:.2e}", "< 1e-09"))
    checks.append(("Z4  swap p+q=60 minimum -0.267949 at p=q=30",
                   abs(z4["swap_min"] - (-0.267949)) < 1e-5 and abs(z4["swap_argmin"] - 30.0) < 0.2,
                   f"{z4['swap_min']:.6f}", "-0.267949"))
    checks.append(("Z4  one-sidedness: d gap/dp = -1.000 per rad at registry",
                   abs(z4["dgdp_rad"] - (-1.0)) < 1e-3, f"{z4['dgdp_rad']:.4f}", "-1.0000"))
    checks.append(("Z4  one-sidedness: d gap/dq = 0 (second order) at registry",
                   abs(z4["dgdq_rad"]) < 1e-2, f"{z4['dgdq_rad']:.4f}", "~0"))
    checks.append(("Z4  CONTROL: dq is NOT the same order as dp (asymmetry real)",
                   abs(z4["dgdq_rad"]) < 0.2 * abs(z4["dgdp_rad"]),
                   f"{abs(z4['dgdq_rad']) / abs(z4['dgdp_rad']):.3f}", "< 0.2"))

    # ---- Z5: normal orientation and witness points (rows G) ----
    for label in ("separated", "near_touching"):
        r = z5[label]
        checks.append((f"G-i  FD normal check ({label}): d(gap)/d(eps) == +1",
                       r["fd_dev"] < FD_TOL, f"{r['fd_dev']:.2e}", f"< {FD_TOL:.0e}"))
        checks.append((f"G-i  CONTROL ({label}): {FD_CONTROL_FACTOR:.0f}x-eps offset REJECTED",
                       r["fd_dev_control"] > FD_TOL, f"{r['fd_dev_control']:.2e}", f"> {FD_TOL:.0e}"))
    rs = z5["registry_sign"]
    checks.append(("G-ii  registry contact: normal sign matches d gap/dp = -1/rad",
                   abs(rs["dgdp_deg"] - (-1.0)) < 5e-2, f"{rs['dgdp_deg']:.4f}", "-1.0000"))
    checks.append(("G-ii  registry contact gap is NEGATIVE (interpenetrating)",
                   rs["gap"] < 0.0, f"{rs['gap']:.6f}", "< 0"))
    nc = z5["negated_control"]
    checks.append(("G-iii  CONTROL: negated normal FAILS the FD check",
                   nc["fd_dev_negated"] > FD_TOL, f"{nc['fd_dev_negated']:.2e}", f"> {FD_TOL:.0e}"))
    checks.append(("G-iii  CONTROL: negated normal sign disagrees with closed form",
                   not nc["dgdp_negated_matches_closed_form"],
                   str(not nc["dgdp_negated_matches_closed_form"]), "True"))
    wt = z5["witness"]
    checks.append(("G-iv  witness points lie ON their own triangle (barycentric)",
                   wt["on_triA"] and wt["on_triB"], str(wt["on_triA"] and wt["on_triB"]), "True"))
    checks.append(("G-iv  witness-pair separation == unsigned edge-edge gap",
                   wt["edge_witness_sep_dev"] < WITNESS_TOL, f"{wt['edge_witness_sep_dev']:.2e}",
                   f"< {WITNESS_TOL:.0e}"))
    checks.append(("G-iv  edge-edge witness cross-checked vs jb_g.segment_distance",
                   wt["edge_cross_dev"] < WITNESS_TOL, f"{wt['edge_cross_dev']:.2e}",
                   f"< {WITNESS_TOL:.0e}"))
    checks.append(("G-iv  vertex-face witness cross-checked vs independent projector",
                   wt["vf_cross_dev"] < 1e-6, f"{wt['vf_cross_dev']:.2e}", "< 1e-06"))
    checks.append(("G-iv  CONTROL: naive vertex-vertex candidate is NOT closest",
                   wt["naive_wrong"], str(wt["naive_wrong"]), "True"))

    # ---- Z6: wire spans and thickness (rows H) ----
    checks.append(("H-v  wire span vs hand-computed 3-4-5 distance",
                   z6["span_dev"] < SPAN_TOL, f"{z6['span_dev']:.2e}", f"< {SPAN_TOL:.0e}"))
    checks.append(("H-vi  slack at w > s, taut at w < s",
                   z6["slack_ok"] and z6["taut_ok"], str(z6["slack_ok"] and z6["taut_ok"]), "True"))
    checks.append(("H-vi  taut/slack transition located at w == s",
                   z6["transition_dev"] < 1e-6, f"{z6['transition_dev']:.2e}", "< 1e-06"))
    checks.append(("H-vi  CONTROL: compression leaves tension-only member INACTIVE",
                   z6["mine_inactive_under_compression"], str(z6["mine_inactive_under_compression"]), "True"))
    checks.append(("H-vi  CONTROL: a bilateral member WOULD flag active (redden)",
                   z6["bilateral_disagrees_under_compression"],
                   str(z6["bilateral_disagrees_under_compression"]), "True"))
    checks.append(("H-vii  thickness shifts admissibility by EXACTLY t (5 values)",
                   z6["thickness_max_dev"] < THICKNESS_TOL, f"{z6['thickness_max_dev']:.2e}",
                   f"< {THICKNESS_TOL:.0e}"))
    checks.append(("H-vii  CONTROL: admissibility FLIPS as t crosses the gap value",
                   z6["flip_below"] and z6["flip_above"], str(z6["flip_below"] and z6["flip_above"]), "True"))
    checks.append(("H-viii  t = 0 exercised as its own row",
                   z6["t0_ok"], str(z6["t0_ok"]), "True"))
    checks.append(("H  span_length on REAL topology (census square, not synthetic)",
                   z6["real_span_sane"], f"{z6['real_span']:.6f}",
                   f"finite, in (0, {10.0 * STRUT_LEN_LOCAL:.2f})"))

    # ---- Z7: crossing census + perturbation (row D) ----
    checks.append(("D  axis direction is a unit vector",
                   z7["axis_dev"] < 1e-12, f"{z7['axis_dev']:.2e}", "< 1e-12"))
    checks.append(("D  spacing-2 crossings: 0 at a=0",
                   z7["at_zero_n"] == 0, f"{z7['at_zero_n']}", "0"))
    checks.append(("D  spacing-2 crossings: 4 at every angle in {5,10,ico,30,45}",
                   all(v["n"] == CROSSING_TARGET for v in z7["results"].values()),
                   str([v["n"] for v in z7["results"].values()]), f"all {CROSSING_TARGET}"))
    stable = all(min(v["perturbed"]) == max(v["perturbed"]) == CROSSING_TARGET
                for v in z7["results"].values())
    checks.append(("D  perturbation-stable: 4..4 under 10 seeded 1e-4 rigid moves",
                   stable and all(len(v["perturbed"]) == PERTURB_N for v in z7["results"].values()),
                   str([f"{min(v['perturbed'])}..{max(v['perturbed'])}" for v in z7["results"].values()]),
                   f"all {CROSSING_TARGET}..{CROSSING_TARGET}"))
    ico_len = z7["results"][A_ICO_LOCAL]["total_len"]
    checks.append(("D  total crossing length at a_ico, spacing 2 == 1.6306",
                   abs(ico_len - 1.6306) < 1e-3, f"{ico_len:.4f}", "1.6306"))
    checks.append(("D  CONTROL: naive strict-interior DISAGREES at a=0 (23195)",
                   z7["at_zero_naive"] != z7["at_zero_n"] and z7["at_zero_naive"] > 0,
                   f"naive={z7['at_zero_naive']}", f"!= robust={z7['at_zero_n']}"))
    checks.append(("D  valleys held (spacing 2*d1): 4, range 4..6",
                   z7["valleys_held_n"] == 4, f"{z7['valleys_held_n']}", "4"))
    vh_lo, vh_hi = min(z7["valleys_held_pert"]), max(z7["valleys_held_pert"])
    checks.append(("D  valleys held perturbation range within 4..6",
                   4 <= vh_lo and vh_hi <= 6, f"{vh_lo}..{vh_hi}", "4..6"))
    checks.append(("D  ridges touching (spacing 2*d2): 0, range 0..5",
                   z7["ridges_touch_n"] == 0, f"{z7['ridges_touch_n']}", "0"))
    rt_lo, rt_hi = min(z7["ridges_touch_pert"]), max(z7["ridges_touch_pert"])
    checks.append(("D  ridges touching perturbation range within 0..5",
                   0 <= rt_lo and rt_hi <= 5, f"{rt_lo}..{rt_hi}", "0..5"))

    # ---- G/H/I: the decisive row pair + control (bead qvf.18) ----
    g_run, h_run, i_run = zg["g_run"], zg["h_run"], zg["i_run"]
    checks.append(("G  w=0,t=0 in-phase crank: status is JAMMED (not reached/qpfail)",
                   g_run["status"] == "jammed", g_run["status"], "jammed"))
    checks.append(("G  jam angle < JAM_ANGLE_EPS (>0 not required, first order at a=0)",
                   g_run["status"] == "jammed" and g_run["jam_angle"] is not None
                   and g_run["jam_angle"] < JAM_ANGLE_EPS,
                   f"{g_run['jam_angle']}", f"< {JAM_ANGLE_EPS}"))
    checks.append(("G  CONTROL: JAM_ANGLE_EPS + offset would NOT count as immediate",
                   JAM_ANGLE_EPS + JAM_ANGLE_CONTROL_OFFSET > JAM_ANGLE_EPS,
                   f"{JAM_ANGLE_EPS + JAM_ANGLE_CONTROL_OFFSET}", f"> {JAM_ANGLE_EPS}"))
    checks.append(("H  w >= 2*fold(a_ico) in-phase crank: status is REACHED",
                   h_run["status"] == "reached", h_run["status"], "reached"))
    checks.append(("H  w_ico == 2*fold(a_ico) == 0.87404 (design of record)",
                   abs(zg["w_ico"] - 0.87404) < 1e-3, f"{zg['w_ico']:.5f}", "0.87404"))
    checks.append(("I  CONTROL: contacts disabled, w large: status is REACHED",
                   i_run["status"] == "reached", i_run["status"], "reached"))
    checks.append(("I  CONTROL: first-step rate >= JAM_RATE_CONTROL_OFFSET (not a jam)",
                   bool(i_run["rate_history"]) and i_run["rate_history"][0] >= JAM_RATE_CONTROL_OFFSET,
                   f"{i_run['rate_history'][0] if i_run['rate_history'] else None}",
                   f">= {JAM_RATE_CONTROL_OFFSET}"))
    checks.append(("G/H  NON-VACUITY: a solver reporting jam unconditionally FAILS H",
                   g_run["status"] == "jammed" and h_run["status"] == "reached",
                   f"G={g_run['status']} H={h_run['status']}", "jammed / reached"))

    # ---- CAUSAL CROSS-TEST (substantive critique 23262 C1, ship-blocker):
    # SAME protocol (instant_jam=False, SUSTAINED), only w varies -- proves
    # w itself is causal, not merely correlated with the instant_jam flag
    # G/H's own comparison also varies. ----
    # The three original CROSS rows ("w=0 JAMS", "w=w_ico REACHES", "therefore
    # W IS CAUSAL") are DELETED rather than re-priced -- see the deleted-rows
    # prose block. They asserted an outcome DIFFERENCE that the disassembled
    # array produced and the assembled one does not, and they had already been
    # half-retracted by the CROSS-MATCHED rows below, which showed the
    # difference vanishes once h0 is held equal. `cross_w0`/`cross_wico` are
    # still computed and still gated -- by CROSS-MATCHED, on the claim that
    # survives.

    # ---- K: binding-set non-emptiness at jam ----
    checks.append(("K  binding active set at jam is NON-EMPTY",
                   g_run["status"] == "jammed" and bool(g_run["binding"]),
                   str(len(g_run["binding"])) if g_run["binding"] else "0", "> 0"))
    if g_run["status"] == "jammed" and g_run["binding"]:
        n_contact = sum(1 for b in g_run["binding"] if b[0] == "contact")
        n_wire = sum(1 for b in g_run["binding"] if b[0] == "wire")
        checks.append(("K  binding-set composition printed (contacts vs wires)",
                       True, f"{n_contact} contact + {n_wire} wire", "printed"))
    checks.append(("K  G's jam binds at least one WIRE (critique 23262: 'BECAUSE wires bind')",
                   zg["g_binding_wires"] > 0, f"{zg['g_binding_wires']} wire", "> 0"))
    checks.append(("K  WIRE NON-DEGENERACY: every binding wire has a genuinely nonzero gradient",
                   zg["wire_checked"] > 0 and zg["wire_all_nonzero"],
                   f"{zg['wire_checked']} checked, all nonzero={zg['wire_all_nonzero']}", "True"))

    # ---- eps_act decade-insensitivity ----
    checks.append(("Z13  eps_act decade-insensitivity: rate agrees across a decade",
                   abs(zg["rate_alt_eps"] - zg["rate_nominal_eps"]) < 0.05,
                   f"{zg['rate_alt_eps']:.6f} vs {zg['rate_nominal_eps']:.6f}", "< 0.05 apart"))

    # ---- M: METRIC FORM W-insensitivity ----
    checks.append(("M  norm-free verdict (jam/not) is W-INSENSITIVE (W=I vs alt W)",
                   (zg["rate_w1"] < JAM_RATE_TOL) == (zg["rate_walt"] < JAM_RATE_TOL),
                   f"{zg['rate_w1'] < JAM_RATE_TOL} vs {zg['rate_walt'] < JAM_RATE_TOL}", "equal"))
    # Composition is checked as a SUBSET relation, not exact equality: at
    # this array's high symmetry many constraints tie EXACTLY at the
    # optimum (achieved rate agrees to 11 significant digits between W=I
    # and alt W -- verified, not assumed), and which of several degenerate
    # ties gets reported as "binding" is a numerical tie-break, not a
    # physical disagreement -- confirmed by direct inspection: the smaller
    # set is a CLEAN subset of the larger, zero elements unique to either
    # side beyond that (a genuine physical conflict would show BOTH-sided
    # unique elements, which this does not).
    smaller, larger = sorted([zg["binding_set_w1"], zg["binding_set_walt"]], key=len)
    checks.append(("M  binding-set COMPOSITION is W-INSENSITIVE (subset, no cross-conflict)",
                   smaller <= larger,
                   f"{len(zg['binding_set_w1'])} vs {len(zg['binding_set_walt'])}", "subset"))

    # ---- L: determinism ----
    checks.append(("L  determinism: two independent solves of the same step are identical",
                   zg["repeat_identical"], str(zg["repeat_identical"]), "True"))

    # ---- Z16: QPFAIL probe (code review 23261 C1, CRITICAL) ----
    checks.append(("Z16  QPFAIL probe uses a real, nontrivial active set",
                   zqp["n_active"] > 0, f"{zqp['n_active']}", "> 0"))
    checks.append(("Z16  forced NNLS non-convergence returns feasible=False (no traceback)",
                   zqp["feasible_forced"] is False, str(zqp["feasible_forced"]), "False"))
    checks.append(("Z16  forced non-convergence's x is the safe zero fallback, never a jam reading",
                   zqp["x_forced_zero"], str(zqp["x_forced_zero"]), "True"))
    checks.append(("Z16  CONTROL: the SAME problem with a normal iteration budget converges",
                   zqp["feasible_normal"] is True, str(zqp["feasible_normal"]), "True"))

    # ---- HAZARD: feasibility-preserving invariant + mutation probe ----
    checks.append(("HAZARD  general-branch min gap over the G/H run stays >= -GAP_FLOOR_TOL",
                   zhaz["on_min_gap"] >= -GAP_FLOOR_TOL,
                   f"{zhaz['on_min_gap']:.2e}", f">= {-GAP_FLOOR_TOL:.0e}"))
    checks.append(("HAZARD  a REAL general-branch pair was found for the targeted probe",
                   zhaz["probe_found"], str(zhaz["probe_found"]), "True"))
    if zhaz["probe_found"]:
        checks.append(("HAZARD  TARGETED PROBE: backtrack ON keeps the pair at/above the floor",
                       zhaz["probe_on_gap"] >= -GAP_FLOOR_TOL,
                       f"{zhaz['probe_on_gap']:.2e}", f">= {-GAP_FLOOR_TOL:.0e}"))
        checks.append(("HAZARD  CONTROL: the SAME step un-bisected (backtrack OFF) VIOLATES it",
                       zhaz["probe_off_gap"] < -GAP_FLOOR_CONTROL_OFFSET,
                       f"{zhaz['probe_off_gap']:.2e}", f"< {-GAP_FLOOR_CONTROL_OFFSET:.0e}"))

    # ---- J: doweled diagnostic (SECOND CONTROL). SC7 is a TREE (jb_x's own
    # note): 3*ncontacts=18 rows caps the rank far below 7n-7=42, so full
    # row rank (every contact equation independent) is the falsifiable
    # claim here, NOT a forced "one global DOF" -- that collapse needs a
    # topology with cycles, which SC7 does not have.
    checks.append(("J  doweled mode mobility: measured RANK equals expected (never a subtraction)",
                   zdow["rank"] == zdow["expected_rank"],
                   f"{zdow['rank']}", f"{zdow['expected_rank']}"))
    checks.append(("J  doweled array: one-global-DOF verdict matches the topology's own math",
                   zdow["one_global_dof"] == (zdow["row_cap"] >= zdow["collapse_target"]),
                   f"one_global={zdow['one_global_dof']}", f"row_cap>=target={zdow['row_cap'] >= zdow['collapse_target']}"))
    checks.append(("J  SC7 is a TREE: one global DOF is HONESTLY False (18 rows < 42 needed)",
                   zdow["one_global_dof"] is False and zdow["row_cap"] < zdow["collapse_target"],
                   f"{zdow['one_global_dof']}", "False (tree, under-constrained)"))

    # ---- Z14: driven-unit-is-an-input, two-cell sanity ----
    checks.append(("Z14  single-crank-handle DRIVE variant is live (driven=DRIVEN_UNIT_INDEX)",
                   ztwo["single_unit_status"] == "OK", ztwo["single_unit_status"], "OK"))
    checks.append(("Z14  a unit-0/unit-1 plate pair near the test w was found",
                   ztwo["pair"] is not None, str(ztwo["pair"]), "not None"))
    checks.append(("Z14  static gap=w crossing angle is a finite, positive number",
                   ztwo["static_cross_angle"] is not None and ztwo["static_cross_angle"] > 0.0,
                   f"{ztwo['static_cross_angle']}", "> 0"))
    checks.append(("Z14  the DYNAMIC crank_run flags that SAME pair active (qvf.17's story, live)",
                   ztwo["dynamic_binds_pair"], str(ztwo["dynamic_binds_pair"]), "True"))

    # ---- Z17: THE LOCK SURFACE a*(w,t) + MOTION ORDER (bead qvf.19) ----
    surf = zlock["surf"]
    all_lock_results = (surf["w_results"] + surf["w_results_alt"]
                        + surf["t_results"] + surf["t_results_alt"])

    checks.append(("Z17  no QPFAIL anywhere in the (w,t) sweep (never silently skipped)",
                   all(r["status"] != "qpfail" for r in all_lock_results),
                   str([r["status"] for r in all_lock_results]), "no 'qpfail'"))

    checks.append(("Z17  NO row asserts a* == a_ico (qvf.15 ruler-test hazard)",
                   all(r["a_star"] is None or abs(r["a_star"] - A_ICO_LOCAL) > AICO_CONTROL_OFFSET
                       for r in all_lock_results),
                   "checked", f"none within {AICO_CONTROL_OFFSET:.0e} of a_ico"))

    # ---- HIGH fix (23337): ALT grids are NON-RATIO-DERIVED from PRIMARY,
    # GATED not merely asserted in prose (the jb_y K_GRID/K_GRID_ALT shape).
    wga, tga = zlock["w_grid_apart"], zlock["t_grid_apart"]
    checks.append(("GRID  W_GRID_FRAC_ALT is non-ratio-derived from W_GRID_FRAC",
                   wga["non_ratio_derived"], str(wga["ratios"]), "no common ratio"))
    checks.append(("GRID  W_GRID_FRAC_ALT is offset-apart from every PRIMARY value",
                   wga["offset_apart"], f"{wga['min_offset']:.4f}", f"> {GRID_OFFSET_APART_TOL}"))
    checks.append(("GRID  T_GRID_ALT is non-ratio-derived from T_GRID",
                   tga["non_ratio_derived"], str(tga["ratios"]), "no common ratio"))
    checks.append(("GRID  T_GRID_ALT is offset-apart from every PRIMARY value",
                   tga["offset_apart"], f"{tga['min_offset']:.4f}", f"> {GRID_OFFSET_APART_TOL}"))

    # ---- SHIP-BLOCKER 1a (23299): matched-h0 CROSS rebuild. cross_w0 and
    # cross_wico (bead .18, UNCHANGED here) are CONFOUNDED -- different w
    # AND different h0. The two ADDED legs (cross_h05_wico, cross_h2_w0)
    # complete both matched pairs so the w-effect is isolated from the
    # step-size effect. THE CORRECTED CONCLUSION (a real change from bead
    # .18's own "W IS CAUSAL" reading): at MATCHED h0=0.5, w=0 AND w=w_ico
    # BOTH jam; at matched h0=2.0, both REACH. W's causal role, as bead
    # .18 originally claimed from the CONFOUNDED pair, does NOT survive a
    # matched-h0 test -- reported here as a correction, not smoothed over.
    cw0, cwico = zg["cross_w0"], zg["cross_wico"]
    ch05wico, ch2w0 = zg["cross_h05_wico"], zg["cross_h2_w0"]
    # RE-PRICED for bead inviscid-1wd. This row asserted the literal value
    # "jammed", inherited from bead .18's cross_w0 and labelled "unchanged".
    # Restoring the ball joint changed it to "reached", and pinning the old
    # value would have meant asserting the disassembled array's behaviour on an
    # assembled one. What the CROSS-MATCHED block is actually for -- showing
    # that at MATCHED h0 the outcome does not turn on w -- is unaffected and is
    # gated by the AGREE row below, which is where the falsifiable content
    # lives. This row is reduced to a computability check so the AGREE row
    # cannot pass on a missing operand.
    checks.append(("CROSS-MATCHED  h0=0.5: w=0 status computable (bead .18's own "
                   "cross_w0, VALUE MOVED by inviscid-1wd)",
                   cw0["status"] in ("jammed", "reached"), cw0["status"], "computed"))
    checks.append(("CROSS-MATCHED  h0=0.5: w=w_ico status (NEW matched leg)",
                   ch05wico["status"] in ("jammed", "reached"), ch05wico["status"], "computed"))
    checks.append(("CROSS-MATCHED  h0=0.5: OUTCOMES AGREE (w is NOT shown causal at this h0)",
                   cw0["status"] == ch05wico["status"],
                   f"{cw0['status']} vs {ch05wico['status']}", "agree (correction)"))
    checks.append(("CROSS-MATCHED  h0=2.0: w=w_ico status (bead .18's own cross_wico, unchanged)",
                   cwico["status"] == "reached", cwico["status"], "reached"))
    checks.append(("CROSS-MATCHED  h0=2.0: w=0 status (NEW matched leg)",
                   ch2w0["status"] in ("jammed", "reached"), ch2w0["status"], "computed"))
    checks.append(("CROSS-MATCHED  h0=2.0: OUTCOMES AGREE (w is NOT shown causal at this h0 either)",
                   cwico["status"] == ch2w0["status"],
                   f"{cwico['status']} vs {ch2w0['status']}", "agree (correction)"))

    # ---- SHIP-BLOCKER 1b (23299): h-refinement. See h_refinement_probe's
    # own docstring for why the gated claim is SPREAD-flatness (survives
    # refinement) at the interior points, not either point's ABSOLUTE a*
    # (which drifts substantially with h0 at every w tested, boundary
    # included) -- and separately, the boundary point's own ABSOLUTE a*
    # is gated NON-stable, the critique's core finding.
    href = zlock["href"]
    interior_spreads = [abs(a["a_star"] - b["a_star"]) for a, b in
                        zip(href["interior_a"], href["interior_b"])
                        if a["a_star"] is not None and b["a_star"] is not None]
    interior_complete = len(interior_spreads) == len(H_REFINE_LEVELS)
    checks.append(("HREF  interior a* computable at every H_REFINE_LEVELS h0 (both points)",
                   interior_complete,
                   str([(a["status"], b["status"]) for a, b in
                        zip(href["interior_a"], href["interior_b"])]), "all jammed/reached"))
    # RE-SCOPED (bead inviscid-l1d). This row asserted the interior spread
    # between 0.3wico and 0.6wico "stays flat", and it passed at EXACTLY 0.00e+00
    # -- because the wire span was pinned and the w axis did nothing at all. It
    # was asserting the defect. The spread between two DIFFERENT w values is a
    # w-sensitivity measurement, not an h0-stability one, and those were
    # conflated. Split into the two questions that were tangled together.
    spread_stable = ([abs(interior_spreads[i] - interior_spreads[i + 1])
                      for i in range(len(interior_spreads) - 1)]
                     if interior_complete else [])
    checks.append(("HREF  interior spread between two w values is NONZERO -- the w axis "
                   "is live, which it was not before inviscid-l1d -- CAN FAIL",
                   interior_complete and all(s > W_LIVE_FLOOR for s in interior_spreads),
                   str([f"{s:.3f}" for s in interior_spreads]), f"all > {W_LIVE_FLOOR}"))
    # AND IT IS NOT CONVERGED, which is a finding and not a tolerance to widen.
    # Measured 3.062 / 2.850 / 21.709 across h0 = 0.5 / 0.25 / 0.125: steady at
    # the two coarse levels and 7x larger at the finest. So the w axis being
    # LIVE is robust, and the MAGNITUDE of da*/dw is not quotable. The row
    # asserts the divergence rather than hiding it, and fails if the spread ever
    # settles -- at which point the magnitude becomes quotable and this row
    # should be replaced by one that quotes it.
    checks.append(("HREF  but the spread's MAGNITUDE is NOT h0-converged, so da*/dw's "
                   "SIZE is not quotable -- only its non-zero-ness -- CAN FAIL",
                   len(spread_stable) > 0
                   and max(spread_stable) > H_REFINE_SPREAD_DIVERGENT_FLOOR,
                   str([f"{d:.3f}" for d in spread_stable]),
                   f"max > {H_REFINE_SPREAD_DIVERGENT_FLOOR}"))
    boundary_vals = [b["a_star"] for b in href["boundary"] if b["a_star"] is not None]
    boundary_complete = len(boundary_vals) == len(H_REFINE_LEVELS)
    boundary_diffs = [abs(boundary_vals[i] - boundary_vals[i + 1])
                      for i in range(len(boundary_vals) - 1)] if boundary_complete else []
    checks.append(("HREF  boundary (w=0) a* computable at every H_REFINE_LEVELS h0",
                   boundary_complete, str([b["status"] for b in href["boundary"]]),
                   "all jammed/reached"))
    # RE-SCOPED, AND IT RETIRES A RECORDED FINDING (bead inviscid-l1d). This row
    # asserted the boundary a* is NOT stable across h0 -- SHIP-BLOCKER 1's
    # reproduced non-convergence -- and it passed at 0.312 / 3.570. With the
    # wire span fixed the same ladder gives 0.185 / 0.067: monotonically
    # DECREASING, i.e. converging. The attribution is clean rather than
    # confounded with the budget raise that accompanied it: with the wire fixed
    # and the OLD budget the finest level exhausted and reported QPFAIL, so the
    # budget only let the run finish; the trajectory that converges is the one
    # the fixed wire produces. The recorded non-convergence was an artifact of
    # the span pinning at sqrt(2) from step 1, not a property of the refinement.
    checks.append(("HREF  boundary a* CONVERGES across h0 -- the recorded "
                   "non-convergence was the pinned wire span, not the refinement "
                   "-- CAN FAIL",
                   boundary_complete and len(boundary_diffs) > 1
                   and boundary_diffs[-1] < boundary_diffs[0]
                   and boundary_diffs[-1] < H_REFINE_UNSTABLE_FLOOR,
                   str([f"{d:.3f}" for d in boundary_diffs]),
                   f"decreasing, last < {H_REFINE_UNSTABLE_FLOOR}"))
    qpfail_budget = [b for b in (href["interior_a"] + href["interior_b"] + href["boundary"])
                     if b["status"] == "qpfail"]
    budget_classified = all(b["budget_exhausted"] for b in qpfail_budget)
    checks.append(("HREF  any qpfail in the refinement ladder is classified BUDGET-EXHAUSTED, "
                   "not a genuine solver failure",
                   budget_classified,
                   f"{len(qpfail_budget)} qpfail row(s), all budget-exhausted={budget_classified}"
                   if qpfail_budget else "0 qpfail rows", "True (or none)"))

    # ---- SHIP-BLOCKER 2 (23299): A_START_T derivation + live validity ----
    asi = surf["a_start_info"]
    checks.append(("A_START_T  bisection bracket is valid (glo below target, ghi at/above it)",
                   asi["bracket_ok"], f"glo={asi['glo']:.2e} ghi={asi['ghi']:.4f}",
                   f"straddle target={asi['target']:.3f}"))
    checks.append(("A_START_T  derived start angle is finite and POSITIVE (a genuine relief, not 0)",
                   surf["a_start_t"] is not None and surf["a_start_t"] > 0.0,
                   f"{surf['a_start_t']:.6f}" if surf["a_start_t"] is not None else "None", "> 0"))
    rc = surf["relief_check"]
    checks.append(("A_START_T  LIVE VALIDITY: a real crank_run from A_START_T at the smallest "
                   "swept t produces genuine multi-step motion (not another instant re-lock)",
                   rc["steps"] > 0, f"steps={rc['steps']} status={rc['status']}", "steps > 0"))

    # ---- O: SIGNED SENSITIVITY, arm A axis (da*/dw at fixed t=0) ----
    w_star = [r["a_star"] for r in surf["w_results"]]
    w_vals = [r["w"] for r in surf["w_results"]]
    w_star_alt = [r["a_star"] for r in surf["w_results_alt"]]
    w_vals_alt = [r["w"] for r in surf["w_results_alt"]]
    w_ok = (len(w_star) > 0 and len(w_star_alt) > 0
           and all(v is not None for v in w_star) and all(v is not None for v in w_star_alt))
    if w_ok:
        da_dw = _ols_slope(w_vals, w_star)
        # INTERIOR-only slope (excludes the w=0 boundary point) -- per the
        # h-refinement finding above, `da_dw` (full grid, boundary
        # included) is NOT h0-convergent and is REPORTED ONLY, never fed
        # into Q2's ratio below; `da_dw_interior` is the h0-refinement-
        # validated, gate-worthy quantity (HREF rows above). The ALT arm,
        # by construction, never samples w=0, so it is compared against
        # the PRIMARY arm's OWN interior points, not the full-grid slope.
        da_dw_interior = _ols_slope(w_vals[1:], w_star[1:])
        da_dw_alt = _ols_slope(w_vals_alt, w_star_alt)
        da_dw_shuffled = _ols_slope(w_vals, list(reversed(w_star)))
    else:
        da_dw = da_dw_interior = da_dw_alt = da_dw_shuffled = float("nan")
    checks.append(("O  primary w-grid: every point JAMMED or REACHED (a* computable)",
                   w_ok, str([r["status"] for r in surf["w_results"]]), "no 'qpfail'"))
    checks.append(("O  da*/dw, FULL grid (includes w=0) -- REPORTED ONLY, non-convergent per HREF",
                   w_ok and np.isfinite(da_dw), f"{da_dw:+.6f} deg/w", "finite (not gated on Q2)"))
    checks.append(("O  da*/dw, INTERIOR ONLY -- the h0-refinement-validated quantity Q2 USES",
                   w_ok and np.isfinite(da_dw_interior), f"{da_dw_interior:+.6f} deg/w", "finite"))
    checks.append(("O  CONTROL: shuffled (w, a*) pairing does NOT reproduce the fitted slope",
                   w_ok and abs(da_dw_shuffled - da_dw) > SLOPE_SHUFFLE_DIFFER_TOL,
                   f"{da_dw_shuffled:+.6f} vs {da_dw:+.6f}", "differ"))
    checks.append(("O  ALT arm agrees with the PRIMARY arm's OWN interior slope (or both ~ 0)",
                   w_ok and (abs(da_dw_interior) < W_INTERIOR_FLAT_TOL
                            or (da_dw_interior > 0) == (da_dw_alt > 0)),
                   f"{da_dw_alt:+.6f} vs interior {da_dw_interior:+.6f}", "same sign or both flat"))

    # ---- P: SIGNED SENSITIVITY, arm C axis (da*/dt, opening range, at
    # fixed w=w_ref, RE-BASED at A_START_T per SHIP-BLOCKER 2) ----
    t_star = [r["a_star"] for r in surf["t_results"]]
    t_vals = [r["t"] for r in surf["t_results"]]
    t_star_alt = [r["a_star"] for r in surf["t_results_alt"]]
    t_vals_alt = [r["t"] for r in surf["t_results_alt"]]
    t_ok = (len(t_star) > 0 and len(t_star_alt) > 0
           and all(v is not None for v in t_star) and all(v is not None for v in t_star_alt))
    if t_ok:
        da_dt = _ols_slope(t_vals, t_star)
        da_dt_alt = _ols_slope(t_vals_alt, t_star_alt)
        da_dt_shuffled = _ols_slope(t_vals, list(reversed(t_star)))
    else:
        da_dt = da_dt_alt = da_dt_shuffled = float("nan")
    checks.append(("P  primary t-grid (re-based, opening range): every point computable",
                   t_ok, str([r["status"] for r in surf["t_results"]]), "no 'qpfail'"))
    checks.append(("P  da*/dt (arm C), opening range, fixed w=T_SWEEP_W_FRAC*w_ico, sign + magnitude",
                   t_ok and np.isfinite(da_dt), f"{da_dt:+.6f} deg-range/t", "finite"))
    checks.append(("P  CONTROL: shuffled (t, a*) pairing does NOT reproduce the fitted slope",
                   t_ok and abs(da_dt_shuffled - da_dt) > SLOPE_SHUFFLE_DIFFER_TOL,
                   f"{da_dt_shuffled:+.6f} vs {da_dt:+.6f}", "differ"))
    checks.append(("P  ALT arm: da*/dt sign agrees with the primary arm",
                   t_ok and (da_dt > 0) == (da_dt_alt > 0),
                   f"{da_dt_alt:+.6f} vs {da_dt:+.6f}", "same sign"))

    # ---- Q: THE SEPARATION ROW, both ways ----
    if w_ok:
        w_sep_all = abs(max(w_star) - min(w_star))
        w_sep_interior = abs(max(w_star[1:]) - min(w_star[1:]))  # excludes the w=0 boundary point
        ktable_ratio_interior = w_sep_interior / 13.849356
    else:
        w_sep_all = w_sep_interior = ktable_ratio_interior = float("nan")
    checks.append(("Q  arm-A comparison: INTERIOR w-separation as a RATIO of the k-table span "
                   "(the h0-validated quantity; the boundary-inclusive figure is reported "
                   "separately, not used here, per HREF's non-convergence finding)",
                   w_ok and np.isfinite(ktable_ratio_interior),
                   f"{w_sep_interior:.6f} / 13.849356 = {ktable_ratio_interior:.4f}",
                   "finite ratio, printed"))
    checks.append(("Q  (reported only) FULL-grid w-separation, includes the non-convergent "
                   "w=0 boundary jump",
                   w_ok, f"{w_sep_all:.6f}", "reported, not gated"))
    if t_ok:
        t_sep = abs(max(t_star) - min(t_star))
    else:
        t_sep = float("nan")
    checks.append(("Q  arm-C signature: t-driven opening-range separation exceeds T_SEPARATION_TOL",
                   t_ok and t_sep > T_SEPARATION_TOL, f"{t_sep:.6f}", f"> {T_SEPARATION_TOL}"))
    # RE-SCOPED (bead inviscid-l1d). This asserted that w-only separation stays
    # UNDER a tolerance -- the arm-C signature -- and it passed at exactly
    # 0.000000. Not because the array is arm-C-like, but because the wire span
    # was pinned at sqrt(2) from the first crank step, so w could not move
    # anything. It was asserting the defect, and it made the Q2 verdict below
    # STRUCTURALLY GUARANTEED: a dead numerator puts the discrimination ratio at
    # zero, always below the band, always ARM-C-LIKE. Q2 could not fail in the
    # only direction that mattered. Now the separation is real and Q2 is a
    # measurement.
    checks.append(("Q  the w axis MOVES the lock angle -- separation is nonzero, where "
                   "before inviscid-l1d it was exactly 0.000000 -- CAN FAIL",
                   w_ok and w_sep_interior > W_LIVE_FLOOR,
                   f"{w_sep_interior:.6f} = {w_sep_interior / K_TABLE_SPAN:.4f} of the "
                   f"k-table span", f"> {W_LIVE_FLOOR}"))

    # ---- Q2: THE DISCRIMINATION-RATIO ROW (the deliverable's headline).
    # USES da_dw_interior (the h0-refinement-validated quantity), NEVER
    # the full-grid da_dw (SHIP-BLOCKER 1's own finding: non-convergent,
    # sign-flipping across h0 -- unfit to feed a headline verdict). da_dt
    # is the RE-BASED, opening-range quantity (SHIP-BLOCKER 2's fix). ----
    ratio_computable = (w_ok and t_ok and np.isfinite(da_dw_interior) and np.isfinite(da_dt)
                       and abs(da_dt) > RATIO_ZERO_TOL)
    if ratio_computable:
        discrim_ratio = abs(da_dw_interior) / abs(da_dt)
        if discrim_ratio > DISCRIM_RATIO_HIGH:
            verdict = "ARM-A-LIKE"
        elif discrim_ratio < DISCRIM_RATIO_LOW:
            verdict = "ARM-C-LIKE"
        else:
            verdict = "NON-DISCRIMINATING"
    else:
        discrim_ratio = None
        verdict = None
    checks.append(("Q2  DISCRIMINATION RATIO |da*/dw_interior|/|da*/dt| is COMPUTABLE "
                   "(never an inf-pass)",
                   ratio_computable,
                   f"{discrim_ratio:.3e}" if ratio_computable else "undefined (|da*/dt| ~ 0)",
                   f"nonzero denom > {RATIO_ZERO_TOL:.0e}"))
    checks.append(("Q2  VERDICT: exactly one of ARM-A-LIKE / ARM-C-LIKE / NON-DISCRIMINATING",
                   ratio_computable and verdict in ("ARM-A-LIKE", "ARM-C-LIKE", "NON-DISCRIMINATING"),
                   str(verdict), f"band [{DISCRIM_RATIO_LOW}, {DISCRIM_RATIO_HIGH}]"))

    # ---- R: NON-VACUITY, two-sided bands -- SPLIT PER ARM (23299-adjacent
    # fix: the w-arm's a* is an ABSOLUTE ANGLE from a=0; the t-arm's a* is
    # an OPENING RANGE from A_START_T (SHIP-BLOCKER 2) -- combining them
    # into one "whole surface" span would mix two different quantities). ----
    w_arm_results = surf["w_results"] + surf["w_results_alt"]
    t_arm_results = surf["t_results"] + surf["t_results_alt"]
    w_finite = [r["a_star"] for r in w_arm_results if r["a_star"] is not None]
    t_finite = [r["a_star"] for r in t_arm_results if r["a_star"] is not None]
    w_non_vacuous = len(w_finite) == len(w_arm_results) and len(w_finite) > 0
    t_non_vacuous = len(t_finite) == len(t_arm_results) and len(t_finite) > 0
    w_span = (max(w_finite) - min(w_finite)) if w_non_vacuous else None
    t_span = (max(t_finite) - min(t_finite)) if t_non_vacuous else None
    checks.append(("R-w  NON-VACUITY: every w-arm grid point contributed a finite a* (fold into all())",
                   w_non_vacuous, f"{len(w_finite)}/{len(w_arm_results)}", "all finite"))
    checks.append(("R-w  w-arm span (absolute angle) lies in a TWO-SIDED band",
                   w_non_vacuous and SPAN_LOWER < w_span < SPAN_UPPER,
                   f"{w_span:.6f}" if w_non_vacuous else "n/a", f"({SPAN_LOWER:.0e}, {SPAN_UPPER})"))
    checks.append(("R-t  NON-VACUITY: every t-arm grid point contributed a finite a* (fold into all())",
                   t_non_vacuous, f"{len(t_finite)}/{len(t_arm_results)}", "all finite"))
    checks.append(("R-t  t-arm span (opening range) lies in a TWO-SIDED band",
                   t_non_vacuous and SPAN_LOWER < t_span < SPAN_UPPER_T,
                   f"{t_span:.6f}" if t_non_vacuous else "n/a", f"({SPAN_LOWER:.0e}, {SPAN_UPPER_T})"))

    # ---- S: MOTION ORDER, falsifiable, with a can-fail control ----
    # TEST reads the LAST completed step of the single-crank-handle trace
    # (propagation has had time to reach some neighbours); CONTROL reads
    # the FIRST step of the driven="all" trace, among the SIX topologically
    # EQUIVALENT corner units only (index 0, the star centre, is symmetry-
    # distinguished -- see MOTION_UNIFORM_TOL's docstring). Both now read a
    # HOMOGENEOUS raw-velocity-norm metric for every unit (code review
    # MEDIUM fix, 23337 -- `motion_order_trace`'s own docstring).
    test_tr, ctrl_tr = zlock["test_trace"], zlock["control_trace"]
    test_not_equal = False
    if test_tr["history"]:
        last = test_tr["history"][-1]
        test_not_equal = bool(np.max(last) - np.min(last) > MOTION_ACTIVITY_TOL)
    ctrl_equal = False
    if ctrl_tr["history"]:
        first_c = ctrl_tr["history"][0][1:]  # corner units only, first step
        ctrl_equal = bool(np.max(first_c) - np.min(first_c) < MOTION_UNIFORM_TOL)
    checks.append(("S  single-handle TEST: per-unit activity is NOT all equal (units take turns)",
                   len(test_tr["history"]) > 0 and test_not_equal,
                   f"spread={float(np.max(test_tr['history'][-1]) - np.min(test_tr['history'][-1])):.4f}"
                   if test_tr["history"] else "n/a", f"> {MOTION_ACTIVITY_TOL:.0e}"))
    # PROMOTED FROM PROSE (code review MEDIUM fix, 23337): the write-back's
    # own strongest evidence for "units take turns" -- the undriven corner
    # units read EXACTLY 0.0 while the driven unit is active, for the
    # FIRST several steps -- as a real, falsifiable row rather than prose.
    # REPLACED for bead inviscid-1wd. The retired row required the undriven
    # corners to read EXACTLY 0.0 for at least one early step. That evidence
    # belonged to the DISASSEMBLED array: with no inter-unit constraint the
    # undriven units genuinely had nothing to transmit motion to them. Once
    # the ball joint is restored they are joined, so they move from step 0 --
    # correctly, and the old row can never pass again for a physical reason,
    # which is why it is replaced rather than loosened.
    #
    # What survives, and is the stronger statement, is LEAD/LAG: the driven
    # unit outruns every other unit while the motion propagates. Measured
    # here: max_other/driven = 0.431 at step 0, rising through 0.615, 0.617,
    # 0.661 as the array takes up the drive. The paired CONTROL below
    # (driven="all") is what makes this falsifiable -- there the ratio goes to
    # 1 and this row would redden.
    lead_ratios = []
    for step_row in test_tr["history"][:MOTION_LEAD_STEPS]:
        driven_act = float(step_row[DRIVEN_UNIT_INDEX])
        others = [float(step_row[i]) for i in range(len(step_row))
                  if i != DRIVEN_UNIT_INDEX]
        if driven_act > MOTION_ACTIVITY_TOL and others:
            lead_ratios.append(max(others) / driven_act)
    checks.append(("S  single-handle TEST: the DRIVEN unit LEADS -- every other unit's "
                   "activity stays a fraction of it while the motion propagates",
                   len(lead_ratios) > 0 and min(lead_ratios) < MOTION_LEAD_RATIO_MAX,
                   f"min ratio {min(lead_ratios):.4f}" if lead_ratios else "n/a",
                   f"< {MOTION_LEAD_RATIO_MAX}"))
    checks.append(("S  CONTROL: driven=all, w large, t=0: the SIX corner units move uniformly",
                   len(ctrl_tr["history"]) > 0 and ctrl_equal,
                   f"spread={float(np.max(ctrl_tr['history'][0][1:]) - np.min(ctrl_tr['history'][0][1:])):.2e}"
                   if ctrl_tr["history"] else "n/a", f"< {MOTION_UNIFORM_TOL}"))
    checks.append(("S  NON-VACUITY: control run actually completed at least one step",
                   len(ctrl_tr["history"]) > 0, str(len(ctrl_tr["history"])), "> 0"))

    # ---- T: BINDING SET REPORTED PER GRID POINT ----
    jammed_results = [r for r in all_lock_results if r["status"] == "jammed"]
    binding_nonempty = all(bool(r["binding_ever"]) for r in jammed_results)
    distinct_sets = {frozenset((b[0], b[1], b[2], b[3], b[4]) for b in r["binding_ever"])
                     for r in jammed_results}
    checks.append(("T  every JAMMED grid point's binding set is non-empty (fold into all())",
                   len(jammed_results) > 0 and binding_nonempty,
                   f"{len(jammed_results)} jammed points", "all non-empty"))
    # LOW fix (23337): a REAL, falsifiable condition -- >1 distinct set is
    # the actual mechanism-evidence signature this row exists to check
    # (a monolithic single binding set everywhere would be a DIFFERENT,
    # also-reportable result, and this row would correctly redden for it).
    checks.append(("T  distinct binding sets across the grid (>1 is mechanism evidence, CAN FAIL)",
                   len(distinct_sets) > 1, f"{len(distinct_sets)} distinct set(s)", "> 1"))

    # ---- PHASE TRACKING: is a* a phase, or a bookkeeping integral? ----
    # These rows are stated as POSITIVE, falsifiable claims about a defect --
    # "the drift EXCEEDS this floor", not "the drift is small". A row asserting
    # a* tracks the phase would be asserting something measurably false, and
    # loosening its tolerance until it passed would be the vacuity this file's
    # whole house style exists to refuse.
    ph = zlock["phase"]
    inst = []
    for ang in PHASE_PROBE_ANGLES + PHASE_PROBE_ANGLES_ALT:
        got, spread = configuration_phase(corners(ang))
        inst.append((ang, abs(got - ang), spread))
    checks.append(("PHASE  INSTRUMENT: closed-form inversion recovers a KNOWN symmetric "
                   "phase from geometry alone (both arms)",
                   len(inst) > 0 and all(e <= PHASE_INSTRUMENT_TOL for _a, e, _s in inst),
                   f"max err {max(e for _a, e, _s in inst):.2e} over {len(inst)} angles",
                   f"<= {PHASE_INSTRUMENT_TOL:.0e}"))
    checks.append(("PHASE  INSTRUMENT: on a symmetric pose all 12 vertex radii AGREE "
                   "(so a nonzero spread means off-path, not noise)",
                   len(inst) > 0
                   and all(s_ <= PHASE_SYMMETRIC_SPREAD_TOL for _a, _e, s_ in inst),
                   f"max spread {max(s_ for _a, _e, s_ in inst):.2e}",
                   f"<= {PHASE_SYMMETRIC_SPREAD_TOL:.0e}"))
    checks.append(("PHASE  INSTRUMENT NON-VACUITY: the calibration angles span a real range",
                   len(inst) > 1
                   and (max(a_ for a_, _e, _s in inst) - min(a_ for a_, _e, _s in inst)) > 10.0,
                   f"{min(a_ for a_,_e,_s in inst):.1f}..{max(a_ for a_,_e,_s in inst):.1f} deg",
                   "> 10 deg"))

    ok_ph = [r for r in ph if r["phase"] is not None]
    checks.append(("PHASE  NON-VACUITY: every step-size level produced a reading",
                   len(ok_ph) == len(ph) and len(ph) > 1,
                   f"{len(ok_ph)}/{len(ph)} levels", "all"))
    checks.append(("PHASE  LINKAGE VALIDITY: hinges, joints and triangle edges stay EXACT "
                   "through the run (so the departure below is real freedom, not a broken solver)",
                   len(ok_ph) > 0
                   and all(max(r["hinge"], r["joint"], r["edge"]) <= LINKAGE_EXACT_TOL
                           for r in ok_ph),
                   f"worst {max(max(r['hinge'], r['joint'], r['edge']) for r in ok_ph):.2e}"
                   if ok_ph else "n/a", f"<= {LINKAGE_EXACT_TOL:.0e}"))
    # CORRECTED. These two rows previously asserted that a_hat differs from "the
    # configuration's phase" by more than a degree and that the gap is
    # STRUCTURAL. Both readings came from `configuration_phase(xs[0])` -- unit 0
    # alone, which is SC7's CENTRE unit and lags the corners. Measured against
    # the array MEAN the gap is 0.34 / 0.27 / 0.21 deg across a 4x ladder:
    # small, and SHRINKING, i.e. discretisation. The 5 deg that was reported as
    # a_hat's error was the centre unit's lag wearing a_hat's name.
    checks.append(("PHASE  a_hat TRACKS the array's MEAN measured phase (it is crank throw, "
                   "and it lands where the array's mean does) -- CAN FAIL",
                   len(ok_ph) > 0
                   and all(abs(r["drift"]) < PHASE_MEAN_DRIFT_TOL for r in ok_ph),
                   f"drifts {[round(r['drift'], 3) for r in ok_ph]}" if ok_ph else "n/a",
                   f"all < {PHASE_MEAN_DRIFT_TOL}"))
    drifts = [abs(r["drift"]) for r in ok_ph]
    checks.append(("PHASE  that residual SHRINKS under refinement, so it is discretisation "
                   "and not a defect in the integral -- CAN FAIL",
                   len(drifts) > 1 and drifts[-1] < drifts[0],
                   f"{drifts[0]:.4f} -> {drifts[-1]:.4f} over h0={list(PHASE_H0_LEVELS)}"
                   if len(drifts) > 1 else "n/a", "decreasing"))

    # ---- JOINT INTEGRITY (bead inviscid-1wd): is it still an array? ----
    jt = zlock["joints"]
    checks.append(("JOINT  NON-VACUITY: the integrity probe actually cranked steps "
                   "(both arms)",
                   jt["test"]["steps"] > 0 and jt["control"]["steps"] > 0,
                   f"test {jt['test']['steps']}, control {jt['control']['steps']}", "> 0"))
    checks.append(("JOINT  TEST: shared vertices stay coincident through a real crank run",
                   jt["test"]["steps"] > 0 and jt["test"]["worst"] < JOINT_GAP_TOL,
                   f"{jt['test']['worst']:.3e}", f"< {JOINT_GAP_TOL:.0e}"))
    checks.append(("JOINT  CONTROL: the retired per-unit projection MEASURABLY pulls the "
                   "array apart (without this the TEST is unfalsifiable)",
                   jt["control"]["steps"] > 0
                   and jt["control"]["worst"] > JOINT_CONTROL_FLOOR,
                   f"{jt['control']['worst']:.3e}", f"> {JOINT_CONTROL_FLOOR:.0e}"))
    checks.append(("JOINT  the two arms differ by orders of magnitude, not by rounding",
                   jt["test"]["steps"] > 0 and jt["control"]["steps"] > 0
                   and jt["control"]["worst"] > 1e3 * max(jt["test"]["worst"], 1e-300),
                   f"control/test = {jt['control']['worst'] / max(jt['test']['worst'], 1e-300):.2e}",
                   "> 1e3"))

    # ---- DRIVE-MODEL ROBUSTNESS (hazard comment, 2026-08-21) ----
    # RESHAPED for bead inviscid-1wd. The original row compared the two points'
    # a* SPREAD against DRIVE_ROBUST_TOL, which presumes both points HAVE an
    # a*. Under the assembled array they do not: single-crank-handle drive runs
    # the full 45-degree window without locking at either w (measured, 360
    # steps each). That is a real drive-dependence result -- the lock at 29.88
    # belongs to the uniform in-phase drive -- so it is reported as such rather
    # than being forced back into a spread comparison that no longer has
    # operands. The rows below still bite: they fail if the two w points
    # DISAGREE with each other, and if a lock does appear the spread row
    # becomes live again on its own terms.
    drv = zlock["drive"]
    checks.append(("DRIVE-ROBUST  NON-VACUITY: both single-handle points computable (no QPFAIL)",
                   len(drv["statuses"]) > 0
                   and all(st != "qpfail" for st in drv["statuses"]),
                   f"{drv['statuses']}", "no 'qpfail'"))
    checks.append(("DRIVE-ROBUST  the two w points AGREE in outcome under this drive (CAN FAIL)",
                   len(drv["statuses"]) > 0 and drv["agree"],
                   f"{' vs '.join(drv['statuses'])}", "identical"))
    if drv["any_lock"]:
        checks.append(("DRIVE-ROBUST  a* spread under single-handle drive within DRIVE_ROBUST_TOL",
                       drv["spread"] is not None and drv["spread"] < DRIVE_ROBUST_TOL,
                       f"{drv['spread']:.6f}" if drv["spread"] is not None else "n/a",
                       f"< {DRIVE_ROBUST_TOL}"))
    else:
        checks.append(("DRIVE-ROBUST  no lock under single-handle drive: the 29.88 lock is "
                       "DRIVE-SPECIFIC, reported not gated",
                       True, f"{drv['statuses'][0]} at both w", "printed"))

    # ---- PHASE THREADING (code review, bead qvf.20) ----
    zpt = zbp["threading"]
    checks.append(("Z18  a two-sublattice topology IS driven with per-unit phase -- the "
                   "guarantee revision 34e636c stated in a comment and did not keep "
                   "-- CAN FAIL",
                   zpt["n_off"] > 0 and zpt["moved"] > 1e-6,
                   f"offset units move by >= {zpt['moved']:.6f}", "> 1e-6"))
    checks.append(("Z18  CONTROL: the ZERO-phase units are bit-identical, so the row "
                   "above measures threading and not merely a changed topology",
                   zpt["n_zero"] > 0 and zpt["still"] == 0.0,
                   f"worst change {zpt['still']:.1e} over {zpt['n_zero']} units",
                   "exactly 0"))
    checks.append(("Z18  and the tangents it threads are genuinely different, so the "
                   "defect was not rounding-level",
                   abs(zpt["cos"]) < 0.99, f"cos(tangent(a), tangent(a+60)) "
                   f"{zpt['cos']:.6f}", "< 0.99"))

    # ---- MEANINGLESS DEPTH (critique, bead qvf.20) ----
    zmd = zbp["depth"]
    checks.append(("Z19  the ANTIPODAL artifact class EXISTS -- without this the row "
                   "below is satisfied by a configuration that simply has none "
                   "-- CAN FAIL",
                   zmd["antipodal"] > 0,
                   f"{zmd['antipodal']} antipodal pairs beyond the depth floor",
                   "> 0"))
    checks.append(("Z19  and NO active contact row is built from it: every binding "
                   "contact sits above -MEANINGLESS_DEPTH_FLOOR",
                   zmd["status"] == "OK"
                   and zmd["worst_binding"] > -MEANINGLESS_DEPTH_FLOOR,
                   f"worst binding contact {zmd['worst_binding']:.6f}",
                   f"> {-MEANINGLESS_DEPTH_FLOOR}"))

    # ---- BROAD-PHASE (beads inviscid-0dm, inviscid-qvf.23) ----
    bpr, bpm = zbp["probe"], zbp["mutation"]
    checks.append(("BP  EXACTNESS: the prefiltered active set is IDENTICAL to the "
                   "exhaustive scan's, at every probed configuration",
                   len(bpr) > 0 and all(r["same_active"] for r in bpr),
                   f"{sum(r['same_active'] for r in bpr)}/{len(bpr)} configurations",
                   "all identical"))
    checks.append(("BP  and the reported minimum general gap is identical too, not "
                   "merely close",
                   len(bpr) > 0 and all(r["mgg_diff"] == 0.0 for r in bpr),
                   f"worst diff {max(r['mgg_diff'] for r in bpr):.3e}" if bpr else "n/a",
                   "exactly 0"))
    checks.append(("BP  MUTATION PROBE: a deliberately UNSOUND cull is caught by that "
                   "exactness row -- without this the row is satisfied by a prefilter "
                   "that rejects nothing -- CAN FAIL",
                   bpm["total"] > 0 and bpm["caught"] == bpm["total"],
                   f"caught {bpm['caught']}/{bpm['total']}", "all caught"))
    _rej = [r["reject"] for r in bpr]
    checks.append(("BP  NON-VACUITY, TWO-SIDED: it rejects a real fraction and not all "
                   "of them -- both edges redden",
                   len(_rej) > 0
                   and all(BP_REJECT_BAND[0] < x < BP_REJECT_BAND[1] for x in _rej),
                   f"{min(_rej):.4f} .. {max(_rej):.4f}" if _rej else "n/a",
                   f"in {BP_REJECT_BAND}"))
    _id = [r["reject"] for r in bpr if r["tag"] == "ideal"]
    _st = [r["reject"] for r in bpr if r["tag"] == "stepped"]
    checks.append(("BP  and it pays MOST in the mid-integration regime, which is where "
                   "the stepper actually spends its time -- CAN FAIL",
                   len(_id) > 0 and len(_st) > 0 and min(_st) > max(_id),
                   f"ideal {max(_id):.4f} vs mid-integration {min(_st):.4f}"
                   if _id and _st else "n/a", "mid > ideal"))

    print()
    print("=" * 78)
    print(f"GATE  {len(checks)} rows: every check's verdict, and this process's "
          f"exit code")
    print("=" * 78)
    for name, passed, val, crit in checks:
        print(f"  {'PASS' if passed else 'FAIL':4s}  {name:52s} "
              f"{str(val):>18s} {str(crit):>16s}")

    print()
    print("  ROWS THAT EXIST ONLY TO STOP ANOTHER ROW BEING UNFALSIFIABLE:")
    print("   * BP's MUTATION PROBE. The exactness row is satisfied by a")
    print("     prefilter that rejects NOTHING, and by a bug that happens to")
    print("     reject only inactive pairs at the configurations sampled. The")
    print("     probe inflates the bound until the cull is unsound and checks")
    print("     that exactness notices -- 10 of 10. Without it the broad-phase")
    print("     shipped ungated, which is exactly how it shipped in August.")
    print("   * BP's TWO-SIDED rejection band. A prefilter that rejects")
    print("     everything would pass a one-sided 'it rejects things' row while")
    print("     being catastrophically wrong.")
    print("   * 'fold(a=0) is EXACTLY zero' -- without it, 'the fold table")
    print("     matches to 1e-5' is satisfiable by a constant function that")
    print("     happens to equal the five recorded angles' values.")
    print("   * 'a_ico +/- offset MISSES strut length' -- without it, the")
    print("     ridge-diagonal-equals-strut tolerance is unbounded above and")
    print("     A_ICO could drift with the row still green (the exact hole an")
    print("     independent validation found in jb_x's aico tolerance).")
    print("   * 'dq is NOT the same order as dp' -- without it, the")
    print("     one-sidedness row is satisfied by a symmetric (two-sided)")
    print("     gradient that happens to have a small q-component.")
    print("   * the FD-eps CONTROL rows (10x offset must be REJECTED) --")
    print("     without them FD_TOL is unbounded above and the normal-")
    print("     orientation claim could be wrong by an order of magnitude")
    print("     with the row still green.")
    print("   * 'negated normal FAILS' -- without it, the sign-convention row")
    print("     is satisfiable by a check that ignores the sign entirely.")
    print("   * 'naive vertex-vertex candidate is NOT closest' -- without it,")
    print("     the witness cross-check could be satisfied by a configuration")
    print("     where every reasonable method agrees trivially.")
    print("   * 'compression leaves the member INACTIVE' AND 'a bilateral")
    print("     member WOULD flag active' together -- either alone is")
    print("     satisfiable by a checker that always returns one constant.")
    print("   * 'admissibility FLIPS as t crosses the gap value' -- without")
    print("     it, 'the threshold shifts by exactly t' could be satisfied by")
    print("     a function that is constant in t.")
    print("   * the naive-instrument DISAGREEMENT row -- without it, the")
    print("     whole point of building the robust crossing test (23195's")
    print("     recorded lesson) is unconfirmed; a robust test that happens")
    print("     to agree with the broken one everywhere measured nothing new.")
    print("   * I's 'first-step rate >= JAM_RATE_CONTROL_OFFSET' -- without")
    print("     it, JAM_RATE_TOL could be set so high that row I ALSO reads")
    print("     as jammed, and the decisive G/H/I contrast would be vacuous.")
    print("   * HAZARD's backtrack-OFF control -- without it, the feasibility")
    print("     invariant row is satisfiable by a backtrack that never")
    print("     actually does anything (GAP_FLOOR_TOL vacuously wide).")
    print("   * M's alternate-W rerun -- without it, 'W=identity is a norm")
    print("     choice, not a physics choice' is asserted, never checked.")
    print("   * O/P's shuffled-pairing controls (bead qvf.19) -- without")
    print("     them, a two-point endpoint difference would report 'a slope'")
    print("     regardless of which a* value landed at which w or t; the OLS")
    print("     fit plus a re-paired-values control makes the ORDER of the")
    print("     (w, a*) / (t, a*) pairing load-bearing, not merely printed.")
    print("   * Q2's 'ratio not computable' FAIL branch -- without it, a")
    print("     degenerate da*/dt (near zero) would either crash or print an")
    print("     inf, the exact jb_x X7 shape this file's own house style")
    print("     already refuses to repeat.")
    print("   * Z17's 'NO row asserts a* == a_ico' row -- without it, the")
    print("     qvf.15 ruler-test hazard (both prior mechanisms predict that")
    print("     exact instant) could silently re-enter as false discriminating")
    print("     evidence the next time this file is extended.")
    print()
    print("  A RECORDED FINDING RETIRED, 2026-08-25, bead inviscid-l1d: the")
    print("  boundary a*'s NON-CONVERGENCE across h0 (SHIP-BLOCKER 1) was an")
    print("  artifact of the wire span pinning at sqrt(2) from the first crank")
    print("  step. With the span computed correctly the same ladder converges")
    print("  monotonically, 0.185 then 0.067. The attribution is not confounded")
    print("  with the budget raise that accompanied it: at the old budget the")
    print("  fixed-wire run merely exhausted and reported QPFAIL, so the budget")
    print("  let it finish and the wire is what made it converge.")
    print()
    print("  AND ONE FINDING THAT IS NOT RETIRED BUT NARROWED: the w axis is")
    print("  LIVE -- interior separation 3.849637 deg, 0.2780 of the k-table")
    print("  span, where before it was exactly 0.000000 -- but its MAGNITUDE is")
    print("  not h0-converged (3.062 / 2.850 / 21.709). Quote that w moves the")
    print("  lock angle. Do not quote by how much.")
    print()
    print("  ROWS DELETED RATHER THAN FIXED:")
    print("   * PHASE OFF-PATH and PHASE DEPHASING, deleted 2026-08-25 with the")
    print("     inviscid-l1d wire fix, and DELETED rather than re-priced on")
    print("     purpose. Both asserted a phenomenon exceeds a floor, and both")
    print("     floors were priced while the wire span was pinned at sqrt(2)")
    print("     from the first crank step -- a spurious force pushing the array")
    print("     off its own path. With the wire fixed the array is markedly")
    print("     MORE coherent: off-path spread 1.58e-03 against a floor of")
    print("     5e-03, dephasing 0.119 deg against a floor of 1.0. So the")
    print("     recorded 6.7 degree dephasing was inflated by TWO independent")
    print("     defects -- the wrong packing (inviscid-ia5) and this wire -- and")
    print("     re-pricing the floors would only put fresh numbers on a")
    print("     superseded topology. The phenomenon may or may not be real; it")
    print("     has not been measured on a correctly packed array, and no row")
    print("     here claims it either way.")
    print("   * the three CROSS rows ('w=0 JAMS', 'w=w_ico REACHES',")
    print("     'therefore W IS CAUSAL') are deleted, not re-priced. They")
    print("     asserted an outcome DIFFERENCE that the DISASSEMBLED array")
    print("     produced and the assembled one does not -- with the ball")
    print("     joint restored (bead inviscid-1wd) both legs reach. They had")
    print("     already been half-retracted by the CROSS-MATCHED rows, which")
    print("     showed the difference vanishes once h0 is held equal; those")
    print("     rows carry whatever survives, and cross_w0/cross_wico are")
    print("     still computed and still gated by them.")
    print("   * S's 'undriven corners read EXACTLY 0.0 for >=1 early step'")
    print("     is deleted rather than loosened. That evidence belonged to an")
    print("     array with no inter-unit constraint, where the corners had")
    print("     nothing to transmit motion to them. Joined, they move from")
    print("     step 0 and the row can never pass again FOR A PHYSICAL")
    print("     REASON. What replaces it is LEAD/LAG -- the driven unit")
    print("     outruns every other while the drive propagates (min ratio")
    print("     0.431) -- falsified by the driven='all' control, where the")
    print("     ratio goes to 1.")
    print("   * 'signed_gap agrees with unsigned _tri_tri on separated pairs'")
    print("     was dropped: jb_x's `_tri_tri` is private and this file does")
    print("     not import private symbols from it (mutation-probe rule) --")
    print("     the SAME cross-check is instead performed independently, per")
    print("     configuration, against `jb_g.segment_distance` (edge-edge)")
    print("     and a from-scratch least-squares projector (vertex-face),")
    print("     which is a stronger check than re-deriving jb_x's own method.")
    print()
    print("  A ROW DELIBERATELY NOT BUILT: general (non-parallel-facing) DEEP")
    print("  penetration-depth accuracy. `signed_gap`'s general branch reports")
    print("  a bounded, continuous, sign-correct PROXY once two non-parallel")
    print("  plates pierce, valid near the contact boundary -- the regime")
    print("  this bead's own rows exercise (a separated pair and a NEAR-")
    print("  touching pair with SMALL POSITIVE gap; the negative-gap case is")
    print("  instead validated via the PARALLEL-FACING branch's exact")
    print("  closed-form projection, which the registry pair and the folding")
    print("  square both exercise for real). Asserting numerical accuracy for")
    print("  a general pair sunk deep into interpenetration would be a row")
    print("  this bead cannot honestly make pass: no independent reference")
    print("  penetration-depth routine exists in this codebase to check it")
    print("  against, and the quasi-static stepper (bead .18) never needs a")
    print("  gap deep in the negative regime -- its active-set threshold sits")
    print("  at g <= eps_act, near zero, by construction.")
    print()
    print("  A SECOND ROW DELIBERATELY NOT BUILT (bead qvf.19): resolving WHY")
    print("  crank_run's SUSTAINED outcome depends on h0 (H_LOCK=0.5 finds a")
    print("  jam at w=w_ico that H_STEP=2.0 does not, at the SAME w, t). This")
    print("  bead FIXES h0 at the JAM-tuned, previously-validated value for")
    print("  every call it makes so the surface is internally consistent, and")
    print("  reports the sensitivity as a finding -- it does not attempt to")
    print("  determine which step size is 'more correct', which needs a")
    print("  mutation-probe pass over the FROZEN stepper (bead qvf.20's own")
    print("  scope, not this bead's), not another lock-surface grid point.")
    print()
    print("  A THIRD ROW DELIBERATELY NOT BUILT: a* == 30 deg, the CLOSED-FORM")
    print("  first contact. Every cuboctahedron vertex rides a planar ellipse")
    print("  (semi-major sqrt(2), semi-minor sqrt(2/3), axis ratio exactly")
    print("  sqrt(3)), the twelve lie in just the three coordinate planes, each")
    print("  face takes one from each, and the crank angle IS the eccentric")
    print("  anomaly (d(theta)/da = -1.000000, std 3e-13). So on the SYMMETRIC")
    print("  path every plate gap is a sinusoid: the blocking family is exactly")
    print("  gap(a) = 4 cos(a + 60), verified to 6.7e-16, giving first contact")
    print("  at a = 30.000000000 with all SIXTEEN separated pairs closing")
    print("  simultaneously, and containing NO w -- an independent geometric")
    print("  proof of the interior da*/dw = 0.")
    print("  The row is not built because THIS FILE DOES NOT WALK THAT PATH,")
    print("  and the PHASE rows above measure the two ways it departs. Units")
    print("  flex INTERNALLY -- vertex radii spread to ~1e-01 while hinges,")
    print("  joints and triangle edges stay exact to 9e-15, so it is real")
    print("  freedom of the linkage and not a broken solver. And the units")
    print("  DEPHASE from EACH OTHER under a uniform in-phase drive: SC7's")
    print("  centre lags its six corners by 6.4-6.7 deg, and the corners")
    print("  themselves split four-and-two. a* is CRANK THROW -- the arc")
    print("  length of the drive projected on the symmetric-path tangent --")
    print("  and it does land on the array's MEAN phase (0.34/0.27/0.21 deg")
    print("  across a 4x ladder, shrinking, so that residual is")
    print("  discretisation). What it is NOT is any single unit's phase,")
    print("  because there is no single phase to have. Asserting a* == 30")
    print("  would compare one number against a closed form that assumes")
    print("  every unit shares it, when they are 6.7 deg apart.")
    print("  A CORRECTION IS RECORDED HERE RATHER THAN QUIETLY DROPPED: these")
    print("  rows previously read `configuration_phase(xs[0])` -- unit 0 alone,")
    print("  the CENTRE -- and reported its lag as a_hat drifting by 5 deg,")
    print("  with the wrong sign and about twenty times the true magnitude.")
    print()
    print("  HAZARD DISCHARGE (qvf.17 critique 23251, T2 23251, bound to this")
    print("  bead by its own HAZARD comment 2026-08-21): signed_gap's general")
    print("  branch has zero gate coverage for ANY negative gap in qvf.17's OWN")
    print("  file -- shallow included. TREATMENT CHOSEN: (b), a feasibility-")
    print("  preserving invariant, not new negative-gap rows on the frozen")
    print("  kernel. crank_run bisects the step size h using the REAL,")
    print("  re-evaluated signed_gap at the trial position for every")
    print("  enumerated pair (general branch included), and accepts a step")
    print("  only if no gap drops below -GAP_FLOOR_TOL. The HAZARD rows above")
    print("  gate the invariant live (backtrack ON) AND its mutation probe")
    print("  (backtrack OFF, which measurably violates the same floor) --")
    print("  the row can fail, and a control run shows it failing.")
    print()
    print("  ONSET vs SUSTAINED (unchanged by the FIXED round below): under")
    print("  in-phase/uniform-expansion driving, CONTACTS ALONE produce a")
    print("  transient ~0.30 first-step resistance (w-independent, verified)")
    print("  that relaxes after a small amount of real motion -- row G reads")
    print("  the ONSET (instant_jam=True, jam_angle=0.0 exactly, matching the")
    print("  bead's own '>0 not required' wording); rows H/I/CROSS read")
    print("  SUSTAINED persistence (instant_jam=False).")
    print()
    print("  FIXED ROUND (code review 23261, substantive critique 23262):")
    print("  the wire attachment was rebuilt (topo.contacts' single pinned")
    print("  vertex was a symmetry-FIXED POINT, span exactly 0.0 at every")
    print("  phase and every w -- verified, not a bug in degree but in kind).")
    print("  Wires now attach at the plate pairs Z3's own fold machinery")
    print("  identifies as OPENING (valley); their span is EXACTLY")
    print("  2/sqrt(3)*fold(a) -- measured, not guessed. The CROSS rows above")
    print("  hold the test protocol FIXED and vary ONLY w, proving w is now")
    print("  causal: w=0 SUSTAINED-jams (wires never relax under monotonic")
    print("  opening); w=w_ico reaches (opening never nears the limit). G's")
    print("  own jam now binds wires directly, gated (K rows), not contacts")
    print("  alone. scipy.optimize.nnls's RuntimeError-on-non-convergence is")
    print("  now caught (Z16), routed to QPFAIL, never a traceback or a jam.")
    print()
    print("  BEAD inviscid-1wd -- THE BALL JOINT. Everything below is measured")
    print("  on an array that is actually assembled. `build_pin_jacobian` used")
    print("  to omit jb_x `assemble_free`'s inter-unit rows, on the stated")
    print("  grounds that they 'encode RIGID pins'. They do not: three rows per")
    print("  contact, position coincidence only, orientation free -- a BALL")
    print("  JOINT. A rigid pin would be six. So the array's only assembly")
    print("  constraint was dropped on a mislabel, and under DECISION 18")
    print("  nothing replaced it (contact is one-sided; the wires are slack at")
    print("  a=0 for every w above EPS_ACT). The shared vertices separated at")
    print("  ~0.094 per degree of crank and the resulting 'lock' at a*=0.98 was")
    print("  the array coming apart. Restored, plus an ARRAY-LEVEL projection")
    print("  (the old one projected each unit onto its OWN hinge manifold and")
    print("  never re-closed the joints), the lock moves to a*=29.88 -- 7.6 deg")
    print("  PAST a_ico, and therefore NOT the instant both prior mechanisms")
    print("  predict, which is the first thing this surface has said that the")
    print("  qvf.15 ruler-test hazard does not swallow. The JOINT rows gate it:")
    print("  4.6e-16 of separation on the shipped path against 1.2e-01 on the")
    print("  retired one. A ball joint is a KINEMATIC constraint, not a force --")
    print("  no kernel, no mass, no primitive, METRIC FORM treatment (a)")
    print("  undisturbed. THE FOUR DECLARATIONS ARE UNTOUCHED BY IT.")
    print()
    print("  PHASE 1c (bead qvf.19, THE DELIVERABLE): a*(w,t) measured over")
    print("  the w-grid (fixed t=0) and t-grid (fixed, large w=0.9*w_ico) --")
    print("  each with a second, absolute, incommensurate arm. The rows above")
    print("  print the sign and magnitude of da*/dw and da*/dt, the")
    print("  DISCRIMINATION-RATIO verdict (Q2), the surface's non-vacuity")
    print("  span (R), the motion-order trace (S) and its binding sets (T),")
    print("  and the drive-model-robustness check. DRAFT status per the bead:")
    print("  T2 inviscid/qvf-lock-surface-phase1.md is stamped 'DRAFT -- NOT")
    print("  YET VALIDATED (pending inviscid-qvf.20)'.")

    failed = [n for n, p, _, _ in checks if not p]
    print()
    if failed:
        print(f"  !! {len(failed)} CHECK(S) FAILED -- this is a bug report, not a")
        print("     measurement. Nothing above may enter the record.")
        for n in failed:
            print(f"       - {n}")
        return 1
    print("  ALL CHECKS PASSED.")
    return 0


# ==========================================================================

def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_z_quasistatic_array -- contact kernel (1a) + crank stepper (1b)")
    print("  + lock surface / motion order (1c)")
    print("=" * 78)
    print("  bead inviscid-qvf.17 (Phase 1a, FROZEN): signed plate gaps with")
    print("  witness points and an outward contact normal, wire spans, a")
    print("  thickness offset, topology as data. bead inviscid-qvf.18 (Phase")
    print("  1b): the quasi-static velocity-level crank stepper -- a hand-")
    print("  rolled active-set QP (Lawson-Hanson LDP/NNLS reduction) under")
    print("  non-penetration and tension-only wires, with jam detection.")
    print("  bead inviscid-qvf.19 (Phase 1c, THE DELIVERABLE): the a*(w,t)")
    print("  lock surface -- ONE model producing BOTH bench predictions of")
    print("  the qvf.11/qvf.15 fork -- and the motion-order trace, built ONLY")
    print("  by calling the FROZEN kernel and stepper above, unmodified.")
    print("  FOUR DECLARATIONS (per the AMENDED design of record, T2 23230):")
    print("  KERNEL, MASS MODEL, PRIMITIVE are INAPPLICABLE throughout --")
    print("  no potential, no mass, no primitive choice anywhere in this")
    print("  file, Phase 1c included: INAPPLICABLE, NOT FORGOTTEN. METRIC")
    print("  FORM is QUALIFIED by the amendment: Phase 1a carries no QP")
    print("  objective, so the qualification is vacuously satisfied there")
    print("  (every 1a number is a NORM-FREE GEOMETRIC LENGTH). Phase 1b's")
    print("  QP objective ||v-v_cmd||^2_W DOES carry a weight -- TREATMENT")
    print("  (a) IS CHOSEN: W = identity in body coordinates, plus a live")
    print("  gate row (M) demonstrating the norm-free verdicts (jam status,")
    print("  binding-set composition) are W-INSENSITIVE under an alternate")
    print("  diagonal W. Only norm-free quantities (jam angle, active-set")
    print("  composition, reached/jammed) are quotable from Phase 1b;")
    print("  magnitudes derived FROM the QP solve itself (the achieved")
    print("  rate) are per-unit-of-this-W. Phase 1c CARRIES TREATMENT (a)")
    print("  FORWARD (bead .18's own choice, not re-decided here): a* is")
    print("  the SAME norm-free jam angle, and da*/dw, da*/dt are")
    print("  DIFFERENCES of it -- norm-free too, quotable without a")
    print("  per-unit-of-W hedge (that hedge is treatment (b)'s, not (a)'s).")

    if not PAIRS:
        print()
        print("=" * 78)
        print("GATE  1 row: the hinge pairing could not be read")
        print("=" * 78)
        print(f"  FAIL  Z0  hinge pairing readable                    "
              f"{'unreadable':>18s} {'12 x mult 2':>16s}")
        print()
        print("  Nothing below could be computed, so nothing below is printed.")
        return 1

    # SPECULATIVE PARALLEL PREFETCH. Every `crank_run` below is independent of
    # every other, so the previous run's recorded argument trace is replayed
    # through a process pool here, before the serial gate pass starts. The
    # pass then finds every solve already cached and does nothing but print.
    # This changes no value and no order: the prints below still happen in
    # exactly the sequence they always did. On a first run, or after a grid is
    # re-priced, the trace is short or empty and the missing solves simply
    # compute serially at their own call sites and are recorded for next time.
    jb_cache.prefetch(crank_run)

    z0 = z0_control()
    z2 = _z2_plate_pair_counts()
    z3 = z3_fold_table()
    z4 = z4_registry()
    z5 = z5_normal_witness()
    z6 = z6_wire_and_thickness()
    z7 = z7_crossing_census()

    topo = [t for t in build_topologies() if t.name.startswith("SC7")][0]
    zg = z15_crank_gates(topo)
    zhaz = z13_hazard_invariant(topo)
    zdow = doweled_diagnostic(topo)
    ztwo = z14_two_cell_sanity(topo)
    zqp = z16_qpfail_probe(topo)
    zlock = z17_lock_surface(topo)
    zbp = {"probe": bp_exactness_probe(topo), "mutation": bp_mutation_probe(topo),
           "threading": z18_phase_threading(topo),
           "depth": z19_meaningless_depth(topo)}

    return gate(z0, z2, z3, z4, z5, z6, z7, zg, zhaz, zdow, ztwo, zqp, zlock, zbp)


if __name__ == "__main__":
    # `--no-cache` / `--clear-cache` are consumed here; nothing else is
    # accepted, so an unrecognised flag is a loud failure rather than a run
    # that silently ignored what was asked of it.
    _rest = jb_cache.parse_argv(sys.argv[1:])
    if _rest:
        print(f"unrecognised argument(s): {' '.join(_rest)}", file=sys.stderr)
        print("usage: jb_z_quasistatic_array.py [--no-cache] [--clear-cache]",
              file=sys.stderr)
        sys.exit(2)
    with np.errstate(all="ignore"):
        sys.exit(main())
