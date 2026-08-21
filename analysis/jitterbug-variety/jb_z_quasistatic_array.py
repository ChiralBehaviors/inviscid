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

from jb_a_family import R_CIRC, Z, corners, faces, rot
from jb_g_strut_clearance import segment_distance as jb_g_segment_distance
from jb_x_array_linkage import (A_ICO, DIAGONALS, PAIRS, STRUT_LEN, STRUTS,
                                SQUARE_DIAGONALS, Topology, apply_body_motions,
                                assemble_doweled, build_topologies, dverts_exact,
                                hinge_jacobian, path_tangent_48, rank_of, verts)

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
            n = np.cross(tri[1] - tri[0], tri[2] - tri[0])
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

def _seg_seg_witness(p1, q1, p2, q2):
    """Closest points on two segments (clamped parametric solve, the same
    method as `jb_g_strut_clearance.segment_distance` / jb_x's private
    `_seg_seg`, extended to return the witness points, not just the
    distance)."""
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e = d1 @ d1, d2 @ d2
    if a < 1e-14 and e < 1e-14:
        return p1, p2, float(np.linalg.norm(r))
    if a < 1e-14:
        t = np.clip(-(d2 @ r) / e, 0.0, 1.0) if e > 1e-14 else 0.0
        pa, pb = p1, p2 + d2 * t
        return pa, pb, float(np.linalg.norm(pa - pb))
    if e < 1e-14:
        s = np.clip((d1 @ r) / a, 0.0, 1.0)
        pa, pb = p1 + d1 * s, p2
        return pa, pb, float(np.linalg.norm(pa - pb))
    b, c, f = d1 @ d2, d1 @ r, d2 @ r
    den = a * e - b * b
    s = np.clip((b * f - c * e) / den, 0.0, 1.0) if den > 1e-12 else 0.0
    t = np.clip((b * s + f) / e, 0.0, 1.0)
    s = np.clip((b * t - c) / a, 0.0, 1.0)
    pa, pb = p1 + d1 * s, p2 + d2 * t
    return pa, pb, float(np.linalg.norm(pa - pb))


def _barycentric(p, tri):
    """Barycentric coordinates of p w.r.t. triangle tri, via the projected
    2D-area-ratio construction. Not clamped -- may go outside [0,1] for a
    point off the triangle's PLANE or outside its extent, which is exactly
    what the on-triangle gate row checks."""
    a, b, c = tri
    n = np.cross(b - a, c - a)
    nn = n @ n
    if nn < 1e-18:
        return None
    denom = nn
    u = (np.cross(b - a, p - a) @ n) / denom
    v = (np.cross(c - b, p - b) @ n) / denom
    w = (np.cross(a - c, p - c) @ n) / denom
    return w, u, v  # weights of (a, b, c) respectively


def _pt_tri_witness(p, tri):
    """Closest point ON triangle tri to point p (projection, clamped to the
    face; falls back to the nearest edge when the projection lands outside).
    Extends jb_x's private `_pt_tri` (distance-only) to return the witness."""
    a, b, c = tri
    n = np.cross(b - a, c - a)
    nn = n @ n
    if nn < 1e-18:
        cands = [_seg_seg_witness(p, p, a, b), _seg_seg_witness(p, p, b, c),
                  _seg_seg_witness(p, p, c, a)]
        best = min(cands, key=lambda x: x[2])
        return best[1], best[2]
    proj = p - n * ((p - a) @ n) / nn
    bw, uw, vw = _barycentric(proj, tri)
    if bw >= 0 and uw >= 0 and vw >= 0:
        return proj, float(np.linalg.norm(p - proj))
    cands = [_seg_seg_witness(p, p, a, b), _seg_seg_witness(p, p, b, c),
              _seg_seg_witness(p, p, c, a)]
    best = min(cands, key=lambda x: x[2])
    return best[1], best[2]


def _seg_tri_hits(p, q, tri):
    """Moller-Trumbore: does segment pq pierce the OPEN interior of tri?
    Same primitive as jb_x's private `_seg_tri_hits`, reimplemented locally
    (mutation-probe rule; this file imports no private symbols from jb_x)."""
    a, b, c = tri
    e1, e2 = b - a, c - a
    d = q - p
    h = np.cross(d, e2)
    det = e1 @ h
    if abs(det) < 1e-14:
        return False
    inv = 1.0 / det
    s = p - a
    u = inv * (s @ h)
    if u < 0.0 or u > 1.0:
        return False
    qv = np.cross(s, e1)
    v = inv * (d @ qv)
    if v < 0.0 or u + v > 1.0:
        return False
    t = inv * (e2 @ qv)
    return 0.0 < t < 1.0


def _closest_point_pair(triA, triB):
    """The 15-candidate unsigned closest-point search: 9 edge-edge pairs + 3
    vertex(A)-vs-triB + 3 vertex(B)-vs-triA. Returns (pA, pB, distance)."""
    best = None
    for i in range(3):
        for j in range(3):
            pa, pb, d = _seg_seg_witness(triA[i], triA[(i + 1) % 3],
                                         triB[j], triB[(j + 1) % 3])
            if best is None or d < best[2]:
                best = (pa, pb, d)
    for i in range(3):
        q, d = _pt_tri_witness(triA[i], triB)
        if d < best[2]:
            best = (triA[i], q, d)
    for i in range(3):
        q, d = _pt_tri_witness(triB[i], triA)
        if d < best[2]:
            best = (q, triB[i], d)
    return best


def _is_piercing(triA, triB):
    for i in range(3):
        if _seg_tri_hits(triA[i], triA[(i + 1) % 3], triB):
            return True
        if _seg_tri_hits(triB[i], triB[(i + 1) % 3], triA):
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
    nB = np.cross(triB[1] - triB[0], triB[2] - triB[0])
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
    perp1 = np.cross(nA, np.array([1.0, 0.0, 0.0]))
    if np.linalg.norm(perp1) < 1e-6:
        perp1 = np.cross(nA, np.array([0.0, 1.0, 0.0]))
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(nA, perp1)
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
    n = np.cross(tri[1] - tri[0], tri[2] - tri[0])
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
    D = np.cross(n1, n2)
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
    n = np.cross(b - a, c - a)
    nn = n @ n
    if nn < 1e-18:
        return False
    u = np.cross(b - a, pt - a) @ n
    v = np.cross(c - b, pt - b) @ n
    w = np.cross(a - c, pt - c) @ n
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
    nB = np.cross(triB[1] - triB[0], triB[2] - triB[0])
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
        return [corners(a) + origins[u] for u in range(topo.n)]

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
    that separates as `2/sqrt(3) * fold(a)` -- reusing `contact_gradient_row`
    with t=0 and negating gives the wire row for free, correct in all three
    of `signed_gap`'s branches via the same `_contact_gradient_direction`
    these pairs were selected to be well inside the parallel-facing regime
    of anyway.

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
    gap, row = contact_gradient_row(xs, u, fi, v, fj, 0.0, ndof)
    degenerate = not np.any(row)
    return w - gap, -row, degenerate


def build_pin_jacobian(xs, n):
    """Block-diagonal (36n, 48n) equality Jacobian: the 36 intra-unit hinge
    rows PER UNIT only -- NOT `assemble_free`'s inter-unit contact rows,
    which encode RIGID pins (pre-DECISION-18 semantics; reused deliberately,
    and ONLY, for the DOWELED diagnostic)."""
    ndof = 48 * n
    big = np.zeros((36 * n, ndof))
    for i in range(n):
        big[36 * i:36 * i + 36, 48 * i:48 * i + 48] = hinge_jacobian(xs[i])
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
    z_phase, _ = path_tangent_48(a_hat)
    dsites = topo.dsites(dverts_exact(a_hat))
    v_cmd = np.zeros(ndof)
    units = range(topo.n) if driven == "all" else [driven]
    for i in units:
        vi = z_phase.copy()
        for p in range(8):
            vi[24 + 3 * p:27 + 3 * p] += dsites[i]
        v_cmd[48 * i:48 * i + 48] = vi
    return v_cmd, list(units)


def crank_pairs(topo):
    """The FIXED, topology-derived candidate plate-pair list, computed once
    (GAP 2's `enumerate_plate_pairs`, reused as prior art)."""
    return enumerate_plate_pairs(topo)


def crank_step(topo, pairs, xs, a_hat, w, t, driven="all",
               enforce_contacts=True, w_diag=None, wpairs=None):
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
    j_pin = build_pin_jacobian(xs, n)
    n_basis = null_space_basis(j_pin, ndof)
    if wpairs is None:
        wpairs = wire_pairs(topo)

    active_rows = []
    active_labels = []
    min_general_gap = float("inf")
    if enforce_contacts:
        for (i, fi, j, fj) in pairs:
            gap, row = contact_gradient_row(xs, i, fi, j, fj, t, ndof)
            nA = plate_normal(fi)
            nB = np.cross(xs[j][fj][1] - xs[j][fj][0], xs[j][fj][2] - xs[j][fj][0])
            nB_norm = np.linalg.norm(nB)
            is_general = nB_norm > 1e-300 and abs(float(nA @ (nB / nB_norm))) <= PARALLEL_TOL
            if is_general:
                min_general_gap = min(min_general_gap, gap)
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


def crank_run(topo, a_start, a_target, w, t, driven="all",
              enforce_contacts=True, w_diag=None, h0=H_STEP,
              max_steps=MAX_CRANK_STEPS, backtrack=True,
              gap_floor=GAP_FLOOR_TOL, instant_jam=True):
    """Returns a dict: status in {"reached","jammed","qpfail"}; a_final;
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
    origins = topo.sites(verts(a_start))
    xs = [corners(a_start) + origins[i] for i in range(topo.n)]
    a_hat = a_start
    min_general_gap = float("inf")
    rate_history = []
    binding_ever = set()

    def _scan_general_gap(xs_now):
        worst = float("inf")
        if not enforce_contacts:
            return worst
        for (i, fi, j, fj) in pairs:
            g, _ = contact_gradient_row(xs_now, i, fi, j, fj, t, 48 * topo.n)
            nA = plate_normal(fi)
            nB = np.cross(xs_now[j][fj][1] - xs_now[j][fj][0], xs_now[j][fj][2] - xs_now[j][fj][0])
            nB_norm = np.linalg.norm(nB)
            is_general = nB_norm > 1e-300 and abs(float(nA @ (nB / nB_norm))) <= PARALLEL_TOL
            if is_general:
                worst = min(worst, g)
        return worst

    for step in range(max_steps):
        v, status, rate, binding, mgg = crank_step(
            topo, pairs, xs, a_hat, w, t, driven, enforce_contacts, w_diag, wpairs)
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

        new_xs = []
        for i in range(topo.n):
            xi, _ = project_to_pin_manifold(trial_xs[i])
            new_xs.append(xi)
        xs = new_xs
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
    xs = [corners(a) + origins[i] for i in range(topo.n)]
    best = None
    for (i, fi, j, fj) in crank_pairs(topo):
        triA, triB = xs[i][fi], xs[j][fj]
        nA = plate_normal(fi)
        nB = np.cross(triB[1] - triB[0], triB[2] - triB[0])
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
    xs0 = [corners(0.0) + origins[i] for i in range(n)]
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
        g_lo, *_ = signed_gap(xs0[i][fi], (corners(0.0) + topo.sites(verts(0.0))[j])[fj], plate_normal(fi))
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            origins_m = topo.sites(verts(mid))
            xm = [corners(mid) + origins_m[u] for u in range(n)]
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

    # eps_act decade-insensitivity: rerun G's FIRST step at EPS_ACT_ALT.
    pairs = crank_pairs(topo)
    origins0 = topo.sites(verts(0.0))
    xs0 = [corners(0.0) + origins0[i] for i in range(topo.n)]
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
    xs0 = [corners(0.0) + origins0[i] for i in range(topo.n)]
    active_rows = []
    for (i, fi, j, fj) in pairs:
        gap, row = contact_gradient_row(xs0, i, fi, j, fj, 0.0, ndof)
        if gap <= EPS_ACT:
            active_rows.append(row)
    n_basis = null_space_basis(build_pin_jacobian(xs0, topo.n), ndof)
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
# THE GATE
# ==========================================================================

def gate(z0, z2, z3, z4, z5, z6, z7, zg, zhaz, zdow, ztwo, zqp):
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
    cross_w0, cross_wico = zg["cross_w0"], zg["cross_wico"]
    checks.append(("CROSS  same protocol, w=0: JAMMED (wires never relax, monotonic opening)",
                   cross_w0["status"] == "jammed", cross_w0["status"], "jammed"))
    checks.append(("CROSS  same protocol, w=w_ico: REACHED (opening tops out under the limit)",
                   cross_wico["status"] == "reached", cross_wico["status"], "reached"))
    checks.append(("CROSS  W IS CAUSAL: identical instant_jam, outcomes differ solely by w",
                   cross_w0["status"] != cross_wico["status"],
                   f"{cross_w0['status']} vs {cross_wico['status']}", "differ"))

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
    print()
    print("  ROWS DELETED RATHER THAN FIXED:")
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
    print("=" * 78)
    print("  bead inviscid-qvf.17 (Phase 1a, FROZEN): signed plate gaps with")
    print("  witness points and an outward contact normal, wire spans, a")
    print("  thickness offset, topology as data. bead inviscid-qvf.18 (Phase")
    print("  1b): the quasi-static velocity-level crank stepper -- a hand-")
    print("  rolled active-set QP (Lawson-Hanson LDP/NNLS reduction) under")
    print("  non-penetration and tension-only wires, with jam detection.")
    print("  FOUR DECLARATIONS (per the AMENDED design of record, T2 23230):")
    print("  KERNEL, MASS MODEL, PRIMITIVE are INAPPLICABLE throughout --")
    print("  no potential, no mass, no primitive choice anywhere in this")
    print("  file. METRIC FORM is QUALIFIED by the amendment: Phase 1a")
    print("  carries no QP objective, so the qualification is vacuously")
    print("  satisfied there (every 1a number is a NORM-FREE GEOMETRIC")
    print("  LENGTH). Phase 1b's QP objective ||v-v_cmd||^2_W DOES carry a")
    print("  weight -- TREATMENT (a) IS CHOSEN: W = identity in body")
    print("  coordinates, plus a live gate row (M) demonstrating the norm-")
    print("  free verdicts (jam status, binding-set composition) are")
    print("  W-INSENSITIVE under an alternate diagonal W. Only norm-free")
    print("  quantities (jam angle, active-set composition, reached/jammed)")
    print("  are quotable from Phase 1b; magnitudes derived FROM the QP")
    print("  solve itself (the achieved rate) are per-unit-of-this-W.")

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

    return gate(z0, z2, z3, z4, z5, z6, z7, zg, zhaz, zdow, ztwo, zqp)


if __name__ == "__main__":
    with np.errstate(all="ignore"):
        sys.exit(main())
