"""Step Y: DEPHASING. What happens to the array lock when the units are NOT in
phase with one another.

Bead `inviscid-qvf.14`. This is the first out-of-phase array model in the
project: `jb_x_array_linkage.py` built the assembled linkage and measured it
entirely IN PHASE, and the bead's whole content is the transverse direction that
file never entered.

THE OBSERVATION BEING MODELLED, and where it came from
------------------------------------------------------
The project owner reports that his physical array freezes ONLY if every unit
reaches the icosahedral phase SIMULTANEOUSLY; given a fractional difference
between units, motion progresses. His reading: simultaneity across every unit is
something an EXTERNAL COORDINATOR imposes. Eight cells by hand, standing outside
the system -- not 10^36 units inside a medium, which cannot agree on simultaneity
without exchanging information, which takes time, which is the thing being
modelled. Hence the reframing: the lock is not an obstacle to explain away, it is
evidence that the IN-PHASE MODE IS NOT PHYSICALLY REALISABLE AT SCALE, so
whatever motion the medium has must carry a PHASE GRADIENT.

THE QUESTION AS THE BEAD POSES IT
---------------------------------
Perturb the in-phase family TRANSVERSALLY and determine whether the binding
span's excess over its member is relieved at FIRST order (the lock is measure
zero) or only at SECOND (a threshold, hence a minimum wavevector, hence a
maximum wavelength).

WHAT THIS FILE ANSWERS, IN ONE PARAGRAPH
----------------------------------------
BOTH bead outcomes are refuted, and for the same reason: the relief question is
FIRST ORDER in every case measured, but its SIGN depends on which member binds
and on whether the cluster carries self-stress. There is no second-order regime
anywhere, so there is NO THRESHOLD, no minimum wavevector and no maximum
wavelength -- the excess is homogeneous of degree one in the dephasing amplitude
to seven decimals on two incommensurate ladders. The INTRA-unit binder (the six
folding diagonals, a = 22.238756093) is NOT relieved: the minimum over the
traceless sup-sphere of the one-sided directional derivative is strictly
POSITIVE and equals |s'(a_ico)|/(N-1) exactly, so the in-phase locus is a strict
transverse minimum with a CORNER (Danskin/Clarke), "the derivative is zero" is
FALSE and "the derivative does not exist" is TRUE. The INTER-unit binder (the
array-induced chord, a = 24.119490) is relieved at first order in CUBE8-M and
not in CUBE27-M, and the discriminator is the STATE OF SELF-STRESS that appears
once a cluster carries BOTH a cycle and a unit with four or more contacts (the
measured discriminator: cycles alone do not do it and coordination alone does not
either, and each has a witness in the run). The relieving pattern exists as a
phase field in both clusters and stops being KINEMATICALLY ADMISSIBLE in the
larger one. In the star, SQUARE4 and CUBE8-R that same chord is not held by the
linkage at all -- the assembly's own mechanisms relieve it with no dephasing
whatever -- which is a correction to jb_x, not a result about dephasing. Static
and path AGREE on every row, so the bead's proposed resolution is empty AT FIRST
ORDER, which is the only order two LPs over the same tangent cone could have
settled; what resolves the tension instead is that the lock is measure zero AS A
SET (Y4) and unrelieved AS A FUNCTIONAL (Y3) -- two different objects, neither
refuting the other. WHICH of them the owner's observation is about is a READING
and is labelled as one in Y4(b) and Y6(5), not a measurement. And the largest
scaling result here is one the first version of this file stopped one cube short
of: ADMISSIBLE DEPHASING SATURATES AT A FIXED DIMENSION independent of array
size (27, 64 and 125 units all give the same number), so the admissible FRACTION
goes to zero as the array grows.

THE MODEL, DECLARED
-------------------
DEPHASING IS INTRINSICALLY A DOWELED-MODEL CONCEPT and this file uses the doweled
model only. In the FREE model of jb_x a unit has 48 body variables and moves on a
six-dimensional internal variety, on which "the phase" is not a coordinate at
all: only the symmetric one-parameter family has a phase, and the DOWEL is what
restricts a unit to it. So a unit here is EXACTLY {rigid placement} x {phase}:
six placement variables and one phase, the same seven variables jb_x's doweled
array carries, and the nonlinear map is exact rather than linearised.

    unit i at phase a_i, placement (w_i, t_i):  x = R(w_i) v_k(a_i) + t_i
    contact (i,k,j,l):  R(w_i) v_k(a_i) + t_i - R(w_j) v_l(a_j) - t_j = 0

Nothing else is imposed. The lattice spacing is never fixed by hand; the unit
origins are an initial guess for the solver and an answer only where the solver
confirms them.

THE MECHANISM CORRECTION, which changes the answer and is not in jb_x
---------------------------------------------------------------------
A span of the assembled array IS NOT A FUNCTION OF THE PHASES. At fixed phases
the vertex-jointed assembly still has internal MECHANISMS -- ker of the placement
Jacobian, minus the six global rigid motions -- and a mechanism moves inter-unit
spans without moving any phase. jb_x's span enumeration evaluates every span at
the PURE-TRANSLATE reference placement, which is one point of that kernel, so its
inter-unit taut angles are placement-conditional. This file therefore optimises
over the mechanisms as well as over the dephasing, which is the only question a
physical array answers, and reports the mechanism-only relief separately as its
own row (Y3c). For the SIX-AROUND-ONE STAR that row is decisive: the star is a
tree, its mechanisms alone relieve the inter-unit chord at rate -5.79e-01 per
unit placement rate, and its inter-unit lock is an artefact of holding the
reference placement rather than a property of the linkage.

The intra-unit diagonals are immune to this correction by construction -- a
diagonal joins two vertices of ONE rigid unit, so its length is a function of
that unit's phase alone and no placement can touch it. That asymmetry is why the
two binders answer the question differently, and it is structural, not numerical.

THE TIED-ORBIT HAZARD, handled rather than avoided
--------------------------------------------------
At the in-phase locus the equivalent spans are TIED: the binding span is one
member of an orbit of equal-length spans, so a max-over-spans statistic is a
MIN-MAX at a symmetric point and naive differentiation is wrong. The recorded
rule (memo of record, RETIRED witness-selection hypothesis) is that smoothness
breaks when a PROPER SUBSET of a tied orbit is selected AND its members have
different derivatives. Both halves are MEASURED here (Y2): the intra-unit orbit
has 6N members and exactly N distinct gradients -- the six diagonals inside one
unit share a gradient to 0.0e+00, and gradients across units are orthogonal --
so the corner is produced by the tie ACROSS units and not by the tie within one.
Nothing here differentiates a max; every first-order statement is a LINEAR
PROGRAM over the one-sided cone, which is exact for a max of linear functions,
and the directional sampling is a confirmation of the LP rather than the method.

FOUR DECLARATIONS: INAPPLICABLE, NOT FORGOTTEN
-----------------------------------------------
The epic's four standing declarations -- interaction KERNEL, MASS MODEL,
PRIMITIVE (vertex or strut-midpoint), METRIC FORM (section or momentum-free) --
are DYNAMICAL, and every quantity in this file is KINEMATIC: ranks and singular
values of constraint Jacobians, span lengths and their derivatives with respect
to angle, one-sided directional derivatives of a max-of-lengths, dimensions of
linear subspaces, and linear-programme optima in those same units. No energy, no
mass, no time, no frequency and no DYNAMICAL metric -- no kinetic form, no
section, no momentum map -- appears anywhere, and the one place a reader might
expect one, the "wavevector" of Y5, is a lattice phase increment per lattice
step, a pure number, with no dispersion relation and no speed attached to it.
KERNEL, MASS MODEL and PRIMITIVE are therefore INAPPLICABLE outright.

METRIC FORM IS INAPPLICABLE IN ITS DYNAMICAL SENSE AND NOT ABSENT AS A
NORMALISATION, and the difference matters enough to state rather than let a
reader discover. Every first-order MAGNITUDE reported here is a rate per unit of
a chosen norm: for D+ that norm is the SUP-NORM ON THE PHASE RATES IN DEGREES,
and in `mechanism_only_relief` it is the sup-norm on a MIXED rotation/translation
placement rate whose components do not share a unit. So the three mechanism-relief
magnitudes are NOT comparable with one another in any unit, and the 1/(N-1)
flattening of Y3(b) is a statement about the SUP-NORM BUDGET and about no other.
What the file's verdicts actually rest on is SIGNS -- whether a traceless u with
D+(u) < 0 exists is invariant under any positive-homogeneous renormalisation of
u -- and every verdict drawn anywhere below is a sign or a dimension, never a
magnitude. Magnitudes are printed because they are what was measured; the reader
is told here, once, that they carry a normalisation and not a physical scale.
This paragraph exists because a critic on jb_x found that the rule can lapse
silently, and a later reader must be able to tell "checked, does not bite" from
"overlooked" -- and because a critic on THIS file found the previous version of
this paragraph claiming, wrongly, that no metric appeared at all.

DECISION 16 governs interference: struts may pass through one another. This file
attaches NO admissibility verdict to any overlap and does not scan for one; the
one place a clearance question could arise (a dephased unit's plates against its
neighbour's) is declared OUT OF SCOPE below rather than measured and left
unlabelled.

TOPOLOGY, STATED RATHER THAN ASSUMED
------------------------------------
Measured for the topologies jb_x already carries: SINGLE-VERTEX contacts, on the
M basis (generators 0,1,3, the one with a diagonal generator pair) and the R
basis (0,1,2, which has none). The IVM honeycomb's face-sharing alternative is a
DIFFERENT constraint set and is NOT measured here; jb_x records that a
face-bonded array has exactly one admissible configuration, which would make the
dephasing question empty, and settling that is not this bead. Topology
sensitivity is reported where it is cheap: the two bases give different
admissible-dephasing dimensions at the same cluster size (Y1c), so the answer IS
topology-sensitive and the M basis is the more constrained of the two.

WHAT THIS FILE DOES NOT SETTLE
------------------------------
* Which members the owner's array actually carries. That is a BUILD FACT, not a
  geometric one, and jb_x already records it as unresolved. This file therefore
  measures BOTH candidate binders and reports that they answer differently,
  rather than picking one.
* Anything beyond first order in the dephasing amplitude for the SIGN of the
  relief. Second order is measured only to establish that it is not needed
  (the first-order term never vanishes), not to characterise it.
* Any dynamics. There is no equation of motion here, so "propagation" means
  a one-parameter family of admissible configurations, not a wave with a speed.
* Clearance between dephased neighbours. Out of scope, and unmeasured.
* Whether the mechanisms this file optimises over are present in the owner's
  rig. A physical wire may be a tension-only member that removes some of them.
* WHETHER THE DOWEL ITSELF IS THE OWNER'S RIG. This is the deferred P1 that
  jb_x's record carries and that this file INHERITS: the doweled model here
  drives a unit along the TRUE PATH TANGENT of the symmetric family, while the
  owner's physical rig rides the shared-vertex ELLIPSE. In the FREE model a unit
  moves on a six-dimensional internal variety on which the phase is not a
  coordinate, so an intra-unit diagonal's length is a function of the whole
  internal state and NOT of a phase -- and a non-symmetric internal deformation
  could then relieve that diagonal WITH NO DEPHASING AT ALL. Every statement
  below about the intra-unit wall surviving is therefore DOWEL-CONDITIONAL, and
  that is the single most consequential unmodelled thing here, ahead of the
  tension-only-member question above. Settling it needs the free model, which is
  jb_x's object and not this bead's.
* THE BULK (Y1e) IS AT FIXED LATTICE PERIOD. The Bloch ansatz carries a rigid
  cell translation and no homogeneous lattice strain, so it describes strict
  plane waves only. The uniform (k = 0) dephasing IS executable on every finite
  cluster measured -- by an exactly affine translation field with a NON-SCALAR
  linear part, i.e. by straining the lattice -- and that motion is outside the
  ansatz by construction, not obstructed by it. Both halves are measured in
  Y1(e); nothing here is inferred from the ansatz alone.

CONVENTIONS INHERITED FROM THIS DIRECTORY
-----------------------------------------
Deterministic and byte-identical across runs; exit code from the gate table; no
raise inside a swept loop; a check whose non-vacuity is printed prose rather than
an assertion cannot fail; a guard band is constrained from ABOVE as well as
below; every sweep grid has a SECOND, ABSOLUTE, INCOMMENSURATE arm because a
commensurate grid misses what lands on it -- three recorded instances in this
directory, one of them in the file this one extends. Run from the repository root
with python3.
"""
import sys

import numpy as np

import jb_cache

#: Importable name of THIS module -- a literal, because `__name__` is
#: "__main__" under direct execution and a prefetch worker in a fresh
#: interpreter must be able to re-import it by name.
_MODULE = "jb_y_dephasing"
from scipy.optimize import brentq, linprog

from jb_x_array_linkage import (A_ICO, DIAGONALS, STRUT_LEN, Topology, _classes,
                                _class_positions, _hat, build_topologies,
                                diagonal_generator_pairs, dverts_exact, verts,
                                PAIRS)
from jb_a_family import rot

# ==========================================================================
# LOCAL CONSTANTS
#
# Every constant a mutation probe needs to reach is defined HERE, locally, and
# nothing below reads a threshold from another module. jb_v's post-mortem
# records a probe that mutated a name the file never defined locally, so the
# mutation never applied and a clean exit was read as confirmation.
# ==========================================================================

#: The intra-unit lock angle: where the six folding square diagonals reach strut
#: length. Imported from jb_x as A_ICO and RE-DERIVED here by root-finding on
#: this file's own span function, so it is a cross-check rather than an input.
A_ICO_RECORD = A_ICO

#: The array-induced inter-unit lock angle recorded in T2 `qvf.11-array-linkage`
#: (RETRACTION R1): the shortest INSTALLABLE inter-unit chord reaches strut
#: length here, ABOVE the icosahedral phase, so in an assembled array carrying
#: that chord it is this angle and not A_ICO that binds. A number from OUTSIDE
#: this file; re-derived below and compared inside a two-sided band.
A_CHORD_RECORD = 24.119490

#: R1's span/strut ratio for that chord at the icosahedral phase.
CHORD_RATIO_RECORD = 1.0705

#: R1's counts of inter-unit spans taut at A_CHORD_RECORD, per cluster. Numbers
#: from outside this file, re-derived and compared.
CHORD_COUNT_RECORD = {"SC7 star (six-around-one)": 8,
                      "CUBE8-M (60/60/90 basis)": 8,
                      "CUBE8-R (60/60/60 basis)": 6,
                      "CUBE27-M": 48}

#: Rank tolerance, relative to the largest singular value. Same convention and
#: same value as jb_x, defined locally so a probe can reach it.
RANK_RTOL = 1e-10

#: A span is ACTIVE (taut) at the lock angle when its length is within this of
#: the member length. Bounded from ABOVE and BELOW by measurement in Y0d: the
#: largest active |excess| must be under it and the smallest INACTIVE |excess|
#: must be over it, and both ends are gated.
ACTIVE_TOL = 1e-7

#: Newton settle threshold for the exact nonlinear dephased solve.
SOLVE_TOL = 1e-13

#: The truncation ladder for the Gauss-Newton step's least-squares solve.
#: A SELF-STRESSED cluster has a rank-deficient contact Jacobian, and near the
#: in-phase locus the deficiency is APPROXIMATE rather than exact: the singular
#: values that vanish at zero dephasing sit at O(eps) instead. `lstsq` with the
#: default cutoff keeps them, and the resulting step is those directions'
#: numerical noise amplified by their reciprocal -- measured on CUBE27-M at
#: eps = 1e-2 along an admissible affine ramp, |dz| came out at 3.6e+01 for a
#: residual of 1.8e-03, no line-search step reduced anything, and the solve
#: returned its own seed. Escalating the cutoff until a step is accepted
#: truncates exactly those directions; the first entry is the untruncated solve,
#: so a well-conditioned cluster takes the same step it always did and its
#: numbers are unchanged. THIS DEFECT WAS INVISIBLE UNTIL THE SOLVE RESIDUAL WAS
#: GATED: `_amplitude_ladder` computed the worst residual and no caller read it.
LSTSQ_RCONDS = (None, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)

#: Consecutive non-improving Gauss-Newton iterations after which the solve stops
#: and returns its residual for the gate to read. A configuration the solver
#: cannot close must cost bounded time and arrive as a NUMBER, not as a hang.
#: This one is a COST guard and not a correctness guard: no gate row distinguishes
#: it from a larger value, and the record says so rather than implying otherwise.
SOLVE_STALL = 3

#: The cluster, direction and amplitude at which the truncated step is DECISIVE,
#: so that LSTSQ_RCONDS is exercised by a row rather than merely believed. A
#: (2,2,4) box carries self-stress; along an affine ramp at this amplitude the
#: untruncated least-squares step stalls two decades above the seed's residual
#: while the truncated one closes to 4e-14. Y0(g) asserts BOTH halves -- a fix
#: whose absence changes no row is a fix nothing tests, which is the shape this
#: directory keeps finding.
RCOND_BOX = (2, 2, 4)
RCOND_EPS = 0.1

#: The amplitude ladder for the order fit. GEOMETRIC WITH AN IRRATIONAL RATIO,
#: so no two rungs are commensurate with one another and none can land on a
#: feature of the function by construction. Golden ratio, so the ladder is as
#: far from any rational ratio as a ladder can be.
EPS_TOP = 1.0
EPS_RATIO = 2.0 / (1.0 + np.sqrt(5.0))       # 1/phi = 0.6180339887...
EPS_RUNGS = 40

#: The SECOND amplitude ladder, absolute and incommensurate with the first, used
#: as the step-independence arm of the order fit. Written as a literal and NOT
#: as a multiple of EPS_RATIO: an earlier file in this directory wrote its second
#: arm as a RATIO of the first, so coarsening one coarsened both in lock-step and
#: the row was structurally unable to see a change.
EPS_TOP_ALT = 0.7071067811865476
EPS_RATIO_ALT = 0.5436890126920763           # 1/(1+1/phi+1/phi^2+...) style; irrational

#: The fit window for the order, bounded at BOTH ends by absolute numbers. Below
#: LADDER_FLOOR the excess is at the level of the Gauss-Newton solve's own
#: residual; above LADDER_CEIL the nonlinear terms are visible and a first-order
#: fit would be measuring something it does not claim. A window open at either
#: end is a fit that cannot fail: open at the bottom it fits solver noise, open
#: at the top it reports curvature as order.
LADDER_FLOOR = 1e-9
LADDER_CEIL = 1e-6

#: Number of seeded quasi-random transverse directions in the corner spread.
#: The spread CONFIRMS the linear programme; it never replaces it.
N_SPREAD = 32
SPREAD_SEED = 20260819

#: Wavevector grids for Y5. Both are INCOMMENSURATE with the lattice: the offset
#: is irrational so no sample lands on a zone centre, a zone boundary or any
#: rational fraction of them. The second grid exists because a threshold sitting
#: on a grid point is invisible to a single grid and a FINER grid of the same
#: family does not help -- only a differently offset one does.
K_GRID = 37
K_OFFSET = 0.1234567891011121
K_GRID_ALT = 23
K_OFFSET_ALT = 0.3819660112501051            # 1/phi^2, irrational, unrelated to the above

#: Boxes used for the admissible-dephasing scaling table. Chosen to bracket the
#: transition: chains of three lengths (no cycle at all), the 2-thick boxes with
#: cycles but no highly-coordinated unit, the first boxes carrying both, and
#: then THREE CUBES OF INCREASING SIZE, because the dimension SATURATES and a
#: table stopping at 3x3x3 cannot see that it does.
SCALING_BOXES = ((2,), (5,), (9,), (2, 2), (2, 4), (4, 2), (3, 3), (4, 4),
                 (2, 2, 2), (2, 2, 3), (2, 2, 4), (2, 3, 3),
                 (3, 3, 3), (4, 4, 4), (5, 5, 5))

#: The cubes whose admissible-dephasing dimension is compared directly, to test
#: whether it grows with N or saturates. Must be at least three sizes: two
#: agreeing numbers are a coincidence, three are a plateau.
SATURATION_BOXES = ((3, 3, 3), (4, 4, 4), (5, 5, 5))

#: Minimum coordination -- number of contacts at one unit -- above which a unit
#: is called HIGHLY COORDINATED in the Y1(b) discriminator. Four, because three
#: is what a corner of a three-generator box carries and four is what its first
#: non-corner carries. The discriminator is MEASURED against the whole table,
#: with a witness on each side, rather than asserted from this number.
COORD_INTERIOR = 4

#: Below this, the Bloch phase column is IDENTICALLY ZERO and the dephasing
#: amplitude is unconstrained for the trivial reason. Bounded from both sides in
#: the gate: the zone corner's ||b|| must fall under it and the nearest control
#: point's must stay over it, so a floor raised until it swallowed real
#: obstructions would redden the second row.
BLOCH_B_FLOOR = 1e-12

#: Angles at which the admissible-dephasing dimension is re-measured, to test
#: whether it is a property of the topology or of the phase. Deliberately
#: spread and deliberately not a uniform grid.
ADM_ANGLES = (5.0, A_ICO_RECORD, 24.119490, 40.0, 55.3)

TOL = {
    "solve": SOLVE_TOL,
    "inphase": 1e-13,        # dephased solver at zero dephasing vs jb_x's reference
    "fd_jacobian": 1e-7,     # analytic (Jr,Jp) vs central differences of the exact residual
    "achord": 1e-6,          # re-derived chord angle vs the recorded value
    "aico": 1e-8,            # re-derived icosahedral angle vs the recorded value
    "tie": 1e-12,            # equality of tied gradients within one unit
    "homog": 1e-5,           # relative spread of |E|/eps down the amplitude ladder
    "lp": 1e-9,              # sign threshold for a linear-programme optimum
}
# An "adm" tolerance was defined here and read by nothing. A constant a
# mutation probe can reach but no row consults is worse than absent: it looks
# like a guard. The admissibility residual is instead bounded by
# BLOCH_B_FLOOR, which two gate rows read from both sides.

#: The DELIBERATE OFFSET that bounds TOL["achord"] FROM ABOVE. The re-derived
#: chord angle must agree with A_CHORD_RECORD to within TOL["achord"] AND must
#: DISAGREE, at that same tolerance, with a value offset by this much. Without
#: the second half the tolerance is unbounded above and the headline angle can
#: drift with the whole gate green -- the worst hole an independent validation
#: found in jb_x, reproduced here as a control rather than re-learned.
ACHORD_CONTROL_OFFSET = 1e-3

#: The quantisation of A_CHORD_RECORD: six decimal places, so half a unit in the
#: last place. A tolerance under this would be testing the record's rounding.
ACHORD_RECORD_QUANTUM = 5e-7

#: The same pair for A_ICO, whose record carries nine decimals.
AICO_CONTROL_OFFSET = 1e-3
AICO_RECORD_QUANTUM = 5e-10

#: The two Y0 FOUNDATION tolerances were bounded from BELOW by their own
#: measured deviations and from ABOVE by nothing at all: an independent
#: validation loosened each to 1.0 and the whole gate stayed green. These are
#: their control offsets, and each gets the same TWO-ROW idiom the chord and the
#: icosahedral angle already carry -- the measurement must PASS at the tolerance
#: and a deliberately offset reference must FAIL at that same tolerance.
#: INPHASE: the reference class positions displaced by this much in one
#: coordinate. FD: the analytic Jacobian's prediction scaled by (1 + this).
INPHASE_CONTROL_OFFSET = 1e-9
FD_CONTROL_OFFSET = 1e-3

#: The margin by which the SELF-STRESSED clusters' measured min D+ must EXCEED
#: the no-self-stress derivation |s'|/(N-1). The previous form was
#: `v > p - 1e-12`, which ACCEPTS v == p -- exactly the enumerator the prose
#: says the row catches, and a probe replacing that arm with True reddened
#: nothing. Bounded from ABOVE in the gate by the measured separation, so it
#: cannot be raised until it excludes the truth either.
DERIVED_MARGIN = 1e-3

#: The finite dephasing excursion of Y3(g). The bead's Outcome 2 is a FINITE
#: phase difference, and a ladder confined to |E| in [1e-9, 1e-6] -- 1e-7 to
#: 1e-4 DEGREES -- cannot see a turning point at a degree. These rungs run to
#: FOUR DEGREES on the exact nonlinear model, and every rung's closure residual
#: is gated: an excess read off a configuration that did not close is not a
#: measurement, which is the hole the discarded `res_max` was hiding.
EXCURSION_RUNGS = (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0)

#: The amplitude at which a KNOWN-NON-CLOSING direction is solved, so that the
#: closure gate above has a demonstrated failure mode rather than a claimed one.
EXCURSION_CTRL_EPS = 0.05

#: Floor for the population-gap band of ACTIVE_TOL. DERIVED, not picked: an
#: excess is the difference of two lengths of order STRUT_LEN, so one unit in
#: the last place of such a length is where the population stops carrying
#: information and a handful of spans sit at EXACTLY zero. Without a floor those
#: exact zeros open an unbounded log gap at the bottom and the widest-gap search
#: finds it instead of the real one. The gate asserts that the band's lower end
#: is STRICTLY ABOVE this floor and that the gap it brackets is many decades
#: wide, so the floor is a discard rule and not the answer.
EXCESS_FLOOR = float(np.finfo(float).eps) * STRUT_LEN

#: The band bracketing ACTIVE_TOL must be at least this many decades wide. A
#: population gap of half a decade is not a separation between "taut" and "not
#: taut"; the measured one here is over fourteen.
EXCESS_GAP_DECADES = 6.0


# ==========================================================================
# THE DEPHASED ARRAY
# ==========================================================================

def _rodrigues(w):
    """Exact rotation from a rotation vector (radians), no small-angle step."""
    th = float(np.linalg.norm(w))
    return np.eye(3) if th < 1e-14 else rot(w / th, np.degrees(th))


def solve_dephased(topo, phases, maxit=400, tol=SOLVE_TOL, rconds=None):
    """Exact placements for PRESCRIBED per-unit phases.

    Gauss-Newton on the exact nonlinear contact residual, seeded at the
    pure-translate placement built from the MEAN phase. Returns
    (z, residual_norm, iterations); never raises, so a topology that cannot
    close reaches the gate as a number rather than as a traceback.

    THE STEP IS TRUNCATED WHEN IT HAS TO BE. On a self-stressed cluster the
    Jacobian's rank deficiency is exact at zero dephasing and APPROXIMATE just
    off it, so the vanishing singular values sit at O(eps) and the untruncated
    least-squares step is their noise divided by them. See LSTSQ_RCONDS: the
    ladder starts at the untruncated solve, so a well-conditioned cluster takes
    exactly the step it took before and reproduces its previous numbers to the
    bit, and only a solve that would otherwise stall pays for the escalation.
    A solve that stalls anyway stops after SOLVE_STALL non-improving iterations
    and returns its residual, which every caller now gates.
    """
    n = topo.n
    V = [verts(a) for a in phases]
    seed = topo.sites(verts(float(np.mean(phases))))
    z = np.zeros(6 * n)
    z[3 * n:] = np.asarray(seed, float).reshape(-1)
    cs = topo.contacts

    def resid(zz):
        w = zz[:3 * n].reshape(n, 3)
        t = zz[3 * n:].reshape(n, 3)
        if not cs:
            return np.zeros(0)
        return np.concatenate([
            _rodrigues(w[i]) @ V[i][k] + t[i] - _rodrigues(w[j]) @ V[j][l] - t[j]
            for (i, k, j, l) in cs])

    def jac(zz):
        w = zz[:3 * n].reshape(n, 3)
        j = np.zeros((3 * len(cs), 6 * n))
        for m, (i, k, l_j, l) in enumerate(cs):
            pi = _rodrigues(w[i]) @ V[i][k]
            pj = _rodrigues(w[l_j]) @ V[l_j][l]
            j[3 * m:3 * m + 3, 3 * i:3 * i + 3] = -_hat(pi)
            j[3 * m:3 * m + 3, 3 * l_j:3 * l_j + 3] = +_hat(pj)
            j[3 * m:3 * m + 3, 3 * n + 3 * i:3 * n + 3 * i + 3] = np.eye(3)
            j[3 * m:3 * m + 3, 3 * n + 3 * l_j:3 * n + 3 * l_j + 3] = -np.eye(3)
        return j

    r = resid(z)
    if r.size == 0:
        return z, 0.0, 0
    it = 0
    stall = 0
    for it in range(1, maxit + 1):
        nr = float(np.linalg.norm(r))
        if not np.isfinite(nr) or nr < tol:
            break
        j = jac(z)
        moved = False
        for rcond in (LSTSQ_RCONDS if rconds is None else rconds):
            dz, *_ = np.linalg.lstsq(j, -r, rcond=rcond)
            for lam in (1.0, 0.5, 0.25, 0.1, 0.03):
                z2 = z + lam * dz
                r2 = resid(z2)
                n2 = float(np.linalg.norm(r2))
                if np.isfinite(n2) and n2 < nr:
                    z, r, moved = z2, r2, True
                    break
            if moved:
                break
        if not moved:
            break
        stall = stall + 1 if float(np.linalg.norm(r)) > 0.999 * nr else 0
        if stall >= SOLVE_STALL:
            break
    return z, float(np.linalg.norm(r)), it


def class_points(topo, members, phases, z):
    """Position of every wired point class in a dephased configuration."""
    n = topo.n
    w = z[:3 * n].reshape(n, 3)
    t = z[3 * n:].reshape(n, 3)
    V = [verts(a) for a in phases]
    return np.array([_rodrigues(w[mem[0][0]]) @ V[mem[0][0]][mem[0][1]]
                     + t[mem[0][0]] for mem in members])


def blocks(topo, a):
    """(Jr, Jp): the linearised contact Jacobian in PLACEMENT and in PHASE.

    Rows are the 3 components of each contact. Jr's columns are the 3N rotation
    rates then the 3N translation rates; Jp's are the N phase rates, in DEGREES,
    to match `dverts_exact`. Derived from the exact residual above, and checked
    against central differences of that residual in Y0b.
    """
    v = verts(a)
    dv = dverts_exact(a)
    n, cs = topo.n, topo.contacts
    jr = np.zeros((3 * len(cs), 6 * n))
    jp = np.zeros((3 * len(cs), n))
    for m, (i, k, j, l) in enumerate(cs):
        jr[3 * m:3 * m + 3, 3 * i:3 * i + 3] = -_hat(v[k])
        jr[3 * m:3 * m + 3, 3 * j:3 * j + 3] = +_hat(v[l])
        jr[3 * m:3 * m + 3, 3 * n + 3 * i:3 * n + 3 * i + 3] = np.eye(3)
        jr[3 * m:3 * m + 3, 3 * n + 3 * j:3 * n + 3 * j + 3] = -np.eye(3)
        jp[3 * m:3 * m + 3, i] += dv[k]
        jp[3 * m:3 * m + 3, j] -= dv[l]
    return jr, jp


def rank_of(m):
    """(rank, singular values). RANK is reported, never a subtraction."""
    if m.size == 0:
        return 0, np.zeros(1)
    s = np.linalg.svd(m, compute_uv=False)
    if not np.isfinite(s).all() or s[0] <= 0.0:
        return 0, s
    return int((s > s[0] * RANK_RTOL).sum()), s


def admissible_dephasing(topo, a):
    """The subspace of phase perturbations the linkage can actually execute.

    A dephasing rate `u` is ADMISSIBLE when some placement rate solves
    Jr z + Jp u = 0 exactly, that is, when Jp u lies in the range of Jr. The
    obstruction is therefore the STATES OF SELF-STRESS of the placement problem
    -- the left null space of Jr -- acting on the phase columns. Returns a dict
    with the counts and an orthonormal basis of the admissible subspace.
    """
    jr, jp = blocks(topo, a)
    rows = jr.shape[0]
    if rows == 0:
        return dict(rows=0, rank=0, selfstress=0, obstruct=0,
                    dim=topo.n, basis=np.eye(topo.n))
    rr, _ = rank_of(jr)
    u, _, _ = np.linalg.svd(jr, full_matrices=True)
    m = u[:, rr:].T @ jp
    if m.size == 0:
        return dict(rows=rows, rank=rr, selfstress=rows - rr, obstruct=0,
                    dim=topo.n, basis=np.eye(topo.n))
    ra, _ = rank_of(m)
    _, _, vt = np.linalg.svd(m)
    return dict(rows=rows, rank=rr, selfstress=rows - rr, obstruct=ra,
                dim=topo.n - ra, basis=vt[ra:].T)


def contact_graph(topo):
    """(cycle rank, max coordination, number of units at or above it).

    The two candidate discriminators for the dephasing obstruction, computed on
    the UNIT ADJACENCY GRAPH alone -- no geometry, no phase, no Jacobian -- so
    that Y1(b) can test each against the measured self-stress instead of
    asserting one. The cycle rank is the first Betti number E - N + C over the
    connected components; the coordination of a unit is its number of contacts.
    """
    n = topo.n
    deg = np.zeros(n, dtype=int)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (i, _, j, _) in topo.contacts:
        deg[i] += 1
        deg[j] += 1
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    comps = len({find(x) for x in range(n)})
    betti = len(topo.contacts) - n + comps
    mx = int(deg.max()) if n else 0
    return betti, mx, int((deg >= COORD_INTERIOR).sum())


def _population_gap(values, floor=EXCESS_FLOOR):
    """(below, above): the two values bracketing the WIDEST log gap in a
    population, and nothing else.

    The same idiom `_spectral_gap` applies to a singular-value spectrum, applied
    here to a population of span EXCESSES. It exists because the previous
    ACTIVE_TOL band was bracketed by "largest excess KEPT" and "smallest excess
    DISCARDED" -- both computed WITH ACTIVE_TOL -- so the band moved with the
    constant it was guarding, and an independent validation measured the row
    passing over more than three hundred decades while the measurement it
    certifies breaks fifteen decades sooner. That is verbatim the RANK_RTOL
    defect this file records fixing, left unfixed one constant along. The
    population here is filtered by NO tolerance at all, so the interval cannot
    move when ACTIVE_TOL moves.
    """
    v = np.sort(np.asarray([x for x in values if np.isfinite(x)], float))
    if v.size < 2:
        return 0.0, float("inf")
    lv = np.log10(np.clip(v, floor, None))
    d = lv[1:] - lv[:-1]
    i = int(np.argmax(d))
    return float(max(v[i], floor)), float(v[i + 1])


# ==========================================================================
# SPANS AND THEIR ACTIVE SETS
# ==========================================================================

def span_length(topo, members, a, p, q):
    """One span at the in-phase, pure-translate reference configuration."""
    pts = _class_positions(a, topo, members)
    return float(np.linalg.norm(pts[p] - pts[q]))


def active_set(topo, members, a, kind, tol=ACTIVE_TOL):
    """The MEMBERS of a candidate build, as ((unit_p, vertex_p),
    (unit_q, vertex_q)) pairs, together with the two numbers that bound
    ACTIVE_TOL from below and above.

    THE TWO BAND NUMBERS ARE COMPUTED WITH NO TOLERANCE. They are the pair
    bracketing the widest log gap in the population of candidate excesses --
    every decreasing intra-unit pair for "intra", every decreasing disjoint span
    for "inter" -- so they are a property of the geometry at this angle and
    nothing else. An earlier version returned "largest excess kept" and
    "smallest excess discarded", both computed WITH `tol`, which made the band
    move with the constant it was guarding: see `_population_gap`.

    A FAITHFUL REVERT OF THAT EARLIER VERSION MUST TOUCH BOTH RETURN SITES --
    the "intra" branch's `return tuple(pairs), worst, gap` above and the
    "inter" branch's below carried the SAME tol-dependent rule, and a revert
    that only restores one of them is not a demonstration of the fix. Verified
    directly (both sites reverted together, swept from tol=1e-300 to tol=1.0
    at this file's own re-derived angles): the pre-fix band check was True
    over the ENTIRE range tested, not merely late by some number of decades --
    it never reddened. A single-branch revert can look like a demonstration
    (it does redden across part of that range) without being a faithful one.

    A member is written in UNIT-AND-VERTEX form rather than as a pair of wired
    point classes, because a class can contain several units and `members[c][0]`
    then picks an arbitrary one. At the in-phase configuration every
    representative of a class sits at the same point, so the LENGTH does not
    care -- but the PHASE ATTRIBUTION does, and attributing a span's rate to the
    wrong unit is exactly the error that would make the tied-orbit analysis
    meaningless. An earlier version of this function returned class pairs and
    reported 129 "intra-unit" members on an eight-unit cluster that has 48.

    kind == "intra": built STRUCTURALLY from `DIAGONALS`, six per unit, which is
        jb_x's independently derived set -- not filtered out of a crossing
        search. Each is verified to sit at strut length and to be decreasing,
        and the verification is what bounds the tolerance.
    kind == "inter": class pairs whose unit sets are DISJOINT, so no
        representative choice can turn one into an intra-unit span.
    """
    _, units = _classes(topo)
    if kind == "intra":
        v0 = verts(a)
        v1 = verts(a + 1e-2)
        pairs, pop = [], []
        for (k, l) in DIAGONALS:
            d0 = float(np.linalg.norm(v0[k] - v0[l]))
            d1 = float(np.linalg.norm(v1[k] - v1[l]))
            if not (d1 < d0):
                continue
            e = abs(d0 - STRUT_LEN)
            pop.append(e)
            if e < tol:
                for i in range(topo.n):
                    pairs.append(((i, k), (i, l)))
        # The rest of the population: every OTHER intra-unit pair, in or out of
        # the active set, tolerance-free. Exact duplicates of a vertex pair
        # (length identically the strut) are not candidates and are dropped.
        for k in range(12):
            for l in range(k + 1, 12):
                if (k, l) in DIAGONALS:
                    continue
                e = abs(float(np.linalg.norm(v0[k] - v0[l])) - STRUT_LEN)
                if e > 1e-14:
                    pop.append(e)
        worst, gap = _population_gap(pop)
        return tuple(pairs), worst, gap

    p0 = _class_positions(a, topo, members)
    p1 = _class_positions(a + 1e-2, topo, members)
    nc = len(members)
    iu = np.triu_indices(nc, 1)
    d0 = np.linalg.norm(p0[iu[0]] - p0[iu[1]], axis=-1)
    d1 = np.linalg.norm(p1[iu[0]] - p1[iu[1]], axis=-1)
    exc = np.abs(d0 - STRUT_LEN)
    dec = d1 < d0
    disjoint = np.array([not (units[int(i)] & units[int(j)])
                         for i, j in zip(*iu)])
    hit = (exc < tol) & dec & disjoint
    pairs = []
    for c in np.nonzero(hit)[0]:
        i, j = int(iu[0][c]), int(iu[1][c])
        pairs.append((members[i][0], members[j][0]))
    worst, gap = _population_gap(list(exc[dec & disjoint]))
    return tuple(pairs), worst, gap


def member_length(phases, z, mem, topo):
    """One member's length in a dephased, solved configuration."""
    n = topo.n
    w = z[:3 * n].reshape(n, 3)
    t = z[3 * n:].reshape(n, 3)
    (ip, kp), (iq, kq) = mem
    p = _rodrigues(w[ip]) @ verts(phases[ip])[kp] + t[ip]
    q = _rodrigues(w[iq]) @ verts(phases[iq])[kq] + t[iq]
    return float(np.linalg.norm(p - q))


def span_rate_rows(topo, pairs, a):
    """(Gz, Gu): d(member length)/d(placement rate) and /d(phase rate).

    Exact first derivatives at the in-phase configuration. For an INTRA-unit
    member both endpoints ride the same rigid body, so the two placement
    contributions cancel and the whole Gz row is zero -- but that is COMPUTED
    here by the same expression as every other row and ASSERTED in Y3c, never
    special-cased, because the row's meaning depends on which case it is in and
    a hand-written zero would be a claim rather than a measurement.
    """
    v = verts(a)
    dv = dverts_exact(a)
    n = topo.n
    gz = np.zeros((len(pairs), 6 * n))
    gu = np.zeros((len(pairs), n))
    for r, ((ip, kp), (iq, kq)) in enumerate(pairs):
        sites = topo.sites(v)
        e = (v[kp] + sites[ip]) - (v[kq] + sites[iq])
        ne = float(np.linalg.norm(e))
        if ne < 1e-14:
            continue
        e = e / ne
        for (i, k), sgn in (((ip, kp), 1.0), ((iq, kq), -1.0)):
            gz[r, 3 * i:3 * i + 3] += sgn * (e @ (-_hat(v[k])))
            gz[r, 3 * n + 3 * i:3 * n + 3 * i + 3] += sgn * e
            gu[r, i] += sgn * float(e @ dv[k])
    return gz, gu


# ==========================================================================
# THE FIRST-ORDER PROGRAMMES
#
# Every first-order statement below is a LINEAR PROGRAMME, not a finite
# difference of a max. The excess is a MAX of linear functions of the rates, so
# its one-sided directional derivative in direction u is exactly
# max_m (Gz_m z + Gu_m u) minimised over the placement rates z that keep the
# contacts closed. That is linear programming and it is exact; differentiating
# it is what the tied-orbit hazard forbids.
# ==========================================================================

class _LPFail:
    """The shape a failed solve returns, so no caller has to branch on type.

    Status 4 is scipy's "numerical difficulties"; it is reused here for "the
    solver refused the inputs", which is a different thing from infeasible and
    from unbounded and must not be silently read as either.
    """
    status = 4
    fun = float("inf")
    x = None


def _safe_linprog(*args, **kwargs):
    """linprog that never raises.

    A malformed programme is a BUG REPORT, and it must arrive as a red gate row
    rather than as a traceback: a traceback destroys the verdict table and
    leaves a reader unable to tell a broken build from a measured result. Two
    mutation probes produced exactly that here -- one of them by making an
    equality block and its right-hand side disagree in length -- and neither
    reached the table.
    """
    try:
        return linprog(*args, **kwargs)
    except Exception:
        return _LPFail()


def _lp_core(gz, gu, jr, jp, n, extra_eq=(), extra_b=(), obj_u=None,
             ubound=1.0):
    """One LP over (placement rate z, dephasing rate u, slack s).

    Constraints always present: contacts stay closed (Jr z + Jp u = 0) and the
    dephasing is TRACELESS (sum u = 0), so nothing here can buy relief by simply
    changing the mean angle -- which is the whole question, and an LP that
    forgot the traceless row would answer it by uniform contraction and report a
    spurious relief. Returns (value, u, status).
    """
    nz = 6 * n
    m = gz.shape[0]
    nv = nz + n + 1
    a_eq = [np.hstack([jr, jp, np.zeros((jr.shape[0], 1))])] if jr.size else []
    b_eq = [np.zeros(jr.shape[0])] if jr.size else []
    row = np.zeros(nv)
    row[nz:nz + n] = 1.0
    a_eq.append(row)
    b_eq.append(np.zeros(1))
    for e, b in zip(extra_eq, extra_b):
        a_eq.append(np.asarray(e, float).reshape(1, -1))
        b_eq.append(np.atleast_1d(np.asarray(b, float)))
    a_ub = np.hstack([gz, gu, -np.ones((m, 1))])
    b_ub = np.zeros(m)
    c = np.zeros(nv)
    if obj_u is None:
        c[-1] = 1.0
    else:
        c[nz:nz + n] = obj_u
    bounds = [(None, None)] * nz + [(-ubound, ubound)] * n + [(None, None)]
    res = _safe_linprog(c, A_ub=a_ub, b_ub=b_ub,
                  A_eq=np.vstack(a_eq), b_eq=np.concatenate(b_eq),
                  bounds=bounds, method="highs")
    # STATUS 3 IS NOT AN ERROR AND MUST NOT BE LUMPED WITH STATUS 2. Unbounded
    # means the objective runs to minus infinity, that is, the excess can be
    # relieved at an unlimited rate -- which happens exactly when a MECHANISM
    # relieves it, since the placement rates are unbounded by design. Reading
    # that as "no answer" and printing +inf, as an earlier version did, inverts
    # the sign of the result on every tree topology.
    # The STATUS is what carries that distinction to the caller; this function
    # returns +inf for every non-optimal status and lets the caller key on the
    # code. An earlier version returned -inf here as well, and a mutation probe
    # showed the value was never read -- dead code wearing the appearance of a
    # decision.
    if res.status != 0:
        return float("inf"), None, int(res.status)
    return float(res.fun), np.asarray(res.x[nz:nz + n]).copy(), 0


def min_directional_derivative(topo, pairs, a):
    """min over |u|_inf = 1, sum u = 0, of the one-sided derivative of the excess.

    EXACT. The excess is positively homogeneous of degree one in the rates, so
    its minimum over the sup-SPHERE is attained at a point with some coordinate
    at +/-1; pinning each coordinate to each sign in turn and taking the best of
    the 2N programmes is therefore the exact minimum over the sphere and not a
    sample of it. A minimum over the sup-BALL would be zero by taking u = 0 and
    would answer nothing, which is the trap this shape exists to avoid.
    Returns (value, argmin, n_infeasible, n_unbounded). An UNBOUNDED programme
    contributes -inf to the minimum, which is the honest reading: the excess can
    be relieved at unlimited rate because a mechanism relieves it.
    """
    jr, jp = blocks(topo, a)
    gz, gu = span_rate_rows(topo, pairs, a)
    n = topo.n
    best, arg, bad, unb = float("inf"), None, 0, 0
    for p in range(n):
        for sgn in (1.0, -1.0):
            e = np.zeros(6 * n + n + 1)
            e[6 * n + p] = 1.0
            val, u, st = _lp_core(gz, gu, jr, jp, n, (e,), (sgn,))
            if st == 3:
                unb += 1
                best = -float("inf")
                continue
            if st != 0 or u is None:
                bad += 1
                continue
            if val < best:
                best, arg = val, u
    return best, arg, bad, unb


def directional_derivative(topo, pairs, a, u):
    """The one-sided derivative in ONE prescribed direction: min over placement
    rates of the max span rate. Exact, and the number the spread confirms."""
    jr, jp = blocks(topo, a)
    gz, gu = span_rate_rows(topo, pairs, a)
    n = topo.n
    eq, bb = [], []
    for i in range(n):
        e = np.zeros(6 * n + n + 1)
        e[6 * n + i] = 1.0
        eq.append(e)
        bb.append(float(u[i]))
    val, _, st = _lp_core(gz, gu, jr, jp, n, eq, bb,
                          ubound=max(1.0, float(np.abs(u).max()) * 1.0001))
    return val, st


def mechanism_only_relief(topo, pairs, a):
    """u = 0: can the MECHANISMS alone relieve the excess, with no dephasing?

    Normalised on the sup-sphere of the PLACEMENT rate, so the answer is a rate
    per unit placement rate. Strictly negative means the reference placement,
    and not the linkage, was holding that span taut -- which is exactly the
    correction jb_x's placement-conditional span enumeration needs.
    """
    jr, _ = blocks(topo, a)
    gz, _ = span_rate_rows(topo, pairs, a)
    nz = 6 * topo.n
    m = gz.shape[0]
    if m == 0:
        return float("inf"), 0
    best, bad = float("inf"), 0
    for p in range(nz):
        for sgn in (1.0, -1.0):
            a_eq = [np.hstack([jr, np.zeros((jr.shape[0], 1))])] if jr.size else []
            b_eq = [np.zeros(jr.shape[0])] if jr.size else []
            e = np.zeros(nz + 1)
            e[p] = 1.0
            a_eq.append(e.reshape(1, -1))
            b_eq.append(np.array([sgn]))
            c = np.zeros(nz + 1)
            c[-1] = 1.0
            res = _safe_linprog(c, A_ub=np.hstack([gz, -np.ones((m, 1))]),
                          b_ub=np.zeros(m), A_eq=np.vstack(a_eq),
                          b_eq=np.concatenate(b_eq),
                          bounds=[(-1.0, 1.0)] * nz + [(None, None)],
                          method="highs")
            if res.status != 0:
                bad += 1
                continue
            best = min(best, float(res.fun))
    return best, bad


# ==========================================================================
# Y0  CONTROL
# ==========================================================================

def y0_control():
    """Nothing downstream means anything if these fail."""
    print()
    print("=" * 78)
    print("Y0  CONTROL: the dephased model reduces to jb_x's in-phase model")
    print("=" * 78)
    out = {}
    topos = {t.name: t for t in build_topologies()}
    tt = topos["CUBE8-M (60/60/90 basis)"]
    mem, _ = _classes(tt)

    # (a) at ZERO dephasing the exact solver must reproduce the reference
    ph = np.full(tt.n, A_ICO_RECORD)
    z, res, its = solve_dephased(tt, ph)
    pts = class_points(tt, mem, ph, z)
    ref = _class_positions(A_ICO_RECORD, tt, mem)
    dev = float(np.abs(pts - ref).max())
    out["inphase_dev"] = dev
    out["inphase_ok"] = res < TOL["solve"] and dev < TOL["inphase"]
    # THE CONTROL THAT BOUNDS TOL[inphase] FROM ABOVE. Without it the tolerance
    # is bounded only below, by its own measured deviation, and an independent
    # validation loosened it to 1.0 with the gate green -- on one of the two
    # rows the whole file rests on. The reference is displaced by a stated
    # offset and the SAME comparison must now REJECT it.
    ref_ctrl = ref.copy()
    ref_ctrl[:, 0] += INPHASE_CONTROL_OFFSET
    out["inphase_ctrl"] = float(np.abs(pts - ref_ctrl).max())
    print(f"  (a) zero dephasing, CUBE8-M: contact residual {res:.2e} in {its} "
          f"iterations,")
    print(f"      class positions vs jb_x's reference placement  {dev:.2e}")
    print(f"      CONTROL: the same comparison against a reference displaced by "
          f"{INPHASE_CONTROL_OFFSET:.0e}")
    print(f"      must be REJECTED at the same tolerance: "
          f"{out['inphase_ctrl']:.2e}")

    # (b) analytic (Jr,Jp) vs central differences of the EXACT residual
    rngl = np.random.default_rng(SPREAD_SEED)
    jr, jp = blocks(tt, A_ICO_RECORD)
    n = tt.n
    zdir = rngl.standard_normal(6 * n) * 1e-3
    udir = rngl.standard_normal(n) * 1e-3
    udir -= udir.mean()
    h = 1e-6

    def exact_res(scale):
        phs = A_ICO_RECORD + scale * udir
        zz = np.zeros(6 * n)
        zz[3 * n:] = np.asarray(tt.sites(verts(A_ICO_RECORD)), float).reshape(-1)
        zz = zz + scale * zdir
        w = zz[:3 * n].reshape(n, 3)
        t = zz[3 * n:].reshape(n, 3)
        vv = [verts(x) for x in phs]
        return np.concatenate([
            _rodrigues(w[i]) @ vv[i][k] + t[i] - _rodrigues(w[j]) @ vv[j][l] - t[j]
            for (i, k, j, l) in tt.contacts])

    fd = (exact_res(h) - exact_res(-h)) / (2 * h)
    an = jr @ zdir + jp @ udir
    fdev = float(np.abs(fd - an).max())
    out["fd"] = fdev
    out["fd_ok"] = fdev < TOL["fd_jacobian"]
    # The same control for the OTHER unbounded-above foundation tolerance: the
    # analytic prediction scaled by a stated relative offset must FAIL the
    # comparison it passes unscaled.
    out["fd_ctrl"] = float(np.abs(fd - an * (1.0 + FD_CONTROL_OFFSET)).max())
    print(f"  (b) analytic (Jr,Jp) vs central differences of the exact nonlinear")
    print(f"      residual, along a seeded mixed direction   {fdev:.2e}")
    print(f"      CONTROL: the analytic prediction scaled by "
          f"(1 + {FD_CONTROL_OFFSET:.0e}) must be REJECTED at the same")
    print(f"      tolerance: {out['fd_ctrl']:.2e}")

    # (c) the two lock angles, RE-DERIVED here from this file's own spans
    dia = DIAGONALS[0] if DIAGONALS else None
    aico = float("nan")
    if dia is not None:
        def fdiag(x):
            v = verts(x)
            return float(np.linalg.norm(v[dia[0]] - v[dia[1]])) - STRUT_LEN
        try:
            aico = brentq(fdiag, 5.0, 40.0, xtol=1e-13)
        except Exception:
            # nan, not a raise; the nan then fails the comparison row below,
            # which is where a reader should learn about it.
            aico = float("nan")
    out["aico"] = aico
    out["aico_dev"] = abs(aico - A_ICO_RECORD) if np.isfinite(aico) else float("inf")
    out["aico_ctrl"] = (abs(aico - (A_ICO_RECORD + AICO_CONTROL_OFFSET))
                        if np.isfinite(aico) else 0.0)

    achord, nch, nfail = _rederive_chord(topos["CUBE8-M (60/60/90 basis)"], mem)
    out["chord_nfail"] = nfail
    out["achord"] = achord
    out["achord_dev"] = (abs(achord - A_CHORD_RECORD) if np.isfinite(achord)
                         else float("inf"))
    out["achord_ctrl"] = (abs(achord - (A_CHORD_RECORD + ACHORD_CONTROL_OFFSET))
                          if np.isfinite(achord) else 0.0)
    print(f"  (c) icosahedral angle re-derived from the folding diagonal   "
          f"{aico:.9f}")
    print(f"      recorded {A_ICO_RECORD:.9f}, deviation {out['aico_dev']:.2e}")
    print(f"      inter-unit chord angle re-derived on CUBE8-M             "
          f"{achord:.9f}")
    print(f"      recorded {A_CHORD_RECORD:.6f}, deviation "
          f"{out['achord_dev']:.2e}, from {nch} chords")

    # (d) the ACTIVE_TOL band, from the POPULATION GAP and therefore computed
    # without reference to ACTIVE_TOL itself. See `_population_gap`.
    pr_i, worst_i, gap_i = active_set(tt, mem, aico, "intra")
    pr_x, worst_x, gap_x = active_set(tt, mem, achord, "inter")
    out["n_intra"] = len(pr_i)
    out["tol_lo"] = max(worst_i, worst_x)
    out["tol_hi"] = min(gap_i, gap_x)
    out["tol_ok"] = out["tol_lo"] < ACTIVE_TOL < out["tol_hi"]
    # ... and the FLOOR the population gap clips at must not be what sets the
    # lower end, or the band is a statement about double precision rather than
    # about this geometry. Both halves are asserted: the lower end is a real
    # member of the population and strictly above the clip, and the gap it
    # brackets is many decades wide rather than a rounding step.
    out["tol_lo_not_floor"] = out["tol_lo"] > EXCESS_FLOOR
    out["tol_gap_decades"] = (
        float(np.log10(out["tol_hi"] / out["tol_lo"]))
        if out["tol_lo"] > 0 and np.isfinite(out["tol_hi"]) else 0.0)
    out["tol_gap_wide"] = out["tol_gap_decades"] > EXCESS_GAP_DECADES
    print(f"  (d) active sets on CUBE8-M: {len(pr_i)} intra at the icosahedral "
          f"phase,")
    print(f"      {len(pr_x)} inter at the chord angle; ACTIVE_TOL band "
          f"{out['tol_lo']:.2e} < {ACTIVE_TOL:.0e} < {out['tol_hi']:.2e}")
    print(f"      (the band is the widest log gap in the population of "
          f"candidate excesses,")
    print(f"      filtered by NO tolerance, so it does not move when "
          f"ACTIVE_TOL moves. Its")
    print(f"      lower end {out['tol_lo']:.2e} is the diagonals' own residual, "
          f"far above the")
    print(f"      double-precision floor {EXCESS_FLOOR:.0e}, so the floor is "
          f"not what bounds it.)")

    # (d2) R1's OTHER recorded number for this chord: its length at the
    # icosahedral phase, as a multiple of the strut. Re-derived from the same
    # pair the angle came from. Without this the chord is identified by its
    # taut ANGLE alone, and two different spans could share an angle.
    ratio = float("nan")
    if pr_x:
        (ip, kp), (iq, kq) = pr_x[0]
        v = verts(aico)
        st_sites = tt.sites(v)
        ratio = float(np.linalg.norm((v[kp] + st_sites[ip])
                                     - (v[kq] + st_sites[iq]))) / STRUT_LEN
    out["chord_ratio"] = ratio
    out["chord_ratio_dev"] = (abs(ratio - CHORD_RATIO_RECORD)
                              if np.isfinite(ratio) else float("inf"))
    out["chord_ratio_ok"] = out["chord_ratio_dev"] < 5e-5
    print(f"      the same chord's length at the icosahedral phase is "
          f"{ratio:.6f} x strut,")
    print(f"      recorded {CHORD_RATIO_RECORD} (four decimals), deviation "
          f"{out['chord_ratio_dev']:.2e}")

    # (e) the recorded per-cluster chord counts, re-derived
    counts = {}
    for name, want in CHORD_COUNT_RECORD.items():
        t2 = topos.get(name)
        if t2 is None:
            counts[name] = (-1, want)
            continue
        m2, _ = _classes(t2)
        p2, _, _ = active_set(t2, m2, achord, "inter")
        counts[name] = (len(p2), want)
    out["counts"] = counts
    out["counts_ok"] = bool(counts) and all(g == w for g, w in counts.values())
    # NON-VACUITY for the T2 cross-check. `all()` over an empty dict is True: a
    # probe that emptied CHORD_COUNT_RECORD reddened NOTHING, so the whole
    # 8/8/6/48 comparison with the recorded table was deletable with the gate
    # green. This is the `all()`-over-empty shape the sibling file was audited
    # for, present here and unguarded until now.
    out["counts_n"] = len(counts)
    out["counts_nonvacuous"] = (len(counts) >= 4
                                and all(w > 0 for _, w in counts.values())
                                and all(g >= 0 for g, _ in counts.values()))
    print("  (e) inter-unit spans taut at the chord angle, re-derived vs T2 R1:")
    for name, (got, want) in counts.items():
        print(f"      {name:34s} {got:4d}   recorded {want:4d}")

    # (g) THE TRUNCATED GAUSS-NEWTON STEP, exercised rather than assumed.
    # A self-stressed cluster's contact Jacobian is exactly rank-deficient at
    # zero dephasing and APPROXIMATELY so just off it, so the vanishing singular
    # values sit at O(eps) and the untruncated least-squares step is their noise
    # divided by them. Both halves are asserted, in the same two-row idiom the
    # chord and the icosahedral angle carry: the shipped solver must CLOSE this
    # configuration and the untruncated one must FAIL to.
    trc = Topology("rcond", "box", (0, 1, 3)[:len(RCOND_BOX)], box=RCOND_BOX)
    rsites = np.array(trc.lattice_sites, float)
    ru = rsites[:, 0] - rsites[:, 0].mean()
    rmx = float(np.abs(ru).max())
    if rmx > 1e-12:
        rph = aico + RCOND_EPS * (ru / rmx)
        _, r_lad, _ = solve_dephased(trc, rph)
        _, r_pln, _ = solve_dephased(trc, rph, rconds=(None,))
    else:
        r_lad, r_pln = float("inf"), 0.0
    out["rcond_ladder_res"] = r_lad
    out["rcond_plain_res"] = r_pln
    print(f"  (g) truncated Gauss-Newton step on {RCOND_BOX} at eps = "
          f"{RCOND_EPS:g} deg along an")
    print(f"      affine ramp: the shipped escalating cutoff closes to "
          f"{r_lad:.2e},")
    print(f"      the UNTRUNCATED solve stalls at {r_pln:.2e}. Both halves are "
          f"gated, so")
    print("      removing the escalation reddens a row instead of passing "
          "unnoticed.")

    # (f) six intra-unit diagonals per unit, and the count is per-unit exactly
    out["per_unit_ok"] = (len(pr_i) == 6 * tt.n) and len(DIAGONALS) == 6
    print(f"  (f) intra actives = 6 per unit: {len(pr_i)} = 6 x {tt.n}? "
          f"{len(pr_i) == 6 * tt.n};  DIAGONALS carries {len(DIAGONALS)}")
    return out


def _rederive_chord(topo, members):
    """The inter-unit chord angle, by root-finding rather than by a literal.

    Brackets around the RECORDED value are deliberately wide (a whole degree
    either side) so that a wrong record shows up as disagreement rather than as
    a bracket that quietly excludes the truth. Returns (angle, n_chords), with
    nan and 0 if nothing crosses -- never a raise.
    """
    lo, hi = A_CHORD_RECORD - 1.0, A_CHORD_RECORD + 1.0
    _, units = _classes(topo)
    p0 = _class_positions(A_CHORD_RECORD, topo, members)
    nc = len(members)
    iu = np.triu_indices(nc, 1)
    d0 = np.linalg.norm(p0[iu[0]] - p0[iu[1]], axis=-1)
    near = np.nonzero(np.abs(d0 - STRUT_LEN) < 1e-3)[0]
    roots, nfail = [], 0
    for c in near:
        i, j = int(iu[0][c]), int(iu[1][c])
        if units[i] & units[j]:
            continue

        def f(x):
            return span_length(topo, members, x, i, j) - STRUT_LEN

        try:
            if f(lo) * f(hi) > 0.0:
                continue
            roots.append(brentq(f, lo, hi, xtol=1e-13))
        except Exception:
            # A DROPPED ROOT, never a raise -- but COUNTED. An excepting
            # refinement that is silently continued is data loss with no trace,
            # so the count comes back and is gated at zero, and the swallow
            # cannot grow quietly. Same discipline as jb_x's crossing search.
            nfail += 1
            continue
    if not roots:
        return float("nan"), 0, nfail
    return float(np.median(roots)), len(roots), nfail


# ==========================================================================
# Y1  THE EQUALITY OBSTRUCTION: which dephasings the linkage can execute at all
# ==========================================================================

def y1_admissible():
    print()
    print("=" * 78)
    print("Y1  WHICH DEPHASINGS ARE KINEMATICALLY ADMISSIBLE AT ALL")
    print("=" * 78)
    print("  Before asking whether dephasing RELIEVES anything, ask whether the")
    print("  linkage can execute it. A dephasing rate u is admissible when some")
    print("  placement rate solves Jr z + Jp u = 0 exactly. The obstruction is")
    print("  the STATES OF SELF-STRESS of the placement problem acting on the")
    print("  phase columns, so a cluster with no self-stress can dephase freely")
    print("  and one with self-stress cannot. This is a NEW constraint: jb_x")
    print("  never varied a phase, so it could not see it.")
    out = {}

    print()
    print("  (a) the named topologies, at the icosahedral phase")
    print(f"      {'topology':34s} {'N':>4s} {'rows':>5s} {'rank':>5s} "
          f"{'self':>5s} {'obstr':>6s} {'dephase-dim':>12s}")
    per = {}
    for t in build_topologies():
        if t.name.startswith("FCC13"):
            continue
        r = admissible_dephasing(t, A_ICO_RECORD)
        per[t.name] = r
        print(f"      {t.name:34s} {t.n:4d} {r['rows']:5d} {r['rank']:5d} "
              f"{r['selfstress']:5d} {r['obstruct']:6d} "
              f"{r['dim']:6d}/{t.n:<5d}")
    out["per"] = per

    print()
    print("  (b) box scaling on the M basis: where the obstruction switches on,")
    print("      and WHAT switches it on. Two candidate graph properties are")
    print("      computed alongside -- the CYCLE RANK and the largest number of")
    print("      contacts at any one unit -- so the discriminator is MEASURED")
    print("      against the table rather than asserted from one reading of it.")
    print(f"      {'box':12s} {'N':>5s} {'cycles':>7s} {'maxcoord':>9s} "
          f"{'selfstress':>11s} {'dephase-dim':>13s}")
    scal, graph = {}, {}
    for box in SCALING_BOXES:
        gens = (0, 1, 3)[:len(box)]
        t = Topology("scal", "box", gens, box=box)
        r = admissible_dephasing(t, A_ICO_RECORD)
        betti, mx, nhi = contact_graph(t)
        scal[box] = (t.n, r["selfstress"], r["dim"])
        graph[box] = (betti, mx, nhi)
        print(f"      {str(box):12s} {t.n:5d} {betti:7d} {mx:9d} "
              f"{r['selfstress']:11d} {r['dim']:6d}/{t.n:<6d}")
    out["scal"] = scal
    out["graph"] = graph
    free = [b for b, (n, ss, d) in scal.items() if d == n]
    obstructed = [b for b, (n, ss, d) in scal.items() if d < n]
    out["free_boxes"] = free
    out["obstructed_boxes"] = obstructed
    out["selfstress_is_the_discriminator"] = all(
        (ss == 0) == (d == n) for (n, ss, d) in scal.values())
    print(f"      UNOBSTRUCTED: {free}")
    print(f"      OBSTRUCTED  : {obstructed}")
    print("      Every unobstructed box has ZERO states of self-stress and every")
    print("      obstructed one has at least one, so the self-stress is the")
    print("      whole discriminator -- BUT read the direction that carries")
    print("      information. `selfstress == 0 => full dimension` is a CODE")
    print("      BRANCH: with no self-stress `admissible_dephasing` returns")
    print("      early with dim = N, so those rows are an identity and not a")
    print("      measurement. The informative half is the OTHER one -- that the")
    print("      obstruction rank equals the self-stress count exactly, on every")
    print("      stressed box -- and that is what the gate asserts separately.")
    print()
    print("      AND IT IS NOT A CYCLE PROPERTY. A previous version of this")
    print("      paragraph said it was; the table refutes it one box further")
    print("      out. Both graph conditions are NECESSARY and NEITHER alone is")
    print("      SUFFICIENT, and each has a witness in this run:")
    cyc_no_ss = [b for b in scal
                 if graph[b][0] > 0 and scal[b][1] == 0]
    coord_no_ss = [t.name for t in build_topologies()
                   if not t.name.startswith("FCC13")
                   and contact_graph(t)[1] >= COORD_INTERIOR
                   and admissible_dephasing(t, A_ICO_RECORD)["selfstress"] == 0]
    out["cycles_without_selfstress"] = cyc_no_ss
    out["coord_without_selfstress"] = coord_no_ss
    print(f"        CYCLES WITHOUT SELF-STRESS: {cyc_no_ss}")
    print(f"        COORDINATION >= {COORD_INTERIOR} WITHOUT SELF-STRESS: "
          f"{coord_no_ss}")
    print("      The star has a unit joined to six neighbours and is a TREE; the")
    print("      2-thick boxes have cycles in abundance and no unit joined to")
    print("      more than three. Neither carries self-stress. What does is the")
    print("      CONJUNCTION -- a cycle AND a unit with at least four contacts,")
    print("      which is to say a cluster that has stopped being all boundary.")
    named = {}
    for t in build_topologies():
        if t.name.startswith("FCC13"):
            continue
        named[t.name] = (contact_graph(t),
                         admissible_dephasing(t, A_ICO_RECORD)["selfstress"])
    conj = [(graph[b][0] > 0 and graph[b][1] >= COORD_INTERIOR)
            == (scal[b][1] > 0) for b in scal]
    conj += [(g[0] > 0 and g[1] >= COORD_INTERIOR) == (ss > 0)
             for (g, ss) in named.values()]
    out["conjunction_hits"] = int(sum(conj))
    out["conjunction_n"] = len(conj)
    out["conjunction_ok"] = bool(conj) and all(conj)
    print(f"      MEASURED on all {len(conj)} clusters in this file -- the "
          f"{len(scal)} boxes above")
    print(f"      and the {len(named)} named topologies: "
          f"(cycle AND coordination >= {COORD_INTERIOR}) agrees with")
    print(f"      (self-stress > 0) on {int(sum(conj))} of {len(conj)}.")
    print("      Reported as a measured correlation over this table. No")
    print("      derivation is offered and none is claimed.")
    print("      2x2x2 is the largest box here with unobstructed dephasing, and")
    print("      2x2x2 is the size of the owner's array. That is a coincidence")
    print("      this file NOTICES and does not lean on: nothing measured here")
    print("      shows his eight cells are an M-basis 2x2x2 box.")
    # THE OBSTRUCTION RANK vs THE SELF-STRESS COUNT, on the STRESSED boxes only
    # -- the half of the discriminator that is NOT a code branch and therefore
    # the half that carries information.
    obs = {}
    for b in scal:
        if scal[b][1] == 0:
            continue
        t = Topology("s", "box", (0, 1, 3)[:len(b)], box=b)
        r = admissible_dephasing(t, A_ICO_RECORD)
        obs[b] = (r["selfstress"], r["obstruct"], t.n)
    out["obstruct"] = obs
    out["n_stressed"] = len(obs)
    out["obstruct_positive"] = bool(obs) and all(o > 0 for _, o, _ in
                                                 obs.values())
    out["obstruct_eq_selfstress"] = [b for b, (ss, o, _) in obs.items()
                                     if o == ss]
    print(f"      {'box':12s} {'selfstress':>11s} {'obstruction rank':>17s}")
    for b, (ss, o, _) in obs.items():
        print(f"      {str(b):12s} {ss:11d} {o:17d}")
    print(f"      The obstruction rank is STRICTLY POSITIVE on all "
          f"{len(obs)} stressed boxes,")
    print(f"      and equals the self-stress count exactly on "
          f"{len(out['obstruct_eq_selfstress'])} of them; on the")
    print("      largest it is smaller, because the phase columns cannot")
    print("      obstruct more directions than the phase space has. That is")
    print("      measured, and it is the same saturation (c1) reports.")

    print()
    print("  (c1) DOES THE DIMENSION GROW WITH THE ARRAY? It does not: it")
    print("       SATURATES. The previous version of this table stopped at")
    print("       3x3x3 and could not see that.")
    print(f"       {'cube':10s} {'N':>5s} {'selfstress':>11s} "
          f"{'dephase-dim':>13s} {'min retention':>14s} {'CONTROL max':>14s}")
    sat = {}
    for box in SATURATION_BOXES:
        t = Topology("sat", "box", (0, 1, 3), box=box)
        r = admissible_dephasing(t, A_ICO_RECORD)
        sites = np.array(t.lattice_sites, float)
        chi = (-1.0) ** sites.sum(axis=1)
        cand = (np.ones(t.n), sites[:, 0], sites[:, 1], sites[:, 2],
                chi, chi * sites[:, 0], chi * sites[:, 1], chi * sites[:, 2])
        # THE CONTROL, without which "retained at 1.000000000" is satisfied by
        # a projector that retains everything -- i.e. by no obstruction at all.
        # A probe replacing this projector with the identity reddened NOTHING
        # until this column existed. Same companion (d) already carries at
        # 3x3x3, applied at every saturation size.
        cctl = ((-1.0) ** sites[:, 0], (-1.0) ** sites[:, 1],
                (-1.0) ** sites[:, 2])
        proj = r["basis"] @ r["basis"].T
        mret = min(float(np.dot(w, proj @ w) / np.dot(w, w)) for w in cand)
        xret = max(float(np.dot(w, proj @ w) / np.dot(w, w)) for w in cctl)
        sat[box] = (t.n, r["selfstress"], r["dim"], mret, xret)
        print(f"       {str(box):10s} {t.n:5d} {r['selfstress']:11d} "
              f"{r['dim']:6d}/{t.n:<6d} {mret:14.9f} {xret:14.9f}")
    out["sat"] = sat
    dims = {v[2] for v in sat.values()}
    out["sat_dim"] = sorted(dims)
    out["sat_constant"] = len(dims) == 1 and len(sat) >= 3
    out["sat_n_grows"] = len({v[0] for v in sat.values()}) == len(sat)
    out["sat_min_ret"] = min((v[3] for v in sat.values()), default=0.0)
    out["sat_ret_ok"] = bool(sat) and out["sat_min_ret"] > 1.0 - 1e-9
    out["sat_ctrl_max"] = max((v[4] for v in sat.values()), default=1.0)
    out["sat_ctrl_ok"] = bool(sat) and out["sat_ctrl_max"] < 0.9
    print("       The dimension is CONSTANT while N runs 27 -> 64 -> 125, and")
    print("       all eight of AFFINE x PARITY are retained at 1.000000000 at")
    print("       every size. ADMISSIBLE DEPHASING IS A FIXED FINITE-DIMENSIONAL")
    print("       AFFINE-TIMES-PARITY SPACE INDEPENDENT OF ARRAY SIZE, so the")
    print("       admissible FRACTION goes to zero as the array grows. That is")
    print("       the strongest scaling statement in this file and it is on the")
    print("       bead's own question: an array large enough to be a medium has")
    print("       essentially no dephasing available to it that the linkage can")
    print("       execute, and what it does have is a uniform gradient and an")
    print("       alternation and nothing between them.")
    print("       SCOPE: three cubes on the M basis at the icosahedral phase.")
    print("       Other bases, other cell shapes and other angles are not swept")
    print("       here; (c) below varies the angle and the basis at fixed size.")

    print()
    print("  (c) topology and angle dependence of the dephasing dimension")
    print(f"      {'basis':6s} {'box':10s} " + " ".join(f"{a:>9.4f}" for a in ADM_ANGLES))
    angdep = {}
    for gens, tag in (((0, 1, 3), "M"), ((0, 1, 2), "R")):
        for box in ((3, 3, 3), (2, 3, 3)):
            t = Topology("ad", "box", gens, box=box)
            dims = []
            for a in ADM_ANGLES:
                dims.append(admissible_dephasing(t, a)["dim"])
            angdep[(tag, box)] = tuple(dims)
            print(f"      {tag:6s} {str(box):10s} "
                  + " ".join(f"{d:9d}" for d in dims))
    out["angdep"] = angdep
    out["angle_independent"] = bool(angdep) and all(
        len(set(v)) == 1 for v in angdep.values())
    out["basis_sensitive"] = any(
        angdep[("M", b)] != angdep[("R", b)] for b in ((3, 3, 3), (2, 3, 3)))
    # NON-VACUITY for "angle-independent". `all(len(set(v)) == 1 ...)` over a
    # single angle, or over the same angle twice, is True for any routine
    # whatever: probes that cut ADM_ANGLES to one angle and to two IDENTICAL
    # angles each reddened NOTHING. The file added exactly this companion for
    # the direction spread ("over a spread of at least 8 directions") and did
    # not add it here. The SPREAD is asserted as well as the count, because
    # three angles a millidegree apart are one angle with rounding.
    out["n_adm_angles"] = len(set(ADM_ANGLES))
    out["adm_angle_spread"] = (max(ADM_ANGLES) - min(ADM_ANGLES)
                               if ADM_ANGLES else 0.0)
    out["adm_angles_ok"] = (len(set(ADM_ANGLES)) >= 3
                            and out["adm_angle_spread"] > 10.0)
    print("      The dimension does NOT move with the phase and DOES move with")
    print("      the basis, so it is a property of the contact topology. The M")
    print("      basis -- the one carrying the diagonal generator pair, hence")
    print("      the taut inter-unit chord -- is the more constrained of the two.")

    print()
    print("  (d) what the admissible dephasings ARE, on a 3x3x3 M cluster")
    t = Topology("id", "box", (0, 1, 3), box=(3, 3, 3))
    r = admissible_dephasing(t, A_ICO_RECORD)
    sites = np.array(t.lattice_sites, float)
    chi = (-1.0) ** sites.sum(axis=1)
    cand = {"1": np.ones(t.n), "s0": sites[:, 0], "s1": sites[:, 1],
            "s2": sites[:, 2], "chi": chi, "chi*s0": chi * sites[:, 0],
            "chi*s1": chi * sites[:, 1], "chi*s2": chi * sites[:, 2]}
    proj = r["basis"] @ r["basis"].T
    ret = {}
    for nm, w in cand.items():
        ret[nm] = float(np.dot(w, proj @ w) / np.dot(w, w))
    # a control: a pattern that must NOT be fully retained
    ctrl = {}
    for nm, kk in (("k=(pi,0,0)", (1, 0, 0)), ("k=(0,pi,0)", (0, 1, 0)),
                   ("k=(pi,pi,0)", (1, 1, 0))):
        w = (-1.0) ** (kk[0] * sites[:, 0] + kk[1] * sites[:, 1]
                       + kk[2] * sites[:, 2])
        ctrl[nm] = float(np.dot(w, proj @ w) / np.dot(w, w))
    out["retention"] = ret
    out["retention_ctrl"] = ctrl
    out["affine_parity_ok"] = all(abs(v - 1.0) < 1e-9 for v in ret.values())
    out["ctrl_ok"] = bool(ctrl) and all(v < 0.9 for v in ctrl.values())
    # THE RANK_RTOL BAND, MEASURED AND INDEPENDENT OF RANK_RTOL. Every dimension
    # in this section is a rank, and a rank is a threshold on a singular-value
    # ratio; the threshold has to lie inside the spectrum's own GAP. An earlier
    # version bracketed it by "largest ratio discarded" and "smallest kept",
    # which are both computed WITH RANK_RTOL -- so the band moved with the
    # constant it was guarding and a probe that loosened RANK_RTOL by nine
    # decades left the row green. That is the shape jb_x's post-mortem names: a
    # guard band that improves when you widen it is not a guard. The widest gap
    # in the spectrum is a property of the matrix alone.
    jr333, jp333 = blocks(t, A_ICO_RECORD)
    rr333, _ = rank_of(jr333)
    u333, _, _ = np.linalg.svd(jr333, full_matrices=True)
    lo1, hi1 = _spectral_gap(jr333)
    lo2, hi2 = _spectral_gap(u333[:, rr333:].T @ jp333)
    out["rtol_lo"] = max(lo1, lo2)
    out["rtol_hi"] = min(hi1, hi2)
    out["rtol_ok"] = out["rtol_lo"] < RANK_RTOL < out["rtol_hi"]
    print(f"      RANK_RTOL band from the SPECTRAL GAP of both matrices this")
    print(f"      section ranks: {out['rtol_lo']:.2e} < {RANK_RTOL:.0e} < "
          f"{out['rtol_hi']:.2e}")
    print("      (the gap is a property of the matrices, so this band does not")
    print("      move when RANK_RTOL does.)")
    print("      retained fraction of a candidate pattern in the admissible space")
    for nm, v in ret.items():
        print(f"        {nm:10s} {v:.9f}")
    print("      CONTROLS that must NOT be retained (single-direction parity):")
    for nm, v in ctrl.items():
        print(f"        {nm:12s} {v:.9f}")
    print(f"      All eight of AFFINE x PARITY are retained exactly; the space")
    print(f"      has dimension {r['dim']}, so eight of it is identified and four")
    print("      are not. The identified part is the physically interesting part:")
    print("      a UNIFORM PHASE GRADIENT along any lattice direction is")
    print("      admissible, and so is the zone-corner alternation, and nothing")
    print("      of intermediate wavelength is.")

    print()
    print("  (e) the same statement in the BULK, by a Bloch calculation, AT")
    print("      FIXED LATTICE PERIOD -- and that scope is load-bearing, not")
    print("      decoration. Read (e3) before reading anything into (e1).")
    print("      One unit per cell, dephasing rate d exp(i k.s), translation")
    print("      rate t exp(i k.s). The contact equation collapses to")
    print("          -(1+z_m) [v_m]x w  +  (1-z_m) t  +  (1+z_m) v'_m d  =  0,")
    print("      z_m = exp(i k_m), for the three generators m. A plane-wave")
    print("      dephasing is admissible exactly when the phase column lies in")
    print("      the range of the placement block. That is three COMPLEX")
    print("      conditions, i.e. SIX REAL conditions, on a THREE-REAL torus, so")
    print("      the generic expectation is that the solution set is EMPTY --")
    print("      not 'isolated points', which an earlier version of this")
    print("      paragraph said and which invites a reader to conclude the scan")
    print("      was merely too coarse. One point survives for a DERIVABLE")
    print("      reason and not by a coincidence of counting: at k = (pi,pi,pi)")
    print("      every (1+z_m) vanishes, the phase column is identically zero")
    print("      and d is free.")
    print()
    print("      (e1) the scan")
    scan = _bloch_scan(K_GRID, K_OFFSET)
    scan2 = _bloch_scan(K_GRID_ALT, K_OFFSET_ALT)
    out["bloch_min"] = scan["min"]
    out["bloch_min_alt"] = scan2["min"]
    out["bloch_degen"] = scan["degenerate"] + scan2["degenerate"]
    cr, cnb = _bloch_residual(np.array([np.pi] * 3))
    xr, xnb = _bloch_residual(np.array([np.pi, np.pi, 3.0]))
    out["bloch_corner_nb"] = cnb
    out["bloch_corner_ctrl"] = xr
    out["bloch_ctrl_nb"] = xnb
    print(f"      incommensurate {K_GRID}^3 scan: min RELATIVE residual "
          f"{scan['min']:.6e}")
    print(f"        at k/pi = {np.round(scan['argmin'] / np.pi, 5)}, and "
          f"{scan['degenerate']} grid points had ||b|| under the floor")
    print(f"      second grid {K_GRID_ALT}^3, different irrational offset: min "
          f"{scan2['min']:.6e}, {scan2['degenerate']} degenerate")
    print(f"      k = (pi,pi,pi) exactly:          rel {cr:.3e}, ||b|| "
          f"{cnb:.3e}  -> FREE, and free because b vanishes")
    print(f"      k = (pi,pi,3.0) as the control:  rel {xr:.3e}, ||b|| "
          f"{xnb:.3e}  -> obstructed")
    # THE TWO GRIDS MUST BE DIFFERENT GRIDS. "Measured on two incommensurate
    # k-grids" is worth nothing if the second is the first under another name,
    # and no row can see that from the two minima alone -- they would simply
    # agree. Probes that made the alternative grid IDENTICAL to the first, and
    # that reduced both to 5^3, each reddened NOTHING. The amplitude ladders
    # carry exactly this companion ("the two ladders share no rung"); the
    # k-grids lacked its analogue until now.
    out["kgrid_distinct"] = (K_GRID != K_GRID_ALT
                             and np.gcd(int(K_GRID), int(K_GRID_ALT)) == 1
                             and abs(K_OFFSET - K_OFFSET_ALT) > 1e-3)
    out["kgrid_min"] = min(int(K_GRID), int(K_GRID_ALT))
    out["kgrid_resolved"] = out["kgrid_min"] >= 15
    print(f"      the two grids are coprime ({K_GRID}, {K_GRID_ALT}), their "
          f"offsets differ by")
    print(f"      {abs(K_OFFSET - K_OFFSET_ALT):.4f}, and the coarser carries "
          f"{min(K_GRID, K_GRID_ALT)}^3 points.")

    print()
    print("      (e2) THE CONTROL THE SCAN NEVER RAN: k = 0. The scan offsets")
    print("      guarantee no grid point lands on it, so the one wavevector")
    print("      whose answer is independently known was the one never")
    print("      evaluated. Evaluated here.")
    zr, znb = _bloch_residual(np.zeros(3))
    out["bloch_k0"] = zr
    out["bloch_k0_nb"] = znb
    print(f"      k = (0,0,0) at FIXED lattice period:  rel {zr:.6e}, ||b|| "
          f"{znb:.3e}")
    print("      -> INADMISSIBLE. And a uniform dephasing is the most obviously")
    print("      executable motion there is, which is the tell.")

    print()
    print("      (e3) WHAT THE ANSATZ LEAVES OUT, AND WHY (e1) IS SCOPED. The")
    print("      only translation freedom above is t exp(i k.s), which at k = 0")
    print("      is a RIGID translation of every cell by the same vector. There")
    print("      is no unknown for a HOMOGENEOUS LATTICE STRAIN -- a field")
    print("      t(s) = F s, which is not a plane wave of any k and so lives")
    print("      outside the ansatz by construction. Measured on the finite")
    print("      clusters, uniform dephasing IS executable and IS executed by")
    print("      exactly such a field:")
    print(f"      {'cluster':10s} {'N':>5s} {'lstsq residual':>15s} "
          f"{'affine-fit rel':>15s} {'|F - (tr F/3) I|':>17s}")
    aff = {}
    for box in ((2, 2, 2), (3, 3, 3), (4, 4, 4)):
        t = Topology("aff", "box", (0, 1, 3), box=box)
        jr, jp = blocks(t, A_ICO_RECORD)
        n = t.n
        b = -jp @ np.ones(n)
        zz, _, _, _ = np.linalg.lstsq(jr, b, rcond=None)
        lsres = float(np.linalg.norm(jr @ zz - b))
        tv = zz[3 * n:].reshape(n, 3)
        sites = np.array(t.lattice_sites, float)
        des = np.hstack([sites, np.ones((n, 1))])
        coef, _, _, _ = np.linalg.lstsq(des, tv, rcond=None)
        rel = (float(np.linalg.norm(tv - des @ coef))
               / max(float(np.linalg.norm(tv)), 1e-300))
        fmat = coef[:3].T
        dev = float(np.abs(fmat - np.trace(fmat) / 3.0 * np.eye(3)).max())
        aff[box] = (n, lsres, rel, dev)
        print(f"      {str(box):10s} {n:5d} {lsres:15.2e} {rel:15.2e} "
              f"{dev:17.3e}")
    out["affine"] = aff
    out["affine_exec_ok"] = bool(aff) and all(v[1] < 1e-9 for v in aff.values())
    out["affine_fit_ok"] = bool(aff) and all(v[2] < 1e-9 for v in aff.values())
    out["affine_nonscalar"] = min((v[3] for v in aff.values()), default=0.0)
    out["affine_nonscalar_ok"] = bool(aff) and out["affine_nonscalar"] > 1e-3
    out["affine_worst_exec"] = max((v[1] for v in aff.values()), default=1.0)
    out["affine_worst_fit"] = max((v[2] for v in aff.values()), default=1.0)
    print("      The translation field is affine to the last bit, and its")
    print("      linear part is NOT a multiple of the identity, so it is a")
    print("      genuine STRAIN and not a rescaling of the lattice. THAT is the")
    print("      freedom (e1) has no unknown for.")
    print()
    print("      SO THE RECONCILIATION IS THE LATTICE PERIOD, NOT THE BOUNDARY.")
    print("      An earlier version of this paragraph said a finite cluster is")
    print("      more permissive because its boundary relaxes the self-stress.")
    print("      That is REFUTED by the table above: the uniform pattern is")
    print("      retained at 1.000000000 at 3x3x3, 4x4x4 and 5x5x5 alike, so")
    print("      the retention does not decay with size and the boundary is not")
    print("      what is doing it. What (e1) measures is that AT FIXED LATTICE")
    print("      PERIOD the only admissible bulk plane wave is k = (pi,pi,pi).")
    print("      The uniform mode is admissible too, and admissible only by")
    print("      straining the lattice -- which this file's own MODEL")
    print("      DECLARATION says is never fixed by hand, so the fixed-period")
    print("      scope is a property of THIS CALCULATION and not of the model.")
    return out


def _spectral_gap(m, floor=1e-18):
    """(below, above): the ratio values bracketing the WIDEST log gap.

    Computed from the singular values alone, with no rank threshold anywhere,
    so the interval it returns cannot move when a rank tolerance moves. The
    spectrum is padded at the bottom with `floor` so that a full-rank matrix
    still reports the gap between its smallest ratio and numerical zero rather
    than an empty interval.
    """
    if m.size == 0:
        return 0.0, 1.0
    s = np.linalg.svd(m, compute_uv=False)
    if not s.size or not np.isfinite(s).all() or s[0] <= 0.0:
        return 0.0, 1.0
    r = np.concatenate([s / s[0], [floor]])
    r = np.clip(r, floor, None)
    lg = np.log10(r)
    d = lg[:-1] - lg[1:]
    if not d.size:
        return 0.0, 1.0
    i = int(np.argmax(d))
    return float(r[i + 1]), float(r[i])


def _bloch_residual(k, a=None, gens=(0, 1, 3)):
    """(relative residual, ||b||) for a plane-wave dephasing of wavevector k.

    Admissible means the phase column b lies in the range of the placement
    block A. The RELATIVE residual ||(I-P)b|| / ||b|| is the scale-free measure
    and is the one reported, with ||b|| returned alongside because it degenerates
    at the zone corner: there every (1+z_m) vanishes, b is IDENTICALLY ZERO, and
    the mode is free for the trivial reason rather than by cancellation. An
    absolute residual conflates the two -- it is small both where b lies in the
    range and where b is merely small -- and an earlier version of this row used
    exactly that and reported its minimum NEAR the corner rather than AT it.
    Callers must read both numbers: free means ||b|| below the floor OR the
    ratio at zero.
    """
    aa = A_ICO_RECORD if a is None else a
    v = verts(aa)
    dv = dverts_exact(aa)
    amat = np.zeros((9, 6), dtype=complex)
    b = np.zeros(9, dtype=complex)
    for m, p in enumerate(gens):
        z = np.exp(1j * k[m])
        amat[3 * m:3 * m + 3, 0:3] = -(1 + z) * _hat(v[p])
        amat[3 * m:3 * m + 3, 3:6] = (1 - z) * np.eye(3)
        b[3 * m:3 * m + 3] = -(1 + z) * dv[p]
    nb = float(np.linalg.norm(b))
    u, s, _ = np.linalg.svd(amat)
    r = int((s > s[0] * RANK_RTOL).sum()) if s.size and s[0] > 0 else 0
    res = float(np.linalg.norm(u[:, r:].conj().T @ b))
    if nb < BLOCH_B_FLOOR:
        return 0.0, nb
    return res / nb, nb


@jb_cache.memoize(_MODULE)
def _bloch_scan(grid, offset):
    """MEMOISED (jb_cache): 62,823 `_bloch_residual` evaluations across two
    calls, 43.9s of a 37.7s-wall gate -- the whole cost of this file. Pure in
    (grid, offset) and in the module constants the residual reads, so it is
    cached on the transitive source closure of all of them and computed in
    parallel ahead of the serial pass. `--no-cache` bypasses both and must
    print identically.

    Incommensurate scan of the three-torus. Never raises.

    Reports the minimum RELATIVE residual over grid points whose ||b|| is above
    the floor, and separately how many points fell below it -- which must be
    zero on an incommensurate grid, since the only such point is the zone
    corner and no incommensurate grid contains it. That count is the row's own
    check that the grid really is off-lattice.
    """
    best, arg, degen = float("inf"), np.zeros(3), 0
    for i0 in range(grid):
        for i1 in range(grid):
            for i2 in range(grid):
                k = 2 * np.pi * np.array([(i0 + offset) / grid,
                                          (i1 + 2 * offset) / grid,
                                          (i2 + 3 * offset) / grid])
                r, nb = _bloch_residual(k)
                if nb < BLOCH_B_FLOOR:
                    degen += 1
                    continue
                if np.isfinite(r) and r < best:
                    best, arg = r, k
    return dict(min=best, argmin=arg, degenerate=degen)


# ==========================================================================
# Y2  THE TIED ORBIT
# ==========================================================================

def y2_tied_orbit(y0):
    print()
    print("=" * 78)
    print("Y2  THE TIED ORBIT, and which half of Danskin's condition holds")
    print("=" * 78)
    print("  The recorded rule: smoothness breaks when a PROPER SUBSET of a tied")
    print("  orbit is selected AND its members have different derivatives. Both")
    print("  halves are measured, because a tie with EQUAL derivatives is smooth")
    print("  and would give a gradient, not a corner.")
    out = {}
    topos = {t.name: t for t in build_topologies()}
    tt = topos["CUBE8-M (60/60/90 basis)"]
    mem, _ = _classes(tt)
    aico, achord = y0["aico"], y0["achord"]

    for tag, a, kind in (("intra", aico, "intra"), ("inter", achord, "inter")):
        pr, _, _ = active_set(tt, mem, a, kind)
        if not pr:
            out[tag] = dict(n=0, ngrad=0, spread=float("nan"),
                            within=float("nan"))
            print(f"  {tag}: no active members")
            continue
        gz, gu = span_rate_rows(tt, pr, a)
        # distinct gradient directions among the tied members
        g = np.hstack([gz, gu])
        uniq = []
        for row in g:
            if not any(np.allclose(row, w, atol=TOL["tie"], rtol=0.0)
                       for w in uniq):
                uniq.append(row)
        # spread WITHIN one unit (intra only): the six diagonals of unit 0
        within = float("nan")
        if kind == "intra":
            same = [g[r] for r, ((ip, _), (iq, _)) in enumerate(pr)
                    if ip == 0 and iq == 0]
            if len(same) > 1:
                arr = np.array(same)
                within = float(np.abs(arr - arr[0]).max())
        out[tag] = dict(n=len(pr), ngrad=len(uniq),
                        spread=float(np.abs(g).max()), within=within)
        print(f"  {tag:5s} at a = {a:.9f}: {len(pr)} tied members, "
              f"{len(uniq)} distinct gradients")
        if kind == "intra":
            print(f"        the six diagonals inside ONE unit agree to "
                  f"{within:.2e} -- the tie WITHIN a unit is smooth")
            print(f"        and the corner comes from the tie ACROSS units, which")
            print(f"        is exactly {len(uniq)} = N distinct gradients.")
    # THE TIE TOLERANCE, BOUNDED FROM ABOVE BY MEASUREMENT. Members inside one
    # unit agree to `within`; members in different units differ by at least
    # `across`. TOL["tie"] must sit strictly between, or "distinct gradients"
    # counts whatever the constant happens to allow.
    pr_i, _, _ = active_set(tt, mem, aico, "intra")
    across = float("inf")
    if pr_i:
        gz_i, gu_i = span_rate_rows(tt, pr_i, aico)
        g_i = np.hstack([gz_i, gu_i])
        for r1, ((ip, _), _) in enumerate(pr_i):
            for r2, ((jp_, _), _) in enumerate(pr_i):
                if ip < jp_:
                    across = min(across,
                                 float(np.abs(g_i[r1] - g_i[r2]).max()))
    out["tie_across"] = across
    out["tie_band_ok"] = (out.get("intra", {}).get("within", float("inf"))
                          < TOL["tie"] < across)
    print(f"  TOL[tie] band: within-unit {out.get('intra', {}).get('within', float('nan')):.1e}"
          f" < {TOL['tie']:.0e} < across-unit {across:.1e}")
    out["intra_ngrad_is_n"] = out.get("intra", {}).get("ngrad", -1) == tt.n
    out["intra_within_tied"] = (np.isfinite(out.get("intra", {}).get("within",
                                                                    np.nan))
                                and out["intra"]["within"] < TOL["tie"])
    out["n_units"] = tt.n
    return out


# ==========================================================================
# Y3  THE ORDER OF RELIEF -- the bead's question
# ==========================================================================

def y3_order(y0):
    print()
    print("=" * 78)
    print("Y3  ORDER OF RELIEF of the binding span's excess under dephasing")
    print("=" * 78)
    aico, achord = y0["aico"], y0["achord"]
    topos = {t.name: t for t in build_topologies()}
    names = ["N2 (one contact)", "CHAIN5", "SQUARE4 (one 4-cycle)",
             "SC7 star (six-around-one)", "CUBE8-M (60/60/90 basis)",
             "CUBE8-R (60/60/60 basis)", "CUBE27-M"]
    out = {"rows": {}}

    print()
    print("  (a) the EXACT first-order optimum. For each cluster and each")
    print("      candidate binder, the minimum over the traceless sup-sphere of")
    print("      the one-sided directional derivative of the max excess, with")
    print("      the placement rates free (mechanisms included) and the contacts")
    print("      held closed. NEGATIVE means relief at first order.")
    print(f"      {'topology':30s} {'binder':6s} {'m':>4s} {'min D+':>16s} "
          f"{'verdict':>12s}")
    for nm in names:
        t = topos.get(nm)
        if t is None:
            continue
        mem, _ = _classes(t)
        for tag, a, kind in (("intra", aico, "intra"), ("inter", achord, "inter")):
            pr, _, _ = active_set(t, mem, a, kind)
            if not pr:
                out["rows"][(nm, tag)] = dict(n=0, val=float("nan"), bad=0)
                print(f"      {nm:30s} {tag:6s} {0:4d} {'no member':>16s} "
                      f"{'-':>12s}")
                continue
            val, arg, bad, unb = min_directional_derivative(t, pr, a)
            if unb:
                verdict = "UNBOUNDED"
            elif val < -TOL["lp"]:
                verdict = "RELIEVED"
            elif val > TOL["lp"]:
                verdict = "worsened"
            else:
                verdict = "flat"
            out["rows"][(nm, tag)] = dict(n=len(pr), val=val, bad=bad, unb=unb,
                                          selfstress=admissible_dephasing(
                                              t, a)["selfstress"],
                                          arg=None if arg is None
                                          else np.asarray(arg).copy())
            print(f"      {nm:30s} {tag:6s} {len(pr):4d} {val:+16.9e} "
                  f"{verdict:>12s}")
    print("      UNBOUNDED is a RESULT, not a failure: the programme runs to")
    print("      minus infinity because the placement rates are unbounded and a")
    print("      MECHANISM relieves the span, so the span was never held by the")
    print("      linkage. Those rows are the ones (c) explains, and they are the")
    print("      correction this file makes to jb_x.")

    print()
    print("  (b) the intra-unit optimum, DERIVED and then compared. A folding")
    print("      diagonal is a function of ONE unit's phase, so its rate is")
    print("      s'(a) u_i with the same s' for every unit and every diagonal;")
    print("      the max over units of s' u_i, minimised over traceless u with")
    print("      |u|_inf = 1, is |s'| / (N-1), attained by ONE unit contracting")
    print("      while the other N-1 expand by 1/(N-1) each.")
    sp, spdev = _diag_slope(aico)
    out["diag_slope"] = sp
    out["diag_slope_dev"] = spdev
    # BOUNDED FROM BOTH SIDES. An upper bound alone is satisfied by a routine
    # that reports its own deviation as zero -- and a probe that did exactly
    # that reddened nothing. Six central differences over six geometrically
    # distinct vertex pairs agreeing to the LAST BIT is not a measurement, it is
    # a routine not measuring.
    out["diag_slope_agree"] = 0.0 < spdev < 1e-8 and abs(sp) > 1e-3
    print(f"      s'(a_ico) measured on all six diagonals: {sp:.9f} per degree")
    print(f"      the six agree to {spdev:.2e}, which is the central")
    print(f"      difference's own truncation level and NOT zero. That number is")
    print("      gated: 'the same on all six diagonals' was printed and ungated")
    print("      until an audit found the function averaging rather than")
    print("      asserting. The magnitude of s' is gated too, so a routine")
    print("      returning six zeros cannot pass by agreeing perfectly.")
    der = {}
    for nm in names:
        t = topos.get(nm)
        if t is None or t.n < 2:
            continue
        r = out["rows"].get((nm, "intra"))
        if not r or r["n"] == 0 or not np.isfinite(r["val"]):
            continue
        pred = abs(sp) / (t.n - 1)
        adm = admissible_dephasing(t, aico)
        der[nm] = (r["val"], pred, adm["selfstress"])
        flag = "==" if abs(r["val"] - pred) < 1e-9 else ">"
        print(f"      {nm:30s} measured {r['val']:.9e}  {flag} derived "
              f"{pred:.9e}   selfstress {adm['selfstress']}")
    out["derived"] = der
    # THE INEQUALITY ARM IS NOW STRICT AND ITS MARGIN IS BOUNDED FROM BOTH
    # SIDES. It used to read `v > p - 1e-12`, which ACCEPTS v == p -- and an
    # enumerator that ignored the obstruction returns exactly v == p, so the arm
    # the prose says catches that enumerator could not catch it. A probe
    # replacing the whole arm with `True` reddened nothing. The margin must now
    # be cleared strictly, and the gate bounds DERIVED_MARGIN from above by the
    # separation actually measured, so it cannot be raised until it excludes the
    # truth either.
    out["derived_ok"] = bool(der) and all(
        (abs(v - p) < 1e-9) if ss == 0 else (v > p + DERIVED_MARGIN)
        for v, p, ss in der.values())
    seps = [v - p for v, p, ss in der.values() if ss > 0]
    out["derived_sep"] = min(seps) if seps else float("inf")
    out["derived_margin_ok"] = (TOL["lp"] < DERIVED_MARGIN
                                < out["derived_sep"])
    out["derived_n_stressed"] = len(seps)
    out["derived_n_free"] = sum(1 for v, p, ss in der.values() if ss == 0)
    print("      Equality holds exactly where the cluster has NO self-stress and")
    print("      STRICT inequality where it has some, because the obstruction")
    print("      removes the cheapest pattern. That is a two-sided check: an")
    print("      enumerator that ignored the obstruction would return exactly")
    print("      the derived value and fail the second half, and one that")
    print("      ignored the derivation would fail the first.")
    print(f"      the strict margin is {DERIVED_MARGIN:.0e}, inside its own band")
    print(f"      {TOL['lp']:.0e} < {DERIVED_MARGIN:.0e} < "
          f"{out['derived_sep']:.3e} (the measured separation), and BOTH")
    print(f"      classes are occupied: {out['derived_n_free']} clusters with no")
    print(f"      self-stress, {out['derived_n_stressed']} with some.")

    print()
    print("  (c) the MECHANISM-ONLY row: can the linkage relieve the span with")
    print("      NO dephasing at all, by moving the units on their joints? This")
    print("      is the correction jb_x's placement-conditional span enumeration")
    print("      needs, and for the star it changes the answer.")
    mech = {}
    for nm in names:
        t = topos.get(nm)
        if t is None:
            continue
        mem, _ = _classes(t)
        for tag, a, kind in (("intra", aico, "intra"), ("inter", achord, "inter")):
            pr, _, _ = active_set(t, mem, a, kind)
            if not pr:
                continue
            gz, _ = span_rate_rows(t, pr, a)
            gzmax = float(np.abs(gz).max())
            val, bad = mechanism_only_relief(t, pr, a)
            mech[(nm, tag)] = (val, gzmax, bad)
            note = ("STRUCTURAL ZERO (intra spans do not see placement)"
                    if gzmax < 1e-14 else
                    ("relieved by mechanisms alone" if val < -1e-9
                     else "mechanisms cannot relieve"))
            print(f"      {nm:30s} {tag:6s} {val:+13.6e}  |Gz|max "
                  f"{gzmax:.2e}  {note}")
    out["mech"] = mech
    out["intra_gz_zero"] = all(g < 1e-14 for (_, tg), (_, g, _) in mech.items()
                               if tg == "intra")
    out["star_inter_mech"] = mech.get(("SC7 star (six-around-one)", "inter"),
                                      (float("nan"), 0.0, 0))[0]

    print()
    print("  (d) THE CORNER. For a spread of transverse directions, the one-sided")
    print("      derivative in BOTH +u and -u. If both are positive the two-sided")
    print("      derivative cannot exist, and 'the derivative is zero' and 'the")
    print("      derivative does not exist' are then DIFFERENT ANSWERS and this")
    print("      is the second one.")
    tt = topos["CUBE8-M (60/60/90 basis)"]
    mem, _ = _classes(tt)
    rng = np.random.default_rng(SPREAD_SEED)
    dirs = []
    for _ in range(N_SPREAD):
        u = rng.standard_normal(tt.n)
        u -= u.mean()
        mx = float(np.abs(u).max())
        if mx > 0:
            dirs.append(u / mx)
    corner = {}
    for tag, a, kind in (("intra", aico, "intra"), ("inter", achord, "inter")):
        pr, _, _ = active_set(tt, mem, a, kind)
        if not pr:
            continue
        vals = []
        for u in dirs:
            vp, sp1 = directional_derivative(tt, pr, a, u)
            vm, sp2 = directional_derivative(tt, pr, a, -u)
            if sp1 or sp2 or not (np.isfinite(vp) and np.isfinite(vm)):
                continue
            vals.append((vp, vm))
        if not vals:
            corner[tag] = dict(n=0)
            continue
        arr = np.array(vals)
        both_pos = int(np.count_nonzero((arr[:, 0] > 0) & (arr[:, 1] > 0)))
        anti = float(np.abs(arr[:, 0] + arr[:, 1]).min())
        corner[tag] = dict(n=len(vals), both_pos=both_pos,
                           min_plus=float(arr[:, 0].min()),
                           min_minus=float(arr[:, 1].min()),
                           min_sum=anti,
                           max_sum=float(np.abs(arr[:, 0] + arr[:, 1]).max()))
        print(f"      {tag:6s}: {len(vals)} directions, both-sides-positive in "
              f"{both_pos}; smallest D+(u) {arr[:, 0].min():+.6e},")
        print(f"              smallest D+(-u) {arr[:, 1].min():+.6e}; "
              f"min |D+(u)+D+(-u)| = {anti:.6e}")
        print("              (an odd function would give D+(u)+D+(-u) = 0 for")
        print("               every direction, so a small minimum here is the")
        print("               falsifier for 'this is really a gradient'.)")
    out["corner"] = corner

    print()
    print("  (e) THE ORDER, by amplitude ladder. E(eps u) along the ladder,")
    print("      geometric with an irrational ratio so no rung is commensurate")
    print("      with another. A first-order excess gives E/eps CONSTANT; a")
    print("      second-order one gives E/eps proportional to eps.")
    order = {}
    for tag, a, kind in (("intra", aico, "intra"), ("inter", achord, "inter")):
        pr, _, _ = active_set(tt, mem, a, kind)
        if not pr:
            continue
        base = out["rows"].get(("CUBE8-M (60/60/90 basis)", tag), {}).get("arg")
        cand = [("LP argmin", base)] if base is not None else []
        # GUARDED. `dirs` is built from a seeded generator and its length is
        # N_SPREAD, but indexing it unconditionally makes the whole file depend
        # on a constant a mutation probe can reach -- and a probe that set
        # N_SPREAD to 1 produced an IndexError that destroyed the verdict table
        # instead of a red row. This is the sentinel shape the directory's
        # post-mortems name: a list assigned under one condition and subscripted
        # under none.
        for idx in (0, 1):
            if len(dirs) > idx:
                cand.append((f"seeded dir {idx}", dirs[idx]))
        for label, u in cand:
            if u is None:
                continue
            res = _amplitude_ladder(tt, pr, a, np.asarray(u, float),
                                    EPS_TOP, EPS_RATIO, EPS_RUNGS)
            res2 = _amplitude_ladder(tt, pr, a, np.asarray(u, float),
                                     EPS_TOP_ALT, EPS_RATIO_ALT, EPS_RUNGS)
            order[(tag, label)] = (res, res2)
            print(f"      {tag:6s} {label:14s} slope {res['slope']:.9f} over "
                  f"{res['nfit']} rungs, |E|/eps spread {res['spread']:.3e}, "
                  f"sign {res['sign']:+.0f}")
            print(f"      {'':6s} {'second ladder':14s} slope "
                  f"{res2['slope']:.9f} over {res2['nfit']} rungs, spread "
                  f"{res2['spread']:.3e}, sign {res2['sign']:+.0f}")
            print(f"      {'':6s} {'closure':14s} worst Gauss-Newton residual "
                  f"down both ladders {max(res['res'], res2['res']):.2e}")
    out["order"] = order
    # THREE QUANTITIES THAT WERE COMPUTED AND NEVER READ, NOW GATED. Each backs
    # a sentence the verdict prints. `one_sign` backs "the sign is constant down
    # every ladder"; `nfit` backs "over N rungs"; and `res` -- the worst
    # Gauss-Newton closure residual down a ladder -- backs the whole ladder,
    # because an excess read off a configuration that did not close is not a
    # measurement of anything. That last one turned out to be hiding a live
    # solver defect the moment it was read: see LSTSQ_RCONDS.
    out["order_one_sign"] = bool(order) and all(
        r1["one_sign"] and r2["one_sign"] for r1, r2 in order.values())
    out["order_nfit_min"] = min((min(r1["nfit"], r2["nfit"])
                                 for r1, r2 in order.values()), default=0)
    out["order_res_max"] = max((max(r1["res"], r2["res"])
                                for r1, r2 in order.values()), default=float("inf"))
    out["order_n"] = len(order)
    # THE TWO LADDERS MUST SHARE NO RUNG. "Measured on two ladders" is worth
    # nothing if the second is the first under another name, and no gate row
    # can see that from the results alone -- the spreads would simply agree.
    # Measured directly on the rung sets, in log space so the comparison is
    # scale-free.
    l1 = np.log(np.array([EPS_TOP * EPS_RATIO ** j for j in range(EPS_RUNGS)]))
    l2 = np.log(np.array([EPS_TOP_ALT * EPS_RATIO_ALT ** j
                          for j in range(EPS_RUNGS)]))
    out["ladder_sep"] = float(np.abs(l1[:, None] - l2[None, :]).min())
    print(f"      the two ladders share no rung: closest pair differs by "
          f"{out['ladder_sep']:.4f} in log-amplitude")
    print("      The SIGN is constant down every ladder, which is part of the")
    print("      first-order claim rather than decoration: a term of even order")
    print("      would flip nothing but a term of higher odd order competing")
    print("      with the linear one would move the crossover into the window.")

    print()
    print("  (f) CONTROL on the admissibility constraint inside the programme.")
    print("      Y1 measures that a single-direction parity pattern is NOT an")
    print("      admissible dephasing of a 3x3x3 cluster and that an affine ramp")
    print("      IS. If the contact rows really bind inside the linear programme")
    print("      then the first must come back INFEASIBLE and the second must")
    print("      SOLVE. Without this row nothing checks that the equality block")
    print("      is doing any work at all -- every other Y3 row is satisfiable")
    print("      by a programme that dropped it.")
    t27 = topos.get("CUBE27-M")
    ctrl = {}
    if t27 is not None:
        m27, _ = _classes(t27)
        pr27, _, _ = active_set(t27, m27, achord, "inter")
        sites = np.array(t27.lattice_sites, float)
        pats = {"single-direction parity (INADMISSIBLE)":
                (-1.0) ** sites[:, 0],
                "affine ramp along g0 (ADMISSIBLE)": sites[:, 0].copy(),
                "zone-corner parity (ADMISSIBLE)":
                (-1.0) ** sites.sum(axis=1)}
        for lab, w in pats.items():
            w = w - w.mean()
            mx = float(np.abs(w).max())
            if mx < 1e-12 or not pr27:
                continue
            w = w / mx
            val, st = directional_derivative(t27, pr27, achord, w)
            ctrl[lab] = (val, st)
            print(f"      {lab:42s} status {st}  D+ {val:+.6e}")
    out["adm_ctrl"] = ctrl
    out["adm_ctrl_ok"] = (
        ctrl.get("single-direction parity (INADMISSIBLE)", (0, 0))[1] != 0
        and ctrl.get("affine ramp along g0 (ADMISSIBLE)", (0, 1))[1] == 0
        and ctrl.get("zone-corner parity (ADMISSIBLE)", (0, 1))[1] == 0)

    print()
    print("  (g) THE FINITE EXCURSION, because the bead's Outcome 2 is a FINITE")
    print("      phase difference and the ladder above is not. The ladder fits")
    print(f"      |E| in [{LADDER_FLOOR:.0e}, {LADDER_CEIL:.0e}], which for this")
    print("      excess is a dephasing of order 1e-7 to 1e-4 DEGREES -- three to")
    print("      six decades below the question. A strictly positive FIRST-ORDER")
    print("      slope does not by itself preclude the excess turning around at")
    print("      finite amplitude, so 'no threshold' does not follow from 'first")
    print("      order everywhere' and must be measured where it could fail.")
    print("      CUBE27-M is that place: the one cluster here where relief fails")
    print("      at first order is the only one where a threshold could live.")
    print("      Rungs run to FOUR DEGREES on the exact nonlinear model and")
    print("      EVERY RUNG'S CLOSURE RESIDUAL IS GATED. An excess read off a")
    print("      configuration that did not close is not a measurement -- and")
    print("      that is not hypothetical here: two other directions tried")
    print("      during this work returned finite excesses at residuals up to")
    print("      6e-01, non-monotone and meaningless, and the control below")
    print("      exhibits one of them so the closure gate has a demonstrated")
    print("      failure mode rather than a claimed one.")
    exc = {}
    if t27 is not None:
        m27, _ = _classes(t27)
        sites = np.array(t27.lattice_sites, float)
        chi = (-1.0) ** sites.sum(axis=1)
        u = chi - chi.mean()
        mx = float(np.abs(u).max())
        u = u / mx if mx > 1e-12 else u
        for tag, a, kind in (("intra", aico, "intra"),
                             ("inter", achord, "inter")):
            prx, _, _ = active_set(t27, m27, a, kind)
            if not prx:
                continue
            vals, vals_min, res = [], [], []
            for eps in EXCURSION_RUNGS:
                ph = a + eps * u
                z, r, _ = solve_dephased(t27, ph)
                res.append(r)
                # THE `max` HERE IS THE WORST-CASE SPAN AND IS OBJECT-LEVEL
                # LOAD-BEARING, not decoration a probe corpus happened to
                # never reach: `min` in its place (best-case span among the
                # SAME active members, SAME solved configuration) is a
                # falsifier for `exc_positive` below that mutates the object
                # rather than TOL["lp"]. It is installed as a gated CONTROL,
                # not narrated -- `exc_falsifier_ok` and its row below, which
                # follow this section's own `exc_ctrl_res` idiom. Computed
                # here, in the same loop and from the same `z`, so it costs
                # no extra solve.
                lengths = [member_length(ph, z, mm, t27) - STRUT_LEN
                          for mm in prx]
                vals.append(max(lengths))
                vals_min.append(min(lengths))
            arr = np.array(vals, float)
            arr_min = np.array(vals_min, float)
            exc[tag] = dict(vals=arr, vals_min=arr_min,
                            res=float(max(res, default=1.0)), n=len(vals))
            print(f"      {tag:6s} on CUBE27-M along the zone-corner "
                  f"alternation, |u|_inf = 1:")
            print("        eps(deg) " + " ".join(f"{e:>11g}"
                                                 for e in EXCURSION_RUNGS))
            print("        excess   " + " ".join(f"{v:>+11.4e}" for v in arr))
            print(f"        worst closure residual over the {len(vals)} rungs: "
                  f"{max(res, default=1.0):.2e}")
    out["excursion"] = exc
    # THE COMPARISON ITSELF LIVES IN THE GATE, not here. A boolean precomputed
    # in a section hides its own threshold from the verdict table and gives a
    # probe a second place to loosen it; the gate row below reads the NUMBER and
    # compares it against SOLVE_TOL, which is already reached by two probes.
    out["exc_res_max"] = max((v["res"] for v in exc.values()),
                             default=float("inf"))
    out["exc_positive"] = bool(exc) and all(float(v["vals"].min()) > TOL["lp"]
                                            for v in exc.values())
    out["exc_monotone"] = bool(exc) and all(
        bool(np.all(np.diff(v["vals"]) > 0.0)) for v in exc.values())
    out["exc_top"] = max(EXCURSION_RUNGS) if EXCURSION_RUNGS else 0.0
    out["exc_n"] = min((v["n"] for v in exc.values()), default=0)
    # THE FALSIFIER FOR `exc_positive`, INSTALLED AS A GATED CONTROL rather
    # than narrated. `vals_min` (collected above, same loop, same solve) picks
    # the BEST-case span among the same active members instead of the
    # worst-case one `vals` reads -- the same object the real measurement
    # reads, aggregated the other way. If it comes out NEGATIVE at every
    # rung, `exc_positive` is demonstrated to be sensitive to which object the
    # geometry hands it, rather than to a threshold a probe could just move.
    out["exc_falsifier_max"] = (
        max(float(v["vals_min"].max()) for v in exc.values())
        if exc else float("inf"))
    out["exc_falsifier_ok"] = bool(exc) and all(
        float(v["vals_min"].max()) < -TOL["lp"] for v in exc.values())
    print(f"      CONTROL: the SAME solves, aggregated by MIN instead of MAX")
    print(f"      over the active set, go NEGATIVE at every rung -- least")
    print(f"      negative case {out['exc_falsifier_max']:+.4e}, versus "
          f"TOL[lp] {TOL['lp']:.0e}. Not")
    print("      all active-set members relieve together, so selecting the")
    print("      wrong one flips the sign the row above asserts: the")
    print("      measurement is object-sensitive, not a constant-detector.")
    # THE CONTROL: a direction that does NOT close, solved at one amplitude, so
    # that the closure gate above is demonstrated rather than asserted.
    out["exc_ctrl_res"] = 0.0
    if t27 is not None:
        sites = np.array(t27.lattice_sites, float)
        w = sites[:, 0] - sites[:, 0].mean()
        mw = float(np.abs(w).max())
        if mw > 1e-12:
            _, rr, _ = solve_dephased(
                t27, achord + EXCURSION_CTRL_EPS * (w / mw))
            out["exc_ctrl_res"] = rr
    print(f"      CONTROL: the affine ramp along g0 at eps = "
          f"{EXCURSION_CTRL_EPS:g} deg closes only to")
    print(f"      {out['exc_ctrl_res']:.2e}, far above SOLVE_TOL "
          f"{SOLVE_TOL:.0e}. It is a first-order")
    print("      admissible direction that does NOT integrate to a finite")
    print("      motion at this size, which is a scope statement about Y1's")
    print("      affine patterns and is reported rather than smoothed over.")
    print("      SO: along the direction that DOES close, the excess is")
    print("      strictly positive and strictly increasing all the way to four")
    print("      degrees. NO TURNING POINT, hence no threshold, on that")
    print("      direction. This is NOT a statement about every direction at")
    print("      finite amplitude, and the file does not make one.")
    return out


def _diag_slope(a, h=1e-6):
    """(mean d(folding diagonal)/da, largest deviation from that mean).

    THE SECOND NUMBER IS THE POINT. The docstring here used to promise "the
    assertion that all six agree" and the function averaged the six and asserted
    nothing -- so "the same on all six diagonals", which the verdict prints, was
    never gated. It is now: the deviation comes back and the gate reads it. It
    sits at the central difference's own truncation level, not at zero, which is
    why it is reported as a number rather than compared to zero.
    """
    if not DIAGONALS:
        return float("nan"), float("inf")
    vals = []
    for (k, l) in DIAGONALS:
        vp = verts(a + h)
        vm = verts(a - h)
        vals.append((float(np.linalg.norm(vp[k] - vp[l]))
                     - float(np.linalg.norm(vm[k] - vm[l]))) / (2 * h))
    arr = np.array(vals, float)
    return float(arr.mean()), float(np.abs(arr - arr.mean()).max())


def _amplitude_ladder(topo, pairs, a, u, top, ratio, rungs):
    """|E(eps u)| down an irrational-ratio ladder, on the EXACT nonlinear model.

    THE ABSOLUTE VALUE MATTERS. Along a RELIEVING direction E is negative, and a
    fit restricted to positive E has nothing to fit -- an earlier version
    returned nan for exactly the direction the whole file is about. The sign is
    reported separately and is constant down the ladder, which is itself part of
    the first-order claim.

    The fit window is bounded at BOTH ends by absolute, stated numbers:
    |E| must clear LADDER_FLOOR, which is where the Gauss-Newton solve's own
    residual starts to show, and must stay under LADDER_CEIL, above which the
    nonlinear terms are visible and a first-order fit would be measuring
    something it does not claim. A window open at the top would let curvature at
    eps = 1 pollute a statement about eps -> 0; a window open at the bottom would
    fit the solver's noise.
    """
    eps, e_vals, res_max = [], [], 0.0
    for j in range(rungs):
        s = top * ratio ** j
        ph = a + s * u
        z, r, _ = solve_dephased(topo, ph)
        res_max = max(res_max, r)
        best = -np.inf
        for mm in pairs:
            best = max(best, member_length(ph, z, mm, topo) - STRUT_LEN)
        eps.append(s)
        e_vals.append(best)
    eps = np.array(eps)
    e_vals = np.array(e_vals)
    mag = np.abs(e_vals)
    good = np.isfinite(mag) & (mag > LADDER_FLOOR) & (mag < LADDER_CEIL)
    signs = np.sign(e_vals[good]) if np.any(good) else np.zeros(0)
    one_sign = bool(signs.size) and bool(np.all(signs == signs[0]))
    if int(np.count_nonzero(good)) < 4:
        return dict(slope=float("nan"), spread=float("nan"), nfit=0,
                    res=res_max, sign=0.0, one_sign=False)
    lx = np.log(eps[good])
    ly = np.log(mag[good])
    slope = float(np.polyfit(lx, ly, 1)[0])
    rv = mag[good] / eps[good]
    spread = float((rv.max() - rv.min()) / max(abs(rv.mean()), 1e-300))
    return dict(slope=slope, spread=spread, nfit=int(np.count_nonzero(good)),
                res=res_max, sign=float(signs[0]), one_sign=one_sign)


# ==========================================================================
# Y4  THE STATIC QUESTION AND THE PATH QUESTION
# ==========================================================================

def y4_static_vs_path(y0):
    print()
    print("=" * 78)
    print("Y4  STATIC versus PATH: two questions the bead says may disagree")
    print("=" * 78)
    print("  STATIC: can the array EXPAND -- can the mean phase decrease -- once")
    print("          it is dephased?")
    print("  PATH  : can units TAKE TURNS -- is there a motion at FIXED mean")
    print("          phase in which some units expand while others contract?")
    print("  Both are first-order linear programmes over the same cone, so they")
    print("  are comparable rather than merely analogous.")
    out = {}
    aico, achord = y0["aico"], y0["achord"]
    topos = {t.name: t for t in build_topologies()}
    names = ["CUBE8-M (60/60/90 basis)", "CUBE27-M", "SC7 star (six-around-one)"]
    print()
    print("  BOTH ARE FIRST-ORDER PROGRAMMES AT THE SAME POINT, so neither can")
    print("  distinguish a static relaxation from a path in a stronger sense")
    print("  than the tangent cone allows. Where the answer below is 'empty',")
    print("  it is EMPTY AT FIRST ORDER, and that qualifier is not droppable.")
    print()
    print(f"      {'topology':30s} {'binder':6s} {'STATIC min mean(u)':>20s} "
          f"{'PATH min u_i':>14s} {'':>6s} {'agree':>7s}")
    rows = {}
    for nm in names:
        t = topos.get(nm)
        if t is None:
            continue
        mem, _ = _classes(t)
        for tag, a, kind in (("intra", aico, "intra"), ("inter", achord, "inter")):
            pr, _, _ = active_set(t, mem, a, kind)
            if not pr:
                continue
            st = _static_lp(t, pr, a)
            pa = _path_lp(t, pr, a)
            agree = (st < -TOL["lp"]) == (pa < -TOL["lp"])
            # MARK THE VARIABLE BOUND. `-1.000000e+00` in the PATH column is
            # the LOWER BOUND on u_i, not a computed optimum interior to the
            # feasible set -- the programme wanted to go further and the box
            # stopped it. That exact number is the signature of the vacuous-
            # slack bug this file's own history records, and printing it
            # unmarked next to a genuine optimum invites the reader to make
            # that mistake a second time.
            atb = "BOUND" if abs(abs(pa) - 1.0) < 1e-12 else ""
            rows[(nm, tag)] = (st, pa, agree, bool(atb))
            print(f"      {nm:30s} {tag:6s} {st:+20.9e} {pa:+14.6e} "
                  f"{atb:>6s} {str(agree):>7s}")
    out["rows"] = rows
    out["all_agree"] = bool(rows) and all(v[2] for v in rows.values())
    out["n_disagree"] = sum(1 for v in rows.values() if not v[2])
    out["n_at_bound"] = sum(1 for v in rows.values() if v[3])
    print("      BOUND marks a PATH value sitting on the variable bound |u_i|")
    print("      <= 1 rather than at an interior optimum: the programme is")
    print(f"      reporting that it ran out of box, on {out['n_at_bound']} of "
          f"{len(rows)} rows.")
    # The INTRA rows have a DERIVATION: the constraint is s'(a) u_i <= 0 with
    # s' < 0, hence u_i >= 0 for every unit; the mean is then >= 0 with equality
    # only at u = 0, and the smallest single coordinate is 0. So both
    # programmes must return exactly zero, and that is what is gated.
    iv = [(st, pa) for (_, tg), (st, pa, _, _) in rows.items()
          if tg == "intra"]
    out["intra_worst"] = (max(max(abs(a), abs(b)) for a, b in iv) if iv
                          else float("inf"))
    out["intra_zero_ok"] = bool(iv) and out["intra_worst"] < 1e-9
    xv = [min(st, pa) for (_, tg), (st, pa, _, _) in rows.items()
          if tg == "inter" and np.isfinite(st) and np.isfinite(pa)]
    out["inter_min"] = min(xv) if xv else float("inf")
    out["inter_neg"] = bool(xv) and out["inter_min"] < -1e-9
    print()
    if out["all_agree"]:
        print("  MEASURED: the two questions AGREE on every row, so the")
        print("  distinction the bead proposed is EMPTY AT FIRST ORDER -- and")
        print("  the qualifier is the whole content of the sentence, because")
        print("  both programmes are LPs over the same tangent cone at the same")
        print("  point and could not have distinguished a relaxation from a")
        print("  path even if the numbers had come out differently. What is")
        print("  refuted is the bead's FIRST-ORDER form of the distinction.")
        print("  The reason is structural rather than numerical: both programmes ask")
        print("  whether the traceless cone {all member rates <= 0} is more than")
        print("  the origin, and the static programme's extra freedom -- letting")
        print("  the mean fall -- is worth nothing, because a falling mean")
        print("  lengthens every active span by construction. That is a")
        print("  REFUTATION of the proposed resolution, not a confirmation.")
    else:
        print(f"  MEASURED: the two questions DISAGREE on {out['n_disagree']} of")
        print(f"  {len(rows)} rows, so the distinction the bead proposed is REAL.")
        print("  Where they differ, the static programme finds no net expansion")
        print("  while the path programme finds a unit free to expand at fixed")
        print("  mean phase: the array can redistribute what it cannot gain.")
    print("  Either way the sentence above is COMPUTED from the table, not")
    print("  written in advance -- an earlier version asserted agreement in")
    print("  prose while the table it sat under said otherwise.")

    print()
    print("  (b) WHAT IS NOT EMPTY, and what resolves the tension instead. The")
    print("      lock is measure zero AS A SET and unrelieved AS A FUNCTIONAL.")
    print("      Those are different objects and both are true.")
    print("      Below: the feasible dephasing polytope at mean-phase slack eta,")
    print("      for the intra-unit binder, where feasibility is exactly")
    print("      a_i >= a_lock for every unit. Its sup-radius is eta and it")
    print("      COLLAPSES TO A POINT at eta = 0 -- so the fully locked state is")
    print("      one point of an (N-1)-dimensional dephasing space.")
    tt = topos["CUBE8-M (60/60/90 basis)"]
    mem, _ = _classes(tt)
    pr, _, _ = active_set(tt, mem, aico, "intra")
    radii = {}
    for eta in (0.0, 1e-6, 1e-3, 1e-1, 1.0):
        radii[eta] = _feasible_radius(tt, aico, eta)
    out["radii"] = radii
    out["radius_linear"] = all(
        abs(radii[e] - e) < 1e-6 * max(e, 1e-6) for e in (1e-3, 1e-1, 1.0))
    out["radius_zero_at_zero"] = radii[0.0] < 1e-12
    for eta, r in radii.items():
        print(f"        eta = {eta:<10g} sup-radius of the feasible dephasing "
              f"set = {r:.9g}")
    print("      The radius is eta exactly. So the arrest is total only when")
    print("      every unit arrives at the wall at once, and the set of ways to")
    print("      do that is a single point.")
    print()
    print("      WHAT THAT IS AND IS NOT, separated on purpose. WHAT IS")
    print("      MEASURED: on CUBE8-M, for the intra-unit binder, the feasible")
    print("      set's sup-radius equals the slack to fourteen decimals and")
    print("      collapses at zero slack. WHAT IS NEAR-TAUTOLOGOUS: feasibility")
    print("      for that binder is a_i >= a_lock for every i, so the feasible")
    print("      set is the positive orthant intersected with the mean-zero")
    print("      hyperplane, which meets it only at the origin -- radius(eta) =")
    print("      eta is that fact, and the bisection is confirming a")
    print("      derivation. WHAT IS NOT MEASURED HERE: the same radius for the")
    print("      INTER-unit binder, and the same radius on a cluster WITH")
    print("      self-stress. CUBE8-M is the one cluster in this file with")
    print("      none, and on CUBE27-M the worst pattern this bisection assumes")
    print("      need not be an admissible dephasing at all -- the Y1")
    print("      obstruction is not applied here. The set statement is")
    print("      therefore scoped to the intra binder on an unobstructed")
    print("      cluster.")
    print("      AND WHAT IS INTERPRETATION, NOT MEASUREMENT: that the owner's")
    print("      observation is about the SET rather than the FUNCTIONAL. The")
    print("      two objects genuinely are different and neither refutes the")
    print("      other -- that much is established. But he said motion")
    print("      PROGRESSES given a fractional difference, and 'does motion")
    print("      progress' is the FUNCTIONAL/PATH side, where this file's intra")
    print("      answer is that it does not. Assigning his observation to the")
    print("      set half is a reading, offered as one, and it is not what")
    print("      makes both claims true. What makes both true is that they are")
    print("      claims about different objects.")

    print()
    print("  (c) TAKING TURNS, made quantitative. At slack eta a travelling")
    print("      pattern of amplitude eta is feasible and its trough units sit")
    print("      exactly at the wall while the rest have slack. How many units")
    print("      are at the wall AT ONE INSTANT depends on whether the pattern")
    print("      is commensurate with the lattice: a commensurate wave arrests a")
    print("      whole plane of units together, an incommensurate one arrests")
    print("      them one at a time. The grid discipline this directory learned")
    print("      three times over is a physical statement here, not a method")
    print("      rule.")
    turns, tlo, thi = _taking_turns()
    out["turns"] = turns
    for lab, (kk, nwall) in turns.items():
        print(f"        {lab:34s} k = {kk:.6f}, units at the wall: {nwall}")
    out["turns_ok"] = (turns["commensurate  k = 2pi/4"][1] > 1
                       and turns["incommensurate k = 2pi/phi^2"][1] == 1)
    out["turns_lo"] = tlo
    out["turns_hi"] = thi
    # STRICT ON BOTH SIDES, matching the other three two-sided bands in this
    # file (tol_ok, rtol_ok, tie_band_ok all use `<` on both ends, not `<=`).
    # Verified this does not move the shipped verdict: `tlo` is bit-exact 0.0
    # (the trough units sit at exactly the minimum in double precision, so
    # `tied.max()` in `_taking_turns` can only ever be 0.0), and TURNS_TOL is
    # 1e-9 -- nowhere near that boundary. Aligning to `<` matches the band's
    # own docstring ("TURNS_TOL must EXCEED the largest residual") and the
    # criterion column's "0.0e+00..3.1e-04" reading, which is inclusive of
    # neither endpoint.
    out["turns_band_ok"] = tlo < TURNS_TOL < thi
    print(f"      TURNS_TOL band, measured: units counted as tied agree to "
          f"{tlo:.1e}")
    print(f"      and the nearest unit that is NOT tied is {thi:.3e} away, so "
          f"{TURNS_TOL:.0e}")
    print("      sits inside a band computed without reference to it. 'Exactly")
    print("      one unit' is therefore a measurement.")
    print("      SCOPE, stated because the sentence above reads bigger than it")
    print("      is: this is cos(k s) on 64 integers -- the equidistribution of")
    print("      an irrational rotation. No jitterbug, no contact and no")
    print("      Jacobian enters it. It says what a travelling phase pattern")
    print("      does on a lattice, not what THIS linkage does.")
    return out


def _static_lp(topo, pairs, a):
    """min mean(u) subject to every active member rate <= 0, contacts closed.

    NOTE the traceless row is DROPPED here and only here -- this is the one
    programme that is allowed to move the mean, since that is its question.
    Zero means the array cannot expand at all; negative means it can.
    """
    jr, jp = blocks(topo, a)
    gz, gu = span_rate_rows(topo, pairs, a)
    n = topo.n
    nz = 6 * n
    nv = nz + n
    a_eq = np.hstack([jr, jp]) if jr.size else np.zeros((0, nv))
    c = np.zeros(nv)
    c[nz:] = 1.0 / n
    res = _safe_linprog(c, A_ub=np.hstack([gz, gu]), b_ub=np.zeros(gz.shape[0]),
                  A_eq=a_eq if a_eq.size else None,
                  b_eq=np.zeros(a_eq.shape[0]) if a_eq.size else None,
                  bounds=[(None, None)] * nz + [(-1.0, 1.0)] * n,
                  method="highs")
    if res.status == 3:
        return -float("inf")
    return float(res.fun) if res.status == 0 else float("inf")


def _path_lp(topo, pairs, a):
    """min over i and over the cone of u_i, at FIXED mean phase.

    Negative means some unit may expand while the mean is held, that is, units
    can take turns. The traceless row is present, so nothing here is bought by
    moving the mean.

    THE MEMBER ROWS ARE `rate <= 0`, WITH NO SLACK VARIABLE. An earlier version
    reused the min-max programme, which carries a free slack `s` and rows
    `rate <= s`; with `s` absent from the objective it runs to +infinity and
    every member row is satisfied vacuously. That version reported -1.0 -- the
    variable bound -- on every row of the table, including the intra-unit rows
    where the constraint is a product of half-lines and the true answer is 0.
    """
    jr, jp = blocks(topo, a)
    gz, gu = span_rate_rows(topo, pairs, a)
    n = topo.n
    nz = 6 * n
    nv = nz + n
    a_eq = [np.hstack([jr, jp])] if jr.size else []
    b_eq = [np.zeros(jr.shape[0])] if jr.size else []
    row = np.zeros(nv)
    row[nz:] = 1.0
    a_eq.append(row.reshape(1, -1))
    b_eq.append(np.zeros(1))
    best = float("inf")
    for i in range(n):
        c = np.zeros(nv)
        c[nz + i] = 1.0
        res = _safe_linprog(c, A_ub=np.hstack([gz, gu]), b_ub=np.zeros(gz.shape[0]),
                      A_eq=np.vstack(a_eq), b_eq=np.concatenate(b_eq),
                      bounds=[(None, None)] * nz + [(-1.0, 1.0)] * n,
                      method="highs")
        if res.status == 3:
            return -float("inf")
        if res.status == 0 and res.fun < best:
            best = float(res.fun)
    return best


def _feasible_radius(topo, a_lock, eta):
    """Largest r with EVERY traceless |delta|_inf <= r feasible at a_lock + eta.

    The worst traceless pattern of sup-norm r is the one that drives a single
    unit down by the full r; every other unit is then above the wall, so the
    excess is decided by that one unit and the radius is found by bisection on
    the EXACT nonlinear intra-unit excess. Bisection is bracketed above by an
    absolute 10 degrees so a runaway cannot loop, and the bracket is widened
    geometrically rather than assumed.
    """
    n = topo.n
    if n < 2:
        return float("inf")

    def worst(r):
        d = np.full(n, r / (n - 1))
        d[0] = -r
        ph = a_lock + eta + d
        best = -np.inf
        for i in range(n):
            vv = verts(ph[i])
            for (k, l) in DIAGONALS:
                best = max(best, float(np.linalg.norm(vv[k] - vv[l]))
                           - STRUT_LEN)
        return best

    if worst(0.0) > 1e-12:
        return 0.0
    lo, hi = 0.0, 1e-9
    while hi < 10.0 and worst(hi) <= 0.0:
        lo, hi = hi, hi * 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if worst(mid) <= 0.0:
            lo = mid
        else:
            hi = mid
    return lo


TURNS_TOL = 1e-9


def _taking_turns(n=64):
    """How many units of a chain sit at the trough of a dephasing wave, and the
    two-sided band that TURNS_TOL has to sit inside.

    The trough is the wall in a driven expansion, so this counts how many units
    arrest at the same instant. The docstring here used to claim the tie
    tolerance was "bounded from BOTH sides by the measured gap"; the gaps were
    sorted into a variable nothing read and the tolerance was a bare literal.
    They are returned now, and the band is a gate row: TURNS_TOL must exceed the
    largest residual among the units counted AS tied and fall under the smallest
    nonzero gap to a unit that is not, so 'exactly one unit' is a measurement.
    """
    s = np.arange(n, dtype=float)
    out = {}
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    lo, hi = 0.0, float("inf")
    for lab, k in (("commensurate  k = 2pi/4", 2 * np.pi / 4.0),
                   ("commensurate  k = 2pi/8", 2 * np.pi / 8.0),
                   ("incommensurate k = 2pi/phi^2", 2 * np.pi / (phi ** 2)),
                   ("incommensurate k = 2pi/pi", 2.0)):
        d = np.cos(k * s)
        gaps = np.sort(d - d.min())
        # THE SPLIT IS TOLERANCE-FREE: the units at the trough sit at EXACTLY
        # the minimum in double precision, so the two populations are "zero" and
        # "positive" and neither is defined by TURNS_TOL. A band bracketed at
        # TURNS_TOL would move with the constant it guards, which is the defect
        # this file records fixing for RANK_RTOL and fixes again for ACTIVE_TOL.
        tied = gaps[gaps <= 0.0]
        rest = gaps[gaps > 0.0]
        if tied.size:
            lo = max(lo, float(tied.max()))
        if rest.size:
            hi = min(hi, float(rest.min()))
        out[lab] = (k, int(np.count_nonzero(gaps < TURNS_TOL)))
    return out, lo, hi


# ==========================================================================
# Y5  AMPLITUDE AND WAVEVECTOR
# ==========================================================================

def y5_amplitude_wavevector(y0, y3):
    print()
    print("=" * 78)
    print("Y5  THE (AMPLITUDE, WAVEVECTOR) RELATION -- and why there is no")
    print("    threshold to put a curve through")
    print("=" * 78)
    print("  The bead expects a threshold on the neighbour phase difference,")
    print("  which for a plane wave scales like amplitude x wavevector and would")
    print("  therefore be a CURVE in (A, |k|) rather than a single k_min. That")
    print("  expectation presumes SECOND-ORDER relief. Y3 measures FIRST order")
    print("  everywhere, and a first-order excess is exactly homogeneous of")
    print("  degree one in A, so E(A, k) = A f(k) with no threshold in A at any")
    print("  k. The relation is a PRODUCT, not a curve, and the only thing left")
    print("  to report is f(k).")
    out = {}
    achord = y0["achord"]
    topos = {t.name: t for t in build_topologies()}

    print()
    print("  (a) homogeneity, over the whole amplitude ladder")
    homog = {}
    for key, (r1, r2) in y3["order"].items():
        homog[key] = (r1["slope"], r1["spread"], r2["slope"], r2["spread"])
    out["homog"] = homog
    ok = [abs(s - 1.0) < 1e-3 and sp < TOL["homog"]
          for s, sp, s2, sp2 in homog.values()]
    out["homog_ok"] = bool(ok) and all(ok)
    for key, (s, sp, s2, sp2) in homog.items():
        print(f"      {str(key):32s} slope {s:.9f} / {s2:.9f}   "
              f"E/eps spread {sp:.2e} / {sp2:.2e}")
    print("      Slope 1 and a spread at the solver's noise floor, on two")
    print("      incommensurate ladders. There is no amplitude at which the")
    print("      behaviour changes, so THERE IS NO THRESHOLD, hence no minimum")
    print("      wavevector and NO MAXIMUM WAVELENGTH. Outcome 2 of the bead is")
    print("      refuted by measurement, not set aside.")

    print()
    print("  (b) f(k) for the INTRA-unit binder: A DERIVATION, AND TWO GATE ROWS")
    print("      DELETED. A folding diagonal is a function of ONE unit's phase,")
    print("      so for a dephasing A cos(k s_i + phi) the intra rate at unit i")
    print("      is s'(a) A cos(k s_i + phi) and")
    print("          f(k) = sup_phi max_i s'(a) cos(k s_i + phi) = |s'(a)|,")
    print("      because the supremum over phi is attained by choosing phi to")
    print("      put SOME unit exactly at the trough, and a unit at the trough")
    print("      contributes s'(a) * (-1) = |s'(a)| whatever k is. There is no")
    print("      k left in the answer.")
    print()
    print("      TWO ROWS THAT USED TO SIT HERE HAVE BEEN REMOVED, and the")
    print("      removal is the point. They asserted that f(k) is flat and that")
    print("      it equals |s'(a_ico)|, and they were computed by a routine that")
    print("      took the sup over the candidate offsets {pi - k s_j} -- among")
    print("      which j = i gives cos(pi) exactly. The routine therefore")
    print("      returned |s'(a)| IDENTICALLY, for every k, every grid, every")
    print("      offset and every angle: the measured relative spread was")
    print("      exactly 0.0e+00 on both k-grids, which reads as strong evidence")
    print("      and is the tell that nothing was being measured. NO STATE OF")
    print("      THE ARRAY CAN REDDEN SUCH A ROW. It can be reddened only by")
    print("      mutating the code, which makes it a change-detector wearing a")
    print("      correctness label; this directory has retired one of those")
    print("      before (the rank-2 skew row of the sibling file) and the")
    print("      precedent is to state the derivation instead. Stated above.")
    print()
    print("      NOTE FOR THE RECORD: the two reviewers of this file SPLIT on")
    print("      these rows. The independent check suite reddened them by code")
    print("      mutation and counted them as covered; the critique held that no")
    print("      physical configuration can redden them and that they therefore")
    print("      say nothing about the geometry. Both readings are correct about")
    print("      what they checked. They are resolved the critique's way here,")
    print("      because a row that only a code edit can falsify measures the")
    print("      code and not the object, and the conclusion it carried is")
    print("      available for free as a derivation.")
    print(f"      THE CONCLUSION IS UNCHANGED: |s'(a_ico)| = "
          f"{abs(y3['diag_slope']):.9f} and the intra-unit")
    print("      binder carries NO wavevector information whatsoever.")

    print()
    print("  (c) f(k) for the INTER-unit binder. Here k does enter, because the")
    print("      chord's two ends sit on different units. The wavevectors a")
    print("      2x2x2 cluster carries are the eight zone points, and BOTH SIGNS")
    print("      of each are evaluated: D+ is positively homogeneous but NOT")
    print("      odd, so D+(w) and D+(-w) are different numbers and scanning one")
    print("      sign only can miss the optimum by its whole magnitude, which is")
    print("      exactly what an earlier version of this row did.")
    tt = topos["CUBE8-M (60/60/90 basis)"]
    mem, _ = _classes(tt)
    pr, _, _ = active_set(tt, mem, achord, "inter")
    zone = {}
    if pr:
        sites = np.array(tt.lattice_sites, float)
        for kx in range(2):
            for ky in range(2):
                for kz in range(2):
                    w = (-1.0) ** (kx * sites[:, 0] + ky * sites[:, 1]
                                   + kz * sites[:, 2])
                    w = w - w.mean()
                    mx = float(np.abs(w).max())
                    if mx < 1e-12:
                        continue
                    w = w / mx
                    for sgn in (+1.0, -1.0):
                        val, st = directional_derivative(tt, pr, achord, sgn * w)
                        zone[(kx, ky, kz, int(sgn))] = (val if st == 0
                                                        else float("inf"))
    out["zone"] = zone
    for kk, v in sorted(zone.items()):
        print(f"      k/pi = {kk[:3]} sign {kk[3]:+d}   D+ = {v:+.9e}")
    best = (min(zone.items(), key=lambda kv: kv[1]) if zone
            else (None, float("nan")))
    out["zone_best"] = best
    lp_best = y3["rows"].get(("CUBE8-M (60/60/90 basis)", "inter"),
                             {}).get("val", float("nan"))
    out["zone_matches_lp"] = (best[0] is not None and np.isfinite(lp_best)
                              and abs(best[1] - lp_best) < 1e-7)
    print(f"      best zone direction {best[0]}, D+ = {best[1]:+.9e};")
    print(f"      the LP optimum over the whole traceless sphere is "
          f"{lp_best:+.9e}")
    # WHICH generators the winning alternation runs along, compared against
    # jb_x's independently derived DIAGONAL GENERATOR PAIR. This is the row that
    # turns "some zone point wins" into a structural statement: the pair that
    # CREATES the taut chord is the pair the relieving wave alternates across,
    # and jb_x computes that pair from the a = 0 cuboctahedron with no reference
    # to anything here.
    dgp = diagonal_generator_pairs(tt.gens)
    out["dgp"] = dgp
    if best[0] is not None:
        won = tuple(i for i in range(3) if best[0][i] == 1)
        out["zone_axes"] = won
        out["zone_is_dgp"] = bool(dgp) and won == tuple(sorted(dgp[0]))
        print(f"      the winning alternation runs along generators {won};")
        print(f"      jb_x's diagonal generator pair for this basis is "
              f"{dgp[0] if dgp else None}")
    else:
        out["zone_axes"] = ()
        out["zone_is_dgp"] = False
    if out["zone_matches_lp"]:
        print("      -- the SAME number, so the optimal relieving dephasing IS a")
        print("      plane wave, and it is the alternation across the DIAGONAL")
        print("      GENERATOR PAIR, the very pair that creates the taut chord.")
        print("      Its wavelength is TWO lattice steps in each of those two")
        print("      directions: the SHORTEST the lattice carries, not the")
        print("      longest. That is a MINIMUM-wavelength statement, and it is")
        print("      the opposite of the maximum wavelength the bead expected.")
    else:
        print("      -- DIFFERENT numbers, so the optimal relieving dephasing is")
        print("      NOT a plane wave of this lattice and no wavelength can be")
        print("      attached to it. Reported as measured.")
    print()
    print("  (d) and the sting, WITH ITS SCOPE ATTACHED. Y1(e) measures the")
    print("      admissible bulk plane wave AT FIXED LATTICE PERIOD to be")
    print("      k = (pi,pi,pi) alone, and Y3 measures the 3x3x3 cluster -- the")
    print("      first with an interior -- to have a POSITIVE minimum. So the")
    print("      relieving pattern and the bulk-admissible plane wave are")
    print("      different waves, and the cluster sizes where they coincide are")
    print("      exactly the ones with no self-stress.")
    print("      THE SCOPE IS NOT DECORATION. Y1(e3) measures that the UNIFORM")
    print("      (k = 0) dephasing is executable on every finite cluster, by an")
    print("      exactly affine translation field with a non-scalar linear part")
    print("      -- a homogeneous lattice STRAIN, which is not a plane wave of")
    print("      any k and which the Bloch ansatz has no unknown for. So the")
    print("      statement 'k = (pi,pi,pi) is the only bulk-admissible wave' is")
    print("      true of STRICT PLANE WAVES AT FIXED PERIOD and is NOT a")
    print("      statement that nothing else is admissible in the bulk. The")
    print("      sting survives that scoping -- the relieving wave is k/pi =")
    print("      (1,0,1) and the uniform mode relieves nothing, being uniform --")
    print("      but it survives with the scope stated rather than without it.")
    return out



# ==========================================================================
# Y6  VERDICT
# ==========================================================================

def y6_verdict(y0, y1, y2, y3, y4, y5):
    print()
    print("=" * 78)
    print("Y6  VERDICT")
    print("=" * 78)
    print("  THE ORDER OF RELIEF IS FIRST, EVERYWHERE, AND ITS SIGN IS WHAT")
    print("  VARIES. Neither of the bead's two outcomes is what happened.")
    print()
    print("  1. THE INTRA-UNIT BINDER IS NOT RELIEVED, AT ANY ORDER.")
    print("     The six folding diagonals of a unit are functions of THAT unit's")
    print("     phase alone, so the constraint set is a product of half-lines and")
    print("     no dephasing can trade one unit's slack for another's. The")
    print("     minimum over the traceless sup-sphere of the one-sided")
    print(f"     directional derivative is STRICTLY POSITIVE, and equals")
    print(f"     |s'(a_ico)|/(N-1) = {abs(y3['diag_slope']):.9f}/(N-1) exactly in")
    print("     every cluster with no self-stress. The in-phase locus is a")
    print("     strict transverse minimum with a CORNER.")
    print()
    print("     WHICH ANSWER THE BEAD ASKED FOR: 'the derivative is ZERO' is")
    print("     FALSE. 'The derivative DOES NOT EXIST' is TRUE. D+(u) and D+(-u)")
    print("     are both strictly positive in every direction measured, so no")
    print("     linear functional can agree with both and there is no gradient.")
    print("     The Danskin condition is met in the measured form: the tied orbit")
    print("     has 6N members and N distinct gradients, the six inside one unit")
    print("     agreeing exactly, so the corner comes from the tie ACROSS units.")
    print()
    print("     BUT THE CORNER FLATTENS -- IN THE SUP-NORM BUDGET AND IN NO")
    print("     OTHER. The optimal cost falls as 1/(N-1), so in a large array")
    print("     the transverse slope out of the in-phase state is arbitrarily")
    print("     small while remaining strictly positive. Read the scope: the")
    print("     minimiser is one unit at +1 and the other N-1 at -1/(N-1), the")
    print("     least wave-like field available, and min over {sum u = 0,")
    print("     |u|_inf = 1} of max_i s' u_i is |s'|/(N-1) for ANY monotone s --")
    print("     no topology, no contact and no jitterbug enter, which is exactly")
    print("     why it is exact. Under an l2 budget the same argument gives")
    print("     |s'|/sqrt(N(N-1)); under a fixed phase difference PER LATTICE")
    print("     STEP, which is the normalisation a wave actually carries, there")
    print("     is no N-dependence at all. So this is the defensible form of the")
    print("     owner's 'easy for eight, impossible for 10^36' ONLY for a")
    print("     sup-norm budget on the dephasing: not that the lock dissolves,")
    print("     but that its transverse stiffness per unit of that budget")
    print("     vanishes with size. The genuinely informative row of the same")
    print("     table is CUBE27-M, where self-stress lifts min D+ strictly")
    print("     ABOVE the derived value -- that one is geometry.")
    print()
    print("  2. THE INTER-UNIT BINDER: RELIEVED AT FIRST ORDER IN CUBE8-M, NOT")
    print("     IN CUBE27-M, AND HELD BY THE LINKAGE IN ONLY TWO OF THE SIX")
    print("     CLUSTERS AT ALL.")
    # float("nan") rather than None: a missing key must print as a number the
    # reader can see is missing, not raise a TypeError inside the verdict and
    # take the whole gate table with it.
    v8 = y3["rows"].get(("CUBE8-M (60/60/90 basis)", "inter"),
                        {}).get("val", float("nan"))
    v27 = y3["rows"].get(("CUBE27-M", "inter"), {}).get("val", float("nan"))
    v8 = float("nan") if v8 is None else v8
    v27 = float("nan") if v27 is None else v27
    # The self-stress counts are READ FROM THE MEASUREMENT, not typed in. The
    # "15" here was a hardcoded literal in the prose, ungated, in a verdict
    # whose whole point is that the number is the discriminator.
    ss8 = y3["rows"].get(("CUBE8-M (60/60/90 basis)", "inter"),
                         {}).get("selfstress", -1)
    ss27 = y3["rows"].get(("CUBE27-M", "inter"), {}).get("selfstress", -1)
    print(f"     CUBE8-M  ({ss8:2d} self-stresses):  min D+ = {v8:+.9e}  "
          f"RELIEVED")
    print(f"     CUBE27-M ({ss27:2d} self-stresses):  min D+ = {v27:+.9e}  "
          f"NOT relieved")
    print("     In SQUARE4, the SC7 star and CUBE8-R the same chord is relieved")
    print("     by the assembly's own MECHANISMS with no dephasing at all, so")
    print("     there is no lock there to ask the question about -- see the")
    print("     correction below. Only the M-basis cubes hold it.")
    print("     Where it IS held, the relieving pattern exists as a phase field in")
    print("     both and is not KINEMATICALLY ADMISSIBLE in the larger one: the")
    print("     states of self-stress that appear once a cluster has an interior")
    print("     are exactly the functionals that forbid it. So the array-induced")
    print("     lock is measure zero in a cluster small enough to have no bulk and")
    print("     is NOT measure zero once it has one. THIS RUNS THE OPPOSITE WAY")
    print("     TO THE OWNER'S SCALING INTUITION and is reported as measured.")
    print("     It is also the narrower of the two locks: the intra-unit wall at")
    print(f"     {A_ICO_RECORD:.9f} survives in every cluster and at every "
          f"size")
    print("     MEASURED HERE, so no amount of dephasing takes the array below")
    print("     it. TWO SCOPES ON THAT SENTENCE, both of which the file used to")
    print("     leave off. (i) 'every topology' is wrong: the face-sharing IVM")
    print("     alternative is a different constraint set and is declared")
    print("     unmeasured. (ii) IT IS DOWEL-CONDITIONAL. The wall is a")
    print("     statement about a unit whose configuration IS a phase. In the")
    print("     free model a unit rides a six-dimensional internal variety on")
    print("     which the diagonal's length is a function of the whole internal")
    print("     state, so a non-symmetric internal deformation could relieve")
    print("     that diagonal WITH NO DEPHASING AT ALL, and the owner's rig")
    print("     rides the shared-vertex ELLIPSE rather than this file's true")
    print("     path tangent. That gap is jb_x's open P1, inherited here")
    print("     unresolved, and it sits above every other scope limit in this")
    print("     file including the tension-only-member question.")
    print()
    print("  3. NO THRESHOLD, HENCE NO MAXIMUM WAVELENGTH -- AND HOW FAR OUT")
    print("     THAT IS ACTUALLY MEASURED.")
    print("     The excess is homogeneous of degree one in the dephasing")
    print("     amplitude on two incommensurate ladders, over the whole range")
    print("     where it stands above the solver's floor. That window is |E| in")
    print("     [1e-9, 1e-6], i.e. dephasing of order 1e-7 to 1e-4 DEGREES,")
    print("     while the bead's Outcome 2 is a FINITE excursion -- so 'no")
    print("     threshold' does NOT follow from 'first order everywhere' and is")
    print("     not left to. For the INTRA binder it follows from monotonicity")
    print("     of s. For the INTER binder on CUBE27-M -- the only cluster where")
    print("     relief fails, hence the only place a threshold could live --")
    print("     Y3(g) carries it out to FOUR DEGREES on the exact nonlinear")
    print("     model along the zone-corner alternation, with every rung's")
    print("     closure residual gated, and the excess is strictly positive and")
    print("     strictly increasing at every rung. Two other directions tried")
    print("     did not close at finite amplitude and are reported as not")
    print("     measured rather than counted as agreeing.")
    print("     Nowhere measured does the behaviour change with amplitude, so")
    print("     the (amplitude, wavevector)")
    print("     relation is the product A f(k) and not a curve. f(k) is exactly")
    print("     FLAT for the intra-unit binder -- DERIVED in Y5(b), not")
    print("     measured, because the only routine that could measure it")
    print("     returned the answer identically and the two rows that read it")
    print("     have been deleted -- so that binder carries no wavevector")
    print("     information at all; and for the inter-unit binder")
    print(f"     its optimum is the zone point k/pi = {y5['zone_best'][0][:3] if y5['zone_best'][0] else None},")
    print(f"     an alternation across generators {y5['zone_axes']}, which is")
    print(f"     jb_x's DIAGONAL GENERATOR PAIR {y5['dgp'][0] if y5['dgp'] else None} -- the very pair")
    print("     that creates the taut chord. Wavelength TWO lattice steps along")
    print("     each of those two directions: a MINIMUM wavelength statement, the")
    print("     opposite of what the bead anticipated. And D+ at that point is")
    print("     SIGN-SENSITIVE: one sign of the alternation relieves by exactly")
    print("     as much as the other worsens, which is the antisymmetry a")
    print("     dephasing wave has and a uniform angle change cannot.")
    print()
    print("  4. STATIC AND PATH AGREE AT FIRST ORDER. The bead's proposed")
    print("     resolution -- that the static question and the taking-turns")
    print("     question might differ -- is EMPTY AT FIRST ORDER, measured on")
    print("     every row of Y4(a). The qualifier is load-bearing and is not")
    print("     droppable: both are linear programmes over the SAME tangent cone")
    print("     at the SAME point, so neither could have distinguished a static")
    print("     relaxation from a path in any stronger sense whatever the")
    print("     numbers had been. Both ask whether the traceless cone {all")
    print("     member rates <= 0} is more than the origin, and letting the mean")
    print("     fall buys nothing because a falling mean lengthens every active")
    print("     span by construction. Of the six rows, three intra ones are")
    print("     derivable to zero and the CUBE27-M inter row is 0/0, so two rows")
    print("     carry the information.")
    print()
    print("  5. WHAT DOES RESOLVE THE TENSION, and it is not the static/path")
    print("     split. Two different objects were being compared. AS A")
    print("     FUNCTIONAL the excess has a strict positive-slope corner at the")
    print("     in-phase locus: dephasing makes the worst span worse. AS A SET")
    print("     the fully arrested state is a single point -- the feasible")
    print("     dephasing polytope at mean-phase slack eta has sup-radius exactly")
    print("     eta and collapses to a point at eta = 0. BOTH ARE TRUE, and THAT")
    print("     is the resolution: they are claims about different objects and")
    print("     neither refutes the other, so the appearance that they did was a")
    print("     category error.")
    print("     WHERE THE MEASUREMENT ENDS AND THE READING BEGINS. Measured: the")
    print("     radius equals the slack, on CUBE8-M, for the intra binder, on")
    print("     the one cluster here with no self-stress -- and it is close to a")
    print("     tautology, since feasibility there is a_i >= a_lock for every i")
    print("     and the positive orthant meets the mean-zero hyperplane only at")
    print("     the origin. Not measured: the same radius for the inter binder,")
    print("     or on a cluster where the Y1 obstruction bites. INTERPRETATION,")
    print("     offered as such: that the owner's observation is the SET")
    print("     statement. He said motion PROGRESSES given a fractional")
    print("     difference, and 'does motion progress' is the functional/path")
    print("     side, where this file's intra answer is that it does not. The")
    print("     assignment is a reading. It is not what makes both claims true.")
    print()
    print("  5b. THE STRONGEST SCALING STATEMENT HERE, AND IT IS NEW.")
    print("      ADMISSIBLE DEPHASING SATURATES. Measured on three cubes (Y1c1):")
    sat = y1.get("sat", {})
    for box, (nn, _, dd, mr, _x) in sat.items():
        print(f"        {str(box):10s} N = {nn:4d}   dephasing dimension "
              f"{dd:3d}   min retention {mr:.9f}")
    print("      The dimension does not grow with N -- it is the same number at")
    print("      27, 64 and 125 units -- and all eight of AFFINE x PARITY are")
    print("      retained exactly at every size. So admissible dephasing is a")
    print("      FIXED FINITE-DIMENSIONAL AFFINE-TIMES-PARITY SPACE INDEPENDENT")
    print("      OF ARRAY SIZE, and the admissible FRACTION goes to zero as the")
    print("      array grows. That is directly on the bead's question: a medium")
    print("      large enough to be a medium has essentially no dephasing its")
    print("      linkage can execute, and what it does have is a uniform")
    print("      gradient and a zone-corner alternation and nothing between.")
    print("      SCOPE: three cubes, M basis, icosahedral phase.")
    print("      This also corrects a sentence Y1(b) used to print. The")
    print("      obstruction is NOT a cycle property: four boxes here carry")
    print("      independent cycles and no self-stress, and the star carries a")
    print("      six-fold-coordinated unit and no self-stress. It takes BOTH --")
    print("      a cycle and a unit with at least four contacts -- and that")
    print("      conjunction agrees with the measured self-stress on every")
    print(f"      cluster in this file ({y1.get('conjunction_hits', 0)} of "
          f"{y1.get('conjunction_n', 0)}). Reported as a measured")
    print("      correlation over that table; no derivation is claimed.")
    print()
    print("  6. WHAT PROPAGATION THIS SUPPORTS, STATED NARROWLY. An arrest that")
    print("     arrives unit by unit rather than all at once is a FRONT, and Y4c")
    print("     measures that a commensurate dephasing arrests a whole plane at")
    print("     one instant while an incommensurate one arrests exactly one unit")
    print("     at a time. That is a kinematic statement about the ORDER of")
    print("     arrival. It is not a wave: there is no equation of motion in this")
    print("     file, no speed, and no dispersion relation, and the four")
    print("     declarations are inapplicable precisely because none of that is")
    print("     here. 'The motion must propagate' is NOT established. What is")
    print("     established is that uniform arrest is a single point and")
    print("     non-uniform arrest is everything else.")
    print()
    print("  CORRECTIONS THIS FILE MAKES TO jb_x")
    print("  -----------------------------------")
    print("  jb_x's inter-unit taut angles are PLACEMENT-CONDITIONAL. Every span")
    print("  there is evaluated at the pure-translate reference placement, but at")
    print("  fixed phases the assembly still has internal mechanisms, and a")
    print("  mechanism moves inter-unit spans. Measured here (Y3c): for the")
    print("  SIX-AROUND-ONE STAR the mechanisms alone relieve the inter-unit")
    print(f"  chord at {y3['star_inter_mech']:+.6e} per unit placement rate, and")
    print("  SQUARE4 and CUBE8-R do the same, so in three of the five clusters")
    print("  that carry the chord its taut angle is an artefact of holding the")
    print("  reference placement rather than a property of the linkage. Only the")
    print(f"  M-basis cubes hold it against their own mechanisms. THE "
          f"{A_CHORD_RECORD:.6f}")
    print("  ANGLE THEREFORE NEEDS A TOPOLOGY ATTACHED TO IT, and T2 R1's table")
    print(f"  of per-cluster counts -- reproduced exactly here, "
          f"{'/'.join(str(g) for g, _ in y0['counts'].values())} -- counts")
    print("  spans that are taut at a placement, not members that bind.")
    print(f"  The intra-unit angle {A_ICO_RECORD:.9f} is untouched by this "
          f"correction --")
    print("  a diagonal joins two vertices of ONE rigid unit and no placement can")
    print("  reach it -- which is measured, not asserted, in the same row: the")
    print("  placement block of every intra row is identically zero.")


# ==========================================================================
# THE GATE
# ==========================================================================

def gate(y0, y1, y2, y3, y4, y5):
    """Every check's verdict in one table, and this process's exit code.

    ON WHAT IS IN THE SECTION DICTS AND NOT READ HERE. Six keys carry RAW
    MEASUREMENT TABLES -- `per`, `graph`, `obstruct`, `affine`, `derived`,
    `excursion` -- which are printed in full in their own sections and are
    reduced to gated scalars before they reach this function. They are the
    section's record, not decisions dressed as checks, and that is the
    distinction the recurring defect in this directory turns on: a value
    computed AS IF it were a guard and then read by nothing is the defect; a
    table printed for a reader and summarised into a gated number is not. Every
    scalar that looked like a guard and was not -- `bad_free`, `bad_stressed`,
    `gaps`, `res_max`, `order_one_sign`, `order_nfit_min` -- is read below.
    """
    rows = y3["rows"]
    r8 = rows.get(("CUBE8-M (60/60/90 basis)", "inter"), {})
    r27 = rows.get(("CUBE27-M", "inter"), {})
    i8 = rows.get(("CUBE8-M (60/60/90 basis)", "intra"), {})
    corner_i = y3["corner"].get("intra", {})
    corner_x = y3["corner"].get("inter", {})
    homog_vals = list(y5["homog"].values())
    slopes = [s for s, sp, s2, sp2 in homog_vals] + \
             [s2 for s, sp, s2, sp2 in homog_vals]
    spreads = [sp for s, sp, s2, sp2 in homog_vals] + \
              [sp2 for s, sp, s2, sp2 in homog_vals]
    nbad = sum(v.get("bad", 0) for v in rows.values())
    scal = y1["scal"]
    n_free = len(y1["free_boxes"])
    n_obs = len(y1["obstructed_boxes"])
    zone = y5["zone"]
    zone_min = min(zone.values()) if zone else float("nan")
    zone_max = max(zone.values()) if zone else float("nan")
    # FOUR AGGREGATES THAT USED TO BE WRITTEN BARE. `min` and `max` over an
    # empty iterable RAISE, and an index into a dict key that is not there
    # raises too -- and a raise inside `gate` destroys the verdict table, which
    # is the one failure mode this file's whole discipline exists to prevent.
    # None is reachable from the current data; all four are reachable if an
    # upstream section returns nothing, which is exactly when a reader most
    # needs the table. Each now defaults to a value that FAILS RED.
    decisive = [abs(x) for x in (i8.get("val", 0.0), r8.get("val", 0.0),
                                 r27.get("val", 0.0)) if np.isfinite(x)]
    decisive_min = min(decisive) if decisive else 0.0
    angdep_vals = list(y1["angdep"].values())
    angdep_worst = max((len(set(v)) for v in angdep_vals), default=0)
    angdep_m333 = y1["angdep"].get(("M", (3, 3, 3)), (-1,))[0]
    angdep_r333 = y1["angdep"].get(("R", (3, 3, 3)), (-1,))[0]
    ret_vals = list(y1["retention"].values())
    ret_min = min(ret_vals) if ret_vals else 0.0
    ctrl_vals = list(y1["retention_ctrl"].values())
    ctrl_max = max(ctrl_vals) if ctrl_vals else 1.0
    mech_intra = [g for (_, tg_), (_, g, _) in y3["mech"].items()
                  if tg_ == "intra"]
    mech_intra_max = max(mech_intra) if mech_intra else 1.0
    # UNBOUNDED must coincide with MECHANISM RELIEF, on every row, both ways.
    # That is the cross-check that the status-3 reading is the right one and
    # not a way of hiding a solver failure as a result.
    unb_vs_mech = []
    for key, v in rows.items():
        if v.get("n", 0) == 0:
            continue
        mv = y3["mech"].get(key)
        if mv is None:
            continue
        unb_vs_mech.append(bool(v.get("unb", 0)) == bool(mv[0] < -1e-9))
    n_unb = sum(1 for v in rows.values() if v.get("unb", 0))
    # Infeasible pins are legitimate ONLY where the cluster carries self-stress:
    # with no self-stress every dephasing is executable and every pin must
    # solve. BOTH HALVES ARE NOW ACTUALLY ASSERTED. They were computed here
    # under this same comment and read by no row at all -- prose claiming a
    # guard that did not exist, inside the gate, which is the worst place for
    # it. The rows are `Y3 no infeasible pin where the cluster is FREE` and the
    # non-vacuity companion below it.
    bad_free = sum(v.get("bad", 0) for v in rows.values()
                   if v.get("selfstress", 1) == 0)
    bad_stressed = sum(v.get("bad", 0) for v in rows.values()
                       if v.get("selfstress", 0) > 0)
    n_free_rows = sum(1 for v in rows.values()
                      if v.get("n", 0) and v.get("selfstress", 1) == 0)
    n_stressed_rows = sum(1 for v in rows.values()
                          if v.get("n", 0) and v.get("selfstress", 0) > 0)

    checks = [
        ("Y0  zero dephasing reproduces jb_x's reference", y0["inphase_ok"],
         f"{y0['inphase_dev']:.2e}", f"< {TOL['inphase']:.0e}"),
        # TOL[inphase] and TOL[fd_jacobian] were the two FOUNDATION tolerances
        # bounded from below by their own deviations and from above by nothing:
        # an independent validation loosened each to 1.0 with the gate green.
        # Each now carries the achord/aico two-row idiom -- a control that must
        # be REJECTED, and a band the tolerance must sit inside.
        ("Y0  CONTROL: rejects a reference offset by 1e-9",
         y0["inphase_ctrl"] > TOL["inphase"], f"{y0['inphase_ctrl']:.2e}",
         f"> {TOL['inphase']:.0e}"),
        ("Y0  ... and TOL[inphase] is inside its measured band",
         y0["inphase_dev"] <= TOL["inphase"] < INPHASE_CONTROL_OFFSET,
         f"{TOL['inphase']:.0e}",
         f"{y0['inphase_dev']:.1e}..{INPHASE_CONTROL_OFFSET:.0e}"),
        ("Y0  analytic (Jr,Jp) vs exact-residual FD", y0["fd_ok"],
         f"{y0['fd']:.2e}", f"< {TOL['fd_jacobian']:.0e}"),
        ("Y0  CONTROL: rejects that Jacobian scaled by 1+1e-3",
         y0["fd_ctrl"] > TOL["fd_jacobian"], f"{y0['fd_ctrl']:.2e}",
         f"> {TOL['fd_jacobian']:.0e}"),
        ("Y0  ... and TOL[fd_jacobian] is inside that band",
         y0["fd"] < TOL["fd_jacobian"] < y0["fd_ctrl"],
         f"{TOL['fd_jacobian']:.0e}", f"{y0['fd']:.1e}..{y0['fd_ctrl']:.1e}"),
        ("Y0  a_ico re-derived vs the record", y0["aico_dev"] < TOL["aico"],
         f"{y0['aico_dev']:.2e}", f"< {TOL['aico']:.0e}"),
        ("Y0  CONTROL: rejects a_ico offset by 1e-3",
         y0["aico_ctrl"] > TOL["aico"], f"{y0['aico_ctrl']:.2e}",
         f"> {TOL['aico']:.0e}"),
        ("Y0  ... and TOL[aico] is inside its derived band",
         AICO_RECORD_QUANTUM < TOL["aico"] < AICO_CONTROL_OFFSET,
         f"{TOL['aico']:.0e}", "5e-10..1e-03"),
        ("Y0  chord angle re-derived vs T2 R1's record",
         y0["achord_dev"] < TOL["achord"], f"{y0['achord_dev']:.2e}",
         f"< {TOL['achord']:.0e}"),
        ("Y0  CONTROL: rejects that angle offset by 1e-3",
         y0["achord_ctrl"] > TOL["achord"], f"{y0['achord_ctrl']:.2e}",
         f"> {TOL['achord']:.0e}"),
        ("Y0  ... and TOL[achord] is inside its derived band",
         ACHORD_RECORD_QUANTUM < TOL["achord"] < ACHORD_CONTROL_OFFSET,
         f"{TOL['achord']:.0e}", "5e-07..1e-03"),
        # THE BAND IS NOW COMPUTED WITHOUT REFERENCE TO ACTIVE_TOL. The old
        # one was bracketed by "largest excess kept" and "smallest discarded",
        # both computed WITH the tolerance, so it moved with what it guarded
        # and passed over three hundred decades. Two companions: the lower end
        # must be a real member of the population rather than the clip, and
        # the gap it brackets must be decades wide rather than a rounding step.
        ("Y0  ACTIVE_TOL inside its MEASURED two-sided band", y0["tol_ok"],
         f"{ACTIVE_TOL:.0e}", f"{y0['tol_lo']:.1e}..{y0['tol_hi']:.1e}"),
        ("Y0  ... its lower end is data, not the precision floor",
         y0["tol_lo_not_floor"], f"{y0['tol_lo']:.2e}",
         f"> {EXCESS_FLOOR:.1e}"),
        ("Y0  ... and the population gap is decades wide",
         y0["tol_gap_wide"], f"{y0['tol_gap_decades']:.1f}",
         f"> {EXCESS_GAP_DECADES:.0f}"),
        ("Y0  intra actives are exactly 6 per unit", y0["per_unit_ok"],
         f"{y0['n_intra']}", "6N"),
        # LSTSQ_RCONDS EXERCISED, BOTH WAYS. Without these two rows the solver
        # fix is a change no check distinguishes from its absence -- which is
        # how a fix becomes folklore.
        ("Y0  truncated GN step CLOSES a self-stressed box",
         0.0 < y0["rcond_ladder_res"] < SOLVE_TOL,
         f"{y0['rcond_ladder_res']:.2e}", f"0 < .. < {SOLVE_TOL:.0e}"),
        ("Y0  CONTROL: the UNTRUNCATED step does not",
         y0["rcond_plain_res"] > SOLVE_TOL, f"{y0['rcond_plain_res']:.2e}",
         f"> {SOLVE_TOL:.0e}"),
        ("Y0  no root-finder failure was swallowed",
         y0["chord_nfail"] == 0, f"{y0['chord_nfail']}", "0"),
        ("Y0  chord length at a_ico vs T2 R1's 1.0705x strut",
         y0["chord_ratio_ok"], f"{y0['chord_ratio']:.6f}",
         f"{CHORD_RATIO_RECORD} +/- 5e-5"),
        ("Y0  per-cluster chord counts vs T2 R1", y0["counts_ok"],
         "/".join(str(g) for g, _ in y0["counts"].values()),
         "/".join(str(w) for _, w in y0["counts"].values())),
        # NON-VACUITY. `all()` over an empty dict is True: a probe that emptied
        # CHORD_COUNT_RECORD reddened NOTHING and the whole T2 cross-check was
        # deletable with the gate green.
        ("Y0  ... over >= 4 clusters with positive counts",
         y0["counts_nonvacuous"], f"{y0['counts_n']}", ">= 4"),

        # The value column must carry a COMPUTED number. A row that prints the
        # literal "True" next to its own predicate tells a reader nothing when
        # it goes red, and three rows here did exactly that until a mutation
        # probe produced the line "FAIL ... True  True".
        ("Y1  self-stress is the whole discriminator",
         y1["selfstress_is_the_discriminator"],
         f"{sum(1 for (n, ss, d) in scal.values() if (ss == 0) == (d == n))}"
         f"/{len(scal)}", f"{len(scal)}/{len(scal)}"),
        ("Y1  ... and BOTH classes are non-empty", n_free > 0 and n_obs > 0,
         f"{n_free}/{n_obs}", "> 0 both"),
        # THE HALF OF THAT ROW THAT IS NOT A CODE BRANCH. `selfstress == 0 =>
        # full dimension` is an early return in `admissible_dephasing`, so
        # those rows are an identity. The obstruction rank being STRICTLY
        # POSITIVE on every stressed box is the measurement.
        ("Y1  ... obstruction rank > 0 on every STRESSED box",
         y1["obstruct_positive"], f"{y1['n_stressed']}", "> 0 rows, all"),
        # NOT A CYCLE PROPERTY, and not a coordination property either. Each
        # witness list must be NON-EMPTY: they are what makes the conjunction
        # a claim rather than a restatement.
        ("Y1  cycles alone do NOT imply self-stress (witness)",
         len(y1["cycles_without_selfstress"]) > 0,
         f"{len(y1['cycles_without_selfstress'])}", "> 0"),
        ("Y1  coordination alone does NOT either (witness)",
         len(y1["coord_without_selfstress"]) > 0,
         f"{len(y1['coord_without_selfstress'])}", "> 0"),
        ("Y1  (cycle AND coord>=4) == self-stress, all clusters",
         y1["conjunction_ok"],
         f"{y1['conjunction_hits']}/{y1['conjunction_n']}", "all"),
        # THE SATURATION. Three sizes, because two agreeing numbers are a
        # coincidence; and N must actually differ across them, or "constant in
        # N" is a statement about one cluster measured three times.
        ("Y1  dephasing dimension SATURATES over 3 cube sizes",
         y1["sat_constant"], f"{y1['sat_dim']}", "one value"),
        ("Y1  ... and N really grows across them (non-vacuity)",
         y1["sat_n_grows"] and len(y1["sat"]) >= 3,
         f"{sorted(v[0] for v in y1['sat'].values())}", "3 distinct"),
        ("Y1  ... affine x parity retained 1.0 at every size",
         y1["sat_ret_ok"], f"{y1['sat_min_ret']:.9f}", "1.000000000"),
        # ... against a CONTROL at every size. Without it the retention row is
        # satisfied by a projector that retains everything, and a probe
        # replacing it with the identity reddened nothing.
        ("Y1  ... CONTROL: single-direction parity is NOT, at any",
         y1["sat_ctrl_ok"], f"{y1['sat_ctrl_max']:.4f}", "< 0.9"),
        ("Y1  a 9-unit CHAIN dephases as freely as a pair",
         scal.get((9,), (0, -1, -1))[2] == scal.get((9,), (0, 0, 0))[0],
         f"{scal.get((9,), (0, 0, 0))[2]}/9", "9/9"),
        ("Y1  2x2x2 unobstructed, 3x3x3 obstructed",
         scal.get((2, 2, 2), (0, 1, 0))[1] == 0
         and scal.get((3, 3, 3), (0, 0, 0))[1] > 0,
         f"{scal.get((2,2,2),(0,-1,0))[1]}/{scal.get((3,3,3),(0,-1,0))[1]}",
         "0 / > 0"),
        ("Y1  dephasing dimension is angle-independent",
         y1["angle_independent"], f"{angdep_worst}", "1"),
        ("Y1  ... over >= 3 DISTINCT angles spanning > 10 deg",
         y1["adm_angles_ok"],
         f"{y1['n_adm_angles']}/{y1['adm_angle_spread']:.1f}", ">=3 / >10"),
        ("Y1  ... and IS basis-sensitive (non-vacuity)",
         y1["basis_sensitive"], f"{angdep_m333}v{angdep_r333}", "differ"),
        ("Y1  affine x parity fully admissible on 3x3x3",
         y1["affine_parity_ok"], f"{ret_min:.9f}", "1.000000000"),
        ("Y1  CONTROL: single-direction parity is NOT",
         y1["ctrl_ok"], f"{ctrl_max:.4f}", "< 0.9"),
        ("Y1  Bloch: no bulk plane wave but the zone corner",
         y1["bloch_min"] > 0.1 and y1["bloch_min_alt"] > 0.1,
         f"{min(y1['bloch_min'], y1['bloch_min_alt']):.3e}", "> 0.1"),
        ("Y1  ... and both grids are off-lattice (0 degenerate)",
         y1["bloch_degen"] == 0, f"{y1['bloch_degen']}", "0"),
        # THE TWO GRIDS MUST BE TWO GRIDS. Probes that made the second grid
        # identical to the first, and that cut both to 5^3, each reddened
        # NOTHING. The amplitude ladders carry this companion; the k-grids did
        # not.
        ("Y1  ... the two k-grids are coprime and offset apart",
         y1["kgrid_distinct"], f"{K_GRID}/{K_GRID_ALT}", "coprime"),
        ("Y1  ... and each resolves at least 15^3 points",
         y1["kgrid_resolved"], f"{y1['kgrid_min']}", ">= 15"),
        ("Y1  ... the zone corner IS free (non-vacuity)",
         y1["bloch_corner_nb"] < BLOCH_B_FLOOR,
         f"{y1['bloch_corner_nb']:.1e}", f"< {BLOCH_B_FLOOR:.0e}"),
        ("Y1  ... and a nearby control is NOT",
         y1["bloch_ctrl_nb"] > BLOCH_B_FLOOR and y1["bloch_corner_ctrl"] > 0.1,
         f"{y1['bloch_corner_ctrl']:.2e}", "> 0.1"),

        # THE SCOPE OF THE BLOCH SCAN, ASSERTED RATHER THAN ONLY PRINTED. The
        # ansatz has no unknown for a homogeneous lattice strain, so k = 0 must
        # read INADMISSIBLE in it -- while the finite clusters execute exactly
        # that motion, by an exactly affine translation field with a NON-SCALAR
        # linear part. Both sides are gated, because the previous
        # reconciliation ("the boundary relaxes the self-stress") was prose and
        # was wrong: retention does not decay with size.
        ("Y1  k=0 IS inadmissible at fixed lattice period",
         y1["bloch_k0"] > 0.1, f"{y1['bloch_k0']:.3e}", "> 0.1"),
        # ... and NOT for the trivial reason. The zone corner is free because
        # its phase column VANISHES; a k = 0 row reporting a large relative
        # residual on a vanishing column would be meaningless in the same way.
        # This is the companion the corner rows already carry, applied to k = 0.
        ("Y1  ... on a phase column that does NOT vanish",
         y1["bloch_k0_nb"] > BLOCH_B_FLOOR, f"{y1['bloch_k0_nb']:.1e}",
         f"> {BLOCH_B_FLOOR:.0e}"),
        ("Y1  ... yet finite clusters DO execute uniform u",
         y1["affine_exec_ok"],
         f"{y1['affine_worst_exec']:.1e}", "< 1e-09"),
        ("Y1  ... by an EXACTLY AFFINE translation field",
         y1["affine_fit_ok"],
         f"{y1['affine_worst_fit']:.1e}", "< 1e-09"),
        ("Y1  ... whose linear part is NOT scalar (a real strain)",
         y1["affine_nonscalar_ok"], f"{y1['affine_nonscalar']:.2e}",
         "> 1e-3"),
        ("Y1  RANK_RTOL inside its MEASURED two-sided band", y1["rtol_ok"],
         f"{RANK_RTOL:.0e}", f"{y1['rtol_lo']:.1e}..{y1['rtol_hi']:.1e}"),
        ("Y2  TOL[tie] inside its MEASURED two-sided band",
         y2["tie_band_ok"], f"{TOL['tie']:.0e}",
         f"{y2['intra']['within']:.1e}..{y2['tie_across']:.1e}"),
        ("Y2  intra orbit: N distinct gradients", y2["intra_ngrad_is_n"],
         f"{y2['intra']['ngrad']}", f"{y2['n_units']}"),
        ("Y2  ... and the six inside one unit are TIED",
         y2["intra_within_tied"], f"{y2['intra']['within']:.1e}",
         f"< {TOL['tie']:.0e}"),

        ("Y3  no infeasible pin anywhere in the sweep",
         nbad == 0, f"{nbad}", "0"),
        # THE TWO HALVES THE COMMENT ABOVE `bad_free` ALREADY CLAIMED WERE
        # ASSERTED, AND WERE NOT. With no self-stress every dephasing is
        # executable, so every pin must solve; where there IS self-stress an
        # infeasible pin is a legitimate result. Both classes must be occupied
        # or the first half is a statement about no rows at all.
        ("Y3  ... none of them on a cluster with NO self-stress",
         bad_free == 0, f"{bad_free}", "0"),
        ("Y3  ... and both classes of row occur (non-vacuity)",
         n_free_rows > 0 and n_stressed_rows > 0,
         f"{n_free_rows}/{n_stressed_rows}", "> 0 both"),
        # AND THE SPLIT MUST BE EXHAUSTIVE. Without this, "no infeasible pin on
        # a free cluster" is green whenever the classification silently drops
        # rows -- `bad_free` would count nothing and pass. `bad_stressed` is
        # the other half of the partition and it is READ here; it was computed
        # and discarded in the version this fix pass started from, under a
        # comment claiming both halves were asserted.
        # THIS ROW IS A DATA-INTEGRITY GUARD, NOT A MATHEMATICAL IDENTITY. It
        # does not hold unconditionally by construction: `bad_free` and
        # `bad_stressed` classify each row by `v.get("selfstress", 1) == 0`
        # and `v.get("selfstress", 0) > 0` respectively, so a row dict that is
        # missing the `selfstress` key falls through BOTH filters -- its `bad`
        # count (if any) lands in neither bucket, while `nbad` above counts it
        # regardless of that key. Such a row would make this equality FAIL,
        # which is the intended catch: completeness of the free/stressed split
        # over every pin record, not an algebraic tautology.
        ("Y3  ... and that split accounts for every pin",
         bad_free + bad_stressed == nbad,
         f"{bad_free}+{bad_stressed}", f"= {nbad}"),
        ("Y3  CONTROL: an INADMISSIBLE direction IS infeasible",
         y3["adm_ctrl_ok"],
         str(y3["adm_ctrl"].get("single-direction parity (INADMISSIBLE)",
                                (0, -1))[1]), "!= 0"),
        ("Y3  UNBOUNDED coincides with mechanism relief, both ways",
         bool(unb_vs_mech) and all(unb_vs_mech),
         f"{sum(unb_vs_mech)}/{len(unb_vs_mech)}", "all"),
        ("Y3  ... and UNBOUNDED is reported as -inf, not +inf",
         all(v.get("val") == -float("inf") for v in rows.values()
             if v.get("unb", 0)),
         f"{sorted({v.get('val') for v in rows.values() if v.get('unb', 0)})}",
         "[-inf]"),
        ("Y3  ... and BOTH cases occur (non-vacuity)",
         0 < n_unb < len([v for v in rows.values() if v.get('n', 0)]),
         f"{n_unb}", "0 < n < rows"),
        ("Y3  INTRA min D+ is STRICTLY POSITIVE",
         np.isfinite(i8.get("val", np.nan)) and i8.get("val", -1) > TOL["lp"],
         f"{i8.get('val', float('nan')):+.6e}", f"> {TOL['lp']:.0e}"),
        ("Y3  ... = |s'|/(N-1), derived then measured", y3["derived_ok"],
         f"{abs(y3['diag_slope']) / 7:.9e}",
         f"{i8.get('val', float('nan')):.9e}"),
        # THE INEQUALITY ARM. It used to read `v > p - 1e-12`, which ACCEPTS
        # v == p -- exactly what the enumerator the prose says it catches would
        # return, and a probe replacing the arm with True reddened nothing. The
        # margin is now strict and bounded from ABOVE by the separation
        # actually measured, so it cannot be raised until it excludes the truth.
        ("Y3  ... DERIVED_MARGIN inside its measured band",
         y3["derived_margin_ok"], f"{DERIVED_MARGIN:.0e}",
         f"{TOL['lp']:.0e}..{y3['derived_sep']:.1e}"),
        ("Y3  ... and both classes of cluster occur (non-vacuity)",
         y3["derived_n_free"] > 0 and y3["derived_n_stressed"] > 0,
         f"{y3['derived_n_free']}/{y3['derived_n_stressed']}", "> 0 both"),
        # "the same on all six diagonals" was printed by the verdict and gated
        # by nothing: `_diag_slope` averaged the six under a docstring
        # promising an assertion. The deviation is returned and read now, and
        # the slope's own magnitude is gated so six zeros cannot pass by
        # agreeing perfectly.
        ("Y3  the six folding diagonals share one slope",
         y3["diag_slope_agree"], f"{y3['diag_slope_dev']:.2e}",
         "< 1e-08 (& > 0)"),
        ("Y3  INTER on CUBE8-M is RELIEVED at first order",
         np.isfinite(r8.get("val", np.nan)) and r8.get("val", 1) < -TOL["lp"],
         f"{r8.get('val', float('nan')):+.6e}", f"< -{TOL['lp']:.0e}"),
        ("Y3  INTER on CUBE27-M is NOT (self-stress bites)",
         np.isfinite(r27.get("val", np.nan)) and r27.get("val", -1) > TOL["lp"],
         f"{r27.get('val', float('nan')):+.6e}", f"> {TOL['lp']:.0e}"),
        # THE FINITE EXCURSION. "No threshold" does not follow from "first
        # order everywhere": the ladder window is three to six decades below a
        # finite phase difference, and CUBE27-M is the only cluster where a
        # threshold could live. Closure is gated FIRST, because an excess read
        # off a configuration that did not close measures nothing.
        ("Y3  finite excursion on CUBE27-M CLOSES at every rung",
         0.0 < y3["exc_res_max"] < SOLVE_TOL, f"{y3['exc_res_max']:.2e}",
         f"0 < .. < {SOLVE_TOL:.0e}"),
        ("Y3  ... CONTROL: a non-closing direction is CAUGHT",
         y3["exc_ctrl_res"] > SOLVE_TOL, f"{y3['exc_ctrl_res']:.2e}",
         f"> {SOLVE_TOL:.0e}"),
        ("Y3  ... excess stays POSITIVE out to 4 degrees",
         y3["exc_positive"], f"{y3['exc_top']:g} deg", f"> {TOL['lp']:.0e}"),
        # THE OBJECT-MUTATING FALSIFIER FOR THE ROW ABOVE, GATED. Without
        # this, the row's only coverage moves TOL["lp"] -- a change-detector
        # shape. This CONTROL swaps which active-set member each rung reads
        # (min instead of max, same solve) and asserts the result is
        # NEGATIVE, so the row is proven sensitive to the geometry rather
        # than to its own threshold. Mirrors `exc_ctrl_res`'s idiom above.
        ("Y3  ... CONTROL: the min-swapped falsifier IS negative",
         y3["exc_falsifier_ok"], f"{y3['exc_falsifier_max']:+.3e}",
         f"< -{TOL['lp']:.0e}"),
        ("Y3  ... and strictly INCREASES: no turning point",
         y3["exc_monotone"], f"{y3['exc_n']}", "rungs, all up"),
        ("Y3  ... over at least 6 rungs reaching >= 1 degree",
         y3["exc_n"] >= 6 and y3["exc_top"] >= 1.0,
         f"{y3['exc_n']}/{y3['exc_top']:g}", ">= 6 / >= 1"),
        ("Y3  intra spans do not see placement (Gz == 0)",
         y3["intra_gz_zero"], f"{mech_intra_max:.1e}", "< 1e-14"),
        ("Y3  star's inter chord IS mechanism-relieved",
         np.isfinite(y3["star_inter_mech"]) and y3["star_inter_mech"] < -1e-3,
         f"{y3['star_inter_mech']:+.3e}", "< -1e-3"),
        ("Y3  CORNER: D+(u) and D+(-u) both > 0, all dirs",
         corner_i.get("n", 0) > 0
         and corner_i.get("both_pos", -1) == corner_i.get("n", 0),
         f"{corner_i.get('both_pos', -1)}/{corner_i.get('n', 0)}", "all"),
        # The spread's SIZE is part of its meaning: "positive in every direction
        # measured" over one direction is a statement about one direction. A
        # probe that reduced N_SPREAD to 1 left the corner rows green.
        ("Y3  ... over a spread of at least 8 directions",
         corner_i.get("n", 0) >= 8 and corner_x.get("n", 0) >= 8,
         f"{min(corner_i.get('n', 0), corner_x.get('n', 0))}", ">= 8"),
        ("Y3  ... so no odd functional fits (non-vacuity)",
         corner_i.get("min_sum", -1.0) > TOL["lp"],
         f"{corner_i.get('min_sum', float('nan')):.2e}", f"> {TOL['lp']:.0e}"),
        ("Y3  ... and the INTER corner is two-sided too",
         corner_x.get("n", 0) > 0 and corner_x.get("min_sum", -1.0) > TOL["lp"],
         f"{corner_x.get('min_sum', float('nan')):.2e}", f"> {TOL['lp']:.0e}"),

        # NOT "they agree" -- whether they agree is the MEASUREMENT, and a gate
        # row asserting the answer would have to be rewritten whenever the
        # answer changed, which is the definition of a row that proves nothing.
        # What IS gated is the DERIVED value on the rows where a derivation
        # exists, plus non-vacuity.
        ("Y4  intra rows: STATIC and PATH both exactly 0 (derived)",
         y4["intra_zero_ok"], f"{y4['intra_worst']:.1e}", "< 1e-09"),
        ("Y4  ... and some INTER row is strictly negative",
         y4["inter_neg"], f"{y4['inter_min']:+.3e}", "< -1e-09"),
        ("Y4  the comparison is non-empty", len(y4["rows"]) >= 4,
         f"{len(y4['rows'])}", ">= 4"),
        ("Y4  feasible dephasing radius is exactly eta",
         y4["radius_linear"], f"{y4['radii'][1e-1]:.9f}", "0.100000000"),
        ("Y4  ... and collapses to a POINT at eta = 0",
         y4["radius_zero_at_zero"], f"{y4['radii'][0.0]:.1e}", "< 1e-12"),
        ("Y4  incommensurate wave arrests ONE unit at a time",
         y4["turns_ok"],
         f"{y4['turns']['incommensurate k = 2pi/phi^2'][1]}", "1"),
        ("Y4  ... and a commensurate one arrests a plane",
         y4["turns"]["commensurate  k = 2pi/4"][1] > 1,
         f"{y4['turns']['commensurate  k = 2pi/4'][1]}", "> 1"),
        # The docstring promised a tolerance "bounded from BOTH sides by the
        # measured gap" and the gaps went into a variable nothing read. They
        # are read now, and the split that produces them is tolerance-free.
        ("Y4  TURNS_TOL inside its MEASURED two-sided band",
         y4["turns_band_ok"], f"{TURNS_TOL:.0e}",
         f"{y4['turns_lo']:.1e}..{y4['turns_hi']:.1e}"),

        ("Y5  excess is HOMOGENEOUS DEGREE ONE in amplitude",
         bool(slopes) and max(abs(s - 1.0) for s in slopes) < 1e-3,
         f"{max(abs(s - 1.0) for s in slopes):.2e}" if slopes else "none",
         "< 1e-3"),
        ("Y5  ... on BOTH incommensurate ladders", y5["homog_ok"],
         f"{max(spreads):.2e}" if spreads else "none",
         f"< {TOL['homog']:.0e}"),
        ("Y5  the optimal dephasing IS a zone-corner wave",
         y5["zone_matches_lp"],
         f"{y5['zone_best'][1]:+.6e}" if y5["zone_best"][0] is not None
         else "none", f"{r8.get('val', float('nan')):+.6e}"),
        ("Y5  ... and the zone points DISAGREE (non-vacuity)",
         np.isfinite(zone_min) and np.isfinite(zone_max)
         and zone_max - zone_min > 1e-3, f"{zone_max - zone_min:.2e}",
         "> 1e-3"),
        ("Y5  ... and it alternates across jb_x's DIAGONAL pair",
         y5["zone_is_dgp"], f"{y5['zone_axes']}",
         f"{tuple(sorted(y5['dgp'][0])) if y5['dgp'] else None}"),
        ("Y3  TOL[lp] is under every decisive |value|",
         TOL["lp"] < decisive_min, f"{TOL['lp']:.0e}",
         f"< {decisive_min:.1e}"),
        ("Y5  the two ladders share no rung (independence)",
         y3["ladder_sep"] > 1e-3, f"{y3['ladder_sep']:.4f}", "> 1e-3"),
        # THREE LADDER QUANTITIES THAT WERE COMPUTED AND NEVER READ. Each backs
        # a sentence the verdict prints. The residual one was hiding a live
        # solver defect: see LSTSQ_RCONDS.
        # BOUNDED FROM BOTH SIDES, for the same reason as the diagonal-slope
        # row: a worst Gauss-Newton residual of EXACTLY zero over 240 solves is
        # not convergence, it is a routine reporting a constant, and a probe
        # that forced it reddened nothing until the lower bound existed.
        ("Y5  EVERY ladder rung closed to SOLVE_TOL",
         0.0 < y3["order_res_max"] < SOLVE_TOL,
         f"{y3['order_res_max']:.2e}", f"0 < .. < {SOLVE_TOL:.0e}"),
        ("Y5  ... the sign is constant down every ladder",
         y3["order_one_sign"] and y3["order_n"] > 0,
         f"{y3['order_n']}", "all, n > 0"),
        ("Y5  ... over at least 10 rungs in the fit window",
         y3["order_nfit_min"] >= 10, f"{y3['order_nfit_min']}", ">= 10"),
        ("Y5  TOL[homog] inside its band (measured .. 2nd order)",
         (max(spreads) if spreads else 1.0) < TOL["homog"] < 1e-3,
         f"{TOL['homog']:.0e}",
         f"{max(spreads):.1e}..1e-03" if spreads else "none"),
    ]

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
    print("   * 'BOTH classes are non-empty' -- without it, 'self-stress is the")
    print("     discriminator' is satisfied by a table with no obstructed box in")
    print("     it, which is exactly how the equivalent row in jb_x was vacuous.")
    print("   * 'the zone corner IS free' and 'a nearby control is NOT' -- the")
    print("     Bloch minimum being large is satisfiable by a residual that is")
    print("     large everywhere, including where the mode is provably free.")
    print("   * 'single-direction parity is NOT retained' -- without it, the")
    print("     affine-times-parity row is satisfied by a projector that retains")
    print("     everything, i.e. by no obstruction at all.")
    print("   * 'no odd functional fits' -- without it, 'both one-sided")
    print("     derivatives are positive' is satisfiable by a smooth function")
    print("     with a positive gradient measured only in the uphill half.")
    print("   * 'the zone points DISAGREE' -- without it, 'the optimum is a zone")
    print("     corner' is satisfied by a functional that is constant in k, which")
    print("     is precisely the intra-unit case and would be the wrong reading.")
    print("   * 'the STATIC/PATH comparison is non-empty' -- 'they agree' over")
    print("     zero rows is True.")
    print("   * 'IS basis-sensitive' -- without it, 'angle-independent' is")
    print("     satisfied by a routine returning a constant.")
    print("   * '>= 3 DISTINCT angles spanning > 10 deg' -- 'the dimension does")
    print("     not move with the phase' over ONE angle, or over the same angle")
    print("     twice, is True for any routine whatever.")
    print("   * 'the two k-grids are coprime and offset apart' and 'each")
    print("     resolves at least 15^3' -- 'measured on two incommensurate")
    print("     grids' is worth nothing if the second is the first renamed.")
    print("   * '>= 4 clusters with positive counts' -- all() over an empty")
    print("     record is True, and the whole T2 cross-check was deletable.")
    print("   * 'N really grows across the three cubes' -- 'the dimension is")
    print("     constant in N' over one cluster measured three times is not a")
    print("     statement about N.")
    print("   * 'cycles alone do NOT imply self-stress' and 'coordination alone")
    print("     does NOT either' -- without both witnesses the conjunction row")
    print("     is satisfied by either condition being the whole story.")
    print("   * 'a non-closing direction is CAUGHT' -- without it, 'every rung")
    print("     closed' is satisfiable by a residual gate nothing can trip.")
    print("   * 'the min-swapped falsifier IS negative' -- without it, 'excess")
    print("     stays POSITIVE' is covered only by moving TOL[lp], a change-")
    print("     detector shape; this asserts an object mutation (which active-")
    print("     set member each rung reads) genuinely flips the sign.")
    print("   * 'the UNTRUNCATED step does not' -- without it the least-squares")
    print("     truncation is a change no check distinguishes from its absence.")
    print("   * 'single-direction parity is NOT retained, at any size' -- the")
    print("     saturation retention row is otherwise satisfied by a projector")
    print("     that retains everything, i.e. by no obstruction at all.")
    print("   * the ZERO lower bounds on the ladder residual and on the")
    print("     diagonal-slope deviation -- a diagnostic reporting EXACTLY zero")
    print("     is a routine not measuring, and an upper bound alone accepts it.")
    print("   * 'both classes of row occur' under the infeasible-pin rows and")
    print("     under DERIVED_MARGIN -- each of those halves is about a class")
    print("     of cluster, and a class with no members is not a check.")
    print()
    print("  TWO ROWS DELETED RATHER THAN FIXED. 'f(k) is FLAT for the intra")
    print("  binder' and '... and equals |s'(a_ico)|' were driven by a routine")
    print("  that returned |s'(a)| IDENTICALLY -- for every k, grid, offset and")
    print("  angle -- so no state of the array could redden them and only a code")
    print("  edit could. A row falsifiable only by mutating the code measures")
    print("  the code. The claim they carried is a two-line derivation and is")
    print("  stated in Y5(b) instead. Precedent: the rank-2 skew row of the")
    print("  sibling file, retired for the same reason.")
    print()
    print("  A ROW DELIBERATELY NOT BUILT. 'The dephased array closes' is not a")
    print("  gate row, because for a tree it is an identity of trees and for a")
    print("  cycle it is FALSE in general -- that falsity is the Y1 finding, and")
    print("  a row asserting closure would have had to be scoped to the")
    print("  unobstructed clusters, where it could not fail. The closure residual")
    print("  is instead a MEASUREMENT inside Y0(a) and an input to Y1.")

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
    print("jb_y_dephasing -- the first OUT-OF-PHASE array model")
    print("=" * 78)
    print("  bead inviscid-qvf.14. Units joined at SINGLE VERTICES; each unit")
    print("  is a rigid placement times a phase (the DOWELED model, the only one")
    print("  in which 'phase' is a variable at all). The four standing")
    print("  declarations are INAPPLICABLE: nothing here is dynamical. METRIC")
    print("  FORM carries a qualification stated in the module docstring --")
    print("  every MAGNITUDE below is a rate per unit of a chosen norm (the")
    print("  sup-norm on phase rates in degrees, and in the mechanism row a")
    print("  sup-norm on a mixed rotation/translation rate), so magnitudes are")
    print("  not comparable across clusters in any unit. Every VERDICT drawn is")
    print("  a sign or a dimension, and those are norm-free.")

    if not PAIRS:
        print()
        print("=" * 78)
        print("GATE  1 row: the hinge pairing could not be read")
        print("=" * 78)
        print(f"  FAIL  Y0  hinge pairing readable                    "
              f"{'unreadable':>18s} {'12 x mult 2':>16s}")
        print()
        print("  Nothing below could be computed, so nothing below is printed --")
        print("  and this arrives as a FAIL ROW rather than as a traceback,")
        print("  because a traceback destroys the verdict table and leaves a")
        print("  reader unable to tell a broken build from a measured result.")
        return 1

    # SPECULATIVE PARALLEL PREFETCH: the two scans are independent.
    jb_cache.prefetch(_bloch_scan)

    y0 = y0_control()
    y1 = y1_admissible()
    y2 = y2_tied_orbit(y0)
    y3 = y3_order(y0)
    y4 = y4_static_vs_path(y0)
    y5 = y5_amplitude_wavevector(y0, y3)
    y6_verdict(y0, y1, y2, y3, y4, y5)
    return gate(y0, y1, y2, y3, y4, y5)


if __name__ == "__main__":
    # numpy's matmul on this platform's BLAS raises spurious divide/overflow
    # warnings whose text names lines that do no division at all. They are
    # suppressed so the output stays byte-identical across runs. Every place
    # where a non-finite value could actually MATTER carries an explicit
    # np.isfinite check that breaks the iteration instead of propagating, and
    # every quantity in the gate is compared against a finite threshold, so a
    # genuine nan reaches the table as a FAIL rather than as a warning nobody
    # reads.
    # `--no-cache` / `--clear-cache` are consumed here; anything else is a
    # loud failure rather than a run that silently ignored what was asked.
    _rest = jb_cache.parse_argv(sys.argv[1:])
    if _rest:
        print(f"unrecognised argument(s): {' '.join(_rest)}", file=sys.stderr)
        print("usage: jb_y_dephasing.py [--no-cache] [--clear-cache]",
              file=sys.stderr)
        sys.exit(2)
    with np.errstate(all="ignore"):
        sys.exit(main())
