"""Step X: the FIRST MULTI-UNIT jitterbug model in this project.

Everything before this file is ONE unit. This is the ARRAY: several jitterbug
units joined to their neighbours by SINGLE-VERTEX CONTACTS, which is the
join the project owner chose physically (from four options) for the wood-and-
wire array whose behaviour opened bead `inviscid-qvf.11`.

THE OBSERVATION BEING MODELLED
------------------------------
With every unit at the same phase `a` (IN PHASE), at the icosahedral phase
a = 22.238756093 the physical array LOCKS: it cannot EXPAND, it can still
CONTRACT. One-sided.

THE DIRECTIONAL FACT THAT SHAPES THE WHOLE FILE
-----------------------------------------------
Expansion is `a` DECREASING (22.24 -> 0, toward the vector equilibrium). Over
that range the circumradius RISES (0.951*L -> 1.000*L) and the folding square
diagonal RISES (1.000*L -> 1.414*L) while the struts stay constant. Expansion
LENGTHENS spans. So whatever blocks expansion must be something that MUST GET
LONGER AND CANNOT: a TENSION-ONLY member. A collision blocks APPROACH, not
separation, and is the wrong shape for this observation.

That is why the PRIMARY DELIVERABLE of this file is a SPAN ENUMERATION and not
a collision scan. Section X5 enumerates every span in the assembled array whose
length increases as `a` decreases, and ranks them by the angle at which an
inextensible member of strut length across that span would go taut.

WHAT IS MODELLED, AND WHAT IS NOT
---------------------------------
Only the VERTEX IDENTIFICATIONS are imposed. The lattice spacing is NEVER fixed
by hand -- fixing it would beg the question this file exists to answer. Each
unit carries 6 internal degrees of freedom (measured, R2 of the memo: 48 body
DOF - 36 hinge constraints - 6 rigid motions) plus 6 rigid-placement DOF; each
inter-unit identification is 3 scalar equality constraints. The solve decides
which placements are consistent.

Two mechanical models, both required by the bead:

  FREE array   -- 8 rigid triangles per unit, 48 variables per unit, 36 intra-
                  unit hinge equations per unit, 3 equations per contact.
  DOWELED array-- the owner's rig puts a DOWEL through each triangle's centre
                  forcing prescribed guide paths, which mechanically enforces
                  the symmetric 1-DOF motion per unit. Modelled EXACTLY as
                  restricting each unit's admissible motion to the 7-dimensional
                  span {6 rigid placements} + {symmetric path tangent}. So a
                  doweled unit has 7 variables, not 48, and the dowels remove
                  five of the six internal DOF by construction.

CONTROL, non-negotiable and run first (X0): this machinery applied to N = 1 must
reproduce 6 internal DOF and RANK 36, and must reproduce sigma_36 = 0.5987572 at
a = 60 -- a number recorded from the Java `JitterbugLinkage` before this file
existed. If X0 fails, nothing downstream means anything.

TOPOLOGY IS A PARAMETER, NOT AN ASSUMPTION
------------------------------------------
An earlier draft of this work assumed FCC twelve-around-one with no source.
Fuller's own array is SIX AROUND ONE -- 784.30 "Six icosahedra may be arrayed
around a nuclear icosahedron in a true XYZ-coordinate model", with 784.20/.40/.41
on unlimited periodic arrays and 724.31 on the jitterbug inside the tensegrity
icosahedron. Six-around-one is what physical builders use. Both are built here,
along with the small cycles that discriminate between them, and the file reports
what each one does rather than choosing in advance.

CAUTION carried from the corpus sweep and NOT flattened: "icosahedra cannot
aggregate" is VERIFIED for CLOSEST-PACKED SPHERES (415.50, 419.37, 987.065,
1011.35-.38, 1052.43) and REFUTED for PHYSICAL ARRAYS (784.20-.41). Different
objects. Nothing in this file is a sphere-packing claim.

WHAT THIS FILE MEASURES (the bead's questions, in order)
--------------------------------------------------------
Q1 X2  Does an in-phase configuration exist for every `a`, or does the family
       terminate? Swept, per topology, with a non-degeneracy guard.
Q2 X4  If it terminates or degenerates, is a = 22.238756093 distinguished?
Q3 X4  Is the obstruction one-sided? A rank drop is two-sided; one-sidedness
       needs a boundary/cusp or an INEQUALITY.
Q4 X3  RANK of the assembled constraint Jacobian vs `a`, free and doweled, per
       topology. RANK is reported, never a subtraction.
Q5 X6  Chirality frustration -- PRE-REGISTERED before measurement.
Q6 X7  Inter-unit clearance vs `a` (demoted; measurement only, no admissibility
       verdict -- DECISION 16 permits interference in the model).

X5 is the primary deliverable and does not correspond to a numbered question:
the span enumeration nobody could run before, because no array model existed.

CONVENTIONS INHERITED FROM THIS DIRECTORY
-----------------------------------------
Deterministic and byte-identical across runs; exit code from the gate table; no
raise inside a swept loop (a traceback destroys the verdict table -- seven of
them came out of one sentinel bug in jb_v); a check whose non-vacuity is printed
prose rather than an assertion cannot fail; a guard band must be constrained from
ABOVE as well as below. Run from the repository root with python3.

THE FOUR DECLARATIONS, and why none of them is invoked here
-----------------------------------------------------------
The bead states that the epic's four standing declarations -- interaction
KERNEL, MASS MODEL, PRIMITIVE, METRIC FORM -- "still apply to any number this
bead produces". They are DYNAMICAL declarations, and this file produces no
dynamical quantity: every number in it is KINEMATIC (ranks and singular values
of constraint Jacobians, taut angles, span lengths, clearances). No energy, no
mass, no time, no frequency, no metric on configuration space appears anywhere.
So all four are INAPPLICABLE rather than forgotten, and this sentence exists so
that the next reader can tell those two apart.

WHAT WAS PRE-REGISTERED, AND WHERE THE PRE-REGISTRATION LIVES
--------------------------------------------------------------
Q5 (chirality) was pre-registered EXTERNALLY, in the bead: `bd show
inviscid-qvf.11`. That is the record that does the work, because it is
timestamped outside this file. X6 restates the prediction for the reader's
convenience; a file vouching for its own chronology is not evidence.

DOCUMENTATION NOTE: the directory README was being edited by another agent while
this file was written, so the README entry for this script was DEFERRED. This
docstring carries it.

REVISION NOTE. This file was revised after a substantive critique and an
independent gate validation, both of which found claims here that exceeded what
was measured. Three headline sentences were RETRACTED rather than softened --
"the lock angle is unchanged by assembly" (X5c/X8, refuted by this file's own
enumeration once it is ranked by span length instead of by taut angle), "the
extra blocker belongs to an incomplete cluster / a physical array has
boundaries" (refuted by a holed cluster that has boundaries everywhere and by
CUBE27-M which has boundaries and no extra blocker), and "full row rank ...
no cusp anywhere on the swept interval" (which no row asserted). The sections
that carry the retractions say so in place.
"""
import itertools as it
import sys

import numpy as np
from scipy.optimize import brentq

from jb_a_family import corners, cluster, rot, L_EDGE, Z, faces

# ==========================================================================
# LOCAL CONSTANTS
#
# Every constant a mutation probe needs to reach is defined HERE, locally.
# jb_v's post-mortem records a probe that mutated a name the file never
# defined locally, so the mutation never applied and a clean exit was read as
# confirmation.
# ==========================================================================

#: The icosahedral phase. Recorded independently in the memo of record and in
#: `jb_a_family.__main__`; RE-DERIVED in X5 from this file's own span
#: enumeration, so it is a cross-check here rather than an input.
A_ICO = 22.238756093

#: Strut length. The brief specifies members "of strut length"; this is that
#: length, and it is the member length used for the X5 ranking.
STRUT_LEN = L_EDGE

#: Where the hinge combinatorics are read. Any angle whose configuration is 12
#: shared vertices of multiplicity 2 will do; 30 degrees is generic. Same
#: convention and same value as Java `JitterbugLinkage.PAIRING_PROBE_DEG`. The
#: pairing is read ONCE and HELD: the linkage does not re-wire itself at the
#: merge angles 60/120/240/300, where the twelve vertices become six of
#: multiplicity 4.
PAIRING_PROBE_DEG = 30.0

#: Rank tolerance, relative to the largest singular value.
RANK_RTOL = 1e-10

#: Newton/Gauss-Newton settle threshold for "a solution exists".
SOLVE_TOL = 1e-12

#: Non-degeneracy guard for the array solve, as a fraction of the STRUT LENGTH.
#: Single-vertex identification does NOT prevent two units from occupying the
#: same place -- excluding that needs an INEQUALITY, not an equality constraint
#: -- and an unguarded solver reports the total collapse (every unit centre at
#: one point) as an exact success. The threshold is ABSOLUTE rather than
#: relative to the a = 0 spacing, because the array genuinely CONTRACTS as `a`
#: rises and a spacing-relative floor would reject legitimate contraction as
#: collapse. Constrained from ABOVE as well as below in the gate: the guard
#: must reject at least one exact-residual solution somewhere in the sweep, or
#: it is a guard that never fires.
COLLAPSE_FRAC = 0.05

#: Angle grid for the span enumeration and the sweeps. Endpoints avoided: a = 0
#: and a = 60 are cone points of the hull (memo C5) and a = 60 is the vertex
#: merge, so the grid stays strictly inside.
SWEEP_LO, SWEEP_HI, SWEEP_STEP = 0.10, 59.90, 0.10

#: Fine step used as the SECOND arm of X5d's step-independence comparison. An
#: ABSOLUTE LITERAL, deliberately, and incommensurate with SWEEP_STEP. An
#: earlier version wrote this as `SWEEP_STEP * 0.37`, which made the comparison
#: a RATIO: coarsening SWEEP_STEP coarsened both arms in lock-step and the row
#: could not see it. A step-independence row whose two arms move together is
#: structurally invisible to any mutation of the step.
SWEEP_STEP_ALT = 0.037

#: Coarse grid for the two-pass span search: pass 1 finds which spans STRADDLE
#: the member length at all, pass 2 refines only those.
SPAN_COARSE_STEP = 0.5

#: The SECOND coarse step, absolute and incommensurate with SPAN_COARSE_STEP.
#: This exists because the coarse grid was the one grid in the file that was
#: never swept, and it turned out to decide a headline result. The pre-filter
#: dropped a straddling span whenever the COARSE grid happened to land on an
#: angle where that span's length is exactly zero, which is grid LANDING, not
#: grid fineness: 0.11/0.25/0.37/0.5/0.7/1.0/1.5 all gave 194 crossings for the
#: star and 0.05/0.13 gave 190, losing the sub-icosahedral angle entirely. The
#: filter is fixed (see `_span_crossings`) and the two-step comparison below is
#: what keeps it fixed. 0.13 is chosen because it is one of the two steps that
#: exhibited the bug, so the row is a REGRESSION test and not a formality.
SPAN_COARSE_ALT = 0.13

#: A span whose length varies by less than this over the whole window is
#: CONSTANT and is excluded from the taut-angle search. Struts are the obvious
#: case, but the assembled array has more of them, and they matter: a span
#: constant AT the member length oscillates around it at the last bit and a
#: naive sign-change scan reports dozens of spurious crossings. An earlier
#: draft of this file did exactly that and reported 2019 crossings where there
#: are 194. Constant spans are counted and reported, not silently dropped --
#: they are real spans, they simply cannot LENGTHEN and so cannot block
#: expansion.
#: ONE tolerance, used by BOTH constant-span filters (`_span_crossings` for the
#: assembled array and X5a for the single unit). X5a previously carried its own
#: hard-coded 1e-12, which is two thresholds for one phenomenon. The value is
#: bounded from ABOVE AND BELOW by measurement rather than chosen: X5a gates
#: `largest constant range < CONST_TOL < smallest varying range`, and the two
#: measured ends are 1.3e-15 and 3.1e-01, so the band is fourteen decades wide
#: and the constant is inside it with room to spare -- but it is now a MEASURED
#: band and a widening past 0.3 turns the row red.
CONST_TOL = 1e-9

TOL = {
    "fd_jacobian": 1e-9,       # analytic Jacobian vs central differences
    "hinge_residual": 1e-14,   # symmetric family satisfies its own hinges
    "solve": SOLVE_TOL,
    "ktable": 1e-6,            # agreement with the independently recorded k-table
    "aico": 1e-8,              # re-derived a_ico vs the recorded value
    "antipode": 1e-12,         # central symmetry of the 12 shared vertices
    "mirror": 1e-12,           # mirror image is a PROPER rotation of the unit
    "sigma36": 1e-7,           # sigma_36 at a = 60 vs the Java-era recorded value
}

#: The DELIBERATE OFFSET that bounds `TOL["aico"]` FROM ABOVE, and the worst
#: hole an independent gate validation found in the first version of this file:
#: `TOL["aico"]` was constrained from below by the measurement and from above by
#: NOTHING. Loosening it from 1e-8 to 1.0 left the output BYTE-IDENTICAL, and
#: with A_ICO also driven to 22.9 the whole gate still exited 0 -- the headline
#: number of the entire deliverable could be wrong by 0.66 degrees with every
#: row green, because A_ICO reaches the gate through exactly ONE row governed by
#: exactly ONE unbounded tolerance.
#: Closed by a CONTROL comparison: the re-derived angle must agree with A_ICO to
#: within `TOL["aico"]` AND must DISAGREE, at that same tolerance, with a value
#: offset by this much. That forces `TOL["aico"] < AICO_CONTROL_OFFSET`.
#: The band is bounded at both ends by things outside the choice: BELOW by the
#: recorded value's own quantisation (A_ICO is given to nine decimals, so
#: 5e-10 is the finest agreement the RECORD can support and a tolerance under
#: that would be testing the record's rounding), ABOVE by this offset. Both
#: ends are gated.
AICO_CONTROL_OFFSET = 1e-3

#: The quantisation of the recorded A_ICO: nine decimal places, so half a unit
#: in the last place. `TOL["aico"]` below this would be testing the record's
#: rounding rather than this file's agreement with it.
AICO_RECORD_QUANTUM = 5e-10

#: sigma_36 of the single-unit hinge Jacobian at a = 60, recorded in the memo of
#: record ("rank is constant 36 through it, sigma_36 = 0.5987572") from the Java
#: implementation. A number entering from OUTSIDE this file.
SIGMA36_AT_60 = 0.5987572

#: The k-table from T2 `nuclear-commensurability-measured-fuller-461.md`,
#: produced by a DIFFERENT script: the angle at which a tension-only member of
#: length k*L spanning a folding square diagonal goes taut. Re-derived in X5
#: from this file's own enumeration and compared. Numbers from outside.
K_TABLE = ((0.90, 26.555073), (0.95, 24.426009), (1.00, 22.238756),
           (1.05, 19.984782), (1.10, 17.653723), (1.20, 12.705717))


# ==========================================================================
# THE SINGLE UNIT, RE-EXPRESSED FOR ASSEMBLY
# ==========================================================================

def _read_pairing():
    """The 12 hinges as ((faceA, cornerA), (faceB, cornerB)), read once.

    Deliberately the same construction as Java `JitterbugLinkage`'s static
    block, including its two structural assertions. Returned in vertex-label
    order, which is the order `cluster` assigns by first appearance over the
    fixed [face][corner] scan -- deterministic, and (measured in X1) constant
    in `a` for every generic angle.
    """
    c = cluster(corners(PAIRING_PROBE_DEG), 1e-7)
    reps, mult, labels = c
    if len(reps) != 12 or sorted(set(mult.tolist())) != [2]:
        return None
    lab = labels.reshape(8, 3)
    slots = [[] for _ in range(12)]
    for i in range(8):
        for j in range(3):
            slots[lab[i, j]].append((i, j))
    return tuple(tuple(s) for s in slots), lab


_PAIRING = _read_pairing()
PAIRS = _PAIRING[0] if _PAIRING else ()
_LAB = _PAIRING[1] if _PAIRING else None

#: The 24 struts as unordered pairs of vertex labels. Constant length in `a` by
#: rigidity, so X5 must exclude them from the taut-angle ranking -- and the
#: count 24 is ASSERTED there rather than assumed.
STRUTS = (frozenset(frozenset((int(_LAB[i, j]), int(_LAB[i, (j + 1) % 3])))
                    for i in range(8) for j in range(3))
          if _LAB is not None else frozenset())


def verts(a):
    """The 12 shared vertices at phase `a`, in the held label order."""
    x = corners(a)
    return np.array([x[f][c] for (f, c), _ in PAIRS])


def dverts_exact(a):
    """d(verts)/da in DEGREES, in closed form.

    x_{f,j}(a) = R(u_f, sigma_f (a - 60)) (v - c_f) + u_f Z cos a, so

        dx/da = (pi/180) [ sigma_f  u_f x x_{f,j}  -  u_f Z sin a ]

    because u x u = 0 kills the translation's contribution to the cross
    product. Verified against central differences in X0; a sign error in the
    sigma term is exactly what that check exists to catch.
    """
    x = corners(a)
    fs = faces()
    out = np.empty((12, 3))
    k = np.pi / 180.0
    for v, ((f, c), _) in enumerate(PAIRS):
        _, _, u, sigma = fs[f]
        out[v] = k * (sigma * np.cross(u, x[f][c])
                      - u * Z * np.sin(np.radians(a)))
    return out


def _hat(r):
    return np.array([[0.0, -r[2], r[1]],
                     [r[2], 0.0, -r[0]],
                     [-r[1], r[0], 0.0]])


def _centroids(x):
    return x.mean(axis=1)


def hinge_residual(x):
    """36 components: for each hinge, the vector between its two corners."""
    return np.concatenate([x[fa][ca] - x[fb][cb] for (fa, ca), (fb, cb) in PAIRS])


def hinge_jacobian(x):
    """The 36x48 intra-unit constraint Jacobian.

    Columns 0..23 are the eight bodies' angular blocks, 24..47 the
    translational ones -- the same packing as the Java, so the two are directly
    comparable.
    """
    cen = _centroids(x)
    j = np.zeros((36, 48))
    for h, ((fa, ca), (fb, cb)) in enumerate(PAIRS):
        ra, rb = x[fa][ca] - cen[fa], x[fb][cb] - cen[fb]
        j[3 * h:3 * h + 3, 3 * fa:3 * fa + 3] += -_hat(ra)
        j[3 * h:3 * h + 3, 3 * fb:3 * fb + 3] -= -_hat(rb)
        j[3 * h:3 * h + 3, 24 + 3 * fa:27 + 3 * fa] += np.eye(3)
        j[3 * h:3 * h + 3, 24 + 3 * fb:27 + 3 * fb] -= np.eye(3)
    return j


def position_jacobian_row(x, f, c):
    """The 3x48 derivative of corner (f, c)'s position wrt this unit's motions."""
    r = x[f][c] - _centroids(x)[f]
    m = np.zeros((3, 48))
    m[:, 3 * f:3 * f + 3] = -_hat(r)
    m[:, 24 + 3 * f:27 + 3 * f] = np.eye(3)
    return m


def apply_body_motions(x0, z):
    """Exact nonlinear body motion: each body rotates about its own centroid.

    Rotation vectors are exponentiated, so this is the map whose derivative at
    z = 0 is `hinge_jacobian`, and X0 checks exactly that.
    """
    cen = _centroids(x0)
    out = np.empty_like(x0)
    for i in range(8):
        w = z[3 * i:3 * i + 3]
        th = float(np.linalg.norm(w))
        r = np.eye(3) if th < 1e-14 else rot(w / th, np.degrees(th))
        out[i] = (r @ (x0[i] - cen[i]).T).T + cen[i] + z[24 + 3 * i:27 + 3 * i]
    return out


def global_rigid_motions(x):
    """The six global rigid motions of one unit as 48-vectors."""
    cen = _centroids(x)
    g = np.zeros((6, 48))
    for d in range(3):
        for i in range(8):
            g[d][24 + 3 * i + d] = 1.0
    for d in range(3):
        e = np.zeros(3)
        e[d] = 1.0
        for i in range(8):
            g[3 + d][3 * i + d] = 1.0
            g[3 + d][24 + 3 * i:27 + 3 * i] = np.cross(e, cen[i])
    return g


def path_tangent_48(a, h=1e-6):
    """The symmetric path direction as a 48-vector of body motions.

    Built from the exact corner derivative by least squares against the
    position Jacobian: the path moves each rigid triangle, so its corner
    velocity field IS in the range of that Jacobian, and the residual of the
    fit is reported in X0 as the check that it is.
    """
    x = corners(a)
    p = np.vstack([position_jacobian_row(x, i, j).reshape(3, 48)
                   for i in range(8) for j in range(3)])
    dx = ((corners(a + h) - corners(a - h)) / (2 * h)).reshape(72)
    z, *_ = np.linalg.lstsq(p, dx, rcond=None)
    return z, float(np.linalg.norm(p @ z - dx))


def rank_of(m):
    """(rank, singular values). RANK is what this file reports, never a
    subtraction -- naive DOF counting has already failed here once (memo R2)."""
    if m.size == 0:
        return 0, np.zeros(1)
    s = np.linalg.svd(m, compute_uv=False)
    return int((s > s[0] * RANK_RTOL).sum()), s


# ==========================================================================
# TOPOLOGY
# ==========================================================================

def antipode_permutation(a):
    """For each shared vertex, the index of its antipode, plus the worst miss.

    Returns (None, inf) when the pairing is unreadable, so that a bad probe
    angle reaches the gate as a FAIL row rather than as an import-time
    traceback.
    """
    if not PAIRS:
        return None, float("inf")
    v = verts(a)
    d = np.linalg.norm(v[:, None, :] + v[None, :, :], axis=-1)
    return d.argmin(axis=1), float(d.min(axis=1).max())


_ANTI_RAW, _ANTI_DEV = antipode_permutation(PAIRING_PROBE_DEG)
ANTI = tuple(int(k) for k in _ANTI_RAW) if _ANTI_RAW is not None else ()

#: The six antipodal pairs of shared vertices, sorted for determinism.
ANTIPODAL_PAIRS = (tuple(sorted({tuple(sorted((k, ANTI[k]))) for k in range(12)}))
                   if ANTI else ())


class Topology:
    """A contact topology: named sites and a held list of vertex identifications.

    `sites(v)` maps the 12-vertex array at phase `a` to the N unit ORIGINS of
    the reference (pure-translate) placement. That placement is a STARTING
    POINT for the solver and an ANSWER only where the solver confirms it; the
    contacts are what is imposed.
    """

    def __init__(self, name, kind, gens, box=None, contacts=None, note="",
                 holes=()):
        self.name, self.kind, self.gens, self.box = name, kind, tuple(gens), box
        self.note = note
        self.holes = tuple(sorted(holes))
        if kind == "box":
            # `holes` removes lattice sites, giving a cluster with a VACANCY.
            # That is not a cosmetic variation: X5c shows the vacancy, and not
            # the boundary, is what decides whether the array carries an extra
            # inter-unit member at the icosahedral phase.
            self.lattice_sites = tuple(s for s in sorted(
                it.product(*[range(b) for b in box])) if s not in self.holes)
            idx = {s: i for i, s in enumerate(self.lattice_sites)}
            cs = []
            for s in self.lattice_sites:
                for m in range(len(gens)):
                    t = list(s)
                    t[m] += 1
                    t = tuple(t)
                    if t in idx:
                        cs.append((idx[s], gens[m], idx[t], ANTI[gens[m]]))
            self.contacts = tuple(cs)
        elif kind == "star":
            # centre plus one unit at +g and one at -g for each generator
            self.lattice_sites = None
            cs = []
            n = 1
            self._star = []
            for m, k in enumerate(gens):
                self._star.append((k, +1))
                cs.append((0, k, n, ANTI[k]))
                n += 1
                self._star.append((k, -1))
                cs.append((0, ANTI[k], n, k))
                n += 1
            self.contacts = tuple(cs)
            self.n = n
        elif kind == "explicit":
            self.lattice_sites = None
            self.contacts = tuple(contacts)
        if kind == "box":
            self.n = len(self.lattice_sites)
        elif kind == "explicit":
            self.n = max(max(c[0], c[2]) for c in self.contacts) + 1

    def dsites(self, dv):
        """d(sites)/da. The lattice BREATHES: the contact vectors 2*v_k are
        functions of the phase, so the unit origins move as `a` moves. This is
        the derivative the in-phase mode of X3b is built from, and it is why
        the lattice spacing is never fixed by hand anywhere in this file."""
        return self.sites(dv)

    def sites(self, v):
        g = np.array([2.0 * v[k] for k in self.gens])
        if self.kind == "box":
            return np.array([np.array(s, float) @ g for s in self.lattice_sites])
        if self.kind == "star":
            out = [np.zeros(3)]
            for k, sgn in self._star:
                out.append(sgn * 2.0 * v[k])
            return np.array(out)
        # explicit: FCC-style, origin plus 2*v_k for every k
        return np.vstack([np.zeros(3)] + [2.0 * v[k] for k in range(12)])


def _fcc13_contacts(v0):
    """Twelve-around-one, read at a = 0 where it closes, then HELD.

    Centre-to-neighbour: 12. Neighbour-to-neighbour: every ordered pair of
    neighbour directions whose difference is again a contact vector. Read at
    a = 0 and never re-read, for the same reason the hinge pairing is never
    re-read at a merge angle.
    """
    cs = [(0, k, 1 + k, ANTI[k]) for k in range(12)]
    for i in range(12):
        for j in range(i + 1, 12):
            d = 2.0 * v0[j] - 2.0 * v0[i]
            for m in range(12):
                if np.linalg.norm(d - 2.0 * v0[m]) < 1e-9:
                    cs.append((1 + i, m, 1 + j, ANTI[m]))
    return tuple(cs)


def build_topologies():
    v0 = verts(0.0)
    fcc = Topology("FCC13 (twelve-around-one)", "explicit", range(12),
                   contacts=_fcc13_contacts(v0),
                   note="Fuller's closest-packed reading; NOT 784.30")
    return (
        Topology("N1 (control)", "box", [0, 1, 3], box=(1, 1, 1),
                 note="the mandated control: must give rank 36, 6 internal DOF"),
        Topology("N2 (one contact)", "box", [0], box=(2,),
                 note="the minimum array"),
        Topology("CHAIN5", "box", [0], box=(5,), note="a path, no cycle"),
        Topology("SQUARE4 (one 4-cycle)", "box", [0, 1], box=(2, 2),
                 note="smallest cycle in a two-generator lattice"),
        Topology("SC7 star (six-around-one)", "star", [0, 1, 3],
                 note="Fuller 784.30 coordination, as a TREE"),
        Topology("CUBE8-M (60/60/90 basis)", "box", [0, 1, 3], box=(2, 2, 2),
                 note="six-around-one WITH cycles"),
        Topology("CUBE8-R (60/60/60 basis)", "box", [0, 1, 2], box=(2, 2, 2),
                 note="the other lattice class"),
        Topology("CUBE27-M", "box", [0, 1, 3], box=(3, 3, 3),
                 note="full coordination for the interior unit"),
        fcc,
    )


# ==========================================================================
# ASSEMBLY: the two mechanical models
# ==========================================================================

def assemble_free(a, topo, rots=None):
    """The FREE array's constraint Jacobian.

    48 variables per unit; 36 intra-unit hinge rows per unit; 3 rows per
    contact. `rots` is the per-unit orientation; None means all identity
    (the pure-translate placement).
    """
    n = topo.n
    x0 = corners(a)
    if rots is None:
        rots = [np.eye(3)] * n
    xs = [np.einsum("pq,ijq->ijp", rots[i], x0) for i in range(n)]
    big = np.zeros((36 * n + 3 * len(topo.contacts), 48 * n))
    for i in range(n):
        big[36 * i:36 * i + 36, 48 * i:48 * i + 48] = hinge_jacobian(xs[i])
    for e, (i, k, j, l) in enumerate(topo.contacts):
        (fa, ca), (fb, cb) = PAIRS[k][0], PAIRS[l][0]
        r = 36 * n + 3 * e
        big[r:r + 3, 48 * i:48 * i + 48] += position_jacobian_row(xs[i], fa, ca)
        big[r:r + 3, 48 * j:48 * j + 48] -= position_jacobian_row(xs[j], fb, cb)
    return big


def assemble_doweled(a, topo, rots=None):
    """The DOWELED array's constraint Jacobian.

    7 variables per unit -- (omega, tau, adot) -- and 3 rows per contact. The
    intra-unit hinge rows are absent because the dowel restricts motion to a
    subspace on which they are identically satisfied: the dowels ARE the
    symmetric 1-DOF sector, built in wood.
    """
    n = topo.n
    v, dv = verts(a), dverts_exact(a)
    if rots is None:
        rots = [np.eye(3)] * n
    big = np.zeros((3 * len(topo.contacts), 7 * n))
    for e, (i, k, j, l) in enumerate(topo.contacts):
        r = 3 * e
        big[r:r + 3, 7 * i:7 * i + 3] += -_hat(rots[i] @ v[k])
        big[r:r + 3, 7 * i + 3:7 * i + 6] += np.eye(3)
        big[r:r + 3, 7 * i + 6] += rots[i] @ dv[k]
        big[r:r + 3, 7 * j:7 * j + 3] -= -_hat(rots[j] @ v[l])
        big[r:r + 3, 7 * j + 3:7 * j + 6] -= np.eye(3)
        big[r:r + 3, 7 * j + 6] -= rots[j] @ dv[l]
    return big


def _rodrigues(w):
    th = float(np.linalg.norm(w))
    return np.eye(3) if th < 1e-14 else rot(w / th, np.degrees(th))


def solve_inphase(a, topo, seed, maxit=400):
    """Solve for rigid PLACEMENTS at a common phase `a`.

    Unknowns: (omega_i, t_i) per unit, unit 0 pinned (that is the global gauge,
    removed by zeroing its columns rather than by deleting variables, so the
    packing stays uniform). Equations: 3 per contact. Levenberg-damped
    Gauss-Newton with backtracking; no raise on any path.

    Returns (residual_norm, min_centre_separation, iterations).
    """
    n = topo.n
    v = verts(a)
    t0 = topo.sites(v)
    w = np.zeros((n, 3))
    t = t0.copy()
    if seed is not None:
        rng = np.random.default_rng(seed)
        w = rng.standard_normal((n, 3)) * 0.35
        t = t0 + rng.standard_normal((n, 3)) * 0.25
    w[0] = 0.0
    t[0] = 0.0
    cs = topo.contacts

    def pack(w_, t_):
        q = [_rodrigues(x) for x in w_]
        return q, [q[i] @ v.T + t_[i][:, None] for i in range(n)]

    q, pos = pack(w, t)
    res = np.concatenate([pos[i][:, k] - pos[j][:, l] for (i, k, j, l) in cs])
    used = 0
    for used in range(1, maxit + 1):
        jm = np.zeros((3 * len(cs), 6 * n))
        for e, (i, k, j, l) in enumerate(cs):
            r = 3 * e
            jm[r:r + 3, 6 * i:6 * i + 3] += -_hat(q[i] @ v[k])
            jm[r:r + 3, 6 * i + 3:6 * i + 6] += np.eye(3)
            jm[r:r + 3, 6 * j:6 * j + 3] -= -_hat(q[j] @ v[l])
            jm[r:r + 3, 6 * j + 3:6 * j + 6] -= np.eye(3)
        jm[:, 0:6] = 0.0
        with np.errstate(all="ignore"):
            dz = np.linalg.solve(jm.T @ jm + 1e-9 * np.eye(6 * n), -jm.T @ res)
        if not np.all(np.isfinite(dz)):
            break            # a stalled step, not a raise: the sweep continues
        dz = dz.reshape(n, 6)
        step, improved = 1.0, False
        for _ in range(40):
            w2, t2 = w + step * dz[:, 0:3], t + step * dz[:, 3:6]
            q2, pos2 = pack(w2, t2)
            r2 = np.concatenate([pos2[i][:, k] - pos2[j][:, l] for (i, k, j, l) in cs])
            if np.linalg.norm(r2) < np.linalg.norm(res):
                improved = True
                break
            step *= 0.5
        if not improved:
            break
        w, t, q, pos, res = w2, t2, q2, pos2, r2
        if np.linalg.norm(res) < SOLVE_TOL * 1e-2:
            break
    cen = np.array([pos[i].mean(axis=1) for i in range(n)])
    sep = min((float(np.linalg.norm(cen[i] - cen[j]))
               for i in range(n) for j in range(i + 1, n)), default=np.inf)
    return float(np.linalg.norm(res)), sep, used


def best_inphase(a, topo, seeds=(None, 0, 1, 2, 3, 4)):
    """Best over several starts, with the collapse guard applied SEPARATELY.

    Reported as two numbers, never one: the best residual found, and the best
    residual found among NON-DEGENERATE configurations. An unguarded solver
    reports the total collapse as an exact success, and that branch is real --
    it appears above roughly a = 47 for the twelve-around-one topology.
    """
    floor = COLLAPSE_FRAC * STRUT_LEN
    out = [solve_inphase(a, topo, s) for s in seeds]
    best = min(out, key=lambda r: r[0])
    nd = [r for r in out if r[1] > floor]
    bestnd = min(nd, key=lambda r: r[0]) if nd else (np.inf, 0.0, 0)
    fired = any(r[0] < SOLVE_TOL and r[1] <= floor for r in out)
    return best, bestnd, floor, fired


# ==========================================================================
# X0  CONTROL
# ==========================================================================

def x0_control():
    print("=" * 78)
    print("X0  CONTROL -- the array machinery applied to N = 1")
    print("=" * 78)
    print("  Mandated by the brief and run FIRST: if this fails, the array code")
    print("  is wrong and nothing downstream means anything. Three independent")
    print("  things are checked -- the analytic Jacobian against finite")
    print("  differences of the EXACT nonlinear residual, the rank, and one")
    print("  number recorded from the Java before this file existed.")
    print()
    ok_pair = (_PAIRING is not None and len(PAIRS) == 12)
    print(f"  hinge pairing read at a = {PAIRING_PROBE_DEG}: 12 hinges of "
          f"multiplicity 2?  {ok_pair}")
    print(f"  24 struts identified as constant-length pairs?  {len(STRUTS) == 24}"
          f"   (count {len(STRUTS)})")
    print()
    print(f"  {'a':>12s} {'rank':>6s} {'nullity':>8s} {'FD dev':>11s} "
          f"{'|C(x)|':>11s} {'rigid leak':>11s} {'tangent fit':>12s}")
    rows = []
    rng = np.random.default_rng(0)
    for a in (0.0, 5.0, A_ICO, 45.0, 60.0, 75.0):
        x = corners(a)
        j = hinge_jacobian(x)
        rk, s = rank_of(j)
        z = rng.standard_normal(48) * 1e-6
        fd = (hinge_residual(apply_body_motions(x, z))
              - hinge_residual(apply_body_motions(x, -z))) / 2.0
        dev = float(np.abs(fd - j @ z).max())
        cres = float(np.abs(hinge_residual(x)).max())
        leak = float(np.abs(j @ global_rigid_motions(x).T).max())
        tz, tfit = path_tangent_48(a)
        tleak = float(np.abs(j @ tz).max())
        rows.append((a, rk, dev, cres, leak, tfit, tleak, s))
        print(f"  {a:12.7f} {rk:6d} {48 - rk:8d} {dev:11.3e} {cres:11.3e} "
              f"{leak:11.3e} {tfit:12.3e}")
    print()
    print("  NULLITY 12 = 6 global rigid motions + 6 INTERNAL degrees of freedom.")
    print("  That decomposition is not asserted from the number: the rigid-motion")
    print("  leak column shows the six global motions are IN the kernel, so the")
    print("  remaining six are internal. Memo R2 records why the count has to be")
    print("  measured -- the naive count from 12 shared vertices gives 12, and only")
    print("  the RANK settled it.")
    print()
    s60 = rank_of(hinge_jacobian(corners(60.0)))[1]
    print(f"  sigma_36 at a = 60 : {s60[35]:.7f}   recorded (Java era): "
          f"{SIGMA36_AT_60:.7f}   dev {abs(s60[35] - SIGMA36_AT_60):.2e}")
    print("  a = 60 is where the twelve shared vertices merge into six of")
    print("  multiplicity 4. Rank is CONSTANT 36 straight through it (memo R3);")
    print("  a multiplicity jump is not a rank drop, and 'multiplicity jumps =>")
    print("  singular' was measured false in this project once already.")
    print()
    dv_e = dverts_exact(30.0)
    h = 1e-6
    dv_f = (verts(30.0 + h) - verts(30.0 - h)) / (2 * h)
    dvdev = float(np.abs(dv_e - dv_f).max())
    print(f"  closed-form dV/da vs central differences at a = 30: {dvdev:.3e}")
    print("  This is the derivative the DOWELED model runs on. A sign error in")
    print("  the sigma term would leave every doweled rank below unchanged in")
    print("  shape and wrong in fact, so it is checked here and not inferred.")
    print()
    dow = assemble_doweled(A_ICO, build_topologies()[0])
    print(f"  DOWELED N = 1: variables {dow.shape[1]}, constraint rows "
          f"{dow.shape[0]}")
    print("  7 = 6 rigid placements + 1 phase, and the dowels remove five of the")
    print("  six internal DOF BY CONSTRUCTION. Note what that is and is not: it")
    print("  is the PARAMETERISATION of the doweled model, not a measurement of")
    print("  it. N = 1 has NO contacts, so this Jacobian is 0 x 7 and any")
    print("  implementation whatever would report nullity 7 - 0 = 7. An earlier")
    print("  version GATED that subtraction, which put an arithmetic identity in")
    print("  the verdict table dressed as a control. The row is gone; what is")
    print("  gated instead is the SHAPE (0 rows, 7 columns, which a change to the")
    print("  parameterisation would break) and, in X3c, that all seven variables")
    print("  actually reach the constraints once there ARE contacts.")

    ranks_ok = all(r[1] == 36 for r in rows)
    fd_ok = max(r[2] for r in rows) < TOL["fd_jacobian"]
    cres_ok = max(r[3] for r in rows) < TOL["hinge_residual"]
    leak_ok = max(r[4] for r in rows) < 1e-9
    tan_ok = max(r[5] for r in rows) < 1e-6 and max(r[6] for r in rows) < 1e-7
    s36_ok = abs(s60[35] - SIGMA36_AT_60) < TOL["sigma36"]
    return dict(ranks_ok=ranks_ok, fd=max(r[2] for r in rows), fd_ok=fd_ok,
                cres=max(r[3] for r in rows), cres_ok=cres_ok,
                leak=max(r[4] for r in rows), leak_ok=leak_ok,
                tan=max(max(r[5] for r in rows), max(r[6] for r in rows)),
                tan_ok=tan_ok,
                s36=float(s60[35]), s36_ok=s36_ok, dvdev=dvdev,
                dv_ok=dvdev < 1e-6, pair_ok=ok_pair and len(STRUTS) == 24,
                dow_rows=int(dow.shape[0]), dow_vars=int(dow.shape[1]))


# ==========================================================================
# X1  TOPOLOGY -- what single-vertex contact actually admits
# ==========================================================================

def x1_topology():
    print()
    print("=" * 78)
    print("X1  TOPOLOGY: what a single-vertex contact admits, measured")
    print("=" * 78)
    print("  The owner joined units at SINGLE VERTICES. That is a BALL JOINT: 3")
    print("  scalar equations, no constraint whatever on relative rotation. The")
    print("  question this section settles is which unit placements can bring a")
    print("  vertex of one unit onto a vertex of another, and the answer is not")
    print("  free -- it is forced by a symmetry of the vertex set.")
    print()
    print("  CENTRAL SYMMETRY of the 12 shared vertices, and the antipode map:")
    perms, devs = [], []
    for a in (0.0, 5.0, A_ICO, 45.0, 59.0, 70.0):
        p, d = antipode_permutation(a)
        perms.append(tuple(int(x) for x in p))
        devs.append(d)
        print(f"    a = {a:10.6f}   worst |p_i + p_anti(i)| = {d:.3e}   "
              f"perm = {list(p)}")
    perm_stable = len(set(perms)) == 1
    print(f"  antipode permutation constant in a: {perm_stable}   "
          f"worst deviation {max(devs):.3e}")
    invol = all(ANTI[ANTI[k]] == k for k in range(12)) if ANTI else False
    fpf = all(ANTI[k] != k for k in range(12)) if ANTI else False
    print(f"  ANTI is an INVOLUTION: {invol}    and FIXED-POINT-FREE: {fpf}")
    print("  Both gated, and they are the non-vacuity of the two rows above")
    print("  them. Without them the pair is satisfied by the TRIVIAL SELF-")
    print("  PAIRING: an independent validation flipped the sign inside")
    print("  `antipode_permutation` so that argmin returned each vertex itself,")
    print("  and 'central symmetry' then measured |v - v| = 0.0, printed")
    print("  0.00e+00 and PASSED, while 'permutation constant in a' printed")
    print("  True. Five other rows went red; neither of these did. A row whose")
    print("  statistic is zero for the degenerate answer needs a companion that")
    print("  rejects the degenerate answer.")
    print()
    print("  CONSEQUENCE, and it removes the modelling freedom the brief warned")
    print("  about: |p - q| over all vertex pairs is maximised exactly by the")
    print("  antipodal pair, so if two units are placed with their centres")
    print("  2R apart along a vertex direction, the ONLY possible coincidence is")
    print("  vertex k of one with vertex anti(k) of the other. A translate array")
    print("  has no choice of which vertices meet.")
    print()
    v0 = verts(0.0)
    print(f"  the six antipodal pairs, and the contact vector 2*v_k they define:")
    for (k, l) in ANTIPODAL_PAIRS:
        print(f"    ({k:2d},{l:2d})  v_k(a=0) = {np.round(v0[k], 6)}   "
              f"|2 v_k| = {2 * np.linalg.norm(v0[k]):.9f}")
    print()
    print("  SIX-AROUND-ONE means picking THREE of those six pairs as lattice")
    print("  generators (Fuller 784.30, 'a true XYZ-coordinate model'). Of the")
    print("  20 triples:")
    classes = {}
    for tri in it.combinations(range(6), 3):
        ks = [ANTIPODAL_PAIRS[t][0] for t in tri]
        g = np.array([2 * v0[k] for k in ks])
        det = abs(float(np.linalg.det(g)))
        ang = tuple(sorted(round(float(np.degrees(np.arccos(
            abs(g[i] @ g[j]) / (np.linalg.norm(g[i]) * np.linalg.norm(g[j]))))), 4)
            for i, j in it.combinations(range(3), 2)))
        classes.setdefault((round(det, 6), ang), []).append(ks)
    n_ortho = 0
    for (det, ang), v in sorted(classes.items()):
        tag = "DEGENERATE (coplanar)" if det < 1e-9 else "independent"
        if all(abs(x - 90.0) < 1e-6 for x in ang):
            n_ortho += len(v)
        print(f"    |det| = {det:9.6f}  generator angles {ang}  "
              f"count {len(v):2d}  {tag}   e.g. ks = {v[0]}")
    print(f"    MUTUALLY ORTHOGONAL triples among all 20: {n_ortho}")
    print("  NOTE, and it is why no triple is 'the' XYZ one: no three shared")
    print("  vertices of a cuboctahedron are mutually PERPENDICULAR, so a")
    print("  vertex-contact array cannot be literally orthogonal. Fuller's XYZ in")
    print("  784.30 is about the tensegrity icosahedron's strut directions, not")
    print("  about vertex contacts. Both independent classes are built below.")
    print()
    print("  TWELVE-AROUND-ONE is a different demand. It needs the 12 contact")
    print("  vectors to be DIFFERENCE-CLOSED -- 2v_j - 2v_i must again be a")
    print("  contact vector -- because neighbours then touch each other too.")
    print("  Measured against the 48 ordered adjacent pairs read at a = 0 and")
    print("  HELD (the same discipline as the hinge pairing):")
    adj = []
    for i in range(12):
        for j in range(12):
            if i == j:
                continue
            d = 2 * v0[j] - 2 * v0[i]
            for m in range(12):
                if np.linalg.norm(d - 2 * v0[m]) < 1e-9:
                    adj.append((i, j, m))
                    break
    print(f"    adjacent ordered pairs at a = 0: {len(adj)}")
    defects = []
    for a in (0.0, 0.25, 0.5, 1.0, 2.0, 10.0, A_ICO, 40.0):
        v = verts(a)
        d = max(float(np.linalg.norm(2 * v[j] - 2 * v[i] - 2 * v[m]))
                for (i, j, m) in adj)
        defects.append((a, d))
        print(f"    a = {a:10.6f}   worst difference-closure defect {d:.6e}")
    lin = [(d / a) for a, d in defects if a > 0]
    print(f"    defect/a over a in (0, 2]: "
          f"{', '.join(f'{d / a:.6f}' for a, d in defects if 0 < a <= 2.0)}")
    print("    LINEAR in a from zero. The twelve-around-one array is not")
    print("    obstructed AT some angle -- it is obstructed everywhere except")
    print("    the vector equilibrium, and the obstruction opens at first order.")
    return dict(perm_stable=perm_stable, perm_dev=max(devs), n_ortho=n_ortho,
                anti_invol=invol, anti_fpf=fpf,
                n_triples=sum(len(v) for v in classes.values()),
                n_adjacent=len(adj), defect0=defects[0][1],
                defect_ico=[d for a, d in defects if a == A_ICO][0],
                lin_spread=(max(lin[:4]) - min(lin[:4])) if len(lin) >= 4 else np.inf)


# ==========================================================================
# X2  Q1 -- does an in-phase configuration exist for every a?
# ==========================================================================

def x2_existence(topos):
    print()
    print("=" * 78)
    print("X2  Q1: does an IN-PHASE configuration exist for every a?")
    print("=" * 78)
    print("  Only the vertex identifications are imposed. The lattice spacing is")
    print("  NOT fixed -- fixing it would beg the question. Every unit is at the")
    print("  same phase a; the solver looks for rigid placements that close every")
    print("  contact. Reported per topology: the best residual found, AND the best")
    print("  found among NON-DEGENERATE configurations, because single-vertex")
    print("  identification does not forbid two units occupying one place and an")
    print("  unguarded solver reports that collapse as an exact success.")
    print()
    print(f"  collapse guard: min centre separation must exceed "
          f"{COLLAPSE_FRAC * STRUT_LEN:.6f} = {COLLAPSE_FRAC} * strut length.")
    print()
    sweep = (0.0, 0.25, 1.0, 5.0, 10.0, A_ICO, 30.0, 40.0, 50.0, 55.0)
    out = {}
    guard_fired = False
    for topo in topos:
        if len(topo.contacts) == 0:
            print(f"  {topo.name}: no contacts (single unit). Skipped.")
            continue
        print(f"  {topo.name}  --  N = {topo.n}, contacts = {len(topo.contacts)}, "
              f"equations = {3 * len(topo.contacts)}, placement unknowns = "
              f"{6 * topo.n - 6}")
        print(f"    {topo.note}")
        print(f"    {'a':>11s} {'best resid':>12s} {'sep':>9s} "
              f"{'best NON-DEGEN':>15s} {'sep':>9s} {'exists?':>9s}")
        rows = []
        for a in sweep:
            b, bnd, floor, fired = best_inphase(a, topo)
            guard_fired = guard_fired or fired
            ex = bnd[0] < SOLVE_TOL
            rows.append((a, b[0], b[1], bnd[0], bnd[1], ex))
            snd = "-" if not np.isfinite(bnd[0]) else f"{bnd[0]:15.3e}"
            ssep = "-" if not np.isfinite(bnd[0]) else f"{bnd[1]:9.4f}"
            print(f"    {a:11.6f} {b[0]:12.3e} {b[1]:9.4f} {snd:>15s} "
                  f"{ssep:>9s} {'YES' if ex else 'no':>9s}")
        out[topo.name] = rows
        allex = all(r[5] for r in rows)
        someex = any(r[5] for r in rows)
        print(f"    exists at EVERY swept a: {allex}     at a = 0 only: "
              f"{rows[0][5] and not any(r[5] for r in rows[1:])}")
        if not allex and someex:
            print("    ^ the family TERMINATES. Where, and whether the icosahedral")
            print("      phase is distinguished, is X4's question.")
        print()
    print(f"  DID THE COLLAPSE GUARD EVER FIRE? {guard_fired}")
    print("  Gated. A guard that never rejects anything is not a guard, and this")
    print("  one has real work to do: above roughly a = 47 the twelve-around-one")
    print("  solver finds an EXACT solution (residual at machine precision) in")
    print("  which all thirteen units sit at one point. Single-vertex")
    print("  identification does not forbid it. Excluding it takes an inequality.")
    return out, guard_fired


# ==========================================================================
# X3  Q4 -- RANK of the assembled Jacobian, free and doweled
# ==========================================================================

def x3_rank(topos, exist):
    print()
    print("=" * 78)
    print("X3  Q4: RANK of the assembled constraint Jacobian, free and doweled")
    print("=" * 78)
    print("  RANK is reported. Never a subtraction: the naive count has already")
    print("  failed in this project (memo R2 -- 12 shared vertices 'therefore'")
    print("  1 DOF; truth 6, settled only by rank). The smallest singular value")
    print("  is printed beside it, because a rank DROP needs sigma_min -> 0 and a")
    print("  table of integers alone cannot show that it does not.")
    print()
    print("  Evaluated only where X2 found a solution -- a Jacobian at a")
    print("  configuration that does not satisfy the constraints is not a")
    print("  statement about the variety.")
    print()
    fine = np.round(np.arange(A_ICO - 2.0, A_ICO + 2.0001, 0.25), 6)
    out = {}
    for topo in topos:
        rows = exist.get(topo.name)
        solvable = topo.name == "N1 (control)" or (
            rows is not None and all(r[5] for r in rows))
        tag = "" if solvable else "   (NO in-phase solution off a = 0 -- rank " \
                                  "reported at the pure-translate placement " \
                                  "ONLY as a diagnostic, NOT as a variety fact)"
        print(f"  {topo.name}: N = {topo.n}, contacts = {len(topo.contacts)}{tag}")
        print(f"    {'a':>11s} | {'FREE rank':>10s} {'of rows':>8s} "
              f"{'nullity':>8s} {'sigma_min':>11s} | {'DOW rank':>9s} "
              f"{'of rows':>8s} {'nullity':>8s} {'sigma_min':>11s}")
        rr = []
        for a in (1.0, 10.0, A_ICO, 30.0, 45.0, 55.0):
            jf = assemble_free(a, topo)
            rf, sf = rank_of(jf)
            jd = assemble_doweled(a, topo)
            rd, sd = rank_of(jd)
            smf = float(sf[rf - 1]) if rf > 0 else np.nan
            smd = float(sd[rd - 1]) if rd > 0 else np.nan
            # DERIVED numerical-rank margin: how many times the floating-point
            # noise floor sigma_max * eps * sqrt(min(m,n)) does the smallest
            # surviving singular value exceed? A genuine rank drop drives this
            # to O(1). Nothing here is fitted.
            mf = (smf / (sf[0] * np.finfo(float).eps
                         * np.sqrt(min(jf.shape))) if rf > 0 else np.nan)
            md = (smd / (sd[0] * np.finfo(float).eps
                         * np.sqrt(min(jd.shape))) if rd > 0 else np.nan)
            rr.append((a, rf, jf.shape[0], 48 * topo.n - rf, smf,
                       rd, jd.shape[0], 7 * topo.n - rd, smd, mf, md))
            print(f"    {a:11.6f} | {rf:10d} {jf.shape[0]:8d} "
                  f"{48 * topo.n - rf:8d} {smf:11.6f} | {rd:9d} {jd.shape[0]:8d} "
                  f"{7 * topo.n - rd:8d} {smd if np.isfinite(smd) else 0.0:11.6f}")
        rk_const = len({r[1] for r in rr}) == 1 and len({r[5] for r in rr}) == 1
        # FULL ROW RANK, read off the SAME tuple entries the table printed.
        # Those entries -- `jf.shape[0]` and `jd.shape[0]` -- were written,
        # printed, and never read again. "Every assembled Jacobian is of FULL
        # ROW RANK" was the central Q4 claim and NO gate row asserted it; the
        # claim lived in a printed column nothing compared against. Measured on
        # this file: a nine-decade loosening of RANK_RTOL turns SEVENTEEN of the
        # forty-eight gated (topology, angle) rows rank-deficient, across FIVE
        # topologies, where the baseline has ZERO -- and before this row existed
        # the gate stayed green through all of it. That is not a no-op, it is an
        # UNASSERTED DETECTION: the program saw the mutation, printed it three
        # times, and asserted none of it. (An independent validation reported
        # twenty-four rows for the same probe; seventeen is what this sha
        # counts excluding FCC13, whose doweled Jacobian is rank-deficient at
        # baseline too and is excluded from the gate for that reason. The
        # DIRECTION and the topology count agree; the row count does not, and
        # the number stated here is the one this file reproduces.)
        full_row = all(r[1] == r[2] and r[5] == r[6] for r in rr)
        print(f"    rank constant over the coarse sweep: {rk_const}")
        print(f"    FULL ROW RANK at every coarse angle, free AND doweled: "
              f"{full_row}"
              f"{'' if full_row else '   <- rank < rows somewhere above'}")
        fine_free = fine_dow = rk_const
        if solvable and len(topo.contacts) > 0:
            fr = [rank_of(assemble_free(a, topo))[0] for a in fine]
            dr = [rank_of(assemble_doweled(a, topo))[0] for a in fine]
            print(f"    FINE sweep across a_ico +/- 2 deg at 0.25 deg: FREE rank "
                  f"{'constant ' + str(fr[0]) if len(set(fr)) == 1 else set(fr)}"
                  f", DOWELED rank "
                  f"{'constant ' + str(dr[0]) if len(set(dr)) == 1 else set(dr)}")
            fine_free, fine_dow = len(set(fr)) == 1, len(set(dr)) == 1
        out[topo.name] = dict(rr=rr, fine_free=fine_free, fine_dow=fine_dow,
                              coarse_const=rk_const, full_row=full_row,
                              solvable=solvable)
        print()
    allrr = [r for v in out.values() for r in v["rr"]]
    sm_free = min((r[4] for r in allrr if np.isfinite(r[4])),
                  default=float("nan"))
    sm_dow = min((r[8] for r in allrr if np.isfinite(r[8]) and r[8] > 0.0),
                 default=float("nan"))
    worst_free = min((r[9] for r in allrr if np.isfinite(r[9])),
                     default=float("nan"))
    worst_dow = min((r[10] for r in allrr if np.isfinite(r[10]) and r[10] > 0.0),
                    default=float("nan"))
    print("  THE RANK'S OWN MARGIN, over every topology and angle above:")
    print(f"    smallest surviving singular value: FREE {sm_free:.9f}, "
          f"DOWELED {sm_dow:.9f}")
    print(f"    as a multiple of the numerical noise floor "
          f"sigma_max * eps * sqrt(min(m,n)):")
    print(f"      FREE {worst_free:.3e} x,  DOWELED {worst_dow:.3e} x")
    print("  The DENOMINATOR is derived; the 1e6 gate on the multiple is not.")
    print("  What this statistic IS: a distance from the floating-point noise")
    print("  floor at the tolerance in force. What it is NOT, and was recorded")
    print("  as being: a detector of tolerance loosening. CORRECTION TO THE")
    print("  RECORD -- an earlier note said the nine-decade RANK_RTOL probe was")
    print("  'closed by a derived rank margin'. IT WAS NOT. The numerator is")
    print("  sigma[rank-1] and `rank` is itself computed with RANK_RTOL, so")
    print("  loosening the tolerance CUTS more singular values, makes")
    print("  sigma[rank-1] a LARGER one, and the margin RISES: measured")
    print("  1.049e+12 -> 1.337e+13 for FREE, an order of magnitude in the WRONG")
    print("  DIRECTION. A guard band that improves when you widen it is not a")
    print("  guard. What actually catches that probe is the FULL ROW RANK row")
    print("  and the COARSE-SWEEP CONSTANCY row above, both gated now; the")
    print("  margin is retained as a conditioning statistic and labelled as one.")
    print("  The doweled margin is the smaller of the two because at a = 55 the")
    print("  M-type lattice's own generators are approaching each other (the")
    print("  twelve vertices merge in pairs at a = 60), which is a real")
    print("  degeneracy of that basis and not of the model.")
    print()
    print("  READING, and stated no wider than it is sampled. Rank is constant")
    print("  and full over SIX coarse angles plus a 0.25-degree fine sweep on")
    print("  a_ico +/- 2, and constant rank on a neighbourhood means the solution")
    print("  set is locally a smooth manifold THERE: no branch point and no")
    print("  boundary on the sampled neighbourhood. Nothing here samples the")
    print("  whole interval and nothing here tests for a CUSP directly -- an")
    print("  earlier version claimed 'no cusp anywhere on the swept interval',")
    print("  which no row asserts, and that is withdrawn. And note what a rank")
    print("  DROP would give even if one occurred: a TWO-SIDED degeneracy. It")
    print("  could not by itself produce an obstruction that forbids expansion")
    print("  and permits contraction.")
    print()
    print("  " + "-" * 74)
    print("  X3b  THE IN-PHASE MODE IS A KERNEL DIRECTION -- the doweled model's")
    print("       whole content, asserted rather than assumed")
    print("  " + "-" * 74)
    print("  The array's own jitterbug motion is: every unit advances its phase")
    print("  at the same rate, no unit rotates, and each unit's origin follows")
    print("  the BREATHING lattice, tau_i = d(site_i)/da. If that vector is not")
    print("  in the kernel of the assembled doweled Jacobian, the array cannot")
    print("  perform the motion the owner performs by hand, and every existence")
    print("  result above would be about some other mechanism.")
    print()
    print("  This check exists because a mutation probe ZEROED the doweled")
    print("  Jacobian's phase column -- deleting the dowel's single degree of")
    print("  freedom outright -- and NOT ONE gate row noticed. The rank rows")
    print("  could not: the contact Jacobian is already of full ROW rank, so")
    print("  removing a column changes no rank at all.")
    print()
    print(f"    {'topology':34s} {'|J z|':>12s} {'|breathing|':>12s} "
          f"{'|J z| / |z|':>13s}")
    inphase, infeasible = {}, {}
    breathe_min = np.inf
    for topo in topos:
        if len(topo.contacts) == 0:
            continue
        rows = exist.get(topo.name)
        if rows is None:
            continue
        live = all(r[5] for r in rows)
        worst = 0.0
        jz_worst = 0.0
        br_here = np.inf
        for a in (5.0, A_ICO, 40.0):
            dv = dverts_exact(a)
            dt = topo.dsites(dv)
            z = np.zeros(7 * topo.n)
            for i in range(topo.n):
                z[7 * i + 3:7 * i + 6] = dt[i]
                z[7 * i + 6] = 1.0
            jz = float(np.linalg.norm(assemble_doweled(a, topo) @ z))
            zn = float(np.linalg.norm(z))
            worst = max(worst, jz / zn)
            jz_worst = max(jz_worst, jz)
            # The BREATHING part: the translation block alone. This, and not
            # |z|, is the non-vacuity that matters. |z| >= sqrt(n) by
            # construction because every unit carries 1.0 in its phase slot, so
            # a row asserting |z| > 1 cannot fail for any topology built here --
            # and the earlier version also read it from the LAST topology
            # (FCC13) rather than the minimum, because the accumulator was reset
            # inside the loop. The breathing norm CAN be zero: it is zero
            # exactly when the lattice is held rigid while the units move, which
            # is the modelling error this whole file exists to avoid ("the
            # lattice spacing is never fixed by hand").
            br_here = min(br_here,
                          float(np.linalg.norm(z.reshape(topo.n, 7)[:, 3:6])))
        breathe_min = min(breathe_min, br_here)
        print(f"    {topo.name:34s} {jz_worst:12.3e} {br_here:12.3e} "
              f"{worst:13.3e}   {'' if live else '<- NO in-phase solution'}")
        (inphase if live else infeasible)[topo.name] = worst
    worst_inphase = max(inphase.values()) if inphase else float("inf")
    best_infeasible = min(infeasible.values()) if infeasible else 0.0
    breathe_floor = float(breathe_min)
    print(f"  worst relative residual among SOLVABLE topologies: "
          f"{worst_inphase:.3e}")
    print(f"  best  relative residual among INFEASIBLE ones:     "
          f"{best_infeasible:.3e}")
    print(f"  smallest BREATHING norm over every topology and angle: "
          f"{breathe_floor:.3e}")
    print("  The second line is the non-vacuity of the first, and it is not a")
    print("  contrivance: twelve-around-one has no in-phase solution off the")
    print("  vector equilibrium, so the in-phase mode is NOT in its kernel, and")
    print("  the same check that passes at 1e-17 for every feasible topology")
    print("  returns O(0.1) for the infeasible one. A build in which the mode")
    print("  were trivially annihilated would fail THAT row.")
    print("  The third line is the other non-vacuity, and it replaces a row that")
    print("  could not fail: the mode's breathing component is bounded below, so")
    print("  the mode is not a pure phase advance with a frozen lattice.")
    print()
    print("  " + "-" * 74)
    print("  X3c  THE ROTATIONAL COLUMNS, WHICH X3b CANNOT SEE")
    print("  " + "-" * 74)
    print("  X3b closed a hole of the form 'a deleted COLUMN is invisible to a")
    print("  full-ROW-rank matrix'. It closed the INSTANCE and not the CLASS: an")
    print("  independent validation then zeroed each unit's three ROTATION")
    print("  columns -- asserting that a unit may rotate freely with no effect at")
    print("  its joints, three of its seven doweled variables gone -- and not one")
    print("  gate row moved. X3b cannot see it, because the in-phase mode it")
    print("  tests has NO rotation component, so those columns are multiplied by")
    print("  zero either way. Three of the seven doweled variables per unit were")
    print("  undefended.")
    print("  Closed here by a motion that USES them: a GLOBAL RIGID ROTATION of")
    print("  the whole assembly, omega_i = omega for every unit with each origin")
    print("  carried along as tau_i = omega x s_i. That is a rigid motion of the")
    print("  assembled array, so it must lie in the kernel; and it is nonzero in")
    print("  exactly the columns X3b leaves untouched.")
    print()
    print(f"    {'topology':34s} {'|J z|/|z|':>12s} {'CONTROL':>12s} "
          f"{'cols zeroed':>12s}")
    om = np.array([0.3, -0.5, 0.8])
    # SENTINELS, chosen so that "nothing ran" FAILS rather than passes. `rot_ok`
    # starts NEGATIVE and is gated from below as well as above: a max()
    # accumulator initialised at 0.0 would report 0.0 for an empty loop and sail
    # through a "< 1e-9" row. The inf-initialised minima fail their "> 0.1"
    # rows only because the gate tests isfinite first; both halves are needed.
    rot_ok, rot_ctrl, rot_bad, colmin = -1.0, np.inf, np.inf, np.inf
    for topo in topos:
        if len(topo.contacts) == 0:
            continue
        rows = exist.get(topo.name)
        if rows is None:
            continue
        live = all(r[5] for r in rows)
        v = verts(A_ICO)
        s = topo.sites(v)
        z = np.zeros(7 * topo.n)
        for i in range(topo.n):
            z[7 * i:7 * i + 3] = om
            z[7 * i + 3:7 * i + 6] = np.cross(om, s[i])
        jm = assemble_doweled(A_ICO, topo)
        zn = float(np.linalg.norm(z))
        r_ok = float(np.linalg.norm(jm @ z)) / zn
        # CONTROL: the same rotation WITHOUT carrying the origins. Not a rigid
        # motion, must NOT be in the kernel. Without this the row would be
        # satisfied by a Jacobian that annihilates everything.
        z2 = z.copy()
        z2.reshape(topo.n, 7)[:, 3:6] = 0.0
        r_ct = float(np.linalg.norm(jm @ z2)) / float(np.linalg.norm(z2))
        # And the mutation itself, run in-line so the defence is demonstrated
        # rather than described: zero the rotation columns and re-measure.
        jm2 = jm.copy()
        for i in range(topo.n):
            jm2[:, 7 * i:7 * i + 3] = 0.0
        r_bd = float(np.linalg.norm(jm2 @ z)) / zn
        # the cheap general companion: no doweled column may be identically zero
        colmin = min(colmin, float(np.abs(jm).max(axis=0).min()))
        print(f"    {topo.name:34s} {r_ok:12.3e} {r_ct:12.3e} {r_bd:12.3e}"
              f"   {'' if live else '<- NO in-phase solution'}")
        if live:
            rot_ok = max(rot_ok, r_ok)
            rot_ctrl = min(rot_ctrl, r_ct)
            rot_bad = min(rot_bad, r_bd)
    print(f"  worst over solvable topologies: in kernel {rot_ok:.3e}; control "
          f"{rot_ctrl:.3e}; with the rotation columns zeroed {rot_bad:.3e}")
    print(f"  smallest column-infinity-norm over every doweled Jacobian: "
          f"{colmin:.3e}")
    print("  Three gated rows, and the third is deliberately the cheap general")
    print("  one: any doweled column that is identically zero is a variable the")
    print("  model claims to have and does not. It is weaker than the kernel")
    print("  membership above and it catches the whole deletion class, not one")
    print("  instance of it -- which is the lesson the first fix missed.")
    return (out, worst_free, worst_dow, worst_inphase, breathe_floor,
            best_infeasible, rot_ok, rot_ctrl, rot_bad, colmin)


# ==========================================================================
# X4  Q2 and Q3 -- is a_ico distinguished, and is anything one-sided?
# ==========================================================================

def x4_distinguished(topos, exist):
    print()
    print("=" * 78)
    print("X4  Q2/Q3: is a = 22.238756093 distinguished, and is anything ONE-SIDED?")
    print("=" * 78)
    print("  Q2 asks where the in-phase family terminates or degenerates and")
    print("  whether the icosahedral phase is that place. Q3 asks whether any")
    print("  obstruction found is one-sided -- and one-sidedness cannot come from")
    print("  a rank drop, which is symmetric in the two directions. It needs a")
    print("  boundary/cusp of the variety, or an INEQUALITY.")
    print()
    print("  TERMINATION, per topology, from X2's sweep:")
    term = {}
    for topo in topos:
        rows = exist.get(topo.name)
        if rows is None:
            continue
        good = [r[0] for r in rows if r[5]]
        bad = [r[0] for r in rows if not r[5]]
        if not bad:
            verdict = "exists at every swept a -- no termination"
        elif good == [0.0]:
            verdict = "exists ONLY at a = 0 (the vector equilibrium)"
        else:
            verdict = f"exists on {good}, fails on {bad}"
        term[topo.name] = (tuple(good), tuple(bad))
        print(f"    {topo.name:34s}  {verdict}")
    print()
    print("  TWO-SIDED CONTINUATION TEST at the icosahedral phase. If the family")
    print("  had a boundary there, the branch would continue in ONE direction")
    print("  only. Solved on both sides at three step sizes:")
    print("  Split by whether the topology has a solution AT the icosahedral")
    print("  phase at all -- for one that does not, 'does the branch continue")
    print("  both ways' is not a question about one-sidedness, and mixing the")
    print("  two into one verdict would let an ABSENT family read as a boundary.")
    print(f"    {'topology':34s} {'delta':>8s} {'resid a-d':>11s} "
          f"{'resid a+d':>11s} {'verdict':>14s}")
    twoside, absent_both = {}, {}
    for topo in topos:
        if len(topo.contacts) == 0:
            continue
        rows = exist.get(topo.name)
        if rows is None:
            continue
        here = [r for r in rows if r[0] == A_ICO]
        live = bool(here) and here[0][5]
        oks = []
        for d in (1.0, 0.1, 0.01):
            lo = best_inphase(A_ICO - d, topo)[1][0]
            hi = best_inphase(A_ICO + d, topo)[1][0]
            if live:
                v = (lo < SOLVE_TOL) and (hi < SOLVE_TOL)
                tag = "BOTH SIDES" if v else "ONE SIDE ONLY"
            else:
                v = (lo >= SOLVE_TOL) and (hi >= SOLVE_TOL)
                tag = "absent BOTH" if v else "absent ONE side"
            oks.append(v)
            print(f"    {topo.name:34s} {d:8.2f} {lo:11.3e} {hi:11.3e} "
                  f"{tag:>14s}")
        (twoside if live else absent_both)[topo.name] = all(oks)
    print()
    print("  The twelve-around-one row is the informative one in the second")
    print("  group: its obstruction is present on BOTH sides of the icosahedral")
    print("  phase, so it is a TWO-SIDED absence and not a one-sided stop. It is")
    print("  also not AT that phase -- X2 puts it everywhere off a = 0.")
    print()
    print("  VERDICT ON Q2/Q3 FROM EQUALITY CONSTRAINTS ALONE: recorded in the")
    print("  gate. The exact angle, wherever anything happens, is printed to nine")
    print("  decimals and is never rounded toward 22.238756093.")
    print()
    print("  WHERE ONE-SIDEDNESS CAN COME FROM, stated before X5 measures it:")
    print("  an inequality. A tension-only member of length m across a span s")
    print("  imposes s(a) <= m, not s(a) = m. Where s is DECREASING in a, that")
    print("  inequality reads a >= a*, a HALF-LINE: expansion blocked, contraction")
    print("  free. Fuller's only genuinely one-sided constraint (711.32, a")
    print("  tensegrity sphere that 'cannot get bigger than its discretely")
    print("  designed dimensions' and 'can only relax inwardly') has exactly that")
    print("  form, and it is attached to a TENSION NETWORK. Meanwhile 541.19 has")
    print("  the bare linkage 'turn around at the negative tetrahedron to reexpand")
    print("  therefrom', and 905.40/905.55/938.13 call the icosahedral phase")
    print("  UNSTABLE, not locked. The corpus points away from a variety boundary.")
    return term, twoside, absent_both


# ==========================================================================
# X5  THE PRIMARY DELIVERABLE -- span enumeration and taut-angle ranking
# ==========================================================================

def folding_diagonals():
    """The six folding square diagonals, derived STRUCTURALLY, not read off.

    At a = 0 the twelve shared vertices are a cuboctahedron whose six square
    faces have diagonals of length sqrt(2) * strut. Half of those twelve
    diagonals SHORTEN as `a` rises and half LENGTHEN; the shortening six are
    the ones that reach strut length and make the icosahedron's 30 equal
    edges. Selected here by measurement of the a = 0 length and of the sign of
    the derivative, so the set is independent of the crossing search that X5a
    runs -- and the two are then compared, which is the point of deriving it
    twice.
    """
    if not PAIRS:
        return (), ()
    v0, v1, v2 = verts(0.0), verts(1.0), verts(2.0)
    target = np.sqrt(2.0) * STRUT_LEN
    at_target, out = [], []
    for k in range(12):
        for l in range(k + 1, 12):
            if frozenset((k, l)) in STRUTS:
                continue
            if abs(float(np.linalg.norm(v0[k] - v0[l])) - target) > 1e-9:
                continue
            at_target.append((k, l))
            if np.linalg.norm(v2[k] - v2[l]) < np.linalg.norm(v1[k] - v1[l]):
                out.append((k, l))
    return tuple(sorted(out)), tuple(sorted(at_target))


DIAGONALS, SQUARE_DIAGONALS = folding_diagonals()

#: The distance between the CENTRES of two units in contact, at the icosahedral
#: phase: |2 v_k|, the same for all twelve shared vertices. DERIVED, not chosen.
#: X5c uses it as the INSTALLABILITY scale -- a span shorter than this at a_ico
#: runs between two points closer together than two contacting unit centres, so
#: a member across it is a plausible build member.
CONTACT_SPACING_ICO = (2.0 * float(np.linalg.norm(verts(A_ICO)[0]))
                       if PAIRS else float("nan"))


def diagonal_generator_pairs(gens):
    """Which pairs of a basis's generator VERTICES are folding diagonals.

    This is the whole of the basis dependence in X5c. When generators g_i and
    g_j are an antipodal-diagonal pair, the neighbour at 2v_i carries vertex j
    at v_j + 2v_i and the neighbour at 2v_j carries vertex i at v_i + 2v_j, so
    the chord between them is v_j - v_i -- a folding diagonal, reaching strut
    length at exactly the icosahedral phase. The M basis (0,1,3) has exactly one
    such pair; the R basis (0,1,2) has none.
    """
    out = []
    for i in range(len(gens)):
        for j in range(i + 1, len(gens)):
            k, l = int(gens[i]), int(gens[j])
            if (min(k, l), max(k, l)) in DIAGONALS:
                out.append((i, j))
    return tuple(out)


def mediating_site_audit(box):
    """In a FULL rectangular box, is the mediating site ALWAYS present?

    Exhaustive over every site and every pair of generator directions: whenever
    c+e_i and c+e_j are both in the box, is c+e_i+e_j? Returns (pairs, missing).

    This exists because the gate row "the full lattices have 0 spans taut at
    a_ico" was reported as a MEASUREMENT and is in fact forced by the shape of a
    box: the mediating unit is present in every case, so the row can report
    nothing but zero. Running the audit turns that from an unexamined zero into
    a stated identity -- and X5c then gates a cluster where the site CAN be
    absent, which is a row that can actually come back nonzero.
    """
    sites = set(it.product(*[range(b) for b in box]))
    npair = miss = 0
    d = len(box)
    for c in sites:
        for i in range(d):
            for j in range(i + 1, d):
                a, b = list(c), list(c)
                a[i] += 1
                b[j] += 1
                m = list(c)
                m[i] += 1
                m[j] += 1
                if tuple(a) in sites and tuple(b) in sites:
                    npair += 1
                    if tuple(m) not in sites:
                        miss += 1
    return npair, miss


def vacancy_topologies():
    """Clusters built to make the X5c inter-unit row FALSIFIABLE.

    Every full rectangular box contains its mediating sites (see
    `mediating_site_audit`), so a row reading "the full lattices have 0" cannot
    come back with anything else. These four can:

      L-TROMINO-M   three units at c, c+e_i, c+e_j on the M basis's DIAGONAL
                    generator pair, with the mediating site c+e_i+e_j ABSENT.
      QUAD-M        the same three plus the mediating unit. The completion.
      HOLED8-M      CUBE8-M with the single mediating site (1,0,1) removed --
                    a bounded, cycle-carrying, seven-unit cluster.
      L-TROMINO-R   the same L on the R basis, whose generator pairs include no
                    diagonal. The basis control.

    All four are FINITE and have BOUNDARIES EVERYWHERE, which is the point:
    CUBE27-M also has boundaries and reports zero, so "boundary" was never the
    discriminator. The vacancy is.
    """
    m, r = (0, 1, 3), (0, 1, 2)
    dp = diagonal_generator_pairs(m)
    if not dp:
        return ()
    i, j = dp[0]
    def cell(*axes):
        c = [0, 0, 0]
        for ax in axes:
            c[ax] = 1
        return tuple(c)
    keep_l = {cell(), cell(i), cell(j)}
    keep_q = keep_l | {cell(i, j)}
    box = (2, 2, 2)
    allsites = set(it.product(range(2), repeat=3))
    return (
        Topology("L-TROMINO-M (mediator ABSENT)", "box", m, box=box,
                 holes=allsites - keep_l,
                 note="three units on the diagonal generator pair, no mediator"),
        Topology("QUAD-M (mediator PRESENT)", "box", m, box=box,
                 holes=allsites - keep_q,
                 note="the same three units plus the mediating unit"),
        Topology("HOLED8-M (one vacancy)", "box", m, box=box,
                 holes={cell(i, j)},
                 note="CUBE8-M minus the single mediating site"),
        Topology("L-TROMINO-R (basis control)", "box", r, box=box,
                 holes=allsites - keep_l,
                 note="the same shape on a basis with no diagonal generator pair"),
    )


def _classes(topo):
    """Union-find the contact identifications: the array's WIRED POINTS.

    Two vertices joined by a contact ARE one point, so a span must be measured
    between point CLASSES and not between (unit, vertex) labels. Skipping this
    double-counts: a span from unit A's vertex k to unit B's vertex l, where
    the endpoint in B is identified with a vertex of A, is not an inter-unit
    span at all -- it is one of A's own spans seen from the far side of a
    shared vertex. An earlier draft of this section reported sixteen such
    shadows as new inter-unit blockers.
    """
    n = topo.n
    par = list(range(12 * n))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for (i, k, j, l) in topo.contacts:
        ra, rb = find(12 * i + k), find(12 * j + l)
        if ra != rb:
            par[ra] = rb
    roots = sorted({find(x) for x in range(12 * n)})
    ridx = {r: m for m, r in enumerate(roots)}
    members = [[] for _ in roots]
    for x in range(12 * n):
        members[ridx[find(x)]].append((x // 12, x % 12))
    units = [frozenset(u for u, _ in mem) for mem in members]
    return members, units


def _class_positions(a, topo, members):
    v = verts(a)
    s = topo.sites(v)
    return np.array([v[mem[0][1]] + s[mem[0][0]] for mem in members])


def _span_crossings(topo, member, lo=SWEEP_LO, hi=SWEEP_HI,
                    coarse=SPAN_COARSE_STEP, fine=SWEEP_STEP):
    """Every span of the assembled array that reaches `member` while LENGTHENING
    under expansion (that is, while DECREASING in `a`).

    Two passes so that a 27-unit cluster is affordable: pass 1 keeps only the
    running min/max per span and selects the straddlers; pass 2 stores just
    those and refines each sign change with brentq. No raise on any path.

    THE STRADDLE FILTER, and a bug it carried. Pass 1 must exclude spans that
    are IDENTICALLY ZERO -- pairs of wired points that coincide at every angle.
    It used to do that with `lo_arr > 1e-9`, where `lo_arr` is the minimum over
    the COARSE grid. That is the wrong test: it also throws away a span that is
    TRANSIENTLY zero, passing through zero at one angle and growing on either
    side, and it does so only when the coarse grid happens to LAND on that
    angle. Measured on the star: coarse 0.11/0.25/0.37/0.5/0.7/1.0/1.5 gave 194
    crossings and coarse 0.05/0.13 gave 190, the four missing spans being
    exactly the ones carrying the sub-icosahedral taut angle -- a headline
    result that existed or not depending on where the grid points fell. This is
    the THIRD instance in this file of "a commensurate grid misses what lands on
    it", after the constant-span oscillation and the root at exactly a = 30, and
    the first one in a filter rather than in a root test.
    The correct test is the one already used for `n_const`: a span that is
    identically zero does not VARY, so `varies` excludes it on its own and the
    `lo_arr` term was never needed. Removing it makes the count 194 at every
    coarse step tried from 0.019 to 1.5. X5d gates that with a second,
    incommensurate coarse step.
    """
    members, units = _classes(topo)
    nc = len(members)
    iu = np.triu_indices(nc, 1)
    cg = np.round(np.arange(lo, hi + 1e-9, coarse), 6)
    lo_arr = np.full(iu[0].size, np.inf)
    hi_arr = np.full(iu[0].size, -np.inf)
    max_slope = 0.0
    prev = None
    for a in cg:
        p = _class_positions(a, topo, members)
        d = np.linalg.norm(p[iu[0]] - p[iu[1]], axis=-1)
        lo_arr = np.minimum(lo_arr, d)
        hi_arr = np.maximum(hi_arr, d)
        if prev is not None:
            max_slope = max(max_slope, float(np.abs(d - prev).max()) / coarse)
        prev = d
    rng_arr = hi_arr - lo_arr
    varies = rng_arr > CONST_TOL
    n_const = int(np.count_nonzero(~varies & (lo_arr > 1e-9)))
    n_const_at_member = int(np.count_nonzero(
        ~varies & (np.abs(lo_arr - member) < CONST_TOL)))
    # `varies` alone: see the docstring. A span identically zero cannot vary, so
    # the old `lo_arr > 1e-9` term added nothing except a grid-dependent loss.
    strad = np.nonzero((lo_arr < member) & (hi_arr > member) & varies)[0]
    fg = np.round(np.arange(lo, hi + 1e-9, fine), 6)
    # CONST_TOL's two-sided band, MEASURED on this cluster rather than assumed:
    # the largest range among spans called constant, and the smallest among
    # spans called varying. The tolerance must lie strictly between them, and
    # the gap is what makes the classification robust.
    c_hi = float(rng_arr[~varies].max()) if np.any(~varies) else 0.0
    v_lo = float(rng_arr[varies].min()) if np.any(varies) else float("inf")
    stats = dict(nc=nc, npairs=int(iu[0].size), nstrad=int(strad.size),
                 max_slope=max_slope, coarse_risk=max_slope * coarse,
                 n_const=n_const, n_const_at_member=n_const_at_member,
                 const_hi=c_hi, vary_lo=v_lo, n_brentq_fail=0, n_interval=0)
    if strad.size == 0:
        return [], stats
    col = np.empty((fg.size, strad.size))
    for t, a in enumerate(fg):
        p = _class_positions(a, topo, members)
        col[t] = np.linalg.norm(p[iu[0][strad]] - p[iu[1][strad]], axis=-1)

    def dist(a, pi, qi):
        p = _class_positions(a, topo, members)
        return float(np.linalg.norm(p[pi] - p[qi]))

    rows, n_interval, n_fail = [], 0, 0
    for c in range(strad.size):
        pi, qi = int(iu[0][strad[c]]), int(iu[1][strad[c]])
        arr = col[:, c]
        here = []
        for t in range(fg.size - 1):
            # HALF-OPEN sign convention: a value of exactly `member` at a grid
            # point counts as NEGATIVE, so a root sitting on a grid point is
            # attributed to exactly one bracket and is neither missed nor
            # double counted. The obvious test (product < 0) SKIPS it, and the
            # star cluster has a real crossing at exactly a = 30, so the step
            # sweep in X5d disagreed with itself until this was fixed.
            if (arr[t] > member) == (arr[t + 1] > member):
                continue
            try:
                astar = brentq(lambda a: dist(a, pi, qi) - member,
                               fg[t], fg[t + 1], xtol=1e-13)
            except Exception:
                # A dropped row, never a raise -- but COUNTED. An excepting
                # refinement that is silently `continue`d is data loss with no
                # trace; the count is returned and gated at zero, so the
                # swallow cannot grow quietly.
                n_fail += 1
                continue
            # DECREASING in a == lengthening under expansion == the member
            # blocks EXPANSION (a >= a*). INCREASING blocks CONTRACTION.
            here.append((float(astar), "dec" if arr[t + 1] < arr[t] else "inc"))
        if any(d == "dec" for _, d in here) and any(d == "inc" for _, d in here):
            n_interval += 1
        shared = bool(units[pi] & units[qi])
        for astar, direc in here:
            rows.append((astar, pi, qi, "intra" if shared else "inter", direc))
    # Rounded angle first in the sort key: two runs at different step sizes
    # bracket the same root differently and brentq lands on it at the last bit,
    # which would otherwise reorder rows that are equal to nine decimals and
    # make the printed order look step-dependent when the values are not.
    rows.sort(key=lambda r: (-round(r[0], 9), r[1], r[2]))
    stats["n_interval"] = n_interval
    stats["n_brentq_fail"] = n_fail
    return rows, stats


def x5_spans(topos):
    print()
    print("=" * 78)
    print("X5  THE PRIMARY DELIVERABLE: every span that LENGTHENS under expansion")
    print("=" * 78)
    print("  Expansion is a DECREASING. A span whose length increases as a")
    print("  decreases is a span an inextensible member can hold, and a member")
    print("  across it imposes span(a) <= member: the half-line a >= a*.")
    print("  Enumerated here intra-unit and -- for the first time, because no")
    print("  array model existed -- INTER-UNIT, and ranked by the angle at which")
    print("  a member of STRUT LENGTH goes taut.")
    print()
    print(f"  member length for the ranking: {STRUT_LEN:.12f} = the strut")
    print(f"  angle window [{SWEEP_LO}, {SWEEP_HI}], coarse {SPAN_COARSE_STEP}, "
          f"fine {SWEEP_STEP}")
    print()

    # ---- X5a intra-unit, topology independent ---------------------------
    print("  " + "-" * 74)
    print("  X5a  INTRA-UNIT spans: 66 vertex pairs, 24 of them struts")
    print("  " + "-" * 74)
    grid = np.round(np.arange(SWEEP_LO, SWEEP_HI + 1e-9, SWEEP_STEP), 6)
    vg = [verts(a) for a in grid]
    n_const, pairs, arrs = 0, [], []
    const_hi, vary_lo = 0.0, float("inf")
    for k in range(12):
        for l in range(k + 1, 12):
            arr = np.array([float(np.linalg.norm(v[k] - v[l])) for v in vg])
            rng = float(arr.max() - arr.min())
            # ONE named tolerance, shared with `_span_crossings`. This test used
            # to carry its own hard-coded 1e-12 -- two thresholds for one
            # phenomenon, and the second one unnamed, so no probe could reach it.
            if rng < CONST_TOL:
                n_const += 1
                const_hi = max(const_hi, rng)
                continue
            vary_lo = min(vary_lo, rng)
            pairs.append((k, l))
            arrs.append(arr)
    print(f"  spans of CONSTANT length: {n_const}   (must be exactly the 24")
    print("  struts -- gated, because a labelling error there silently changes")
    print("  which spans the ranking scans)")
    print(f"  CONST_TOL = {CONST_TOL:.0e} is bounded from BOTH sides by this")
    print(f"  measurement, not chosen: largest CONSTANT range {const_hi:.3e} <")
    print(f"  {CONST_TOL:.0e} < smallest VARYING range {vary_lo:.3e}. Both ends")
    print("  are gated, so the tolerance can be neither tightened into the noise")
    print("  nor widened until it swallows a real span. A guard band constrained")
    print("  from below only is a guard band that can be widened for free.")
    n_intra_fail = [0]

    def intra_cross(member):
        out = []
        for (k, l), arr in zip(pairs, arrs):
            for t in range(len(grid) - 1):
                if (arr[t] > member) == (arr[t + 1] > member):
                    continue          # half-open convention, see _span_crossings
                if arr[t + 1] >= arr[t]:
                    continue
                try:
                    astar = brentq(
                        lambda a, k=k, l=l: float(np.linalg.norm(
                            verts(a)[k] - verts(a)[l])) - member,
                        grid[t], grid[t + 1], xtol=1e-13)
                except Exception:
                    n_intra_fail[0] += 1   # counted, gated at zero, never a raise
                    continue
                out.append((float(astar), (k, l)))
        out.sort(key=lambda r: (-round(r[0], 9), r[1]))
        return out

    intra_rows = intra_cross(STRUT_LEN)
    print(f"  intra-unit spans reaching strut length while lengthening under")
    print(f"  expansion: {len(intra_rows)}")
    for astar, (k, l) in intra_rows:
        print(f"    a* = {astar:.9f}   vertices ({k:2d},{l:2d})   "
              f"is a derived folding diagonal: {(k, l) in DIAGONALS}")
    distinct = sorted({round(r[0], 9) for r in intra_rows}, reverse=True)
    print(f"  distinct taut angles: {distinct}")
    print()
    print(f"  THE SAME SIX, DERIVED A SECOND WAY without any crossing search:")
    print(f"    cuboctahedron square-face diagonals at a = 0 (length "
          f"{np.sqrt(2) * STRUT_LEN:.9f}): {len(SQUARE_DIAGONALS)}")
    print(f"    of those, the ones that SHORTEN as a rises: {len(DIAGONALS)}")
    print(f"    {list(DIAGONALS)}")
    same = tuple(sorted(k for _, k in intra_rows)) == tuple(sorted(DIAGONALS))
    print(f"    identical to the crossing search's answer: {same}")
    print("  Two derivations, one comparison. The crossing search could have")
    print("  returned a different set -- 36 non-strut pairs were scanned and any")
    print("  of them could have crossed. The ABSENCE of competitors is the")
    print("  measurement here, and it is what the coordinator's 'ruler test'")
    print("  could not supply: that test compared the diagonals to the struts at")
    print("  the lock, which is true by construction (omnitriangulation and")
    print("  diagonal-reaches-strut-length are the same instant) and therefore")
    print("  discriminated nothing. It confirmed the ANGLE only.")
    ico_from_spans = distinct[0] if distinct else float("nan")
    aico_dev = abs(ico_from_spans - A_ICO)
    aico_ctrl = abs(ico_from_spans - (A_ICO + AICO_CONTROL_OFFSET))
    print()
    print(f"  A_ICO, AGREEMENT AND CONTROL. The re-derived angle is")
    print(f"  {ico_from_spans:.9f}; the recorded A_ICO is {A_ICO:.9f}; the")
    print(f"  deviation is {aico_dev:.3e}, gated below TOL['aico'] = "
          f"{TOL['aico']:.0e}.")
    print(f"  AND THE CONTROL, which is what bounds that tolerance from ABOVE:")
    print(f"  the same comparison against a value offset by "
          f"{AICO_CONTROL_OFFSET:.0e} deg")
    print(f"  gives {aico_ctrl:.3e} and must EXCEED the tolerance. Without this")
    print("  row TOL['aico'] was unbounded above: an independent validation")
    print("  loosened it from 1e-8 to 1.0 and the whole output stayed BYTE-")
    print("  IDENTICAL, and with A_ICO also moved to 22.9 the gate still exited")
    print("  0 -- the headline number of this deliverable could be wrong by 0.66")
    print("  degrees with every row green, because A_ICO reached the gate")
    print("  through one row governed by one unbounded tolerance. The band is")
    print(f"  now closed at both ends: {AICO_RECORD_QUANTUM:.0e} (the recorded")
    print(f"  value's own nine-decimal quantum) < TOL['aico'] < "
          f"{AICO_CONTROL_OFFSET:.0e}.")
    print()

    # ---- X5b member-length sensitivity ----------------------------------
    print("  " + "-" * 74)
    print("  X5b  MEMBER LENGTH SENSITIVITY -- prediction (A), made quantitative")
    print("  " + "-" * 74)
    print("  If the lock is a tension-only member, its angle MOVES with the")
    print("  member's length. Column 2 is the DIAGONAL's taut angle, comparable")
    print("  with a table produced by a different script (T2 22682) -- numbers")
    print("  entering from outside. Column 4 is the highest taut angle over ALL")
    print("  intra-unit spans, which is what would actually bind if every vertex")
    print("  pair carried a member of that length.")
    print(f"    {'k':>6s} {'diagonal (here)':>18s} {'T2 22682':>12s} "
          f"{'dev':>10s} {'highest intra':>15s} {'binder':>10s}")
    ktab = []
    for k, recorded in K_TABLE:
        rows = intra_cross(k * STRUT_LEN)
        diag = [r for r in rows if r[1] in DIAGONALS]
        got = diag[0][0] if diag else float("nan")
        top = rows[0][0] if rows else float("nan")
        dev = abs(got - recorded)
        binder = "diagonal" if abs(top - got) < 1e-9 else "OTHER SPAN"
        ktab.append((k, got, recorded, dev, top, binder))
        print(f"    {k:6.2f} {got:18.9f} {recorded:12.6f} {dev:10.2e} "
              f"{top:15.9f} {binder:>10s}")
    ktab_ok = all(np.isfinite(d) and d < TOL["ktable"] for _, _, _, d, _, _ in ktab)
    kspan = max(g for _, g, _, _, _, _ in ktab) - min(g for _, g, _, _, _, _ in ktab)
    print(f"  agreement with the recorded table within {TOL['ktable']:.0e}: "
          f"{ktab_ok}")
    print(f"  the diagonal's taut angle MOVES by {kspan:.6f} deg over k in")
    print("  [0.90, 1.20]. That span is the non-vacuity of the row above it: an")
    print("  enumerator returning a constant would satisfy 'agrees at k = 1.00'")
    print("  and nothing else. Both are gated.")
    print()
    print("  A CORRECTION TO THE RECORDED k-TABLE'S SCOPE, found here. The table")
    print("  gives the DIAGONAL's angle, and for k <= 1.00 the diagonal is also")
    print("  the binding span. For larger k it is NOT: another intra-unit span")
    print("  goes taut at a higher angle, so a build using over-long members")
    print("  everywhere would stop before the diagonals ever came into play. The")
    print("  'binder' column above says which is which. The recorded table is")
    print("  correct for what it measures; it is its scope that needed narrowing.")
    print()
    print("  AND A SCOPE LIMIT ON THE CORRECTION ITSELF, which is the same kind")
    print("  of narrowing applied to this file: the 'binder' column ranges over")
    print("  INTRA-UNIT spans only. X5c finds inter-unit spans that go taut")
    print("  ABOVE the icosahedral phase on short, installable chords, so in an")
    print("  ASSEMBLED array the binder at any k may be an inter-unit span that")
    print("  this column never considered. The falsifier below therefore")
    print("  predicts the DIAGONAL's swing, which is the right prediction only")
    print("  for a build whose binding member is a diagonal -- exactly the build")
    print("  fact this file cannot settle.")
    print()

    # ---- X5c the assembled array ----------------------------------------
    print("  " + "-" * 74)
    print("  X5c  THE ASSEMBLED ARRAY: spans between IDENTIFIED wired points")
    print("  " + "-" * 74)
    print("  Contacts are unioned first, so each wired point is counted once and")
    print("  a span is measured between point CLASSES. Without that, a span whose")
    print("  far endpoint is identified into the near unit is reported as a new")
    print("  inter-unit blocker when it is that unit's own diagonal seen from the")
    print("  other side of a shared vertex.")
    print()
    print("  THE INSTALLABILITY SCALE, derived and not chosen: two units in")
    print("  contact have their centres 2|v| apart, and every one of the twelve")
    print("  shared vertices has the same norm, so that is one number. At the")
    print(f"  icosahedral phase it is {CONTACT_SPACING_ICO:.9f} = "
          f"{CONTACT_SPACING_ICO / STRUT_LEN:.4f} x strut. A span SHORTER than")
    print("  that at a_ico runs between points closer together than two")
    print("  contacting unit centres, so a member across it is a plausible build")
    print("  member; a span several strut lengths long is not. This is the axis")
    print("  the first version of this section failed to rank by.")
    print()
    want = ("SC7 star (six-around-one)", "CUBE8-M (60/60/90 basis)",
            "CUBE8-R (60/60/60 basis)", "CUBE27-M")
    inter = {}
    for topo in topos:
        if topo.name not in want:
            continue
        members, units = _classes(topo)
        rows, st = _span_crossings(topo, STRUT_LEN)
        p_ico = _class_positions(A_ICO, topo, members)
        p_one = _class_positions(1.0, topo, members)
        dec = [r for r in rows if r[4] == "dec"]
        ints = [r for r in dec if r[3] == "inter"]
        intras = [r for r in dec if r[3] == "intra"]
        at_ico = [r for r in ints if abs(r[0] - ico_from_spans) < TOL["aico"]]
        below = [r for r in ints if r[0] < ico_from_spans - TOL["aico"]]
        intra_angles = sorted({round(r[0], 9) for r in intras})
        print(f"  {topo.name}: units {topo.n}, contacts {len(topo.contacts)}, "
              f"wired points {st['nc']}")
        print(f"    spans {st['npairs']}; CONSTANT in a: {st['n_const']} "
              f"({st['n_const_at_member']} of them sitting exactly AT strut")
        print(f"    length, so they can never lengthen and never go taut); "
              f"straddling: {st['nstrad']}")
        print(f"    crossings that BLOCK EXPANSION: {len(dec)} "
              f"({len(intras)} intra / {len(ints)} inter)")
        print(f"    spans whose member gives an INTERVAL rather than a "
              f"half-line: {st['n_interval']}")
        print(f"    intra-unit: expected {6 * topo.n} (6 diagonals x "
              f"{topo.n} units), found {len(intras)}, all at {intra_angles}")
        # ---- THE RANKING, CORRECTED. The first version of this section sorted
        # by DESCENDING a* and printed the top six, which are systematically the
        # LONGEST spans -- the ones a builder cannot wire -- and then argued
        # from their unbuildability that assembly does not move the lock. The
        # short competitors were invisible to that presentation, and two of them
        # are installable and bind ABOVE the icosahedral phase. The table now
        # ranks by SPAN AT a_ico, ascending, which is the axis that decides
        # whether a member can be installed across the span at all. The taut
        # angle is still printed; it is no longer the sort key.
        above = [r for r in ints if r[0] > ico_from_spans + TOL["aico"]]
        by_len = sorted(
            ((float(np.linalg.norm(p_ico[r[1]] - p_ico[r[2]])), r)
             for r in above), key=lambda t: (round(t[0], 9), -round(t[1][0], 9)))
        print(f"    inter-unit crossings ABOVE the icosahedral phase: "
              f"{len(above)}, ranked by SPAN AT a_ico (ascending):")
        print(f"    {'rank':>5s} {'span at a_ico':>15s} {'/strut':>8s} "
              f"{'a* (deg)':>16s} {'span at a=1':>12s}  units")
        shown = 0
        seen_key = set()
        for slen, (astar, pi, qi, kind, _) in by_len:
            key = (round(slen, 9), round(astar, 9))
            if key in seen_key:
                continue
            seen_key.add(key)
            shown += 1
            if shown > 6:
                break
            print(f"    {shown:5d} {slen:15.9f} {slen / STRUT_LEN:8.4f} "
                  f"{astar:16.9f} "
                  f"{float(np.linalg.norm(p_one[pi] - p_one[qi])):12.6f}  "
                  f"{sorted(units[pi])}/{sorted(units[qi])}")
        bind_len = by_len[0][0] if by_len else float("nan")
        bind_a = by_len[0][1][0] if by_len else float("nan")
        n_short = sum(1 for slen, _ in by_len if slen < CONTACT_SPACING_ICO)
        print(f"    SHORTEST such span at a_ico: {bind_len:.9f} = "
              f"{bind_len / STRUT_LEN:.4f} x strut, going taut at "
              f"a* = {bind_a:.9f}")
        print(f"    how many of the {len(above)} are SHORTER at a_ico than the "
              f"contact spacing {CONTACT_SPACING_ICO:.6f}: {n_short}")
        print(f"    (a span shorter at a_ico than the distance between two")
        print(f"     contacting unit centres is a NEAR-NEIGHBOUR chord: a wire")
        print(f"     across it is installable, and it binds at a* > a_ico.)")
        lowest = min((r[0] for r in ints), default=float("nan"))
        print(f"    HIGHEST taut angle {dec[0][0] if dec else float('nan'):.9f}"
              f";  LOWEST inter-unit {lowest:.9f}")
        print(f"    inter-unit spans taut AT the icosahedral phase: "
              f"{len(at_ico)}    strictly BELOW it: {len(below)}")
        for astar, pi, qi, _, _ in at_ico + below:
            print(f"      a* = {astar:.9f}  units {sorted(units[pi])} and "
                  f"{sorted(units[qi])}  span at a_ico "
                  f"{float(np.linalg.norm(p_ico[pi] - p_ico[qi])):.9f}")
        print(f"    coarse-step risk: largest |d span / d a| "
              f"{st['max_slope']:.4f} per deg, x {SPAN_COARSE_STEP} deg =")
        print(f"    {st['coarse_risk']:.4f} of a strut -- PRICED here, and")
        print(f"    BOUNDED by X5d's two-incommensurate-coarse-step row, which")
        print(f"    is what a price with no bound attached was missing")
        print(f"    refinements that failed to converge (gated at zero): "
              f"{st['n_brentq_fail']}")
        inter[topo.name] = dict(top=(dec[0][0] if dec else float("nan")),
                                n_inter=len(ints), n_intra=len(intras),
                                n_at_ico=len(at_ico), n_below=len(below),
                                below=sorted({round(r[0], 9) for r in below}),
                                lowest=lowest, intra_angles=intra_angles,
                                n_interval=st["n_interval"], st=st,
                                n_above=len(above), bind_len=bind_len,
                                bind_a=bind_a, n_short=n_short,
                                n_fail=st["n_brentq_fail"])
        print()

    star = inter.get("SC7 star (six-around-one)", {})
    lat = {k: v for k, v in inter.items() if k != "SC7 star (six-around-one)"}
    print("  WHAT X5c FOUND. Three things, none of them visible from one unit.")
    print()
    print("  (1) THE INTRA-UNIT PICTURE SURVIVES ASSEMBLY UNCHANGED. Every")
    print("      cluster shows exactly 6 blocking spans per unit, all at the")
    print("      icosahedral phase and nowhere else. Assembling units neither")
    print("      adds nor removes a unit's own diagonals.")
    print()
    print(f"  (2) THE STAR HAS {star.get('n_at_ico', 0)} GENUINE INTER-UNIT SPANS")
    print("      TAUT AT EXACTLY THE ICOSAHEDRAL PHASE -- endpoints in different")
    print("      units with no unit in common. The arithmetic is exact: the")
    print("      neighbour at 2 v_i carries vertex j at v_j + 2 v_i, the")
    print("      neighbour at 2 v_j carries vertex i at v_i + 2 v_j, and the")
    print("      difference is v_j - v_i. So that chord IS the intra-unit span")
    print("      (i, j), and when (i, j) is a folding diagonal it reaches strut")
    print("      length at exactly the icosahedral phase. It exists only when a")
    print("      GENERATOR PAIR is itself a diagonal pair -- the M-type basis has")
    print("      one such pair, the R-type basis none.")
    print("      AND IT DISAPPEARS WHEN THE MEDIATING UNIT AT 2 v_i + 2 v_j IS")
    print("      PRESENT: both endpoints are then identified into it and the")
    print("      chord is that unit's own diagonal. The full boxes show zero:")
    for name, v in sorted(lat.items()):
        print(f"        {name:28s} inter-unit spans at a_ico: "
              f"{v.get('n_at_ico', -1)}")
    print()
    print("      RETRACTION. The first version of this file read those zeros as")
    print("      'the extra blocker belongs to an INCOMPLETE cluster -- a")
    print("      physical array's boundary -- and a finite wooden array has")
    print("      boundaries.' THAT IS WRONG, in two separate ways.")
    print("      FIRST, the zeros are FORCED. In any full rectangular box, if")
    print("      c+e_i and c+e_j are present then so is c+e_i+e_j, so the")
    print("      mediating site is never missing and the row cannot report")
    print("      anything but zero. Exhaustively, over every site and every")
    print("      generator pair:")
    for box in ((2, 2, 2), (3, 3, 3)):
        npair, miss = mediating_site_audit(box)
        print(f"        box {box}: co-neighbour pairs {npair:3d}, "
              f"mediating site MISSING {miss}")
    print("      SECOND, boundaries are not the discriminator. CUBE27-M is a")
    print("      finite cluster with boundaries on every face and reports ZERO.")
    print("      What decides it is a VACANCY at 2 v_i + 2 v_j. Four clusters")
    print("      built to show that, each of them finite and bounded everywhere:")
    vac = {}
    for topo in vacancy_topologies():
        mem_v, un_v = _classes(topo)
        rows_v, st_v = _span_crossings(topo, STRUT_LEN)
        dec_v = [r for r in rows_v if r[4] == "dec" and r[3] == "inter"]
        at_v = [r for r in dec_v if abs(r[0] - ico_from_spans) < TOL["aico"]]
        vac[topo.name] = dict(n=topo.n, n_at_ico=len(at_v),
                              n_fail=st_v["n_brentq_fail"])
        print(f"        {topo.name:34s} units {topo.n:2d}   "
              f"inter-unit spans at a_ico: {len(at_v)}")
    print("      L-TROMINO-M has a boundary everywhere and ONE extra member;")
    print("      QUAD-M is the SAME cluster with the mediating unit added back")
    print("      and has NONE; HOLED8-M is CUBE8-M with that one site removed")
    print("      and has ONE where CUBE8-M has zero. The R-basis control has")
    print("      none, because no R generator pair is a diagonal pair.")
    print("      SO: THE DISCRIMINATOR IS A VACANCY, NOT A BOUNDARY. A solid")
    print("      rectangular build carries no extra member at all; a sparse or")
    print("      HOLED build does. That is the opposite of what the retracted")
    print("      sentence would have told a builder.")
    print()
    print("  (3) A NEW ANGLE, BELOW THE ICOSAHEDRAL PHASE, THAT NO SINGLE UNIT")
    print("      HAS. Inter-unit spans go taut at angles strictly below 22.24:")
    for name, v in sorted(inter.items()):
        print(f"        {name:28s} count {v.get('n_below', -1):3d}   "
              f"angles {v.get('below', [])}")
    print("      Basis-dependent again: present in every M-basis cluster,")
    print("      absent from the R-basis one. It does NOT bind, because the")
    print("      array's allowed set is the intersection of half-lines a >= a*")
    print("      and the diagonals' 22.238756093 is the larger bound. It would")
    print("      bind in a build with no diagonal members.")
    print()
    print("  AND A CAUTION THE SINGLE-UNIT PICTURE DOES NOT PREPARE YOU FOR:")
    print("  inter-unit spans are NOT MONOTONE in a. The sub-icosahedral one")
    print("  falls to a minimum and rises again, so a member across it goes taut")
    print("  at TWO angles and confines a to an INTERVAL, not a half-line -- it")
    print("  blocks contraction as well as expansion. The count of such spans is")
    print("  printed per cluster above. The intra-unit diagonals are monotone on")
    print("  [0, 60] (1.414 L at the vector equilibrium down to 0 at the")
    print("  octahedron), which is exactly why THEY give the one-sided behaviour")
    print("  the owner reports and these do not.")
    print()
    print("  (4) AND THE RETRACTION THAT MATTERS MOST. The first version of this")
    print("      file concluded: 'THE LOCK ANGLE IS UNCHANGED BY ASSEMBLY. Every")
    print("      span taut at the icosahedral phase has the same length function,")
    print("      so the array MULTIPLIES the members holding the lock without")
    print("      MOVING it.' THAT IS FALSE, and this file's own enumeration is")
    print("      what refutes it. The defence offered was that a strut-length")
    print("      member across a span SEVERAL STRUT LENGTHS long is not")
    print("      installable -- true of the rows that version printed, because it")
    print("      sorted by DESCENDING a* and printed six rows, which selects")
    print("      exactly the long spans. Ranked instead by span at a_ico:")
    print(f"      {'cluster':30s} {'shortest span@a_ico':>20s} {'/strut':>8s} "
          f"{'binds at a*':>13s}")
    for name, v in sorted(inter.items()):
        print(f"      {name:30s} {v.get('bind_len', float('nan')):20.9f} "
              f"{v.get('bind_len', float('nan')) / STRUT_LEN:8.4f} "
              f"{v.get('bind_a', float('nan')):13.6f}")
    print("      Those spans are SHORTER at a_ico than the distance between two")
    print("      contacting unit centres. They are installable, and a")
    print("      strut-length member across one imposes a >= a* with a* ABOVE")
    print("      the icosahedral phase, which is the LARGER lower bound, so it")
    print("      BINDS and the diagonals never engage.")
    print("      CORRECTED STATEMENT: the array does NOT merely multiply members")
    print("      at the same angle -- it introduces EARLIER-BINDING ones. Which")
    print("      angle an assembled array actually locks at depends on WHICH of")
    print("      these chords the build wires, and that is a build fact this")
    print("      file cannot settle. What the file can say is that 22.238756093")
    print("      is the lock angle of a build carrying the six intra-unit")
    print("      diagonals AND NOT these shorter inter-unit chords.")
    print()
    print("  " + "-" * 74)
    print("  X5d  GRID INDEPENDENCE -- both grids are choices, so both are gated")
    print("  " + "-" * 74)
    print("  TWO grids select data here and they need different treatment.")
    print("  THE FINE GRID only selects the BRACKET that brentq then refines, so")
    print("  it cannot move an angle, only lose or double-count a root. THE")
    print("  COARSE GRID decides which spans are LOOKED AT AT ALL, so it can")
    print("  delete a finding outright -- and it did: see `_span_crossings`.")
    print("  Both arms below use ABSOLUTE literals incommensurate with the")
    print("  defaults. An earlier version wrote the fine comparison step as")
    print("  `SWEEP_STEP * 0.37`, which made the row a RATIO test: coarsening")
    print("  SWEEP_STEP coarsened both arms together and the row printed the")
    print("  identical string at 0.0259, 0.037, 0.07, 0.10, 0.74 and 2.0. A")
    print("  step-independence row whose two arms move in lock-step is")
    print("  structurally blind to any mutation of the step.")
    star_topo = [t for t in topos if t.name.startswith("SC7")][0]
    step_dev, step_counts = float("inf"), (0, 0)
    ra, _ = _span_crossings(star_topo, STRUT_LEN, fine=SWEEP_STEP)
    rb, _ = _span_crossings(star_topo, STRUT_LEN, fine=SWEEP_STEP_ALT)
    aa = sorted(round(r[0], 9) for r in ra if r[4] == "dec")
    bb = sorted(round(r[0], 9) for r in rb if r[4] == "dec")
    step_counts = (len(aa), len(bb))
    if len(aa) == len(bb) and aa:
        step_dev = max(abs(x - y) for x, y in zip(aa, bb))
    print()
    print(f"    FINE   step {SWEEP_STEP} : {len(aa)} blocking crossings;  "
          f"step {SWEEP_STEP_ALT} (absolute): {len(bb)}")
    print(f"    largest disagreement between the two angle lists: "
          f"{step_dev:.3e} deg")
    print("    The two lists first disagreed by EIGHT crossings, all of them a")
    print("    genuine inter-unit span reaching strut length at EXACTLY a = 30,")
    print("    a grid point of the default step. The obvious sign test (product")
    print("    of the two endpoint values < 0) skips a root that lands on a grid")
    print("    point; the enumeration now uses a half-open convention instead.")
    rc, _ = _span_crossings(star_topo, STRUT_LEN, coarse=SPAN_COARSE_ALT)
    cc = sorted(round(r[0], 9) for r in rc if r[4] == "dec")
    coarse_counts = (len(aa), len(cc))
    coarse_dev = float("inf")
    if len(aa) == len(cc) and aa:
        coarse_dev = max(abs(x - y) for x, y in zip(aa, cc))
    print()
    print(f"    COARSE step {SPAN_COARSE_STEP} : {len(aa)} blocking crossings;  "
          f"step {SPAN_COARSE_ALT} (absolute, incommensurate): {len(cc)}")
    print(f"    largest disagreement between the two angle lists: "
          f"{coarse_dev:.3e} deg")
    print("    THIS ARM IS THE NEW ONE, and it is a regression test for a real")
    print("    loss. With the old zero-span pre-filter these two counts were 194")
    print("    and 190: the coarse grid at 0.13 LANDS on the angle where four")
    print("    spans pass through length zero, the filter threw them away, and")
    print("    the sub-icosahedral taut angle vanished from every M cluster. A")
    print("    FINER grid did not fix it (0.05 lost them too) -- only an")
    print("    INCOMMENSURATE one exposes it, which is the same lesson the fine")
    print("    arm above teaches, applied to the grid that can lose data rather")
    print("    than to the grid that cannot. The coarse-step risk printed per")
    print("    cluster above was a PRICE with no bound attached; this is the")
    print("    bound.")
    print()
    print("  ON THE LONG SPANS. A member of strut length across a span several")
    print("  strut lengths long at the icosahedral phase is not installable, and")
    print("  the 'span at a_ico' column is what says so. That observation is")
    print("  still true and it is still worth printing -- what it does NOT")
    print("  support, and was previously made to support, is the conclusion that")
    print("  assembly leaves the lock angle alone. See (4).")
    return dict(intra_rows=intra_rows, n_const=n_const, ico=ico_from_spans,
                ktab=ktab, ktab_ok=ktab_ok, kspan=kspan, inter=inter,
                distinct=distinct, diag_same=same,
                n_square_diag=len(SQUARE_DIAGONALS), n_diag=len(DIAGONALS),
                step_dev=step_dev, step_counts=step_counts,
                coarse_dev=coarse_dev, coarse_counts=coarse_counts,
                vac=vac, aico_dev=aico_dev, aico_ctrl=aico_ctrl,
                const_hi=const_hi, vary_lo=vary_lo,
                n_intra_fail=n_intra_fail[0],
                contact_spacing=CONTACT_SPACING_ICO)


# ==========================================================================
# X6  Q5 -- CHIRALITY FRUSTRATION (pre-registered)
# ==========================================================================

def x6_chirality(topos):
    print()
    print("=" * 78)
    print("X6  Q5: CHIRALITY FRUSTRATION -- pre-registered, then measured")
    print("=" * 78)
    print("  THE PRE-REGISTRATION IS EXTERNAL, and that is the point of it:")
    print("  Q5 is written into the bead, `bd show inviscid-qvf.11`, timestamped")
    print("  outside this file and before it existed. What follows restates the")
    print("  question and the prediction for the reader; a file asserting its own")
    print("  chronology ('written before any number below was produced', as an")
    print("  earlier version put it here) is not evidence of anything, and that")
    print("  claim is withdrawn in favour of the citation.")
    print("  Chirality inside a unit is FORCED, not")
    print("  chosen: sigma = sx*sy*sz, and uniform sigma = +1 breaks vertex")
    print("  sharing outright (24 distinct corners instead of 12). Contracting")
    print("  units twist. The frustration question is whether neighbours meeting")
    print("  at a shared vertex need compatible twist, and whether a consistent")
    print("  assignment exists on the lattice.")
    print()
    print("  THE PREDICTION MADE IN ADVANCE: not frustrated, for two independent")
    print("  reasons -- (i) a single-vertex contact is a BALL JOINT and imposes")
    print("  nothing at all on relative rotation, so there is no compatibility")
    print("  condition to satisfy; (ii) the mirror image of a unit is expected to")
    print("  be a PROPER rotation of that same unit, because the shared-vertex")
    print("  set is centrally symmetric, in which case handedness is not even an")
    print("  independent label. Both are measured below and either could fail.")
    print()
    print("  (i) DOES A CONTACT CONSTRAIN RELATIVE ROTATION? The doweled contact")
    print("      Jacobian's columns for unit i's rotation block are -[R v_k]_x,")
    print("      which has RANK 2, not 3 -- the component about the contact")
    print("      direction is missing. Measured over the two units of N2:")
    for a in (5.0, A_ICO, 45.0):
        v = verts(a)
        blk = -_hat(v[0])
        r, s = rank_of(blk)
        print(f"      a = {a:10.6f}   rank of the rotation block = {r}   "
              f"singular values {np.round(s, 6)}")
    print("      Rank 2 of 3 at every angle: one whole rotational direction per")
    print("      unit per contact is unconstrained. A ball joint transmits no")
    print("      twist, so no twist-compatibility condition exists to frustrate.")
    print("      AND THIS IS A DERIVATION, NOT A MEASUREMENT, so it carries no")
    print("      gate row. `rank_of(-_hat(u))` is 2 for ANY nonzero u -- a 3x3")
    print("      skew-symmetric matrix has rank 2 by construction and can never")
    print("      return 3. The earlier version gated it, which put a row in the")
    print("      table that no state of this array could redden. The printed")
    print("      values above are the evidence that the code computes what the")
    print("      derivation says; the derivation is what carries the claim.")
    print()
    print("  (ii) IS THE MIRROR IMAGE A PROPER ROTATION OF THE UNIT?")
    print("       DERIVED FIRST: the shared-vertex set V is centrally symmetric,")
    print("       so (-I) V = V as sets. For the mirror M = diag(-1,1,1),")
    print("       M V = M (-I) V = diag(1,-1,-1) V, and diag(1,-1,-1) has")
    print("       determinant +1. So the mirror is a ROTATION of the original.")
    print("       MEASURED, because a derivation is not a measurement:")
    mm = np.diag([-1.0, 1.0, 1.0])
    dd = np.diag([1.0, -1.0, -1.0])
    mdevs = []
    for a in (5.0, A_ICO, 45.0, 55.0):
        v = verts(a)
        a_set = (mm @ v.T).T
        b_set = (dd @ v.T).T
        d = np.linalg.norm(a_set[:, None, :] - b_set[None, :, :], axis=-1)
        haus = max(float(d.min(axis=1).max()), float(d.min(axis=0).max()))
        mdevs.append(haus)
        print(f"       a = {a:10.6f}   Hausdorff( M V , diag(1,-1,-1) V ) = "
              f"{haus:.3e}   det diag(1,-1,-1) = {np.linalg.det(dd):+.1f}")
    print("       CONTROL that can fail: the same comparison against a")
    print("       deliberately WRONG rotation must NOT come back zero.")
    ctrl = []
    wrong = rot(np.array([0.0, 0.0, 1.0]), 17.0)
    for a in (5.0, A_ICO, 45.0, 55.0):
        v = verts(a)
        a_set = (mm @ v.T).T
        b_set = (wrong @ v.T).T
        d = np.linalg.norm(a_set[:, None, :] - b_set[None, :, :], axis=-1)
        ctrl.append(max(float(d.min(axis=1).max()), float(d.min(axis=0).max())))
    print(f"       control Hausdorff values: "
          f"{', '.join(f'{c:.4f}' for c in ctrl)}   (must all be large)")
    print("       AND A LIMIT ON WHAT THE FIRST ROW PINS. The vertex set admits")
    print("       more than one rotation carrying M V onto itself, so replacing")
    print("       diag(1,-1,-1) with diag(-1,1,-1) also passes at 3.3e-16. The")
    print("       row therefore establishes that the mirror IS a rotation, not")
    print("       that it is THAT rotation; the control above is what makes it")
    print("       falsifiable at all, and the determinant is what makes it a")
    print("       PROPER one. Stated because a row that passes under a wrong")
    print("       matrix should say so rather than be read as pinning the matrix.")
    print()
    print("  (iii) MIXED-CHIRALITY ARRAYS, and this leg is ENTAILED by (ii)")
    print("        rather than independent of it. If the mirror unit is a proper")
    print("        rotation of the unit, then a 'mirrored' unit IS a rotated")
    print("        unit, and every sign pattern is the all-same placement problem")
    print("        re-seeded -- so it must solve. It is run anyway, as a")
    print("        consistency check on the implementation of that entailment,")
    print("        and it is labelled as such in the gate. Presenting it as a")
    print("        third independent leg, as an earlier version did, overstates")
    print("        the evidence by one leg.")
    sq = [t for t in topos if t.name.startswith("SQUARE4")][0]
    worst_pat, worst_res = None, -1.0
    n_ok = 0
    for bits in range(2 ** sq.n):
        signs = [1 if (bits >> i) & 1 == 0 else -1 for i in range(sq.n)]
        res = _solve_mixed(A_ICO, sq, signs)
        if res < SOLVE_TOL:
            n_ok += 1
        if res > worst_res:
            worst_res, worst_pat = res, tuple(signs)
    print(f"        patterns solved: {n_ok} of {2 ** sq.n}   worst residual "
          f"{worst_res:.3e} at pattern {worst_pat}")
    print()
    print("  VERDICT ON Q5, stated plainly as the brief demands: chirality")
    print("  frustration is NOT the mechanism, and the reason is stronger than")
    print("  'a consistent assignment exists'. There is no assignment to make.")
    print("  A single-vertex contact transmits no twist, and the mirror unit is")
    print("  the same unit rotated, so 'handedness' is not a degree of freedom of")
    print("  this array at all.")
    return dict(mirror_dev=max(mdevs), ctrl_min=min(ctrl),
                mixed_ok=(n_ok == 2 ** sq.n), mixed_worst=worst_res,
                rot_block_rank=rank_of(-_hat(verts(A_ICO)[0]))[0])


def _solve_mixed(a, topo, signs):
    """Solve the placement problem with per-unit handedness `signs`.

    A unit with sign -1 uses the MIRRORED vertex set. No raise on any path.
    """
    n = topo.n
    v = verts(a)
    mm = np.diag([-1.0, 1.0, 1.0])
    vs = [v if s > 0 else (mm @ v.T).T for s in signs]
    t = topo.sites(v).copy()
    w = np.zeros((n, 3))
    rng = np.random.default_rng(7)
    w = rng.standard_normal((n, 3)) * 0.2
    t = t + rng.standard_normal((n, 3)) * 0.1
    w[0] = 0.0
    t[0] = 0.0
    cs = topo.contacts
    res = None
    for _ in range(400):
        q = [_rodrigues(x) for x in w]
        pos = [q[i] @ vs[i].T + t[i][:, None] for i in range(n)]
        res = np.concatenate([pos[i][:, k] - pos[j][:, l] for (i, k, j, l) in cs])
        jm = np.zeros((3 * len(cs), 6 * n))
        for e, (i, k, j, l) in enumerate(cs):
            r = 3 * e
            jm[r:r + 3, 6 * i:6 * i + 3] += -_hat(q[i] @ vs[i][k])
            jm[r:r + 3, 6 * i + 3:6 * i + 6] += np.eye(3)
            jm[r:r + 3, 6 * j:6 * j + 3] -= -_hat(q[j] @ vs[j][l])
            jm[r:r + 3, 6 * j + 3:6 * j + 6] -= np.eye(3)
        jm[:, 0:6] = 0.0
        with np.errstate(all="ignore"):
            dz = np.linalg.solve(jm.T @ jm + 1e-9 * np.eye(6 * n), -jm.T @ res)
        if not np.all(np.isfinite(dz)):
            break
        dz = dz.reshape(n, 6)
        step, improved = 1.0, False
        for _ in range(40):
            w2, t2 = w + step * dz[:, 0:3], t + step * dz[:, 3:6]
            q2 = [_rodrigues(x) for x in w2]
            pos2 = [q2[i] @ vs[i].T + t2[i][:, None] for i in range(n)]
            r2 = np.concatenate([pos2[i][:, k] - pos2[j][:, l]
                                 for (i, k, j, l) in cs])
            if np.linalg.norm(r2) < np.linalg.norm(res):
                improved = True
                break
            step *= 0.5
        if not improved:
            break
        w, t, res = w2, t2, r2
        if np.linalg.norm(res) < SOLVE_TOL * 1e-2:
            break
    return float(np.linalg.norm(res))


# ==========================================================================
# X7  Q6 -- INTER-UNIT CLEARANCE (measurement only)
# ==========================================================================

def _seg_seg(p1, q1, p2, q2):
    """Shortest distance between two segments. Lifted unchanged in method from
    jb_g_strut_clearance, whose routine is EXACT (4.4e-16 against an independent
    reference), not approximate."""
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e = d1 @ d1, d2 @ d2
    if a < 1e-14 or e < 1e-14:
        return float(np.linalg.norm(r))
    b, c, f = d1 @ d2, d1 @ r, d2 @ r
    den = a * e - b * b
    s = np.clip((b * f - c * e) / den, 0.0, 1.0) if den > 1e-12 else 0.0
    t = np.clip((b * s + f) / e, 0.0, 1.0)
    s = np.clip((b * t - c) / a, 0.0, 1.0)
    return float(np.linalg.norm((p1 + d1 * s) - (p2 + d2 * t)))


def _pt_tri(p, tri):
    """Distance from a point to a triangle (projection, clamped to the face)."""
    a, b, c = tri
    n = np.cross(b - a, c - a)
    nn = n @ n
    if nn < 1e-18:
        return min(_seg_seg(p, p, a, b), _seg_seg(p, p, b, c), _seg_seg(p, p, c, a))
    proj = p - n * ((p - a) @ n) / nn
    u = np.cross(b - a, proj - a) @ n
    v = np.cross(c - b, proj - b) @ n
    w = np.cross(a - c, proj - c) @ n
    if u >= 0 and v >= 0 and w >= 0:
        return float(abs((p - a) @ n) / np.sqrt(nn))
    return min(_seg_seg(p, p, a, b), _seg_seg(p, p, b, c), _seg_seg(p, p, c, a))


def _seg_tri_hits(p, q, tri):
    """Moller-Trumbore: does the segment pq pierce the triangle interior?"""
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


def _tri_tri(t1, t2):
    for i in range(3):
        if _seg_tri_hits(t1[i], t1[(i + 1) % 3], t2):
            return 0.0
        if _seg_tri_hits(t2[i], t2[(i + 1) % 3], t1):
            return 0.0
    best = np.inf
    for i in range(3):
        for j in range(3):
            best = min(best, _seg_seg(t1[i], t1[(i + 1) % 3],
                                      t2[j], t2[(j + 1) % 3]))
    for i in range(3):
        best = min(best, _pt_tri(t1[i], t2), _pt_tri(t2[i], t1))
    return float(best)


def _seg_tri(p, q, tri):
    """Shortest distance from segment pq to triangle `tri`. Exact for disjoint
    inputs; returns 0 when the segment pierces the face."""
    if _seg_tri_hits(p, q, tri):
        return 0.0
    best = min(_seg_seg(p, q, tri[i], tri[(i + 1) % 3]) for i in range(3))
    return float(min(best, _pt_tri(p, tri), _pt_tri(q, tri)))


def x7_clearance():
    print()
    print("=" * 78)
    print("X7  Q6: INTER-UNIT CLEARANCE vs a -- MEASUREMENT ONLY")
    print("=" * 78)
    print("  Demoted by the directional deduction: a collision blocks APPROACH,")
    print("  not separation, and cannot forbid expansion while permitting")
    print("  contraction. Measured anyway because it is cheap and has never been")
    print("  computed. NO ADMISSIBILITY VERDICT IS ATTACHED. DECISION 16 permits")
    print("  interference in the model; the owner's physical plates cannot")
    print("  interpenetrate; measurement and verdict stay apart.")
    print()
    print("  Two units at the M-type contact 2*v_0. Plate-to-plate (filled")
    print("  triangles, with a Moller-Trumbore piercing test), strut-to-plate")
    print("  (segment against filled triangle), and strut-to-strut (no faces at")
    print("  all, the reading Fuller's own jitterbug supports). Triangles sharing")
    print("  the contact point are EXCLUDED: they meet there by construction and")
    print("  would report 0 at every angle.")
    print()
    print(f"    {'a':>10s} {'plate-plate':>13s} {'strut-plate':>13s} "
          f"{'strut-strut':>13s} {'state':>16s}")
    grid = [1.0, 5.0, 8.0, 10.0, 12.0, 15.0, 18.0, 20.0, A_ICO, 25.0, 30.0,
            35.0, 40.0, 45.0, 50.0, 55.0, 58.0]
    rows = []
    for a in grid:
        x = corners(a)
        v = verts(a)
        off = 2.0 * v[0]
        tri_a = [x[i] for i in range(8)]
        tri_b = [x[i] + off for i in range(8)]
        cpt = v[0]
        touch_a = [any(np.linalg.norm(p - cpt) < 1e-7 for p in t) for t in tri_a]
        touch_b = [any(np.linalg.norm(p - cpt) < 1e-7 for p in t) for t in tri_b]
        pp, sp, ss = np.inf, np.inf, np.inf
        for i in range(8):
            for j in range(8):
                if touch_a[i] and touch_b[j]:
                    continue
                pp = min(pp, _tri_tri(tri_a[i], tri_b[j]))
                for e in range(3):
                    sp = min(sp,
                             _seg_tri(tri_a[i][e], tri_a[i][(e + 1) % 3], tri_b[j]),
                             _seg_tri(tri_b[j][e], tri_b[j][(e + 1) % 3], tri_a[i]))
                    for f in range(3):
                        ss = min(ss, _seg_seg(tri_a[i][e], tri_a[i][(e + 1) % 3],
                                              tri_b[j][f], tri_b[j][(f + 1) % 3]))
        state = "PLATES OVERLAP" if pp < 1e-9 else ""
        rows.append((a, float(pp), float(sp), float(ss)))
        print(f"    {a:10.6f} {pp:13.9f} {sp:13.9f} {ss:13.9f} {state:>16s}")
    contact = [r[0] for r in rows if r[1] < 1e-9]
    fc = contact[0] if contact else float("nan")
    pmin = min(r[1] for r in rows)
    pmax = max(r[1] for r in rows)
    ico_gap = [r[1] for r in rows if r[0] == A_ICO][0]
    tail = [r for r in rows if r[0] >= 12.0]
    monotone = all(tail[i + 1][1] < tail[i][1] for i in range(len(tail) - 1))
    faces = max(abs(r[1] - r[3]) for r in rows)
    # THE RATIO CANNOT DIVIDE BY ZERO. `pmin` is exactly 0.0 the moment plates
    # touch -- and plate contact is precisely the condition the "no contact" row
    # exists to detect, so an unguarded ratio meant that the ONLY mutation
    # capable of reddening that row destroyed the verdict table instead. It did:
    # moving the two units to 1.4 * v[0] raised ZeroDivisionError at this line,
    # exit 1 with ZERO gate rows printed. That is the same class as the pairing
    # traceback the file already closed, in the row that exists to prevent
    # vacuity. An infinite ratio is the honest reading (the statistic spans from
    # zero) and it keeps the table intact so the CONTACT row can go red.
    ratio = (pmax / pmin) if pmin > 0.0 else float("inf")
    rstr = "inf (contact)" if not np.isfinite(ratio) else f"{ratio:.3f}"
    print()
    print(f"  inter-unit plate contact anywhere on the grid: "
          f"{'NONE' if not contact else f'first at {fc:.6f} deg'}")
    print(f"  clearance at the icosahedral phase: {ico_gap:.9f}   "
          f"range over the grid {pmin:.9f} .. {pmax:.9f}  "
          f"(ratio {rstr})")
    print(f"  strictly DECREASING on [12, 58]: {monotone}")
    print(f"  |plate-plate  -  strut-strut| over the grid: {faces:.3e}")
    print()
    print("  THE FACES NEVER WIN HERE, which is a ONE-SIDED statement and is")
    print("  stated as one. `_tri_tri` minimises over the same nine edge-edge")
    print("  terms that `ss` minimises, plus point-triangle terms, so pp <= ss")
    print("  holds by construction and the only content of the comparison is")
    print("  that no vertex-face approach ever undercuts the edge-edge one --")
    print("  the figures are bit-identical, not agreeing to a tolerance. An")
    print("  earlier wording said 'plate-plate == strut-strut', which reads as a")
    print("  two-sided agreement the arithmetic cannot supply. jb_g found the")
    print("  same for INTRA-unit clearance; this is the array analogue, and it")
    print("  had to be checked separately.")
    print()
    print("  EXPECTED RESULT, and it is the one obtained: inter-unit clearance is")
    print("  nowhere near zero at 22.24, has NO local feature there, and falls")
    print("  monotonically through it. Reported plainly, not strained to fit.")
    print("  Three gate rows: no contact anywhere on the grid; monotone through")
    print("  the icosahedral phase, so no feature; and the statistic actually")
    print("  MOVES over the sweep -- without that third row, 'no feature' would")
    print("  be satisfied by a measurement that returned a constant.")
    print()
    print("  For contrast, INTRA-unit interference (Java `Interpenetration`,")
    print("  edge-through-face) is 0 for a <= 60 and 24 throughout 60 < a < 120,")
    print("  so intra-unit interference starts at the octahedron and cannot cause")
    print("  a lock at 22.24 either.")
    return dict(rows=rows, first_contact=fc, ico_gap=ico_gap,
                no_contact=(not contact), monotone=monotone,
                ratio=ratio, faces_irrelevant=faces, pmin=pmin)




# ==========================================================================
# X8  VERDICT
# ==========================================================================

def x8_verdict(r0, r1, term, twoside, ranks, r5, r6, r7):
    print()
    print("=" * 78)
    print("X8  VERDICT")
    print("=" * 78)
    print("  Q1  IN-PHASE EXISTENCE. Topology-dependent, and the dependence is")
    print("      the result, not a nuisance:")
    for name, (good, bad) in sorted(term.items()):
        if not bad:
            print(f"        {name:34s} EXISTS AT EVERY SWEPT a")
        elif tuple(good) == (0.0,):
            print(f"        {name:34s} EXISTS ONLY AT a = 0")
        else:
            print(f"        {name:34s} exists on {good}")
    print("      SIX-AROUND-ONE (Fuller 784.30, and every cubic cluster built on")
    print("      it) closes at every angle, because a translate lattice's contact")
    print("      cycles close by the commutativity of translation. TWELVE-AROUND-")
    print("      ONE does not, and its obstruction is AT THE VECTOR EQUILIBRIUM,")
    print("      opening LINEARLY in a -- it is not a lock at any angle, it is an")
    print("      isolated point at a = 0.")
    print()
    print("  Q2  IS 22.238756093 DISTINGUISHED BY THE EQUALITY CONSTRAINTS? NO.")
    print("      Nothing in X2, X3 or X4 happens there. The twelve-around-one")
    print("      family terminates immediately off a = 0, nowhere near it; every")
    print("      six-around-one topology is smooth straight through it. The angle")
    print("      is printed everywhere it appears and is never rounded toward.")
    print()
    print("  Q3  IS THE OBSTRUCTION ONE-SIDED? NOT FROM EQUALITY CONSTRAINTS.")
    print("      Continuation solves on BOTH sides of the icosahedral phase at")
    print("      every step size tried, for every topology with a solution. A")
    print("      rank drop would in any case be two-sided. One-sidedness in this")
    print("      model can only come from an INEQUALITY, and X5 exhibits the")
    print("      inequality: span(a) <= member, on a span that is DECREASING in a,")
    print("      which is the half-line a >= a*. Expansion blocked, contraction")
    print("      free -- the observed asymmetry, with the right shape.")
    print()
    print("  Q4  RANK. Reported in X3 per topology, free and doweled, at six")
    print("      coarse angles and on a 0.25-degree sweep of a_ico +/- 2. On")
    print("      those samples the rank is CONSTANT and FULL ROW RANK wherever a")
    print("      solution exists -- both now gated, the second having previously")
    print("      been asserted by no row at all. Constant rank on a neighbourhood")
    print("      makes the solution set locally a smooth manifold THERE: no")
    print("      branch point and no boundary on the sampled neighbourhood.")
    print("      'No cusp anywhere on the swept interval' was claimed by an")
    print("      earlier version and is WITHDRAWN: nothing here tests for a cusp,")
    print("      and six angles plus a four-degree window are not an interval.")
    print()
    print("  Q5  CHIRALITY. Not the mechanism, and not because an assignment")
    print("      happens to exist: there is nothing to assign. A single-vertex")
    print("      contact is a ball joint whose rotation block has rank 2 of 3, so")
    print("      no twist is transmitted; and the mirror image of a unit is that")
    print("      unit rotated, so handedness is not a label the array carries.")
    print("      The rank-2 half is a DERIVATION about skew-symmetric 3x3")
    print("      matrices, the mirror half is a derivation MEASURED against a")
    print("      control that can fire, and the mixed-chirality sweep is entailed")
    print("      by the mirror result rather than independent of it.")
    print()
    print("  Q6  CLEARANCE. Measured, featureless at 22.24, no verdict attached.")
    print()
    print("  THE PRIMARY DELIVERABLE, X5. Every span in the assembled array that")
    print("  lengthens under expansion is enumerated. Of the 66 intra-unit vertex")
    print("  pairs, exactly SIX reach strut length while lengthening, and all six")
    print(f"  do so at {r5['ico']:.9f} degrees -- the folding square diagonals,")
    print("  derived twice and agreeing.")
    print()
    print("  THE ARRAY-LEVEL FINDINGS, which contradict this file's own first")
    print("  draft and are the reason the enumeration was worth running:")
    star = r5["inter"].get("SC7 star (six-around-one)", {})
    c8m = r5["inter"].get("CUBE8-M (60/60/90 basis)", {})
    c8r = r5["inter"].get("CUBE8-R (60/60/60 basis)", {})
    c27 = r5["inter"].get("CUBE27-M", {})
    print(f"    (i) the SIX-AROUND-ONE STAR has {star.get('n_at_ico', 0)} GENUINE")
    print("        inter-unit spans -- endpoints in different units with no unit")
    print("        in common -- taut at EXACTLY the icosahedral phase, running")
    print("        between two different neighbours of the centre. Their length")
    print("        function is v_j - v_i for two generator vertices: a folding")
    print("        diagonal when the generator pair is itself a diagonal pair.")
    print("        The M basis has one such pair; the R basis has none.")
    print(f"        The FULL boxes have {c8m.get('n_at_ico', -1)} (CUBE8-M), "
          f"{c8r.get('n_at_ico', -1)} (CUBE8-R), {c27.get('n_at_ico', -1)} "
          f"(CUBE27-M) -- but")
    print("        those zeros are FORCED: a full rectangular box always")
    print("        contains the mediating site (X5c audits it exhaustively), so")
    print("        the count cannot come back anything else. The DISCRIMINATOR")
    print("        IS A VACANCY AT 2 v_i + 2 v_j, NOT A BOUNDARY. X5c's holed")
    print("        clusters show it: L-TROMINO-M has 1, the same cluster with")
    print("        the mediating unit added back has 0, CUBE8-M minus that one")
    print("        site has 1 where CUBE8-M has 0. CUBE27-M is finite WITH")
    print("        boundaries on every face and has 0, which is what refutes the")
    print("        earlier reading. 'The extra blocker belongs to an INCOMPLETE")
    print("        cluster -- a physical wooden array has boundaries' is")
    print("        WITHDRAWN; it points a builder the wrong way. A solid")
    print("        rectangular build carries no extra member; a HOLED or sparse")
    print("        one does.")
    below = star.get("below", [])
    print(f"    (ii) a NEW taut angle BELOW the icosahedral phase: "
          f"{below[0] if below else float('nan'):.9f} degrees, the same number")
    print("        in every M-basis cluster and absent from the R basis. It does")
    print("        NOT bind -- the allowed set is an intersection of half-lines")
    print("        a >= a* and 22.238756093 is the larger bound -- but it is a")
    print("        real array-level number that no single unit has, and it would")
    print("        bind in a build carrying no diagonal members.")
    print("    (iii) inter-unit spans are NOT MONOTONE in a. A member across the")
    print("        sub-icosahedral one goes taut at TWO angles and confines a to")
    print("        an INTERVAL, blocking contraction as well as expansion. The")
    print("        intra-unit diagonals ARE monotone on [0, 60], which is exactly")
    print("        why they and not these give the one-sided behaviour observed.")
    print()
    print("    (iv) AND THE ARRAY DOES MOVE THE LOCK. An earlier version of this")
    print("        verdict said 'THE LOCK ANGLE IS UNCHANGED BY ASSEMBLY -- the")
    print("        array MULTIPLIES the members holding the lock without MOVING")
    print("        it'. That is RETRACTED. Ranked by span length at a_ico rather")
    print("        than by taut angle, every cluster has short inter-unit chords")
    print(f"        that go taut ABOVE the icosahedral phase: the shortest is")
    print(f"        {star.get('bind_len', float('nan')):.6f} = "
          f"{star.get('bind_len', float('nan')) / STRUT_LEN:.4f} x strut at "
          f"a_ico, binding at")
    print(f"        a* = {star.get('bind_a', float('nan')):.6f}, and it is")
    print("        shorter at a_ico than the distance between two contacting")
    print("        unit centres, so a member across it is installable. Because")
    print("        the allowed set is an intersection of half-lines a >= a*, the")
    print("        LARGER bound wins: such a member BINDS and the diagonals")
    print("        never engage. The array introduces EARLIER-BINDING members,")
    print("        not merely more members at the same angle. Which one governs")
    print("        a given build is a BUILD FACT.")
    print()
    print("  VERDICT ON THE BEAD'S FORK, in the order the evidence supports it:")
    print("  (B) IS REFUTED FOR THE IDEALISED EQUALITY-CONSTRAINED VARIETY over")
    print("  every topology built here -- constant, full row rank, two-sided")
    print("  continuation, nothing whatever happening at 22.238756093. HENCE (A)")
    print("  BY ELIMINATION. But THE MECHANISM WITHIN (A) IS UNRESOLVED, and")
    print("  that is not a caveat below the verdict, it is part of the verdict:")
    print("   * the bead's (A) is the WIRE, and this file shows the owner's")
    print("     inter-unit wire sits AT the contact, where the span is")
    print("     IDENTICALLY ZERO. A member of any length across a zero-length")
    print("     span constrains only coincidence. On this file's own")
    print("     measurement the inter-unit wires CANNOT be the blocker.")
    print("   * the only spans taut at exactly a_ico are the six INTRA-unit")
    print("     folding diagonals, and the described build -- solid triangles,")
    print("     vertex wires, centre dowels -- contains no member across them.")
    print("   * (A) also covers the DOWEL GUIDE, and the dowel branch is")
    print("     UNTESTED here. The doweled model in this file uses the TRUE PATH")
    print("     TANGENT; the owner's rig uses an ELLIPTICAL guide. Whether those")
    print("     are the same curve is bead P1 and is not addressed here at all,")
    print("     so the one thing that could still make this a genuine boundary")
    print("     of the REALISED mechanism is precisely what is deferred.")
    print("  What IS established without qualification: the equality-constrained")
    print("  variety is smooth and two-sided there, and an inequality reproduces")
    print("  the observation at the observed angle with no fitting -- IF the")
    print("  build carries a member across one of the enumerated spans. So the")
    print("  right correction to (A)'s stated consequence stands ('no consequence")
    print("  for the EQUALITY-constrained variety, fully determined for any model")
    print("  carrying tension-only members'), and the identification of WHICH")
    print("  member does not.")
    print()
    print("  WHAT WOULD FALSIFY THIS, and it is the owner's P0 experiment, which")
    print("  remains cheaper and more decisive than anything above: change the")
    print("  member lengths and see whether the lock angle MOVES. X5b predicts,")
    print(f"  with no free parameters, that it moves by {r5['kspan']:.6f} degrees")
    print("  across k in [0.90, 1.20], and gives the angle for each k.")
    print()
    print("  WHAT THIS FILE DOES NOT SETTLE. Whether the owner's build actually")
    print("  contains members spanning the six folding diagonals is a BUILD FACT")
    print("  and cannot be measured from geometry. If it does not, the mechanism")
    print("  is still open and the next candidate is the dowel guide curve (bead")
    print("  P1: is the true vertex path exactly an ellipse?), which is not")
    print("  addressed here at all.")
    print()
    print("  SCOPE LIMITS, stated so that nothing above is read wider than it was")
    print("  measured:")
    print("   * EVERYTHING HERE IS IN PHASE -- every unit at the same a. An")
    print("     OUT-OF-PHASE array is not modelled at all, and that is precisely")
    print("     the object the epic's propagation question needs. What this file")
    print("     supplies toward it is the assembled constraint machinery, not an")
    print("     answer.")
    print("   * FINITE CLUSTERS ONLY, up to 27 units. Statements about 'the")
    print("     periodic array' are inferences from the largest clusters built,")
    print("     not from a periodic calculation.")
    print("   * NON-EXISTENCE IS EVIDENCE, NOT PROOF. Where a topology is")
    print("     reported to have no in-phase solution, that is a damped")
    print("     Gauss-Newton failing from six starts, one of them the exact")
    print("     translate placement. It is a strong negative and it is not a")
    print("     theorem.")
    print("   * AND THE MIRROR OF THAT BULLET, which the first version omitted:")
    print("     EXISTENCE HERE IS BY CONSTRUCTION, NOT BY DISCOVERY. 'Six-")
    print("     around-one closes at EVERY swept a' is reported from a solver")
    print("     whose DEFAULT SEED is the pure-translate placement, and that")
    print("     placement satisfies every contact identically by central")
    print("     symmetry: vertex k of the unit at s sits at s + v_k, and vertex")
    print("     anti(k) of the unit at s + 2 v_k sits at s + 2 v_k - v_k =")
    print("     s + v_k. The 1e-16 residuals are the seed. The structural reason")
    print("     is the result; the solver number is a consistency check on it.")
    print("     The discipline has to be symmetric or it is not discipline.")
    print("   * UNIQUENESS WAS NEVER ASKED AND IS NOT ANSWERED. Nothing here")
    print("     tests whether the translate branch is the ONLY in-phase branch,")
    print("     or whether it is isolated. 'An in-phase configuration exists' is")
    print("     the trivial half of that question.")
    print("   * THE FOUR DECLARATIONS (kernel, mass model, primitive, metric")
    print("     form) are INAPPLICABLE here, not forgotten: every quantity in")
    print("     this file is kinematic. See the module docstring.")
    print("   * MEMBERS ARE IDEAL AND UNIFORM. Real wire has slack, stretch and")
    print("     knots, none of which is modelled. The k-sweep is the only")
    print("     concession to member length being uncertain.")
    print("   * NO ADMISSIBILITY VERDICT ANYWHERE. DECISION 16 permits")
    print("     interference in the model; X7 measures clearance and stops.")
    print("   * The span enumeration is over the wired points of the clusters")
    print("     built. Spans between units further apart than those clusters")
    print("     reach are not enumerated.")


# ==========================================================================
# THE GATE
# ==========================================================================

def gate(r0, r1, guard_fired, ranks, sig_free, sig_dow, inph, breathe,
         inph_bad, rot_ok, rot_ctrl, rot_bad, colmin,
         term, twoside, absent, r5, r6, r7):
    """Every check's verdict in one table, and this process's exit code.

    Every number here is COMPUTED from what was passed in. A verdict that is
    printed but not asserted is a verdict nobody can break.
    """
    sixaround = [n for n in term if n.startswith(("CUBE", "SC7", "SQUARE", "N2",
                                                 "CHAIN"))]
    six_all = all(not term[n][1] for n in sixaround)
    fcc = [n for n in term if n.startswith("FCC13")]
    fcc_only0 = all(tuple(term[n][0]) == (0.0,) for n in fcc) if fcc else False
    two_all = all(twoside.values()) if twoside else False
    live_ranks = {k: v for k, v in ranks.items() if not k.startswith("FCC13")}
    # `all()` over an empty dict is True. Every one of the three rank rows below
    # would then pass on a build that measured no topology at all, so the
    # non-emptiness is part of each predicate rather than assumed.
    rank_const = bool(live_ranks) and all(v["fine_free"] and v["fine_dow"]
                                          for v in live_ranks.values())
    # The COARSE sweep was computed, printed, and DISCARDED for every solvable
    # topology -- and it is the sweep that flips to False under a loosened
    # RANK_RTOL. Gated separately from the fine one, because they answer
    # different questions (a four-degree window versus the whole range).
    rank_const_coarse = bool(live_ranks) and all(v["coarse_const"]
                                                 for v in live_ranks.values())
    # FULL ROW RANK, the central Q4 claim, previously asserted by no row.
    full_row = bool(live_ranks) and all(v["full_row"]
                                        for v in live_ranks.values())
    n_full = sum(1 for v in live_ranks.values() if v["full_row"])
    n_intra = len(r5["intra_rows"])
    intra_one_angle = len(r5["distinct"]) == 1
    intra_at_ico = (abs(r5["ico"] - A_ICO) < TOL["aico"]
                    if np.isfinite(r5["ico"]) else False)
    inter_nonempty = all(v["n_inter"] > 0 for v in r5["inter"].values())
    sizes = {t.name: t.n for t in build_topologies()}
    per_unit_ok = all(v["n_intra"] == 6 * sizes[name]
                      and v["intra_angles"] == [round(r5["ico"], 9)]
                      for name, v in r5["inter"].items())
    star = r5["inter"].get("SC7 star (six-around-one)", {})
    cube8m = r5["inter"].get("CUBE8-M (60/60/90 basis)", {})
    cube8r = r5["inter"].get("CUBE8-R (60/60/60 basis)", {})
    cube27 = r5["inter"].get("CUBE27-M", {})
    star_at_ico = star.get("n_at_ico", 0) > 0
    lattice_clean = (cube8m.get("n_at_ico", -1) == 0
                     and cube8r.get("n_at_ico", -1) == 0
                     and cube27.get("n_at_ico", -1) == 0)
    # C1: the BINDING competitor. Over every cluster, the shortest inter-unit
    # span at a_ico among crossings ABOVE a_ico -- the one a builder could
    # actually install -- and the angle it binds at. Both gated, because the
    # retracted headline ("the lock angle is unchanged by assembly") is exactly
    # the claim these two numbers refute, and a retraction that nothing asserts
    # is as unfalsifiable as the claim it replaced.
    bl = [v.get("bind_len", float("nan")) for v in r5["inter"].values()]
    ba = [v.get("bind_a", float("nan")) for v in r5["inter"].values()]
    bind_installable = (bool(bl) and all(np.isfinite(x) for x in bl)
                        and max(bl) < r5["contact_spacing"])
    bind_above = (bool(ba) and all(np.isfinite(x) for x in ba)
                  and min(ba) > r5["ico"] + TOL["aico"])
    bind_worst = max(bl) / STRUT_LEN if bind_installable else float("nan")
    bind_lowest_a = min(ba) if ba and all(np.isfinite(x) for x in ba) else float("nan")
    n_short_all = min((v.get("n_short", -1) for v in r5["inter"].values()),
                      default=-1)
    # C3: the VACANCY clusters. The mediating site absent gives 1, present 0.
    vac = r5.get("vac", {})
    def _vn(pref):
        for k, v in vac.items():
            if k.startswith(pref):
                return v["n_at_ico"]
        return -1
    vac_absent = _vn("L-TROMINO-M")
    vac_present = _vn("QUAD-M")
    vac_holed = _vn("HOLED8-M")
    vac_rbasis = _vn("L-TROMINO-R")
    vacancy_ok = (vac_absent == 1 and vac_present == 0 and vac_holed == 1
                  and vac_rbasis == 0)
    # `miss == 0` over ZERO pairs would also be zero, so the pair count is part
    # of the predicate: the audit has to have looked at something.
    _aud = [mediating_site_audit(b) for b in ((2, 2, 2), (3, 3, 3))]
    box_forced = all(n > 0 and m == 0 for n, m in _aud)
    box_pairs = sum(n for n, _ in _aud)
    n_fail_total = (sum(v.get("n_fail", 0) for v in r5["inter"].values())
                    + sum(v.get("n_fail", 0) for v in vac.values())
                    + r5["n_intra_fail"])
    # the sub-icosahedral inter-unit angle: present in every M basis, absent
    # from the R basis, and the SAME number in all three M clusters
    mbelow = [star.get("below", []), cube8m.get("below", []),
              cube27.get("below", [])]
    m_below_ok = (all(len(b) == 1 for b in mbelow)
                  and len({b[0] for b in mbelow}) == 1
                  and mbelow[0][0] < r5["ico"])
    r_below_ok = cube8r.get("n_below", -1) == 0
    a_below = mbelow[0][0] if m_below_ok else float("nan")
    clearance_far = r7["no_contact"] or abs(r7["first_contact"] - A_ICO) > 5.0

    checks = [
        ("X0  hinge pairing: 12 hinges, 24 struts", r0["pair_ok"],
         str(r0["pair_ok"]), "True"),
        ("X0  CONTROL: N=1 free rank is 36 at every a", r0["ranks_ok"],
         str(r0["ranks_ok"]), "True"),
        ("X0  N=1 doweled model shape is 0 rows x 7 vars",
         r0["dow_rows"] == 0 and r0["dow_vars"] == 7,
         f"{r0['dow_rows']}x{r0['dow_vars']}", "0x7"),
        ("X0  analytic Jacobian vs exact-residual FD", r0["fd_ok"],
         f"{r0['fd']:.2e}", f"< {TOL['fd_jacobian']:.0e}"),
        ("X0  symmetric family satisfies its own hinges", r0["cres_ok"],
         f"{r0['cres']:.2e}", f"< {TOL['hinge_residual']:.0e}"),
        ("X0  6 global rigid motions lie in the kernel", r0["leak_ok"],
         f"{r0['leak']:.2e}", "< 1e-9"),
        ("X0  path tangent is a body-motion field, in ker", r0["tan_ok"],
         f"{r0['tan']:.2e}", "< 1e-6"),
        ("X0  sigma_36 at a=60 vs the Java-era record", r0["s36_ok"],
         f"{r0['s36']:.7f}", f"{SIGMA36_AT_60:.7f}"),
        ("X0  closed-form dV/da vs central differences", r0["dv_ok"],
         f"{r0['dvdev']:.2e}", "< 1e-6"),
        ("X1  antipode permutation constant in a", r1["perm_stable"],
         str(r1["perm_stable"]), "True"),
        ("X1  central symmetry of the 12 vertices", r1["perm_dev"] < TOL["antipode"],
         f"{r1['perm_dev']:.2e}", f"< {TOL['antipode']:.0e}"),
        ("X1  ... and ANTI is an involution (non-vacuity)", r1["anti_invol"],
         str(r1["anti_invol"]), "True"),
        ("X1  ... and FIXED-POINT-FREE (kills self-pairing)", r1["anti_fpf"],
         str(r1["anti_fpf"]), "True"),
        ("X1  NO orthogonal generator triple exists (784.30)",
         r1["n_ortho"] == 0 and r1["n_triples"] == 20,
         f"{r1['n_ortho']} of {r1['n_triples']}", "0 of 20"),
        ("X1  48 difference-closed pairs at a=0, exactly",
         r1["n_adjacent"] == 48, str(r1["n_adjacent"]), "48"),
        ("X1  ... and closure holds there", r1["defect0"] < 1e-12,
         f"{r1['defect0']:.2e}", "< 1e-12"),
        ("X1  ... and FAILS at the icosahedral phase", r1["defect_ico"] > 1.0,
         f"{r1['defect_ico']:.3f}", "> 1.0"),
        ("X2  six-around-one closes at EVERY swept a", six_all,
         str(six_all), "True"),
        ("X2  twelve-around-one closes ONLY at a = 0", fcc_only0,
         str(fcc_only0), "True"),
        ("X3  rank constant on the FINE sweep, free AND dow", rank_const,
         str(rank_const), "True"),
        ("X3  ... and on the COARSE sweep (was discarded)", rank_const_coarse,
         str(rank_const_coarse), "True"),
        ("X3  FULL ROW RANK everywhere solvable (was ungated)", full_row,
         f"{n_full}/{len(live_ranks)}", "all"),
        ("X3  conditioning margin over the noise floor, FREE",
         np.isfinite(sig_free) and sig_free > 1e6, f"{sig_free:.3e} x", "> 1e6 x"),
        ("X3  ... and DOWELED", np.isfinite(sig_dow) and sig_dow > 1e6,
         f"{sig_dow:.3e} x", "> 1e6 x"),
        ("X3b the IN-PHASE mode is in the doweled kernel",
         np.isfinite(inph) and inph < 1e-9, f"{inph:.2e}", "< 1e-9"),
        ("X3b ... and its BREATHING part is nonzero",
         np.isfinite(breathe) and breathe > 1e-3, f"{breathe:.3e}", "> 1e-3"),
        ("X3b ... and the SAME check fails where infeasible",
         np.isfinite(inph_bad) and inph_bad > 1e-3, f"{inph_bad:.2e}", "> 1e-3"),
        ("X3c global rigid ROTATION is in the doweled kernel",
         0.0 <= rot_ok < 1e-9, f"{rot_ok:.2e}", "0 <= v < 1e-9"),
        ("X3c ... control: rotation w/o transport is NOT",
         np.isfinite(rot_ctrl) and rot_ctrl > 0.1, f"{rot_ctrl:.3f}", "> 0.1"),
        ("X3c ... and zeroing the rot columns BREAKS it",
         np.isfinite(rot_bad) and rot_bad > 0.1, f"{rot_bad:.3f}", "> 0.1"),
        ("X3c no doweled column is identically zero",
         np.isfinite(colmin) and colmin > 1e-6, f"{colmin:.3e}", "> 1e-6"),
        ("X4  where solvable, continuation works BOTH ways", two_all,
         f"{sum(twoside.values())}/{len(twoside)}", "all"),
        ("X4  the 12-around-1 absence is TWO-SIDED too",
         bool(absent) and all(absent.values()),
         f"{sum(absent.values())}/{len(absent)}", "all"),
        ("X5  exactly 24 constant-length intra-unit spans",
         r5["n_const"] == 24, str(r5["n_const"]), "24"),
        ("X5  CONST_TOL above the largest CONSTANT range",
         r5["n_const"] > 0 and r5["const_hi"] < CONST_TOL,
         f"{r5['const_hi']:.2e}", f"< {CONST_TOL:.0e}"),
        ("X5  ... and below the smallest VARYING one (bound ABOVE)",
         np.isfinite(r5["vary_lo"]) and r5["vary_lo"] > CONST_TOL,
         f"{r5['vary_lo']:.2e}", f"> {CONST_TOL:.0e}"),
        ("X5  every brentq refinement converged (no silent loss)",
         n_fail_total == 0, str(n_fail_total), "0"),
        ("X5  exactly 6 intra-unit spans reach strut length",
         n_intra == 6, str(n_intra), "6"),
        ("X5  ... at a SINGLE angle", intra_one_angle,
         str(len(r5["distinct"])), "1"),
        ("X5  ... which is the icosahedral phase", intra_at_ico,
         f"{r5['ico']:.9f}", f"{A_ICO:.9f}"),
        ("X5  ... CONTROL: rejects A_ICO offset by 1e-3 deg",
         np.isfinite(r5["aico_ctrl"]) and r5["aico_ctrl"] > TOL["aico"],
         f"{r5['aico_ctrl']:.2e}", f"> {TOL['aico']:.0e}"),
        ("X5  ... and TOL[aico] is inside its derived band",
         AICO_RECORD_QUANTUM < TOL["aico"] < AICO_CONTROL_OFFSET,
         f"{TOL['aico']:.0e}",
         f"{AICO_RECORD_QUANTUM:.0e}..{AICO_CONTROL_OFFSET:.0e}"),
        ("X5  12 square diagonals at a=0, 6 of them shorten",
         r5["n_square_diag"] == 12 and r5["n_diag"] == 6,
         f"{r5['n_square_diag']}/{r5['n_diag']}", "12/6"),
        ("X5  ... and that derived set IS the crossing set",
         r5["diag_same"], str(r5["diag_same"]), "True"),
        ("X5b k-table vs T2 22682 (outside source)", r5["ktab_ok"],
         f"{max(d for _, _, _, d, _, _ in r5['ktab']):.2e}",
         f"< {TOL['ktable']:.0e}"),
        ("X5b ... and the angle MOVES with k (non-vacuity)",
         r5["kspan"] > 5.0, f"{r5['kspan']:.4f}", "> 5.0"),
        ("X5c every cluster: 6 intra crossings/unit, at a_ico",
         per_unit_ok, str(per_unit_ok), "True"),
        ("X5c inter-unit spans reaching strut length exist",
         inter_nonempty, str(inter_nonempty), "True"),
        ("X5c STAR has inter-unit spans taut AT a_ico",
         star_at_ico, str(star.get("n_at_ico", -1)), "> 0"),
        ("X5c full boxes contain every mediating site (forced)", box_forced,
         f"0 of {box_pairs}", "0 missing"),
        ("X5c ... so their 0/0/0 is an identity, printed not gated",
         lattice_clean,
         f"{cube8m.get('n_at_ico', -1)}/{cube8r.get('n_at_ico', -1)}/"
         f"{cube27.get('n_at_ico', -1)}", "0/0/0"),
        ("X5c VACANCY, not boundary: absent 1, present 0", vacancy_ok,
         f"{vac_absent}/{vac_present}/{vac_holed}/{vac_rbasis}", "1/0/1/0"),
        ("X5c BINDING competitor is INSTALLABLE at a_ico",
         bind_installable, f"{bind_worst:.4f} x strut",
         f"< {r5['contact_spacing'] / STRUT_LEN:.4f}"),
        ("X5c ... and it binds ABOVE a_ico (lock angle MOVES)",
         bind_above, f"{bind_lowest_a:.6f}", f"> {r5['ico']:.6f}"),
        ("X5c ... in every cluster (count of short competitors)",
         n_short_all > 0, str(n_short_all), "> 0"),
        ("X5c ONE sub-icosa angle, same in all M clusters", m_below_ok,
         f"{a_below:.9f}", "one, < a_ico"),
        ("X5c ... and NONE in the R basis (basis-dependent)", r_below_ok,
         str(cube8r.get("n_below", -1)), "0"),
        ("X5d FINE-step independent (absolute second step)",
         r5["step_counts"][0] == r5["step_counts"][1] and r5["step_dev"] < 1e-9,
         f"{r5['step_dev']:.1e}/{r5['step_counts'][0]}v"
         f"{r5['step_counts'][1]}", "< 1e-9, equal"),
        ("X5d COARSE-step independent (the grid that can LOSE)",
         r5["coarse_counts"][0] == r5["coarse_counts"][1]
         and r5["coarse_dev"] < 1e-9,
         f"{r5['coarse_dev']:.1e}/{r5['coarse_counts'][0]}v"
         f"{r5['coarse_counts'][1]}", "< 1e-9, equal"),
        ("X6  mirror unit == diag(1,-1,-1) * unit",
         r6["mirror_dev"] < TOL["mirror"], f"{r6['mirror_dev']:.2e}",
         f"< {TOL['mirror']:.0e}"),
        ("X6  ... control against a WRONG rotation is large",
         r6["ctrl_min"] > 0.1, f"{r6['ctrl_min']:.4f}", "> 0.1"),
        ("X6  mixed chirality solves (ENTAILED by the row above)",
         r6["mixed_ok"], f"{r6['mixed_worst']:.2e}", f"< {SOLVE_TOL:.0e}"),
        ("X7  plate contact absent, or >5 deg from a_ico", clearance_far,
         "none" if r7["no_contact"] else f"{r7['first_contact']:.4f}",
         "none or >5"),
        ("X7  clearance monotone through a_ico (no feature)",
         r7["monotone"], str(r7["monotone"]), "True"),
        ("X7  ... and the statistic MOVES (non-vacuity)",
         r7["ratio"] > 5.0,
         "inf" if not np.isfinite(r7["ratio"]) else f"{r7['ratio']:.3f}", "> 5"),
        ("X7  faces never win: plate-plate <= strut-strut, bitwise",
         r7["faces_irrelevant"] == 0.0, f"{r7['faces_irrelevant']:.2e}",
         "exactly 0"),
        ("GUARD  the collapse guard actually FIRED somewhere",
         guard_fired, str(guard_fired), "True"),
    ]
    print()
    print("=" * 78)
    print(f"GATE  {len(checks)} rows: every check's verdict, and this process's "
          f"exit code")
    print("=" * 78)
    for name, passed, val, crit in checks:
        print(f"  {'PASS' if passed else 'FAIL':4s}  {name:54s} "
              f"{val:>18s} {crit:>16s}")
    print()
    print("  ROWS THAT EXIST TO STOP A CHECK FROM BEING UNFALSIFIABLE:")
    print("   * 'the angle MOVES with k' -- without it, agreement with the")
    print("     recorded k-table at k = 1.00 would be satisfied by an enumerator")
    print("     that returned a constant.")
    print("   * 'control against a WRONG rotation is large' -- without it, the")
    print("     mirror row would be satisfied by a comparison that always")
    print("     returns zero.")
    print("   * 'exactly 24 constant-length spans' -- without it, a labelling")
    print("     error that mislabelled struts would silently change which spans")
    print("     the ranking scans.")
    print("   * 'twelve-around-one closes ONLY at a = 0' -- the companion to")
    print("     'six-around-one closes at every a'. Either alone is satisfiable")
    print("     by a solver that always succeeds or always fails.")
    print("   * 'ANTI is an involution / fixed-point-free' -- without them the")
    print("     two central-symmetry rows above are satisfied by the TRIVIAL")
    print("     SELF-PAIRING, which measures |v - v| = 0 and passes.")
    print("   * 'the control rejects A_ICO offset by 1e-3' -- without it")
    print("     TOL['aico'] is unbounded above and the headline angle can be")
    print("     wrong by two thirds of a degree with the whole gate green.")
    print("   * 'rotation w/o transport is NOT in the kernel' and 'zeroing the")
    print("     rotation columns BREAKS it' -- without them X3c is satisfied by")
    print("     a Jacobian that annihilates everything, and the column-deletion")
    print("     class stays open in the three variables X3b cannot see.")
    print("   * 'CONST_TOL below the smallest VARYING range' -- the ABOVE half")
    print("     of a band whose BELOW half was the only one gated.")
    print("   * 'COARSE-step independent' -- the grid that can delete a finding")
    print("     outright, and did, while the gated grid could not.")
    print()
    print("  ROWS DELIBERATELY REMOVED IN REVISION, because they could not fail:")
    print("   * 'N=1 doweled nullity is 7'. N=1 has no contacts, so the doweled")
    print("     Jacobian is 0 x 7 and `rank_of` returns 0 for an empty matrix:")
    print("     the row evaluated 7 - 0 == 7, an arithmetic identity, with the")
    print("     word 'nullity 7' hardcoded in the printed line beside it. The")
    print("     SHAPE is gated instead, and X3c gates that the seven variables")
    print("     reach the constraints once contacts exist.")
    print("   * 'contact rotation block has rank 2, not 3'. `rank_of(-_hat(u))`")
    print("     is 2 for any nonzero u; a 3x3 skew-symmetric matrix cannot have")
    print("     rank 3. It is a property of `_hat`, not of this array, and it is")
    print("     now stated in X6 as the derivation it always was.")
    print("   * 'the in-phase mode is not the zero vector' (|z| > 1). Every unit")
    print("     carries 1.0 in its phase slot, so |z| >= sqrt(n) >= sqrt(2) for")
    print("     any topology built here, and the value read was the LAST")
    print("     topology's rather than the minimum. Replaced by the BREATHING")
    print("     norm, which is zero exactly when the lattice is held rigid.")
    print("   * 'the FULL lattices have none there' is retained but RELABELLED:")
    print("     it is an identity of rectangular boxes (audited above), not a")
    print("     measurement, and the falsifiable version is the VACANCY row.")
    print()
    print("  AND ONE ROW LABELLED AS ENTAILED RATHER THAN INDEPENDENT: 'mixed")
    print("  chirality solves' follows from 'mirror unit == a rotation of the")
    print("  unit' -- a mirrored unit IS a rotated unit, so every sign pattern is")
    print("  the all-same problem re-seeded. It is a consistency check on the")
    print("  implementation, not a third leg of the Q5 argument.")
    print()
    print("  A ROW DELIBERATELY NOT BUILT. The coordinator's 'ruler test' (at the")
    print("  lock, are the six non-strut spans the same length as the struts?)")
    print("  COULD NOT FAIL: omnitriangulation and diagonal-reaches-strut-length")
    print("  are the same instant by construction, so both candidate mechanisms")
    print("  predicted a lock at exactly the icosahedral phase and the test")
    print("  discriminated nothing. It confirmed the ANGLE only. No gate row of")
    print("  that shape appears above; the X5 rows assert the ABSENCE of")
    print("  competing spans, which is a different and falsifiable claim.")

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
    print("jb_x_array_linkage -- the first multi-unit jitterbug model")
    print("=" * 78)
    print("  bead inviscid-qvf.11. Units joined at SINGLE VERTICES, in phase.")
    print("  Only the vertex identifications are imposed; the lattice spacing is")
    print("  never fixed by hand.")

    if _PAIRING is None:
        print()
        print("=" * 78)
        print(f"GATE  1 row: the hinge pairing could not be read at "
              f"a = {PAIRING_PROBE_DEG}")
        print("=" * 78)
        print(f"  FAIL  X0  hinge pairing readable at the probe angle       "
              f"{'unreadable':>18s} {'12 x mult 2':>16s}")
        print()
        print("  The pairing must be read where the configuration is 12 shared")
        print("  vertices of multiplicity 2. At a = 60, 120, 240, 300 the twelve")
        print("  merge into six of multiplicity 4 and there is no such pairing to")
        print("  read. Nothing below could be computed, so nothing below is")
        print("  printed -- and this arrives as a FAIL ROW rather than as a")
        print("  traceback, because a traceback destroys the verdict table and")
        print("  leaves a reader unable to tell a broken build from a measured")
        print("  result. (This branch exists because a mutation probe took the")
        print("  probe angle to 60 and got a traceback.)")
        return 1

    topos = build_topologies()
    r0 = x0_control()
    r1 = x1_topology()
    x2, guard_fired = x2_existence(topos)
    (ranks, sig_free, sig_dow, inph, breathe, inph_bad,
     rot_ok, rot_ctrl, rot_bad, colmin) = x3_rank(topos, x2)
    term, twoside, absent = x4_distinguished(topos, x2)
    r5 = x5_spans(topos)
    r6 = x6_chirality(topos)
    r7 = x7_clearance()
    x8_verdict(r0, r1, term, twoside, ranks, r5, r6, r7)
    return gate(r0, r1, guard_fired, ranks, sig_free, sig_dow, inph,
                breathe, inph_bad, rot_ok, rot_ctrl, rot_bad, colmin,
                term, twoside, absent, r5, r6, r7)


if __name__ == "__main__":
    # numpy's matmul on this platform's BLAS raises spurious divide/overflow
    # warnings whose text names lines that do no division at all. They are
    # suppressed so the output stays byte-identical across runs. Every place
    # where a non-finite value could actually MATTER carries an explicit
    # np.isfinite check that breaks the iteration instead of propagating, and
    # every quantity in the gate is compared against a finite threshold, so a
    # genuine nan reaches the table as a FAIL rather than as a warning nobody
    # reads.
    with np.errstate(all="ignore"):
        sys.exit(main())
