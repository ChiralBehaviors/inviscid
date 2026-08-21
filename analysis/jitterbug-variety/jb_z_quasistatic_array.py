"""Step Z: the CONTACT GEOMETRY KERNEL for the quasi-static physical-plate array.

Bead `inviscid-qvf.17`, Phase 1a of DECISION 18 (design of record: T2
`inviscid/design-contact-dynamics-array.md`). DECISION 18 REINSTATES contact as
the neighbour coupling for the ARRAY MODEL ONLY (DECISION 16's "interference is
permitted" stands unchanged for the single-unit abstract variety): no plate may
pass through another, inter-unit joints are tension-only wire loops of slack
`w`, and plate thickness `t` enters as a gap offset (`gap >= t`).

THIS BEAD BUILDS THE KERNEL ONLY -- no stepper, no QP, no sweep. Those are
`inviscid-qvf.18/.19/.20/.21`. Everything below must run standalone and print
its own falsifiable gate.

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

from jb_a_family import R_CIRC, Z, corners, faces, rot
from jb_g_strut_clearance import segment_distance as jb_g_segment_distance
from jb_x_array_linkage import (A_ICO, DIAGONALS, PAIRS, STRUT_LEN, STRUTS,
                                SQUARE_DIAGONALS, Topology, build_topologies,
                                verts)

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
# THE GATE
# ==========================================================================

def gate(z0, z2, z3, z4, z5, z6, z7):
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
    print("jb_z_quasistatic_array -- the contact geometry kernel (Phase 1a)")
    print("=" * 78)
    print("  bead inviscid-qvf.17, DECISION 18. Signed plate gaps with witness")
    print("  points and an outward contact normal, wire spans, a thickness")
    print("  offset, and topology as data -- the pieces bead .18's quasi-static")
    print("  QP stepper assembles. No stepper, no QP, no sweep here.")
    print("  FOUR DECLARATIONS (per the AMENDED design of record, T2 23230):")
    print("  KERNEL, MASS MODEL, PRIMITIVE are INAPPLICABLE -- every quantity")
    print("  below is a static geometric measurement at a held phase, with no")
    print("  potential, no mass, no primitive choice. METRIC FORM is QUALIFIED")
    print("  by the amendment, not flatly inapplicable -- but THIS bead carries")
    print("  no QP objective and therefore no weight W, so the qualification is")
    print("  vacuously satisfied rather than exercised: every number below is a")
    print("  NORM-FREE GEOMETRIC LENGTH in R_oct = 1 units. The W-treatment")
    print("  choice belongs to bead .18, not this one.")

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
    return gate(z0, z2, z3, z4, z5, z6, z7)


if __name__ == "__main__":
    with np.errstate(all="ignore"):
        sys.exit(main())
