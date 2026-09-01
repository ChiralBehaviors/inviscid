"""joint_twist -- G2: the relative-rotation lift across BOTH co-located
vertex joints, over one full overdrive cycle. The answer to section 4's
"if twist accumulates anywhere it is at the vertex joints": it does not
accumulate there either. Every relative loop lifts to +q.

THE QUESTION (note section 5 G2, amended by 8.4; bead inviscid-v3q).
Section 4 argued the face welds cannot hold twist (turn = a - b + 60 = 0)
and pointed at the vertex joints as the place a cycle's 2 pi would
physically sit. Section 8.4 added: "the vertex joint" is not singular --
one point carries TWO independent joints, {O,A} and {B,D} (vertex_point
R1/R2), and both must be measured. G3 (braid class of the joint swap) was
dropped as ill-posed -- the exchange is collinear, there is no transverse
framing to braid in (8.4) -- and folds into this gate: the well-posed
content of "does the swap twist anything" is the lift class of the
relative rotation across each joint, read here.

WHAT IS MEASURED HERE. The four cells of the vertex_point configuration
(O, A, B, D around the vertex at L*(1,1,0)), driven through the full
360-degree coherent cycle at G1's sampling grid. The eight plate corners
on the vertex are grouped by vertex_point's own tools; each involved
plate's rotation path is fit from its matched corners and lifted with the
G0 instrument. For every cross pair (a plate of one cell, a plate of its
joint partner) the RELATIVE path R2(a) R1(a)^T is closed and lifted.
Controls: the absolute lifts (must reproduce G1's -q right here), pairs
NOT tied by any joint, and a cell-void face-weld pair on the ring.

  R1  The setup is the record's and the absolute fact holds here: the
      eight corners group {O,A} and {B,D} at reference, and each of the
      eight involved plates' absolute path lifts to -q on this patch.
  R2  ACROSS BOTH JOINTS THE RELATIVE LOOP IS TRIVIAL: all four cross
      pairs of {O,A} and all four of {B,D} close in SO(3), wander to a
      FULL 180 degrees mid-cycle, and every one lifts to +q. The twist
      one cycle would deposit in a joint cancels pairwise.
  R3  The cancellation is the composition law, not joint geometry:
      sign(relative) = sign(abs) * sign(abs) for every pair measured,
      including pairs no joint ties -- ANY two plates cancel. The Z2 is
      a fact about ABSOLUTE orientation histories only; no pairwise
      relative history carries it.
  R4  The contrast at face welds: across a cell-void face weld the
      relative history is CONSTANT through the whole cycle (excursion at
      the arccos floor) -- section 5's turn law extended through all 360
      degrees. Vertex joints wander to 180 degrees and still lift +q.

WHAT THIS MEANS FOR THE PROPOSAL (sections 4, 8.2, 8.6). Stronger than
8.2's "point joints cannot store the Z2": there is no Z2 in any relative
history TO store. A redesigned, tether-capable joint (the spinor-linkage
route) tied plate-to-plate would return UNTWISTED after one cycle; the -q
is visible only against an external frame -- a tether to the lab, not to
a neighbour. The SU(2)-medium proposal's "one cycle leaves every joint
2 pi-twisted" is measured false in the relative sense that a joint could
ever feel.

T2: [23865]. Ref: su2_boundary_conditions.md section 8 (read before 1-7).
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model import cluster as CL
from analysis.model import plates as Z
from analysis.model.first_principles import geometry as G
from analysis.model.first_principles.vertex_point import NAMES, SITES4, corners_at_vertex
from analysis.model.su2 import lift as SU
from analysis.model.su2.plate_holonomy import coherent_X, grid

RING = [(0, 0, 0), (1, 1, 1), (2, 2, 0), (1, 1, -1)]


def vertex_groups(at, a_ref=-30.0):
    """The two coincidence groups of the eight tracked corners, as (cell, face)
    plate ids, read at a generic angle -- vertex_point R2's grouping with the
    member plates kept."""
    L = RC.lattice_constant(a_ref)
    C = Z.corners(a_ref)
    pts = [C[f][c] + L * np.array(SITES4[k], float) for (k, f, c) in at]
    cls = []
    for i, p in enumerate(pts):
        for g in cls:
            if np.linalg.norm(pts[g[0]] - p) < 1e-9:
                g.append(i)
                break
        else:
            cls.append([i])
    assert len(cls) == 2 and all(len(g) == 4 for g in cls)
    return [sorted({(at[i][0], at[i][1]) for i in g}) for g in cls]


def plate_paths(plates, degs):
    """Each plate's rotation path relative to the cycle start, fit from its
    three matched corners under the model's spacing at every angle."""
    C0 = {kf: Z.corners(degs[0])[kf[1]]
          + RC.lattice_constant(degs[0]) * np.array(SITES4[kf[0]], float) for kf in plates}
    R = {kf: np.empty((len(degs), 3, 3)) for kf in plates}
    for n, a in enumerate(degs):
        C, L = Z.corners(a), RC.lattice_constant(a)
        for kf in plates:
            R[kf][n] = G.kabsch(C0[kf], C[kf[1]] + L * np.array(SITES4[kf[0]], float))[0]
    return R


def rel(R2, R1):
    return np.einsum("nij,nkj->nik", R2, R1)


def max_angle(Rs):
    tr = np.einsum("nii->n", Rs)
    return float(np.degrees(np.arccos(np.clip((tr - 1.0) / 2.0, -1.0, 1.0))).max())


def gate():
    checks = []
    A = checks.append
    degs = grid()
    at = corners_at_vertex()
    gOA, gBD = sorted(vertex_groups(at), key=lambda g: g[0][0])
    R = plate_paths(sorted(set(gOA) | set(gBD)), degs)

    def name(kf):
        return f"{NAMES[kf[0]]}{kf[1]}"

    def cross(g):
        k1 = min(k for k, _f in g)
        return [(p, q) for p in g if p[0] == k1 for q in g if q[0] != k1]

    # ---- R1: the setup is the record's; the absolute fact holds here ----------
    absh = {kf: SU.holonomy(R[kf], max_step_deg=1.0) for kf in R}
    grp_ok = ({NAMES[k] for k, _ in gOA} == {"O", "A"} and
              {NAMES[k] for k, _ in gBD} == {"B", "D"})
    A(("R1  THE SETUP IS THE RECORD'S AND THE ABSOLUTE FACT HOLDS HERE: the eight "
       "corners on the vertex group {O,A} and {B,D} at reference (vertex_point R2's "
       "permanent pairing), and each of the eight involved plates' absolute "
       "orientation path over the cycle lifts to -q on this four-cell patch -- G1's "
       "measurement reproduced on the drive this gate rides.",
       grp_ok and all(v == -1 for v in absh.values()),
       "groups {" + ",".join(name(kf) for kf in gOA) + "} {"
       + ",".join(name(kf) for kf in gBD) + "}; abs signs "
       + " ".join(f"{name(kf)}:{v:+d}" for kf, v in absh.items()),
       "{O,A},{B,D}; all -1"))

    # ---- R2: across both joints the relative loop is trivial ------------------
    got2, ok2 = {}, True
    worst_step = 0.0
    for jn, g in (("OA", gOA), ("BD", gBD)):
        rows = []
        for p, q in cross(g):
            Rr = rel(R[q], R[p])
            qs, ws = SU.lift(Rr, max_step_deg=2.0)
            worst_step = max(worst_step, ws)
            h = SU.holonomy(Rr, max_step_deg=2.0)
            ok2 &= (h == +1)
            rows.append(f"{name(p)}x{name(q)}:{h:+d}@{max_angle(Rr):.0f}deg")
        got2[jn] = " ".join(rows)
    ok2 &= worst_step <= 1.0 + 1e-6
    A(("R2  ACROSS BOTH JOINTS THE RELATIVE LOOP IS TRIVIAL: every cross pair of "
       "{O,A} and of {B,D} closes in SO(3), wanders to a FULL 180 degrees mid-cycle, "
       "and lifts to +q. One cycle deposits NO net twist in either co-located joint: "
       "what it would wind, it unwinds, pairwise.",
       ok2, got2 | {"worst rel step": f"{worst_step:.2f} deg"},
       "all +1, excursions ~180, step <= 1.0"))

    # ---- R3: the composition law, not joint geometry --------------------------
    controls = [(gOA[0], gBD[0]), (gOA[0], gBD[2]), (gOA[2], gBD[1])]
    got3, ok3 = [], True
    for p, q in [pq for g in (gOA, gBD) for pq in cross(g)] + controls:
        h = SU.holonomy(rel(R[q], R[p]), max_step_deg=2.0)
        ok3 &= (h == absh[p] * absh[q])
        got3.append(f"{name(p)}x{name(q)}:{h:+d}")
    A(("R3  THE CANCELLATION IS THE COMPOSITION LAW, NOT JOINT GEOMETRY: "
       "sign(relative) = sign(abs) * sign(abs) for every pair measured, including "
       "pairs NO joint ties. Any two plates cancel; the Z2 is a fact about absolute "
       "orientation histories only, and no pairwise relative history carries it.",
       ok3, " ".join(got3) + " (last three: unjoined controls)",
       "sign(rel) = product of abs signs, all pairs"))

    # ---- R4: the contrast at face welds ---------------------------------------
    fa = CL.face_along(np.array([1.0, 1.0, 1.0]))
    fb = CL.face_along(np.array([-1.0, -1.0, -1.0]))
    X0 = coherent_X(RING, degs[0])
    Ca0, Cb0 = X0[0][G.TRI[fa]], X0[1][G.TRI[fb]]
    Rw = np.empty((len(degs), 2, 3, 3))
    for n, a in enumerate(degs):
        X = coherent_X(RING, a)
        Rw[n, 0] = G.kabsch(Ca0, X[0][G.TRI[fa]])[0]
        Rw[n, 1] = G.kabsch(Cb0, X[1][G.TRI[fb]])[0]
    exc = max_angle(rel(Rw[:, 1], Rw[:, 0]))
    A(("R4  THE CONTRAST AT FACE WELDS: across a cell-void face weld on the ring the "
       "relative history is CONSTANT through the whole 360 degrees -- section 5's "
       "turn law (a - b + 60 = 0) extended through the full cycle; the excursion is "
       "the arccos conditioning floor, not motion. The vertex joints' relative paths "
       "wander to 180 degrees and STILL lift trivially.",
       exc < 1e-4,
       f"face-weld relative excursion {exc:.1e} deg over the cycle",
       "< 1e-4 deg"))
    return checks


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("joint_twist -- G2: the relative-rotation lift across BOTH co-located "
          "vertex joints is +q")
    print("=" * 78)
    checks = gate()
    bad = 0
    for name, ok, got, want in checks:
        bad += 0 if ok else 1
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")
    print()
    print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
    print("   * SECTION 4'S HOPE IS MEASURED FALSE IN THE ONLY SENSE A JOINT COULD")
    print("     FEEL: the relative loop across every vertex-joint pair is trivial.")
    print("     Even a redesigned, tether-capable joint (the spinor-linkage route,")
    print("     section 8.2) tied plate-to-plate would return UNTWISTED each cycle;")
    print("     the -q exists only against an external frame.")
    print("   * G3 IS SUBSUMED: the braid class was ill-posed (collinear exchange,")
    print("     no transverse framing, 8.4); its well-posed remnant is this relative")
    print("     lift, and it reads +q.")
    print("   * NOTHING HERE SELECTS A BOUNDARY CONDITION: closure for finite")
    print("     simulations remains the G5/G6 lane.")
    print()
    print("  ALL CHECKS PASSED." if not bad
          else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
