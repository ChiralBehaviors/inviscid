"""physical_arc -- G7: the physical sixty-degree arc never exercises the
orientation holonomy. The -q of G1 is a fact about the analytic
continuation, bounded away from every physical motion by the crossing
walls. Section 6 bullet 2 made executable (the gate 8.4 said was missing).

THE QUESTION (note sections 6, 8.5; bead inviscid-1lc). G1 measured the
360-degree overdrive cycle as the nontrivial loop of SO(3) for every
plate. But the medium cannot drive that cycle: the physical range is the
sixty degrees between the two octahedra, [-60, 0]. Does ANY closed
physical loop -- the breathe, a Bloch mode, any boundary condition --
exercise the Z2 at all? Expected: no. This gate measures the expectation
instead of asserting it.

WHAT IS MEASURED HERE, on the record's three patches under the coherent
drive (the G1 machinery throughout):

  R1  THE PHYSICAL RANGE IS AN ARC (OR A BOX OF ARCS) WITH DISTINCT
      ENDS: the closed patches (ring, hc15) have constraint nullity 7 --
      six rigid motions plus the breathe, ONE internal freedom. The free
      block has nullity 15, and the gate IDENTIFIES the surplus: its
      eight weld-degree-1 corner VEs each hang by one shared face and
      carry one independent fold (face_to_face's two-independent-folds
      law at the block's corners) -- internal gamma rank 9, rank 1 on
      non-corner bodies (the breathe), rank 8 on the corners. Every
      freedom is a fold ARC, so the accessible space is a box of
      intervals either way, and the ends differ in orientation: every
      plate turns exactly 60 degrees about its own axis from -60 to 0.
  R2  THE ARC IS CROSSING-FREE AND ITS CONTINUATIONS ARE NOT: zero strut
      crossings at every half-degree of [-60, 0] on every patch; within
      ONE degree beyond either end the count is already in the dozens to
      hundreds, and mid-passage in the hundreds to thousands (overdrive
      R4's walls, located). The 2 pi loop is physically untraversable.
  R3  THE CANONICAL PHYSICAL LOOP READS +q: the breathe -60 -> 0 -> -60
      lifts to +q for every plate of every patch, with the orientation
      excursion peaking at exactly 60 degrees and closing at zero net.
      Closed patches retrace the one arc; the block's extra corner folds
      are arcs too, so its accessible space is a simply connected box
      and no loop anywhere in it can wind.
  R4  OSCILLATION NEVER LIFTS, TRAVERSAL DOES: out-and-back spin loops of
      amplitude 30, 90, 170, even 350 degrees all lift +q; the monotone
      360 lifts -q. What separates -q from +q is MONOTONE TRAVERSAL, not
      amplitude -- and monotone traversal is exactly what the crossing
      walls deny the medium. A Bloch mode of any amplitude oscillates,
      so no mode and no boundary condition exercises the holonomy.

WHAT THIS SETTLES (section 8.5's strongest objection, gated). The Z2 of
G1 is real and it is a property of the PARAMETRISATION: the fold circle
of the analytic continuation, not the configuration path of any physical
motion. Together with G2 (no relative history carries the sign), the
SU(2)-medium proposal has no physical carrier in this model: closure
selection for finite simulations rests entirely on the G5/G6 lane.

T2: [23865]. Ref: su2_boundary_conditions.md section 8 (read before 1-7).
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model.first_principles import geometry as G
from analysis.model.first_principles.overdrive import PATCHES
from analysis.model.su2 import lift as SU
from analysis.model.su2.plate_holonomy import coherent_X, drive_rotations

ARC = np.arange(-60.0, 0.0 + 1e-9, 0.5)
OUTSIDE = (-90.0, -61.0, 1.0, 30.0)


def freedoms(sites, gc=-30.0):
    """(nullity, deg-1 body count, internal gamma ranks): the constraint null
    space, with its internal part (rigid motions projected out) read on the
    per-body fold coordinates -- total rank, rank off the weld-degree-1
    bodies, rank on them."""
    asm, deg = RC.honeycomb(sites, gc=gc)
    q = asm.q0()
    ctr, R, gam, B = asm.frames(q)
    C = asm.constraint_jacobian(asm.cell_jacobians(ctr, R, B))
    _, s, Vt = np.linalg.svd(C)
    null = Vt[(s > s[0] * 1e-8).sum():]
    Gq, _ = np.linalg.qr(np.asarray(asm.globals(ctr)).T)
    P = null.T - Gq @ (Gq.T @ null.T)
    Ui, si, _ = np.linalg.svd(P, full_matrices=False)
    internal = Ui[:, si > si[0] * 1e-8].T
    gpart = internal[:, 6::7]
    hang = [i for i in range(len(sites)) if deg[i] == 1]
    rest = [i for i in range(len(sites)) if deg[i] != 1]
    return (null.shape[0], len(hang),
            int(np.linalg.matrix_rank(gpart, tol=1e-8)),
            int(np.linalg.matrix_rank(gpart[:, rest], tol=1e-8)),
            int(np.linalg.matrix_rank(gpart[:, hang], tol=1e-8)) if hang else 0)


def ncross(sites, a):
    P, Q, owner = G.assembly_struts(coherent_X(sites, a))
    return G.crossings(P, Q, owner)


def max_angle(Rs):
    tr = np.einsum("nii->n", Rs)
    return float(np.degrees(np.arccos(np.clip((tr - 1.0) / 2.0, -1.0, 1.0))).max())


def gate():
    checks = []
    A = checks.append
    out_back = np.concatenate([ARC, ARC[::-1][1:]])

    # ---- R1: an arc (or a box of arcs), distinct ends -------------------------
    nul = {n: freedoms(s) for n, s in PATCHES.items()}
    end_angles = {}
    Rob = {}
    for n, sites in PATCHES.items():
        R = drive_rotations(sites, out_back)
        Rob[n] = R
        i0 = len(ARC) - 1                                    # the a = 0 frame
        tr = np.einsum("pii->p", R[:, i0])
        ang = np.degrees(np.arccos(np.clip((tr - 1.0) / 2.0, -1.0, 1.0)))
        end_angles[n] = (float(ang.min()), float(ang.max()))
    ok1 = (nul["ring"] == (7, 0, 1, 1, 0) and nul["hc15"] == (7, 0, 1, 1, 0)
           and nul["block"] == (15, 8, 9, 1, 8)
           and all(abs(lo - 60) < 1e-9 and abs(hi - 60) < 1e-9
                   for lo, hi in end_angles.values()))
    A(("R1  THE PHYSICAL RANGE IS AN ARC (OR A BOX OF ARCS) WITH DISTINCT ENDS: the "
       "closed patches have nullity 7 -- six rigid motions plus the breathe, ONE "
       "internal freedom -- while the free block's nullity 15 is IDENTIFIED: its eight "
       "weld-degree-1 corner VEs each hang by one shared face and carry one "
       "independent fold (face_to_face's law at the block's corners; internal gamma "
       "rank 9, rank 1 off the corners, 8 on them). Every freedom is a fold arc; and "
       "the ends differ in orientation, every plate turning exactly 60 degrees about "
       "its own axis from -60 to 0.",
       ok1,
       "(nullity, deg-1 bodies, gamma ranks total/rest/hanging) " + str(nul)
       + "; end-to-end plate turn min/max "
       + str({n: f"{lo:.6f}/{hi:.6f}" for n, (lo, hi) in end_angles.items()}),
       "ring, hc15 (7,0,1,1,0); block (15,8,9,1,8); all turns exactly 60"))

    # ---- R2: crossing-free arc, walled continuations --------------------------
    worst_arc = {n: max(ncross(s, a) for a in ARC) for n, s in PATCHES.items()}
    walls = {n: {a: ncross(s, a) for a in OUTSIDE} for n, s in PATCHES.items()}
    ok2 = (all(v == 0 for v in worst_arc.values())
           and all(c > 0 for w in walls.values() for c in w.values()))
    A(("R2  THE ARC IS CROSSING-FREE AND ITS CONTINUATIONS ARE NOT: zero strut "
       "crossings at every half-degree of [-60, 0] on every patch; within ONE degree "
       "beyond either end the count is already in the dozens to hundreds, mid-passage "
       "in the hundreds to thousands. The walls overdrive R4 gated, located at the "
       "arc's ends: the 2 pi loop is physically untraversable.",
       ok2,
       f"max on arc {worst_arc}; outside " + str(walls),
       "all 0 on arc; all > 0 outside"))

    # ---- R3: the canonical physical loop reads +q -----------------------------
    got3, ok3 = {}, True
    for n, R in Rob.items():
        h = np.array([SU.holonomy(R[p], max_step_deg=1.0) for p in range(R.shape[0])])
        exc = max(max_angle(R[p]) for p in range(R.shape[0]))
        ok3 &= bool((h == +1).all()) and abs(exc - 60) < 1e-9
        got3[n] = f"{int((h == +1).sum())}/{len(h)} +q, excursion {exc:.6f}"
    A(("R3  THE CANONICAL PHYSICAL LOOP READS +q: the breathe -60 -> 0 -> -60 lifts "
       "to +q for every plate of every patch, the orientation excursion peaking at "
       "exactly 60 degrees and closing at zero net. On the closed patches every "
       "physical loop is a retrace of this one; on the free block the added corner "
       "folds are arcs too (R1), the accessible space a simply connected box, so no "
       "loop anywhere in it can wind.",
       ok3, got3, "all +q, excursion 60"))

    # ---- R4: oscillation never lifts; traversal does --------------------------
    osc = {}
    for amp in (30.0, 90.0, 170.0, 350.0):
        degs = np.concatenate([np.arange(0.0, amp + 1e-9, 0.5),
                               np.arange(amp, 0.0 - 1e-9, -0.5)])
        osc[amp] = SU.holonomy(SU.path("skew", degs), max_step_deg=1.0)
    mono = SU.holonomy(SU.path("skew", np.arange(0.0, 360.0 + 1e-9, 0.5)),
                       max_step_deg=1.0)
    A(("R4  OSCILLATION NEVER LIFTS, TRAVERSAL DOES: out-and-back spin loops of "
       "amplitude 30, 90, 170, even 350 degrees all lift +q; the monotone 360 lifts "
       "-q. What separates the signs is MONOTONE TRAVERSAL, not amplitude -- and "
       "traversal is exactly what the crossing walls deny the medium. A Bloch mode of "
       "any amplitude oscillates: no mode and no boundary condition exercises the "
       "holonomy.",
       all(v == +1 for v in osc.values()) and mono == -1,
       ", ".join(f"A={a:.0f}:{v:+d}" for a, v in osc.items()) + f"; monotone 360:{mono:+d}",
       "all +1; -1"))
    return checks


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("physical_arc -- G7: the physical sixty degrees never exercises the "
          "orientation holonomy")
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
    print("   * SECTION 8.5'S STRONGEST OBJECTION IS NOW A GATED FACT: the Z2 of G1")
    print("     is a property of the fold circle's analytic continuation. No physical")
    print("     motion of the medium -- breathe, mode, or boundary condition -- can")
    print("     traverse the loop that carries the sign.")
    print("   * WITH G2 (no relative history carries it), the SU(2)-medium proposal")
    print("     has NO physical carrier in this model. Closure selection for finite")
    print("     simulations rests entirely on the G5/G6 lane.")
    print("   * THE KINEMATIC FACT STANDS: the continuation IS the nontrivial loop")
    print("     (G1). Nothing here retracts it; this gate bounds where it lives.")
    print()
    print("  ALL CHECKS PASSED." if not bad
          else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
