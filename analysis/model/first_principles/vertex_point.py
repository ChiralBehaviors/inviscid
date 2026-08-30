"""vertex_point -- one vertex of a VE in the single-covering lattice: four
cells meet there, two permanent joints share the point, and the tied array
folds one way from the VE because the other way drives those two joints
head-on through each other.

THE QUESTION (owner, 2026-08-30). "make the 8 octa around 1 ve. in this
config, what is the answer you get for the ve" (single covering ONLY).
Then "i believe the correct answer is 4 ... 4 triangles joined at the
vertex." Then, of a first animation that showed the tie sets swapping across
the VE: "at the 20 degree mark, the vertexes all have 4 triangles. these
triangles *never* separate ... that is incorrect." He was right; the error
was a mirrored spacing on the a > 0 frames. Everything below uses the
model's own spacing at every fold.

WHAT IS MEASURED HERE (27 VEs, the all-even sites of {-2,0,2}^3,
`honeycomb_single`; the centre's vertex at L*(1,1,0) with cells O=(0,0,0),
A=(2,0,0), B=(0,2,0), D=(2,2,0)).
  R1  At the VE every one of the centre's twelve vertices is met by FOUR
      cells -- itself, two axis neighbours, the face-diagonal neighbour --
      so eight plate corners sit on the point, and all eight voids around
      the VE are exact octahedra with strut-length edges. The model's held
      ties join four of the eight: one axis neighbour per vertex.
  R2  The joints are permanent: with the model's spacing at every a, the
      eight corners are {O+A} and {B+D} at every a except 0, where all
      eight coincide.
  R3  The tied 27-cell block walks THROUGH the VE on its weld manifold,
      both ways (nullity 7 at -10, -2.5 and 0; residual 1e-16), with its
      spacing following `lattice_constant(a)` on both sides.
  R4  What stops it is the two joints: they separate along one line, the
      same distance at -a and +a, in exactly opposite directions -- past the
      VE they have passed head-on through each other -- and the struts
      hanging off them cross from +0.5 on. On the physical side they open
      a gap.

So on the single covering the array folds from the VE one way, the way its
ties were laid, and the far side is two rubber joints trying to occupy one
point. Seen from the void, ring.py: at the VE the void is at ITS octahedron.

T2: [23789]. Pages: "One Vertex, Four Cells", "Two Joints, One Point"
(pages/export_vertex.py, pages/export_joints.py).
"""
from __future__ import annotations

import itertools as it
import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model import cell as IC
from analysis.model import plates as Z
from analysis.model.first_principles import geometry as G

SITES27 = list(it.product((-2, 0, 2), repeat=3))
SITES4 = [(0, 0, 0), (2, 0, 0), (0, 2, 0), (2, 2, 0)]
NAMES = ["O", "A", "B", "D"]


def corners_at_vertex():
    """The eight (cell, plate, corner) that sit on the vertex L*(1,1,0) at a = 0."""
    L0 = RC.lattice_constant(0.0)
    vtx = L0 * np.array([1.0, 1.0, 0.0])
    C0 = Z.corners(0.0)
    at = [(k, f, c) for k, s in enumerate(SITES4) for f in range(8) for c in range(3)
          if np.linalg.norm(C0[f][c] + L0 * np.array(s, float) - vtx) < 1e-9]
    assert len(at) == 8
    return at


def groups_at(a, at):
    """Coincidence groups of the tracked corners with the model's spacing at a."""
    L = RC.lattice_constant(a)
    C = Z.corners(a)
    cells = [C + L * np.array(s, float) for s in SITES4]
    pts = [cells[k][f][c] for (k, f, c) in at]
    cls = []
    for i, p in enumerate(pts):
        for g in cls:
            if np.linalg.norm(pts[g[0]] - p) < 1e-9:
                g.append(i)
                break
        else:
            cls.append([i])
    return [tuple(sorted({NAMES[at[i][0]] for i in g})) for g in cls], cells


def block(gc):
    asm, _ = RC.honeycomb_single(SITES27, gc=gc)
    return asm


def internal_tangent(asm, q):
    ctr, R, gam, B = asm.frames(q)
    C = asm.constraint_jacobian(asm.cell_jacobians(ctr, R, B))
    _, s, Vt = np.linalg.svd(C)
    null = Vt[(s > s[0] * 1e-8).sum():]
    Gq, _ = np.linalg.qr(np.asarray(asm.globals(ctr)).T)
    P = null.T - Gq @ (Gq.T @ null.T)
    u = np.linalg.svd(P)[0][:, 0]
    return u / abs(u.reshape(-1, 7)[:, 6].mean())


def gate():
    checks, out = [], {}
    A = checks.append
    L0 = RC.lattice_constant(0.0)

    # ---- R1: the census at the VE -----------------------------------------------
    asm, deg = RC.honeycomb_single(SITES27, gc=0.0)
    X = asm.positions(asm.q0())
    origin = SITES27.index((0, 0, 0))
    per_vertex = []
    for i in range(IC.NV):
        v = X[origin][i]
        cells = {k for k in range(asm.N) for j in range(IC.NV) if np.linalg.norm(X[k][j] - v) < 1e-9}
        per_vertex.append(len(cells))
    held = sum(len(p) for k, l, p in asm.welds if origin in (k, l))
    allpts = X.reshape(-1, 3)
    octa = 0
    for s in it.product((-1, 1), repeat=3):
        c = L0 * np.array(s, float)
        verts = [c + L0 * e for e in np.eye(3)] + [c - L0 * e for e in np.eye(3)]
        present = all(np.min(np.linalg.norm(allpts - w, axis=1)) < 1e-9 for w in verts)
        edges = sorted(np.linalg.norm(p - q) for p, q in it.combinations(verts, 2))[:12]
        octa += int(present and np.allclose(edges, IC.EL))
    out["R1"] = (sorted(set(per_vertex)), held, deg[origin], octa)
    A(("R1  AT THE VE EVERY VERTEX OF A CELL IS MET BY FOUR CELLS -- itself, two axis "
       "neighbours, the face-diagonal neighbour between them -- so EIGHT plate corners sit "
       "on one point; the eight voids around it are exact octahedra with strut edges; "
       "the model's held ties join four of the eight (one axis neighbour per vertex, "
       "twelve pairs, degree six).",
       set(per_vertex) == {4} and held == 12 and deg[origin] == 6 and octa == 8,
       f"cells per vertex {sorted(set(per_vertex))}, held pairs {held}, degree {deg[origin]}, octahedral voids {octa}/8",
       "{4}, 12, 6, 8"))

    # ---- R2: the joints are permanent ------------------------------------------
    at = corners_at_vertex()
    grp = {a: groups_at(a, at)[0] for a in (-60.0, -30.0, -10.0, -2.5, 0.0, 2.5, 10.0)}
    ref = sorted(grp[-30.0])
    out["R2"] = grp
    A(("R2  THE JOINTS ARE PERMANENT. With the model's spacing at every a the eight "
       "corners are {O+A} and {B+D} at every fold except a = 0, where all eight coincide; "
       "the pairing never changes across the VE.",
       all(sorted(grp[a]) == ref for a in grp if a != 0.0) and ref == [("A", "O"), ("B", "D")]
       and grp[0.0] == [("A", "B", "D", "O")],
       {k: v for k, v in grp.items()}, "{O,A},{B,D} everywhere but 0; all four at 0"))

    # ---- R3: the tied block walks through the VE ---------------------------------
    nul = {}
    for gc in (-10.0, -2.5, 0.0):
        b = block(gc)
        q = b.q0()
        ctr, R, gam, B = b.frames(q)
        C = b.constraint_jacobian(b.cell_jacobians(ctr, R, B))
        nul[gc] = C.shape[1] - RC.rank_of(C)[0]
    b0 = block(0.0)
    q0 = b0.q0()
    u = internal_tangent(b0, q0)
    i0, i1 = SITES27.index((0, 0, 0)), SITES27.index((2, 0, 0))
    walked = {}
    for sign in (-1, +1):
        q = q0.copy()
        for _ in range(3):
            q = RC.walk(b0, q, u, sign * np.radians(2.5))
        g = b0.frames(q)[2]
        c = q.reshape(-1, 8)[:, 0:3]
        sep = np.linalg.norm(c[i1] - c[i0]) / 2.0
        walked[sign] = (float(g.mean()), float(g.std()), float(sep), float(RC.lattice_constant(g.mean())),
                        float(np.abs(b0.weld_residual(q)).max()))
    out["R3"] = (nul, walked)
    ok3 = (all(v == 7 for v in nul.values())
           and all(w[1] < 1e-9 and abs(w[2] - w[3]) < 1e-6 and w[4] < 1e-12 for w in walked.values())
           and (walked[-1][0] < -5 < 5 < walked[+1][0] or walked[+1][0] < -5 < 5 < walked[-1][0]))
    A(("R3  THE WELDS DO NOT STOP IT: nullity 7 (6 rigid + the breathe) at -10, -2.5 and 0, "
       "and a finite walk from the VE along the breathe tangent converges on BOTH sides "
       "(residual 1e-16, all cells at one fold), the spacing following lattice_constant(a) "
       "continued smoothly through zero.",
       ok3, {k: f"fold {v[0]:+.2f} (spread {v[1]:.0e}), spacing {v[2]:.5f} vs law {v[3]:.5f}, residual {v[4]:.0e}"
             for k, v in walked.items()} | {"nullity": nul},
       "7 everywhere; both walks past +-5 deg with spacing = law"))

    # ---- R4: the two joints, head-on ----------------------------------------------
    def joints_sep(a):
        _, cells = groups_at(a, at)
        pO = [cells[k][f][c] for (k, f, c) in at if k == 0][0]
        pB = [cells[k][f][c] for (k, f, c) in at if k == 2][0]
        return pB - pO
    dm, dp = joints_sep(-7.5), joints_sep(7.5)
    cosang = float(dm @ dp / (np.linalg.norm(dm) * np.linalg.norm(dp)))
    cross = {}
    for a in (-2.5, 2.5):
        _, cells = groups_at(a, at)
        P, Q, owner = [], [], []
        for k in range(4):
            for f, p, q in G.cell_struts(cells[k]):
                P.append(p)
                Q.append(q)
                owner.append(k * 8 + f)
        cross[a] = G.crossings(P, Q, owner)
    out["R4"] = (np.linalg.norm(dm), np.linalg.norm(dp), cosang, cross)
    A(("R4  WHAT STOPS IT IS THE TWO JOINTS ON THE POINT: they separate along one line, "
       "the same distance at -a and +a, in exactly opposite directions (cos = -1) -- past "
       "the VE they have passed head-on through each other -- and the struts hanging off "
       "them cross from the first frame past zero, while on the physical side they open a "
       "gap and nothing crosses.",
       abs(np.linalg.norm(dm) - np.linalg.norm(dp)) < 1e-9 and abs(cosang + 1.0) < 1e-9
       and cross[-2.5] == 0 and cross[2.5] > 0,
       f"joints {np.linalg.norm(dm):.4f} apart at -7.5 and {np.linalg.norm(dp):.4f} at +7.5, cos {cosang:+.4f}; "
       f"crossings at -2.5: {cross[-2.5]}, at +2.5: {cross[2.5]}", "equal; -1; 0; > 0"))
    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("vertex_point -- four cells on one vertex, two permanent joints, and the far "
              "side of the VE")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * '4 TRIANGLES JOINED AT THE VERTEX' IS ONE JOINT: two corners of a cell")
        print("     (its own hinge) plus two of its axis partner. The point holds TWO.")
        print("   * THE TIED ARRAY FOLDS ONE WAY FROM THE VE, the way its ties were laid.")
        print("     Tie O to its x-neighbour and it folds one way; to its y-neighbour, the")
        print("     other. Both chiralities exist; each array has one, chosen when tied.")
        print("   * WHETHER THE RIG'S JOINTS SLIP PAST EACH OTHER is a rig fact the model")
        print("     cannot supply: its joints have no size and no give.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
