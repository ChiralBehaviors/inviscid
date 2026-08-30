"""face_to_face -- a VE and an octahedron cell welded on one shared triangle:
their folds are independent, the octahedron is a passenger that rides the
shared plate, and each shared joint carries six corners.

THE QUESTION (owner, 2026-08-30). "now do 2. 1 VE and an adjacent octa cell.
notice anything?" -- and, on seeing it: "there's your problem ;) every
octahedral 'void' is really just a VE that expands. in your animation, the
octa 'void' just rotates with the face." Then: "add a third octa around the
ve."

WHAT IS MEASURED HERE.
  R1  `assembly.cluster(gc=0, sites=[(1,1,1)])`: 14 dof, one three-pair
      weld of rank 6, nullity 8 = 6 rigid + TWO internal. The folds are
      independent. Each shared point carries 2 corners of the VE and 4 of
      the octahedron (the octahedron's two merged joints).
  R2  The whole (a, b) family: the octahedron's pose is the rigid motion
      that puts its welded plate on the VE's (a fit to the three pairs),
      exact to 1e-15, and equal to the cluster builder's on b = a + 60. Its
      turn about the shared axis is a - b + 60 exactly -- zero only on that
      line -- and the centres are ZC (cos a + cos b) apart.
  R3  Nowhere in the family do the two cells' struts cross.
  R4  Two octahedra on adjacent faces of the VE, tied to each other at the
      vertex point they share: still two passengers. Drive the VE either way
      and both stay at 60.00; the tie is satisfied automatically because
      that point is the VE's own vertex.

So an octahedron hung on a VE face does not expand -- there is nothing to
make it. What makes a void expand is being bounded by plates of DIFFERENT
VEs in a closed ring: ring.py.

T2: [23791]. Page: "Face to Face" (pages/export_face_to_face.py).
"""
from __future__ import annotations

import itertools as it
import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model import cell as IC
from analysis.model import plates as Z
from analysis.model.first_principles import geometry as G

AX = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)


def family_pose(a, b, pairs):
    """R, t carrying the octahedron's body corners (fold b) onto the VE's (fold a)."""
    P1 = np.array([IC.cell_verts(a, np.zeros(3))[i] for i, _ in pairs])
    P2 = np.array([IC.cell_verts(b, np.zeros(3))[j] for _, j in pairs])
    R, t = G.kabsch(P2, P1)
    return R, t, float(np.abs((R @ P2.T).T + t - P1).max())


def two_passengers(tie=True):
    sites = [(1, 1, 1), (1, 1, -1)]
    asm = RC.cluster(gc=0.0, sites=sites)
    if not tie:
        return asm, []
    probe = RC.cluster(gc=-30.0, sites=sites)
    Xr = probe.positions(probe.q0())
    pairs, seen = [], []
    for i in range(IC.NV):
        v = Xr[1][i]
        if any(np.linalg.norm(v - u) < 1e-9 for u in seen):
            continue
        for j in range(IC.NV):
            if np.linalg.norm(v - Xr[2][j]) < 1e-9:
                pairs.append((i, j))
                seen.append(v)
                break
    return RC.Assembly(asm.gam0, asm.ctr0, list(asm.welds) + [(1, 2, pairs)]), pairs


def drive(asm, cell, sign, steps=3, deg=10.0):
    q = asm.q0()
    u = np.zeros((asm.N, 7))
    u[cell, 6] = 1.0
    folds = []
    for _ in range(steps):
        uu = np.asarray(asm.project_velocity(q, u, momentum=False)[0]).reshape(-1)
        gd = uu.reshape(-1, 7)[:, 6]
        uu = uu / abs(gd[cell])
        q = RC.walk(asm, q, uu, sign * np.radians(deg))
        folds.append(asm.frames(q)[2].copy())
    return folds, float(np.abs(asm.weld_residual(q)).max())


def gate():
    checks, out = [], {}
    A = checks.append

    # ---- R1: two cells, two folds --------------------------------------------
    asm = RC.cluster(gc=0.0, sites=[(1, 1, 1)])
    q0 = asm.q0()
    ctr, R, gam, B = asm.frames(q0)
    C = asm.constraint_jacobian(asm.cell_jacobians(ctr, R, B))
    rank, _ = RC.rank_of(C)
    (_, _, pairs) = asm.welds[0]
    X = asm.positions(q0)
    C1 = Z.corners(0.0)
    C2 = Z.corners(60.0) + ctr[1]
    census = []
    for i, _ in pairs:
        p = X[0][i]
        n1 = sum(1 for f in range(8) for c in range(3) if np.linalg.norm(C1[f][c] - p) < 1e-9)
        n2 = sum(1 for f in range(8) for c in range(3) if np.linalg.norm(C2[f][c] - p) < 1e-9)
        census.append((n1, n2))
    out["R1"] = (C.shape, rank, census)
    A(("R1  A VE AND AN OCTAHEDRON CELL WELDED ON ONE TRIANGLE HAVE TWO INTERNAL "
       "FREEDOMS: 14 dof, one three-pair weld of rank 6, nullity 8 = 6 rigid + 2. The "
       "folds are independent. Each shared joint carries 2 corners of the VE and 4 of "
       "the octahedron -- an octahedron's vertex is two of its joints on one point.",
       C.shape == (9, 14) and rank == 6 and census == [(2, 4)] * 3,
       f"rows {C.shape[0]}, dof {C.shape[1]}, rank {rank}, nullity {C.shape[1] - rank}; corners per shared point {census}",
       "9, 14, 6, 8; [(2, 4)] x 3"))

    # ---- R2: the whole family ---------------------------------------------------
    ang = np.arange(-60.0, 60.0 + 1e-9, 10.0)
    res_max, turn_err, dist_err, builder_err, kinds = 0.0, 0.0, 0.0, 0.0, 0
    for a in ang:
        for b in ang:
            Rm, t, res = family_pose(a, b, pairs)
            res_max = max(res_max, res)
            turn = G.rot_angle_about(Rm, AX)
            want = (a - b + 60.0 + 180.0) % 360.0 - 180.0
            turn_err = max(turn_err, abs((turn - want + 180.0) % 360.0 - 180.0))
            dist_err = max(dist_err, abs(np.linalg.norm(t) - IC.ZC * (np.cos(np.radians(a)) + np.cos(np.radians(b)))))
            if abs(b - (a + 60.0)) < 1e-9 and -60.0 <= a <= 0.0:
                asm_ab = RC.cluster(gc=float(a), sites=[(1, 1, 1)])
                builder_err = max(builder_err, float(np.abs(asm_ab.ctr0[1] - t).max()), float(np.abs(Rm - np.eye(3)).max()))
                kinds += 1
    out["R2"] = (res_max, turn_err, dist_err, builder_err, kinds)
    A(("R2  ACROSS THE WHOLE (a, b) FAMILY the octahedron's pose is fixed by the shared "
       "plate alone (fit residual ~1e-15), equals the cluster builder's on b = a + 60, "
       "turns about the shared axis by exactly a - b + 60, and sits ZC (cos a + cos b) "
       "from the VE.",
       res_max < 1e-12 and turn_err < 1e-8 and dist_err < 1e-9 and builder_err < 1e-9 and kinds == 7,
       f"fit residual {res_max:.1e}, turn error {turn_err:.1e}, distance error {dist_err:.1e}, "
       f"builder agreement {builder_err:.1e} on {kinds} angles", "< 1e-12, < 1e-8, < 1e-9, < 1e-9, 7"))

    # ---- R3: they never cross ---------------------------------------------------
    worst, nearest = 0, np.inf
    for a in ang:
        for b in ang:
            Rm, t, _ = family_pose(a, b, pairs)
            X2 = (Rm @ Z.corners(b).reshape(-1, 3).T).T + t
            X2 = X2.reshape(8, 3, 3)
            S1, S2 = G.cell_struts(Z.corners(a)), G.cell_struts(X2)
            P = [s[1] for s in S1] + [s[1] for s in S2]
            Q = [s[2] for s in S1] + [s[2] for s in S2]
            owner = [s[0] for s in S1] + [8 + s[0] for s in S2]
            worst = max(worst, G.crossings(P, Q, owner))
            for (f1, p1, q1) in S1:
                for (f2, p2, q2) in S2:
                    if min(np.linalg.norm(x - y) for x in (p1, q1) for y in (p2, q2)) < 1e-9:
                        continue
                    nearest = min(nearest, G.segdist(p1, q1, p2, q2)[0])
    out["R3"] = (worst, nearest)
    A(("R3  NOWHERE IN THE FAMILY DO THE TWO CELLS' STRUTS CROSS. Drive the VE from -60 "
       "to +60 with the octahedron held anywhere: it rides the plate, turning and "
       "sliding, and never expands and never collides.",
       worst == 0 and nearest > 0.03,
       f"crossings {worst}; nearest foreign struts {nearest:.3f}", "0; > 0.03"))

    # ---- R4: two passengers tied to each other ----------------------------------
    tied, tie_pairs = two_passengers(True)
    free, _ = two_passengers(False)
    nul = {}
    for name, a2 in (("tied", tied), ("free", free)):
        q = a2.q0()
        c, r, g, b = a2.frames(q)
        Cj = a2.constraint_jacobian(a2.cell_jacobians(c, r, b))
        nul[name] = Cj.shape[1] - RC.rank_of(Cj)[0]
    stay, resid = 0.0, 0.0
    for sign in (-1, +1):
        folds, res = drive(tied, 0, sign)
        stay = max(stay, max(abs(f[1] - 60.0) for f in folds), max(abs(f[2] - 60.0) for f in folds))
        resid = max(resid, res)
    out["R4"] = (nul, len(tie_pairs), stay, resid)
    A(("R4  TWO OCTAHEDRA ON ADJACENT FACES, TIED TO EACH OTHER AT THE VERTEX POINT THEY "
       "SHARE, ARE STILL TWO PASSENGERS: nullity 8 with the tie (9 without); drive the VE "
       "either way and both stay at 60.00, welds holding. The tie costs nothing because "
       "that point is the VE's own vertex and both octahedra hang on plates containing it.",
       nul["tied"] == 8 and nul["free"] == 9 and len(tie_pairs) == 2 and stay < 1e-6 and resid < 1e-12,
       f"nullity tied {nul['tied']}, free {nul['free']}; tie pairs {len(tie_pairs)}; "
       f"octahedra's folds move by {stay:.1e}; weld residual {resid:.0e}", "8, 9; 2; < 1e-6; < 1e-12"))
    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("face_to_face -- a VE and an octahedron cell on one shared plate: independent "
              "folds, a passenger")
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
        print("   * ONE SHARED PLATE TIES NOTHING BUT THE PLATE. Two folds, independent;")
        print("     the octahedron turns with the plate by a - b + 60 and never expands.")
        print("   * THE LATTICE'S LAW b = a + 60 IS THE ONE LINE WHERE NEITHER CELL TURNS")
        print("     relative to the other. Cycles force the lattice onto it: ring.py.")
        print("   * HANGING MORE OCTAHEDRA ON THE VE CHANGES NOTHING (R4).")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
