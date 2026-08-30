"""ring -- the smallest closed ring, two VEs and the two voids between them:
one motion, the void's fold is the VE's plus 60, the void expands as the VE
closes, and both ends of the motion are octahedra.

THE QUESTION (owner, 2026-08-30). "every octahedral 'void' is really just a
VE that expands." What is the smallest structure in which that is forced?
Not an octahedron hung on a VE face (face_to_face.py: a passenger), but a
void bounded by plates of DIFFERENT VEs in a closed cycle.

WHAT IS MEASURED HERE (`assembly.honeycomb` on VE O(0,0,0), void U(1,1,1),
VE D(2,2,0), void L(1,1,-1) -- four face welds around one vertex point;
the voids carried as jitterbugs whose plates ARE the VEs' plates).
  R1  Nullity 7 (6 rigid + ONE internal) at every fold but a = 0, where it
      is 8 -- the VE end is a singular point of the ring.
  R2  Started at -10 and walked both ways, the one motion is the exchange:
      void fold - VE fold = 60.00 on every step, no cell rotates, the
      spacing follows `lattice_constant(a)`, welds to 1e-14. It reaches the
      octahedron (a -> -60, voids -> 0) and it runs on past the VE.
  R3  Past the VE it runs on only by taking the void to b > 60 -- folded
      past its own octahedron -- where the void's own struts cross. At every
      fold in the physical range they do not.

So "every void is a VE that expands" is exactly what the ring does, and it
is what makes the VE state an END of the array's motion: there the void is
at ITS octahedron, and an octahedron with fixed joints opens one way
(one_cell.py R4). The cell and its void are one jitterbug's worth of motion
split in two: cell 0 -> -60 while void 60 -> 0. Each body's own fold covers
sixty degrees; the full -60 <-> +60 exists, half of it lived by the void.

T2: [23791]. Pages: "The Exchange", "A Full Cycle" (pages/export_ring.py,
pages/export_cycle.py).
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model.first_principles import geometry as G

SITES = [(0, 0, 0), (1, 1, 1), (2, 2, 0), (1, 1, -1)]
NAMES = ["VE O", "void U", "VE D", "void L"]


def nullity(gc):
    asm, _ = RC.honeycomb(SITES, gc=gc)
    q = asm.q0()
    ctr, R, gam, B = asm.frames(q)
    C = asm.constraint_jacobian(asm.cell_jacobians(ctr, R, B))
    r, gap = RC.rank_of(C)
    return C.shape[1] - r, gap


def walk(asm, sign, steps, deg):
    q = asm.q0()
    _, R0, _, _ = asm.frames(q)
    u = np.zeros((asm.N, 7))
    u[0, 6] = 1.0
    rows = []
    for _ in range(steps):
        uu = np.asarray(asm.project_velocity(q, u, momentum=False)[0]).reshape(-1)
        gd = uu.reshape(-1, 7)[:, 6]
        uu = uu / abs(gd[0])
        q = RC.walk(asm, q, uu, sign * np.radians(deg))
        _, Rq, g, _ = asm.frames(q)
        rot = max(np.degrees(np.arccos(np.clip((np.trace(Rq[k] @ R0[k].T) - 1) / 2, -1, 1))) for k in range(4))
        c = q.reshape(-1, 8)[:, 0:3]
        L = np.linalg.norm(c[2] - c[0]) / np.sqrt(8.0)
        rows.append((g.copy(), rot, float(np.abs(asm.weld_residual(q)).max()), L, RC.lattice_constant(g[0])))
    return rows


def void_crossings(gc):
    asm, _ = RC.honeycomb(SITES, gc=gc)
    X = asm.positions(asm.q0())
    S = [(f, X[1][G.TRI[f][i]], X[1][G.TRI[f][j]]) for f in range(8) for i, j in ((0, 1), (1, 2), (2, 0))]
    return G.crossings([s[1] for s in S], [s[2] for s in S], [s[0] for s in S])


def gate():
    checks, out = [], {}
    A = checks.append

    # ---- R1: one internal freedom, singular at the VE ---------------------------
    nul = {gc: nullity(gc)[0] for gc in (-30.0, -10.0, -2.5, 0.0)}
    out["R1"] = nul
    A(("R1  THE RING HAS ONE INTERNAL FREEDOM: nullity 7 (6 rigid + 1) at -30, -10 and "
       "-2.5, and 8 at exactly a = 0 -- the VE end is a singular point of the ring's "
       "weld manifold.",
       nul[-30.0] == 7 and nul[-10.0] == 7 and nul[-2.5] == 7 and nul[0.0] == 8,
       nul, "7, 7, 7, 8"))

    # ---- R2: the exchange, both ways -------------------------------------------
    asm, _ = RC.honeycomb(SITES, gc=-10.0)
    down = walk(asm, -1, 10, 5.0)
    up = walk(asm, +1, 8, 2.5)
    def worst(rows):
        return (max(abs(g[1] - g[0] - 60.0) for g, *_ in rows),
                max(abs(g[3] - g[2] - 60.0) for g, *_ in rows),
                max(abs(g[2] - g[0]) for g, *_ in rows),
                max(r for _, r, *_ in rows), max(res for _, _, res, *_ in rows),
                max(abs(L - law) for _, _, _, L, law in rows))
    wd, wu = worst(down), worst(up)
    out["R2"] = (down[-1][0], up[-1][0], wd, wu)
    A(("R2  THE ONE MOTION IS THE EXCHANGE: walked from -10 both ways, void fold - VE "
       "fold = 60.00 on every step, both VEs at one fold, NO cell rotates, the spacing "
       "follows lattice_constant(a), welds to 1e-14; it reaches the octahedron (VEs -> -60, "
       "voids -> 0) and runs on past the VE.",
       all(w[0] < 1e-6 and w[1] < 1e-6 and w[2] < 1e-6 and w[3] < 1e-8 and w[4] < 1e-12 and w[5] < 1e-6 for w in (wd, wu))
       and down[-1][0][0] < -59.0 and up[-1][0][0] > 9.0,
       f"toward the octahedron: reaches VE {down[-1][0][0]:+.2f} / void {down[-1][0][1]:+.2f}; "
       f"past the VE: reaches VE {up[-1][0][0]:+.2f} / void {up[-1][0][1]:+.2f}; worst |b-a-60| "
       f"{max(wd[0], wu[0]):.1e}, rotation {max(wd[3], wu[3]):.1e}, residual {max(wd[4], wu[4]):.0e}, "
       f"spacing-law {max(wd[5], wu[5]):.1e}",
       "< -59 and > +9; all < 1e-6 / 1e-8 / 1e-12"))

    # ---- R3: past the VE the void folds past its octahedron ------------------------
    cr = {gc: void_crossings(gc) for gc in (-60.0, -30.0, 0.0, 2.5, 5.0, 10.0)}
    out["R3"] = cr
    A(("R3  PAST THE VE THE VOID IS FOLDED PAST ITS OWN OCTAHEDRON and its struts cross; "
       "at every fold in the physical range they do not. The VE state is an END of the "
       "ring's motion because there the void is at ITS octahedron.",
       all(cr[g] == 0 for g in (-60.0, -30.0, 0.0)) and all(cr[g] > 0 for g in (2.5, 5.0, 10.0)),
       cr, "0 at -60/-30/0; > 0 at +2.5/+5/+10"))
    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("ring -- two VEs and the two voids between them: one motion, the exchange")
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
        print("   * A VOID IS A VE THAT EXPANDS -- when it is bounded by plates of")
        print("     different VEs in a closed cycle. Then the folds are tied, b = a + 60.")
        print("   * BOTH ENDS OF THE TIED MOTION ARE OCTAHEDRA (the voids' at a = 0, the")
        print("     cells' at a = -60), so both are dead ends: the body at the VE at either")
        print("     end nominally has the choice; its partner at the octahedron forecloses")
        print("     it. The tied array's range is one sixty-degree segment.")
        print("   * NOTHING PAST a = 0 IS PHYSICAL: overdrive.py explores it deliberately.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
