"""overdrive -- the ring, and larger patches, driven through a full 360
degrees of the fold: past both dead ends, through the intersections, and
back. An exploration the owner asked for with the physical analogies
suspended; nothing here is the medium.

THE QUESTION (owner, 2026-08-30). "What I'd like to explore with this is
the *overdrive* where we start intersecting in the octas. I'd like to keep
driving that." And: the plates "have sides ... we can track the
orientation." And: "eventually we'll assume that the system drives through
the joint and magically reorients itself if needed." Then: "if we expand
the animated patch, does it look much the same?" and "if we expand the
patch again?"

WHAT IS MEASURED HERE, on three patches -- the four-body ring, the fifteen
bodies of HC15 (a VE, its eight voids, its six axis VEs), and the 3x3x3
block of 27 VEs with its eight interior voids -- at every 30 degrees from
-60 to 300 (= -60 + 360):
  R1  Welds held BY IDENTITY from a = -30 stay satisfied to 1e-15 on every
      frame of every patch. The joints never needed the magic: bodies pass
      through each other and through their own centres; no joint parts.
  R2  The spacing is a pure sinusoid, L(a) = (2/sqrt3) cos(a + 30): largest
      at -30, zero at +60 and +240, negative between -- the lattice turned
      inside out. The physical range [-60, 0] is the sixty degrees near the
      top.
  R3  Sides: every plate's front (fixed to point out of its body at a = 0)
      turns inward when the plate passes through its body's centre --
      voids around +30 and back at +210, cells around +90 and back at
      +270; from 120 to 180 every plate in the array is inside out.
  R4  Strut crossings are confined to the passages: none at any multiple
      of 60, hundreds mid-passage, the count scaling with the patch.
  R5  At +60 and +240 the spread of every centre about the mean is zero:
      the whole patch, cells and voids alike, is ONE coincident octahedron;
      and 300 is -60 again, position for position.

None of it depends on the size of the patch: the coherent motion is
translation-invariant, so an infinite array driven this way is a homothety
to a point and back, mirrored each time, period 360, with the sixty degrees
we can build sitting at the top of the spacing curve between the two
collapses.

T2: [23794] [23795] [23796]. Pages: "Overdrive" x3 (pages/export_overdrive.py).
"""
from __future__ import annotations

import itertools as it
import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model import kinematics as MJ
from analysis.model.first_principles import geometry as G

PATCHES = {
    "ring": [(0, 0, 0), (1, 1, 1), (2, 2, 0), (1, 1, -1)],
    "hc15": [tuple(int(c) for c in s) for s in MJ.hc15_sites()],
    "block": [tuple(int(c) for c in s) for s in RC.brick(5, 5, 5)],
}
ANGLES = [float(a) for a in range(-60, 301, 30)]


def drive(sites):
    """The patch at every angle, welds held by identity from -30."""
    ref, _ = RC.honeycomb(sites, gc=-30.0)
    welds = ref.welds
    kind = np.array([all(c % 2 == 0 for c in s) for s in sites])      # True = cell
    asm0, _ = RC.honeycomb(sites, gc=0.0)
    S = G.front_signs(asm0.positions(asm0.q0()), asm0.ctr0)
    rows = {}
    for a in ANGLES:
        asm, _ = RC.honeycomb(sites, gc=a)
        held = RC.Assembly(asm.gam0, asm.ctr0, welds)
        q = held.q0()
        X = held.positions(q)
        res = float(np.abs(held.weld_residual(q)).max())
        out = G.fronts_out(X, asm.ctr0, S)
        P, Q, owner = G.assembly_struts(X)
        cross = G.crossings(P, Q, owner)
        spread = float(np.max(np.linalg.norm(asm.ctr0 - asm.ctr0.mean(0), axis=1)))
        rows[a] = dict(res=res, L=RC.lattice_constant(a), cells_out=int(out[kind].sum()),
                       voids_out=int(out[~kind].sum()), cross=cross, spread=spread, X=X)
    return rows, int(kind.sum()) * 8, int((~kind).sum()) * 8, len(welds)


def gate():
    checks, out = [], {}
    A = checks.append
    data = {name: drive(sites) for name, sites in PATCHES.items()}
    out["data"] = data

    # ---- R1: welds by identity ---------------------------------------------------
    worst = {n: max(r["res"] for r in rows.values()) for n, (rows, *_ ) in data.items()}
    A(("R1  THE JOINTS NEVER NEED THE MAGIC: welds read once at -30 and held by identity "
       "stay satisfied to 1e-15 through all 360 degrees on every patch. Bodies pass through "
       "each other and through their own centres; no joint parts.",
       all(v < 1e-12 for v in worst.values()), {n: f"{v:.0e}" for n, v in worst.items()}, "all < 1e-12"))

    # ---- R2: the spacing law -------------------------------------------------------
    err = max(abs(RC.lattice_constant(a) - (2 / np.sqrt(3)) * np.cos(np.radians(a + 30.0))) for a in ANGLES)
    L = {a: RC.lattice_constant(a) for a in (-30.0, 0.0, 60.0, 150.0, 240.0)}
    A(("R2  THE SPACING IS A PURE SINUSOID, L(a) = (2/sqrt3) cos(a + 30): 1.155 at -30, 1 at "
       "the two dead ends, ZERO at +60 and +240, -1.155 at +150. The physical range is the "
       "sixty degrees near the top.",
       err < 1e-9 and abs(L[60.0]) < 1e-12 and abs(L[240.0]) < 1e-12 and L[150.0] < -1.15,
       {a: round(v, 4) for a, v in L.items()} | {"law error": f"{err:.0e}"}, "1.1547, 1, 0, -1.1547, 0; < 1e-9"))

    # ---- R3: sides -------------------------------------------------------------------
    ok3, got3 = True, {}
    for n, (rows, nc, nv, _) in data.items():
        for a, r in rows.items():
            cc, cv = np.cos(np.radians(a)), np.cos(np.radians(a + 60.0))
            if abs(cc) > 1e-9:
                ok3 &= (r["cells_out"] == (nc if cc > 0 else 0))
            if abs(cv) > 1e-9:
                ok3 &= (r["voids_out"] == (nv if cv > 0 else 0))
        got3[n] = " ".join(f"{int(a):+d}:{r['cells_out']}/{r['voids_out']}" for a, r in rows.items())
    A(("R3  SIDES FLIP IN TWO WAVES: every plate's front points out of its body while the "
       "plate is on its body's near side of centre (cos fold > 0) and in once it has passed "
       "through -- voids around +30 and back at +210, cells around +90 and back at +270; "
       "from 120 to 180 every plate in the array is inside out.",
       ok3, got3, "cells out iff cos a > 0; voids iff cos(a + 60) > 0"))

    # ---- R4: crossings confined to the passages ---------------------------------------
    ok4, got4 = True, {}
    for n, (rows, *_ ) in data.items():
        for a, r in rows.items():
            ok4 &= (r["cross"] > 0) if (a in (30.0, 90.0, 210.0, 270.0)) else (r["cross"] == 0)
        got4[n] = {int(a): r["cross"] for a, r in rows.items()}
    A(("R4  STRUT CROSSINGS ARE CONFINED TO THE PASSAGES through the collapse (30, 90, "
       "210, 270, where |L| = 0.577): none at any multiple of 60 and none at -30 or 150 "
       "where the spacing is widest; hundreds mid-passage, the count scaling with the patch.",
       ok4, got4, "> 0 only at 30, 90, 210, 270"))

    # ---- R5: the collapse, and the period ----------------------------------------------
    ok5, got5 = True, {}
    for n, (rows, *_ ) in data.items():
        sp = {int(a): round(r["spread"], 4) for a, r in rows.items()}
        period = float(np.abs(rows[300.0]["X"] - rows[-60.0]["X"]).max())
        ok5 &= sp[60] < 1e-9 and sp[240] < 1e-9 and sp[0] > 1.0 and period < 1e-9
        got5[n] = {"spread": sp, "|X(300) - X(-60)|": f"{period:.0e}"}
    A(("R5  AT +60 AND +240 EVERY CENTRE IS ON ONE POINT -- the whole patch is one "
       "coincident octahedron -- and 300 is -60 again, position for position. None of it "
       "depends on the size of the patch.",
       ok5, got5, "spread 0 at 60 and 240; period exact"))
    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("overdrive -- the ring and larger patches driven through 360 degrees")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")
        print("\n  THE THREE PATCHES")
        print(f"   {'patch':6s} {'bodies':>6s} {'welds':>5s}   " + " ".join(f"{int(a):>5d}" for a in ANGLES))
        for n, (rows, nc, nv, nw) in out["data"].items():
            print(f"   {n:6s} {nc // 8 + nv // 8:6d} {nw:5d}   crossings " + " ".join(f"{rows[a]['cross']:5d}" for a in ANGLES))
            print(f"   {'':6s} {'':6s} {'':5s}   spread    " + " ".join(f"{rows[a]['spread']:5.2f}" for a in ANGLES))
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * NOTHING HERE IS THE MEDIUM. It is the model's parametrisation continued")
        print("     analytically through configurations rigid struts cannot occupy.")
        print("   * WHAT IT SHOWS: the joints' kinematics never object; the passage does.")
        print("     An infinite array driven this way is a homothety to a point and back,")
        print("     mirrored each time, period 360.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
