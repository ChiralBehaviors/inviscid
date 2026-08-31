"""plate_holonomy -- G1: the quaternion holonomy of every plate over one
full 360-degree overdrive cycle, on the record's three patches. A
CONFIRMATION gate: the record already predicts -q.

THE QUESTION (bead inviscid-6h9; note sections 8.3, 8.4, 8.6). The record
holds that a plate's spin about its own axis is the fold angle (one_cell
R2), that body frames sit at identity through the coherent drive
(assembly's closed form), and that the drive returns position-for-position
after 360 degrees (overdrive R5). One cycle is therefore 2 pi of plate
spin, and the continuous SU(2) lift of any plate's orientation path must
return with the OPPOSITE sign: q(300) = -q(-60). Section 8.3 calls this a
near-corollary of gated facts; this module makes it a measurement instead
of an inference -- fine-sampled, every plate, both body kinds, three patch
sizes. A +q anywhere means the lift code or the record is wrong, and is to
be investigated, not published.

WHAT IS MEASURED HERE. The drive is `assembly.coherent` -- the medium's
one motion in closed form (gamma_even = a, gamma_odd = a + 60, R = I,
centre = L(a) * site), gated against the constrained model by assembly
R5h -- swept a = -60 .. 300 on a grid of 0.5 degrees, refined to 0.25
through the inside-out region [90, 270] and to 0.05 within 5 degrees of
both collapses (+60, +240). At every frame each plate's rotation relative
to its start pose is fit from its three matched corners (batched Kabsch,
cross-checked against geometry.kabsch); each plate's rotation path is then
lifted to SU(2) by the G0-calibrated instrument (model.su2.lift) with the
step guard tightened to 1 degree.

  R1  The drive here IS the record's drive, and it closes: the closed-form
      positions match the overdrive's held-weld construction at all 13
      gated angles; the batched rotation fit matches geometry.kabsch
      there; every plate's R(300) equals its R(-60).
  R2  The sampling obeys the gate plan mechanically: steps <= 0.5 deg
      everywhere, <= 0.25 through [90, 270], <= 0.05 through the
      collapses; the lift's guard is 1 deg and the worst step it measures
      is 0.5 deg.
  R3  The whole cycle is 2 pi of pure spin about each plate's OWN axis:
      every plate's fit rotation equals rot(u_f, sigma_f * (a + 60)) on
      every frame -- one_cell R2's spin law extended through the full
      cycle, voids running 60 degrees ahead on the same law.
  R4  HOLONOMY: every plate of every body on every patch lifts to -q --
      ring 32/32, hc15 120/120, block 280/280, cells and voids alike.
      q(300) = -q(-60): the record's prediction, confirmed.

WHAT THIS DOES NOT SAY. -q is necessary, NOT sufficient, for any
"SU(2) medium" claim (section 8.4): it is classical Z2 holonomy of a path
in SO(3), the point joints as modeled store nothing (section 8.2), and the
PHYSICAL sixty-degree arc never traverses this loop -- the overdrive
passages are strut-crossing territory (overdrive R4, section 8.5). Whether
the physical arc exercises holonomy at all is G7 (inviscid-1lc).

T2: [23865]. Ref: su2_boundary_conditions.md section 8 (read before 1-7).
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model import plates as Z
from analysis.model.first_principles import geometry as G
from analysis.model.first_principles.overdrive import ANGLES, PATCHES
from analysis.model.su2 import lift as SU

FACES = Z.faces()
A_START, A_END = -60.0, 300.0


def grid():
    """The G1 sampling grid: 0.5 base, 0.25 through [90, 270], 0.05 at the
    collapses -- the bead's sampling requirement, constructed not promised."""
    segs = [np.arange(A_START, A_END + 1e-9, 0.5),
            np.arange(90.0, 270.0 + 1e-9, 0.25)]
    for c in (60.0, 240.0):
        segs.append(np.arange(c - 5.0, c + 5.0 + 1e-9, 0.05))
    return np.unique(np.round(np.concatenate(segs), 6))


def coherent_X(sites, a):
    """The coherent drive's vertex positions at fold angle a: identity body
    frames, centres L(a) * site, corners body(gamma) -- assembly.coherent
    unpacked, with body() evaluated once per distinct gamma."""
    q, _ = RC.coherent(sites, a)
    ctr, _quat, gam = RC.Assembly.unpack(q)
    V = {g: RC.body(g) for g in np.unique(gam)}
    return np.array([ctr[k] + V[gam[k]] for k in range(len(sites))])


def plate_corners(X):
    """(bodies*8, 3, 3): each plate's three matched corners, plate-major."""
    return np.array([X[k][G.TRI[f]] for k in range(X.shape[0]) for f in range(8)])


def kabsch_batch(P0, P):
    """Batched proper-rotation fit R[f, p] with R @ (P0[p]-c0) = P[f, p]-c,
    the same algebra as geometry.kabsch, over frames f and plates p."""
    A = P0 - P0.mean(1, keepdims=True)                      # (p, 3, 3)
    B = P - P.mean(2, keepdims=True)                        # (F, p, 3, 3)
    H = np.einsum("pki,fpkj->fpij", A, B)
    U, _, Vt = np.linalg.svd(H)
    VtT, UT = np.swapaxes(Vt, -1, -2).copy(), np.swapaxes(U, -1, -2)
    VtT[..., :, 2] *= np.sign(np.linalg.det(VtT @ UT))[..., None]
    return VtT @ UT


def drive_rotations(sites, degs):
    """(plates, frames, 3, 3): each plate's rotation relative to A_START."""
    C0 = plate_corners(coherent_X(sites, A_START))
    R = np.empty((len(C0), len(degs), 3, 3))
    for i, a in enumerate(degs):
        R[:, i] = kabsch_batch(C0, plate_corners(coherent_X(sites, a))[None, ...])[0]
    return R


def gate():
    checks = []
    A = checks.append
    degs = grid()

    # ---- R1: this drive is the record's drive, and it closes ------------------
    worst_pos, worst_fit = 0.0, 0.0
    for sites in PATCHES.values():
        C0 = plate_corners(coherent_X(sites, A_START))
        for a in ANGLES:
            asm, _ = RC.honeycomb(sites, gc=a)
            ref, _ = RC.honeycomb(sites, gc=-30.0)
            held = RC.Assembly(asm.gam0, asm.ctr0, ref.welds)
            Xh = held.positions(held.q0())
            Xc = coherent_X(sites, a)
            worst_pos = max(worst_pos, float(np.abs(Xh - Xc).max()))
            Cp = plate_corners(Xc)
            Rb = kabsch_batch(C0, Cp[None, ...])[0]
            for p in (0, len(C0) // 2, len(C0) - 1):
                Rs, _t = G.kabsch(C0[p], Cp[p])
                worst_fit = max(worst_fit, float(np.abs(Rb[p] - Rs).max()))
    Rpaths = {n: drive_rotations(s, degs) for n, s in PATCHES.items()}
    worst_close = max(float(np.abs(R[:, -1] - R[:, 0]).max()) for R in Rpaths.values())
    A(("R1  THE DRIVE HERE IS THE RECORD'S DRIVE, AND IT CLOSES: the closed-form coherent "
       "positions match the overdrive's held-weld construction at all 13 gated angles on "
       "all three patches; the batched rotation fit matches geometry.kabsch there; and "
       "every plate's R(300) equals its R(-60).",
       worst_pos < 1e-12 and worst_fit < 1e-12 and worst_close < 1e-9,
       f"positions {worst_pos:.0e}, fit vs kabsch {worst_fit:.0e}, closure {worst_close:.0e}",
       "< 1e-12, < 1e-12, < 1e-9"))

    # ---- R2: the sampling obeys the gate plan mechanically --------------------
    d = np.diff(degs)
    mid = (degs[:-1] + degs[1:]) / 2
    inside = d[(mid > 90) & (mid < 270)]
    coll = d[(np.abs(mid - 60) < 5) | (np.abs(mid - 240) < 5)]
    _q, worst_step = SU.lift(Rpaths["ring"][0], max_step_deg=1.0)
    A(("R2  THE SAMPLING OBEYS THE GATE PLAN MECHANICALLY: steps at most 0.5 deg "
       "everywhere, 0.25 through the inside-out region [90, 270], 0.05 within 5 deg of "
       "both collapses; the lift guard is tightened to 1 deg and the worst step it "
       "measures is the base step.",
       d.max() <= 0.5 + 1e-9 and inside.max() <= 0.25 + 1e-9 and coll.max() <= 0.05 + 1e-9
       and abs(worst_step - 0.5) < 1e-6,
       f"{len(degs)} frames; max step {d.max():.2f}, inside-out {inside.max():.2f}, "
       f"collapses {coll.max():.3f}; worst lift step {worst_step:.4f} deg",
       "0.50, 0.25, 0.050; 0.5"))

    # ---- R3: the cycle is 2 pi of pure own-axis spin --------------------------
    ana = np.array([[SU.rot(FACES[f][2], FACES[f][3] * (a - A_START)) for a in degs]
                    for f in range(8)])
    worst_spin = {}
    for n, R in Rpaths.items():
        dev = 0.0
        for p in range(R.shape[0]):
            dev = max(dev, float(np.abs(R[p] - ana[p % 8]).max()))
        worst_spin[n] = dev
    A(("R3  THE WHOLE CYCLE IS 2 PI OF PURE SPIN ABOUT EACH PLATE'S OWN AXIS: every "
       "plate's fit rotation equals rot(u_f, sigma_f * (a + 60)) on every frame of every "
       "patch -- one_cell R2's spin law extended through the full cycle, the voids "
       "running sixty degrees ahead on the same law, no plate ever borrowing any other "
       "axis. The net turn is 360 degrees by construction of the match.",
       all(v < 1e-9 for v in worst_spin.values()),
       {n: f"{v:.0e}" for n, v in worst_spin.items()}, "all < 1e-9"))

    # ---- R4: holonomy -q everywhere -------------------------------------------
    got4, ok4 = {}, True
    for n, R in Rpaths.items():
        sites = PATCHES[n]
        kind = np.repeat([all(c % 2 == 0 for c in s) for s in sites], 8)
        h = np.array([SU.holonomy(R[p], max_step_deg=1.0) for p in range(R.shape[0])])
        ok4 &= bool((h == -1).all())
        got4[n] = (f"{int((h == -1).sum())}/{len(h)} -q "
                   f"({int((h[kind] == -1).sum())} cell, {int((h[~kind] == -1).sum())} void)")
    A(("R4  HOLONOMY: EVERY PLATE OF EVERY BODY ON EVERY PATCH LIFTS TO -q. One full "
       "overdrive cycle carries each plate around the nontrivial loop of SO(3): "
       "q(300) = -q(-60), cells and voids alike, invariant in the patch size. The "
       "record's prediction (section 8.3), confirmed as a measurement.",
       ok4, got4, "ring 32/32, hc15 120/120, block 280/280, all kinds"))
    return checks


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("plate_holonomy -- G1: every plate lifts to -q over one 360-degree "
          "overdrive cycle")
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
    print("   * THE KINEMATIC FACT IS NOW MEASURED, NOT INFERRED: one overdrive cycle")
    print("     is the nontrivial loop of SO(3) for every plate, at fine sampling,")
    print("     on three patch sizes. A CONFIRMATION of the record (section 8.3).")
    print("   * -q IS NECESSARY, NOT SUFFICIENT, for any SU(2)-medium claim (section")
    print("     8.4). It is classical Z2 holonomy of a parametrised path; the point")
    print("     joints as modeled store nothing (section 8.2), and no closure choice")
    print("     for finite simulations follows from it -- that is the G5/G6 lane.")
    print("   * THE PHYSICAL ARC IS A DIFFERENT QUESTION: the sixty degrees that can")
    print("     be built never traverses this loop (crossings forbid the passages,")
    print("     overdrive R4). Whether it exercises holonomy at all is G7.")
    print()
    print("  ALL CHECKS PASSED." if not bad
          else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
