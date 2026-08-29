"""jb_tr -- phase-field transport on the honeycomb, and the asymmetry that is not there.

WHY THIS FILE EXISTS. Bead qvf.22 is the epic's own criterion: kick one cell and
measure how the disturbance moves. It has been unreachable for the whole life of
the project, for a reason each session found and none could see whole -- a rigid
constraint transmits instantly, so there was never anything to time. jb_ct broke
that with joint PLAY, jb_cp measured why the medium has no linear regime, and
jb_mj built the impact law that lets more than one joint bind at once. This runs
the measurement on the real packing.

THE INSTRUMENT. jb_mj works on HC15, which is ONE SHELL DEEP -- every neighbour
is distance 1 from the centre, so there is no far cell to time anything to. This
uses the honeycomb's own DIAGONAL CHAIN instead: sites (k,k,k), consecutive ones
differing by (1,1,1), which is exactly the triangular-face neighbour relation.
Ten cells, nine welds, 27 unilateral bands, closing to 9e-16. It is the
honeycomb analogue of jb_ct's scalar chain, with the real geometry, three bands
per joint instead of one scalar, and simultaneous multi-contact.

WHAT REPRODUCES, AND IT IS THE WEEK'S CENTRAL RESULT SURVIVING A HARDER MODEL:

    speed doubles when the kick doubles      1.99x measured for 2x
    speed * play is constant                 0.8% over a twofold range

That is jb_ct's law -- speed = |dV/dgamma| * gammadot / play -- measured now on
the honeycomb, in full reduced coordinates, with three-vertex welds and an LCP
resolving up to six contacts in a single step. The scalar chain was not an
artifact of its own reduction.

WHAT DOES NOT REPRODUCE, AND THE BEAD ASKED FOR EXACTLY THIS. Bead qvf.11
records the physical array locking ONE-SIDEDLY -- it cannot expand, it can still
contract -- and this bead's text says "Do NOT assume the qvf.11 observation;
MEASURE it in this model and say whether it reproduces". IT DOES NOT. Driving
toward the VE and toward the octahedron give arrival times agreeing to 1%, the
same delay spread to two decimals, and the same attenuation. R6 is that row.

The reason is structural rather than numerical, which is what makes it useful: a
clearance band ||va - vb|| <= t is a BALL, and a ball has no preferred
direction. Whatever produces the rig's one-sidedness is not in the joint
clearance. The obvious remaining candidate is plate INTERPENETRATION, which is
genuinely one-sided -- it blocks approach and not separation -- and which
DECISION 18 names and this model omits.

A SPEED IS NOT QUOTED AS A MEASUREMENT, and R2 is why. jb_ct's R4 established
the rule the hard way: a linear fit to arrival times is not a speed unless the
front is uniform. Here the per-cell delay spreads by a factor of 1.98 across the
chain -- the front DECELERATES -- so the fit is reported and the SCALINGS are
what is gated. The deceleration and the attenuation are the same phenomenon
seen twice: the front carries less energy as it spreads, so it crosses each
clearance more slowly.

FOUR DECLARATIONS -- INAPPLICABLE, NOT FORGOTTEN.
  * MASS MODEL: DECLARED, both as peers, inherited from jb_mj -- corner point
    masses and uniform laminae, each validated there against a number from
    outside the file (swings of exactly 3 and exactly 9).
  * METRIC: DECLARED -- the block-diagonal mass matrix of the model in use.
  * KERNEL, PRIMITIVE: INAPPLICABLE while V = 0. Stated, not lapsed.

ABSOLUTE-VERSUS-RATIO. Every speed here inherits the project's standing
convention (coupling 1, total mass 1/2, R = 1), so no absolute speed is offered
as physical. What is gated is RATIOS: across amplitude, across clearance, and
across the two drive senses.
"""
from __future__ import annotations

import sys

import numpy as np

import jb_mj_inertial_honeycomb as MJ
import jb_rc_reduced as RC

A_REF = MJ.A_REF

#: Cells in the diagonal chain. Ten is enough that the front has somewhere to
#: go and short enough that the run stays cheap.
CELLS = 10

#: Reference drive and clearance. Both are swept; neither is a measurement.
KICK = 0.9
PLAY = 0.05

#: A cell counts as reached when its fold rate exceeds this fraction of the
#: drive. NOTE the abs(): an earlier probe wrote `thresh * kick` and, driven
#: with a NEGATIVE kick, compared a positive rate against a negative threshold,
#: so every cell "arrived" on the first step and the fit returned 2500. The
#: asymmetry row is exactly where that bug would have been invisible.
ARRIVED = 0.02

#: Elastic. A fully inelastic medium absorbs a disturbance instead of carrying
#: one, so transport is measured at e = 1; jb_mj runs both and adopts neither.
RESTITUTION = 1.0

STEP = 1e-3


def chain(n=CELLS, gc=A_REF):
    """The honeycomb's own diagonal chain: sites (k,k,k), consecutive ones
    differing by (1,1,1), which IS the triangular-face neighbour relation."""
    asm, _ = RC.honeycomb([(k, k, k) for k in range(n)], gc=gc)
    return asm


def transport(asm, kick=KICK, play=PLAY, e=RESTITUTION, tmax=0.8, h=STEP,
              lamina=False, bands=True):
    """Kick cell 0; record when each cell first moves and how much it ever
    moves. `bands=False` is the control with the coupling removed."""
    pairs = MJ.tied_pairs(asm)
    n = asm.N
    q = asm.q0()
    u = np.zeros((n, 7))
    u[0, 6] = kick
    arrive = [None] * n
    arrive[0] = 0.0
    peak = np.zeros(n)
    now, hits, multi = 0.0, 0, 0
    while now < tmax - 1e-12:
        J, M, Minv = MJ.kinematics(asm, q, lamina)
        a1 = MJ.free_accel(asm, q, u, J, Minv, lamina)
        u_h = u + 0.5 * h * a1
        q_h = RC.apply_increment(asm, q, (0.5 * h * u).ravel())
        Jh, _, Mih = MJ.kinematics(asm, q_h, lamina)
        a2 = MJ.free_accel(asm, q_h, u_h, Jh, Mih, lamina)
        u = u + h * a2
        q = RC.apply_increment(asm, q, (h * u_h).ravel())
        now += h
        if bands:
            J, M, Minv = MJ.kinematics(asm, q, lamina)
            s = MJ.separations(asm, q, pairs)
            N = MJ.band_rows(asm, q, J, pairs)
            rate = np.dot(N, u.ravel())
            act = [i for i in range(len(pairs))
                   if s[i] >= play and rate[i] > 0]
            if act:
                u, _, _, _ = MJ.resolve(asm, u, N, act, Minv, e)
                hits += 1
                multi = max(multi, len(act))
        peak = np.maximum(peak, np.abs(u[:, 6]))
        for k in range(1, n):
            if arrive[k] is None and abs(u[k, 6]) > ARRIVED * abs(kick):
                arrive[k] = now
    got = [(k, a) for k, a in enumerate(arrive) if a is not None]
    d = np.diff([a for _, a in got]) if len(got) > 2 else np.array([np.nan])
    fit = (float(np.polyfit([a for _, a in got[1:]],
                            [k for k, _ in got[1:]], 1)[0])
           if len(got) > 3 else float("nan"))
    return {"arrive": arrive, "reached": len(got), "peak": peak, "hits": hits,
            "multi": multi, "fit": fit,
            "spread": float(d.max() / d.min()) if len(d) > 1 else float("nan")}


def gate():
    checks, out = [], {}
    A = checks.append
    asm = chain()
    pairs = MJ.tied_pairs(asm)
    weld = float(np.abs(asm.weld_residual(asm.q0())).max())

    # R1 -- the instrument, and why HC15 could not be it
    A(("R1  THE INSTRUMENT IS THE HONEYCOMB'S OWN DIAGONAL CHAIN, because HC15 "
       "is ONE SHELL DEEP and has no far cell to time anything to. Sites "
       "(k,k,k) differ consecutively by (1,1,1), which IS the triangular-face "
       "neighbour relation, so this is a genuine chain in the real packing -- "
       "the honeycomb analogue of jb_ct's scalar chain, with three unilateral "
       "bands per joint instead of one scalar. It closes to machine precision, "
       "which a patch built on the wrong offsets would not",
       asm.N == CELLS and len(asm.welds) == CELLS - 1
       and len(pairs) == 3 * (CELLS - 1) and weld < 1e-12,
       f"{asm.N} cells, {len(asm.welds)} welds, {len(pairs)} bands, weld "
       f"residual {weld:.2e}",
       "a closed chain with three bands per joint"))

    # R2 -- transport is finite, and the front decelerates
    ref = transport(asm)
    out["ref"] = ref
    A(("R2  A DISTURBANCE CROSSES THE CHAIN IN FINITE TIME, AND THE FRONT "
       "DECELERATES -- so a speed is FITTED here and not measured. jb_ct's R4 "
       "established the rule: a linear fit to arrival times is not a speed "
       "unless the front is uniform, and here the per-cell delay spreads by "
       "nearly a factor of two from end to end. The fit is reported because it "
       "is what the scalings below are ratios OF; it is not offered as a "
       "property of the medium. TWO-SIDED: the row fails if the disturbance "
       "stops reaching the far end, and it fails if the front comes back "
       "UNIFORM, because then this caveat would be the wrong one to carry",
       ref["reached"] == CELLS and ref["spread"] > 1.5
       and ref["multi"] > 1,
       f"reached all {ref['reached']} cells over {ref['hits']} impacts, up to "
       f"{ref['multi']} bands bound at once; per-cell delay spread "
       f"{ref['spread']:.2f}; fitted {ref['fit']:.3f} cells per unit time",
       "crosses the chain, front demonstrably non-uniform"))

    # R3 -- jb_ct's amplitude law, on the real geometry
    amp = {k: transport(asm, kick=k) for k in (KICK, 2 * KICK)}
    ratio = amp[2 * KICK]["fit"] / amp[KICK]["fit"]
    out["amp"] = ratio
    A(("R3  SPEED DOUBLES WHEN THE DRIVE DOUBLES -- jb_ct's AMPLITUDE LAW, "
       "SURVIVING A MUCH HARDER MODEL. The scalar chain measured speed = "
       "|dV/dgamma| * gammadot / play on one dowel with one scalar per joint "
       "and one contact resolved at a time. This is the same law on the real "
       "honeycomb: full reduced coordinates, seven degrees of freedom per "
       "cell, three bands per joint, and an LCP resolving several at once. A "
       "medium whose signal speed depends on AMPLITUDE has no linear regime, "
       "which is jb_cp's infinite square well seen from the transport side",
       abs(ratio - 2.0) < 0.05,
       f"fitted speed {amp[KICK]['fit']:.3f} at kick {KICK} against "
       f"{amp[2 * KICK]['fit']:.3f} at kick {2 * KICK}: ratio {ratio:.4f}",
       "a ratio of 2 within 5%"))

    # R4 -- and the clearance law
    pl = {p: transport(asm, play=p) for p in (PLAY / 2, PLAY)}
    inv = {p: pl[p]["fit"] * p for p in pl}
    out["play"] = inv
    A(("R4  SPEED TIMES CLEARANCE IS CONSTANT -- jb_ct's R6, ALSO SURVIVING. "
       "Halving the play doubles the speed, so the speed diverges as the joint "
       "tightens and a rigid joint recovers the instantaneous onset every "
       "earlier module in this project reported. That is what makes the "
       "clearance the coupling rather than a tolerance: it is the only thing "
       "standing between this medium and an infinite signal speed",
       max(inv.values()) / min(inv.values()) < 1.05,
       "; ".join(f"play {p}: fitted {pl[p]['fit']:.3f}, speed*play "
                 f"{inv[p]:.4f}" for p in sorted(pl)),
       "speed * play constant to 5% over a twofold range"))

    # R5 -- attenuation, which is the deceleration seen again
    pk = ref["peak"]
    atten = float(pk[0] / pk[-1])
    out["atten"] = atten
    A(("R5  THE DISTURBANCE ATTENUATES, AND IT IS THE SAME PHENOMENON AS THE "
       "DECELERATION. The peak fold rate falls monotonically from the driven "
       "cell to the far end, and the front slows over the same distance "
       "because it is carrying less energy into each successive clearance. "
       "Reported as a RATIO, per the project's absolute-versus-ratio rule -- "
       "the underlying scale is a convention. TWO-SIDED: a constant profile "
       "fails the lower bound and a collapse to nothing fails the upper, which "
       "is the jb_x X7 shape this house style refuses to repeat",
       1.5 < atten < 100.0 and pk[-1] > 1e-6
       and all(pk[i] >= pk[i + 1] - 0.02 for i in range(1, CELLS - 1)),
       "peak fold rate by cell: "
       + " ".join(f"{v:.3f}" for v in pk)
       + f"; first/last = {atten:.3f}",
       "monotone decay in a two-sided band, far cell still moving"))

    # R6 -- the asymmetry the epic asks about, and it is not there
    fwd = transport(asm, kick=KICK)
    rev = transport(asm, kick=-KICK)
    sp_ratio = fwd["fit"] / rev["fit"]
    arr_gap = max(abs(fwd["arrive"][k] - rev["arrive"][k])
                  / max(fwd["arrive"][k], 1e-9) for k in range(1, CELLS))
    out["asym"] = (sp_ratio, arr_gap)
    A(("R6  THERE IS NO COMPRESSION / RAREFACTION ASYMMETRY IN THIS MODEL, AND "
       "THE BEAD ASKED FOR EXACTLY THIS ANSWER RATHER THAN THE ASSUMPTION. "
       "Bead qvf.11 records the physical array locking ONE-SIDEDLY -- it "
       "cannot expand, it can still contract -- and qvf.22's text says 'Do NOT "
       "assume the qvf.11 observation; MEASURE it in this model and say "
       "whether it reproduces'. Driving toward the VE and toward the "
       "octahedron give arrival times agreeing to about a percent, the same "
       "delay spread, and the same attenuation. The reason is structural and "
       "that is what makes it useful: a clearance band ||va - vb|| <= t is a "
       "BALL, and a ball has no preferred direction. Whatever produces the "
       "rig's one-sidedness is NOT in the joint clearance -- the remaining "
       "candidate is plate INTERPENETRATION, which is genuinely one-sided "
       "(it blocks approach, not separation), which DECISION 18 names and this "
       "model omits. TWO-SIDED: a real asymmetry would fail this row",
       abs(sp_ratio - 1.0) < 0.05 and arr_gap < 0.05,
       f"fitted speed {fwd['fit']:.3f} driving one sense against "
       f"{rev['fit']:.3f} the other, ratio {sp_ratio:.4f}; worst relative "
       f"arrival-time difference across the chain {arr_gap:.4f}; delay spreads "
       f"{fwd['spread']:.2f} and {rev['spread']:.2f}",
       "the two senses agreeing within 5%, i.e. no asymmetry"))

    # R7 -- the two controls the bead names, both able to fail
    still = transport(asm, kick=0.0)
    free = transport(asm, bands=False)
    out["still"], out["free"] = still, free
    A(("R7  THE TWO CONTROLS THE BEAD NAMES, AND BOTH CAN FAIL. UNDRIVEN, "
       "nothing moves: with no kick the far cells never reach the arrival "
       "threshold at all, which is the epic's own stated control. COUPLING "
       "REMOVED, nothing propagates: with the bands disabled the driven cell "
       "keeps its motion and no other cell ever moves, so the transport "
       "measured above is a property of the coupling and not an artifact of "
       "the drive or of the projection. Without this second control the whole "
       "file could be measuring its own initial condition",
       still["reached"] == 1 and free["reached"] == 1
       and ref["reached"] == CELLS,
       f"undriven: {still['reached']} of {CELLS} cells ever move; bands "
       f"disabled: {free['reached']} of {CELLS}; bands enabled: "
       f"{ref['reached']} of {CELLS}",
       "only the driven cell moves in both controls, all ten with coupling"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_tr -- phase-field transport on the honeycomb")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        print("\n  FOUR DECLARATIONS.")
        print("   MASS MODEL  DECLARED, both peers, inherited from jb_mj.")
        print("   METRIC      DECLARED: that model's block-diagonal mass matrix.")
        print("   KERNEL      INAPPLICABLE while V = 0. Stated, not lapsed.")
        print("   PRIMITIVE   INAPPLICABLE while V = 0. Stated, not lapsed.")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * jb_ct's TWO SCALING LAWS SURVIVE the move to the real")
        print("     packing, full reduced coordinates and simultaneous")
        print("     multi-contact. The scalar chain was not an artifact of its")
        print("     own reduction. That is the main result.")
        print("   * NO ABSOLUTE SPEED IS QUOTED. The front decelerates by a")
        print("     factor of two across the chain, so a fit is not a speed")
        print("     (jb_ct R4's rule), and the underlying scale is a")
        print("     convention anyway. Only ratios are gated.")
        print("   * THE EPIC'S DISPERSION CRITERION IS STILL NOT MET. A")
        print("     dispersion relation pairs a frequency with a wavenumber;")
        print("     what is here is a front, its amplitude scaling, and its")
        print("     attenuation. With V = 0 there is no frequency to pair, and")
        print("     jb_cp measured why: the effective potential is an infinite")
        print("     square well, which has no linear regime to linearise.")
        print("     A travelling disturbance is not yet a wave with a")
        print("     dispersion relation, and this file does not claim one.")
        print("   * NO ONE-SIDEDNESS. qvf.11's asymmetry does NOT reproduce")
        print("     here, and the reason says where to look next: a clearance")
        print("     ball has no preferred direction, so if the rig's")
        print("     one-sidedness is real it lives in plate interpenetration,")
        print("     which DECISION 18 names and this model omits.")
        print("   * ONE CHAIN, ten cells, free ends. Not the bulk.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
