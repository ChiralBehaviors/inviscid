"""jb_1v -- one covering, the sonic vacuum: the global laws and the solitary
front survive, the driver hands off clean, and jb_sv's local law does not
hold because the fold's decay is recoil, not attenuation.

WHY THIS FILE EXISTS. DECISION 19 (T2 [23673]) reformulated the epic's
criterion into four sonic-vacuum parts and jb_sv met them on a hard-wall
chain: speed ~ A^1, speed ~ 1/play, an attenuation, and a front that does
not spread. Every one of those was measured on `jb_rc.honeycomb`'s body-
diagonal chain, the double covering. DECISION 21 (T2 [23727]) rules the
voids empty. jb_1e re-measured the SMOOTH joint under one covering and the
exponents survived; this file re-measures the HARD WALL -- jb_mj's LCP
pointed at the single covering for the first time. `jb_sv.front_history`
takes an assembly and `MJ.tied_pairs`, so the instrument is unchanged: only
the chain is new, jb_1e's axis chain, sixteen cells, fifteen two-pair welds,
thirty unilateral bands against jb_sv's forty-five.

WHAT SURVIVES.
  * THE GLOBAL LAWS. Speed against kick, exponent 0.994 over a fourfold
    drive; speed against clearance, exponent -0.994 over fourfold; and the
    two collapse onto one number, speed * play / kick = 0.995, spread 0.9%
    across every run (R3). jb_ct's law, jb_tr's R3/R4, on a chain whose
    joint is a different object.
  * THE SOLITARY FRONT. One cell wide at every cell it reaches, from cell 1
    to cell 15 (jb_sv's is two cells wide at cell 1, then one). The
    reformulated criterion's part 4, sharper.

WHAT CHANGES, AND IT IS NOT SUBTLE.
  * THE DRIVER HANDS OFF CLEAN. jb_sv R4 gates, as an honest limitation,
    that its driven cell KEEPS 0.3593 of its 0.9 when the front reaches the
    far end and sheds a train behind the front. Here the driven cell keeps
    0.0009. The front cell carries 99.9% of the chain's kinetic energy at
    cell 1 and 98.9% at cell 15; the wake holds 1.1% after fifteen
    handoffs. That caveat does not carry to the single covering.
  * jb_sv's LOCAL LAW DOES NOT HOLD. jb_sv R2/R3 measure that WITHIN one
    run the front's speed falls by the factor its amplitude falls by
    (1.8767 against 1.9238, p = 0.9880), and read from that that the
    deceleration and the attenuation are one phenomenon. On the single
    chain the transit time per cell is CONSTANT to the integration step --
    0.056 at every one of fifteen cells at kick 0.9 -- while the fold
    amplitude decays from 0.899 to 0.820. Fitted the way jb_sv fits it,
    p_local = 0.00.
  * WHERE THE FOLD GOES. Not into the wake. The decayed fold energy is in
    the front cell's OWN axial translation: its velocity along the chain
    grows monotonically from -0.027 at cell 1 to -0.359 at cell 15, with
    zero transverse velocity and zero angular velocity to 1e-12, and
    KE_f - KE_fold grows 0.003 -> 0.515 while KE_f stays within 1% of the
    total. Each handoff converts a little fold into recoil, inside the
    cell that carries the front; the front's energy and speed are
    untouched. Linear momentum along the chain is conserved to 1e-12
    throughout (R4), so the recoil is balanced by the wake's slow drift the
    other way.

So on the single chain there is neither deceleration nor attenuation to
speak of, and the thing jb_sv measured as one phenomenon is, here, not a
phenomenon. The amplitude jb_sv tracks -- the fold rate at the front -- is
still the right observable for the GLOBAL law (it is what the kick sets),
and is the wrong one for a local law on this chain, because the front cell
trades fold for translation as it goes.

    quantity                        single (axis, here)   double (jb_sv)
    front width, cells               1 from cell 1         2 at cell 1, then 1
    driver left of 0.9               0.0009                0.3593
    front cell's share of KE         0.999 -> 0.989        (a train behind)
    local speed, first/last          1.000                 1.8767
    local fold amplitude, first/last 1.095                 1.9238
    p_local (jb_sv's fit)            0.00                  0.9880
    p_global (kick sweep)            0.994                 1 (jb_tr R3, 1.99x)
    play exponent                    -0.994                -1 (jb_tr R4)
    speed*play/kick                  0.995, spread 0.9%    (decelerates; not quoted)

WHAT IS NOT CLAIMED. A mechanism. jb_1e R2 shows the axis bond's two points
straddle the square opening, so the joint's own geometry depends on the fold
and each bond leaves a hinge; that a two-point elastic impact between such
cells transfers the fold like a cradle is a reading, not a row. The counts
above are what is measured.

SCOPE. One chain, sixteen cells along an axis, free ends. Elastic, e = 1,
as jb_sv (a fully inelastic medium absorbs). Clearance 0.05, exaggerated
as everywhere in the contact line so the LCP is exercised. Point masses via
jb_mj; the lamina peer is not swept. V = 0: the impact law is the whole
dynamics. Every speed inherits the standing convention; ratios and
exponents are the measurements. jb_sv, jb_tr and jb_ct are untouched and
NOT retracted: they measure the double-covered chain, which is what the
double covering connects.
"""
from __future__ import annotations

import sys

import numpy as np

import jb_mj_inertial_honeycomb as MJ
import jb_rc_reduced as RC
import jb_sv_sonic_vacuum as SV

A_REF = MJ.A_REF
CELLS = SV.CELLS
KICK = SV.KICK
PLAY = SV.PLAY
STEP = SV.STEP

#: jb_sv's published numbers, T2 [23673], the instrument check.
SV_CELLS, SV_VR, SV_AR, SV_P, SV_DRIVER = 12, 1.8767, 1.9238, 0.9880, 0.3593

#: The sweeps. Time budgets scale with play / kick because the speed does.
KICKS = (0.45, 0.9, 1.8)
PLAYS = (0.025, 0.05, 0.1)


def chain(ncells=CELLS, gc=A_REF):
    """jb_1e's axis chain: the line the single covering connects."""
    asm, _ = RC.honeycomb_single([(2 * k, 0, 0) for k in range(ncells)], gc=gc)
    return asm


def momentum(asm, J, u):
    """Linear momentum of the whole chain under the point-mass model."""
    p = np.zeros(3)
    for k in range(asm.N):
        v = (J[k] @ u[k]).reshape(-1, 3)
        p += (RC.VMASS[:, None] * v).sum(axis=0)
    return p


def history(asm, kick=KICK, play=PLAY, e=SV.RESTITUTION, tmax=SV.TMAX,
            h=STEP, bands=True):
    """`jb_sv.front_history`, step for step, recording in addition the front
    cell's full velocity state, the chain's energy split, and momentum.

    R2 asserts the (t, cell, amplitude, width, driver) columns agree with
    `SV.front_history` on the same chain to 1e-12, so nothing here is a
    second instrument: it is jb_sv's, with more columns read off it.
    """
    pairs = MJ.tied_pairs(asm)
    n = asm.N
    q = asm.q0()
    u = np.zeros((n, 7))
    u[0, 6] = kick
    now, seen, rec = 0.0, 1, []
    J, M, Minv = MJ.kinematics(asm, q, False)
    e0 = MJ.kinetic(M, u)
    p_worst = 0.0
    while now < tmax - 1e-12:
        J, M, Minv = MJ.kinematics(asm, q, False)
        a1 = MJ.free_accel(asm, q, u, J, Minv, False)
        u_h = u + 0.5 * h * a1
        q_h = RC.apply_increment(asm, q, (0.5 * h * u).ravel())
        Jh, _, Mih = MJ.kinematics(asm, q_h, False)
        a2 = MJ.free_accel(asm, q_h, u_h, Jh, Mih, False)
        u = u + h * a2
        q = RC.apply_increment(asm, q, (h * u_h).ravel())
        now += h
        J, M, Minv = MJ.kinematics(asm, q, False)
        if bands:
            s = MJ.separations(asm, q, pairs)
            N = MJ.band_rows(asm, q, J, pairs)
            rate = np.dot(N, u.ravel())
            act = [i for i in range(len(pairs))
                   if s[i] >= play and rate[i] > 0]
            if act:
                u, _, _, _ = MJ.resolve(asm, u, N, act, Minv, e)
        p_worst = max(p_worst, float(np.abs(momentum(asm, J, u)).max()))
        p = np.abs(u[:, 6])
        hot = np.where(p > SV.FRONT_THRESH * abs(kick))[0]
        if not len(hot):
            continue
        f = int(hot.max())
        if f >= seen:
            seen = f + 1
            j, cells = f, 1
            while (j - 1 >= 0 and p[j - 1] < p[j]
                   and p[j - 1] > SV.SHOULDER * p[f]):
                j -= 1
                cells += 1
            ke = [0.5 * u[k] @ M[k] @ u[k] for k in range(n)]
            rec.append(dict(
                t=now, cell=f, amp=float(p[f]), width=cells,
                driver=float(p[0]), vx=float(u[f, 0]),
                vperp=float(np.hypot(u[f, 1], u[f, 2])),
                w=float(np.linalg.norm(u[f, 3:6])),
                ke_f=float(ke[f]), ke_fold=float(0.5 * M[f][6, 6] * u[f, 6] ** 2),
                ke_wake=float(sum(ke[:f])), ke_tot=float(sum(ke))))
    _, Mf, _ = MJ.kinematics(asm, q, False)
    return rec, dict(E0=e0, E1=MJ.kinetic(Mf, u), p_worst=p_worst)


def local_fit(rec):
    """jb_sv R2/R3's reduction: (v first/last, A first/last, p_local)."""
    t = np.array([r["t"] for r in rec])
    a = np.array([r["amp"] for r in rec])
    v = 1.0 / np.diff(t)
    a_mid = a[1:]
    p = float(np.polyfit(np.log(a_mid), np.log(v), 1)[0])
    return float(v[0] / v[-1]), float(a_mid[0] / a_mid[-1]), p


def speed(rec):
    """Linear fit of cell against time -- quoted only where R4 shows the
    front is uniform, which on this chain it is."""
    t = np.array([r["t"] for r in rec])
    c = np.array([r["cell"] for r in rec], float)
    return float(np.polyfit(t, c, 1)[0])


def gate():
    checks, out = [], {}
    A = checks.append

    # ---- R1: the instrument reproduces jb_sv --------------------------------
    d = SV.front_history(SV.chain())
    dt = [r[0] for r in d]
    da = [r[2] for r in d]
    dv = 1.0 / np.diff(dt)
    d_vr, d_ar = float(dv[0] / dv[-1]), float(da[1] / da[-1])
    d_p = float(np.polyfit(np.log(np.array(da[1:])), np.log(dv), 1)[0])
    d_drv = float(d[-1][4])
    out["R1"] = (len(d), d_vr, d_ar, d_p, d_drv)
    A(("R1  THE INSTRUMENT REPRODUCES jb_sv's PUBLISHED NUMBERS ON jb_sv's OWN "
       "CHAIN BEFORE IT IS POINTED AT THE SINGLE COVERING. Twelve cells "
       "reached, local speed falling 1.8767x while local amplitude falls "
       "1.9238x, p = 0.9880, the driver holding 0.3593 -- T2 [23673] to "
       "every quoted digit. [23746] lesson 4 again: check the instrument "
       "against the published number first. TWO-SIDED: any drift fails",
       len(d) == SV_CELLS and abs(d_vr - SV_VR) < 1e-3
       and abs(d_ar - SV_AR) < 1e-3 and abs(d_p - SV_P) < 1e-3
       and abs(d_drv - SV_DRIVER) < 1e-3,
       f"double chain: {len(d)} cells, v falls {d_vr:.4f}x, A falls "
       f"{d_ar:.4f}x, p = {d_p:.4f}, driver {d_drv:.4f}",
       f"{SV_CELLS} cells, {SV_VR}, {SV_AR}, {SV_P}, {SV_DRIVER}"))

    # ---- the reference run on the single chain, read two ways ---------------
    sgl = chain()
    ref, ref_e = history(sgl)
    sv_ref = SV.front_history(sgl)
    out["ref"], out["ref_e"] = ref, ref_e
    same = (len(sv_ref) == len(ref) and max(
        max(abs(a[0] - b["t"]), abs(a[2] - b["amp"]), abs(a[4] - b["driver"]))
        for a, b in zip(sv_ref, ref)) < 1e-12
        and all(a[1] == b["cell"] and a[3] == b["width"]
                for a, b in zip(sv_ref, ref)))

    # ---- R2: solitary, and a clean handoff ----------------------------------
    widths = sorted({r["width"] for r in ref})
    drv = ref[-1]["driver"]
    share = [r["ke_f"] / r["ke_tot"] for r in ref]
    out["R2"] = (widths, drv, share, same, sgl.N, len(sgl.welds),
                 len(MJ.tied_pairs(sgl)))
    A(("R2  THE FRONT IS ONE CELL WIDE AT EVERY CELL AND THE DRIVER HANDS IT "
       "OFF CLEAN. On jb_1e's axis chain -- sixteen cells, fifteen two-pair "
       "welds, thirty bands against jb_sv's forty-five -- the leading front "
       "is one cell wide from cell 1 to cell 15, where jb_sv's is two wide at "
       "cell 1. And the caveat jb_sv R4 carries as an honest limitation does "
       "not carry here: its driver keeps 0.3593 of 0.9 and sheds a train; "
       "this driver keeps under 0.01 and the front cell holds over 98% of "
       "the chain's kinetic energy at every cell it reaches. The columns "
       "this file adds are read off jb_sv's own stepper: the (t, cell, "
       "amplitude, width, driver) record agrees with SV.front_history on "
       "this chain to 1e-12. TWO-SIDED: a width above one, a driver keeping "
       "more than 1%, or a front cell below 98% of the energy fails",
       same and len(ref) == CELLS - 1 and widths == [1]
       and drv < 0.01 * KICK and min(share) > 0.98,
       f"{len(ref)} cells reached, widths {widths}; driver left "
       f"{drv:.4f} of {KICK}; front cell's share of KE {share[0]:.4f} -> "
       f"{share[-1]:.4f}; agrees with SV.front_history: {same}",
       "width 1 everywhere, driver under 1%, front share over 98%"))

    # ---- R3: the global laws, and their collapse -----------------------------
    runs = {}
    for k in KICKS:
        rec, _ = history(sgl, kick=k, tmax=SV.TMAX * KICK / k + 0.2)
        runs[("kick", k)] = rec
    for pl in PLAYS:
        if pl == PLAY:
            runs[("play", pl)] = runs[("kick", KICK)]
            continue
        rec, _ = history(sgl, play=pl, tmax=SV.TMAX * pl / PLAY + 0.2)
        runs[("play", pl)] = rec
    vk = {k: speed(runs[("kick", k)]) for k in KICKS}
    vp = {pl: speed(runs[("play", pl)]) for pl in PLAYS}
    p_glob = float(np.polyfit(np.log(KICKS), np.log([vk[k] for k in KICKS]),
                              1)[0])
    q_play = float(np.polyfit(np.log(PLAYS), np.log([vp[p] for p in PLAYS]),
                              1)[0])
    coll = [vk[k] * PLAY / k for k in KICKS] + [vp[p] * p / KICK for p in PLAYS]
    coll_spread = max(coll) / min(coll) - 1.0
    complete = all(len(r) == CELLS - 1 for r in runs.values())
    out["R3"] = (vk, vp, p_glob, q_play, coll)
    A(("R3  THE GLOBAL SONIC-VACUUM LAWS HOLD ON THE SINGLE CHAIN, AND THEY "
       "COLLAPSE ONTO ONE NUMBER. Speed against kick over a fourfold drive "
       "gives exponent 1; speed against clearance over fourfold gives -1; "
       "and speed * play / kick is the same to under 1% in every run, which "
       "is jb_ct's law -- speed = |dV/dgamma| gammadot / play -- with the "
       "front reaching the far end in every run. jb_tr R3/R4 measured both "
       "on the double chain (1.99x for 2x; speed*play constant to 0.8%). "
       "Here the front is uniform (R4), so the linear fit IS a speed and the "
       "exponents are read from it. TWO-SIDED: an exponent off by 0.02, a "
       "collapse spread above 1%, or a run that never reached the end fails",
       complete and abs(p_glob - 1.0) < 0.02 and abs(q_play + 1.0) < 0.02
       and coll_spread < 0.01,
       "speed by kick " + ", ".join(f"{k:g}: {vk[k]:.3f}" for k in KICKS)
       + f" (p = {p_glob:.4f}); by play "
       + ", ".join(f"{p:g}: {vp[p]:.3f}" for p in PLAYS)
       + f" (exponent {q_play:.4f}); speed*play/kick "
       + ", ".join(f"{c:.4f}" for c in coll) + f" (spread {100 * coll_spread:.2f}%)",
       "p = 1 and -1 within 0.02, collapse within 1%"))

    # ---- R4: the local law does not hold, and where the fold goes -----------
    vr, ar, p_loc = local_fit(ref)
    vx = [r["vx"] for r in ref]
    rot = max(max(r["vperp"], r["w"]) for r in ref)
    conv = [r["ke_f"] - r["ke_fold"] for r in ref]
    mono_vx = all(abs(vx[i + 1]) > abs(vx[i]) for i in range(len(vx) - 1))
    mono_cv = all(conv[i + 1] > conv[i] for i in range(len(conv) - 1))
    out["R4"] = (vr, ar, p_loc, vx, rot, conv, ref_e["p_worst"])
    A(("R4  jb_sv's LOCAL LAW DOES NOT HOLD HERE, BECAUSE THE FOLD'S DECAY IS "
       "RECOIL INSIDE THE FRONT CELL, NOT ATTENUATION. jb_sv R2/R3 gate that "
       "within one run the speed falls by the factor the amplitude falls by "
       "and fit p = 0.9880 from it. On the single chain the transit time per "
       "cell is constant to the step -- the speed's first/last ratio is 1 "
       "within 1% -- while the fold amplitude at the front decays by over "
       "8%, so jb_sv's fit returns p_local of essentially zero. The decayed "
       "fold energy is measured where it went: the front cell's velocity "
       "ALONG the chain grows in magnitude at every handoff, its transverse "
       "and angular velocities stay at zero to 1e-9, and KE_f - KE_fold "
       "grows monotonically while the cell's total stays within 2% of the "
       "chain's. Linear momentum along the chain is conserved to 1e-9 "
       "throughout, so the recoil is balanced by the wake's slow drift the "
       "other way. TWO-SIDED: a speed that tracked the amplitude (ratio "
       "near 1) would mean jb_sv's law holds and this row fails; a fold "
       "energy that went to the wake, or into rotation, fails it too",
       abs(vr - 1.0) < 0.01 and ar > 1.08 and abs(p_loc) < 0.1
       and mono_vx and mono_cv and rot < 1e-9 and ref_e["p_worst"] < 1e-9,
       f"local speed first/last {vr:.4f} (jb_sv {SV_VR}); fold amplitude "
       f"first/last {ar:.4f} (jb_sv {SV_AR}); p_local = {p_loc:+.4f} (jb_sv "
       f"{SV_P}); front cell vx {vx[0]:+.4f} -> {vx[-1]:+.4f} monotone; "
       f"transverse/angular {rot:.1e}; KE_f - KE_fold {conv[0]:.4f} -> "
       f"{conv[-1]:.4f}; worst |momentum| {ref_e['p_worst']:.1e}",
       "speed constant, amplitude falling, recoil growing, no rotation, "
       "momentum conserved"))

    # ---- R5: energy, and the wake ---------------------------------------------
    wake = [r["ke_wake"] / r["ke_tot"] for r in ref]
    e_drift = abs(ref_e["E1"] - ref_e["E0"]) / ref_e["E0"]
    out["R5"] = (e_drift, wake)
    A(("R5  ENERGY IS CONSERVED AND THE WAKE IS ONE PERCENT OF IT. V = 0 and "
       "e = 1, so the run must conserve kinetic energy to the integrator's "
       "drift; and the reformulated criterion's part 3, attenuation with "
       "distance, is measured here as an energy statement: the cells behind "
       "the front hold under 2% of the chain's energy when the front reaches "
       "the far end, against jb_tr R5's 3.27x amplitude decay over nine "
       "cells and jb_sv's driver keeping 0.36. TWO-SIDED: an energy drift "
       "above 1e-5 means the impacts are not elastic, and a wake above 2% "
       "means the front IS attenuating and R4's reading is wrong",
       e_drift < 1e-5 and max(wake) < 0.02
       and all(wake[i + 1] >= wake[i] for i in range(len(wake) - 1)),
       f"E {ref_e['E0']:.5f} -> {ref_e['E1']:.5f} (drift {e_drift:.1e}); "
       f"wake share {wake[0]:.4f} -> {wake[-1]:.4f}, monotone",
       "drift under 1e-5, wake under 2%"))

    # ---- R6: controls -----------------------------------------------------------
    still, _ = history(sgl, kick=0.0)
    free, _ = history(sgl, bands=False)
    out["R6"] = (len(still), len(free))
    A(("R6  BOTH CONTROLS, ON THIS CHAIN. Undriven, no front forms. Bands "
       "disabled, no front forms either: the driven cell keeps its fold rate "
       "and nothing reaches a neighbour, so every measurement above is a "
       "property of the two-pair contact and not of the drive, the index, or "
       "the hinge freedom",
       len(still) == 0 and len(free) == 0 and len(ref) == CELLS - 1,
       f"undriven: {len(still)} cells reached; bands disabled: {len(free)}; "
       f"bands enabled: {len(ref)}",
       "no front in either control, a full front with coupling"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_1v -- one covering: the sonic vacuum on the axis chain, "
              "against jb_sv")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        print("\n  THE FRONT, CELL BY CELL (kick 0.9, play 0.05, e = 1)")
        print(f"   {'t':>7s} {'cell':>4s} {'amp':>7s} {'width':>5s} "
              f"{'driver':>7s} {'vx':>8s} {'KE_f':>8s} {'KE_fold':>8s} "
              f"{'wake':>7s}")
        for r in out["ref"]:
            print(f"   {r['t']:7.3f} {r['cell']:4d} {r['amp']:7.4f} "
                  f"{r['width']:5d} {r['driver']:7.4f} {r['vx']:8.4f} "
                  f"{r['ke_f']:8.5f} {r['ke_fold']:8.5f} "
                  f"{r['ke_wake'] / r['ke_tot']:7.4f}")

        n_d, d_vr, d_ar, d_p, d_drv = out["R1"]
        vr, ar, p_loc = out["R4"][:3]
        vk, vp, p_glob, q_play, coll = out["R3"]
        print("\n  SINGLE AGAINST DOUBLE")
        print(f"   {'':32s} {'single (axis)':>14s} {'double (jb_sv)':>15s}")
        print(f"   {'front width':32s} {'1 from cell 1':>14s} "
              f"{'2 at 1, then 1':>15s}")
        print(f"   {'driver left of 0.9':32s} "
              f"{out['R2'][1]:>14.4f} {d_drv:>15.4f}")
        print(f"   {'local speed first/last':32s} {vr:>14.4f} {d_vr:>15.4f}")
        print(f"   {'local amplitude first/last':32s} {ar:>14.4f} "
              f"{d_ar:>15.4f}")
        print(f"   {'p_local':32s} {p_loc:>14.4f} {d_p:>15.4f}")
        print(f"   {'p_global (kick sweep)':32s} {p_glob:>14.4f} "
              f"{'1 (jb_tr R3)':>15s}")
        print(f"   {'play exponent':32s} {q_play:>14.4f} "
              f"{'-1 (jb_tr R4)':>15s}")
        print(f"   {'speed*play/kick':32s} {np.mean(coll):>14.4f} "
              f"{'decelerates':>15s}")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * THE SONIC VACUUM SURVIVES THE COVERING: p = 1, the -1")
        print("     clearance law, and a one-cell solitary front.")
        print("   * THE LOCAL LAW DOES NOT. Speed is constant while the fold")
        print("     decays into the front cell's own recoil. jb_sv's 'one")
        print("     phenomenon' is a statement about the double-covered")
        print("     chain; on this one there is neither deceleration nor")
        print("     attenuation to speak of.")
        print("   * THE HANDOFF IS CLEAN. jb_sv R4's caveat is the double")
        print("     chain's, not the medium's.")
        print("   * NO MECHANISM IS CLAIMED for the clean transfer.")
        print("   * jb_sv, jb_tr, jb_ct ARE NOT RETRACTED. Different chain.")
        print("   * ONE CHAIN, e = 1, play 0.05, point masses. Ratios only.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
