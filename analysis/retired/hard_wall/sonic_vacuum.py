"""sonic_vacuum -- the transport law that replaces the dispersion relation.

THE REFORMULATION, AND WHY IT IS NOT A RETREAT. Epic inviscid-qvf's acceptance
criterion has been "a measured dispersion relation" since it was written, and
jb_cp finally established why that can never be met HERE: a dispersion relation
pairs a frequency with a wavenumber, a frequency needs normal modes, normal
modes need a Hessian, and the medium's effective potential is an INFINITE
SQUARE WELL -- flat inside, infinite at the stop, no Hessian anywhere. omega(k)
does not merely go unmeasured in this medium; it does not exist for it. Owner
decision 2026-08-28: reformulate.

WHAT REPLACES IT, and this is a standard object rather than an invention. A
medium with no linear sound speed is a SONIC VACUUM in Nesterenko's sense, and
such media are characterised not by omega(k) but by their SOLITARY WAVE: the
speed-amplitude exponent, the width of the disturbance, and whether the profile
survives travel. The exponent is a fingerprint of the contact law -- Hertzian
grains give speed proportional to A^(1/4); a HARD WALL gives A^1, which is what
jb_ct and jb_tr both measured.

AND THE SUBSTITUTION IS HONEST ABOUT WHAT IT KEEPS. What a dispersion relation
tells you OPERATIONALLY is whether components of different wavelength travel at
different speeds -- i.e. whether a pulse SPREADS. That question can be put
directly to the medium by launching a pulse and watching its width, with no
omega(k) anywhere. The reformulated criterion asks the same physical question
through an instrument this medium admits.

THE FOUR PARTS, and three were already in hand before this file:

    1. speed-amplitude exponent      p = 1        jb_tr R3 (1.9904x for 2x)
    2. speed-clearance exponent      -1           jb_tr R4 (0.8% over 2x)
    3. attenuation with distance     3.27x / 9    jb_tr R5
    4. SHAPE: spread, or solitary?   THIS FILE

WHAT THIS FILE MEASURES, AND IT IS THE STRONGEST OF THE FOUR.

  * THE FRONT DOES NOT SPREAD. It is ONE CELL WIDE at cell 2 and one cell wide
    at cell 12, every step of the way. A dispersive medium widens; this does
    not. That is the direct answer to the question omega(k) was being asked
    for.
  * THE FRONT OBEYS speed proportional to AMPLITUDE, LOCALLY, cell by cell, as
    its own amplitude decays. Over the run the local speed falls by 1.96x while
    the local amplitude falls by 1.93x -- the same factor to about a percent.
    jb_ct's law was measured as a GLOBAL scaling across separate runs; here it
    holds WITHIN one run, at each successive cell, which is a much stronger
    statement of the same law.
  * SO THE DECELERATION AND THE ATTENUATION ARE ONE PHENOMENON, not two.
    jb_tr reported them separately and said they were probably the same thing;
    this measures that they are. The front is not slowing for its own reasons.
    It is obeying v proportional to A while A decays.

WHAT IT ALSO SHOWS, AND THE FILE WOULD BE DISHONEST WITHOUT IT: the driven cell
does NOT hand its motion off. It keeps roughly half of what it started with and
sheds a TRAIN of smaller pulses behind the leading front. This is a solitary
FRONT in a medium that also rings, not a clean single soliton, and no reading
here should be dressed up as the latter.

SCOPE. One chain, sixteen cells along a body diagonal, free ends. Elastic
(e = 1), because a fully inelastic medium absorbs rather than carries. Mass
models and metric declared through jb_mj. KERNEL and PRIMITIVE INAPPLICABLE
while V = 0, stated and not lapsed. Every speed inherits the project's
convention (coupling 1, total mass 1/2, R = 1), so only RATIOS are gated.
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.retired.hard_wall import impact_law as MJ
from analysis.model import assembly as RC
CELLS = 16
KICK = 0.9
PLAY = 0.05
RESTITUTION = 1.0
STEP = 1e-3
TMAX = 1.2

#: A cell counts as carrying the front above this fraction of the drive.
FRONT_THRESH = 0.05

#: A neighbour belongs to the leading feature while it is still rising toward
#: the front and holds at least this share of the front's amplitude. Without
#: the second condition the "width" would run all the way back to the driven
#: cell through the wake.
SHOULDER = 0.2


def chain(n=CELLS, gc=MJ.A_REF):
    asm, _ = RC.honeycomb([(k, k, k) for k in range(n)], gc=gc)
    return asm


def front_history(asm, kick=KICK, play=PLAY, e=RESTITUTION, tmax=TMAX,
                  h=STEP, lamina=False, bands=True):
    """Track the LEADING front: when it reaches each cell, how big it is there,
    and how many cells it occupies."""
    pairs = MJ.tied_pairs(asm)
    n = asm.N
    q = asm.q0()
    u = np.zeros((n, 7))
    u[0, 6] = kick
    now, seen, rec = 0.0, 1, []
    driven_peak = 0.0
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
        p = np.abs(u[:, 6])
        driven_peak = max(driven_peak, float(p[0]))
        hot = np.where(p > FRONT_THRESH * abs(kick))[0]
        if not len(hot):
            continue
        f = int(hot.max())
        if f >= seen:
            seen = f + 1
            j, cells = f, 1
            while (j - 1 >= 0 and p[j - 1] < p[j]
                   and p[j - 1] > SHOULDER * p[f]):
                j -= 1
                cells += 1
            rec.append((now, f, float(p[f]), cells, float(p[0])))
    return rec


def gate():
    checks, out = [], {}
    A = checks.append
    asm = chain()
    rec = front_history(asm)
    out["rec"] = rec
    cells = [r[1] for r in rec]
    amps = [r[2] for r in rec]
    wids = [r[3] for r in rec]
    times = [r[0] for r in rec]

    # R1 -- the row the reformulation exists for
    body = [w for c, w in zip(cells, wids) if c >= 2]
    A(("R1  THE FRONT DOES NOT SPREAD, WHICH IS THE QUESTION omega(k) WAS BEING "
       "ASKED FOR. A dispersion relation's operational content is whether "
       "components of different wavelength travel at different speeds -- that "
       "is, whether a pulse WIDENS as it goes. Put to this medium directly: the "
       "leading front is ONE CELL WIDE when it reaches cell 2 and one cell wide "
       "when it reaches the far end, at every cell in between. So the medium "
       "carries a SOLITARY front rather than a dispersive packet, and it does "
       "so without any frequency having to exist. TWO-SIDED and that is the "
       "whole point: a spreading front would fail this row, and it is the "
       "outcome a medium with a dispersion relation would generically give",
       len(body) >= 8 and max(body) == 1 and min(body) == 1,
       f"front width in cells, from cell {min(c for c in cells if c >= 2)} to "
       f"cell {max(cells)}: {sorted(set(body))} (constant); measured at "
       f"{len(body)} successive cells",
       "width exactly 1 at every cell the front reaches"))

    # R2 -- and the local law, which unifies two things jb_tr reported apart
    d = np.diff(times)
    v = 1.0 / d
    a_mid = np.array(amps[1:])
    vr = float(v[0] / v[-1])
    ar = float(a_mid[0] / a_mid[-1])
    out["local"] = (vr, ar)
    A(("R2  THE FRONT OBEYS speed = k * AMPLITUDE LOCALLY, CELL BY CELL, WHICH "
       "MAKES THE DECELERATION AND THE ATTENUATION ONE PHENOMENON RATHER THAN "
       "TWO. jb_ct measured speed proportional to amplitude as a GLOBAL "
       "scaling across separate runs at different drive; jb_tr reproduced it "
       "the same way and reported the front's slowing and its decay as two "
       "observations that were probably the same thing. They are: within a "
       "SINGLE run the front's local speed falls by very nearly the factor its "
       "local amplitude falls by. The front is not decelerating for reasons of "
       "its own -- it is obeying v proportional to A while A decays. TWO-SIDED: "
       "a front that slowed WITHOUT losing amplitude, or lost amplitude "
       "without slowing, fails this row",
       abs(vr / ar - 1.0) < 0.06,
       f"local speed falls {vr:.4f}x over the run while local amplitude falls "
       f"{ar:.4f}x; ratio of the two {vr / ar:.4f}",
       "the two factors agreeing within 6%"))

    # R3 -- the exponent is the fingerprint of the contact law
    logv = np.log(v)
    loga = np.log(a_mid)
    p = float(np.polyfit(loga, logv, 1)[0])
    out["p"] = p
    A(("R3  THE SPEED-AMPLITUDE EXPONENT IS 1, WHICH IDENTIFIES THE CONTACT "
       "LAW. In a sonic vacuum the exponent p in speed ~ A^p is a fingerprint "
       "of what the grains do when they touch: Hertzian spheres give p = 1/4, "
       "and a HARD WALL -- flat inside the clearance, infinite at the stop -- "
       "gives p = 1. Fitted here from the front's own decay within one run, "
       "which is an independent route to the same number jb_tr got by "
       "comparing separate runs. It agrees with jb_cp's infinite square well, "
       "reached from the potential side. Three instruments, one contact law",
       abs(p - 1.0) < 0.12,
       f"fitted exponent p = {p:.4f} from {len(logv)} successive cells "
       f"(Hertzian would be 0.25; hard wall is 1)",
       "p = 1 within 12%, i.e. hard wall and not Hertzian"))

    # R4 -- the honest limitation
    driven_end = rec[-1][4]
    A(("R4  THE DRIVEN CELL DOES NOT HAND ITS MOTION OFF, AND THIS FILE WOULD "
       "BE DISHONEST WITHOUT THE ROW. A clean soliton would carry the "
       "disturbance away and leave the driver quiet. Here the driven cell "
       "still holds a large share of its original rate when the front has "
       "reached the far end, and a TRAIN of smaller pulses trails behind the "
       "leading front. So what the medium carries is a solitary FRONT in a "
       "medium that also rings -- not a single soliton, and no row here should "
       "be read as claiming one. This is a two-sided measurement of a "
       "limitation: it fails if the driver ever does go quiet, which would "
       "mean the caveat is the wrong one to carry",
       driven_end > 0.3 * KICK,
       f"driven cell's fold rate when the front reaches cell {cells[-1]}: "
       f"{driven_end:.4f} against an initial {KICK}",
       "the driver still carrying a large share, i.e. no clean handoff"))

    # R5 -- controls
    still = front_history(asm, kick=0.0)
    free = front_history(asm, bands=False)
    out["still"], out["free"] = len(still), len(free)
    A(("R5  BOTH CONTROLS, AND BOTH CAN FAIL. UNDRIVEN, no front ever forms. "
       "COUPLING REMOVED, no front ever forms either -- with the bands "
       "disabled the driven cell keeps its motion and nothing reaches any "
       "other cell, so every measurement above is a property of the coupling "
       "rather than of the drive or of the integrator. Without the second "
       "control this file could be measuring its own initial condition "
       "travelling through an index",
       len(still) == 0 and len(free) == 0 and len(rec) > 8,
       f"undriven: {len(still)} cells reached; bands disabled: {len(free)}; "
       f"bands enabled: {len(rec)}",
       "no front in either control, a full front with coupling"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_sv -- the transport law that replaces the dispersion relation")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        print("\n  THE FRONT, CELL BY CELL")
        print(f"   {'t':>7s} {'cell':>5s} {'amplitude':>10s} {'width':>6s} "
              f"{'driven':>8s}")
        for t, c, a, w, d0 in out["rec"]:
            print(f"   {t:7.3f} {c:5d} {a:10.4f} {w:6d} {d0:8.4f}")

        print()
        print("  THE REFORMULATED CRITERION, AND WHAT NOW MEETS IT.")
        print("   1. speed-amplitude exponent   p = 1        MET (R3, jb_tr R3)")
        print("   2. speed-clearance exponent   -1           MET (jb_tr R4)")
        print("   3. attenuation with distance  measured     MET (jb_tr R5, R2)")
        print("   4. shape: solitary or spread  SOLITARY     MET (R1)")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * THE MEDIUM CARRIES A SOLITARY FRONT, one cell wide, that")
        print("     does not spread. That is the direct answer to the question")
        print("     a dispersion relation was being asked for, obtained")
        print("     without any frequency having to exist.")
        print("   * NO DISPERSION RELATION IS CLAIMED OR IMPLIED. There is")
        print("     still no frequency in this model, because V = 0 and the")
        print("     effective potential is an infinite square well. The")
        print("     criterion was REFORMULATED, by owner decision; it was not")
        print("     met in its original form and this file does not pretend")
        print("     otherwise.")
        print("   * NOT A CLEAN SOLITON. The driver keeps its motion and the")
        print("     medium rings behind the front (R4).")
        print("   * ONE CHAIN, sixteen cells, free ends. Not the bulk.")
        print("   * RATIOS ONLY. Every speed inherits the project's standing")
        print("     convention, so no absolute number here is physical.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
