"""jb_cp -- the medium's V is a hard wall, and that is why the eigensolver is idle.

WHY THIS FILE EXISTS. Bead qvf.30 assembles a V lead from four places that never
cited each other: T2 [22331]'s untested note that "a conical minimum is not
automatically illegitimate (contact/impact potentials look like this)"; the
exhaustive failure of the smooth-potential survey jb_h..jb_q; contact becoming
the coupling (jb_ct, T2 [23639]); and InternalMassMetric + GeneralizedEigensolver
sitting BUILT, "CHOICE-FREE SCAFFOLDING", and IDLE.

The bead names the obvious tension: a conical minimum has no Hessian, and
GeneralizedEigensolver wants one. It asks whether the frequencies come from
somewhere else, or whether the expansion point is wrong.

BOTH ANSWERS ARE NO, AND THE THIRD ONE IS THE RESULT.

FIRST, A STALE PREMISE CORRECTED, because the bead and T2 [23643] section III
both carry it: the scaffolding is NOT "waiting for bead qvf.2 to choose a
potential". qvf.2 CLOSED on 2026-08-15 having chosen one -- DECISION 17, the raw
all-pairs kernel on real distances -- with six real frequencies at the VE and
zero zero-modes at every tolerance from 1e-4 to 1e-12. Bead inviscid-6dp, the
scaffolding itself, is closed too. What is true is narrower and more
interesting, and both halves were checked rather than assumed:

  * those six frequencies were computed in PYTHON (attic/jb_s_frequency_spectrum
    and the jb_o / jb_t / jb_u family), never through the Java scaffolding.
    GeneralizedEigensolver still has no caller anywhere outside its own file and
    its own tests, and its javadoc's "once bead inviscid-qvf.2 chooses a
    potential" has been stale for a fortnight.
  * DECISION 17's V is a SINGLE-UNIT potential. What has no V is the ARRAY,
    whose coupling is contact -- and that is the object this file measures.

So the scaffolding is idle for two stacked reasons, neither of them the one on
record: the single-unit spectrum was computed somewhere else, and the array's
potential turns out to be the one shape it cannot consume.

THE INSTRUMENT: THE EXPONENT, NOT THE HESSIAN. For V = k|x|^n the period of a
bounded orbit is T ~ A^(1 - n/2), so

    omega ~ A^(n/2 - 1)

which exists for EVERY n. The Hessian is not what makes a frequency; it is the
n = 2 special case where the frequency stops depending on amplitude. So "no
Hessian" never meant "no frequency" -- it meant "no amplitude-independent
frequency", which is a different and much weaker statement, and it is one this
programme has already measured from the other side as "no linear sound speed"
(jb_ct R5). R1 validates the exponent instrument against four known answers
before it is pointed at anything.

WHAT THE MEDIUM MEASURES. jb_ct's chain is V = 0 inside the band with an elastic
impact at |c| = t. That IS an infinite square well in the constraint coordinate,
and R2 measures its exponents: the impact rate goes as kick^(+1.000000) over a
sixteenfold range and as play^(-1) over an eightfold one, with rate*play/kick
flat to 0.4%. A harmonic minimum would give (0, 0). A CONICAL one would give
(0, -1/2). The medium gives (+1, -1) -- a hard wall, and neither of the others.

    THE EFFECTIVE POTENTIAL IS  V(c) = 0 for |c| < t,  +infinity at |c| = t.

AND IT IS NOT Vol_hull'S CONE. [22331]'s local model at the VE is
V/V_tet = 20 + c|a| - (k/2)a^2 with c = 4 sqrt(3): degree-one, conical, n = 1.
That predicts a width exponent of -1/2. The medium measures -1.03. They are a
factor of two apart, so the cone point is NOT this medium's effective potential
(R4). The bead asked for [22331] to be re-read before discarding again; this is
a measurement instead, and it discards it on evidence rather than on the old
"no Hessian, therefore useless" reasoning -- which R3 shows was never a good
reason.

SO THE EIGENSOLVER IS IDLE BECAUSE THERE IS NOTHING FOR IT TO DO, and that is a
RESULT rather than a gap in the scaffolding. A hard wall has no Hessian anywhere
-- identically flat inside, infinite at the boundary -- so there is no linear
regime to linearise about and no normal-mode spectrum to compute. Bead qvf.2's
premise, that choosing a potential would let InternalMassMetric and
GeneralizedEigensolver produce frequencies, does not survive the potential the
medium actually has.

THE WALLS MOVE (R5). The hard-wall law is EXACT instantaneously -- the first
crossing takes 2t/(|V'| gammadot) to 0.54% -- and then the intervals lengthen,
17% over 61 impacts, because |V'(gamma)| = Z|sin gamma| IS the wall speed and
gamma winds. Predicted from the drift alone: 1.1774 against 1.1693 measured. The
box is rigid in c and has a configuration-dependent width in gamma, which is the
same phase dependence jb_pr measured in the coupling coefficient C(a0) and jb_ct
in the front's uniformity.

WHAT THIS FILE COULD NOT MEASURE, AND WHY IT IS THE INTERESTING PART. A
dispersion relation needs a plane wave. Launch one and EVERY joint reaches its
stop at the same instant -- 7 of 7 at theta = pi, 4 at theta = pi/2, against
exactly 1 for the single-end kick jb_ct was built around (R6). jb_ct's
integrator resolves ONE contact per event, which is correct for a staggered
front and is the wrong impact law for a plane wave: simultaneous multi-contact
needs the LCP that T2 [23230] DECISION 18 specified as Phase 2 and that was
never built (beads qvf.21 / qvf.22). So the dispersion relation is not merely
unmeasured here, it is BLOCKED on a named, specified, unbuilt piece -- and R6 is
the measurement that says so rather than an assertion that it does.

SCOPE.
  * The effective V is read off jb_ct's SCALAR radial model. jb_pr showed that
    reduction understates the joint's coupling by sqrt(3); it does not change
    any exponent here, because every exponent is a ratio.
  * One joint (R2, R5) and an 8-cell chain (R6). Nothing here is the bulk.
  * "Frequency" throughout means the reciprocal of a turning or bounce time. No
    normal mode is computed, because R3 is the finding that there are none.
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.integrate import solve_ivp

import jb_ct_contact_chain as CT

Z = CT.Z
A_REF = -30.0

#: Wall speed at the reference: |V'(-30)| = Z sin(30). The rate at which a
#: joint's clearance is consumed per unit of fold rate.
WALL = Z * np.sin(np.radians(30.0))


def power_omega(n, A, k=1.0, m=1.0):
    """Angular frequency of a bounded orbit in V = k|x|^n, released at rest
    from x = A. MEASURED by timing the return to the turning point, not taken
    from the closed form -- the closed form is what R1 tests it against."""
    def f(_t, y):
        return [y[1], -(k * n * np.abs(y[0]) ** (n - 1) * np.sign(y[0])) / m]

    def turn(_t, y):
        return y[1]
    turn.direction = 1
    turn.terminal = True
    E = k * abs(A) ** n
    tmax = 400.0 * A / max(np.sqrt(2 * E / m), 1e-12)
    s = solve_ivp(f, [0, tmax], [A, 0.0], events=turn, rtol=1e-12, atol=1e-14)
    return 2 * np.pi / float(s.t_events[0][0])


def _step(g, gd, sep0, play, h):
    """One event-driven step of the two-cell joint. Same law as jb_ct.chain --
    bisect to the crossing, elastic impulse along the constraint gradient --
    reproduced here rather than imported because this file needs the impact
    TIMES, which jb_ct only counts."""
    ng, ngd = CT._rk4(g, gd, h)

    def viol(x):
        return abs(CT.radial(x[0]) + CT.radial(x[1]) - sep0) - play

    if viol(ng) <= 0:
        return ng, ngd, h, False
    lo, hi = 0.0, h
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        tg, _ = CT._rk4(g, gd, mid)
        if viol(tg) > 0:
            hi = mid
        else:
            lo = mid
    g, gd = CT._rk4(g, gd, lo)
    c = CT.radial(g[0]) + CT.radial(g[1]) - sep0
    u = CT.radial_d(g[0]) * gd[0] + CT.radial_d(g[1]) * gd[1]
    hit = False
    if (c > 0 and u > 0) or (c < 0 and u < 0):
        lam = -2.0 * u / (CT.radial_d(g[0]) ** 2 / CT.inertia(g[0])
                          + CT.radial_d(g[1]) ** 2 / CT.inertia(g[1]))
        gd[0] += lam * CT.radial_d(g[0]) / CT.inertia(g[0])
        gd[1] += lam * CT.radial_d(g[1]) / CT.inertia(g[1])
        hit = True
    return g, gd, lo, hit


def joint(kick=0.30, play=0.02, tmax=16.0, a0=A_REF, dt=1e-3):
    """One joint, two cells. Returns (impact times, final g, initial g)."""
    g = np.array([np.radians(a0), np.radians(a0 + 60.0)])
    gd = np.array([kick, 0.0])
    g0 = g.copy()
    sep0 = CT.radial(g[0]) + CT.radial(g[1])
    now, ts = 0.0, []
    while now < tmax - 1e-12:
        g, gd, adv, hit = _step(g, gd, sep0, play, min(dt, tmax - now))
        now += adv
        if hit:
            ts.append(now)
    return np.array(ts), g, g0


def simultaneous(gd0, play=0.01, a0=A_REF, dt=1e-3, events=24):
    """How many joints stand AT their stop when an event fires. The integrator
    resolves one; this counts how many it should have resolved."""
    n = len(gd0)
    g = np.array([np.radians(a0 + 60.0 * (k % 2)) for k in range(n)])
    gd = np.array(gd0, dtype=float)
    sep0 = CT.radial(g[0]) + CT.radial(g[1])
    now, seen, out = 0.0, 0, []

    def viol(x):
        return np.abs(CT.radial(x[:-1]) + CT.radial(x[1:]) - sep0) - play

    while seen < events and now < 50.0:
        ng, ngd = CT._rk4(g, gd, dt)
        if not (viol(ng) > 0).any():
            g, gd, now = ng, ngd, now + dt
            continue
        lo, hi = 0.0, dt
        for _ in range(35):
            mid = 0.5 * (lo + hi)
            tg, _ = CT._rk4(g, gd, mid)
            if (viol(tg) > 0).any():
                hi = mid
            else:
                lo = mid
        g, gd = CT._rk4(g, gd, lo)
        now += lo
        out.append(int(np.sum(viol(g) > -1e-9)))
        c = CT.radial(g[:-1]) + CT.radial(g[1:]) - sep0
        k = int(np.argmax(np.abs(c)))
        u = CT.radial_d(g[k]) * gd[k] + CT.radial_d(g[k + 1]) * gd[k + 1]
        if (c[k] > 0 and u > 0) or (c[k] < 0 and u < 0):
            lam = -2.0 * u / (CT.radial_d(g[k]) ** 2 / CT.inertia(g[k])
                              + CT.radial_d(g[k + 1]) ** 2 / CT.inertia(g[k + 1]))
            gd[k] += lam * CT.radial_d(g[k]) / CT.inertia(g[k])
            gd[k + 1] += lam * CT.radial_d(g[k + 1]) / CT.inertia(g[k + 1])
        seen += 1
    return out


def gate():
    checks, out = [], {}
    A = checks.append
    amps = np.array([0.25, 0.5, 1.0, 2.0, 4.0])

    # R1 -- the instrument, on four answers that are known and different
    exps = {}
    for n in (1.0, 2.0, 3.0, 4.0):
        w = np.array([power_omega(n, a) for a in amps])
        exps[n] = float(np.polyfit(np.log(amps), np.log(w), 1)[0])
    out["exps"] = exps
    A(("R1  THE EXPONENT IS THE INSTRUMENT, AND THE HESSIAN IS ONLY ITS n = 2 "
       "SPECIAL CASE. For V = k|x|^n a bounded orbit has omega ~ A^(n/2 - 1), "
       "which is finite and well defined for EVERY n -- the frequency comes "
       "from the turning time, and a Hessian is what you get when that turning "
       "time stops depending on amplitude. Measured by timing the return to "
       "the turning point, never from the closed form, at four exponents that "
       "are known and different: -1/2 (conical), 0 (harmonic), +1/2, +1. "
       "TWO-SIDED and the whole file rests on it: if the harmonic case came "
       "back amplitude-DEPENDENT, or any other came back flat, every verdict "
       "below would be the instrument's and not the medium's",
       max(abs(exps[n] - (n / 2 - 1)) for n in exps) < 1e-4,
       "; ".join(f"n={n:.0f}: measured {exps[n]:+.6f} against {n / 2 - 1:+.1f}"
                 for n in sorted(exps)),
       "all four exponents to 1e-4, harmonic flat and the others not"))

    # R2 -- the medium's own exponents say infinite square well
    kicks = np.array([0.075, 0.15, 0.30, 0.60, 1.20])
    krate = []
    for kk in kicks:
        ts, _, _ = joint(kick=kk, play=0.02, tmax=16.0 * 0.30 / kk)
        krate.append(len(ts) / (16.0 * 0.30 / kk))
    kexp = float(np.polyfit(np.log(kicks), np.log(krate), 1)[0])
    plays = np.array([0.0025, 0.005, 0.01, 0.02])
    prate = []
    for pp in plays:
        ts, _, _ = joint(kick=0.30, play=pp, tmax=16.0)
        prate.append(len(ts) / 16.0)
    pexp = float(np.polyfit(np.log(plays), np.log(prate), 1)[0])
    inv = [r * p / 0.30 for r, p in zip(prate, plays)]
    out["kexp"], out["pexp"] = kexp, pexp
    A(("R2  THE MEDIUM'S EFFECTIVE POTENTIAL IS AN INFINITE SQUARE WELL, and "
       "its exponents are what say so. jb_ct's chain is V = 0 inside the band "
       "with an elastic impact at |c| = t, which IS a hard wall in the "
       "constraint coordinate -- so the impact rate must go as kick^(+1) and "
       "play^(-1), and rate*play/kick must be a constant. Measured over a "
       "SIXTEENFOLD range in kick and an EIGHTFOLD one in play. TWO-SIDED "
       "against the two rivals this file exists to separate: a harmonic "
       "minimum gives (0, 0) and a CONICAL one gives (0, -1/2). Note the "
       "confound this row had to remove -- an earlier sweep scaled tmax with "
       "the play, which let the pair wind further at large clearance and bent "
       "the exponent to -1.17. Holding tmax fixed is what makes the play axis "
       "mean the play",
       abs(kexp - 1.0) < 0.02 and abs(pexp + 1.0) < 0.06
       and max(inv) / min(inv) < 1.02,
       f"kick exponent {kexp:+.6f} (predicts +1); play exponent {pexp:+.6f} "
       f"(predicts -1); rate*play/kick = "
       + ", ".join(f"{v:.4f}" for v in inv),
       "(+1, -1) with the invariant flat to 2%, and neither (0,0) nor (0,-1/2)"))

    # R3 -- so the bead's tension dissolves, and it is measurable that it does
    small = np.array([1e-3, 1e-2, 1e-1, 1.0])
    cone = np.array([power_omega(1.0, a) for a in small])
    A(("R3  \"NO HESSIAN\" NEVER MEANT \"NO FREQUENCY\", WHICH DISSOLVES THE "
       "BEAD'S TENSION WITHOUT MOVING THE EXPANSION POINT. A conical minimum "
       "has a perfectly well defined period at every nonzero amplitude -- "
       "measured here across three decades -- and what it lacks is a LIMIT as "
       "the amplitude goes to zero, so there is no single number to call THE "
       "frequency. That is what the missing Hessian actually costs. T2 [22331] "
       "read it as fatal to Vol_hull ('no linearisation, no frequency, no "
       "dispersion relation'), and the first two thirds of that sentence are "
       "too strong. TWO-SIDED: the frequencies must be finite at every "
       "amplitude AND must diverge as the amplitude shrinks -- a conical well "
       "that produced a finite limit would refute the reading, and one that "
       "produced no orbit at all would refute the row",
       np.all(np.isfinite(cone)) and cone[0] > cone[-1] * 10
       and abs(np.polyfit(np.log(small), np.log(cone), 1)[0] + 0.5) < 1e-4,
       "conical V = |x|, omega at A = 1e-3, 1e-2, 1e-1, 1: "
       + ", ".join(f"{w:.4f}" for w in cone)
       + f"; exponent {np.polyfit(np.log(small), np.log(cone), 1)[0]:+.6f}",
       "finite everywhere, diverging as A^(-1/2), no zero-amplitude limit"))

    # R4 -- and the medium is not the cone. A factor of two, measured.
    A(("R4  THE MEDIUM IS NOT Vol_hull'S CONE, AND THIS DISCARDS IT ON "
       "EVIDENCE RATHER THAN ON THE OLD REASON. [22331]'s local model at the "
       "VE is V/V_tet = 20 + c|a| - (k/2)a^2 with c = 4 sqrt(3): degree one, "
       "conical, n = 1, which R1 measures as a width exponent of -1/2. The "
       "medium measures -1. They differ by a factor of two, so the cone point "
       "is not this medium's effective potential. The bead asked for [22331] "
       "to be re-read under a contact interpretation before being discarded "
       "again; the honest answer is that it can be MEASURED against instead, "
       "and R3 has already removed the reason it was discarded the first time",
       abs(pexp - (-1.0)) < abs(pexp - (-0.5)) and abs(pexp + 0.5) > 0.4,
       f"medium {pexp:+.4f}; hard wall predicts -1 (distance "
       f"{abs(pexp + 1.0):.4f}); cone predicts -1/2 (distance "
       f"{abs(pexp + 0.5):.4f})",
       "nearer the wall than the cone, by more than the gap between them"))

    # R5 -- the walls are rigid in c and move in gamma
    ts, g, g0 = joint()
    gaps = np.diff(ts)
    first = 2 * 0.02 / (WALL * 0.30)
    w0 = abs(Z * np.sin(g0[0])) + abs(Z * np.sin(g0[1]))
    w1 = abs(Z * np.sin(g[0])) + abs(Z * np.sin(g[1]))
    out["gapratio"] = (float(gaps[-1] / gaps[0]), float(w0 / w1))
    A(("R5  THE WALLS ARE RIGID IN THE CONSTRAINT AND MOVE IN THE ANGLE. The "
       "hard-wall law is EXACT instantaneously: the first crossing takes "
       "2t/(|V'| gammadot) to within 0.6%. Then the intervals lengthen, "
       "because |V'(gamma)| = Z|sin gamma| IS the wall speed and gamma winds "
       "with nothing to restore it. Predicting the interval ratio from the "
       "measured drift alone lands within 1% of it, which is what makes this "
       "the mechanism rather than a coincidence -- an earlier attempt "
       "attributed the same drift to V's curvature ACROSS ONE CROSSING, which "
       "is a different and far smaller quantity, and correcting for it changed "
       "nothing. TWO-SIDED: the row fails if the first crossing stops matching "
       "the closed form OR if the drift stops being predicted by the wind",
       abs(gaps[0] / first - 1.0) < 0.01
       and abs((gaps[-1] / gaps[0]) / (w0 / w1) - 1.0) < 0.02,
       f"first interval {gaps[0]:.6f} against 2t/(|V'| gammadot) = "
       f"{first:.6f}; over {len(ts)} impacts the interval grows "
       f"{gaps[-1] / gaps[0]:.4f}x while the wall speed falls {w0:.5f} -> "
       f"{w1:.5f}, predicting {w0 / w1:.4f}",
       "closed form to 1%, and the drift predicted by the wind to 2%"))

    # R6 -- and the dispersion relation is blocked on a named, unbuilt piece
    n = 8
    census = {}
    for label, v in (("plane wave, theta = pi", 0.30 * np.cos(np.pi * np.arange(n))),
                     ("plane wave, theta = pi/2",
                      0.30 * np.cos(np.pi / 2 * np.arange(n))),
                     ("single end kick (jb_ct's own)",
                      np.array([0.30] + [0.0] * (n - 1)))):
        census[label] = simultaneous(v)
    out["census"] = census
    single = census["single end kick (jb_ct's own)"]
    plane = census["plane wave, theta = pi"]
    A(("R6  THE DISPERSION RELATION IS BLOCKED ON THE UNBUILT PHASE 2 IMPACT "
       "LAW, AND THIS IS THE MEASUREMENT THAT SAYS SO RATHER THAN AN "
       "ASSERTION. A dispersion relation needs a plane wave, and a plane wave "
       "on this medium puts EVERY joint at its stop in the same instant -- all "
       "7 of an 8-cell chain at theta = pi, 4 at theta = pi/2 -- against "
       "exactly 1 for the single-end kick jb_ct was built around. jb_ct's "
       "integrator resolves ONE contact per event, which is right for a "
       "staggered front and is the WRONG IMPACT LAW for simultaneous "
       "multi-contact: that needs the LCP T2 [23230] DECISION 18 specified as "
       "Phase 2 (0 <= lambda_N perp grad g . v+ >= 0) and which beads qvf.21 / "
       "qvf.22 never built. TWO-SIDED, and that is the point: the row fails if "
       "jb_ct's own initial condition ever needs more than one contact "
       "resolved, which would mean the existing measurements were wrong too",
       max(single) == 1 and min(plane) == n - 1,
       "; ".join(f"{k}: {v[:6]}" for k, v in census.items()),
       "exactly 1 for a single kick, all n-1 for a plane wave"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_cp -- the medium's V is a hard wall, and the eigensolver is idle")
    print("=" * 78)
    checks, out = gate()
    bad = 0
    for name, ok, got, want in checks:
        bad += 0 if ok else 1
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")

    print("\n  THE ANSWER TO BEAD qvf.30, IN ONE LINE:")
    print("   V(c) = 0 for |c| < t and +infinity at |c| = t. An infinite square")
    print("   well in the constraint coordinate -- not a cone, not a spring.")
    print()
    print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
    print("   * THE EIGENSOLVER IS IDLE BECAUSE THERE IS NOTHING FOR IT TO DO.")
    print("     A hard wall has no Hessian anywhere: identically flat inside,")
    print("     infinite at the boundary. That is a RESULT, not a gap in the")
    print("     scaffolding.")
    print("   * AND NOT FOR THE REASON ON RECORD. Bead qvf.30 and T2 [23643]")
    print("     section III both say the scaffolding waits on bead qvf.2 to")
    print("     choose a potential. qvf.2 closed 2026-08-15 having chosen one")
    print("     (DECISION 17, raw all-pairs kernel), with six real")
    print("     frequencies at the VE. Those were computed in PYTHON; the Java")
    print("     GeneralizedEigensolver still has no caller outside its own")
    print("     tests, and its javadoc's 'once qvf.2 chooses a potential' is a")
    print("     fortnight stale. DECISION 17's V is a SINGLE-UNIT potential --")
    print("     it is the ARRAY that has none, and the array's is this wall.")
    print("   * 'NO HESSIAN' DOES NOT MEAN 'NO FREQUENCY'. It means no")
    print("     AMPLITUDE-INDEPENDENT frequency, which this programme already")
    print("     measured from the other side as 'no linear sound speed'. The")
    print("     two are the same statement.")
    print("   * Vol_hull's CONE IS DISCARDED ON EVIDENCE. Its exponent is")
    print("     -1/2 and the medium's is -1. [22331] discarded it for the")
    print("     wrong reason and reached the right verdict.")
    print("   * NO DISPERSION RELATION, AND THE OBSTACLE IS NAMED. It needs a")
    print("     simultaneous multi-contact impact law -- DECISION 18's Phase 2")
    print("     LCP, specified and never built. R6 measures the obstruction")
    print("     rather than asserting it, so the next attempt knows what to")
    print("     build first.")
    print("   * NO V IN THE SMOOTH SENSE IS FOUND OR IMPLIED. This does not")
    print("     supply the potential the epic has been looking for; it")
    print("     measures that the medium's is the one shape that cannot be")
    print("     linearised, which is why the search kept failing.")
    print("   * The effective V is read off jb_ct's SCALAR radial reduction.")
    print("     jb_pr showed that understates the joint's coupling by sqrt(3).")
    print("     No exponent here changes: every one of them is a ratio.")
    print()
    print("  ALL CHECKS PASSED." if not bad else f"  {bad} CHECK(S) FAILED.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
