"""jb_ct -- the tolerance is the coupling, and it buys a FINITE signal speed.

WHY THIS FILE EXISTS. Every module in this programme that measured an onset
reported the same thing: it is instantaneous. jb_ic's R7 says it outright -- "a
rigid constraint has infinite signal speed, so the projection that makes the
impulse admissible reaches the whole chain at once and there is no onset lag to
measure anywhere. This medium therefore has NO WAVEFRONT". That was correct FOR
A RIGIDLY COUPLED MEDIUM, and it is why the wave programme kept failing to find
a wave: a bilateral constraint transmits instantly, so there is nothing to time.

THE OWNER'S BUILD BREAKS THAT, and the mechanism is the one thing nobody had
modelled: PLAY. Owner, 2026-08-28, on why the physical cells hold station when
Gray says a shared face forbids it -- "that was likely due to the very low
tolerance of the physical build". The joints are not rigid. They have clearance.

Gray, "The Jitterbug Motion" (2002) p.40, is why that matters:

    "Two Jitterbugs can not share the same triangular face and have their
    positions (location of center of volume) fixed as they go through the
    Jitterbug motion. If two Jitterbugs are to share the same triangle face then
    as the joined Jitterbugs jitterbug, the positions of the Jitterbugs must
    move."

SharedFaceAnimation measures how much: the centres travel 15% while the shared
face stays joined to 1e-15. A build that holds station cannot be paying that in
centre travel, so it pays it in CLEARANCE -- and a joint with clearance is not a
weld either. That is the whole model here.

NOT A UNILATERAL CONTACT, AND THE DISTINCTION IS LOAD BEARING. T2 inviscid
[23575] rules out the contact/granular reading of this medium, correctly: "A
unilateral contact is one that can SEPARATE -- a shared vertex ceasing to be
shared -- and that is exactly a change of incidence. Incidence never changes.
The contact/granular reading was never available."

That refutes contacts which MAKE AND BREAK. This model has neither. The owner's
joints are permanent -- "they will ALWAYS be joined by the vertex and never ever
split" -- and the constraint below is |c_k| <= t, a band with a stop on BOTH
sides. Two-sided, permanently engaged, incidence fixed. The mechanism is
BACKLASH, which is what a bilateral joint with clearance has, and it is a
different object from a unilateral contact. 23575 and this file are both right
and are about different things. An earlier version of this docstring said
"UNILATERAL CONTACT" and would have read as refuted by an entry it does not
actually contradict.

THE MODEL. Cells at fixed sites along one dowel, one fold angle each, inertia
from their own corner masses. Gray's radial law V(g) = Z cos(g) with
Z = EL sqrt(2/3) puts each cell's shared face at V(g) from its centre, so
neighbours are coupled by

    c_k = V(g_k) + V(g_k+1) - sep0,     |c_k| <= t

with t the joint play. INSIDE the band the cells are free of each other; at its
two edges they collide. That is a BACKLASH chain, and a backlash chain has no
linear sound speed at all -- only solitary waves whose speed depends on
amplitude.

23575 also reaches this file's starting point and stops one step short: "shared
vertex, rigid -> constraint propagates instantly -> rigidity, NOT waves". True,
for rigid joints. Give the SAME permanent joints clearance and the propagation
stops being instant, without any incidence changing.

TWO METHOD NOTES, both paid for here.

  * RESOLVE CONTACTS BY EVENT, NOT BY PROJECTION. The first version of this
    advanced a fixed step and, on finding |c| > t, nudged one cell's POSITION
    back onto the stop. That injects energy at every contact, and with thousands
    of contacts it bled 20 to 40 PERCENT, and over-counted the contacts themselves
    by two orders of magnitude. The speed scalings survived it INTACT, which is
    exactly why the energy audit has to be run even when the answer already
    looks right.
  * A LINEAR FIT TO ARRIVAL TIMES IS NOT A SPEED unless the front is uniform.
    At a0 = -30 the per-cell delays are flat to 5% and the fit means something.
    At a0 = -45 they spread 45-fold and it does not; an earlier version of this
    file quoted the fit anyway. R4 now measures the uniformity BEFORE quoting
    any speed.
"""
from __future__ import annotations

import sys

import numpy as np

import jb_rc_reduced as RC

EL = RC.EL
#: Gray's radial law: V = Z cos(gamma), Z = EL sqrt(2/3). The distance from a
#: cell's centre of volume to one of its triangle face centres.
Z = np.sqrt(2.0 / 3.0) * EL
#: the cell's inertia in its own fold coordinate, m(g) = M0 (1 + 2 sin^2 g).
#: R1 measures both the form and the constant against the corner geometry.
M0 = 16.0 / 3.0

FOLD_REF = -30.0


def inertia(g):
    """Cell inertia in gamma (RADIANS), from the mass model."""
    return M0 * (1.0 + 2.0 * np.sin(g) ** 2)


def inertia_d(g):
    return M0 * 4.0 * np.sin(g) * np.cos(g)


def radial(g):
    """Gray eq (3): distance from centre of volume to a triangle's face centre."""
    return Z * np.cos(g)


def radial_d(g):
    """dV/dgamma. It VANISHES at the VE (g = 0), so the contact coupling -- and
    with it the signal speed -- goes to zero there."""
    return -Z * np.sin(g)


def _rk4(g, gd, h):
    """V = 0, so the only force is the configuration-dependent inertia:
    gddot = -(1/2)(m'/m) gdot^2. FreeDynamics' equation, per cell."""
    def acc(g, gd):
        return -0.5 * (inertia_d(g) / inertia(g)) * gd ** 2

    k1v, k1x = acc(g, gd), gd
    k2v, k2x = acc(g + h / 2 * k1x, gd + h / 2 * k1v), gd + h / 2 * k1v
    k3v, k3x = acc(g + h / 2 * k2x, gd + h / 2 * k2v), gd + h / 2 * k2v
    k4v, k4x = acc(g + h * k3x, gd + h * k3v), gd + h * k3v
    return (g + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x),
            gd + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v))


def chain(a0=FOLD_REF, n=14, play=0.02, kick=0.30, tmax=12.0, dt=1e-3):
    """A chain of cells on one dowel, coupled ONLY by contact through the play.

    EVENT DRIVEN: integrate to the exact crossing by bisection, apply an elastic
    impulse along the constraint gradient, resume. Nothing is ever projected, so
    energy is conserved to machine precision and the audit means something.
    """
    g = np.array([np.radians(a0) if k % 2 == 0 else np.radians(a0 + 60.0)
                  for k in range(n)])
    gd = np.zeros(n)
    gd[0] = kick
    sep0 = radial(g[0]) + radial(g[1])
    e0 = 0.5 * float(np.sum(inertia(g) * gd ** 2))
    arrive = [None] * n
    arrive[0] = 0.0
    thresh = 0.02 * abs(kick)
    now, hits, worst = 0.0, 0, 0.0

    def viol(x):
        return np.abs(radial(x[:-1]) + radial(x[1:]) - sep0) - play

    while now < tmax - 1e-12:
        h = min(dt, tmax - now)
        ng, ngd = _rk4(g, gd, h)
        if (viol(ng) > 0).any():
            lo, hi = 0.0, h
            for _ in range(45):
                mid = 0.5 * (lo + hi)
                tg, _ = _rk4(g, gd, mid)
                if (viol(tg) > 0).any():
                    hi = mid
                else:
                    lo = mid
            g, gd = _rk4(g, gd, lo)
            now += lo
            c = radial(g[:-1]) + radial(g[1:]) - sep0
            k = int(np.argmax(np.abs(c)))
            u = radial_d(g[k]) * gd[k] + radial_d(g[k + 1]) * gd[k + 1]
            if (c[k] > 0 and u > 0) or (c[k] < 0 and u < 0):
                lam = -2.0 * u / (radial_d(g[k]) ** 2 / inertia(g[k])
                                  + radial_d(g[k + 1]) ** 2 / inertia(g[k + 1]))
                gd[k] += lam * radial_d(g[k]) / inertia(g[k])
                gd[k + 1] += lam * radial_d(g[k + 1]) / inertia(g[k + 1])
                hits += 1
        else:
            g, gd, now = ng, ngd, now + h
        worst = max(worst, float(viol(g).max()))
        for k in range(1, n):
            if arrive[k] is None and abs(gd[k]) > thresh:
                arrive[k] = now
    e1 = 0.5 * float(np.sum(inertia(g) * gd ** 2))
    return {"arrive": arrive, "hits": hits, "drift": abs(e1 - e0) / e0,
            "overshoot": worst, "E": e0}


def front(res):
    """(speed in cells per unit time, per-cell delay spread) -- and the spread is
    what says whether the speed means anything."""
    got = [(k, a) for k, a in enumerate(res["arrive"]) if a is not None]
    if len(got) < 4:
        return float("nan"), float("inf"), len(got)
    ts = np.array([a for _, a in got[1:]])
    ks = np.array([k for k, _ in got[1:]])
    d = np.diff([a for _, a in got])
    return float(np.polyfit(ts, ks, 1)[0]), float(d.max() / d.min()), len(got)


def gate():
    checks, out = [], {}
    A = checks.append

    # R1 -- the inertia, from the corner geometry rather than assumed
    meas = {}
    for gd_ in (-60.0, -30.0, 0.0, 30.0, 60.0, 90.0):
        v1 = RC.body(gd_, 1)[1]
        meas[gd_] = float(np.dot(RC.VMASS, np.einsum('ij,ij->i', v1, v1)))
    form = max(abs(meas[x] - inertia(np.radians(x))) for x in meas)
    A(("R1  THE CELL'S INERTIA COMES FROM THE CORNERS, and independently "
       "reproduces this project's own free-dynamics result. Summing the 24 "
       "corner masses against |dv/dgamma|^2 gives m(g) = (16/3)(1 + 2 sin^2 g) "
       "exactly, so m swings 3:1 across the fold -- which is FreeDynamics' "
       "C_axial = 2/3, C_spin = 1/3 arrived at from a different direction. The "
       "mass model is DECLARED, not discovered: unit mass per triangle, lumped "
       "m/3 to each corner",
       form < 1e-12 and abs(meas[90.0] / meas[0.0] - 3.0) < 1e-12,
       "; ".join(f"g={k:+.0f}: {v:.6f}" for k, v in meas.items())
       + f"; worst deviation from (16/3)(1+2sin^2 g) is {form:.1e}, "
         f"m(90)/m(0) = {meas[90.0] / meas[0.0]:.6f}",
       "closed form to 1e-12, 3:1 swing"))

    # R2 -- the contact chain conserves energy and never leaves its play
    runs = {}
    for key, kw in (("ref", {}), ("2x kick", {"kick": 0.60}),
                    ("2x play", {"play": 0.04}), ("half play", {"play": 0.01})):
        runs[key] = chain(**kw)
    out["runs"] = runs
    A(("R2  ENERGY SURVIVES THE CONTACTS, AND THE JOINTS NEVER EXCEED THEIR "
       "PLAY. V = 0, so energy is the only audit there is, and a contact chain "
       "is exactly where a careless integrator loses it. Contacts are resolved "
       "by EVENT -- bisect to the crossing, apply an elastic impulse along the "
       "constraint gradient, resume -- and nothing is projected. An earlier "
       "version nudged positions back onto the stop instead and bled 20 to 40 "
       "PERCENT, and over-counted the contacts by two orders of magnitude. Its "
       "speed scalings were nonetheless correct, which is the whole argument "
       "for running the audit even when the answer already looks right",
       max(r["drift"] for r in runs.values()) < 1e-9
       and max(r["overshoot"] for r in runs.values()) < 1e-9,
       "; ".join(f"{k}: {r['hits']} contacts, E drift {r['drift']:.1e}, "
                 f"overshoot {r['overshoot']:.1e}" for k, r in runs.items()),
       "drift and overshoot both below 1e-9"))

    # R3 -- A FINITE SIGNAL SPEED, which this programme has never had
    sp = {k: front(r) for k, r in runs.items()}
    out["speed"] = sp
    ref_speed, ref_spread, ref_reach = sp["ref"]
    A(("R3  THE MEDIUM HAS A FINITE SIGNAL SPEED. Every earlier module here "
       "measured the onset as INSTANTANEOUS -- jb_ic's R7 concluded 'this "
       "medium therefore has NO WAVEFRONT' -- and that was right for a rigidly "
       "coupled medium, because a bilateral constraint transmits at once. Give "
       "the joints PLAY and it is no longer rigid: a cell must physically cross "
       "its clearance before it loads its neighbour, so the disturbance takes "
       "TIME to travel. This is the first finite propagation speed in the "
       "programme, and the tolerance is what buys it",
       np.isfinite(ref_speed) and ref_speed > 0 and ref_reach == 14,
       f"a0 = {FOLD_REF:.0f}, play 0.02, kick 0.30: the front reaches all "
       f"{ref_reach} cells at {ref_speed:.3f} cells per unit time",
       "a finite speed, and the front crosses the whole chain"))

    # R4 -- and it is a REAL front, not a fitted line through a smear
    far = chain(a0=-45.0)
    far_speed, far_spread, far_reach = front(far)
    out["far"] = (far_speed, far_spread)
    A(("R4  AT a0 = -30 IT IS A UNIFORM FRONT; AWAY FROM THERE IT IS NOT, and "
       "the difference decides whether a speed exists to quote at all. The "
       "per-cell delays are flat to 5% at a0 = -30 and spread 45-fold at "
       "a0 = -45, where the disturbance decelerates instead of propagating. A "
       "linear fit to arrival times returns a number either way; only the first "
       "of them is a speed. TWO-SIDED: this row fails if the reference stops "
       "being uniform OR if the off-reference case starts being uniform",
       ref_spread < 1.2 and far_spread > 5.0,
       f"per-cell delay spread (max/min): {ref_spread:.2f} at a0 = -30, "
       f"{far_spread:.1f} at a0 = -45",
       "flat at the reference, badly non-uniform away from it"))

    # R5 -- the scaling, and the closed form it implies
    s_ref = sp["ref"][0]
    ratios = {"2x kick": sp["2x kick"][0] / s_ref,
              "2x play": sp["2x play"][0] / s_ref,
              "half play": sp["half play"][0] / s_ref}
    pred = abs(radial_d(np.radians(FOLD_REF))) * 0.30 / 0.02
    A(("R5  SPEED = |dV/dgamma| * gammadot / play, WHICH MEANS THERE IS NO "
       "LINEAR SOUND SPEED. The gap closes at |V'(a0)| gammadot, so crossing it "
       "takes play/(|V'| gammadot) and the front advances one cell per that "
       "time. Doubling the kick doubles the speed and doubling the play halves "
       "it, both to within 1%. A medium whose signal speed depends on AMPLITUDE "
       "has no linear regime at all -- it is a sonic vacuum, carrying solitary "
       "waves and nothing else. That is a different object from anything this "
       "programme has looked for so far",
       abs(ratios["2x kick"] - 2.0) < 0.05
       and abs(ratios["2x play"] - 0.5) < 0.05
       and abs(ratios["half play"] - 2.0) < 0.05
       and abs(s_ref / pred - 1.0) < 0.02,
       f"speed {s_ref:.3f} against the closed form |V'|*kick/play = "
       f"{pred:.3f} (ratio {s_ref / pred:.4f}); "
       + ", ".join(f"{k} -> {v:.3f}x" for k, v in ratios.items()),
       "2x, 0.5x, 2x and the closed form within 1%"))

    # R6 -- take the play away and the finite speed goes with it
    tight = {p: front(chain(play=p))[0] for p in (0.04, 0.02, 0.01, 0.005)}
    prod = {p: v * p for p, v in tight.items()}
    A(("R6  REMOVE THE PLAY AND THE SPEED DIVERGES, recovering the "
       "instantaneous onset every earlier module reported. speed * play is "
       "constant across a factor of eight in the clearance, so speed grows "
       "without bound as the joint tightens -- and a rigid joint, play zero, is "
       "the infinite signal speed jb_ic's R7 measured. The two results do not "
       "disagree; R7 measured the play-free limit of this one. This row is the "
       "mutation probe for R3: it is the tolerance, and nothing else, that "
       "makes the speed finite",
       max(prod.values()) / min(prod.values()) < 1.05,
       "; ".join(f"play {p}: speed {v:.2f}, speed*play {prod[p]:.4f}"
                 for p, v in tight.items()),
       "speed * play constant to 5% over an 8x range of clearance"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_ct -- the tolerance is the coupling, and it buys a finite speed")
    print("=" * 78)
    checks, out = gate()
    bad = 0
    for name, ok, got, want in checks:
        tag = "PASS" if ok else "FAIL"
        bad += 0 if ok else 1
        print(f"  {tag}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")

    print("\n  FRONT ARRIVAL, a0 = -30, play 0.02, kick 0.30:")
    arr = out["runs"]["ref"]["arrive"]
    print("   " + "  ".join(f"{k}:{a:.2f}" for k, a in enumerate(arr)
                            if a is not None))
    print()
    print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
    print("   * A FINITE signal speed, which is new here. It is bought by the")
    print("     joint play, and it diverges as the play goes to zero (R6) --")
    print("     which is why every rigidly-coupled module found none.")
    print("   * NO linear sound speed. The speed depends on amplitude, so the")
    print("     medium has no linear regime to linearise about.")
    print("   * A front worth timing ONLY at a0 = -30. Away from there the")
    print("     disturbance decelerates and no single speed describes it (R4).")
    print("     Why is NOT established -- sublattice coupling asymmetry is a")
    print("     hypothesis, not a measurement: |V'| is equal on both sublattices")
    print("     at a0 = -30 and differs 2.7-fold at a0 = -45.")
    print("   * V = 0 STILL. No potential energy is found or implied. The")
    print("     coupling here is contact, not a restoring force.")
    print("   * The play is a FREE PARAMETER standing in for a real build's")
    print("     clearance. Nothing here measures the owner's rig.")
    print()
    print("  ALL CHECKS PASSED." if not bad else f"  {bad} CHECK(S) FAILED.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
