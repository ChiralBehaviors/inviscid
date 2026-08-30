"""phase_ramp -- play admits a bounded phase gradient, and the bound is the wave.

WHY THIS FILE EXISTS. A wave IS a phase gradient, and this programme has
already measured that a rigid chain forbids one. LineChainAnimation:36-41 and
:48-52:

    "every cell's fold angle must lie in [-60, 60], and a ramp pins p = 0, a
    single frozen configuration with no motion left" ... "There is no way to
    give the far cell a different phase from the near one while keeping the
    faces shared. The chain moves as a unit."

That is the rigid statement of exactly what a wave needs and cannot have.
jb_ct then measured that joint PLAY buys the medium a finite signal speed. The
obvious consequence had never been tested: with play, does a bounded ramp
become admissible? It does, and the bound is not the one the bead predicted.

THE INSTRUMENT, AND WHY IT IS NOT jb_ct's. jb_ct reduces a joint to ONE SCALAR
-- the radial separation V(g_k) + V(g_k+1) against a fixed site spacing. That
is a projection, and it is blind to the channel this question turns on. Two
neighbours sharing a triangular face must agree on where the face IS (radial)
AND on how it is TURNED about its own normal (twist). A phase gradient is
exactly a twist mismatch, so a radial-only model cannot see the constraint that
actually binds it.

So this file measures the joint the way the honeycomb actually holds it: the
worst discrepancy between the three corners of one cell's shared face and the
three corners of its neighbour's, resolved by permutation search, on jb_hc's
own placement. Reference joints close to 6e-16 in all eight directions and both
orientations (R1), so the instrument reproduces the structure before it is
asked anything.

THE RESULT, in the two coordinates the joint has:

    M(s, d) = C * |d| + C * s^2 + higher order

  * d is the DIFFERENTIAL phase (neighbours at 60 + d instead of 60). FIRST
    order. A ramp costs clearance immediately.
  * s is the COMMON phase (both cells shifted together). SECOND order, because
    a = -30 is where the lattice constant is stationary -- the same fact that
    makes the medium's coherent breathe cheap.
  * C is ONE coefficient for both, and it is the cell's own vertex speed
    |dv/dgamma| = sqrt(m(gamma)/8) with m jb_ct's measured inertia. At the
    reference it is EL/sqrt(2) = 1 exactly, and it is MAXIMAL there, so a = -30
    is where the medium is STIFFEST against a phase gradient.

Two bounds follow, and they are different orders of the same clearance:

    per neighbour   d_max = t / C          FIRST order in the play
    across a chain  s_max = sqrt(t / C)    SECOND order in the play

THE CORRECTION THIS FILE OWES THE RECORD. Bead qvf.29 and T2 [23639] predict
~0.70 degrees per neighbour for 1 mm of clearance on a 100 mm edge, from
t/Z with Z = EL sqrt(2/3). That is the RADIAL channel alone. The measured
coefficient is EL/sqrt(2), not EL sqrt(2/3), and the bound is 0.807 degrees
(0.810 at first order, approached from below) -- about 15% larger. The
coherent figure in the same record, ~6.8 degrees, is
CONFIRMED here to three digits. The "roughly 10x" between them is 8.41, and
that number is not a coincidence: s_max / d_max = sqrt(C/t) is also the chain
length at which the two bounds cross (R5).

AND THE CORRECTION IT OWES jb_ct. Same reduction, same reason: jb_ct's coupling
constant is |V'| = Z sin(30) = 0.5774 where the joint's real coefficient is
C = 1. Its speed law speed = |V'| gammadot / play has the RIGHT FORM and a
constant low by exactly sqrt(3). Every scaling in jb_ct's R5 and R6 is
unaffected -- they are ratios -- and its headline speed is not (R7).

WHAT A RAMP IS NOT. A monotone ramp accumulates common mode, so it saturates:
a long chain carries a SMALLER gradient, capped in total at s_max no matter how
many cells it is spread over (R5). A WAVE does not have that problem, because
its phase deviation oscillates instead of accumulating -- which is why R6
measures the admissible amplitude by wavenumber rather than quoting the ramp
bound as though it applied to a wave.

SCOPE, stated rather than discovered later.
  * ONE chain along one body diagonal. A bulk cell has eight triangular-face
    neighbours, not two, and nothing here measures how eight constraints on one
    cell interact.
  * The play is ONE ISOTROPIC clearance t compared against the worst matched-
    corner discrepancy. A real joint's clearance need not be isotropic, and if
    it is not, the binding combination is not this one.
  * V = 0 STILL. Nothing here finds a potential. A static ramp at rest is an
    equilibrium of this model TRIVIALLY -- inside the band there are no forces
    at all -- and that is not evidence of stability (R7).
  * The dynamics in R7 runs in jb_ct's SCALAR model, so its constant carries
    that sqrt(3). The signs and the monotonicity do not.
"""
from __future__ import annotations

import itertools as it
import sys

import numpy as np

from analysis.retired.hard_wall import contact_chain as CT
from analysis.model import plates as ZP
from analysis.retired.strut_springs import honeycomb_waves as HC
from analysis.model import assembly as RC
EL = CT.EL

#: The reference phase: midpoint of the exchange, widest lattice, stationary
#: point of the separation law, and -- measured in R1 -- the maximum of the
#: mismatch coefficient. jb_ct's uniform front lives here too.
A_REF = -30.0

#: The lattice spacing the chain is held at. FIXED, which is the whole model:
#: Gray p.40 says shared faces force the centres to move, the owner's build
#: holds station instead, and the price is paid in clearance.
L_REF = HC.lattice(A_REF)

#: One body diagonal. The chain runs along it, alternating VE and hole cells.
DIR = np.array([1.0, 1.0, 1.0])

#: The owner's rig, as a ratio: 1 mm of clearance on a 100 mm edge.
T_RIG = 0.01 * EL

#: Clearances swept by the gate, an eightfold range, matching jb_ct's.
T_GRID = (0.005, 0.01, 0.02, 0.04)


def mismatch(pa, pb, L=L_REF, d=DIR):
    """Worst corner discrepancy across the shared triangular face.

    Cell at phase `pa` sits at the origin; its neighbour at phase `pb` sits at
    `L*d`. The corner correspondence is resolved by permutation search, exactly
    as jb_hc's `face_pairing` does -- the face is an equilateral triangle, so
    which corner meets which is a bookkeeping question and not a geometric one,
    and choosing it by nearest match is what makes this a mismatch rather than
    a labelling artifact.
    """
    d = np.asarray(d, dtype=float)
    A = ZP.corners(pa)[HC.face_toward(d)]
    B = ZP.corners(pb)[HC.face_toward(-d)] + L * d
    return min(max(float(np.linalg.norm(A[pm[c]] - B[c])) for c in range(3))
               for pm in it.permutations(range(3)))


def coefficient(a0, h=1e-6):
    """dM/d(differential phase) at `a0`, per RADIAN. Central difference."""
    return mismatch(a0 - h / 2.0, a0 + 60.0 + h / 2.0,
                    HC.lattice(a0)) / np.radians(h)


def coefficient_closed(a0):
    """The same, from the cells' OWN analytic vertex velocities.

    A phase change moves each cell's corners at dv/dgamma, so the mismatch of a
    MATCHED corner pair grows at |dv_p/dgamma + dv_q/dgamma| / 2 when the two
    cells split the difference. The sum rather than the difference because the
    two sublattices fold in opposite senses across a shared face. This is an
    independent route -- RC's analytic derivative, not a difference of the
    placement -- so R1's agreement is a cross-check and not a tautology.
    """
    A = ZP.corners(a0)[HC.face_toward(DIR)]
    B = ZP.corners(a0 + 60.0)[HC.face_toward(-DIR)] + HC.lattice(a0) * DIR
    pm = min(it.permutations(range(3)),
             key=lambda q: sum(np.linalg.norm(A[q[c]] - B[c]) for c in range(3)))
    dp = RC.body(a0, 1)[1]
    dq = RC.body(a0 + 60.0, 1)[1]
    ip = [RC.SLOT[(HC.face_toward(DIR), c)] for c in range(3)]
    iq = [RC.SLOT[(HC.face_toward(-DIR), c)] for c in range(3)]
    return max(float(np.linalg.norm(dp[ip[pm[c]]] + dq[iq[c]])) / 2.0
               for c in range(3))


def chain_worst(dev, a0=A_REF, L=L_REF):
    """Worst joint mismatch along a chain whose cell k deviates by `dev[k]`
    degrees from its own sublattice reference."""
    g = [a0 + 60.0 * (k % 2) + dev[k] for k in range(len(dev))]
    return max(mismatch(g[k], g[k + 1], L) for k in range(len(dev) - 1))


def largest(pattern, t, hi, iters=48):
    """Largest amplitude of `pattern` (a callable amplitude -> deviations) that
    keeps every joint inside the play. Bisection; returns 0.0 if even an
    infinitesimal amplitude already violates."""
    if chain_worst(pattern(1e-12)) > t:
        return 0.0
    lo = 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if chain_worst(pattern(mid)) <= t:
            lo = mid
        else:
            hi = mid
    return lo


def ramp(n):
    return lambda d: d * np.arange(n)


def uniform(n):
    return lambda s: s * np.ones(n)


def wave(n, theta):
    return lambda A: A * np.cos(theta * np.arange(n))


def gate():
    checks, out = [], {}
    A = checks.append
    C = coefficient(A_REF)
    out["C"] = C

    # R1 -- the instrument, and the coefficient it exposes
    ref = max(max(mismatch(A_REF, A_REF + 60.0, L_REF, d),
                  mismatch(A_REF + 60.0, A_REF, L_REF, d)) for d in HC.DIRS)
    per_dir = [mismatch(A_REF - 5e-7, A_REF + 60.0 + 5e-7, L_REF, d) / np.radians(1e-6)
               for d in HC.DIRS]
    vspeed = [float(np.linalg.norm(RC.body(a, 1)[1], axis=1).max())
              for a in (-50.0, -45.0, -30.0, -10.0)]
    inert = [float(np.sqrt(CT.inertia(np.radians(a)) / 8.0))
             for a in (-50.0, -45.0, -30.0, -10.0)]
    A(("R1  THE JOINT IS NOT A SCALAR, AND ITS COEFFICIENT IS THE CELL'S OWN "
       "VERTEX SPEED. Two neighbours sharing a triangular face must agree on "
       "where the face sits AND on how it is turned about its normal, and a "
       "phase gradient is exactly a twist mismatch -- so jb_ct's radial-only "
       "reduction is blind to the channel this question turns on. Measured on "
       "jb_hc's own placement, the reference joints close to 6e-16 in all "
       "eight directions and both orientations, and the mismatch grows at "
       "EXACTLY EL/sqrt(2) = 1 per radian of differential phase, the same in "
       "every direction. That constant is not fitted: it is the matched "
       "corner's |dv/dgamma|, which is sqrt(m(gamma)/8) with m jb_ct's own "
       "measured inertia (16/3)(1 + 2 sin^2 g), and it agrees with RC's "
       "analytic derivative to 1e-7",
       ref < 1e-14 and max(abs(c - 1.0) for c in per_dir) < 1e-6
       and max(abs(v - i) for v, i in zip(vspeed, inert)) < 1e-12,
       f"reference joints worst {ref:.2e}; dM/dd per radian across the 8 "
       f"directions {min(per_dir):.7f}..{max(per_dir):.7f} against EL/sqrt(2) "
       f"= {EL / np.sqrt(2):.7f}; |dv/dg| vs sqrt(m/8) at g = -50,-45,-30,-10 "
       + ", ".join(f"{v:.7f}/{i:.7f}" for v, i in zip(vspeed, inert)),
       "joints exact, coefficient EL/sqrt(2) in every direction, and it IS "
       "the vertex speed"))

    # R2 -- and the coefficient is MAXIMAL at the reference
    span = {a: coefficient(a) for a in (-55.0, -45.0, -30.0, -15.0, -5.0)}
    closed = {a: coefficient_closed(a) for a in span}
    A(("R2  a = -30 IS WHERE THE MEDIUM IS STIFFEST AGAINST A PHASE GRADIENT. "
       "The coefficient is symmetric about the reference and MAXIMAL there, so "
       "the admissible gradient is at its SMALLEST exactly where jb_ct's front "
       "is uniform and where the lattice is widest -- three properties of one "
       "point, not three coincidences. Cross-checked against the closed form "
       "|dv_p/dgamma + dv_q/dgamma|/2 on the matched corners, which is an "
       "independent route through RC's analytic derivative rather than a "
       "difference of the same placement. TWO-SIDED: this row fails if the "
       "coefficient stops peaking at the reference OR if the two routes part",
       max(span, key=span.get) == -30.0
       and max(abs(span[a] - closed[a]) for a in span) < 1e-6
       and span[-55.0] < span[-30.0] and span[-5.0] < span[-30.0],
       "; ".join(f"a={a:+.0f}: {span[a]:.7f} (closed {closed[a]:.7f})"
                 for a in sorted(span)),
       "peak at a = -30, both routes agreeing to 1e-6"))

    # R3 -- the ramp bound, and the correction it forces on the record
    bound = {t: largest(ramp(2), t, 20.0) for t in T_GRID}
    lin = {t: np.degrees(t / C) for t in T_GRID}
    short = {t: 1.0 - bound[t] / lin[t] for t in T_GRID}
    rig = largest(ramp(2), T_RIG, 20.0)
    out["rig_ramp"] = rig
    old = np.degrees(T_RIG / (np.sqrt(2.0 / 3.0) * EL))
    A(("R3  PLAY ADMITS A RAMP, AT t/C PER NEIGHBOUR -- and the record's "
       "number for it is 15% low. A rigid chain forbids a phase gradient "
       "outright (LineChainAnimation: 'the chain moves as a unit'); give the "
       "joint clearance and a gradient becomes admissible. Bead qvf.29 and T2 "
       "[23639] predict t/Z with Z = EL sqrt(2/3), which is the RADIAL channel "
       "alone and gives 0.70 degrees per neighbour for 1 mm on a 100 mm edge. "
       "The joint's real coefficient is EL/sqrt(2), because the twist channel "
       "binds too, and the measured answer is 0.807. The bound approaches t/C "
       "FROM BELOW and never reaches it: even a single joint carries a common "
       "mode of d/2, whose second-order cost eats a share of the clearance "
       "proportional to t -- measured at 0.285*t across an eightfold range, "
       "against d/4 = 0.25*t from the two orders alone. TWO-SIDED: this row "
       "fails if the shortfall vanishes (the second order is real) OR if it "
       "stops being proportional to t (then it is not second order)",
       all(s > 0 for s in short.values())
       and max(short[t] / t for t in T_GRID)
       / min(short[t] / t for t in T_GRID) < 1.05
       and abs(rig - 0.807) < 0.002,
       "; ".join(f"t={t}: {bound[t]:.5f} deg, {100 * short[t]:.2f}% under "
                 f"t/C = {lin[t]:.5f}" for t in T_GRID)
       + "; shortfall/t = " + ", ".join(f"{short[t] / t:.3f}" for t in T_GRID)
       + f"; at the rig's 1 mm on 100 mm: {rig:.4f} deg per neighbour "
         f"(first order {np.degrees(T_RIG / C):.4f}), against the record's "
         f"t/Z = {old:.4f}",
       "shortfall positive, proportional to t within 5%, and 0.807 "
       "deg/neighbour at rig scale"))

    # R4 -- the coherent mode is second order, confirming the record
    coh = {t: largest(uniform(2), t, 40.0) for t in T_GRID}
    quad = {t: np.degrees(np.sqrt(t / C)) for t in T_GRID}
    rig_coh = largest(uniform(2), T_RIG, 40.0)
    out["rig_coh"] = rig_coh
    A(("R4  THE COHERENT MODE IS SECOND ORDER, AND THE RECORD'S 6.8 DEGREES IS "
       "CONFIRMED. Shifting both cells together costs nothing at first order, "
       "because a = -30 is the stationary point of the lattice constant -- the "
       "medium breathes for free and carries a disturbance grudgingly, which "
       "is the same asymmetry from the constraint side rather than the "
       "dynamical one. The bound goes as sqrt(t/C), so it is the SQUARE ROOT "
       "of the clearance where the ramp is LINEAR in it, and the two cannot be "
       "compared without saying which order each is",
       max(abs(coh[t] / quad[t] - 1.0) for t in T_GRID) < 0.02
       and abs(rig_coh - 6.81) < 0.05,
       "; ".join(f"t={t}: {coh[t]:.4f} deg against sqrt(t/C) = {quad[t]:.4f}"
                 for t in T_GRID)
       + f"; at rig scale {rig_coh:.3f} deg against the record's ~6.8",
       "sqrt(t/C) to 2%, and 6.81 deg at rig scale"))

    # R5 -- the ramp saturates, and the crossover IS the ratio
    ns = (4, 8, 14, 30, 60)
    sat = {n: largest(ramp(n), 0.02, 20.0) for n in ns}
    total = {n: (n - 1) * sat[n] for n in ns}
    cap = np.degrees(np.sqrt(0.02 / C))
    star = np.sqrt(C / T_RIG)
    out["sat"], out["cap"], out["star"] = sat, cap, star
    A(("R5  A MONOTONE RAMP SATURATES, AND THE CROSSOVER LENGTH IS THE RATIO "
       "BETWEEN THE TWO BOUNDS. A ramp accumulates COMMON mode as it goes, so "
       "the far joints pay the second-order cost even when the per-cell "
       "gradient is small: the total phase a chain can hold is capped at "
       "s_max no matter how many cells it is spread over, and a longer chain "
       "therefore carries a SMALLER gradient, falling as s_max/n. The two "
       "bounds cross at n* = s_max/d_max = sqrt(C/t), which is 8.4 cells at "
       "the rig's clearance -- so the record's 'roughly 10x' between coherent "
       "and differential range and the length at which a ramp stops being "
       "gradient-limited are THE SAME NUMBER. TWO-SIDED: this row fails if "
       "the per-cell bound ever exceeds the local one, or if the total stops "
       "converging on the cap",
       all(sat[a] > sat[b] for a, b in zip(ns, ns[1:]))
       and max(sat.values()) <= np.degrees(0.02 / C)
       and max(total.values()) <= cap
       and total[60] / cap > 0.95,
       "; ".join(f"n={n}: {sat[n]:.4f} deg/cell, total {total[n]:.3f}"
                 for n in ns)
       + f"; local bound {np.degrees(0.02 / C):.4f}, cap sqrt(t/C) = "
         f"{cap:.4f}, n* at rig scale {star:.2f} cells",
       "monotone, under both bounds, and within 5% of the cap by n = 60"))

    # R6 -- the admissible amplitude spectrum: what a WAVE may do
    thetas = (np.pi, np.pi / 2, np.pi / 4, np.pi / 8, np.pi / 16, np.pi / 32)
    spec = {}
    for th in thetas:
        f = np.cos(th * np.arange(64))
        g = float(np.abs(np.diff(f)).max())
        p = float(np.abs(0.5 * (f[:-1] + f[1:])).max())
        meas = largest(wave(64, th), 0.02, 60.0)
        pred = min(np.degrees(0.02 / (C * g)), np.degrees(np.sqrt(0.02 / C)) / p)
        spec[th] = (meas, pred, 2 * np.pi / th)
    out["spec"] = spec
    A(("R6  THE ADMISSIBLE AMPLITUDE SPECTRUM, WHICH IS WHAT A WAVE ACTUALLY "
       "HAS TO LIVE INSIDE. R5's saturation is a property of a MONOTONE ramp, "
       "not of a wave: a wave's phase deviation oscillates instead of "
       "accumulating, so it never runs up the common-mode cost. Imposing "
       "cos(k*theta) and bisecting the amplitude gives a two-branch envelope "
       "-- gradient-limited at short wavelength, amplitude-limited at long -- "
       "with the crossover near 2*pi/sqrt(t/C) cells. This is the closest "
       "thing to a dispersion relation this medium has produced, and it is "
       "NOT one: it says which amplitudes are ADMISSIBLE at each wavenumber, "
       "and with V = 0 there is no frequency to pair them with",
       max(abs(m / p - 1.0) for m, p, _ in spec.values()) < 0.05,
       "; ".join(f"lambda={lam:.0f} cells: {m:.4f} deg (envelope {p:.4f})"
                 for m, p, lam in spec.values()),
       "measured amplitude matches the two-branch envelope to 5%"))

    # R7 -- the two corrections to jb_ct, and the dynamic half of the bead
    pre = {}
    for r in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5):
        res = CT.chain(ramp=r)
        sp, spread, reach = CT.front(res)
        pre[r] = (sp, spread, res["drift"])
    quotable = [r for r in pre if pre[r][1] < 1.2]
    A(("R7  jb_ct'S COUPLING CONSTANT IS LOW BY sqrt(3), AND A CHAIN ALREADY "
       "CARRYING A GRADIENT TRANSMITS SLOWER. Two things, one cause. jb_ct's "
       "|V'| = Z sin(30) = 0.5774 is the radial projection of the coefficient "
       "measured here as 1, so its speed = |V'| gammadot / play has the right "
       "FORM with a constant low by exactly sqrt(3); every ratio in its R5 "
       "and R6 is untouched, its headline speed is not. And driving jb_ct's "
       "own integrator from a RAMPED start -- joints beginning part-way across "
       "their play rather than centred -- the front slows monotonically AND "
       "stops being uniform, so beyond the smallest ramps there is no speed "
       "to quote at all, by R4's own rule. Pre-stress does not stiffen this "
       "medium the way it would a granular chain; it detunes it. Energy holds "
       "at 1e-15 throughout, so this is the mechanism and not the integrator",
       abs(C / (CT.Z * 0.5) - np.sqrt(3.0)) < 1e-6
       and all(pre[a][0] > pre[b][0] for a, b in zip(sorted(pre),
                                                     sorted(pre)[1:]))
       and all(pre[a][1] < pre[b][1] for a, b in zip(sorted(pre),
                                                     sorted(pre)[1:]))
       and max(v[2] for v in pre.values()) < 1e-9,
       f"C/|V'| = {C / (CT.Z * 0.5):.7f} against sqrt(3) = {np.sqrt(3.0):.7f}; "
       + "; ".join(f"ramp {r:.1f}: speed {pre[r][0]:.3f}, spread {pre[r][1]:.2f}"
                   for r in sorted(pre))
       + f"; uniform (spread < 1.2) only at ramp <= {max(quotable):.1f}",
       "constant off by exactly sqrt(3); speed monotone down, spread "
       "monotone up, energy 1e-15"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_pr -- play admits a bounded phase gradient, and the bound is the wave")
    print("=" * 78)
    checks, out = gate()
    bad = 0
    for name, ok, got, want in checks:
        bad += 0 if ok else 1
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")

    print("\n  THE TWO BOUNDS AT THE OWNER'S SCALE (1 mm of clearance on a")
    print("  100 mm edge), both checkable by hand on the physical rig:")
    print(f"   * DIFFERENTIAL, first order:  {out['rig_ramp']:.3f} deg of phase "
          "per neighbour")
    print(f"   * COHERENT,     second order: {out['rig_coh']:.2f} deg, all cells "
          "together")
    print(f"   * they cross at n* = {out['star']:.1f} cells, which is also their "
          "ratio")
    print("   Looser joints carry more, linearly in the differential channel")
    print("   and as a square root in the coherent one -- so the two do NOT")
    print("   scale together, and a rig built to a different tolerance moves")
    print("   them apart.")
    print()
    print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
    print("   * A BOUNDED PHASE GRADIENT IS ADMISSIBLE. The lock")
    print("     LineChainAnimation measured is the play-free limit of this,")
    print("     exactly as jb_ic's instantaneous onset is the play-free limit")
    print("     of jb_ct's finite speed. Same clearance, twice.")
    print("   * The ramp bound is t/C with C = EL/sqrt(2), NOT t/Z. The")
    print("     record's 0.70 deg per neighbour is corrected to 0.807; its")
    print("     6.8 deg coherent figure is confirmed at 6.82.")
    print("   * A MONOTONE RAMP SATURATES; A WAVE DOES NOT. R5's cap is a")
    print("     property of accumulation, and R6 measures what a wave may")
    print("     actually do instead. Quoting R5's bound at a wave would")
    print("     understate it at every wavelength above the crossover.")
    print("   * NOT A DISPERSION RELATION. R6 gives admissible AMPLITUDE by")
    print("     wavenumber. With V = 0 there is no frequency to pair with it,")
    print("     and nothing here supplies one.")
    print("   * V = 0 STILL. A static ramp at rest is an equilibrium of this")
    print("     model TRIVIALLY -- inside the band there are no forces at all")
    print("     -- so 'the ramp persists' is a statement about the model")
    print("     having nothing to relax toward, not about stability.")
    print("   * ONE CHAIN, one body diagonal. A bulk cell has eight")
    print("     triangular-face neighbours and nothing here measures how")
    print("     eight constraints on one cell interact.")
    print("   * The play is ONE ISOTROPIC clearance against the worst matched")
    print("     corner. An anisotropic joint binds on a different combination.")
    print()
    print("  ALL CHECKS PASSED." if not bad else f"  {bad} CHECK(S) FAILED.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
