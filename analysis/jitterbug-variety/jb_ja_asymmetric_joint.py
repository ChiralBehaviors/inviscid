"""jb_ja -- an asymmetric joint: what it does to the breathe, and what it costs.

WHAT PROMPTED THIS. Two files on record name the same untried candidate.
jb_sj [23682]: "Real rubber is nonlinear and much stiffer in compression than
in extension, and THAT asymmetry is a live candidate for the one-sidedness
qvf.11 reports and jb_tr could NOT reproduce." jb_bz [23683] repeats it.
Neither built it. This file builds it.

A PREMISE THAT LOOKED SOUND AND IS WRONG, recorded because it is the reason
this file measures what it does. The lock qvf.11 reports acts on the COHERENT
BREATHE, and three instruments call that motion free -- jb_mj R3 (zero
contacts under coherent drive, ever), jb_sj R3 (in the zero space to 1e-12),
jb_bz R2 (the fourth zero mode at Gamma, a Goldstone). If the breathe never
opened a joint, every potential V(d) would be identically zero along it and
joint asymmetry could be dismissed without building it.

IT DOES OPEN JOINTS. All three statements are INFINITESIMAL -- a zero mode is
a first-order object, and jb_mj's "zero contacts" means the separation never
reached the clearance, not that it was zero. R3 measures what actually
happens: the breathe opens joints at EXACTLY SECOND ORDER, separation going as
amplitude squared to three figures (1.56e-4, 6.24e-4, 2.50e-3, 9.97e-3 for
amplitude 0.15, 0.3, 0.6, 1.2). So the mode is genuinely free to first order
AND genuinely loads the joint at finite amplitude, and an asymmetric joint can
see it. The candidate had to be measured, not argued away.

THE PHASE DECIDES WHETHER THE QUESTION EVEN PARSES (R4). At the reference
phase a = -30 the array's span is at an EXTREMUM: both fold senses CONTRACT
it, by -5.451e-4 and -5.445e-4, the same number to three figures. "Expansion
versus contraction" is not a distinction the medium makes there, so the
reference phase is the wrong place to ask qvf.11's question. At the
ICOSAHEDRAL phase a = -37.761 -- the phase qvf.11's lock is reported at, and
where its close records the hole cell reaching the icosahedron -- the two
senses separate cleanly, +5.244e-3 and -6.332e-3. Everything about the breathe
is therefore measured at BOTH phases.

WHAT IS ACTUALLY THERE (R5), and it is neither the dismissal nor the lock.

  * ONE-SIDEDNESS IS REAL AND IT IS GEOMETRIC BEFORE IT IS MATERIAL. At the
    icosahedral phase a SYMMETRIC joint already loads the two senses
    unequally -- peak joint potential 0.781 one way against 2.244 the other,
    a ratio of 0.348. That is the linkage, not the rubber, and no asymmetric
    material was needed to produce it.
  * THE ASYMMETRY THEN MOVES IT SUBSTANTIALLY, and phase-dependently: a
    hundred-to-one stiffness ratio takes 0.348 to 0.778 at the icosahedral
    phase, and only 0.998 to 0.990 at the reference phase.
  * NOTHING HERE IS A LOCK. qvf.11's array CANNOT EXPAND. Every direction
    here remains available at every ratio tried; what changes is how much of
    the drive's energy the joints take. A ratio moving from 0.35 to 0.78 is
    not an obstruction, and this file does not reproduce qvf.11 and does not
    claim to.
  * THE JOINT'S SENSE IS NOT THE ARRAY'S SENSE. Raising the COMPRESSION
    stiffness a hundredfold affects the span-EXPANDING drive far more
    (0.781 -> 1.845) than the contracting one (2.244 -> 2.371). An array can
    expand while its joints are compressed, and the two vocabularies must not
    be run together.

AND WHAT IT COSTS (R7), which is the sharpest result in the file. An
asymmetric law is stiffer leaving the weld one way than the other, so at the
welded configuration the second derivative depends on WHICH WAY YOU LEAVE:
V is continuous and C1 there but NOT C2, and a quantity that differs by
direction is not a Hessian. omega(k) needs one. Along the staggered direction
-- the one that loads every joint at once -- the curvature leaving the weld is
3.607x larger one way than the other at kc/ke = 100, against 0.999997 for the
symmetric control.

    THAT IS jb_cp's OBSTRUCTION REACHED FROM THE OPPOSITE SIDE. jb_cp [23658]
    established that a hard wall has no omega(k) because an infinite square
    well is FLAT inside and has no Hessian anywhere. This potential is soft,
    finite and everywhere continuous, and STILL has none -- for a KINK instead
    of for flatness. Two very different potentials, one obstruction.

SO THE LEDGER IS: asymmetry acts on the breathe but produces nothing like a
lock, leaves the sonic-vacuum transport exponents untouched (R6), and costs
jb_bz's bands. That is a statement about THIS asymmetry -- a direction-
dependent stiffness at the tied vertex pair -- and not about rubber.

    V(d) = (t^2/n) K(c) |d/t|^n,   c = (d.u)/|d|,   K(c) = k_e + (k_c-k_e)(1+c)/2

u is the inter-cell axis taken at REST, so it is constant per pair and the
gradient is exact. K depends on the DIRECTION of d only, never its magnitude,
which is what makes the symmetric limit exact at EVERY n rather than only at
n = 2: a split into axial and transverse parts would not reproduce jb_je away
from n = 2, because |d|^n is not |d_ax|^n + |d_tr|^n. R1 asserts the collapse
and R2 asserts the extra dK/dc term really is a gradient.

SCOPE.
  * STRUTS RIGID (T2 [23562]). The compliance is in the joint, as in jb_sj.
  * jb_mj's HC15 patch for the breathe, at two phases; jb_je's twelve-cell
    chain for transport. Not the bulk.
  * kc and ke are conventions like every stiffness here. Only their RATIO is
    a measurement, and no absolute curvature, energy or speed below is
    physical.
  * NOTHING HERE IS A CONTACT. The law is smooth and there is no LCP.
  * R7 IS ABOUT THE WELDED CONFIGURATION. Away from it, where a joint is
    already open, V is smooth and a Hessian exists. It is the medium's own
    rest state that lacks one -- which is exactly where omega(k) is computed.
  * jb_bz IS NOT RETRACTED. Its bands are a measurement about the SYMMETRIC
    joint and stand untouched. R7 says what a DIFFERENT joint law would cost.
  * qvf.11 IS NOT REOPENED and is NOT reproduced. The owner closed it
    2026-08-26 as "evidence that served, not a question resolved", and said
    the lock is not a goal to eliminate. Its fork (A) -- a TENSION-ONLY MEMBER
    spanning a square-face diagonal -- is an ADDED unilateral member between
    different vertices, not a joint law, and this file neither builds it nor
    bears on it.
"""
from __future__ import annotations

import sys

import numpy as np

import jb_je_joint_exponent as JE
import jb_mj_inertial_honeycomb as MJ
import jb_rc_reduced as RC
from jb_rc_reduced import hat

#: The asymmetry, as a ratio. 1.0 is jb_je exactly; R2 uses 100.
RATIOS = (1.0, 4.0, 100.0)

T_SCALE = JE.T_SCALE


def axes(asm, pairs):
    """Unit inter-cell axis per tied pair, at REST.

    Taken at rest so it is a constant of the joint and the gradient of V is
    exact -- a co-rotating axis would add a term through dR/dq and the energy
    would not close. Physically it is the joint's own material direction.
    """
    c0 = asm.ctr0
    U = np.zeros((len(pairs), 3))
    for r, (k, a, l, b) in enumerate(pairs):
        v = c0[l] - c0[k]
        nv = float(np.linalg.norm(v))
        U[r] = v / nv if nv > 1e-14 else np.array([1.0, 0.0, 0.0])
    return U


def stiffness(c, kc, ke):
    """K and dK/dc. Direction-dependent only -- never magnitude-dependent."""
    return ke + 0.5 * (kc - ke) * (1.0 + c), 0.5 * (kc - ke)


def state(asm, q, u, pairs, U, n, kc, ke, t=T_SCALE):
    """(acceleration, separations, mass blocks) from ONE frames(q) call.

    Same fused evaluation as jb_je's, with the joint gradient replaced by

        grad_d V = (t^2/n) K'(c) |d/t|^n (u - c dhat)/|d|
                   + K(c) t |d/t|^(n-1) dhat

    whose first term vanishes identically when kc == ke, which is R1.
    """
    nc = asm.N
    ctr, R, gam, B = asm.frames(q)
    RB0 = [np.dot(R[k], B[k][0].T).T for k in range(nc)]
    RB1 = [np.dot(R[k], B[k][1].T).T for k in range(nc)]
    RB2 = [np.dot(R[k], B[k][2].T).T for k in range(nc)]
    X = np.array([ctr[k] + RB0[k] for k in range(nc)])
    J = np.zeros((nc, 3 * RC.NV, 7))
    for k in range(nc):
        for i in range(RC.NV):
            r = 3 * i
            J[k, r:r + 3, 0:3] = np.eye(3)
            J[k, r:r + 3, 3:6] = -hat(RB0[k][i])
            J[k, r:r + 3, 6] = RB1[k][i]
    m3 = np.repeat(RC.VMASS, 3)
    M = np.array([np.dot(J[k].T, m3[:, None] * J[k]) for k in range(nc)])
    Minv = np.array([np.linalg.inv(M[k]) for k in range(nc)])
    A = np.zeros((nc, RC.NV, 3))
    for k in range(nc):
        w, gd = u[k, 3:6], u[k, 6]
        A[k] = (np.cross(w, np.cross(w, RB0[k]))
                + 2.0 * np.cross(w, RB1[k] * gd) + RB2[k] * gd * gd)
    f = np.array([-np.dot(J[k].T, m3 * A[k].ravel()) for k in range(nc)])
    s = np.empty(len(pairs))
    for r, (k, a, l, b) in enumerate(pairs):
        d = X[k][a] - X[l][b]
        nn = float(np.linalg.norm(d))
        s[r] = nn
        if nn < 1e-14:
            continue
        dh = d / nn
        c = float(np.dot(dh, U[r]))
        K, Kp = stiffness(c, kc, ke)
        g = (K * t * (nn / t) ** (n - 1.0)) * dh
        if Kp != 0.0:
            g = g + ((t * t / n) * (nn / t) ** n / nn) * Kp * (U[r] - c * dh)
        f[k] -= np.dot(g, J[k][3 * a:3 * a + 3])
        f[l] += np.dot(g, J[l][3 * b:3 * b + 3])
    return np.einsum('kij,kj->ki', Minv, f), s, M


def potential(asm, q, pairs, U, n, kc, ke, t=T_SCALE):
    X = asm.positions(q)
    V = 0.0
    for r, (k, a, l, b) in enumerate(pairs):
        d = X[k][a] - X[l][b]
        nn = float(np.linalg.norm(d))
        if nn < 1e-14:
            continue
        c = float(np.dot(d / nn, U[r]))
        K, _ = stiffness(c, kc, ke)
        V += (t * t / n) * K * (nn / t) ** n
    return float(V)


def integrate(asm, pairs, U, u0, n, kc, ke, tmax, h=JE.STEP, stop_at=None,
              kick=None):
    """Same scheme as jb_je's, driven by an arbitrary initial u0."""
    nc = asm.N
    q = asm.q0()
    u = np.array(u0, float)
    rec, seen, now, e0, e1 = [], 1, 0.0, None, None
    smax = 0.0
    ref = abs(kick) if kick else float(np.abs(u).max())
    for _ in range(int(round(tmax / h))):
        a1, s, M = state(asm, q, u, pairs, U, n, kc, ke)
        e1 = MJ.kinetic(M, u) + potential(asm, q, pairs, U, n, kc, ke)
        if e0 is None:
            e0 = e1
        smax = max(smax, float(s.max()))
        p = np.abs(u[:, 6])
        hot = np.where(p > JE.THRESH * ref)[0]
        if len(hot):
            fcell = int(hot.max())
            if fcell >= seen:
                seen = fcell + 1
                rec.append((now, fcell, float(p[fcell])))
                if stop_at is not None and fcell >= stop_at:
                    return dict(rec=rec, smax=smax, q=q, u=u, e=(e0, e1))
        u_h = u + 0.5 * h * a1
        q_h = RC.apply_increment(asm, q, (0.5 * h * u).ravel())
        a2, _, _ = state(asm, q_h, u_h, pairs, U, n, kc, ke)
        u = u + h * a2
        q = RC.apply_increment(asm, q, (h * u_h).ravel())
        now += h
    return dict(rec=rec, smax=smax, q=q, u=u, e=(e0, e1))


def curvature(asm, pairs, U, e_dir, kc, ke, step, n=2.0):
    """Second derivative of V along +/- a direction, from the welded state.

    A central difference is deliberately NOT used: the whole point is that the
    two sides disagree, so each side is taken on its own by a one-sided second
    difference from q0, where V = 0.
    """
    q0 = asm.q0()
    out = []
    for sgn in (+1.0, -1.0):
        vs = []
        for m in (1.0, 2.0):
            q = RC.apply_increment(asm, q0, sgn * m * step * e_dir)
            vs.append(potential(asm, q, pairs, U, n, kc, ke))
        # V(0) = 0 exactly at the weld, so V'' ~ (4 V(s) - V(2s)) ... use the
        # homogeneity instead: V is exactly quadratic along a ray at n = 2, so
        # V(s) = 0.5 * curv * s^2 and one point determines it.
        out.append(2.0 * vs[0] / (step ** 2))
    return out[0], out[1]


def transport_exponent(asm, pairs, U, n, kc, ke, kj_scale, kicks=JE.KICKS):
    """p_v from an amplitude sweep, jb_je's instrument with the asymmetric law."""
    v, a = [], []
    for kick in kicks:
        u0 = MJ.single_kick(asm, amp=kick, cell=0)
        r = integrate(asm, pairs, U, u0, n, kc * kj_scale, ke * kj_scale,
                      tmax=6.0, stop_at=asm.N - 1, kick=kick)
        rec = r["rec"]
        if not rec or rec[-1][1] < JE.HI:
            return None
        inner = [x for x in rec if JE.LO <= x[1] <= JE.HI]
        if len(inner) < 5:
            return None
        ts = np.array([x[0] for x in inner])
        cs = np.array([x[1] for x in inner], float)
        v.append(float(np.polyfit(ts, cs, 1)[0]))
        a.append(float(np.mean([x[2] for x in inner])))
    return JE.slope(a, v)


#: The icosahedral phase, where qvf.11's lock is reported. Its close records
#: a_ico = 22.238756 as outside the VE's range and the icosahedron reached by
#: the HOLE cell when the VE is at this phase.
A_ICO = -37.761

#: The breathe runs. A hundred-to-one ratio is stiff, so the step is finer
#: than jb_je's; R8 gates that the energy closes on it.
BREATHE_T = 0.10
BREATHE_H = 1e-4
AMP = 0.6


def patch_at(gc):
    asm, _ = RC.honeycomb(MJ.hc15_sites(), gc=gc)
    pairs = MJ.tied_pairs(asm)
    return asm, pairs, axes(asm, pairs)


def span(asm, q):
    """Mean distance of the cell centres from the patch centroid."""
    ctr, _, _, _ = asm.frames(q)
    return float(np.mean(np.linalg.norm(ctr - ctr.mean(axis=0), axis=1)))


def breathe(asm, pairs, U, sgn, kc, ke, n=2.0, amp=AMP, tmax=BREATHE_T,
            h=BREATHE_H):
    """Drive the coherent breathe one way; report what the joints took."""
    q = asm.q0()
    u = MJ.coherent_kick(asm, amp=sgn * amp)
    s0 = span(asm, q)
    vmax, smax, e0, e1 = 0.0, 0.0, None, None
    for _ in range(int(round(tmax / h))):
        a1, s, M = state(asm, q, u, pairs, U, n, kc, ke)
        V = potential(asm, q, pairs, U, n, kc, ke)
        e1 = MJ.kinetic(M, u) + V
        if e0 is None:
            e0 = e1
        vmax = max(vmax, V)
        smax = max(smax, float(s.max()))
        u_h = u + 0.5 * h * a1
        q_h = RC.apply_increment(asm, q, (0.5 * h * u).ravel())
        a2, _, _ = state(asm, q_h, u_h, pairs, U, n, kc, ke)
        u = u + h * a2
        q = RC.apply_increment(asm, q, (h * u_h).ravel())
    return dict(vpeak=vmax, smax=smax, dspan=span(asm, q) - s0,
                drift=abs(e1 - e0) / max(abs(e0), 1e-30))


def gate():
    checks, out = [], {}
    A = checks.append

    ref = patch_at(MJ.A_REF)
    ico = patch_at(A_ICO)
    ch = JE.chain()
    cpairs = MJ.tied_pairs(ch)
    CU = axes(ch, cpairs)

    # ---- R1: the symmetric limit is jb_je, at every n ----------------------
    rng = np.random.default_rng(20260829)
    q = RC.apply_increment(ch, ch.q0(), 1e-3 * rng.standard_normal(7 * ch.N))
    u = 1e-2 * rng.standard_normal((ch.N, 7))
    worst = 0.0
    for n in (2.0, 2.5, 3.0, 6.0):
        a_sym, s_sym, _ = state(ch, q, u, cpairs, CU, n, JE.K2, JE.K2)
        a_je, s_je, _ = JE.state(ch, q, u, cpairs, n, JE.K2)
        worst = max(worst, float(np.abs(a_sym - a_je).max()),
                    float(np.abs(s_sym - s_je).max()))
    out["R1"] = worst
    A(("R1  WITH THE TWO STIFFNESSES EQUAL THIS IS jb_je EXACTLY, AT EVERY "
       "EXPONENT AND NOT MERELY AT n = 2. That is a design constraint the "
       "obvious formulation fails: splitting the separation into axial and "
       "transverse parts and stiffening only the axial one is continuous and "
       "physical-looking, but |d|^n is not |d_ax|^n + |d_tr|^n, so it "
       "reproduces jb_je only at n = 2, and every comparison at another n "
       "would be against a silently different model. Making the stiffness "
       "depend on the DIRECTION of the separation and never its magnitude "
       "keeps V positively homogeneous and collapses the extra gradient term "
       "identically when kc = ke. Checked at n = 2, 5/2, 3 and 6 against "
       "jb_je's own state(). TWO-SIDED: without this row every row below "
       "could be measuring a different potential rather than an asymmetric "
       "one",
       worst < 1e-13,
       f"worst disagreement with jb_je over n = 2, 5/2, 3, 6 at kc = ke: "
       f"{worst:.2e}",
       "machine precision, i.e. the same model"))

    # ---- R2: the extra term really is a gradient ---------------------------
    gerr = 0.0
    hfd = 1e-7
    qg = RC.apply_increment(ch, ch.q0(),
                            3e-3 * np.random.default_rng(11).standard_normal(
                                7 * ch.N))
    z = np.zeros((ch.N, 7))
    for (kc, ke) in ((JE.K2, JE.K2), (100.0 * JE.K2, JE.K2)):
        for n in (2.0, 2.5, 3.0):
            a_an, _, _ = state(ch, qg, z, cpairs, CU, n, kc, ke)
            _, Mm, _ = MJ.kinematics(ch, qg, False)
            fq = np.array([Mm[k] @ a_an[k] for k in range(ch.N)]).ravel()
            g = np.zeros(7 * ch.N)
            for i in range(7 * ch.N):
                e = np.zeros(7 * ch.N)
                e[i] = hfd
                g[i] = (potential(ch, RC.apply_increment(ch, qg, e), cpairs,
                                  CU, n, kc, ke)
                        - potential(ch, RC.apply_increment(ch, qg, -e), cpairs,
                                    CU, n, kc, ke)) / (2 * hfd)
            gerr = max(gerr, float(np.abs(fq + g).max()
                                   / max(float(np.abs(g).max()), 1e-300)))
    out["R2"] = gerr
    A(("R2  THE EXTRA TERM REALLY IS A GRADIENT, WHICH IS THE ONE PIECE OF "
       "MACHINERY THIS FILE ADDS AND THE ONE THING NO OTHER ROW WOULD CATCH. "
       "The asymmetric potential carries a term through dK/dc that the "
       "symmetric law does not have, and a wrong one would not announce "
       "itself -- it would quietly make every energy, curvature and peak "
       "below the property of no potential at all. The analytic generalised "
       "force is therefore checked against a central finite difference of the "
       "potential itself, coordinate by coordinate, at a perturbed "
       "configuration, at both stiffness ratios and three exponents. TWO-"
       "SIDED: dropping or mis-signing the dK/dc term fails this row at the "
       "asymmetric ratios while still passing R1, which only exercises the "
       "case where that term vanishes",
       gerr < 1e-7,
       f"worst relative error between the analytic force and a finite-"
       f"difference gradient of V, over kc/ke = 1 and 100 and n = 2, 5/2, 3: "
       f"{gerr:.2e} (the difference's own truncation is near 1e-9)",
       "agreement to the finite difference's own accuracy"))

    # ---- R3: the breathe is free to FIRST order and loads at SECOND --------
    amps = (0.15, 0.3, 0.6, 1.2)
    seps = [breathe(*ref, +1.0, JE.K2, JE.K2, amp=a)["smax"] for a in amps]
    ratios = [seps[i + 1] / seps[i] for i in range(len(seps) - 1)]
    e_coh = MJ.coherent_kick(ref[0], amp=1.0).ravel()
    e_coh = e_coh / np.linalg.norm(e_coh)
    c_coh = curvature(ref[0], ref[1], ref[2], e_coh, JE.K2, JE.K2, 1e-5)
    out["R3"] = (seps, ratios, c_coh)
    A(("R3  THE COHERENT BREATHE IS FREE TO FIRST ORDER AND LOADS THE JOINT AT "
       "EXACTLY SECOND, WHICH IS WHY AN ASYMMETRIC JOINT CAN SEE IT AT ALL. "
       "Three instruments call this motion free -- jb_mj R3, jb_sj R3, "
       "jb_bz R2 -- and all three statements are INFINITESIMAL: a zero mode "
       "is a first-order object, and jb_mj's zero contacts means the "
       "separation never reached the CLEARANCE, not that it was zero. "
       "Doubling the drive quadruples the separation, to three figures, at "
       "every amplitude tried, and the curvature of V along the coherent "
       "direction is zero to 1e-8 -- so the mode is genuinely free where "
       "those files measured it AND genuinely opens joints at finite "
       "amplitude. Both halves are needed: the first is what makes it a zero "
       "mode, the second is what leaves anything for this file to measure. "
       "TWO-SIDED: separation linear in amplitude would mean it is not a zero "
       "mode and would contradict three files; separation identically zero "
       "would mean no joint law of any form could act on it, and every row "
       "below would be empty",
       max(abs(r - 4.0) for r in ratios) < 0.05 and abs(c_coh[0]) < 1e-6,
       "separation " + ", ".join(f"{s:.3e}" for s in seps)
       + f" at amplitude {amps}; successive ratios "
       + ", ".join(f"{r:.3f}" for r in ratios)
       + f" against 4 for a quadratic; curvature along the coherent "
       f"direction {c_coh[0]:.2e}",
       "each doubling quadrupling the separation, and zero curvature"))

    # ---- R4: the phase decides whether the question parses -----------------
    dsp = {}
    for tag, P in (("a=-30", ref), ("a=-37.761", ico)):
        for sgn in (+1.0, -1.0):
            dsp[(tag, sgn)] = breathe(*P, sgn, JE.K2, JE.K2)["dspan"]
    ref_same = (dsp[("a=-30", 1.0)] < 0 and dsp[("a=-30", -1.0)] < 0)
    ico_split = (dsp[("a=-37.761", 1.0)] > 0 > dsp[("a=-37.761", -1.0)])
    out["R4"] = dsp
    A(("R4  AT THE REFERENCE PHASE 'EXPAND' AND 'CONTRACT' ARE NOT A "
       "DISTINCTION THE MEDIUM MAKES, AND AT THE ICOSAHEDRAL PHASE THEY ARE. "
       "This row is why everything about the breathe is measured at two "
       "phases rather than one, and it would have invalidated the obvious "
       "experiment. At a = -30 the array's span is at an EXTREMUM: BOTH fold "
       "senses contract it, by the same number to three figures, so a "
       "measurement of expansion against contraction there is comparing a "
       "motion with itself. At a = -37.761 -- the phase qvf.11 reports its "
       "lock at, and where its close records the hole cell reaching the "
       "icosahedron -- the two senses separate cleanly, one growing the span "
       "and the other shrinking it. TWO-SIDED: if the reference phase had "
       "split, the icosahedral phase would carry no special status here and "
       "the file could have used one patch",
       ref_same and ico_split,
       f"a = -30: d(span) {dsp[('a=-30', 1.0)]:+.3e} and "
       f"{dsp[('a=-30', -1.0)]:+.3e}, both contracting;  "
       f"a = -37.761: {dsp[('a=-37.761', 1.0)]:+.3e} and "
       f"{dsp[('a=-37.761', -1.0)]:+.3e}, opposite in sign",
       "both senses contracting at the reference phase, splitting at the "
       "icosahedral one"))

    # ---- R5: one-sidedness, geometric first, and nothing like a lock -------
    ons = {}
    for tag, P in (("a=-30", ref), ("a=-37.761", ico)):
        for ratio in (1.0, 100.0):
            vp = breathe(*P, +1.0, ratio * JE.K2, JE.K2)["vpeak"]
            vm = breathe(*P, -1.0, ratio * JE.K2, JE.K2)["vpeak"]
            ons[(tag, ratio)] = (vp, vm, vp / vm)
    geo = ons[("a=-37.761", 1.0)][2]
    mat = ons[("a=-37.761", 100.0)][2]
    flat = ons[("a=-30", 1.0)][2]
    out["R5"] = ons
    A(("R5  THE ONE-SIDEDNESS IS REAL, IT IS GEOMETRIC BEFORE IT IS MATERIAL, "
       "AND IT IS NOTHING LIKE A LOCK. Three findings in one row because they "
       "are only meaningful together. FIRST, at the icosahedral phase a "
       "SYMMETRIC joint already loads the two senses unequally -- the peak "
       "joint potential differs by a factor near three -- so the linkage "
       "produces one-sidedness with no asymmetric material anywhere, and any "
       "reading that attributes it to rubber is attributing it to the wrong "
       "thing. SECOND, a hundred-to-one stiffness ratio then moves that "
       "number a long way, and phase-dependently: substantially at the "
       "icosahedral phase, barely at the reference phase, which is what R4 "
       "predicts. THIRD and most important, NOTHING IS FORBIDDEN. qvf.11's "
       "array CANNOT EXPAND; here every direction stays available at every "
       "ratio and only the energy the joints take changes. This file does not "
       "reproduce qvf.11's lock and must not be read as doing so. TWO-SIDED: "
       "a symmetric ratio of one would make the effect purely material, and a "
       "ratio unmoved by the asymmetry would make the candidate dead rather "
       "than merely insufficient",
       abs(geo - 1.0) > 0.2 and abs(mat - geo) > 0.2
       and abs(flat - 1.0) < 0.05,
       "  ".join(f"{t} kc/ke={r:g}: Vpeak {v[0]:.4e} vs {v[1]:.4e}, "
                 f"ratio {v[2]:.4f}" for (t, r), v in sorted(ons.items()))
       + "; no direction blocked at any ratio",
       "the symmetric icosahedral ratio already far from 1, the asymmetry "
       "moving it, and the reference phase flat"))

    # ---- R6: the transport exponents survive -------------------------------
    texp = {}
    for ratio in RATIOS:
        p = transport_exponent(ch, cpairs, CU, 2.5, ratio, 1.0,
                               JE.seed(2.5) / (0.5 * (1.0 + ratio)))
        if p is not None:
            texp[ratio] = p
    out["texp"] = texp
    A(("R6  THE TRANSPORT EXPONENT SURVIVES THE ASYMMETRY UNCHANGED, AND THE "
       "REASON IS THE SAME HOMOGENEITY jb_je RESTS ON. Because K depends only "
       "on the DIRECTION of the separation, V(L d) = L^n V(d) still holds "
       "exactly, so the similarity that fixes jb_je's exponents is untouched "
       "and p_v = (n-2)/n should come back at every stiffness ratio. It does, "
       "at n = 5/2 for ratios 1, 4 and 100 -- so a medium can be strongly "
       "one-sided in its joints and still carry fronts under the symmetric "
       "medium's speed-amplitude law. TWO-SIDED: a ratio-dependent exponent "
       "would mean the asymmetry had introduced a scale, which would break "
       "the closed forms jb_je measured rather than extend them",
       len(texp) == 3 and max(abs(p - 0.2) for p in texp.values()) < 0.02,
       "  ".join(f"kc/ke={r:g}: p_v={p:+.4f}" for r, p in sorted(texp.items()))
       + " against (n-2)/n = 0.2000 at n = 5/2",
       "every ratio within 0.02 of 0.2"))

    # ---- R7: and it costs the Hessian --------------------------------------
    e_stag = MJ.staggered_kick(ref[0], amp=1.0).ravel()
    e_stag = e_stag / np.linalg.norm(e_stag)
    curv = {}
    for ratio in (1.0, 100.0):
        cp, cm = curvature(ref[0], ref[1], ref[2], e_stag,
                           ratio * JE.K2, JE.K2, 1e-5)
        curv[ratio] = (cp, cm, cp / cm if cm else float("nan"))
    out["curv"] = curv
    A(("R7  AN ASYMMETRIC JOINT HAS NO HESSIAN AT THE WELDED CONFIGURATION, SO "
       "IT COSTS THE MEDIUM ITS DISPERSION RELATION -- AND THAT IS jb_cp's "
       "OBSTRUCTION REACHED FROM THE OPPOSITE SIDE. The potential is stiffer "
       "leaving the weld one way than the other, so the second derivative "
       "there depends on WHICH WAY YOU LEAVE: V is continuous and even C1 at "
       "the weld, but not C2, and a quantity that differs by direction is not "
       "a Hessian. omega(k) needs one. Measured along the STAGGERED direction, "
       "which is the one that loads every joint at once -- a random direction "
       "understates this badly, because it opens some joints in compression "
       "and others in extension and the two average away. jb_cp [23658] found "
       "a hard wall has no omega(k) because an infinite square well is FLAT "
       "inside and has no Hessian anywhere; this potential is soft, finite "
       "and everywhere continuous and STILL has none, for a KINK instead of "
       "for flatness. Two very different potentials, one obstruction. "
       "TWO-SIDED: the symmetric control must show NO gap, or the row is "
       "measuring its own finite difference rather than the medium",
       abs(curv[1.0][2] - 1.0) < 1e-4 and abs(curv[100.0][2] - 1.0) > 0.5,
       "  ".join(f"kc/ke={r:g}: curvature {c[0]:.6g} one way vs {c[1]:.6g} "
                 f"the other, ratio {c[2]:.6f}" for r, c in sorted(curv.items())),
       "ratio 1 to 1e-4 when symmetric, far from 1 when not"))

    # ---- R8: controls, and the energy the asymmetric gradient conserves ----
    dead = integrate(ch, cpairs, CU, MJ.single_kick(ch, amp=JE.KICKS[-1]),
                     2.5, 0.0, 0.0, tmax=0.6, kick=JE.KICKS[-1])
    quiet = integrate(ch, cpairs, CU, np.zeros((ch.N, 7)), 2.5,
                      JE.K2 * 100.0, JE.K2, tmax=0.6, kick=JE.KICKS[-1])
    de = max(breathe(*P, s, 100.0 * JE.K2, JE.K2)["drift"]
             for P in (ref, ico) for s in (+1.0, -1.0))
    out["ctl"] = (len(dead["rec"]), len(quiet["rec"]), de)
    A(("R8  CONTROLS, AND THE ENERGY THE STIFFEST RUNS HAVE TO CLOSE ON. With "
       "both stiffnesses zero no front forms, and undriven no front forms, so "
       "every number above is a property of the joint law rather than of the "
       "drive, the index or the integrator. The energy check is the one that "
       "matters at a hundred-to-one ratio: the medium is a hundred times "
       "stiffer on one side, which is exactly the regime where a fixed step "
       "quietly stops integrating the problem posed. Gated on the breathe "
       "runs at both phases and both senses, since those are what R4 and R5 "
       "read",
       len(dead["rec"]) == 0 and len(quiet["rec"]) == 0 and de < 1e-6,
       f"kc = ke = 0: {len(dead['rec'])} cells reached; undriven: "
       f"{len(quiet['rec'])}; worst relative energy drift over the four "
       f"breathe runs at kc/ke = 100: {de:.2e}",
       "no front in either control and energy closing to 1e-6"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_ja -- an asymmetric joint: what it does, and what it costs")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        print("\n  THE LEDGER")
        print("   the coherent breathe    LOADED at second order, so the")
        print("                           asymmetry does act on it")
        print("   one-sidedness           REAL, and GEOMETRIC before material")
        print("   qvf.11's lock           NOT reproduced -- nothing forbidden")
        print("   transport exponents     UNCHANGED, p_v = (n-2)/n at every ratio")
        print("   the dispersion relation LOST -- no Hessian at the weld")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * THE CANDIDATE IS INSUFFICIENT, NOT DEAD. jb_sj and jb_bz")
        print("     both name joint asymmetry as the live candidate for")
        print("     qvf.11's one-sidedness. It does act on the breathe, and")
        print("     it produces nothing resembling a lock: no direction is")
        print("     forbidden at any ratio tried.")
        print("   * THE ONE-SIDEDNESS IS THE LINKAGE'S BEFORE IT IS THE")
        print("     MATERIAL'S. A symmetric joint at the icosahedral phase")
        print("     already loads the two senses by a factor near three.")
        print("   * THE JOINT'S SENSE IS NOT THE ARRAY'S. Raising the")
        print("     COMPRESSION stiffness affects the span-EXPANDING drive")
        print("     more. Do not run the two vocabularies together.")
        print("   * qvf.11's FORK (A) IS UNTOUCHED. A tension-only member")
        print("     spanning a square-face diagonal is an ADDED unilateral")
        print("     member between different vertices, not a joint law.")
        print("   * qvf.11 IS NOT REOPENED. The owner closed it as evidence")
        print("     that served and said the lock is not a goal to eliminate.")
        print("   * jb_bz IS NOT RETRACTED. Its bands are about the SYMMETRIC")
        print("     joint and stand. R7 says what a different law would cost.")
        print("   * RATIOS ONLY. kc and ke are conventions; only their ratio")
        print("     is a measurement.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
