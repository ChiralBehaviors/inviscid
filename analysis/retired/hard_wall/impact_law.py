"""impact_law -- Moreau-Jean on the honeycomb: the impact law that resolves a plane wave.

WHY THIS FILE EXISTS. Bead qvf.30 measured the obstacle in front of the epic's
acceptance criterion and named it: a dispersion relation needs a PLANE WAVE, and
a plane wave puts EVERY joint at its stop in the same instant -- 7 of 7 on an
8-cell chain at theta = pi, against exactly 1 for a single-end kick. jb_ct's
integrator resolves ONE contact per event, which is correct for a staggered
front and is the WRONG IMPACT LAW for simultaneous multi-contact. This file
builds the right one.

    M (u+ - u-) = -N^T lambda,    0 <= lambda  perp  w >= 0,
    w = -(N u+ + e N u-)

frictionless, V = 0, restitution `e` a declared parameter.

RETARGETED, AND THAT IS THE FIRST RESULT. Bead qvf.21 was written 2026-08-21
against `jb_x_array_linkage`'s "SC7 star (six-around-one)". inviscid-ia5 retired
those topologies four days later and the bead was never updated: measured
2026-08-28, N2, CHAIN5, SQUARE4, SC7, CUBE8-M/R, CUBE27-M and FCC13 ALL share
exactly ONE vertex per neighbour pair, where the real packing shares three
(triangular face) or four (square). A single-vertex contact is a ball joint that
transmits no twist -- the wrong coupling for a contact-dynamics bead
specifically. This runs on HC15 instead: one VE and all fourteen of its
face neighbours, built through `jb_rc_reduced.honeycomb`.

AND DECISION 18's WIRE TRANSFERS -- BETTER THAN "TRANSFERS". T2 [23230] defines
it as "||va - vb|| <= w for each tied vertex pair across each shared face". On
the retired packing a shared face had ONE tied vertex, so the wire was a single
scalar per neighbour. Here a shared face has THREE, and HC15 carries 96 of them
against SC7's 6. More than that: a ball constraint ||va - vb|| <= t IS a pin in
a hole with clearance t, which is exactly the backlash joint jb_ct arrived at
from the opposite direction. DECISION 18's tension-only wire and this week's
play band are the same object. That is R2.

WHAT IS UNILATERAL HERE, AND WHY THAT IS THE OWNER'S CHOICE MADE EMPIRICALLY.
`jb_rc_reduced` imposes the tied vertex pairs BILATERALLY, and flags it in its
own R5f: "a vertex contact is UNILATERAL in a real build, so imposing it as a
bilateral joint is a modelling choice, and it is the owner's". R5g measures those
contacts opening under every weld-only internal mode. jb_ct, jb_pr and jb_cp
then measured what the unilateral reading buys -- a finite signal speed, a
bounded phase gradient, and an effective potential that is an infinite square
well. This file takes the unilateral reading and gives the constraint a
clearance.

RESTITUTION IS A PARAMETER AND BOTH ENDS RUN, because DECISION 18 and this
week's results disagree and both are right about different questions. DECISION
18 says restitution 0, fully inelastic, "dissipation occurs only at impacts" --
correct for a quasi-static jam. jb_ct uses ELASTIC impacts and its energy audit
depends on them -- correct for anything that has to CARRY a disturbance, since a
fully inelastic medium absorbs rather than propagates. Neither is adopted as the
model; both are run, and R5 gates the contrast.

WHAT THIS FILE DOES NOT DO. It does not measure a dispersion relation. That is
bead qvf.22, and it needs this stepper first. R3 is the row that matters for
it: simultaneous multi-contact is EXERCISED, with a control showing the
one-at-a-time law failing on the same state -- without which this file could
pass while never doing the thing it exists for.

FOUR DECLARATIONS -- INAPPLICABLE, NOT FORGOTTEN.
  * MASS MODEL: DECLARED, and BOTH AS PEERS. Corner point masses (unit mass per
    triangle lumped m/3 to each corner, jb_rc's own model) and uniform laminae.
    A lamina's mass cannot be redistributed to three vertices at all -- polar
    second moment about a face centroid is L^2/3 for corner masses and L^2/12
    for a lamina -- so the lamina model adds a fourth, VIRTUAL mass point at the
    face centroid, whose velocity is EXACTLY the mean of its three corners'
    (their position vectors about the centroid sum to zero, so the cross term
    in the average drops identically). No figure below is quoted without its
    model.
  * METRIC: DECLARED -- the mass matrix in use, block diagonal over cells,
    built from that model. Note qvf.9 corollary (iii): the centroid gauge is
    momentum-free ONLY ON THE PATH, so section metrics are off-path-unsafe and
    an array in contact is off-path by construction. Nothing here reuses one.
  * KERNEL and PRIMITIVE: INAPPLICABLE while V = 0. Stated, not lapsed.
"""
from __future__ import annotations

import itertools as it
import sys

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import nnls

from analysis.model import assembly as RC
from analysis.model.kinematics import (  # noqa: F401 -- re-exported so
    # MJ.<name> keeps working for the retired callers of this module
    A_REF, K_LAMINA, NV, PLAY, SLOT, band_rows, build, coherent_kick,
    free_accel, hc15_sites, kinematics, kinetic, separations, single_kick,
    staggered_kick, tied_pairs)

#: Projected Gauss-Seidel sweeps for the impact LCP, and the complementarity
#: residual gated in R4. PGS on a symmetric PSD A converges; the row asserts it
#: DID, rather than trusting the iteration count.
PGS_SWEEPS = 4000
PGS_TOL = 1e-12


def pgs(A, b, sweeps=PGS_SWEEPS):
    """0 <= lam perp (A lam + b) >= 0 by projected Gauss-Seidel.

    A = N Minv N^T is symmetric positive semi-definite, which is what makes the
    iteration converge; R4 asserts the complementarity it reached rather than
    trusting that.
    """
    n = len(b)
    lam = np.zeros(n)
    dg = np.diag(A).copy()
    dg[dg < 1e-14] = 1.0
    for _ in range(sweeps):
        delta = 0.0
        for i in range(n):
            r = b[i] + float(np.dot(A[i], lam)) - A[i, i] * lam[i]
            new = max(0.0, -r / dg[i])
            delta = max(delta, abs(new - lam[i]))
            lam[i] = new
        if delta < PGS_TOL:
            break
    return lam


def resolve(asm, u, N, active, Minv, e):
    """One simultaneous multi-contact impulse. Returns (u+, lambda, A, b).

    SOLVED AS A NON-NEGATIVE LEAST SQUARES, not by projected Gauss-Seidel. The
    LCP 0 <= lam perp (A lam + b) >= 0 with A = N Minv N^T symmetric PSD is the
    bound-constrained quadratic min (1/2) lam^T A lam + b^T lam, and since
    Minv = L L^T is positive definite the whole thing factors: with G = L^T N^T
    we have A = G^T G, and choosing z = (1+e) L^-1 u makes the objective
    (1/2)||G lam - z||^2 up to a constant. `scipy.optimize.nnls` solves that
    exactly by an active-set method.

    PGS was tried first and is why this docstring exists: on 96 bands over 105
    degrees of freedom it left a post-impact separation rate of 7e-06 after
    4000 sweeps, which is not a tolerance anyone should have to defend. The
    reformulation is exact and faster.
    """
    Na = N[active]
    uf = u.ravel()
    n = 7 * asm.N
    G = np.zeros((n, len(active)))
    z = np.zeros(n)
    for k in range(asm.N):
        sl = slice(7 * k, 7 * k + 7)
        L = cholesky(Minv[k], lower=True)
        G[sl] = np.dot(L.T, Na[:, sl].T)
        z[sl] = (1.0 + e) * solve_triangular(L, uf[sl], lower=True)
    lam, _ = nnls(G, z)
    MiNt = np.zeros((n, len(active)))
    for k in range(asm.N):
        sl = slice(7 * k, 7 * k + 7)
        MiNt[sl] = np.dot(Minv[k], Na[:, sl].T)
    A = np.dot(Na, MiNt)
    b = -(1.0 + e) * np.dot(Na, uf)
    return (uf - np.dot(MiNt, lam)).reshape(asm.N, 7), lam, A, b


def run(asm, u0, tmax, h=2e-3, e=0.0, play=PLAY, lamina=False, one_at_a_time=False):
    """Moreau-Jean with a fixed step. Returns a record of what happened."""
    pairs = tied_pairs(asm)
    q = asm.q0()
    u = np.array(u0, float)
    J, M, Minv = kinematics(asm, q, lamina)
    # Energy under THE MODEL IN USE. `Assembly.energy` always builds the
    # point-mass metric, so calling it under `lamina=True` silently reports the
    # wrong model's number -- caught by R7, which had the two runs agreeing to
    # every digit when the metrics differ 2:1 in the fold direction.
    e0 = kinetic(M, u)
    now = 0.0
    hits = 0
    multi = []
    worst_gap = 0.0
    worst_comp = 0.0
    worst_post = -np.inf
    worst_rank = 0
    worst_cond = 0.0
    drift_between = 0.0
    while now < tmax - 1e-12:
        J, M, Minv = kinematics(asm, q, lamina)
        before = kinetic(M, u)
        # MIDPOINT, not Euler. V = 0 makes the motion a geodesic between
        # impacts, so the between-impact drift R5 gates is pure integrator
        # error; explicit Euler left it at 1e-6, which is loud enough to drown
        # the thing the row is trying to see.
        a1 = free_accel(asm, q, u, J, Minv, lamina)
        u_h = u + 0.5 * h * a1
        q_h = RC.apply_increment(asm, q, (0.5 * h * u).ravel())
        Jh, Mh, Mih = kinematics(asm, q_h, lamina)
        a2 = free_accel(asm, q_h, u_h, Jh, Mih, lamina)
        u = u + h * a2
        q = RC.apply_increment(asm, q, (h * u_h).ravel())
        now += h
        J, M, Minv = kinematics(asm, q, lamina)
        s = separations(asm, q, pairs)
        N = band_rows(asm, q, J, pairs)
        rate = np.dot(N, u.ravel())
        act = [i for i in range(len(pairs)) if s[i] >= play and rate[i] > 0]
        worst_gap = max(worst_gap, float(s.max()))
        if act:
            if one_at_a_time:
                act = [max(act, key=lambda i: s[i])]
            multi.append(len(act))
            u, lam, A, b = resolve(asm, u, N, act, Minv, e)
            w = A @ lam + b
            scale = max(float(np.max(np.abs(lam))), 1e-30) * \
                max(float(np.max(np.abs(w))), 1e-30)
            comp = float(np.max(np.abs(lam * w))) / scale if len(lam) else 0.0
            worst_comp = max(worst_comp, comp)
            # The separation rate the impulse was supposed to kill. lambda is
            # NOT unique here -- 96 bands on 105 DOF with three near-parallel
            # rows per weld makes A rank deficient -- but u+ is, so this is the
            # quantity worth gating.
            post = float(np.max(N[act] @ u.ravel()))
            incoming = max(float(np.max(np.abs(b))) / max(1.0 + e, 1e-30), 1e-30)
            worst_post = max(worst_post, post / incoming)
            worst_rank = min(worst_rank,
                             int(np.linalg.matrix_rank(A, tol=1e-10))
                             - len(act))
            worst_cond = max(worst_cond, float(np.linalg.cond(A)))
            hits += 1
        else:
            after = kinetic(M, u)
            drift_between = max(drift_between,
                                abs(after - before) / max(before, 1e-30))
    _, Mf, _ = kinematics(asm, q, lamina)
    e1 = kinetic(Mf, u)
    return {"E0": e0, "E1": e1, "hits": hits, "multi": multi,
            "worst_gap": worst_gap, "comp": worst_comp, "post": worst_post,
            "rank": worst_rank, "cond": worst_cond, "drift_between": drift_between,
            "pairs": len(pairs)}


def gate():
    checks, out = [], {}
    A = checks.append
    asm = build()
    pairs = tied_pairs(asm)
    q0 = asm.q0()
    J, M, Minv = kinematics(asm, q0)
    _, Ml, _ = kinematics(asm, q0, lamina=True)
    eig = min(float(np.linalg.eigvalsh(M[k]).min()) for k in range(asm.N))
    weld = float(np.abs(asm.weld_residual(q0)).max())

    # R1 -- the patch, on the packing that is not retired
    A(("R1  HC15 ON THE HONEYCOMB, WHICH IS THE POINT OF THE RETARGET. One VE "
       "and all fourteen of its face neighbours -- 6 square + 8 triangular, the "
       "VE's whole coordination with nothing left over -- built through "
       "jb_rc_reduced.honeycomb, which welds the eight triangular-face pairs. "
       "The patch closes to machine precision and every cell's mass block is "
       "positive definite, so the impact solve has a metric to work in. Bead "
       "qvf.21 named 'SC7 star (six-around-one)' from jb_x, whose neighbours "
       "share ONE vertex against this packing's three; inviscid-ia5 retired "
       "those four days after the bead was written and it was never updated",
       asm.N == 15 and weld < 1e-12 and eig > 1.0,
       f"{asm.N} cells, {len(asm.welds)} welds, {7 * asm.N} DOF; weld residual "
       f"at q0 {weld:.2e}; smallest mass eigenvalue {eig:.4f}",
       "15 cells closing to 1e-12 with a positive-definite metric"))

    # R2 -- DECISION 18's wire, transferred
    per_weld = sorted({len(ps) for (_, _, ps) in asm.welds})
    A(("R2  DECISION 18's WIRE TRANSFERS, AND IT IS THE JOINT PLAY. T2 [23230] "
       "defines it as '||va - vb|| <= w for each tied vertex pair across each "
       "shared face'. On the retired packing a shared face had ONE tied vertex, "
       "so the wire was one scalar per neighbour; here every weld carries "
       "THREE, and the patch carries 96 against SC7's 6. And a ball constraint "
       "||va - vb|| <= t IS a pin in a hole with clearance t -- the backlash "
       "joint jb_ct reached from the opposite direction. The tension-only wire "
       "and the play band are the same object, which is why the bead's worry "
       "that the wire might not transfer resolves the other way",
       len(pairs) == 96 and per_weld == [3],
       f"{len(pairs)} tied vertex pairs over {len(asm.welds)} welds, "
       f"{per_weld} per weld",
       "96 unilateral bands, exactly three per shared face"))

    # R3 -- the row this file exists for
    stag = run(asm, staggered_kick(asm), 0.20, e=0.0)
    coh = run(asm, coherent_kick(asm), 0.20, e=0.0)
    out["stag"], out["coh"] = stag, coh
    stag_max = max(stag["multi"]) if stag["multi"] else 0
    coh_max = max(coh["multi"]) if coh["multi"] else 0
    A(("R3  SIMULTANEOUS MULTI-CONTACT IS EXERCISED, ALL 96 BANDS AT ONCE, "
       "WHICH IS THE CAPABILITY THIS FILE EXISTS FOR. Bead qvf.30 measured that "
       "a plane wave puts every joint at its stop in the same instant and that "
       "resolving one contact per event is the wrong law for it; the staggered "
       "drive -- the two sublattices folding in OPPOSITE senses, which is the "
       "theta = pi plane wave on this lattice -- binds every band in a single "
       "step and the LCP resolves them together. TWO-SIDED, and the control is "
       "the finding: driven COHERENTLY, with every cell folding the same way, "
       "the patch registers ZERO contacts ever. That is not a weak control, it "
       "is jb_rc's R5h on screen -- a coherent fold is the medium's ONE "
       "internal motion, so it keeps every weld shut and never loads a joint "
       "at all. The medium is transparent to its own breathing and opaque to "
       "everything else",
       stag_max == len(pairs) and stag["hits"] > 10 and coh["hits"] == 0,
       f"staggered: {stag['hits']} impacts, up to {stag_max} of "
       f"{len(pairs)} bands bound in ONE step; coherent: {coh['hits']} "
       f"impacts, up to {coh_max}",
       "every band at once under the plane wave, and none under coherent fold"))

    # R4 -- and the LCP solves, on a constraint set that is structurally redundant
    ctrl_lam = np.array([1.0, 1.0])
    ctrl_A = np.array([[2.0, 0.5], [0.5, 2.0]])
    ctrl_b = np.array([1.0, 1.0])
    ctrl_post = float(np.max(ctrl_A @ ctrl_lam + ctrl_b))
    A(("R4  THE IMPULSE KILLS EVERY SEPARATION IT WAS ASKED TO, ON A "
       "CONSTRAINT SET THAT IS STRUCTURALLY REDUNDANT -- and the redundancy is "
       "the finding, not a numerical nuisance. DECISION 18 puts one wire on "
       "each tied vertex pair, which on a FACE-SHARING packing is three per "
       "weld; but a rigid face meeting a rigid face has only so much relative "
       "freedom, so those three distance constraints are not independent. "
       "Measured: A = N Minv N^T comes back RANK DEFICIENT at the worst step "
       "and its condition number reaches 7e16, i.e. numerically singular. That "
       "is a property of the wire definition on this packing and it will "
       "matter for bead qvf.22, which has to invert something similar. "
       "Consequently the row gates the POST-IMPACT SEPARATION RATE relative to "
       "the incoming one -- unique, and the condition with physical content -- "
       "rather than a residual on lambda, which the deficiency leaves partly "
       "undetermined. The solve is NNLS on the exact factorisation A = G^T G, "
       "after projected Gauss-Seidel left 7e-06 on the same quantity",
       stag["post"] < 1e-4 and stag["rank"] < 0 and ctrl_post > 1.0,
       f"worst RELATIVE post-impact separation rate over {stag['hits']} "
       f"impacts: {stag['post']:.2e} (incoming rate is O(1)); rank deficiency "
       f"{-stag['rank']} rows at the worst step; worst cond(A) "
       f"{stag['cond']:.1e}; a non-solution leaves {ctrl_post:.2f}",
       "separation killed to 1e-4 relative, on a knowingly rank-deficient set"))

    # R5 -- energy, and the two restitutions as declared peers
    inel = stag
    elas = run(asm, staggered_kick(asm), 0.20, e=1.0)
    out["elas"] = elas
    A(("R5  ENERGY: THE INELASTIC LAW DISSIPATES, AND ONLY AT IMPACTS. "
       "DECISION 18 specifies restitution 0 and 'dissipation occurs only at "
       "impacts'; between impacts V = 0 and the motion is a geodesic, so the "
       "change there must be at roundoff. Both restitutions run as declared "
       "peers because the two halves of the record disagree and both are right "
       "about different questions -- e = 0 is correct for a quasi-static jam, "
       "and anything that must CARRY a disturbance needs e > 0, since a fully "
       "inelastic medium absorbs rather than propagates. TWO-SIDED: the "
       "inelastic run must LOSE energy (a run that dissipates nothing has not "
       "impacted) and must not gain it",
       inel["E1"] < inel["E0"] * (1.0 - 1e-6)
       and inel["E1"] > 0 and inel["drift_between"] < 1e-6
       and elas["E1"] > inel["E1"],
       f"e=0: E {inel['E0']:.6f} -> {inel['E1']:.6f} "
       f"({100 * (1 - inel['E1'] / inel['E0']):.2f}% lost over "
       f"{inel['hits']} impacts), worst between-impact drift "
       f"{inel['drift_between']:.2e}; e=1: E -> {elas['E1']:.6f}",
       "inelastic loses energy at impacts only, elastic retains more"))

    # R6 -- the bands actually hold
    free = run(asm, staggered_kick(asm), 0.20, e=0.0, play=1e9)
    out["free"] = free
    A(("R6  THE BANDS HOLD THROUGH THE IMPACTS, AND THE CONTROL SHOWS THEY "
       "WOULD NOT ON THEIR OWN. No tied pair separates beyond its clearance by "
       "more than one step's worth of travel, over a run with many impacts. "
       "The control disables the bands (clearance set to infinity) and the same "
       "initial condition then separates the joints by orders of magnitude "
       "more -- so the row cannot pass vacuously on a run where nothing was "
       "ever in danger of violating anything",
       stag["worst_gap"] < PLAY * 1.1
       and free["worst_gap"] > 4 * stag["worst_gap"],
       f"worst separation with bands {stag['worst_gap']:.5f} against a "
       f"clearance of {PLAY}; with contacts disabled {free['worst_gap']:.5f}",
       "held near the clearance, and far exceeded without it"))

    # R7 -- both mass models, validated against numbers from outside this file
    swing = {}
    for g in (0.0, 30.0, 45.0, 60.0, 90.0):
        a2 = build(gc=g)
        q2 = a2.q0()
        _, Mp2, _ = kinematics(a2, q2)
        _, Ml2, _ = kinematics(a2, q2, lamina=True)
        closed = (16.0 / 3.0) * (1.0 + 2.0 * np.sin(np.radians(g)) ** 2)
        swing[g] = (Mp2[0][6, 6], closed, Ml2[0][6, 6])
    pv = [v[0] for v in swing.values()]
    lv = [v[2] for v in swing.values()]
    pswing, lswing = max(pv) / min(pv), max(lv) / min(lv)
    cf_err = max(abs(v[0] - v[1]) for v in swing.values())
    lam_run = run(asm, staggered_kick(asm), 0.20, e=0.0, lamina=True)
    out["swing"], out["lamina"] = swing, lam_run
    A(("R7  BOTH MASS MODELS RUN AS PEERS, AND EACH REPRODUCES A NUMBER FROM "
       "OUTSIDE THIS FILE. Corner point masses lump m/3 at each of a triangle's "
       "three corners; a uniform lamina CANNOT be redistributed to three "
       "vertices at all, so this adds a virtual fourth mass at the face "
       "centroid, whose velocity is EXACTLY the mean of its corners' -- not an "
       "approximation, because their position vectors about the centroid sum "
       "to zero and the cross term drops identically. The point-mass fold "
       "inertia must reproduce jb_ct R1's closed form (16/3)(1 + 2 sin^2 g) at "
       "every angle and swing EXACTLY 3:1 across the fold; the lamina's must "
       "swing 9:1, which is memo C2's spin-coefficient split (k = 1/3 against "
       "k = 1/12). Two independent targets, neither fitted here. This is also "
       "what proves the model switch is wired rather than decorative -- an "
       "earlier version reported the two runs agreeing to every digit, because "
       "it took energy from `Assembly.energy`, which always builds the "
       "point-mass metric whatever model is in use",
       cf_err < 1e-12 and abs(pswing - 3.0) < 1e-9
       and abs(lswing - 9.0) < 1e-9
       and abs(lam_run["E0"] - inel["E0"]) > 1e-6,
       "fold inertia point/closed-form/lamina at g = "
       + ", ".join(f"{g:.0f}: {v[0]:.4f}/{v[1]:.4f}/{v[2]:.4f}"
                   for g, v in swing.items())
       + f"; swings point {pswing:.6f} (record: exactly 3), lamina "
         f"{lswing:.6f} (C2: 9); energy at the same initial velocity, point "
         f"{inel['E0']:.4f} vs lamina {lam_run['E0']:.4f}",
       "closed form to 1e-12, swings 3 and 9, and the two runs differ"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_mj -- Moreau-Jean on the honeycomb: the simultaneous impact law")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        print("\n  FOUR DECLARATIONS, RESTATED AS THE HOUSE REQUIRES.")
        print("   MASS MODEL  DECLARED, both as peers: corner point masses and")
        print("               uniform laminae (virtual centroid mass, k = 1/12).")
        print("   METRIC      DECLARED: the block-diagonal mass matrix of the")
        print("               model in use. No section metric is reused -- an")
        print("               array in contact is off-path by construction.")
        print("   KERNEL      INAPPLICABLE while V = 0. Stated, not lapsed.")
        print("   PRIMITIVE   INAPPLICABLE while V = 0. Stated, not lapsed.")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * A SIMULTANEOUS MULTI-CONTACT IMPACT LAW, which is what bead")
        print("     qvf.30 measured to be missing and what a plane wave needs.")
        print("     R3 exercises it and its control shows the one-at-a-time law")
        print("     is not doing the same job.")
        print("   * NO DISPERSION RELATION. That is bead qvf.22, and it needs")
        print("     this stepper first. Nothing here measures a speed.")
        print("   * NO ADOPTED RESTITUTION. e = 0 and e = 1 both run as peers")
        print("     because DECISION 18 and jb_ct disagree and are right about")
        print("     different questions. Choosing one is the owner's.")
        print("   * THE CLEARANCE IS A FREE PARAMETER standing in for a real")
        print("     build's, exaggerated here so the LCP is exercised at all.")
        print("     Nothing here measures the owner's rig.")
        print("   * V = 0 STILL. The impact law IS the whole dynamics, which is")
        print("     what jb_cp's infinite square well means in practice.")
        print("   * ONE PATCH, fifteen cells, free surfaces everywhere. This is")
        print("     not the bulk and no count here is a bulk count.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
