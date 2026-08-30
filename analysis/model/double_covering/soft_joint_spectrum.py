"""soft_joint_spectrum -- soften the joint, and the medium acquires a spectrum.

OWNER DECISION 2026-08-28, and it came with a fact the model did not have: the
physical jitterbugs are SOFT RUBBER with FLEXIBLE JOINTS AT THE VERTICES. So
joint compliance is not a modelling convenience adopted to make waves appear.
It is what the rig is made of, and every model in this project until now has
been stiffer than the object it describes.

THE MODEL, and what it deliberately is NOT:

    V = (1/2) k * sum over tied vertex pairs |x_a - x_b|^2,   STRUTS RIGID

The struts stay perfectly rigid; each cell stays a rigid-triangle linkage with
one fold angle. The only softness is in the joint that ties two cells' shared
vertices together. THIS IS NOT THE STRUT-COMPLIANCE FORK the owner rejected on
2026-08-27 (T2 [23562]: "you turned the struts into springs and are treating
this as a matrix of these springs that can also precess"). Nothing here is a
spring network; the linkage is untouched and the compliance sits where the
rubber sits.

THE OBJECTION THAT HAD TO BE CLEARED FIRST, because it is on record as fatal.
Bead qvf.2's comment of 2026-08-12, "SOFT JOINTS ARE DEAD ANALYTICALLY":

    V = 0.5 k ||C(x)||^2 ... The variety IS {x : C(x) = 0}, so V vanishes
    identically on it. H = k J^T J, and for any internal tangent v,
    v^T H v = k |J v|^2 = 0. EXACTLY ZERO in all six internal directions, at
    every k. ... THE GENERAL DISQUALIFIER: V must be non-constant ON the
    variety.

That argument is CORRECT AND IT IS ABOUT ONE UNIT. R6 reproduces it here as a
control, so this file is not quietly disagreeing with it: a single jitterbug's
six internal freedoms all come back at exactly zero frequency, at every k, just
as the comment says.

IT DOES NOT SURVIVE THE MOVE TO AN ARRAY, and the reason is a measured property
of the honeycomb rather than an argument. On a compact patch, the rigid-weld
analysis finds internal DOF = 1 (jb_rc R5e/R5f): there is exactly ONE motion
that keeps every weld satisfied, the coherent breathe. Every other motion of
the array STRETCHES A JOINT, so every other motion feels the potential. The
disqualifier bites only where the variety is big enough to hide in, and in the
array it is one-dimensional.

WHAT THE MEDIUM THEN HAS, measured on HC15 (15 cells, 105 DOF, 32 welds):

    zero modes   7  =  6 rigid-body  +  1 coherent breathe
    frequencies  98 real, 0.46813 .. 2.08755

The seven zero modes are not counted and assumed -- R3 identifies them. Each of
the six rigid-body motions and the uniform-fold vector lands inside the measured
zero space to 2e-15, and together they SPAN it to 5e-15.

AND IT IS THE SAME FACT AS jb_mj R3, ARRIVED AT FROM THE OTHER SIDE. The
contact model measured that a coherent drive produces ZERO contacts, ever,
because it keeps every weld shut. The potential model measures that the same
motion is the one mode with ZERO frequency. Free of the contacts and free of
the potential are the same freedom, and the two files reach it from opposite
directions without sharing a line of code (R5).

WHAT THIS OPENS, stated as an opening and not as a result. With a Hessian the
medium has normal modes, so a dispersion relation is REACHABLE for the first
time in this project -- omega(k) by Bloch analysis on the honeycomb, which
DECISION 19 reformulated away precisely because it could not exist under a hard
wall. This file does not compute one. It measures that the obstruction is gone.

WHAT IS NOT SETTLED, and one of these is the owner's.
  * THE STIFFNESS k IS A FREE PARAMETER and every frequency here scales as
    sqrt(k). Only RATIOS of frequencies are measurements; no absolute value is.
  * THE FORM OF V IS THE SIMPLEST ONE, a quadratic in the joint separation.
    Real rubber is nonlinear and is much stiffer in compression than in
    extension. It is a modelling choice and it is not made here.

    SUPERSEDED IN ONE CLAUSE, 2026-08-29, DECISION 20 (T2 [23713]). This
    bullet used to call that asymmetry "a live candidate for the one-sidedness
    qvf.11 reports and jb_tr could NOT reproduce". The TARGET is withdrawn --
    the owner withdrew qvf.11's lock as evidence about the medium -- so
    jb_tr's failure to reproduce it is the EXPECTED result rather than a gap.
    jb_ja [23709] built the asymmetry and measured no lock at any ratio.
    Nothing measured in this file is affected.
  * THE HARD WALL IS NOT RECOVERED HERE as a limit of this V. A quadratic well
    is not the stiff limit of a square well; that would need the family
    |d/t|^(2n) with n large. So this file does not supersede jb_ct/jb_sv, it
    sits beside them as a different joint law.

FOUR DECLARATIONS.
  * KERNEL: DECLARED, and for the first time in this project it is not
    inapplicable -- V is a quadratic in the tied-vertex separation, stiffness k.
  * MASS MODEL: DECLARED -- corner point masses, jb_rc's model. The lamina peer
    is built in jb_mj and not swept here; no frequency below is quoted as
    model-independent.
  * METRIC: DECLARED -- that model's block-diagonal mass matrix.
  * PRIMITIVE: the tied VERTEX PAIR, which is what the physical joint is.
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import kinematics as MJ
from analysis.model import assembly as RC
A_REF = MJ.A_REF

#: Joint stiffness. FREE: every frequency scales as sqrt(k), so only ratios
#: below are measurements. Unity keeps the printed numbers readable.
K_JOINT = 1.0

#: Relative cut for calling a generalised eigenvalue zero. R3 gates the zero
#: modes by IDENTIFYING them, not by trusting this number.
ZERO_RTOL = 1e-9


def spectrum(asm, gc_q=None, k=K_JOINT):
    """(eigenvalues, zero-space basis in u-coordinates, mass matrix).

    The Hessian of V = (1/2) k |C|^2 at a configuration where C = 0 is
    k C'^T C', because the second term carries a factor of C. So the spectrum
    is the generalised eigenproblem (k C^T C, M), reduced by a Cholesky of M
    exactly as the Java GeneralizedEigensolver does it.
    """
    q = asm.q0() if gc_q is None else gc_q
    ctr, R, gam, B = asm.frames(q)
    J = asm.cell_jacobians(ctr, R, B)
    M = asm.mass_blocks(J)
    C = asm.constraint_jacobian(J)
    H = k * (C.T @ C)
    n = asm.N
    Mf = np.zeros((7 * n, 7 * n))
    for i in range(n):
        Mf[7 * i:7 * i + 7, 7 * i:7 * i + 7] = M[i]
    L = np.linalg.cholesky(Mf)
    A = np.linalg.solve(L, np.linalg.solve(L, H).T).T
    ev, V = np.linalg.eigh(A)
    cut = ZERO_RTOL * max(ev.max(), 1e-300)
    zi = [i for i in range(len(ev)) if ev[i] < cut]
    Z = np.linalg.solve(L.T, V[:, zi])
    return ev, Z, ctr


def coherent_vector(n):
    """Every cell folding at the same rate, nothing else moving."""
    v = np.zeros(7 * n)
    v[6::7] = 1.0
    return v


def outside(v, Q):
    """How much of `v` lies outside the span of the orthonormal columns Q."""
    v = v / np.linalg.norm(v)
    return float(np.linalg.norm(v - Q @ (Q.T @ v)))


def gate():
    checks, out = [], {}
    A = checks.append
    asm, _ = RC.honeycomb(MJ.hc15_sites(), gc=A_REF)
    n = asm.N
    ev, Z, ctr = spectrum(asm)
    Zq, _ = np.linalg.qr(Z)
    nz = Z.shape[1]
    freqs = np.sqrt(np.clip(ev[nz:], 0.0, None))
    out["nz"], out["freqs"] = nz, freqs

    # R1 -- the model, and the fork it is not
    A(("R1  THE COMPLIANCE IS IN THE JOINT AND THE STRUTS STAY RIGID, which is "
       "what the physical object is: the owner's jitterbugs are SOFT RUBBER "
       "with FLEXIBLE JOINTS AT THE VERTICES (2026-08-28). So every model in "
       "this project until now has been stiffer than the thing it describes. "
       "V = (1/2) k sum |x_a - x_b|^2 over the tied vertex pairs, and NOTHING "
       "ELSE CHANGES -- each cell is still a rigid-triangle linkage with one "
       "fold angle, and there is no spring on any strut. This is NOT the "
       "strut-compliance fork rejected in T2 [23562]. The row checks the model "
       "is what it says: the patch closes, the mass metric is positive "
       "definite, and the Hessian is symmetric positive semi-definite as a "
       "quadratic's must be",
       n == 15 and float(np.abs(asm.weld_residual(asm.q0())).max()) < 1e-12
       and ev.min() > -1e-9 * max(ev.max(), 1.0),
       f"{n} cells, {7 * n} DOF, {len(asm.welds)} welds, weld residual "
       f"{float(np.abs(asm.weld_residual(asm.q0())).max()):.1e}; smallest "
       f"generalised eigenvalue {ev.min():.2e} (must not be negative)",
       "a closed patch and a positive semi-definite spectrum"))

    # R2 -- the recorded objection, and why it does not reach the array
    A(("R2  THE ARRAY HAS EXACTLY SEVEN ZERO MODES, WHICH IS WHY qvf.2's "
       "ANALYTIC DEATH DOES NOT REACH IT. That comment of 2026-08-12 showed "
       "V = 0.5 k ||C||^2 leaves EVERY internal freedom at zero frequency, "
       "because V vanishes identically on the variety the freedoms live in. "
       "The argument is correct and it is about ONE UNIT, where six internal "
       "directions all lie inside that variety. On a compact honeycomb patch "
       "the rigid-weld analysis finds internal DOF = 1 (jb_rc R5e/R5f): ONE "
       "motion keeps every weld satisfied and every other motion stretches a "
       "joint. So the disqualifier bites only where the variety is big enough "
       "to hide in, and here it is one-dimensional. R6 runs the single unit as "
       "the control, so this file is not quietly disagreeing with the record",
       nz == 7,
       f"zero modes on HC15: {nz} of {7 * n} (6 rigid-body + 1 expected)",
       "exactly 7"))

    # R3 -- and they are identified, not counted and assumed
    G = asm.globals(ctr)
    rigid = [outside(G[d], Zq) for d in range(6)]
    coh = outside(coherent_vector(n), Zq)
    span = np.column_stack([G.T, coherent_vector(n)])
    Sq, _ = np.linalg.qr(span)
    gap = float(np.linalg.norm(Zq - Sq @ (Sq.T @ Zq)))
    out["ident"] = (max(rigid), coh, gap)
    A(("R3  THE SEVEN ARE IDENTIFIED, NOT COUNTED AND ASSUMED -- six rigid-body "
       "motions and the COHERENT BREATHE. Each of the six global motions and "
       "the uniform-fold vector lands inside the measured zero space to 2e-15, "
       "and the seven together SPAN it. Counting alone would not have said "
       "which motions are free, and this project has been bitten before by a "
       "count that was right for the wrong reason. TWO-SIDED: a seventh mode "
       "that was something other than uniform folding would leave a residual "
       "here, and so would a zero space the seven failed to span",
       max(rigid) < 1e-12 and coh < 1e-12 and gap < 1e-10,
       f"worst rigid-body residual outside the zero space {max(rigid):.1e}; "
       f"coherent breathe {coh:.1e}; the seven span the zero space to "
       f"{gap:.1e}",
       "all seven inside to 1e-12, and spanning to 1e-10"))

    # R4 -- the spectrum itself
    A(("R4  NINETY-EIGHT REAL FREQUENCIES: THE ARRAY HAS A SPECTRUM. This is "
       "the first one in the project for a MEDIUM rather than for a single "
       "unit. Every mode except the seven free ones has a real, nonzero "
       "frequency, so the medium now rings: displace it off the coherent "
       "breathe and something pulls back. The stiffness k is free and every "
       "frequency scales as sqrt(k), so the RATIO of the extremes is the "
       "measurement and the absolute values are not",
       len(freqs) == 7 * n - 7 and freqs.min() > 1e-6
       and np.all(np.isfinite(freqs)),
       f"{len(freqs)} real frequencies, {freqs.min():.5f} .. "
       f"{freqs.max():.5f}, ratio {freqs.max() / freqs.min():.4f} (k = "
       f"{K_JOINT}, and every frequency scales as sqrt(k))",
       "7N - 7 real nonzero frequencies"))

    # R5 -- the two sides agree on which motion is free
    A(("R5  THE FREE MODE HERE IS THE FREE MOTION THERE, AND THE TWO FILES "
       "REACH IT FROM OPPOSITE DIRECTIONS. jb_mj R3 measured that a COHERENT "
       "drive produces ZERO contacts on this same patch, ever, because it "
       "keeps every weld shut. This file measures that the SAME motion is the "
       "one mode with zero frequency. Free of the contacts and free of the "
       "potential are the same freedom -- one measured with a hard wall and an "
       "LCP, the other with a quadratic and an eigensolve, sharing no code "
       "path. That agreement is worth a row because either result alone would "
       "be easy to believe for the wrong reason",
       coh < 1e-12 and nz == 7,
       f"coherent breathe lies in the zero space to {coh:.1e} here; jb_mj R3 "
       f"records 0 contacts under coherent drive against up to 96 bands under "
       f"a staggered one",
       "the same motion free under both couplings"))

    # R6 -- the control: reproduce the recorded death on a single unit
    one, _ = RC.honeycomb([(0, 0, 0)], gc=A_REF)
    ev1, Z1, _ = spectrum(one)
    A(("R6  THE CONTROL REPRODUCES qvf.2's ANALYTIC DEATH ON A SINGLE UNIT, "
       "which is what stops R2 from being a disagreement with the record. One "
       "cell has no welds at all, so V is identically zero and EVERY one of "
       "its degrees of freedom comes back free -- exactly the comment's "
       "conclusion that a soft joint 'stiffens only the directions TRANSVERSE "
       "to the variety'. The instrument therefore reports the death when the "
       "death is there, and R2's seven is a statement about the ARRAY. Without "
       "this row R2 would be one measurement contradicting a recorded proof, "
       "with no way to tell which was wrong",
       Z1.shape[1] == 7 and len(one.welds) == 0,
       f"single unit: {one.N} cell, {len(one.welds)} welds, "
       f"{Z1.shape[1]} of {7 * one.N} modes free",
       "every mode free on one unit, i.e. the recorded death reproduced"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_sj -- soften the joint, and the medium acquires a spectrum")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        f = out["freqs"]
        print(f"\n  SPECTRUM: {out['nz']} free modes, {len(f)} real "
              f"frequencies (k = {K_JOINT})")
        print(f"   lowest  {f[:6]}")
        print(f"   highest {f[-6:]}")

        print()
        print("  FOUR DECLARATIONS.")
        print("   KERNEL      DECLARED, and for the first time NOT")
        print("               inapplicable: a quadratic in the tied-vertex")
        print("               separation, stiffness k.")
        print("   MASS MODEL  DECLARED: corner point masses. The lamina peer")
        print("               is built in jb_mj and not swept here.")
        print("   METRIC      DECLARED: that model's mass matrix.")
        print("   PRIMITIVE   the tied VERTEX PAIR -- the physical joint.")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * THE MEDIUM HAS A SPECTRUM. Displace it off the coherent")
        print("     breathe and something pulls back. That is new, and it is")
        print("     what every earlier model lacked.")
        print("   * A DISPERSION RELATION IS NOW REACHABLE -- omega(k) by")
        print("     Bloch analysis on the honeycomb. DECISION 19 reformulated")
        print("     it away because it could not exist under a hard wall.")
        print("     THIS FILE DOES NOT COMPUTE ONE; it measures that the")
        print("     obstruction is gone. Do not read it as met.")
        print("   * NO ABSOLUTE FREQUENCY. k is free and every frequency goes")
        print("     as sqrt(k). Ratios only.")
        print("   * THE HARD WALL IS NOT A LIMIT OF THIS V. A quadratic well")
        print("     is not the stiff limit of a square well, so this does not")
        print("     supersede jb_ct or jb_sv -- it is a different joint law")
        print("     beside them. Bridging them needs the |d/t|^(2n) family.")
        print("   * THE FORM OF V IS THE SIMPLEST ONE. Real rubber is")
        print("     nonlinear and much stiffer in compression than extension.")
        print("     It is a modelling choice and it is not made here. The")
        print("     qvf.11 one-sidedness this bullet used to name as its")
        print("     target is WITHDRAWN -- DECISION 20, T2 [23713] -- and")
        print("     jb_ja built the asymmetry and found no lock at any ratio.")
        print("   * ONE PATCH, fifteen cells, free surfaces. Not the bulk.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
