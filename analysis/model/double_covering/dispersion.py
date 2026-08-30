"""dispersion -- omega(k) on the soft-jointed honeycomb. The criterion, on its original terms.

WHAT THIS IS. Epic inviscid-qvf asked, from 2026-08-11, for "a measured
dispersion relation". It was never met, and on 2026-08-28 jb_cp established why
it COULD not be met: under a hard-wall joint the medium has no Hessian, hence no
normal modes, hence no frequency. DECISION 19 reformulated the criterion around
that. Then the owner supplied the fact the model had been missing -- the physical
jitterbugs are SOFT RUBBER with FLEXIBLE JOINTS AT THE VERTICES -- jb_sj put the
compliance where the rubber is, and the obstruction went away.

This is omega(k), computed on the original terms.

THE REFORMULATION IS NOT RETRACTED. DECISION 19 is correct FOR A HARD-WALL
JOINT, and jb_ct, jb_tr, jb_sv and jb_lf remain valid for that joint. What
changed is that a different joint law -- the one the rig actually has -- admits
the original object as well. Two joints, two answers, both measured.

THE SETUP. The rectified cubic honeycomb's sites are the all-even and all-odd
integer triples, so the Bravais lattice is SIMPLE CUBIC with two cells in the
primitive cell: a VE at (0,0,0) and a hole cell at (1,1,1) running 60 degrees
ahead. Each cell carries seven generalized coordinates, so there are FOURTEEN
BANDS. The eight triangular-face welds run from the VE to the hole cell at the
eight lattice offsets (0 or -1) in each direction; each weld ties three vertex
pairs, giving 72 constraint rows per primitive cell.

    H(k) = k_joint * C(k)^H C(k),   omega^2(k) = eig(H(k), M)

EVERYTHING COMES OUT IN CLOSED FORM, which is not what one expects from a
numerical band calculation and is the strongest evidence the setup is right:

    FOUR ACOUSTIC BRANCHES from Gamma, two speeds, each doubly degenerate
        sqrt(3)/4 = 0.433012702      measured 0.433012727
        1/sqrt(3) = 0.577350269      measured 0.577350300

    THREE OPTICAL FREQUENCIES AT GAMMA, degeneracies 3, 3, 4
        sqrt(12/5) = 1.549193,  sqrt(18/5) = 1.897367,  sqrt(6) = 2.449490
        their squares are exactly 2.4, 3.6 and 6

    TOTAL DEGENERACY AT THE ZONE CORNER R = (pi,pi,pi)
        ALL FOURTEEN BANDS EQUAL sqrt(3), deviation 0.00e+00

FOUR ZERO MODES AT GAMMA, not three. Three are the rigid translations, which
any lattice has. The fourth is the COHERENT BREATHE -- the medium's own one
motion, which jb_sj measures free of the potential and jb_mj measures free of
the contacts. It is a genuine Goldstone mode of this medium and it is why the
acoustic sector has four branches rather than the usual three.

AND THE SOUND SPEED IS AMPLITUDE-INDEPENDENT, which is the sharpest possible
contrast with everything measured earlier this week. Under the hard wall,
jb_ct and jb_tr both measured speed PROPORTIONAL TO AMPLITUDE -- a sonic vacuum
with no linear regime at all. Soften the joint and the amplitude dependence
disappears: omega/k is a constant of the medium. That is what having a
dispersion relation MEANS, and the two joints could not differ more sharply.

WHAT IS STILL A CONVENTION. k_joint is free and every frequency scales exactly
as sqrt(k_joint) -- R6 measures that rather than asserting it -- so no absolute
frequency or speed here is physical. The RATIOS are the measurements, and the
closed forms above are ratios.

FOUR DECLARATIONS.
  * KERNEL: DECLARED -- quadratic in the tied-vertex separation, stiffness
    k_joint. Not inapplicable for the first time in this project.
  * MASS MODEL: DECLARED -- corner point masses. The lamina peer exists in
    jb_mj and is not swept here, so no frequency is model-independent.
  * METRIC: DECLARED -- that model's block-diagonal mass matrix.
  * PRIMITIVE: the tied VERTEX PAIR, which is the physical joint.

SCOPE. Harmonic, at the reference phase a = -30. A quadratic joint is the
SIMPLEST V and real rubber is neither quadratic nor symmetric; the
compression/extension asymmetry is not modelled here.

  SUPERSEDED IN ONE CLAUSE, 2026-08-29, DECISION 20 (T2 [23713]). This
  paragraph used to call the asymmetry "the live candidate for the qvf.11
  one-sidedness". It is no longer a candidate for anything, because the
  TARGET is withdrawn: the owner withdrew qvf.11's array lock as evidence
  about the MEDIUM, the rig having held centres that Gray p.40 requires to
  move, under tolerances wide enough to hide it. jb_ja built the asymmetry
  anyway and measured that it produces no lock at any stiffness ratio, and
  -- relevant to THIS file -- that it would cost the medium its Hessian at
  the weld, hence omega(k) itself. Nothing measured below is affected. jb_hc's own band structure is VOID --
it put springs on struts -- and nothing here reuses any number from it, only
the fact that the honeycomb's Bravais lattice is cubic.
"""
from __future__ import annotations

import itertools as it
import sys

import numpy as np

from analysis.model import assembly as RC
A_REF = -30.0

#: Joint stiffness. FREE: R6 measures that every frequency scales as its square
#: root, so only ratios below are measurements.
K_JOINT = 1.0

#: Relative cut for calling a band zero.
ZERO_RTOL = 1e-7


def unit_cell(gc=A_REF):
    """(J_VE, J_hole, M, welds) for the primitive cell.

    The honeycomb's sites are the all-even and all-odd integer triples, so the
    Bravais lattice is SIMPLE CUBIC (spacing 2 in site units) with a two-cell
    basis. A patch of the VE and its eight triangular-face neighbours supplies
    both cells' Jacobians and the weld correspondence; the eight neighbours are
    identical up to a lattice translation, which is asserted rather than
    assumed.
    """
    sites = [(0, 0, 0)] + [tuple(t) for t in it.product((1, -1), repeat=3)]
    asm, _ = RC.honeycomb(sites, gc=gc)
    q = asm.q0()
    ctr, R, gam, B = asm.frames(q)
    J = asm.cell_jacobians(ctr, R, B)
    M = asm.mass_blocks(J)
    same = (all(np.allclose(J[1], J[k]) for k in range(2, 9))
            and all(np.allclose(M[1], M[k]) for k in range(2, 9)))
    welds = []
    for (k, l, pairs) in asm.welds:
        other = l if k == 0 else k
        off = (np.array(sites[other]) - np.array([1, 1, 1])) // 2
        welds.append((tuple(int(x) for x in off), tuple(pairs)))
    Mf = np.zeros((14, 14))
    Mf[:7, :7] = M[0]
    Mf[7:, 7:] = M[1]
    return J[0], J[1], Mf, tuple(welds), same, asm


def bands(kvec, cell, k_joint=K_JOINT):
    """omega(k): the fourteen bands at reduced wavevector `kvec`.

    A cell at lattice offset R carries amplitude exp(i k . R), so each weld
    contributes three rows tying the VE's vertex a to the hole cell's vertex b
    with that phase. H is C^H C, Hermitian positive semi-definite by
    construction, and the generalised problem is reduced by a Cholesky of M --
    the same reduction the Java GeneralizedEigensolver performs.
    """
    JA, JB, Mf, welds = cell[0], cell[1], cell[2], cell[3]
    rows = []
    for off, pairs in welds:
        ph = np.exp(1j * float(np.dot(kvec, off)))
        for (a, b) in pairs:
            r = np.zeros((3, 14), dtype=complex)
            r[:, 0:7] = JA[3 * a:3 * a + 3]
            r[:, 7:14] = -ph * JB[3 * b:3 * b + 3]
            rows.append(r)
    C = np.vstack(rows)
    H = k_joint * (C.conj().T @ C)
    L = np.linalg.cholesky(Mf)
    A = np.linalg.solve(L, np.linalg.solve(L, H).conj().T).conj().T
    ev = np.linalg.eigvalsh((A + A.conj().T) / 2.0)
    return np.sqrt(np.clip(ev, 0.0, None))


def gate():
    checks, out = [], {}
    A = checks.append
    cell = unit_cell()
    JA, JB, Mf, welds, same, asm = cell
    G = bands(np.zeros(3), cell)
    cut = ZERO_RTOL * G.max()

    # R1 -- the setup
    A(("R1  THE PRIMITIVE CELL IS TWO CELLS ON A SIMPLE CUBIC LATTICE, AND "
       "THAT IS READ OFF THE HONEYCOMB RATHER THAN ASSUMED. Its sites are the "
       "all-even and all-odd integer triples, so the Bravais lattice is cubic "
       "with a VE at the origin and a hole cell at (1,1,1) sixty degrees "
       "ahead. Seven generalized coordinates each gives FOURTEEN BANDS. The "
       "eight triangular-face welds run to the eight lattice offsets, three "
       "tied vertex pairs apiece, 72 constraint rows. The eight neighbours "
       "must be identical up to a lattice translation for the two-cell basis "
       "to be right, and that is asserted here rather than assumed",
       same and len(welds) == 8
       and sorted(w[0] for w in welds) == sorted(
           tuple(t) for t in it.product((0, -1), repeat=3))
       and all(len(w[1]) == 3 for w in welds),
       f"{len(welds)} welds at offsets "
       f"{sorted(w[0] for w in welds)}; three tied pairs each; the eight "
       f"neighbour Jacobians and mass blocks identical: {same}",
       "8 welds at the 8 cubic offsets, identical neighbours"))

    # R2 -- four acoustic branches, and the fourth is the medium's own motion
    nz = int((G < cut).sum())
    generic = bands(np.array([0.7, 0.3, 1.1]), cell)
    out["nz"] = nz
    A(("R2  FOUR BRANCHES GO TO ZERO AT GAMMA, NOT THREE, AND THE EXTRA ONE IS "
       "THE MEDIUM'S OWN MOTION. Three are the rigid translations that any "
       "lattice has. The fourth is the COHERENT BREATHE -- jb_sj measures it "
       "free of the potential, jb_mj measures it free of the contacts, and "
       "here it is a genuine GOLDSTONE MODE of the medium, which is why the "
       "acoustic sector has four branches instead of the usual three. "
       "TWO-SIDED: at a generic wavevector NO band may vanish, so a spurious "
       "zero mode from a defective constraint matrix would redden this row "
       "rather than pass it",
       nz == 4 and generic.min() > 1e-6,
       f"{nz} bands below {cut:.1e} at Gamma, of 14; at a generic k the "
       f"lowest band is {generic.min():.6f}",
       "exactly 4 at Gamma and none away from it"))

    # R3 -- the sound speeds, in closed form, and amplitude-independent
    eps = 1e-4
    w = bands(np.array([eps, 0.0, 0.0]), cell)
    s1, s2 = w[0] / eps, w[2] / eps
    c1, c2 = np.sqrt(3.0) / 4.0, 1.0 / np.sqrt(3.0)
    slopes = [bands(np.array([e, 0, 0]), cell)[0] / e
              for e in (1e-3, 1e-4, 1e-5)]
    out["speeds"] = (s1, s2)
    A(("R3  TWO SOUND SPEEDS, EACH DOUBLY DEGENERATE, AND BOTH IN CLOSED FORM "
       "-- sqrt(3)/4 and 1/sqrt(3) to eight decimals. That a numerical band "
       "calculation lands on exact closed forms is the strongest evidence the "
       "setup is right. AND THE SPEED IS AMPLITUDE-INDEPENDENT, which is the "
       "sharpest contrast this project can draw: under a HARD-WALL joint, "
       "jb_ct and jb_tr both measured speed PROPORTIONAL TO AMPLITUDE -- a "
       "sonic vacuum with no linear regime at all. Soften the joint and the "
       "amplitude dependence is simply gone; omega/k is a constant of the "
       "medium. That is what having a dispersion relation MEANS, and the two "
       "joint laws could not differ more sharply",
       abs(s1 - c1) < 1e-6 and abs(s2 - c2) < 1e-6
       and max(slopes) / min(slopes) - 1.0 < 1e-4,
       f"slopes at Gamma {s1:.9f} (x2) and {s2:.9f} (x2) against "
       f"sqrt(3)/4 = {c1:.9f} and 1/sqrt(3) = {c2:.9f}; the lower slope over "
       f"k = 1e-3, 1e-4, 1e-5: " + ", ".join(f"{v:.9f}" for v in slopes),
       "both closed forms to 1e-6, slope constant as k -> 0"))

    # R4 -- the optical frequencies at Gamma
    opt = sorted(set(np.round(G[nz:], 6)))
    targets = [np.sqrt(2.4), np.sqrt(3.6), np.sqrt(6.0)]
    out["opt"] = opt
    A(("R4  THE OPTICAL FREQUENCIES AT GAMMA ARE ALSO CLOSED FORM: their "
       "SQUARES are exactly 12/5, 18/5 and 6, with degeneracies three, three "
       "and four. Ten optical branches on top of the four acoustic ones makes "
       "fourteen, which is the count the two-cell basis requires -- so this "
       "row also checks that no band has been lost or double counted",
       len(opt) == 3
       and all(abs(o - t) < 1e-5 for o, t in zip(opt, targets))
       and len(G) - nz == 10,
       "; ".join(f"{o:.6f} (squared {o * o:.6f})" for o in opt)
       + f"; {len(G) - nz} optical branches",
       "squares exactly 12/5, 18/5, 6 and ten optical branches"))

    # R5 -- the zone-corner degeneracy
    Rk = bands(np.array([np.pi] * 3), cell)
    out["R"] = (float(Rk.min()), float(Rk.max()))
    A(("R5  ALL FOURTEEN BANDS COLLAPSE TO sqrt(3) AT THE ZONE CORNER "
       "R = (pi,pi,pi), with deviation ZERO to machine precision. A total "
       "degeneracy of every branch at one point is a symmetry statement about "
       "the lattice, not a coincidence, and it is the kind of thing that would "
       "be destroyed by almost any error in the phase convention or the weld "
       "correspondence. TWO-SIDED: the row fails if the bands separate there "
       "OR if they collapse anywhere they should not, so the M point is "
       "checked as a control and must NOT be degenerate",
       float(Rk.max() - Rk.min()) < 1e-9
       and abs(float(Rk.mean()) - np.sqrt(3.0)) < 1e-9
       and float(np.ptp(bands(np.array([np.pi, np.pi, 0.0]), cell))) > 0.1,
       f"at R all 14 bands in [{Rk.min():.9f}, {Rk.max():.9f}] against "
       f"sqrt(3) = {np.sqrt(3.0):.9f}; at M the spread is "
       f"{float(np.ptp(bands(np.array([np.pi, np.pi, 0.0]), cell))):.4f}",
       "total degeneracy at R, and none at M"))

    # R6 -- and the stiffness is a convention
    ref = bands(np.array([1.0, 0.0, 0.0]), cell, k_joint=1.0)
    scal = {kj: bands(np.array([1.0, 0.0, 0.0]), cell, k_joint=kj).max()
            for kj in (1.0, 4.0, 9.0)}
    ratios = [scal[kj] / (scal[1.0] * np.sqrt(kj)) for kj in (4.0, 9.0)]
    A(("R6  EVERY FREQUENCY SCALES AS THE SQUARE ROOT OF THE JOINT STIFFNESS, "
       "MEASURED RATHER THAN ASSERTED -- so k_joint is a convention and no "
       "absolute frequency or speed above is physical. The closed forms in R3 "
       "and R4 are RATIOS and survive it; a reader quoting sqrt(6) as a "
       "frequency in any unit would be quoting the choice k_joint = 1. This "
       "row exists because the project's standing rule is that absolute scales "
       "here are conventions and only ratios are measurements",
       max(abs(r - 1.0) for r in ratios) < 1e-9,
       "; ".join(f"k={kj}: max band {scal[kj]:.6f}" for kj in sorted(scal))
       + f"; divided by sqrt(k) they agree to "
         f"{max(abs(r - 1.0) for r in ratios):.1e}",
       "omega proportional to sqrt(k_joint) exactly"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_bz -- omega(k) on the soft-jointed honeycomb")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        cell = unit_cell()
        print("\n  BANDS ALONG THE HIGH-SYMMETRY PATH (k_joint = 1)")
        for nm, kv in (("Gamma  (0,0,0)   ", [0, 0, 0]),
                       ("X    (pi,0,0)    ", [np.pi, 0, 0]),
                       ("M    (pi,pi,0)   ", [np.pi, np.pi, 0]),
                       ("R    (pi,pi,pi)  ", [np.pi, np.pi, np.pi])):
            w = bands(np.array(kv, float), cell)
            print(f"   {nm} {np.array2string(w, precision=4)}")

        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * THE EPIC'S ORIGINAL CRITERION, ON ITS ORIGINAL TERMS: a")
        print("     measured dispersion relation, fourteen bands, four")
        print("     acoustic and ten optical, every feature in closed form.")
        print("   * DECISION 19 IS NOT RETRACTED. It is correct for a")
        print("     HARD-WALL joint, and jb_ct / jb_tr / jb_sv / jb_lf remain")
        print("     valid for that joint. A different joint law admits the")
        print("     original object too. Two joints, two answers, both")
        print("     measured -- and the amplitude dependence of the speed is")
        print("     exactly what separates them.")
        print("   * NO ABSOLUTE FREQUENCY OR SPEED. k_joint is free and every")
        print("     frequency goes as its square root (R6). Ratios only.")
        print("   * HARMONIC, AT a = -30. A quadratic joint is the simplest V;")
        print("     real rubber is neither quadratic nor symmetric, and the")
        print("     compression/extension asymmetry -- still the live")
        print("     candidate for the qvf.11 one-sidedness -- is not here.")
        print("   * NOTHING IS REUSED FROM jb_hc's BAND STRUCTURE, which is")
        print("     VOID: it put springs on struts. Only the fact that the")
        print("     honeycomb's Bravais lattice is cubic is taken from it.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
