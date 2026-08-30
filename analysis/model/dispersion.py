"""dispersion -- one covering: the voids are empty, and what that costs omega(k).

OWNER DECISION 21, 2026-08-29 (T2 [23727]), settled against an animation frame
of the VE-to-octahedron transform: "empty, as this animation clearly shows. in
the jitterbug transform, the ves become octa and the voids become ves. they
exchange places."

THE EXCHANGE IS WHAT MAKES IT UNAMBIGUOUS. Across [0, -60] the SOLID cells run
VE -> octahedron while the VOIDS run octahedron -> VE. The shapes swap and the
OCCUPANCY does not, so ONE sublattice is solid at every phase and the other is
void at every phase. `jb_rc.honeycomb` puts a cell at EVERY site, which fills
the voids; since every triangular face is shared by one even cell and one odd
cell and BOTH carry a plate there, it draws every interior triangle TWICE.

    plates drawn / distinct positions:  HC9 72/64, box r=2 280/216,
    box r=3 728/512 (216 twice, 296 once, and the 296 are the surface)

The ratio goes to 2 in the bulk. R1 measures it on HC9 and shows
`honeycomb_single` draws 64 for 64.

WHAT CHANGES, AND IT IS NOT SUBTLE. With the voids empty a cell's NEAREST
sites are voids, so it welds to its six AXIS neighbours -- second-nearest, at
(+-2,0,0) and permutations -- and each carries TWO coincident vertex pairs
rather than a triangular face's three. One cell per primitive cell, 7 DOF, 3
bonds, 18 constraint rows, against the double covering's 2 cells, 14 DOF, 8
welds, 72 rows.

                        THIS FILE (single)        jb_bz (double)
  bands                 7                         14
  zero modes at Gamma   4                         4
    identified          3 translations + the      same
                        coherent breathe
  optical at Gamma      12/5 (x3)                 12/5 (x3), 18/5 (x3), 6 (x4)
  sound speeds          1/(2 sqrt2), 1/2          sqrt(3)/4, 1/sqrt(3)
  zone corner R         sqrt(3/5) (x3),           ALL 14 at sqrt(3)
                        sqrt(3) (x4)

THE GOLDSTONE SURVIVES, AND SO DOES ONE OPTICAL BRANCH. The four zero modes
and their identification are unchanged, and 12/5 comes back exactly. What
moves is both sound speeds, and what breaks is jb_bz R5's total degeneracy at
the zone corner -- three of the seven bands drop to sqrt(3/5). R4 also
reproduces `jb_rc` R5e/R5f's published "internal DOF = 1, the coherent
breathe" EXACTLY, so the medium's one free motion was never a double-cover
artifact.

A PREDICTION THAT FAILED, RECORDED BECAUSE IT NEARLY DECIDED THIS. Before
building, the fourteen Gamma modes were decomposed into COMMON (both basis
cells alike) and DIFFERENTIAL (cells moving relatively, which separates every
shared plate), and the differential ones -- 12/5 and 6 -- were read as the
double-cover artifacts, 18/5 as legitimate. THE DIRECT CONSTRUCTION CONTRADICTS
THAT ON TWO OF THREE: 12/5 SURVIVES and 18/5 VANISHES. Removing the second
basis cell does not delete modes from a fixed spectrum; it changes the whole
constraint structure, 18 rows on 7 DOF against 72 on 14, so the reduced
spectrum is not a SUBSET of the larger one and no projection inside the larger
model can predict it. Build the reduced model; do not infer it.

THE ARITY IS PHASE DEPENDENT AND THE HELD PAIR SET IS NOT, which is the one
thing that could have sunk this. The axis contact carries FOUR coincident
vertex pairs at a = 0, TWO through the interior and ONE at a = -60, because
the square is an OPENING that closes as the exchange passes through it rather
than a face. R2 measures that the pair set read at a = -30 has zero weld
residual at EVERY phase, so those two are the persistent joints and the extras
at a = 0 are momentary. This is the same reason `jb_w.HONEYCOMB_REF_PHASE`
exists, reached from the other covering.

SCOPE, and it is the important part of this file.
  * NOTHING IS RETRACTED AND NOTHING IS MIGRATED. `jb_rc.honeycomb` is left
    intact and every module still uses it. This file measures what the OTHER
    covering gives so the two can be compared module by module.
  * jb_bz IS NOT WRONG ABOUT ITS OWN COVERING. Its 14 bands are what a
    double-covered lattice has. Which covering the PHYSICAL array has is
    DECISION 21's business, not a defect in jb_bz's arithmetic.
  * TEN MODULES SIT ON THE DOUBLE COVERING -- jb_bz, jb_ct, jb_ja, jb_je,
    jb_lf, jb_mj, jb_pr, jb_sj, jb_sv, jb_tr. The transport and contact line
    has NOT been re-measured under one covering at all, and this file does not
    do it.
  * k_joint IS A CONVENTION as everywhere here; only RATIOS and closed forms
    are measurements.
  * HARMONIC, at a = -30, point masses via jb_rc. No contact anywhere.
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model.double_covering import dispersion as BZ
from analysis.model import kinematics as MJ
from analysis.model import assembly as RC
A_REF = BZ.A_REF
K_JOINT = BZ.K_JOINT

#: The three lattice directions, in LATTICE units (two site units each). With
#: the voids empty the solid sublattice is simple cubic of spacing 2, so these
#: are its primitive vectors and each carries one bond per primitive cell.
LATT = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

HC9 = [(1, 1, 1)] + [(x, y, z) for x in (0, 2) for y in (0, 2) for z in (0, 2)]


def _faces():
    f = {}
    for (fc, c), slot in RC.SLOT.items():
        f.setdefault(fc, [None, None, None])[c] = slot
    return f


def plate_positions(asm):
    """Canonical position of every plate of every cell, for counting cover."""
    X = asm.positions(asm.q0())
    F = _faces()
    g = {}
    for i in range(asm.N):
        for fc in range(8):
            key = tuple(sorted(tuple(round(float(v), 6) for v in X[i][s])
                               for s in F[fc]))
            g.setdefault(key, []).append((i, fc))
    return g


def periodic_cell(gc=A_REF, ref=RC.SINGLE_REF_PHASE):
    """(J, M, bonds) for ONE solid cell and its three axis bonds."""
    sites = [(0, 0, 0)] + [tuple(2 * np.array(e)) for e in LATT]
    asm, _ = RC.honeycomb_single(sites, gc=gc, ref=ref)
    X = asm.positions(asm.q0())
    bonds = []
    for n, e in enumerate(LATT, start=1):
        pr = tuple((a, b) for a in range(RC.NV) for b in range(RC.NV)
                   if np.linalg.norm(X[0][a] - X[n][b]) < 1e-9)
        bonds.append((e, pr))
    ctr, R, gam, B = asm.frames(asm.q0())
    J = asm.cell_jacobians(ctr, R, B)
    return J[0], asm.mass_blocks(J)[0], tuple(bonds)


def bands(kvec, cell, k_joint=K_JOINT):
    """omega(k): the SEVEN bands of the single covering."""
    J, M, bonds = cell
    rows = []
    for e, prs in bonds:
        ph = np.exp(1j * float(np.dot(kvec, e)))
        for (a, b) in prs:
            rows.append(J[3 * a:3 * a + 3] - ph * J[3 * b:3 * b + 3])
    C = np.vstack(rows)
    H = k_joint * (C.conj().T @ C)
    L = np.linalg.cholesky(M)
    A = np.linalg.solve(L, np.linalg.solve(L, H).conj().T).conj().T
    return np.sqrt(np.clip(np.linalg.eigvalsh((A + A.conj().T) / 2.0), 0.0, None))


def gate():
    checks, out = [], {}
    A = checks.append
    cell = periodic_cell()
    J, M, bonds = cell

    # ---- R1: one covering, measured against the double one -----------------
    dbl, _ = RC.honeycomb(HC9, gc=A_REF)
    sgl, _ = RC.honeycomb_single(HC9, gc=A_REF)
    gd, gs = plate_positions(dbl), plate_positions(sgl)
    dd = sum(1 for v in gd.values() if len(v) > 1)
    ds = sum(1 for v in gs.values() if len(v) > 1)
    out["R1"] = (8 * dbl.N, len(gd), dd, 8 * sgl.N, len(gs), ds)
    A(("R1  THE VOIDS ARE EMPTY, SO EVERY TRIANGLE IS DRAWN ONCE -- AND THE "
       "BUILDER THIS REPLACES DRAWS THE INTERIOR ONES TWICE. DECISION 21: "
       "through the exchange the solid cells run VE -> octahedron while the "
       "VOIDS run octahedron -> VE, so the shapes swap and the occupancy does "
       "not; one sublattice is solid at every phase. `jb_rc.honeycomb` puts a "
       "cell at every site and therefore fills the voids, and because each "
       "triangular face is shared by an even and an odd cell that both carry a "
       "plate there, it draws HC9's eight interior triangles twice -- 72 "
       "plates over 64 distinct positions, a ratio that goes to 2 in the bulk. "
       "`honeycomb_single` keeps the all-even sublattice and draws 64 for 64. "
       "TWO-SIDED, and both sides can fail: the double builder must show "
       "duplication and the single one must show none",
       8 * dbl.N == 72 and len(gd) == 64 and dd == 8
       and 8 * sgl.N == 64 and len(gs) == 64 and ds == 0,
       f"double: {8 * dbl.N} plates over {len(gd)} positions, {dd} doubled; "
       f"single: {8 * sgl.N} plates over {len(gs)} positions, {ds} doubled",
       "72/64 with 8 doubled, against 64/64 with none"))

    # ---- R2: the held pair set is phase-stable, the ARITY is not ----------
    res, arity = {}, {}
    for g in (0.0, -15.0, -30.0, -45.0, -60.0):
        a2, _ = RC.honeycomb_single(HC9, gc=g)
        res[g] = float(np.abs(a2.weld_residual(a2.q0())).max())
        b2, _ = RC.honeycomb_single([(0, 0, 0), (2, 0, 0)], gc=g, ref=g)
        arity[g] = len(b2.welds[0][2]) if b2.welds else 0
    out["R2"] = (res, arity)
    A(("R2  THE HELD PAIR SET IS PHASE-STABLE EVEN THOUGH THE ARITY IS NOT, "
       "WHICH IS THE ONE THING THAT COULD HAVE SUNK THIS COVERING. The axis "
       "contact carries FOUR distinct shared points at a = 0, TWO through the "
       "interior and ONE at a = -60 -- the square is an OPENING that closes as "
       "the exchange passes through it, not a face, so it has no stable arity. "
       "If the welds were re-read at each phase the constraint COUNT would "
       "change with the fold, which no model can carry. Read at a = -30 and "
       "HELD, the two pairs have zero weld residual at EVERY phase of the "
       "exchange: they are the PERSISTENT joints and the extras at a = 0 are "
       "momentary. The count is of DISTINCT POINTS and not of label pairs, "
       "which is not pedantry -- at a = -60 the cell's twelve labels occupy "
       "six positions, so FOUR label pairs name ONE joint there, and a builder "
       "that emitted each would quadruple that joint's constraint. These "
       "numbers reproduce `jb_w_honeycomb`'s own square census exactly. "
       "TWO-SIDED: a pair set that drifted would show a residual, and an arity "
       "that did NOT vary would mean the reference phase is unnecessary and "
       "this row measures nothing",
       max(res.values()) < 1e-12 and arity[0.0] == 4
       and arity[-30.0] == 2 and arity[-60.0] == 1,
       "held-set weld residual " + ", ".join(f"a={g:g}: {r:.1e}"
                                             for g, r in sorted(res.items()))
       + "; distinct shared points READ at each phase "
       + ", ".join(f"a={g:g}: {n}" for g, n in sorted(arity.items())),
       "zero residual at every phase, and arity 4 / 2 / 1 matching jb_w"))

    # ---- R3: the layout change is a no-op for three-pair welds -------------
    hc15, _ = RC.honeycomb(MJ.hc15_sites(), gc=A_REF)
    off_ok = all(hc15._woff[r] == 9 * r for r in range(len(hc15.welds)))
    q = hc15.q0()
    ctr, R, gam, B = hc15.frames(q)
    Jh = hc15.cell_jacobians(ctr, R, B)
    Ch = hc15.constraint_jacobian(Jh)
    rank = int(np.linalg.matrix_rank(Ch, tol=1e-9))
    out["R3"] = (off_ok, hc15.nc, rank)
    A(("R3  GENERALISING THE CONSTRAINT LAYOUT DID NOT MOVE THE DOUBLE "
       "COVERING. `Assembly` used to hardcode nine rows per weld and index "
       "them at 9*r, which is right for a triangular face's three pairs and "
       "cannot express the axis weld's two. It is now a running sum. For an "
       "all-three-pair assembly the offsets must come back EXACTLY 9*r and the "
       "row count and rank must be untouched, or every number ten other "
       "modules have published moved when this file was added. Checked on "
       "HC15, the patch jb_sj measures. TWO-SIDED: a layout that shifted by "
       "even one row fails here",
       off_ok and hc15.nc == 9 * len(hc15.welds) and rank == 98,
       f"HC15 offsets all 9*r: {off_ok}; nc = {hc15.nc} for "
       f"{len(hc15.welds)} welds; constraint rank {rank}, nullity "
       f"{7 * hc15.N - rank} (jb_sj publishes 7 zero modes)",
       "offsets 9*r, nc = 9 x welds, rank 98 giving nullity 7"))

    # ---- R4: internal DOF = 1, jb_rc's own published number ----------------
    s9, _ = RC.honeycomb_single(HC9, gc=A_REF)
    q9 = s9.q0()
    c9, R9, g9, B9 = s9.frames(q9)
    C9 = s9.constraint_jacobian(s9.cell_jacobians(c9, R9, B9))
    r9 = int(np.linalg.matrix_rank(C9, tol=1e-9))
    intern = 7 * s9.N - r9 - 6
    out["R4"] = (7 * s9.N, r9, intern)
    A(("R4  THE SINGLE COVERING HAS INTERNAL DOF = 1, WHICH IS jb_rc R5e/R5f's "
       "OWN PUBLISHED NUMBER FOR THIS MEDIUM. That result -- 'exactly ONE "
       "motion keeps every weld satisfied, the coherent breathe' -- is what "
       "jb_sj's whole escape from qvf.2 rests on, and it is reproduced here "
       "with half the cells and a different weld set. So the medium's one free "
       "motion was NEVER an artifact of drawing the triangles twice. TWO-SIDED "
       "in both directions: more internal freedom would mean the axis welds "
       "under-constrain the array, and less would mean they over-constrain it "
       "and the breathe had been welded shut",
       intern == 1,
       f"single covering HC9: {7 * s9.N} DOF, constraint rank {r9}, nullity "
       f"{7 * s9.N - r9} = 6 rigid + {intern} internal",
       "internal DOF exactly 1"))

    # ---- R5: seven bands, four zero modes, and they are identified --------
    G = np.zeros(3)
    w = bands(G, cell)
    nz = int((w < 1e-7).sum())
    rows = []
    for e, prs in bonds:
        for (a, b) in prs:
            rows.append((J[3 * a:3 * a + 3] - J[3 * b:3 * b + 3]).real)
    C0 = np.vstack(rows)
    sv = np.linalg.svd(C0, compute_uv=False)
    Z = np.linalg.svd(C0)[2][(sv > 1e-9).sum():].T
    Q, _ = np.linalg.qr(Z)
    br = np.zeros(7)
    br[6] = 1.0
    res_br = float(np.linalg.norm(br - Q @ (Q.T @ br)))
    res_tr = max(float(np.linalg.norm(np.eye(7)[d] - Q @ (Q.T @ np.eye(7)[d])))
                 for d in range(3))
    opt = float(w[-1])
    out["R5"] = (w, nz, res_br, res_tr)
    A(("R5  SEVEN BANDS, FOUR OF THEM ZERO AT GAMMA, AND THE FOUR ARE "
       "IDENTIFIED RATHER THAN COUNTED. Three lattice translations and the "
       "COHERENT BREATHE, each landing in the measured zero space to machine "
       "precision -- so the Goldstone mode survives the change of covering "
       "intact, which is the single most important thing this file can say "
       "about it. The surviving optical branch is threefold at omega^2 = 12/5 "
       "EXACTLY, the same closed form the double covering has. TWO-SIDED: a "
       "fourth zero mode that was something other than uniform folding would "
       "leave a residual here, and a fifth would mean the axis welds do not "
       "hold the lattice",
       len(w) == 7 and nz == 4 and res_br < 1e-12 and res_tr < 1e-12
       and abs(opt ** 2 - 2.4) < 1e-9,
       f"{len(w)} bands, {nz} zero at Gamma; breathe outside the zero space "
       f"{res_br:.1e}, worst translation {res_tr:.1e}; optical "
       f"omega^2 = {opt ** 2:.12f} against 12/5",
       "7 bands, 4 zero, all four identified, optical exactly 12/5"))

    # ---- R6: the sound speeds move, and to closed forms --------------------
    sl = {}
    for kk in (1e-3, 1e-4):
        s = np.sort(bands(np.array([kk, 0.0, 0.0]), cell))[:4] / kk
        sl[kk] = s
    s4 = sl[1e-4]
    lo, hi = float(np.mean(s4[:2])), float(np.mean(s4[2:]))
    out["R6"] = (sl, lo, hi)
    A(("R6  BOTH SOUND SPEEDS MOVE, AND BOTH LAND ON CLOSED FORMS. The double "
       "covering gives sqrt(3)/4 and 1/sqrt(3); one covering gives "
       "1/(2 sqrt2) and 1/2 -- factors of sqrt(3/2) and 2/sqrt(3). So the "
       "acoustic sector still has two doubly degenerate speeds and a linear "
       "dispersion, and the NUMBERS jb_bz publishes are properties of ITS "
       "covering rather than of the medium. That a numerical slope lands on a "
       "closed form is the same evidence here it was there. TWO-SIDED: a "
       "slope that drifted with k would mean there is no linear regime to "
       "quote",
       abs(lo - 1.0 / (2 * np.sqrt(2))) < 1e-5 and abs(hi - 0.5) < 1e-5,
       f"slopes at k = 1e-4: {np.round(s4, 9)}; means {lo:.9f} against "
       f"1/(2 sqrt2) = {1 / (2 * np.sqrt(2)):.9f} and {hi:.9f} against 1/2 "
       f"(jb_bz: {np.sqrt(3) / 4:.9f} and {1 / np.sqrt(3):.9f})",
       "1/(2 sqrt2) and 1/2 to 1e-5"))

    # ---- R7: the zone-corner degeneracy BREAKS ----------------------------
    wr = bands(np.array([np.pi] * 3), cell)
    wd = BZ.bands(np.array([np.pi] * 3), BZ.unit_cell())
    lo3 = float(np.mean(wr[:3]))
    hi4 = float(np.mean(wr[3:]))
    out["R7"] = (wr, wd)
    A(("R7  jb_bz's TOTAL DEGENERACY AT THE ZONE CORNER BREAKS UNDER ONE "
       "COVERING, AND BREAKS TO ANOTHER CLOSED FORM. jb_bz R5 measures all "
       "fourteen bands collapsing to sqrt(3) at R = (pi,pi,pi) and reads that "
       "as a symmetry statement about the lattice. With the voids empty, four "
       "of the seven bands sit at sqrt(3) and the other THREE drop to "
       "sqrt(3/5). So the collapse was a property of the double covering, not "
       "of the lattice -- and the part that survives is exact. This row is the "
       "sharpest single difference between the two coverings and the easiest "
       "to check against the record. TWO-SIDED: a total collapse here would "
       "mean the covering does not affect the corner at all",
       abs(lo3 - np.sqrt(0.6)) < 1e-9 and abs(hi4 - np.sqrt(3.0)) < 1e-9
       and float(np.ptp(wd)) < 1e-9,
       f"single at R: 3 bands at {lo3:.12f} (sqrt(3/5) = {np.sqrt(0.6):.12f}) "
       f"and 4 at {hi4:.12f} (sqrt(3) = {np.sqrt(3.0):.12f}); double at R: "
       f"all {len(wd)} within {float(np.ptp(wd)):.1e}",
       "3 at sqrt(3/5) and 4 at sqrt(3), against a total collapse"))

    # ---- R8: controls ------------------------------------------------------
    w0 = bands(G, cell, k_joint=0.0)
    w4 = bands(G, cell, k_joint=4.0)
    scale = float(np.max(np.abs(w4 / np.sqrt(4.0) - w))) if len(w) else 9.9
    out["R8"] = (int((w0 < 1e-12).sum()), scale)
    A(("R8  CONTROLS. With the joint stiffness set to zero every band is zero, "
       "so the spectrum is a property of the welds rather than of the "
       "kinematics or the eigensolver. And every frequency scales as "
       "sqrt(k_joint) to machine precision, which is what makes k_joint a "
       "CONVENTION here as everywhere in this project: no absolute frequency "
       "or speed above is physical, and the closed forms are RATIOS that "
       "survive the choice. Without the first control the closed forms could "
       "be arithmetic on an empty matrix",
       int((w0 < 1e-12).sum()) == 7 and scale < 1e-12,
       f"k_joint = 0: {int((w0 < 1e-12).sum())} of 7 bands zero; "
       f"k_joint = 4 divided by sqrt(4) matches k_joint = 1 to {scale:.1e}",
       "all seven zero without the welds, and exact sqrt(k) scaling"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_1c -- one covering: the voids are empty, and what it costs")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        w = out["R5"][0]
        print(f"\n  THE SEVEN BANDS AT GAMMA: {np.round(w, 6)}")
        print("\n  SINGLE AGAINST DOUBLE")
        print(f"   {'':22s} {'single':>22s} {'double (jb_bz)':>26s}")
        print(f"   {'bands':22s} {7:>22d} {14:>26d}")
        print(f"   {'zero modes at Gamma':22s} {4:>22d} {4:>26d}")
        print(f"   {'optical at Gamma':22s} {'12/5 (x3)':>22s} "
              f"{'12/5, 18/5, 6':>26s}")
        print(f"   {'sound speeds':22s} {'1/(2sqrt2), 1/2':>22s} "
              f"{'sqrt3/4, 1/sqrt3':>26s}")
        print(f"   {'zone corner R':22s} {'sqrt(3/5) x3, sqrt3 x4':>22s} "
              f"{'all 14 at sqrt3':>26s}")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * THE GOLDSTONE SURVIVES. Four zero modes, identified, and")
        print("     jb_rc R5e/R5f's internal DOF = 1 reproduced exactly. The")
        print("     medium's one free motion was never a double-cover artifact.")
        print("   * NOTHING IS RETRACTED AND NOTHING IS MIGRATED. jb_rc.honeycomb")
        print("     is intact and all ten modules still use it. jb_bz is not")
        print("     wrong about its own covering; which covering the PHYSICAL")
        print("     array has is DECISION 21's business.")
        print("   * THE TRANSPORT AND CONTACT LINE IS UNMEASURED under one")
        print("     covering -- jb_ct, jb_sv, jb_tr, jb_lf, jb_mj. Whether the")
        print("     sonic-vacuum exponents or jb_je's (n-2)/n survive is NOT")
        print("     known and is not guessed at here.")
        print("   * RATIOS ONLY. k_joint is a convention (R8).")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
