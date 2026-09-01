"""doubled_block -- G5: the double of a finite block (note section 3c),
built, decomposed, and raced against the free block and the torus. The
section-5 question "does the double lose the seven free-boundary modes?"
is answered NO -- and the reason is exact: the double's spectrum is the
free block's plus a pinned-boundary sector, Neumann plus Dirichlet. What
the double actually buys is CONVERGENCE: its spectrum sits about three
times closer to the Bloch bands than the free block's at every size
tried, rivalling the exact periodic reference by side 4.

THE QUESTION (note sections 3c, 5-G5, 8.6; bead inviscid-44u). Close a
block not into a torus but into its DOUBLE: a second copy glued along
the boundary plate-backs, closing every half-weld, leaving no boundary
and no rig. Does that closure lose the free block's seven zero modes and
converge to the Bloch bands fastest?

WHAT IS BUILT HERE. The double of the all-even box of side n: two copies
of the box assembly, coincident sheet by sheet; every interior axis weld
duplicated per sheet; every HALF-WELD -- an axis bond of the box whose
partner site lies outside -- closed by tying the boundary cell's own
bond-pair vertices to the same vertices of its twin (sheet A to sheet
B). The result has every cell at weld-degree six: no boundary. The torus
reference is the wrapped real-space operator: the free Hessian plus the
deleted boundary bonds re-added with periodically wrapped indices and
ghost-bond geometry -- validated against the union of the gated Bloch
bands over the commensurate k-grid before it is used for anything.

  R1  THE CONSTRUCTOR CLOSES EVERY HALF-WELD: 2 n^3 cells, 2 W_box + 6
      n^2 welds, every cell at degree six, weld residual ~1e-15 at the
      reference -- no boundary, no rig, and no ambient obstruction (the
      second sheet rides coincident, the branched-cover picture made a
      constraint system).
  R2  THE SEVEN SURVIVE, AND THE DECOMPOSITION IS EXACT: the double has
      exactly seven zero modes (six rigid + the coherent breathe of both
      sheets, identified and spanning), and its full spectrum is the
      multiset union of the FREE block's (the sheet-symmetric, Neumann
      sector) and a Dirichlet sector (free Hessian + twice the boundary
      vertex pinning) with NO zero modes. The double does not lose the
      free-boundary modes; it keeps them all and adds a pinned copy.
  R3  THE TORUS INSTRUMENT: the wrapped operator's spectrum equals the
      union of the Bloch bands over the commensurate k-grid to 1e-12,
      and its zero space is FOUR -- three translations plus the breathe,
      the rotations measurably outside -- the k = 0 kernel, not the free
      block's seven.
  R4  THE RACE: Kolmogorov distance from each spectrum to the dense-k
      Bloch density of states. The double sits ~3x closer than the free
      block at every side, runs at the exact torus's level from side 3
      on (edging ahead at side 4, just behind at 5), and wins
      even at matched cell count (double of side 4, 128 cells, beats the
      free box of side 5, 125 cells). Closing into the double removes
      the surface pollution the free boundary causes.

WHAT THIS LICENSES AND WHAT IT DOES NOT. The section-3c closure is
spectrally real and cheap to build; as a BOUNDARY CONDITION for finite
simulation it costs 2x the cells of the box it closes and converges like
the periodic reference while keeping the rigid modes a free simulation
needs. Nothing here selects it as THE closure: the screw sector (G6,
mapping-torus form per screw_prerequisites) is unmeasured, and the
choice among double / torus / free is a modeling decision this gate
informs, not makes.

T2: [23865], [23914]. Ref: su2_boundary_conditions.md section 8.
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model import block_spectrum as BS
from analysis.model import dispersion as OC
from analysis.model.double_covering import soft_joint_spectrum as SJ

A_REF = BS.A_REF
AXES6 = ((2, 0, 0), (0, 2, 0), (0, 0, 2), (-2, 0, 0), (0, -2, 0), (0, 0, -2))


def boxsites(side):
    return [(2 * x, 2 * y, 2 * z) for x in range(side)
            for y in range(side) for z in range(side)]


def boundary_halfwelds(side):
    """(cell index, bond-pair vertex ids) for every axis bond of the box
    whose partner site lies outside -- read from a ghost-padded assembly so
    the pair set is the model's own, never re-derived."""
    sites = boxsites(side)
    idx = {s: i for i, s in enumerate(sites)}
    ghosts = []
    for s in sites:
        for e in AXES6:
            t = tuple(int(v) for v in np.array(s) + e)
            if t not in idx and t not in ghosts:
                ghosts.append(t)
    big, _ = RC.honeycomb_single(sites + ghosts, gc=A_REF)
    ntot, out = len(sites), []
    for (k, l, pairs) in big.welds:
        kin, lin = k < ntot, l < ntot
        if kin == lin:
            continue
        out.append((k, tuple(a for (a, b) in pairs)) if kin
                   else (l, tuple(b for (a, b) in pairs)))
    return out


def double(side):
    """(double assembly, free box assembly): two coincident sheets, interior
    welds per sheet, every half-weld closed sheet-to-sheet."""
    asm, _ = RC.honeycomb_single(boxsites(side), gc=A_REF)
    n = asm.N
    welds = list(asm.welds) + [(k + n, l + n, p) for (k, l, p) in asm.welds]
    for (i, vids) in boundary_halfwelds(side):
        welds.append((i, i + n, [(v, v) for v in vids]))
    return RC.Assembly(np.concatenate([asm.gam0] * 2),
                       np.vstack([asm.ctr0] * 2), welds), asm


def _mass_chol(asm, J):
    M = asm.mass_blocks(J)
    Mf = np.zeros((7 * asm.N, 7 * asm.N))
    for i in range(asm.N):
        Mf[7 * i:7 * i + 7, 7 * i:7 * i + 7] = M[i]
    return np.linalg.cholesky(Mf)


def _eigs(H, L):
    A = np.linalg.solve(L, np.linalg.solve(L, H).T).T
    return np.linalg.eigvalsh((A + A.T) / 2.0)


def dirichlet_eigs(free, side):
    """The sheet-antisymmetric sector: the free Hessian plus twice the
    boundary bond-pair vertex pinning."""
    ctr, R, gam, B = free.frames(free.q0())
    J = free.cell_jacobians(ctr, R, B)
    C = free.constraint_jacobian(J)
    H = C.T @ C
    n7 = 7 * free.N
    for (i, vids) in boundary_halfwelds(side):
        for v in vids:
            row = np.zeros((3, n7))
            row[:, 7 * i:7 * i + 7] = J[i][3 * v:3 * v + 3]
            H += 2.0 * row.T @ row
    return _eigs(H, _mass_chol(free, J))


def torus_eigs(side, cell):
    """The wrapped real-space operator: every axis bond present, boundary
    bonds re-added with wrapped indices and ghost-bond geometry."""
    Jc, Mc, bonds = cell
    sites = boxsites(side)
    idx = {s: i for i, s in enumerate(sites)}
    n = len(sites)
    rows = []
    for s in sites:
        i = idx[s]
        for (e, prs) in bonds:
            t = tuple(int(v) for v in (np.array(s) + 2 * np.array(e)) % (2 * side))
            j = idx[t]
            for (a, b) in prs:
                r = np.zeros((3, 7 * n))
                r[:, 7 * i:7 * i + 7] += Jc[3 * a:3 * a + 3]
                r[:, 7 * j:7 * j + 7] -= Jc[3 * b:3 * b + 3]
                rows.append(r)
    C = np.vstack(rows)
    Mf = np.zeros((7 * n, 7 * n))
    for i in range(n):
        Mf[7 * i:7 * i + 7, 7 * i:7 * i + 7] = Mc
    return _eigs(C.T @ C, np.linalg.cholesky(Mf)), 7 * n


def kolmogorov(w, ref):
    w, ref = np.sort(w), np.sort(ref)
    grid = np.unique(np.concatenate([w, ref]))
    cw = np.searchsorted(w, grid, side="right") / len(w)
    cr = np.searchsorted(ref, grid, side="right") / len(ref)
    return float(np.abs(cw - cr).max())


def gate():
    checks = []
    A = checks.append
    cell = OC.periodic_cell()
    sides = (2, 3, 4, 5)

    # ---- R1: the constructor closes every half-weld ---------------------------
    got1, ok1 = {}, True
    doubles = {}
    for n in sides:
        dbl, free = double(n)
        doubles[n] = (dbl, free)
        degs = np.zeros(dbl.N, int)
        for (k, l, _p) in dbl.welds:
            degs[k] += 1
            degs[l] += 1
        res = float(np.abs(dbl.weld_residual(dbl.q0())).max())
        wexp = 2 * len(free.welds) + 6 * n * n
        ok1 &= (dbl.N == 2 * n ** 3 and len(dbl.welds) == wexp
                and set(degs) == {6} and res < 1e-12)
        got1[n] = f"{dbl.N} cells, {len(dbl.welds)} welds, degrees {set(degs)}, res {res:.0e}"
    A(("R1  THE CONSTRUCTOR CLOSES EVERY HALF-WELD: two coincident sheets, interior "
       "welds per sheet, each of the box's 6n^2 boundary half-welds tied to the twin "
       "sheet's same vertices. Every cell of every double sits at weld-degree SIX -- "
       "no boundary, no rig -- and the weld residual is machine zero at reference.",
       ok1, got1, "2n^3 cells, 2W+6n^2 welds, all degree 6, res < 1e-12"))

    # ---- R2: the seven survive; Neumann + Dirichlet exactly -------------------
    got2, ok2 = {}, True
    for n in (2, 3, 4):
        dbl, free = doubles[n]
        nz_d, _f, Zq_d, c_d = BS.spec(dbl)
        rg, co, _e, gap = BS.identify(dbl, Zq_d, c_d)
        ev_d, _Z, _c = SJ.spectrum(dbl)
        ev_f, _Zf, _cf = SJ.spectrum(free)
        ev_dir = dirichlet_eigs(free, n)
        nz_dir = int((ev_dir < 1e-9 * ev_dir.max()).sum())
        mis = float(np.abs(np.sort(ev_d)
                           - np.sort(np.concatenate([ev_f, ev_dir]))).max())
        ok2 &= (nz_d == 7 and rg < 1e-10 and co < 1e-10 and gap < 1e-8
                and nz_dir == 0 and mis < 1e-10)
        got2[n] = (f"zeros {nz_d} (rigid {rg:.0e}, breathe {co:.0e}, span {gap:.0e}); "
                   f"Dirichlet zeros {nz_dir}; union mismatch {mis:.0e}")
    A(("R2  THE SEVEN SURVIVE, AND THE DECOMPOSITION IS EXACT: the double's zero "
       "space is exactly seven -- six rigid plus the coherent breathe of BOTH sheets, "
       "identified and spanning -- and its full spectrum is the multiset union of the "
       "free block's (the sheet-symmetric Neumann sector) and a Dirichlet sector "
       "(boundary vertices pinned, NO zeros). Section 5-G5's 'does the double lose "
       "the seven?' is NO: it keeps the free spectrum whole and adds a pinned copy.",
       ok2, got2, "7 identified; Dirichlet 0; union exact"))

    # ---- R3: the torus instrument ---------------------------------------------
    got3, ok3 = {}, True
    for n in (2, 3, 4):
        ev_t, _dof = torus_eigs(n, cell)
        ks = 2.0 * np.pi * np.arange(n) / n
        bl = np.sort(np.concatenate(
            [OC.bands(np.array([kx, ky, kz]), cell) ** 2
             for kx in ks for ky in ks for kz in ks]))
        mis = float(np.abs(np.sort(ev_t) - bl).max())
        nz_t = int((ev_t < 1e-9 * ev_t.max()).sum())
        ok3 &= (mis < 1e-10 and nz_t == 4)
        got3[n] = f"vs bands {mis:.0e}, zeros {nz_t}"
    z0 = int((OC.bands(np.zeros(3), cell) < 1e-8).sum())
    ok3 &= (z0 == 4)
    A(("R3  THE TORUS INSTRUMENT: the wrapped real-space operator's spectrum equals "
       "the union of the gated Bloch bands over the commensurate k-grid, and its "
       "zero space is FOUR -- the k = 0 kernel (three translations + the breathe; "
       "bands(0) has exactly four zeros; no rotations) -- not the free block's seven.",
       ok3, got3 | {"bands(0) zeros": z0}, "all < 1e-10; 4 zeros everywhere"))

    # ---- R4: the race ----------------------------------------------------------
    NREF = 20
    ks = 2.0 * np.pi * np.arange(NREF) / NREF
    ref = np.concatenate([OC.bands(np.array([kx, ky, kz]), cell)
                          for kx in ks for ky in ks for kz in ks])
    D = {}
    for n in sides:
        dbl, free = doubles[n]
        wf = np.sqrt(np.clip(SJ.spectrum(free)[0], 0, None))
        wd = np.sqrt(np.clip(SJ.spectrum(dbl)[0], 0, None))
        wt = np.sqrt(np.clip(torus_eigs(n, cell)[0], 0, None))
        D[n] = (kolmogorov(wf, ref), kolmogorov(wd, ref), kolmogorov(wt, ref))
    ok4 = (all(D[n][1] < D[n][0] / 2 for n in sides)
           and D[4][1] < D[5][0]          # matched cells: double(4)=128 vs free(5)=125
           and D[4][1] < D[4][2])
    A(("R4  THE RACE, AS KOLMOGOROV DISTANCE TO THE DENSE-k BLOCH DENSITY OF "
       "STATES: the double sits well under HALF the free block's distance at every "
       "side, runs at the exact periodic reference's level from side 3 on (ahead at "
       "4, just behind at 5), and wins at matched "
       "cell count -- the double of side 4 (128 cells) beats the free box of side 5 "
       "(125 cells). Closing into the double removes the free surface's spectral "
       "pollution.",
       ok4,
       {n: f"free {f:.4f} / double {d:.4f} / torus {t:.4f}" for n, (f, d, t) in D.items()},
       "double < free/2 at every side; double(4) < free(5); double(4) < torus(4)"))
    return checks


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("doubled_block -- G5: the double keeps the seven, splits Neumann + "
          "Dirichlet, and out-converges the free block")
    print("=" * 78)
    with np.errstate(all="ignore"):
        checks = gate()
    bad = 0
    for name, ok, got, want in checks:
        bad += 0 if ok else 1
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")
    print()
    print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
    print("   * THE SECTION-3c CLOSURE IS BUILDABLE AND SPECTRALLY REAL: no")
    print("     boundary, seven zero modes kept, spectrum = Neumann + Dirichlet,")
    print("     convergence to the bands at the periodic reference's level for")
    print("     twice the cells of the box it closes (par at matched count).")
    print("   * IT DOES NOT SELECT THE CLOSURE: the screw sector (G6, in the")
    print("     mapping-torus form screw_prerequisites licenses) is unmeasured,")
    print("     and double / torus / free is a modeling decision this gate")
    print("     informs, not makes.")
    print()
    print("  ALL CHECKS PASSED." if not bad
          else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
