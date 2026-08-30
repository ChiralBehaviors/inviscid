"""joint_exponent -- one covering, the joint exponent family re-measured: the chain
must be re-laid, its internal freedom is a different motion, and the
exponents survive anyway.

WHY THIS FILE EXISTS. DECISION 21 (T2 [23727]) rules the octahedral voids
EMPTY, and jb_1c measures what one covering does to omega(k): both sound
speeds move and the zone-corner degeneracy breaks. Ten modules still stand on
`jb_rc.honeycomb`, the double covering. jb_je is the first to re-measure
because it carries this project's only EXTERNAL validation -- the Hertzian
1/4 and 1/5 at n = 5/2 -- and because the similarity argument that makes its
exponents exact is about the JOINT LAW having no scale, not about the lattice.
If that argument is right the exponents survive any geometry; if the geometry
moves them, the argument was wrong. Either way it has to be measured, and it
was not.

THIS IS NOT A PORT, AND R1 MEASURES WHY. jb_je lays its chain on a BODY
DIAGONAL, `[(k, k, k)]`, which alternates VE and Oct. With the voids empty
only the all-even sites are solid, so those twelve sites hold SIX cells at
(0,0,0), (2,2,2), ... whose differences are (2,2,2) -- and the single
covering welds AXIS neighbours only, at (+-2,0,0) and permutations. Handed
jb_je's sites, `honeycomb_single` returns six cells and ZERO welds: the chain
falls into disconnected cells. The chain here is laid along an axis,
`[(2k, 0, 0)]`, where each bond carries TWO coincident vertex pairs across
the square opening rather than a triangular face's three.

THE INTERNAL FREEDOM IS A DIFFERENT MOTION, and this is the finding that
could have moved the exponents. Both chains have 84 DOF, constraint rank 66
and nullity 18 = 6 rigid + 12 internal. The integer is the same and the
mechanics is not:

    double (jb_je)   3-point pins on RIGID triangles, so a pin removes six
                     relative DOF and leaves every fold free: the 12
                     internal motions are 12 INDEPENDENT FOLDS.
    single (here)    the two shared points straddle the square OPENING, so
                     their in-cell distance depends on the fold. Each bond
                     spends one equation coupling adjacent folds and five on
                     relative rigid motion, leaving a HINGE about the line
                     through its two points: 12 folds + 11 hinges - 11
                     couplings = ONE COHERENT FOLD + 11 HINGES.

R2 identifies all eighteen on both sides: on the single chain exactly one
null direction carries any fold and in it every cell folds at the same rate,
the eleven constructed hinge rotations each lie in the null space, and the
6 + 1 + 11 span it exactly; on the double chain all twelve fold patterns are
free. So a fold-rate kick on cell 0 of the single chain is driving a medium
whose free internal motions are hinges, with the fold locked coherent, and
the front that propagates is not the same motion jb_je measured.

THE EXPONENTS SURVIVE, TO THE SAME FIGURES. The sweep, the fits and the rows
are jb_je's own code on the new chain -- `state`, `run`, `measure` and
`exponents` are imported, not rewritten -- so the only thing that changed is
the assembly.

    n        2      5/2      3       4       6       8    ->  wall
    p_v   0.0000  0.1993  0.3289  0.4883  0.6339  0.7008  ->  1    law (n-2)/n
    p_d   0.0000  0.2493  0.4959  0.9877  1.9511  2.9053  ->  diverges

    jb_je (double, body diagonal):
    p_v   0.0007  0.2003  0.3329  0.4908  0.6507  0.7078
    p_d   0.0007  0.2502  0.4996  0.9908  1.9795  2.9280

The Hertzian pair comes back -- 0.2493 against 1/4 and 0.1993 against 1/5 --
and n = 2 is amplitude-independent to four zeros (jb_je's chain gave 0.0007;
this one gives 0.0000, and the speed varies by 0.000% over the eightfold
drive). This is the similarity argument doing what it says: a pure power
law has no intrinsic scale, so displacement-by-L, time-by-L^((2-n)/2)
carries solutions to solutions on ANY assembly, and the exponents are what
that forces. The GEOMETRY sets the speed; the JOINT LAW sets how the speed
scales, and only the second is an exponent.

ONE HONEST DIFFERENCE. At n = 6 and 8 the single chain sits FURTHER below
the law than jb_je's does -- 0.6339 against its 0.6507, 0.7008 against
0.7078 -- and R9 measures that this is the integrator, on this chain: halving
the step moves p_v from 0.6339 to 0.6481 toward 2/3 and cuts the energy
drift from 4.1e-3 to 5.5e-4. The stiffer the joint the more a fixed step
under-resolves it, and the hinge-carried front is evidently the harder one to
resolve. The bias has the same sign as jb_je R7's, so R8 stays one-sided on
it; what is NOT claimed is that the large-n values here are as close to the
law as jb_je's, because they are not.

WHAT THE GEOMETRY DOES SET. At n = 2 and the same k, the single chain's
front runs at 10.5838 cells per unit time against the double chain's 14.2653:
a ratio of 0.7419 in cells, 0.8567 in length (the axis spacing is 2L, the body
diagonal's sqrt(3) L). Slower, in both units. The chain is not the bulk, so
this is not one of jb_1c's sound-speed ratios (2/sqrt6 = 0.8165 and
sqrt3/2 = 0.8660) and no closed form is claimed for it; it is recorded
because k is a convention and a RATIO at fixed k is the one speed statement
that is not.

SCOPE.
  * ONE CHAIN, twelve cells along an axis, free ends. Not the bulk. jb_je's
    chain is left untouched and its numbers are quoted, not re-run, except
    for the one n = 2 speed R10 needs.
  * SMOOTH JOINT, no LCP, no contact anywhere. The contact line (jb_ct,
    jb_sv, jb_tr) is still unmeasured under one covering.
  * THE WALL IS APPROACHED, NOT REACHED, exactly as in jb_je. R8 is gated
    one-sided on the same grounds and R9 re-measures those grounds HERE
    rather than citing them: the step bias is measured on this chain.
  * k IS A CONVENTION. R6 measures the exponent at two couplings 23x apart.
  * POINT-MASS MODEL, jb_rc's, via jb_mj. Harmonic in nothing.
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model.double_covering import joint_exponent as JE
from analysis.model import kinematics as MJ
from analysis.model import assembly as RC
A_REF = MJ.A_REF
CELLS = JE.CELLS
KICKS = JE.KICKS
K2 = JE.K2

#: The single covering's bond is an AXIS step of two site units; the double
#: covering's is a body-diagonal step. Both in units of the lattice constant.
SPACING_SINGLE = 2.0
SPACING_DOUBLE = np.sqrt(3.0)


def chain(ncells=CELLS, gc=A_REF):
    """The chain re-laid along an axis: the only line the single covering
    connects."""
    asm, _ = RC.honeycomb_single([(2 * k, 0, 0) for k in range(ncells)], gc=gc)
    return asm


def diagonal(ncells=CELLS, gc=A_REF):
    """jb_je's own sites handed to the single builder -- what DECISION 21
    leaves of them."""
    asm, _ = RC.honeycomb_single([(k, k, k) for k in range(ncells)], gc=gc)
    return asm


def _fold_rows(asm):
    return [7 * k + 6 for k in range(asm.N)]


def null_space(asm):
    """(C, rank, Z, centres) at q0: Z's columns span the constraint null space."""
    q = asm.q0()
    ctr, R, gam, B = asm.frames(q)
    C = asm.constraint_jacobian(asm.cell_jacobians(ctr, R, B))
    _, s, Vt = np.linalg.svd(C)
    rank = int((s > 1e-9).sum())
    return C, rank, Vt[rank:].T, ctr


def rigid_vectors(asm, ctr):
    """Three translations and three rotations about the origin, in the
    (cdot, omega, gdot) coordinates."""
    vs = []
    for d in range(3):
        e = np.eye(3)[d]
        t = np.zeros((asm.N, 7))
        t[:, 0:3] = e
        vs.append(t.ravel())
        r = np.zeros((asm.N, 7))
        r[:, 0:3] = np.cross(e, ctr)
        r[:, 3:6] = e
        vs.append(r.ravel())
    return vs


def hinge_vectors(asm, ctr):
    """For each bond, the rigid rotation of every cell beyond it about the
    line through the bond's two shared points. Zero fold everywhere."""
    X = asm.positions(asm.q0())
    vs = []
    for (k, l, pairs) in asm.welds:
        P = [X[k][a] for (a, _b) in pairs]
        n = P[1] - P[0]
        n = n / np.linalg.norm(n)
        h = np.zeros((asm.N, 7))
        for j in range(asm.N):
            if j >= l:
                h[j, 0:3] = np.cross(n, ctr[j] - P[0])
                h[j, 3:6] = n
        vs.append(h.ravel())
    return vs


def gate():
    checks, out = [], {}
    A = checks.append
    sgl = chain()
    ps = MJ.tied_pairs(sgl)
    out["cells"], out["pairs"] = sgl.N, len(ps)

    # ---- R1: the chain must be re-laid -------------------------------------
    dg = diagonal()
    ar = sorted({len(w[2]) for w in sgl.welds})
    _, degs = RC.honeycomb_single([(2 * k, 0, 0) for k in range(CELLS)],
                                  gc=A_REF)
    res = {}
    for g in (0.0, -15.0, -30.0, -45.0, -60.0):
        a2 = chain(gc=g)
        res[g] = float(np.abs(a2.weld_residual(a2.q0())).max())
    out["R1"] = (dg.N, len(dg.welds), sgl.N, len(sgl.welds), ar, res)
    A(("R1  jb_je's CHAIN FALLS APART UNDER ONE COVERING AND MUST BE RE-LAID "
       "ALONG AN AXIS, WHERE EVERY BOND IS A TWO-PAIR WELD. jb_je's twelve "
       "body-diagonal sites alternate VE and Oct; with the voids empty six of "
       "them are solid and they differ by (2,2,2), which the single covering "
       "does not weld -- so the builder returns six cells and no welds at "
       "all. Along an axis it returns twelve cells in a line, eleven welds of "
       "exactly two coincident pairs, end cells of degree one and interior "
       "cells of degree two, and the pair set read at a = -30 has zero weld "
       "residual at every phase of the exchange. TWO-SIDED on both chains: a "
       "body diagonal that welded would mean the rule built is not DECISION "
       "21's, and an axis bond with three pairs would mean a triangular face "
       "is being welded where there is only an opening",
       dg.N == 6 and len(dg.welds) == 0 and sgl.N == 12
       and len(sgl.welds) == 11 and ar == [2]
       and sorted(degs) == [1, 1] + [2] * 10 and max(res.values()) < 1e-12,
       f"body diagonal: {dg.N} cells, {len(dg.welds)} welds; axis: {sgl.N} "
       f"cells, {len(sgl.welds)} welds, arity {ar}, degrees "
       f"{sorted(degs)}; held-pair residual "
       + ", ".join(f"a={g:g}: {r:.1e}" for g, r in sorted(res.items())),
       "6 cells / 0 welds against 12 / 11 of arity 2, zero residual at every "
       "phase"))

    # ---- R2: the internal freedom is a different motion --------------------
    C, rank, Z, ctr = null_space(sgl)
    nul = 7 * sgl.N - rank
    Zf = Z[_fold_rows(sgl), :]
    Uf, sf, Vtf = np.linalg.svd(Zf)
    fold_rank = int((sf > 1e-9).sum())
    ones = np.ones(sgl.N) / np.sqrt(sgl.N)
    dir_err = min(float(np.linalg.norm(Uf[:, 0] - ones)),
                  float(np.linalg.norm(Uf[:, 0] + ones)))
    hv = hinge_vectors(sgl, ctr)
    hinge_err = max(float(np.linalg.norm(C @ h) / np.linalg.norm(h))
                    for h in hv)
    coherent = Z @ Vtf[0]
    span = int(np.linalg.matrix_rank(
        np.column_stack(rigid_vectors(sgl, ctr) + [coherent] + hv), tol=1e-9))
    dbl = JE.chain()
    Cd, rankd, Zd, _ = null_space(dbl)
    nuld = 7 * dbl.N - rankd
    fold_rank_d = int((np.linalg.svd(Zd[_fold_rows(dbl), :], compute_uv=False)
                       > 1e-9).sum())
    out["R2"] = (nul, fold_rank, dir_err, hinge_err, span, nuld, fold_rank_d)
    A(("R2  THE SAME NULLITY HIDES A DIFFERENT MECHANISM: ONE COHERENT FOLD "
       "AND ELEVEN HINGES, AGAINST jb_je's TWELVE FREE FOLDS. Both chains "
       "have 84 DOF, rank 66, nullity 18 = 6 rigid + 12 internal. On the "
       "double chain a bond pins three points of a RIGID triangle, so it "
       "removes six relative DOF and leaves both folds alone: all twelve "
       "fold patterns are free. On the single chain the two shared points "
       "straddle the square opening and their in-cell distance depends on "
       "the fold, so each bond spends one equation tying adjacent folds "
       "together and leaves a hinge about the line through its points. "
       "Measured, not counted: exactly ONE null direction carries any fold "
       "and in it every cell folds at the same rate; the eleven constructed "
       "hinge rotations each satisfy every constraint; and 6 + 1 + 11 span "
       "the null space exactly. This is the motion the exponents had to "
       "survive. TWO-SIDED: a second fold direction, a hinge outside the "
       "null space, or a double chain with fewer than twelve free folds all "
       "fail here",
       nul == 18 and fold_rank == 1 and dir_err < 1e-9 and hinge_err < 1e-9
       and span == 18 and nuld == 18 and fold_rank_d == 12,
       f"single: nullity {nul}, fold-carrying directions {fold_rank}, "
       f"coherent to {dir_err:.1e}, worst hinge residual {hinge_err:.1e}, "
       f"6 rigid + 1 fold + {len(hv)} hinges span {span}; double: nullity "
       f"{nuld}, fold-carrying directions {fold_rank_d}",
       "18 = 6 + 1 + 11 on the single chain, 18 = 6 + 12 folds on the double"))

    # ---- R3: the model, on the two-pair layout -----------------------------
    rng = np.random.default_rng(20260830)
    q = RC.apply_increment(sgl, sgl.q0(), 1e-3 * rng.standard_normal(7 * sgl.N))
    u = 1e-2 * rng.standard_normal((sgl.N, 7))
    a_f, s_f, M_f = JE.state(sgl, q, u, ps, 2.0, K2)
    J, M_r, Minv = MJ.kinematics(sgl, q, False)
    a_free = MJ.free_accel(sgl, q, u, J, Minv, False)
    s_r = MJ.separations(sgl, q, ps)
    N = MJ.band_rows(sgl, q, J, ps)
    f_ref = -(N.T @ (K2 * s_r)).reshape(sgl.N, 7)
    a_ref = a_free + np.einsum('kij,kj->ki', Minv, f_ref)
    kin_err = max(float(np.abs(a_f - a_ref).max()),
                  float(np.abs(s_f - s_r).max()))
    Cq = sgl.constraint_jacobian(J)
    g = sgl.weld_residual(q)
    sj_err = float(np.abs(f_ref.ravel() + K2 * (Cq.T @ g)).max())
    out["R3"] = (kin_err, sj_err, sgl.nc)
    A(("R3  jb_je's FUSED EVALUATION IS THE ESTABLISHED KINEMATICS ON THE "
       "TWO-PAIR LAYOUT TOO, AND AT n = 2 ITS FORCE IS THE GRADIENT OF THE "
       "VARIABLE-ARITY CONSTRAINT SET. jb_je R1 asserted both on a "
       "three-pair chain. Here the assembly has 66 constraint rows laid out "
       "by jb_1c's running-sum offsets, and the same two things are "
       "asserted against it at a randomly perturbed configuration: the fast "
       "path agrees with MJ.free_accel, MJ.separations and MJ.band_rows to "
       "1e-13, and the joint force equals -k C^T C where C is THIS "
       "assembly's constraint Jacobian. Without this row the sweep below "
       "could be integrating a force that disagrees with the constraint set "
       "jb_1c gated. TWO-SIDED: a layout that mis-indexed one pair would "
       "show in sj_err at once",
       kin_err < 1e-13 and sj_err < 1e-12 and sgl.nc == 66,
       f"fused vs MJ kinematics {kin_err:.2e}; joint force vs -k C^T C on "
       f"{sgl.nc} rows {sj_err:.2e}",
       "both under 1e-12 on a 66-row constraint set"))

    # ---- the sweep the next four rows read ---------------------------------
    fam = {}
    for n in (2.0, 2.5, 3.0, 4.0):
        e = JE.exponents(sgl, ps, n)
        if e is not None:
            fam[n] = e
    out["fam"] = fam

    # ---- R4: n = 2 is amplitude independent --------------------------------
    e2 = fam.get(2.0)
    vs2 = [r["v"] for r in e2["rows"]] if e2 else []
    spread2 = (max(vs2) / min(vs2) - 1.0) if vs2 else 9.9
    A(("R4  AT n = 2 THE FRONT SPEED IS AMPLITUDE-INDEPENDENT ON THE SINGLE "
       "CHAIN, SO THE p = 0 END OF THE BRIDGE IS STILL THERE. Same drive "
       "range as jb_je R2, eightfold, same threshold, a different internal "
       "mechanism carrying the front, and the speed does not move. "
       "TWO-SIDED: any drift with amplitude fails the row",
       e2 is not None and abs(e2["p_v"]) < 0.01 and spread2 < 0.01,
       f"p_v = {e2['p_v']:+.4f} at n = 2 (law 0); speed varies "
       f"{100 * spread2:.3f}% over {KICKS[-1] / KICKS[0]:.0f}x drive, "
       f"{min(vs2):.4f}..{max(vs2):.4f} cells per unit time"
       if e2 else "n = 2 sweep incomplete",
       "p_v = 0 within 0.01 and under 1% speed variation"))

    # ---- R5: the external check, again -------------------------------------
    e25 = fam.get(2.5)
    A(("R5  THE HERTZIAN PAIR SURVIVES THE CHANGE OF COVERING. jb_je R3 got "
       "0.2502 against the literature's 1/4 and 0.2003 against its 1/5 on a "
       "chain of free folds; here the folds are locked coherent and the "
       "hinges are free, and the same two published numbers come back. "
       "This is the row DECISION 21 could have taken away, and it is the "
       "project's only external check, so it is gated to jb_je's own "
       "tolerance and no looser. TWO-SIDED: miss either and the family's "
       "independence from the lattice is refuted",
       e25 is not None and abs(e25["p_d"] - 0.25) < 0.015
       and abs(e25["p_v"] - 0.2) < 0.015,
       f"n = 5/2: strain exponent {e25['p_d']:.4f} against 1/4, velocity "
       f"exponent {e25['p_v']:.4f} against 1/5 (jb_je: 0.2502, 0.2003)"
       if e25 else "n = 5/2 sweep incomplete",
       "both within 0.015 of the published values"))

    # ---- R6: the family, and the coupling ----------------------------------
    err_v = {n: abs(e["p_v"] - (n - 2) / n) for n, e in fam.items()}
    err_d = {n: abs(e["p_d"] - (n - 2) / 2) for n, e in fam.items()}
    alt = JE.exponents(sgl, ps, 2.5, kj=JE.seed(2.5) / 23.0)
    out["alt"] = alt
    A(("R6  BOTH EXPONENT LAWS HOLD ACROSS n = 2, 5/2, 3, 4 ON THE SINGLE "
       "CHAIN, AND THE EXPONENT IS STILL NOT THE COUPLING. The similarity "
       "argument is assembly-agnostic -- scaling displacement by L and time "
       "by L^((2-n)/2) maps solutions to solutions of ANY assembly under a "
       "pure power-law joint -- and this row is where that claim meets a "
       "second assembly. n = 5/2 is re-measured at a coupling 23x softer, "
       "as jb_je R4 does. TWO-SIDED: an exponent off its law by 0.02, or "
       "one that moved with k, fails",
       len(fam) == 4 and max(err_v.values()) < 0.02
       and max(err_d.values()) < 0.02 and alt is not None
       and abs(alt["p_v"] - 0.2) < 0.02,
       "  ".join(f"n={n:g}: p_v={e['p_v']:+.4f}/{(n - 2) / n:+.4f} "
                 f"p_d={e['p_d']:+.4f}/{(n - 2) / 2:+.4f}"
                 for n, e in sorted(fam.items()))
       + (f"; n=5/2 at k/23 gives p_v={alt['p_v']:+.4f}" if alt
          else "; alt coupling incomplete"),
       "every measured exponent within 0.02 of its law, at both couplings"))

    # ---- R7: the identity ---------------------------------------------------
    ids = {n: 1.0 / e["p_v"] - 1.0 / e["p_d"]
           for n, e in fam.items() if n > 2.0}
    qs = {n: abs(e["q"] - 2.0 / n) for n, e in fam.items()}
    A(("R7  1/p_v - 1/p_d = 1 AND THE STRAIN-AGAINST-DRIVE EXPONENT IS 2/n, "
       "NEITHER OF WHICH ASSUMES A PREDICTED LAW. The identity is the "
       "travelling-wave kinematics v ~ V_s * delta and nothing else; jb_je's "
       "mutation probes showed it orthogonal to the formulas R5 and R6 "
       "compare against. If the hinge freedom had turned the front into "
       "something other than one travelling motion with one strain and one "
       "particle velocity, this is the row that would say so. TWO-SIDED",
       len(ids) == 3 and max(abs(v - 1.0) for v in ids.values()) < 0.05
       and max(qs.values()) < 0.02,
       "  ".join(f"n={n:g}: 1/p_v-1/p_d={v:.4f}"
                 for n, v in sorted(ids.items()))
       + "; strain-vs-drive "
       + "  ".join(f"n={n:g}: {e['q']:.4f}/{2 / n:.4f}"
                   for n, e in sorted(fam.items())),
       "the identity within 0.05 of 1, and 2/n within 0.02"))

    # ---- R8: the bridge -----------------------------------------------------
    big = {}
    for n in (6.0, 8.0):
        e = JE.exponents(sgl, ps, n)
        if e is not None:
            big[n] = e
    out["big"] = big
    series = sorted({**fam, **big}.items())
    pv = [e["p_v"] for _, e in series]
    pd = [e["p_d"] for _, e in series]
    mono = all(pv[i + 1] > pv[i] for i in range(len(pv) - 1))
    A(("R8  p_v STILL CLIMBS TOWARD THE WALL'S p = 1 WHILE p_d DIVERGES. The "
       "bridge jb_je R6 found runs on the single chain too: the velocity "
       "exponent rises monotonically through n = 6 and 8 toward the hard "
       "wall's measured 1, the strain exponent passes 1 at n = 4 and keeps "
       "going. GATED ONE-SIDED, as jb_je's is, because the step biases p_v "
       "downward -- and R9 measures that bias on THIS chain rather than "
       "borrowing jb_je's",
       mono and len(pv) == 6 and pv[-1] > 0.70 and pd[-1] > 2.85
       and abs(pv[0]) < 0.01,
       "  ".join(f"n={n:g}: p_v={e['p_v']:.4f} p_d={e['p_d']:.4f}"
                 for n, e in series)
       + "; wall p = 1 by jb_ct R3 / jb_tr R3 / jb_sv R3",
       "p_v rising monotonically past 0.70 toward 1 while p_d passes 2.85"))

    # ---- R9: what the step costs, here ---------------------------------------
    e_c = big.get(6.0)
    e_f = JE.exponents(sgl, ps, 6.0, h=JE.FINE)
    out["conv"] = (e_c, e_f)
    t_v = 2.0 / 3.0
    ok9 = (e_c is not None and e_f is not None
           and abs(e_f["p_v"] - t_v) < abs(e_c["p_v"] - t_v)
           and e_f["p_v"] > e_c["p_v"] and e_f["dE"] < e_c["dE"])
    A(("R9  THE LARGE-n SHORTFALL IS THE INTEGRATOR ON THIS CHAIN TOO, AND "
       "THE BIAS HAS THE SAME SIGN. Halving the step at n = 6 moves p_v UP "
       "toward 2/3 and drops the energy drift with it. Measured here rather "
       "than cited from jb_je R7 because the medium the step has to resolve "
       "is a different one -- hinges free, fold coherent -- and a bias with "
       "the opposite sign would have voided R8's one-sided gate",
       ok9,
       (f"n = 6 at h = {JE.STEP:.2e}: p_v = {e_c['p_v']:.4f}, drift "
        f"{e_c['dE']:.1e};  at h = {JE.FINE:.2e}: p_v = {e_f['p_v']:.4f}, "
        f"drift {e_f['dE']:.1e}  (law 0.6667)")
       if ok9 or (e_c and e_f) else "convergence run incomplete",
       "the finer step nearer the law, higher, and less drifty"))

    # ---- R10: single against double, at the same k ---------------------------
    pdb = MJ.tied_pairs(dbl)
    md = JE.measure(dbl, pdb, 2.0, JE.seed(2.0), KICKS[2])
    v_s = float(np.mean(vs2)) if vs2 else float("nan")
    v_d = md["v"] if md else float("nan")
    r_cells = v_s / v_d
    r_len = r_cells * SPACING_SINGLE / SPACING_DOUBLE
    worst_s = max(r["resid"] for r in e2["rows"]) if e2 else 9.9
    out["R10"] = (v_s, v_d, r_cells, r_len, worst_s, md["resid"] if md else 9.9)
    A(("R10 THE COVERING SETS THE SPEED AND THE JOINT LAW SETS THE EXPONENT, "
       "AND THE SINGLE CHAIN IS THE SLOWER ONE. At n = 2 and the same k the "
       "two fronts are both uniform and the single chain's runs slower in "
       "cells per unit time and in length, the axis bond being 2L against "
       "the body diagonal's sqrt(3) L. k is a convention so neither speed is "
       "physical; the RATIO at fixed k is the one speed statement that "
       "survives the convention, and it is recorded. The chain is not the "
       "bulk and this is not one of jb_1c's sound-speed ratios; no closed "
       "form is claimed. ONE-SIDED BY NATURE, a direction claim: a single "
       "chain that ran faster in either unit would fail",
       md is not None and worst_s < 0.1 and md["resid"] < 0.1
       and r_cells < 1.0 and r_len < 1.0,
       f"single {v_s:.4f} against double {v_d:.4f} cells per unit time, "
       f"ratio {r_cells:.4f} in cells, {r_len:.4f} in length; front delay "
       f"spread single {worst_s:.3f}, double {md['resid'] if md else 9.9:.3f} "
       f"cells (jb_1c sound-speed ratios 0.8165 and 0.8660, for the reader)",
       "both fronts uniform to 0.1 cell and the single chain slower in both "
       "units"))

    # ---- R11: controls --------------------------------------------------------
    quiet, _ = JE.run(sgl, ps, 3.0, JE.seed(3.0), 0.0, tmax=0.6)
    loose, _ = JE.run(sgl, ps, 3.0, 0.0, KICKS[-1], tmax=0.6)
    out["ctl"] = (len(quiet), len(loose))
    A(("R11 CONTROLS, ON THIS CHAIN. Undriven, no front forms. Coupling "
       "removed, no front forms: the kicked cell keeps its fold rate and "
       "nothing reaches a neighbour, so every exponent above is a property "
       "of the joint law on the two-pair welds and not of the drive, the "
       "index, or the hinge freedom",
       len(quiet) == 0 and len(loose) == 0,
       f"undriven: {len(quiet)} cells reached; k = 0: {len(loose)}; driven "
       f"with coupling: reaches cell {sgl.N - 1}",
       "no front in either control"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_1e -- one covering: the exponent family re-measured on the "
              "axis chain")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        print("\n  THE FAMILY, MEASURED ON THE SINGLE COVERING")
        print(f"   {'n':>5s} {'p_v':>9s} {'(n-2)/n':>9s} {'p_d':>9s} "
              f"{'(n-2)/2':>9s} {'strain':>8s} {'2/n':>7s}")
        for n, e in sorted({**out["fam"], **out.get("big", {})}.items()):
            print(f"   {n:5.2f} {e['p_v']:9.4f} {(n - 2) / n:9.4f} "
                  f"{e['p_d']:9.4f} {(n - 2) / 2:9.4f} {e['q']:8.4f} "
                  f"{2 / n:7.4f}")
        print(f"   {'wall':>5s} {1.0:9.4f} {1.0:9.4f} {'-':>9s} "
              f"{'diverges':>9s}      (jb_ct R3 / jb_tr R3 / jb_sv R3)")
        print("   jb_je, double covering, body diagonal:")
        print("      n      2      5/2      3       4       6       8")
        print("      p_v  0.0007  0.2003  0.3329  0.4908  0.6507  0.7078")
        print("      p_d  0.0007  0.2502  0.4996  0.9908  1.9795  2.9280")

        v_s, v_d, r_c, r_l, _, _ = out["R10"]
        print("\n  SINGLE AGAINST DOUBLE, n = 2, same k")
        print(f"   {'':22s} {'single (axis)':>16s} {'double (diag)':>16s}")
        print(f"   {'internal DOF':22s} {'1 fold + 11 hinge':>16s} "
              f"{'12 folds':>16s}")
        print(f"   {'weld arity':22s} {2:>16d} {3:>16d}")
        print(f"   {'bond length / L':22s} {SPACING_SINGLE:>16.4f} "
              f"{SPACING_DOUBLE:>16.4f}")
        print(f"   {'front, cells/time':22s} {v_s:>16.4f} {v_d:>16.4f}")
        print(f"   {'front, length/time':22s} "
              f"{v_s * SPACING_SINGLE:>16.4f} {v_d * SPACING_DOUBLE:>16.4f}")
        print(f"   {'ratio single/double':22s} {r_c:>16.4f} (cells)  "
              f"{r_l:.4f} (length)")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * THE EXPONENTS ARE PROPERTIES OF THE JOINT LAW, NOT THE")
        print("     COVERING. The Hertzian 1/4 and 1/5 and the bridge")
        print("     p_v = (n-2)/n -> 1 survive a chain whose internal")
        print("     freedom is a different motion. Similarity is why.")
        print("   * THE SPEED IS A PROPERTY OF THE COVERING. jb_je's")
        print("     absolute speeds, like jb_bz's, belong to the double")
        print("     covering; the ratio at fixed k is recorded, unexplained.")
        print("   * jb_je IS NOT RETRACTED and is not re-run. Its chain is")
        print("     what the double covering connects; this one is what the")
        print("     single covering connects. They are different chains.")
        print("   * THE CONTACT LINE IS STILL UNMEASURED under one covering.")
        print("     Nothing here has an LCP in it.")
        print("   * ONE CHAIN, twelve cells, free ends, point masses. k is a")
        print("     convention; the exponents are the measurement.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
