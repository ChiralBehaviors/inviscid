"""block_spectrum -- one covering, the spectrum: jb_sj's patch is a tree with the voids
empty, the closed patches keep seven free modes at every size, and their
spectrum sits under jb_1c's bands and climbs to them.

THE QUESTION. jb_sj [23682] measures, on HC15 -- one VE and all fourteen of
its face neighbours, 15 cells, 32 welds, 105 DOF -- that a soft joint gives
the medium exactly SEVEN zero modes (six rigid + the coherent breathe) and
98 real frequencies, 0.46813 .. 2.08755. That is the file that cleared
qvf.2's "soft joints are dead analytically", and it is the foundation jb_bz
stands on. DECISION 21 (T2 [23727]) rules the voids EMPTY, so HC15 as built
draws its interior triangles twice, and the question the continuation asked
was: what happens to the 98?

THE ANSWER HAS TWO HALVES, AND THE FIRST IS A TRAP THIS PROJECT HAS ALREADY
PAID FOR. Under one covering HC15's fifteen sites keep SEVEN solid cells --
the centre and its six axis neighbours; the eight triangular neighbours are
voids -- and those seven form a STAR: six bonds, no cycle, a tree. Its
spectrum has THIRTEEN zero modes and 36 frequencies. R2 identifies all
thirteen: six rigid, the coherent breathe, and six hinges, one per leaf,
each a rigid rotation of that leaf about the line through its bond's two
shared points. That is the same-sites comparison, and it is a tree against a
cycle -- the confound recorded in [23746] lesson 4 (HC9's 9 against 1), met
again from the other covering. It measures the PATCH, not the medium, and
"the 98 become 36" is not a statement about anything physical.

THE SECOND HALF IS THE MEDIUM. The closed single-covering patches are the
all-even boxes: side 2 is jb_1c R4's cube (8 cells, 12 bonds), and the
series runs to side 6 (216 cells, 540 bonds, 1512 DOF). At EVERY size the
zero space is exactly seven and R3 identifies it as six rigid plus the
coherent breathe, spanning to 1e-10. So jb_sj's seven is what the single
covering has too, on every closed patch tried, and the hinges that a tree
leaves free are killed by the first cycle.

    patch             cells  bonds   DOF  zero  freqs   min      min*side   max     sqrt3-max
    HC15 double(jb_sj)   15     32   105     7     98  0.46813              2.08755   (ceiling sqrt6)
    HC15 single (star)    7      6    49    13     36  0.19762              1.40245
    box 2 (cube)          8     12    56     7     49  0.38730   0.77460    1.46277   0.26928
    box 3                27     54   189     7    182  0.32093   0.96280    1.59386   0.13819
    box 4                64    144   448     7    441  0.24573   0.98290    1.65003   0.08202
    box 5               125    300   875     7    868  0.19850   0.99248    1.67848   0.05357
    box 6               216    540  1512     7   1505  0.16628   0.99769    1.69446   0.03759

THE TOP IS BOUNDED BY jb_1c's BANDS AND CLIMBS TO THEM (R4). A free patch's
Hessian is the periodic operator with the boundary-crossing bonds deleted
and the mass matrix untouched, so by the Rayleigh quotient no finite patch
can ring above the Bloch maximum. jb_1c's seven bands top out at sqrt(3),
at the zone corner R; the boxes' maxima rise monotonically 1.463 -> 1.694
toward it and the gap shrinks at every step. The same theorem is checked on
the OTHER covering as the instrument's control: jb_sj's 2.08755 sits under
jb_bz's Bloch maximum sqrt(6), which is the omega^2 = 6 mode at Gamma --
the one jb_1c shows the single covering does not have. So the two coverings'
finite spectra are each bounded by their own bands, and the ceilings differ
by sqrt(2).

THE BOTTOM IS ACOUSTIC (R5). The lowest nonzero frequency falls as 1/side:
the log-log slope over sides 3..6 is -0.95 and min*side rises 0.963 ->
0.998 toward a bound. There is no gap above the seven: the medium's
lowest motions on a finite patch are long-wavelength sound, as a medium
with two acoustic speeds must have. (That min*side lands near 1.000 is
noted and NOT claimed as a closed form -- k is a convention and this is one
number at one k.)

WHAT THIS DOES AND DOES NOT DO TO jb_sj. jb_sj stands: its seven is
reproduced by the instrument here (R1) and is what the single covering has
on every closed patch. Its 98 frequencies and their range are properties of
a double-covered patch and are not migrated. Its R5 -- the coherent breathe
is free of the contacts too, jb_mj R3 -- is NOT re-measured here because
jb_mj's LCP has not been pointed at the single covering; that is the
contact line's work and it is still open.

SCOPE.
  * HARMONIC, at a = -30, k = 1 (a convention; R6 measures sqrt(k)
    scaling). Point masses via jb_rc. No contact anywhere.
  * FREE SURFACES on every patch. The bulk is jb_1c's; this file measures
    finite patches against it, not instead of it.
  * THE BLOCH MAXIMA are taken on a grid over the irreducible wedge
    0 <= kz <= ky <= kx <= pi, which for both coverings lands exactly on a
    closed form at a symmetry point, so the grid is not the resolution.
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import dispersion as OC
from analysis.model.double_covering import dispersion as BZ
from analysis.model import kinematics as MJ
from analysis.model import assembly as RC
from analysis.model.double_covering import soft_joint_spectrum as SJ
A_REF = MJ.A_REF
SIDES = (2, 3, 4, 5, 6)

#: jb_sj's published HC15 numbers, the instrument check.
SJ_ZERO, SJ_FREQS, SJ_MIN, SJ_MAX = 7, 98, 0.46813, 2.08755


def box(side, gc=A_REF):
    """The all-even box of `side` cells on an edge: the closed patch."""
    sites = [(2 * x, 2 * y, 2 * z) for x in range(side)
             for y in range(side) for z in range(side)]
    asm, _ = RC.honeycomb_single(sites, gc=gc)
    return asm


def spec(asm, k=SJ.K_JOINT):
    """(zero count, frequencies, orthonormal zero-space basis, centres)."""
    ev, Z, ctr = SJ.spectrum(asm, k=k)
    nz = Z.shape[1]
    Zq, _ = np.linalg.qr(Z) if nz else (np.zeros((7 * asm.N, 0)), None)
    return nz, np.sqrt(np.clip(ev[nz:], 0.0, None)), Zq, ctr


def tree_hinges(asm):
    """One rigid rotation per bond about the line through its two shared
    points, carrying the side of the tree that does not contain cell 0."""
    X = asm.positions(asm.q0())
    ctr = asm.ctr0
    adj = {i: set() for i in range(asm.N)}
    for (k, l, _p) in asm.welds:
        adj[k].add(l)
        adj[l].add(k)
    vs = []
    for (k, l, pairs) in asm.welds:
        comp, stack = {l}, [l]
        while stack:
            c = stack.pop()
            for m in adj[c]:
                if {c, m} == {k, l} or m in comp:
                    continue
                comp.add(m)
                stack.append(m)
        P = [X[k][a] for (a, _b) in pairs]
        n = P[1] - P[0]
        n = n / np.linalg.norm(n)
        h = np.zeros((asm.N, 7))
        for j in comp:
            h[j, 0:3] = np.cross(n, ctr[j] - P[0])
            h[j, 3:6] = n
        vs.append(h.ravel())
    return vs


def identify(asm, Zq, ctr, extra=()):
    """Residuals of the rigid, coherent and `extra` vectors outside the zero
    space, and how well they all span it."""
    G = asm.globals(ctr)
    rigid = max(SJ.outside(G[d], Zq) for d in range(6))
    coh = SJ.outside(SJ.coherent_vector(asm.N), Zq)
    ex = max((SJ.outside(v, Zq) for v in extra), default=0.0)
    S = np.column_stack([G.T, SJ.coherent_vector(asm.N)] + list(extra))
    Sq, _ = np.linalg.qr(S)
    gap = float(np.linalg.norm(Zq - Sq @ (Sq.T @ Zq)))
    return rigid, coh, ex, gap


def bloch_max(bands, cell, n=21):
    """Max band frequency over a grid on the wedge 0 <= kz <= ky <= kx <= pi."""
    ks = np.linspace(0.0, np.pi, n)
    best, arg = 0.0, None
    for i, kx in enumerate(ks):
        for j, ky in enumerate(ks[:i + 1]):
            for kz in ks[:j + 1]:
                w = float(bands(np.array([kx, ky, kz]), cell).max())
                if w > best:
                    best, arg = w, (kx, ky, kz)
    return best, arg


def gate():
    checks, out = [], {}
    A = checks.append

    # ---- R1: the instrument reproduces jb_sj ---------------------------------
    d15, _ = RC.honeycomb(MJ.hc15_sites(), gc=A_REF)
    nz_d, f_d, Zq_d, c_d = spec(d15)
    rg_d, co_d, _, gap_d = identify(d15, Zq_d, c_d)
    out["R1"] = (d15.N, len(d15.welds), nz_d, len(f_d), f_d.min(), f_d.max())
    A(("R1  THE INSTRUMENT REPRODUCES jb_sj's PUBLISHED HC15 SPECTRUM BEFORE "
       "IT IS POINTED AT ANYTHING ELSE. Seven zero modes, identified as six "
       "rigid plus the coherent breathe and spanning; 98 real frequencies; "
       "the same extremes to five figures. [23746] lesson 4: validate the "
       "instrument against a published number first, because the tree-vs-"
       "cycle confound looked decisive until the instrument was checked. "
       "TWO-SIDED: any drift from jb_sj's numbers fails",
       d15.N == 15 and len(d15.welds) == 32 and nz_d == SJ_ZERO
       and len(f_d) == SJ_FREQS and abs(f_d.min() - SJ_MIN) < 1e-5
       and abs(f_d.max() - SJ_MAX) < 1e-5 and rg_d < 1e-12 and co_d < 1e-12
       and gap_d < 1e-10,
       f"HC15 double: {d15.N} cells, {len(d15.welds)} welds, {nz_d} zero, "
       f"{len(f_d)} frequencies {f_d.min():.5f}..{f_d.max():.5f}; rigid "
       f"{rg_d:.1e}, coherent {co_d:.1e}, span gap {gap_d:.1e}",
       f"7 zero, 98 frequencies, {SJ_MIN}..{SJ_MAX}, all identified"))

    # ---- R2: the same sites are a tree -----------------------------------------
    s15, _ = RC.honeycomb_single(MJ.hc15_sites(), gc=A_REF)
    nz_s, f_s, Zq_s, c_s = spec(s15)
    hv = tree_hinges(s15)
    rg_s, co_s, hg_s, gap_s = identify(s15, Zq_s, c_s, hv)
    out["R2"] = (s15.N, len(s15.welds), nz_s, len(f_s), f_s.min(), f_s.max())
    A(("R2  HC15's SITES UNDER ONE COVERING ARE A SEVEN-CELL STAR, A TREE, "
       "WITH THIRTEEN FREE MODES -- AND THAT IS THE PATCH TALKING, NOT THE "
       "MEDIUM. The centre keeps its six axis neighbours and loses its eight "
       "triangular ones to the voids; six bonds, no cycle. All thirteen zero "
       "modes are identified rather than counted: six rigid, the coherent "
       "breathe, and six HINGES, each a rigid rotation of one leaf about the "
       "line through its bond's two shared points, and the thirteen span the "
       "zero space. 'The 98 become 36' is the lesson-4 confound -- a tree "
       "against jb_sj's cyclic patch -- and this row exists so nobody quotes "
       "it. TWO-SIDED: a hinge outside the zero space, or a thirteenth mode "
       "the constructed set does not span, fails",
       s15.N == 7 and len(s15.welds) == 6 and nz_s == 13 and len(f_s) == 36
       and rg_s < 1e-12 and co_s < 1e-12 and hg_s < 1e-10 and gap_s < 1e-10,
       f"HC15 single: {s15.N} cells, {len(s15.welds)} welds, {nz_s} zero, "
       f"{len(f_s)} frequencies {f_s.min():.5f}..{f_s.max():.5f}; rigid "
       f"{rg_s:.1e}, coherent {co_s:.1e}, worst hinge {hg_s:.1e}, the "
       f"thirteen span to {gap_s:.1e}",
       "7 cells, 6 welds, 13 = 6 + 1 + 6 identified and spanning"))

    # ---- R3: the closed patches keep seven at every size ----------------------
    boxes = {}
    for side in SIDES:
        b = box(side)
        nz, f, Zq, ctr = spec(b)
        rg, co, _, gap = identify(b, Zq, ctr)
        boxes[side] = dict(N=b.N, welds=len(b.welds), nz=nz, f=f,
                           rigid=rg, coh=co, gap=gap)
    out["boxes"] = boxes
    ok3 = all(v["nz"] == 7 and len(v["f"]) == 7 * v["N"] - 7
              and v["rigid"] < 1e-12 and v["coh"] < 1e-12 and v["gap"] < 1e-10
              and v["f"].min() > 1e-6 for v in boxes.values())
    A(("R3  EVERY CLOSED SINGLE-COVERING PATCH HAS EXACTLY SEVEN FREE MODES, "
       "AND THEY ARE jb_sj's SEVEN. The all-even boxes from side 2 (jb_1c "
       "R4's cube) to side 6 (216 cells, 1512 DOF): at every size the zero "
       "space is six rigid motions plus the coherent breathe, each inside it "
       "to 1e-12 and the seven spanning it to 1e-10, and every other mode "
       "has a real nonzero frequency. The first cycle kills the hinges the "
       "star left free, and no further cycle frees anything. This is the "
       "half of jb_sj that IS about the medium, and it survives the change "
       "of covering at every size tried. TWO-SIDED in both directions: an "
       "eighth zero mode is a hinge that survived the cycles, a sixth is the "
       "breathe welded shut",
       ok3,
       "  ".join(f"side {s}: N={v['N']} welds={v['welds']} zero={v['nz']} "
                 f"freqs={len(v['f'])} gap={v['gap']:.0e}"
                 for s, v in sorted(boxes.items())),
       "seven identified and spanning at every side"))

    # ---- R4: the top sits under the bands and climbs to them ------------------
    bm_s, arg_s = bloch_max(OC.bands, OC.periodic_cell())
    bm_d, arg_d = bloch_max(BZ.bands, BZ.unit_cell())
    tops = [float(boxes[s]["f"].max()) for s in SIDES]
    gaps = [bm_s - t for t in tops]
    out["R4"] = (bm_s, arg_s, bm_d, arg_d, tops, gaps)
    A(("R4  NO FINITE PATCH RINGS ABOVE ITS OWN COVERING'S BANDS, AND THE "
       "SINGLE PATCHES CLIMB MONOTONICALLY TOWARD jb_1c's CEILING. A free "
       "patch's Hessian is the periodic operator with the boundary-crossing "
       "bonds removed and the mass matrix untouched, so its Rayleigh "
       "quotient can only be smaller: the top of every finite spectrum must "
       "lie under the Bloch maximum. jb_1c's bands top out at sqrt(3) at R; "
       "the boxes rise toward it with a gap that shrinks at every step. "
       "CONTROL ON THE OTHER COVERING: jb_sj's 2.08755 sits under jb_bz's "
       "Bloch maximum sqrt(6) -- the omega^2 = 6 Gamma mode the single "
       "covering does not have -- so the two coverings' finite spectra are "
       "each bounded by their own bands and the ceilings differ by sqrt(2). "
       "TWO-SIDED: a patch above its ceiling means the finite Hessian is not "
       "a restriction of the periodic one; a non-monotone climb means the "
       "series is not converging on it",
       abs(bm_s - np.sqrt(3.0)) < 1e-9 and abs(bm_d - np.sqrt(6.0)) < 1e-9
       and all(t <= bm_s + 1e-9 for t in tops)
       and all(tops[i + 1] > tops[i] for i in range(len(tops) - 1))
       and all(gaps[i + 1] < gaps[i] for i in range(len(gaps) - 1))
       and f_d.max() <= bm_d + 1e-9,
       f"single Bloch max {bm_s:.6f} at k={np.round(arg_s, 3)} (sqrt3 = "
       f"{np.sqrt(3.0):.6f}); box tops " + ", ".join(f"{t:.5f}" for t in tops)
       + "; gaps " + ", ".join(f"{g:.5f}" for g in gaps)
       + f"; double Bloch max {bm_d:.6f} at k={np.round(arg_d, 3)} (sqrt6 = "
       f"{np.sqrt(6.0):.6f}) against jb_sj's {f_d.max():.5f}",
       "every top under its ceiling, the single series rising with a "
       "shrinking gap"))

    # ---- R5: the bottom is acoustic ---------------------------------------------
    mins = [float(boxes[s]["f"].min()) for s in SIDES]
    prod = [m * s for m, s in zip(mins, SIDES)]
    fit_sides = [s for s in SIDES if s >= 3]
    fit_mins = [float(boxes[s]["f"].min()) for s in fit_sides]
    slope = float(np.polyfit(np.log(fit_sides), np.log(fit_mins), 1)[0])
    out["R5"] = (mins, prod, slope)
    A(("R5  THE LOWEST NONZERO FREQUENCY FALLS AS 1/SIDE: THERE IS NO GAP "
       "ABOVE THE SEVEN, AND THE BOTTOM OF THE SPECTRUM IS SOUND. A medium "
       "with two acoustic speeds must give a finite patch long-wavelength "
       "modes whose frequency scales inversely with its size. Measured: the "
       "minimum decreases at every step, min*side rises toward a bound, and "
       "the log-log slope over sides 3..6 is within 0.1 of -1. The cube is "
       "excluded from the fit as the one patch with no interior cell. "
       "TWO-SIDED: a floor that did not fall with size would be a gap -- an "
       "optical-only or surface-pinned bottom -- and a slope far from -1 "
       "would mean the lowest mode is not acoustic",
       all(mins[i + 1] < mins[i] for i in range(len(mins) - 1))
       and all(prod[i + 1] > prod[i] for i in range(len(prod) - 1))
       and abs(slope + 1.0) < 0.1,
       "min by side " + ", ".join(f"{m:.5f}" for m in mins)
       + "; min*side " + ", ".join(f"{p:.5f}" for p in prod)
       + f"; log-log slope over sides 3..6 = {slope:.4f}",
       "falling minimum, rising min*side, slope within 0.1 of -1"))

    # ---- R6: controls -------------------------------------------------------------
    cube = box(2)
    nz1, f1, _, _ = spec(cube, k=1.0)
    nz4, f4, _, _ = spec(cube, k=4.0)
    scale = float(np.max(np.abs(f4 / 2.0 - f1)))
    one, _ = RC.honeycomb_single([(0, 0, 0)], gc=A_REF)
    nz_1, f_1, _, _ = spec(one)
    out["R6"] = (scale, nz_1, len(one.welds))
    A(("R6  CONTROLS. Every frequency scales as sqrt(k) to machine precision, "
       "so k is a convention here as everywhere and only ratios and closed "
       "forms above are measurements. And the single-unit control of jb_sj "
       "R6 holds for the single builder too: one cell, no welds, every one "
       "of its seven modes free -- qvf.2's analytic death reproduced where "
       "it belongs, on one unit, by the instrument that finds seven on the "
       "array",
       nz1 == 7 and nz4 == 7 and scale < 1e-12 and nz_1 == 7
       and len(one.welds) == 0,
       f"k = 4 over 2 matches k = 1 to {scale:.1e}; single unit: {one.N} "
       f"cell, {len(one.welds)} welds, {nz_1} of 7 modes free",
       "exact sqrt(k) scaling, and every mode free on one unit"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_1s -- one covering: the spectrum on finite patches, against "
              "jb_sj and jb_1c")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        bm_s = out["R4"][0]
        print("\n  FINITE PATCHES AGAINST THE BANDS")
        print(f"   {'patch':20s} {'cells':>5s} {'bonds':>5s} {'DOF':>5s} "
              f"{'zero':>4s} {'freqs':>5s} {'min':>8s} {'min*side':>9s} "
              f"{'max':>8s} {'ceiling-max':>11s}")
        N, W, nz, nf, lo, hi = out["R1"]
        print(f"   {'HC15 double (jb_sj)':20s} {N:5d} {W:5d} {7 * N:5d} "
              f"{nz:4d} {nf:5d} {lo:8.5f} {'':>9s} {hi:8.5f} "
              f"{np.sqrt(6.0) - hi:11.5f}")
        N, W, nz, nf, lo, hi = out["R2"]
        print(f"   {'HC15 single (star)':20s} {N:5d} {W:5d} {7 * N:5d} "
              f"{nz:4d} {nf:5d} {lo:8.5f} {'':>9s} {hi:8.5f} "
              f"{bm_s - hi:11.5f}")
        for s in SIDES:
            v = out["boxes"][s]
            print(f"   {'box ' + str(s) + (' (cube)' if s == 2 else ''):20s} "
                  f"{v['N']:5d} {v['welds']:5d} {7 * v['N']:5d} {v['nz']:4d} "
                  f"{len(v['f']):5d} {v['f'].min():8.5f} "
                  f"{v['f'].min() * s:9.5f} {v['f'].max():8.5f} "
                  f"{bm_s - v['f'].max():11.5f}")
        print(f"   Bloch ceilings: single sqrt3 = {bm_s:.6f} at R; double "
              f"sqrt6 = {out['R4'][2]:.6f} at Gamma")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * jb_sj's SEVEN IS THE MEDIUM'S. Every closed single-covering")
        print("     patch has six rigid modes plus the coherent breathe and")
        print("     nothing else free. qvf.2's death stays cleared.")
        print("   * jb_sj's 98 ARE THE PATCH'S. They belong to a double-covered")
        print("     HC15, whose sites under one covering are a tree. Do not")
        print("     quote '98 become 36'; that compares a cycle to a star.")
        print("   * THE FINITE SPECTRUM IS BOUNDED BY jb_1c's BANDS and climbs")
        print("     to them; its bottom is acoustic. No closed form is claimed")
        print("     for min*side.")
        print("   * jb_sj R5 IS NOT RE-MEASURED: jb_mj's LCP has not been")
        print("     pointed at the single covering. The contact line is open.")
        print("   * RATIOS ONLY. k is a convention (R6).")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
