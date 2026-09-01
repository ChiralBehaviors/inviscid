"""inside_out_cover -- G4: X(a+180) = -X(a), exactly, on configurations.
Deck-map candidate (b) of the note, pinned. Near-formality by section 8.3;
measured here rather than left as an inference, and the "+ const" of the
conjecture comes out ZERO.

THE QUESTION (note section 3b, bead inviscid-dcz). The spacing law gives
L(a+180) = -L(a) and the far half of the overdrive cycle LOOKS like the
array inverted through its centre; the conjectured deck map for the
inside-out cover is central inversion composed with fold+180. The note
left one item open: "X(a+180) = -X(a) + const has not been verified on
configurations." This module verifies it on the record's three patches --
every vertex label, both body kinds -- and measures where the two sheets
of the would-be cover actually coincide in ambient space.

WHAT IS MEASURED HERE. The drive is the same closed form G1 rides
(assembly.coherent, unpacked by plate_holonomy.coherent_X), swept at 0.5
degrees. Algebra says the identity should be exact: every centre is
L(a) * site and L flips sign; every body corner is R(u, sigma(g-60))(v-c)
+ u Z cos g, the corner offsets lie in the plate plane perpendicular to u,
so advancing the fold 180 degrees rotates each plate half a turn about its
own axis -- which on its own plane IS the central inversion -- and flips
its height. The gate measures it instead of trusting the algebra.

  R1  THE IDENTITY, WITH CONST = 0: on ring, hc15 and block,
      X(a+180) = -X(a) at every vertex label on every 0.5-degree frame,
      and the fitted translation is zero -- even for the UNCENTRED
      patches (ring, block): the inversion centre is the lattice origin,
      not the patch centroid, because every term is odd.
  R2  THE RECORD'S TWO KNOWN FACTS ARE ITS COROLLARIES: the centre level
      reproduces L(a+180) = -L(a) exactly, and the front census
      complements per BODY -- fronts_out(a) + fronts_out(a+180) = 8 for
      every body at every angle where frontness is bounded away from
      zero (overdrive R3's 100% <-> 0% flip, body by body).
  R3  WHERE THE SHEETS TOUCH IN AMBIENT SPACE IS A PROPERTY OF THE
      PATCH, NOT THE DYNAMICS: one body's twelve vertices are a
      centrally symmetric SET at every fold angle, so a site-symmetric
      patch (hc15) has ambiently coincident sheets at EVERY angle --
      only identity tracking separates them -- while the uncentred ring
      and block coincide ONLY at the collapses, where L = 0 wipes the
      centres. Section 3b's "the collapses are where the sheets touch"
      is the uncentred-patch statement.

WHAT THIS DOES NOT SAY. Pinning the map is not constructing the cover:
nothing here builds a quotient of shape space, computes a fundamental
group, or licenses the period-180 language of section 3b as more than a
description of this identity (section 8.1 discipline). The ambient
indistinguishability of a symmetric patch's sheets feeds section 3c's
two-sided-plate reading but adds no degree-of-freedom claim.

T2: [23865]. Ref: su2_boundary_conditions.md section 8 (read before 1-7).
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model.first_principles import geometry as G
from analysis.model.first_principles.overdrive import PATCHES
from analysis.model.su2.plate_holonomy import coherent_X

#: (a, a+180) census pairs where no plate's frontness crosses zero
#: (cells cross at 90/270, voids at 30/210; both members must be clean)
CENSUS_A = [-60.0, -30.0, 0.0, 60.0, 120.0]


def centres(sites, a):
    """The patch's body centres L(a) * site, from the coherent drive."""
    q, _ = RC.coherent(sites, a)
    return RC.Assembly.unpack(q)[0]


def set_mismatch(P, Q):
    """max over p in P of the distance to its nearest q in Q."""
    D = np.linalg.norm(P[:, None, :] - Q[None, :, :], axis=2)
    return float(D.min(axis=1).max())


def gate():
    checks = []
    A = checks.append
    degs = np.arange(-60.0, 120.0 + 1e-9, 0.5)

    # ---- R1: the identity, with const = 0 -------------------------------------
    worst_res, worst_c = {}, {}
    X = {n: {} for n in PATCHES}
    for n, sites in PATCHES.items():
        res = c = 0.0
        for a in degs:
            Xa, Xb = coherent_X(sites, a), coherent_X(sites, a + 180.0)
            X[n][a], X[n][a + 180.0] = Xa, Xb
            S = Xb + Xa                                   # = const if the map holds
            cv = S.reshape(-1, 3).mean(0)
            res = max(res, float(np.abs(S - cv).max()))
            c = max(c, float(np.abs(cv).max()))
        worst_res[n], worst_c[n] = res, c
    A(("R1  THE IDENTITY, WITH CONST = 0: X(a+180) = -X(a) at every vertex label on "
       "every 0.5-degree frame of all three patches, and the fitted translation is "
       "zero even for the UNCENTRED ring and block -- the inversion centre is the "
       "lattice origin, not the patch centroid, because every term is odd in the map. "
       "Deck-map candidate (b) -- central inversion o fold+180 -- pinned on "
       "configurations, label- and kind-preserving, no translation part.",
       all(v < 1e-9 for v in worst_res.values()) and all(v < 1e-9 for v in worst_c.values()),
       {n: f"res {worst_res[n]:.0e}, const {worst_c[n]:.0e}" for n in PATCHES},
       "all < 1e-9"))

    # ---- R2: the record's two known facts as corollaries ----------------------
    lerr = max(abs(RC.lattice_constant(a + 180.0) + RC.lattice_constant(a)) for a in degs)
    worst_census = {}
    for n, sites in PATCHES.items():
        S0 = G.front_signs(coherent_X(sites, 0.0), centres(sites, 0.0))
        bad = 0
        for a in CENSUS_A:
            oa = G.fronts_out(X[n][a], centres(sites, a), S0)
            ob = G.fronts_out(X[n][a + 180.0], centres(sites, a + 180.0), S0)
            bad = max(bad, int(np.abs(oa + ob - 8).max()))
        worst_census[n] = bad
    A(("R2  THE RECORD'S TWO KNOWN FACTS ARE ITS COROLLARIES: the centre level "
       "reproduces L(a+180) = -L(a) exactly, and the front census complements per "
       "BODY -- every body has fronts_out(a) + fronts_out(a+180) = 8 at every clean "
       "census angle: overdrive R3's 100% <-> 0% flip, body by body.",
       lerr < 1e-12 and all(v == 0 for v in worst_census.values()),
       f"L identity {lerr:.0e}; census residue " + str(worst_census),
       "< 1e-12; all 0"))

    # ---- R3: where the sheets touch is a property of the patch ----------------
    body_cs = max(set_mismatch(RC.body(a), -RC.body(a)) for a in degs)
    generic, collapse = (-30.0, 0.0, 17.0, 105.0), 60.0
    touch = {}
    for n in PATCHES:
        gm = min(set_mismatch(X[n][a + 180.0].reshape(-1, 3), X[n][a].reshape(-1, 3))
                 for a in generic)
        cm = set_mismatch(X[n][collapse + 180.0].reshape(-1, 3),
                          X[n][collapse].reshape(-1, 3))
        touch[n] = (gm, cm)
    ok3 = (body_cs < 1e-12 and touch["hc15"][0] < 1e-12 and touch["hc15"][1] < 1e-12
           and touch["ring"][0] > 1e-2 and touch["block"][0] > 1e-2
           and touch["ring"][1] < 1e-12 and touch["block"][1] < 1e-12)
    A(("R3  WHERE THE SHEETS TOUCH IN AMBIENT SPACE IS A PROPERTY OF THE PATCH, NOT "
       "THE DYNAMICS: one body's twelve vertices are a centrally symmetric SET at "
       "every fold angle, so the site-symmetric hc15's two sheets are ambiently "
       "coincident at EVERY angle -- only identity tracking separates them -- while "
       "the uncentred ring and block coincide ONLY at the collapses, where L = 0 "
       "wipes the centres. Section 3b's picture is the uncentred-patch statement.",
       ok3,
       f"body set vs -set {body_cs:.0e}; generic/collapse mismatch "
       + str({n: f"{g:.0e} / {cme:.0e}" for n, (g, cme) in touch.items()}),
       "body < 1e-12; hc15 both < 1e-12; ring, block > 1e-2 generic, < 1e-12 collapse"))
    return checks


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("inside_out_cover -- G4: X(a+180) = -X(a) exactly; deck-map candidate (b) "
          "pinned")
    print("=" * 78)
    checks = gate()
    bad = 0
    for name, ok, got, want in checks:
        bad += 0 if ok else 1
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")
    print()
    print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
    print("   * CANDIDATE (b) IS PINNED: central inversion o fold+180 is an exact,")
    print("     label- and kind-preserving identity of the coherent drive, with NO")
    print("     translation part, on centred and uncentred patches alike.")
    print("   * IT DOES NOT CONSTRUCT THE COVER: no quotient of shape space is")
    print("     built, no fundamental group computed; period-180 language remains")
    print("     a description of this identity (section 8.1 discipline).")
    print("   * SHEET-TOUCHING IS PATCH GEOMETRY: for site-symmetric patches the")
    print("     halves of the cycle are ambiently indistinguishable everywhere --")
    print("     grist for section 3c's two-sided-plate reading, not a new claim.")
    print()
    print("  ALL CHECKS PASSED." if not bad
          else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
