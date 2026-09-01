"""screw_prerequisites -- the two checks the 60-degree-screw quotient owes
before it is even a candidate closure (note sections 8.1, 8.3; bead
inviscid-4g6): does tau^6 reduce to a pure lattice translation on
configurations, and is the action free and properly discontinuous? Both
answered by measurement -- with the first answer sharper than the note's
question: single applications of tau LEAVE the coherent family, and the
note's Bieberbach language stays withdrawn.

THE QUESTION (note section 3a). The ring gated that the void runs sixty
degrees ahead of the cell (b = a + 60), and section 3a inferred a screw
symmetry tau = (translate by (1,1,1)) o (fold +60) o (cell<->void),
proposing a quotient of space by it -- with a flat-3-manifold
classification attached. Section 8.1 objected that fold+60 moves between
fibers of the family (the correct home is a mapping torus), and 8.3 that
nothing gated checks tau^6 = translation, freeness, proper
discontinuity, or injective tiling. This module measures what tau
actually is.

WHAT TAU IS HERE. The only version that is well-defined per body: carry
each body to the site (1,1,1) away (position += L(a) * (1,1,1)) and
advance ITS OWN fold by sixty. gamma_s = a + 60 p(s) with p the site
parity class, and the site shift flips p, so a single application can
match the family's fold book-keeping for at most one parity class at a
time -- the obstruction is structural, not numerical:
60 p + 60 k = 60 (p XOR (k mod 2)) mod 360 has solutions only at
(k=1, p=0), (k=5, p=1), and k = 0 mod 6 for both.

  R1  SINGLE APPLICATIONS LEAVE THE COHERENT FAMILY: body-for-body over
      the physical arc on all three patches, tau^1 carries exactly the
      ex-cells onto void states (ring.py's one-sided law), tau^5 exactly
      the ex-voids onto cell states, tau^2, tau^3, tau^4 carry NOTHING,
      and tau^6 carries everything. The named single-step tau is not a
      symmetry of the family in any composition; only its sixth power
      returns.
  R2  TAU^6 IS A PURE LATTICE TRANSLATION, EXACTLY: fold+360 is the
      identity on bodies, so tau^6(config(a)) = config(a) + 6 L(a)
      (1,1,1), body for body. And tau^3 is (translate by 3 L (1,1,1))
      composed with PER-BODY central inversion -- the G4 law
      body(g+180) = -body(g) read per body: the screw's half-period is
      the inside-out cover's deck map acting body-locally.
  R3  FREE AND PROPERLY DISCONTINUOUS ON THE PHYSICAL FAMILY, DEGENERATE
      AT THE COLLAPSES: the orbit of any body centre is spaced
      sqrt(3) L(a) >= sqrt(3) on the whole arc (L >= 1, minimum exactly
      1 at both octahedron ends), so the Z-action is free with discrete
      orbits -- injective tiling included; at a = 60 and 240 L = 0
      EXACTLY and the entire orbit collapses onto one point. The screw
      quotient is a candidate on the physical family only.

WHAT THIS LICENSES AND WHAT IT DOES NOT. G6 (inviscid-y13) may proceed,
but its identification must be the MAPPING-TORUS form -- fields at
(s + (1,1,1), a) tied to fields at (s, a + 60) -- not a fixed-fiber
spatial screw, which R1 rules out. The Bieberbach/flat-manifold naming
of section 3a stays withdrawn: tau^6 = translation and proper
discontinuity are necessary for it, not sufficient, and no manifold is
constructed here.

T2: [23865], [23914]. Ref: su2_boundary_conditions.md section 8.
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model.first_principles.overdrive import PATCHES

ARC = np.arange(-60.0, 0.0 + 1e-9, 7.5)          # fold samples for body checks
FINE = np.arange(-60.0, 0.0 + 1e-9, 0.5)         # spacing-law samples


def parity(s):
    return 0 if all(c % 2 == 0 for c in s) else 1


def match_table(sites, a):
    """Per iterate k = 1..6: how many bodies land on a family body state.

    A body of parity p carries fold a + 60 p + 60 k to a site of parity
    p XOR (k mod 2), where the family holds fold a + 60 (p XOR (k mod 2));
    the comparison is on the full body shape, not the fold arithmetic.
    """
    counts = {}
    for k in range(1, 7):
        n = 0
        for s in sites:
            p = parity(s)
            got = RC.body(a + 60.0 * p + 60.0 * k)
            want = RC.body(a + 60.0 * (p ^ (k % 2)))
            n += int(float(np.abs(got - want).max()) < 1e-12)
        counts[k] = n
    return counts


def gate():
    checks = []
    A = checks.append

    # ---- R1: single applications leave the family -----------------------------
    ok1, got1 = True, {}
    for n, sites in PATCHES.items():
        nc = sum(1 for s in sites if parity(s) == 0)
        nv = len(sites) - nc
        expect = {1: nc, 2: 0, 3: 0, 4: 0, 5: nv, 6: len(sites)}
        worst = None
        for a in ARC:
            t = match_table(sites, a)
            if t != expect:
                worst = (a, t)
        ok1 &= worst is None
        got1[n] = f"k=1..6 {list(match_table(sites, ARC[0]).values())} at every a" \
            if worst is None else f"MISMATCH at a={worst[0]}: {worst[1]}"
    A(("R1  SINGLE APPLICATIONS LEAVE THE COHERENT FAMILY: over the physical arc on "
       "all three patches, tau^1 carries exactly the ex-cells onto void states "
       "(ring.py's one-sided law), tau^5 exactly the ex-voids onto cell states, "
       "tau^2..4 carry NOTHING, tau^6 carries everything. The parity flip demands "
       "+60 for one class and -60 for the other at once; only the sixth power "
       "returns to the family.",
       ok1, got1 | {"expect": "cells, 0, 0, 0, voids, all"},
       "ring 2,0,0,0,2,4; hc15 7,0,0,0,8,15; block 27,0,0,0,8,35"))

    # ---- R2: tau^6 is a pure lattice translation, exactly ---------------------
    per360 = max(float(np.abs(RC.body(g + 360.0) - RC.body(g)).max())
                 for a in ARC for g in (a, a + 60.0))
    inv180 = max(float(np.abs(RC.body(g + 180.0) + RC.body(g)).max())
                 for a in ARC for g in (a, a + 60.0))
    tr6 = 0.0
    for n, sites in PATCHES.items():
        for a in (-30.0, -7.5):
            L = RC.lattice_constant(a)
            for s in sites:
                p = parity(s)
                img = L * (np.array(s, float) + 6.0) + RC.body(a + 60.0 * p + 360.0)
                fam = L * np.array(s, float) + RC.body(a + 60.0 * p)
                tr6 = max(tr6, float(np.abs(img - (fam + L * 6.0)).max()))
    A(("R2  TAU^6 IS A PURE LATTICE TRANSLATION, EXACTLY: fold+360 is the identity "
       "on bodies, so tau^6(config(a)) = config(a) + 6 L(a) (1,1,1) body for body; "
       "and tau^3 is translate(3 L (1,1,1)) composed with PER-BODY central inversion "
       "(body(g+180) = -body(g), the G4 law body-locally) -- the screw's half-period "
       "is the inside-out deck map.",
       per360 < 1e-12 and inv180 < 1e-12 and tr6 < 1e-12,
       f"|body(g+360)-body(g)| {per360:.0e}; |body(g+180)+body(g)| {inv180:.0e}; "
       f"tau^6 vs translation {tr6:.0e}",
       "all < 1e-12"))

    # ---- R3: free + properly discontinuous on the arc, degenerate at collapses -
    Ls = np.array([RC.lattice_constant(a) for a in FINE])
    ends = (abs(RC.lattice_constant(-60.0) - 1.0), abs(RC.lattice_constant(0.0) - 1.0))
    coll = (abs(RC.lattice_constant(60.0)), abs(RC.lattice_constant(240.0)))
    spacing = float(np.sqrt(3.0) * Ls.min())
    A(("R3  FREE AND PROPERLY DISCONTINUOUS ON THE PHYSICAL FAMILY, DEGENERATE AT "
       "THE COLLAPSES: the tau-orbit of any body centre is spaced sqrt(3) L(a), and "
       "L >= 1 on the whole arc with the minimum exactly 1 at both octahedron ends "
       "-- no fixed bodies, discrete orbits, injective tiling. At a = 60 and 240 "
       "L = 0 EXACTLY: the whole orbit lands on one point. The screw quotient is a "
       "candidate on the physical family only, never through the collapses.",
       Ls.min() >= 1.0 - 1e-12 and max(ends) < 1e-12 and max(coll) < 1e-12
       and spacing > 1.7,
       f"min L on arc {Ls.min():.12f} (ends off by {max(ends):.0e}); orbit spacing "
       f">= {spacing:.6f}; L at collapses {coll[0]:.0e}, {coll[1]:.0e}",
       "min L = 1 at the ends; spacing sqrt(3); collapses 0"))
    return checks


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("screw_prerequisites -- tau^6 = pure translation; free and properly "
          "discontinuous on the arc")
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
    print("   * G6 MAY PROCEED, IN MAPPING-TORUS FORM: identify (s+(1,1,1), a) with")
    print("     (s, a+60) -- a fixed-fiber spatial screw is ruled out by R1; single")
    print("     applications of tau are not symmetries of any one configuration.")
    print("   * THE QUOTIENT'S DOMAIN IS THE PHYSICAL FAMILY: proper discontinuity")
    print("     holds with margin on the arc and fails exactly at the collapses.")
    print("   * THE BIEBERBACH NAMING STAYS WITHDRAWN (8.1/8.3): these checks are")
    print("     necessary for it, not sufficient; no flat manifold is constructed.")
    print()
    print("  ALL CHECKS PASSED." if not bad
          else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
