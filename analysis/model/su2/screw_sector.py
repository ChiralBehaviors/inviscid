"""screw_sector -- G6, resolved in the NEGATIVE and gated: there is no
60-degree-screw Bloch sector on the medium of record, and the original
boundary-condition question is answered by G5's race, to which the screw
adds nothing.

THE QUESTION (note section 5-G6; bead inviscid-y13). Dispersion with the
screw identification along (1,1,1) versus plain periodic. The bead asked
for a screw-quotient assembly constructor. screw_prerequisites (this
lane) had already shown the fixed-fiber spatial screw is not a symmetry
of any configuration; this module measures the rest of the ground the
sector would need, and finds none of it on the medium of record.

WHY THE NEGATIVE IS THE RESULT (owner decision, 2026-09-01: option (a)
of the reformulation recorded on inviscid-y13). Three independent gated
facts close the lane:
  1. single applications of tau leave the coherent family
     (screw_prerequisites R1);
  2. the (1,1,1) translation carries every SOLID site of the single
     covering to an EMPTY void site (measured here) -- under DECISION 21
     the screw has nothing to act on;
  3. the mapping-torus identification that survives constrains only
     full-cycle driven protocols, which the crossing walls exclude
     physically (physical_arc R2).

  R1  THE SCREW LANE IS EMPTY ON THE MEDIUM OF RECORD: odd multiples of
      (1,1,1) put ZERO solid sites on solid sites (k = 1, 3, 5 measured
      on a 4^3 solid box); even multiples recover the full interior
      overlap. The screw's only trace on the solid sublattice is its
      even powers -- plain lattice translations, i.e. the periodic torus
      G5 already raced.
  R2  THE FIXED-FOLD TWISTED SECTOR FAILS AT EVERY GENERIC FOLD AND
      BELONGS TO THE RETIRED COVERING AT THE ONE FOLD WHERE IT EXISTS:
      searching all 48 signed permutations, NO isometry maps the cell
      body onto the void body at a generic fold (0/48 at a = -17), so
      no point operation can repair T(1,1,1) into a symmetry of the
      fixed-fold operator; at EXACTLY the symmetric phase a = -30,
      24 of 48 succeed (12 proper, 12 improper) -- a glide-type pairing
      real at that one fold -- but its home is the both-kinds covering,
      whose odd sites the record rules EMPTY. On the medium there is
      nothing at the target sites at any fold.
  R3  THE BC QUESTION IS ANSWERED WHERE THE MACHINERY WAS BUILT: with
      the screw lane empty, the closure candidates are free / torus /
      double, and the G5 race is re-measured here as the cross-check --
      at side 3 the double and torus sit within 3 percent of each other
      and both under half the free block's distance to the Bloch DOS.
      The remaining choice, double versus torus, is a modeling decision
      G5 informs; no screw sector exists to widen it.

WHAT THIS DOES NOT SAY. The glide pairing at a = -30 is a real fact
about the both-kinds INSTRUMENT covering and is recorded as such (24/48,
counted); building its twisted 14-band dispersion would measure the
retired covering, not the medium, and is deliberately not done. Nothing
here touches the kinematic Z2 (G1): the lane's verdict stands --
correct kinematics, physically inert, no carrier, and now: no spectral
sector either.

T2: [23865], [23914]. Ref: su2_boundary_conditions.md section 8.
"""
from __future__ import annotations

import itertools as it
import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model import dispersion as OC
from analysis.model.double_covering import soft_joint_spectrum as SJ
from analysis.model.su2.doubled_block import double, kolmogorov, torus_eigs


def setmis(P, Q):
    D = np.linalg.norm(P[:, None, :] - Q[None, :, :], axis=2)
    return float(D.min(axis=1).max())


def oh_elements():
    for perm in it.permutations(range(3)):
        for signs in it.product((1, -1), repeat=3):
            M = np.zeros((3, 3))
            for i, p in enumerate(perm):
                M[i, p] = signs[i]
            yield M


def pairing_census(a):
    """How many of the 48 signed permutations map body(a) onto body(a+60),
    split by determinant."""
    Ba, Bb = RC.body(a), RC.body(a + 60.0)
    proper = improper = 0
    for M in oh_elements():
        if setmis((M @ Ba.T).T, Bb) < 1e-9:
            if round(float(np.linalg.det(M))) == 1:
                proper += 1
            else:
                improper += 1
    return proper, improper


def gate():
    checks = []
    A = checks.append

    # ---- R1: the screw lane is empty on the medium of record ------------------
    solid = {(2 * x, 2 * y, 2 * z) for x in range(4) for y in range(4) for z in range(4)}
    overlap = {}
    for k in range(1, 7):
        t = {(s[0] + k, s[1] + k, s[2] + k) for s in solid}
        interior = {s for s in solid if all(c + k < 8 for c in s)}
        overlap[k] = (len(t & solid), len(interior) if k % 2 == 0 else 0)
    ok1 = all(got == exp for got, exp in overlap.values())
    A(("R1  THE SCREW LANE IS EMPTY ON THE MEDIUM OF RECORD: odd multiples of "
       "(1,1,1) put ZERO solid sites of the single covering on solid sites; even "
       "multiples recover the full interior overlap. Under DECISION 21 (voids "
       "EMPTY) the screw identification has nothing to act on, and its only trace "
       "on the solid sublattice is its even powers -- plain translations, the "
       "periodic torus G5 already raced.",
       ok1,
       {k: f"{got}/{exp}" for k, (got, exp) in overlap.items()},
       "0 at k=1,3,5; full interior at k=2,4,6"))

    # ---- R2: no fixed-fold twisted sector; the glide belongs to the retired ----
    generic = {a: pairing_census(a) for a in (-17.0, -45.0, 7.0)}
    sym = pairing_census(-30.0)
    ok2 = (all(p == 0 and i == 0 for p, i in generic.values())
           and sym == (12, 12))
    A(("R2  THE FIXED-FOLD TWISTED SECTOR FAILS AT EVERY GENERIC FOLD AND BELONGS "
       "TO THE RETIRED COVERING AT THE ONE FOLD WHERE IT EXISTS: no signed "
       "permutation maps the cell body onto the void body at generic folds, so no "
       "point operation repairs T(1,1,1) into a symmetry of the fixed-fold "
       "operator; at exactly a = -30 the pairing exists (12 proper + 12 improper "
       "of 48) -- a glide-type symmetry of the BOTH-KINDS covering, whose odd "
       "sites the record rules empty. On the medium: nothing at the target sites, "
       "at any fold.",
       ok2,
       {"generic": {a: f"{p}+{i}" for a, (p, i) in generic.items()},
        "a=-30": f"{sym[0]} proper + {sym[1]} improper"},
       "0+0 generic; 12+12 at -30"))

    # ---- R3: the BC question is answered where the machinery was built --------
    cell = OC.periodic_cell()
    n = 3
    NREF = 20
    ks = 2.0 * np.pi * np.arange(NREF) / NREF
    ref = np.concatenate([OC.bands(np.array([kx, ky, kz]), cell)
                          for kx in ks for ky in ks for kz in ks])
    dbl, free = double(n)
    Df = kolmogorov(np.sqrt(np.clip(SJ.spectrum(free)[0], 0, None)), ref)
    Dd = kolmogorov(np.sqrt(np.clip(SJ.spectrum(dbl)[0], 0, None)), ref)
    Dt = kolmogorov(np.sqrt(np.clip(torus_eigs(n, cell)[0], 0, None)), ref)
    ok3 = (Dd < Df / 2 and Dt < Df / 2 and abs(Dd - Dt) < 0.03)
    A(("R3  THE BC QUESTION IS ANSWERED WHERE THE MACHINERY WAS BUILT: with the "
       "screw lane empty, the closure candidates are free / torus / double, and "
       "the G5 race cross-checks here -- at side 3 the double and torus sit within "
       "three points of each other and both under half the free block's distance "
       "to the Bloch density of states. The remaining choice, double versus "
       "torus, is a modeling decision G5 informs; no screw sector exists to "
       "widen it.",
       ok3,
       f"side 3: free {Df:.4f} / double {Dd:.4f} / torus {Dt:.4f}",
       "double, torus < free/2 and within 0.03 of each other"))
    return checks


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("screw_sector -- G6 resolved in the negative: no screw Bloch sector on "
          "the medium of record")
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
    print("   * G6 IS CLOSED AS A MEASURED NEGATIVE: no screw-quotient constructor")
    print("     can exist for the single covering (nothing at the target sites),")
    print("     no fixed-fold twisted sector exists at any generic fold, and the")
    print("     glide pairing at a = -30 belongs to the retired both-kinds")
    print("     covering (recorded, not pursued).")
    print("   * THE ORIGINAL BC QUESTION LANDS ON G5: free is dominated; double")
    print("     and torus are spectrally equivalent closures at the sizes tried;")
    print("     choosing between them is a modeling decision, now fully informed.")
    print("   * NOTHING HERE RETRACTS THE KINEMATICS: the Z2 of G1 stands, as a")
    print("     property of the analytic continuation only (G2, G7).")
    print()
    print("  ALL CHECKS PASSED." if not bad
          else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
