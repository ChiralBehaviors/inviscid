"""held_correspondences -- read once and held: which holds are labellings, and which are luck.

WHY THIS FILE EXISTS. Bead inviscid-43o filed the THIRD instance of one defect
class and asked for the other held correspondences to be audited against a
distinction nobody had stated precisely:

    is the thing being held a LABELLING, which must not be re-read, or a
    COINCIDENCE, which only holds where it was read?

The idiom itself is not wrong. The hinge pairing genuinely must not be re-read,
and R2 measures exactly why. What makes the idiom dangerous is that both cases
look identical in the code -- a tuple computed once at module load -- and the
one place the distinction was written down, `_fcc13_contacts`' own docstring,
gets it BACKWARDS:

    "Read at a = 0 and never re-read, for the same reason the hinge pairing is
     never re-read at a merge angle."

Measured here, those are OPPOSITE cases. The hinge pairing is invariant at every
generic phase and degenerates only AT the merge angles, so holding it is
correct. FCC13's neighbour-neighbour contacts exist at a = 0 and NOWHERE ELSE,
so holding them is the bug. The sentence that justified the defect is the one
that named the right precedent for the wrong reason.

THE CRITERION IS CHEAP AND MECHANICAL. Re-read the correspondence at a spread of
phases and compare. That is all it takes, and it is the whole content of this
file -- the audit had never been run because nobody had written down what to
run.

WHAT THE AUDIT FOUND, and the good news is the smaller number:

    hinge pairing              LABELLING    invariant at 9 generic phases
    N1, N2, CHAIN5, SQUARE4    labelling    generated from lattice generators
    SC7, CUBE8-M/R, CUBE27-M   labelling    same, invariant to 4e-15
    FCC13 twelve-around-one    COINCIDENCE  only at a = 0

One instance in jb_x, and it is the one already filed. The eight box- and
star-kind topologies build their contacts from lattice GENERATORS rather than
reading them off geometry, which is why they cannot have this defect: there is
no coincidence to mistake for a structure.

AND THE FCC13 DEFECT SPLITS, which the bead did not know. Of its 36 contacts:

  * the 12 CENTRE-to-neighbour contacts are EXACT at every phase, to 4e-16,
    and not by luck -- unit 0's vertex k and unit 1+k's vertex ANTI[k] both
    land at 2 v[k] identically, because v[ANTI[k]] = -v[k] at every fold
    angle. That half is a labelling.
  * all 24 NEIGHBOUR-to-neighbour contacts fail, and they carry the entire
    residual the bead reports.

REGENERATING PER PHASE DOES NOT RESCUE IT (R4). Re-read at any a != 0 and the
neighbour-neighbour search returns NOTHING: 12 contacts, zero of them
neighbour-to-neighbour. The twelve neighbours stop touching each other the
instant the jitterbug moves. So "twelve-around-one closes ONLY at a = 0"
(jb_x's row X2) is a GEOMETRIC FACT and not an artifact of the held list --
but it is a fact about the neighbour-neighbour half specifically, and X2's
wording attributes to the whole structure what is true of that half.

CONSEQUENCE FOR THE THREE OPTIONS bead 43o offers. Regenerating per phase is
well posed but yields a DIFFERENT OBJECT -- a 13-cell star, which SC7 already
covers -- so it would quietly replace the topology rather than fix it.
Deleting discards a legitimate a = 0 structure. The honest option is the third:
FCC13 is an a = 0 topology, said out loud, with the reason attached.

SCOPE. This audits the held sets reachable from jb_x. The honeycomb's contact
list (read at a = -30, deliberately, because the squares fold) is jb_ht's own
T3 row and is not re-measured here; the weld corner maps are jb_ic's R2b and
jb_cl's C2. Those three are named so the reader knows this audit is not
claiming to have covered them.
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.retired.rig_lock import array_linkage as X
#: Generic phases. Deliberately spread, and deliberately including negatives --
#: a correspondence that is invariant on [0, 60] only would still be a
#: coincidence for the medium, whose range is [-60, 0].
GENERIC = (1.0, 10.0, 22.238756093, 30.0, 45.0, 55.0, 59.0, -30.0, -55.0)

#: Where the twelve vertices merge in pairs into six of multiplicity 4. The
#: hinge pairing must NOT be re-read here, which is the whole reason the
#: read-once idiom exists.
MERGE = (60.0, 120.0)


def pairing_at(a):
    """`_read_pairing`'s own construction, at an arbitrary phase.

    Reproduced rather than called because `_read_pairing` hardcodes its probe
    angle -- which is the point: the invariance it asserts had never been
    exercised at any other angle from outside it.
    """
    reps, mult, labels = X.cluster(X.corners(a), 1e-7)
    if len(reps) != 12 or sorted(set(mult.tolist())) != [2]:
        return None
    lab = labels.reshape(8, 3)
    slots = [[] for _ in range(12)]
    for i in range(8):
        for j in range(3):
            slots[lab[i, j]].append((i, j))
    return tuple(tuple(s) for s in slots)


def residual(topo, a, contacts=None):
    """Worst contact residual under the reference (pure-translate) placement."""
    v = X.verts(a)
    o = topo.sites(v)
    cs = topo.contacts if contacts is None else contacts
    if not cs:
        return 0.0
    return max(float(np.linalg.norm((v[k] + o[i]) - (v[l] + o[j])))
               for (i, k, j, l) in cs)


def fcc13():
    return [t for t in X.build_topologies() if t.name.startswith("FCC13")][0]


def gate():
    checks, out = [], {}
    A = checks.append

    # R1 -- the criterion, and it must separate the two known cases
    base = pairing_at(30.0)
    pair_same = all(pairing_at(a) == base for a in GENERIC)
    f = fcc13()
    fcc_worst = max(residual(f, a) for a in GENERIC)
    A(("R1  THE CRITERION, AND IT SEPARATES THE TWO CASES THAT LOOK IDENTICAL "
       "IN THE CODE. Re-read the correspondence at a spread of phases and "
       "compare: a LABELLING is invariant, a COINCIDENCE is not. That is the "
       "whole audit bead 43o asked for, and it had never been run because "
       "nobody had written down what to run -- both cases are a tuple computed "
       "once at module load and are indistinguishable by inspection. "
       "TWO-SIDED, and this row is worthless without both halves: the hinge "
       "pairing must come back INVARIANT and FCC13's contact set must come "
       "back NOT",
       pair_same and fcc_worst > 1e-3,
       f"hinge pairing identical at all {len(GENERIC)} generic phases: "
       f"{pair_same}; FCC13 worst contact residual over the same phases: "
       f"{fcc_worst:.6f}",
       "the criterion says labelling for one and coincidence for the other"))

    # R2 -- the hinge pairing is a labelling, and the merge angles are why
    merged = {a: X.cluster(X.corners(a), 1e-7) for a in MERGE}
    merge_ok = all(len(m[0]) == 6 and sorted(set(m[1].tolist())) == [4]
                   for m in merged.values())
    A(("R2  THE HINGE PAIRING IS A LABELLING, SO THE IDIOM IS CORRECT WHERE IT "
       "WAS INVENTED. It is identical at nine generic phases including negative "
       "ones, and at the merge angles 60 and 120 the twelve vertices collapse "
       "to SIX of multiplicity 4 -- which is exactly the configuration where "
       "re-reading would rewire the linkage, and exactly what the read-once "
       "rule exists to prevent. `_read_pairing` asserts this invariance in its "
       "own docstring and hardcodes its probe angle, so it had never been "
       "exercised at another angle from outside it; this row does that",
       pair_same and merge_ok,
       f"pairing identical at {len(GENERIC)} generic phases; at a = 60 and 120: "
       + "; ".join(f"{len(m[0])} distinct vertices, multiplicities "
                   f"{sorted(set(m[1].tolist()))}" for m in merged.values()),
       "invariant where it is read, degenerate exactly where it must not be"))

    # R3 -- FCC13's defect splits, and only one half is broken
    cn = [c for c in f.contacts if c[0] == 0]
    nn = [c for c in f.contacts if c[0] != 0]
    probe = (0.0, 1.0, 5.0, 22.238756093, 40.0)
    cn_worst = max(residual(f, a, cn) for a in probe)
    nn_worst = max(residual(f, a, nn) for a in probe)
    out["split"] = (cn_worst, nn_worst)
    A(("R3  THE DEFECT SPLITS, AND HALF OF FCC13 IS EXACT AT EVERY PHASE. Of "
       "its 36 contacts the 12 CENTRE-to-neighbour ones hold to 4e-16 "
       "everywhere, and not by luck: unit 0's vertex k and unit 1+k's vertex "
       "ANTI[k] both land at 2 v[k] identically, because v[ANTI[k]] = -v[k] at "
       "every fold angle. The 24 NEIGHBOUR-to-neighbour ones carry the ENTIRE "
       "residual the bead reports -- 0.069809626 at a = 1, 0.348622971 at "
       "a = 5, reproduced here. So 'the contact set is infeasible' is true of "
       "two thirds of it, and the star at its centre was never the problem",
       cn_worst < 1e-12 and nn_worst > 1.0
       and abs(residual(f, 1.0, nn) - 0.069809626) < 1e-8
       and abs(residual(f, 5.0, nn) - 0.348622971) < 1e-8,
       f"centre-to-neighbour worst {cn_worst:.3e} over {len(probe)} phases; "
       f"neighbour-to-neighbour worst {nn_worst:.6f}; the bead's numbers "
       f"reproduce: a=1 -> {residual(f, 1.0, nn):.9f}, "
       f"a=5 -> {residual(f, 5.0, nn):.9f}",
       "centre half exact everywhere, neighbour half carrying all of it"))

    # R4 -- and regenerating does not rescue it
    regen = {}
    for a in (0.0, 1.0, 5.0, 22.238756093, 40.0):
        cs = X._fcc13_contacts(X.verts(a))
        regen[a] = (len(cs), len([c for c in cs if c[0] != 0]),
                    residual(f, a, cs))
    out["regen"] = regen
    A(("R4  REGENERATING PER PHASE DOES NOT RESCUE IT -- THE NEIGHBOURS STOP "
       "TOUCHING EACH OTHER. Re-read the set at any a != 0 and the "
       "neighbour-neighbour search returns NOTHING: 12 contacts, none of them "
       "between neighbours, against 36 with 24 at a = 0. Twelve-around-one's "
       "neighbour contact is a property of the a = 0 cuboctahedron -- whose "
       "vertex set is closed under the differences the search looks for -- and "
       "it does not survive the motion. So jb_x's row X2, 'twelve-around-one "
       "closes ONLY at a = 0', is a GEOMETRIC FACT and not an artifact of the "
       "held list; what its wording gets wrong is attributing to the whole "
       "structure something true only of the neighbour-neighbour half, since "
       "the centre star closes at every phase (R3). TWO-SIDED: a = 0 must "
       "regenerate all 24, and every other phase must regenerate none",
       regen[0.0][1] == 24
       and all(regen[a][1] == 0 for a in regen if a != 0.0)
       and all(regen[a][2] < 1e-12 for a in regen),
       "; ".join(f"a={a:.3f}: {regen[a][0]} contacts, {regen[a][1]} "
                 f"neighbour-neighbour, residual {regen[a][2]:.1e}"
                 for a in sorted(regen)),
       "24 neighbour contacts at a = 0 and none at any other phase"))

    # R5 -- and FCC13 is the only instance in this module
    verdicts = {}
    for t in X.build_topologies():
        verdicts[t.name] = max(residual(t, a) for a in GENERIC)
    out["verdicts"] = verdicts
    bad = [n for n, w in verdicts.items() if w > 1e-9]
    A(("R5  FCC13 IS THE ONLY INSTANCE HERE, AND THE REASON IS STRUCTURAL. "
       "Every other topology builds its contacts from lattice GENERATORS "
       "rather than reading them off geometry, so it has no coincidence to "
       "mistake for a structure -- all eight are invariant to 4e-15 across "
       "nine phases. That is the audit bead 43o asked for, and it comes back "
       "with one instance, already filed. It also says where the class can "
       "recur: wherever a correspondence is DISCOVERED by measuring positions "
       "rather than CONSTRUCTED from the lattice, which is the same "
       "position-versus-identity distinction this project has now been bitten "
       "by five times",
       bad == ["FCC13 (twelve-around-one)"],
       "coincidence-kind held sets: " + str(bad) + "; worst residual among the "
       "rest: "
       + f"{max(w for n, w in verdicts.items() if n not in bad):.2e}",
       "exactly one, and it is the one on the bead"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_hh -- read once and held: which holds are labellings, which luck")
    print("=" * 78)
    checks, out = gate()
    bad = 0
    for name, ok, got, want in checks:
        bad += 0 if ok else 1
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")

    print("\n  THE AUDIT TABLE")
    for n, w in sorted(out["verdicts"].items(), key=lambda kv: -kv[1]):
        kind = "COINCIDENCE" if w > 1e-9 else "labelling"
        print(f"   {n:34s} {w:.2e}  {kind}")
    print(f"   {'hinge pairing':34s} {0.0:.2e}  labelling")

    print()
    print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
    print("   * THE IDIOM IS SOUND AND ITS ONE WRITTEN JUSTIFICATION IS NOT.")
    print("     `_fcc13_contacts` holds its set 'for the same reason the hinge")
    print("     pairing is never re-read at a merge angle'. Measured here,")
    print("     those are OPPOSITE cases: the pairing is invariant and")
    print("     degenerates only AT the merge angles, while FCC13's neighbour")
    print("     contacts exist at a = 0 and nowhere else.")
    print("   * FCC13 IS AN a = 0 TOPOLOGY. Regenerating per phase yields a")
    print("     13-cell STAR, which SC7 already covers, so it would replace")
    print("     the topology rather than repair it. Restricting it and saying")
    print("     why is the honest option of the three bead 43o offers.")
    print("   * NO GATE ROW ANYWHERE DEPENDS ON FCC13 AT a != 0. jb_y skips it")
    print("     in all three of its call sites, jb_x excludes it from")
    print("     `live_ranks`, and the one row that names it (X2) asserts")
    print("     precisely the restriction R4 confirms.")
    print("   * THIS AUDIT IS NOT COMPLETE. It covers the held sets reachable")
    print("     from jb_x. The honeycomb list read at a = -30 is jb_ht's T3")
    print("     row; the weld corner maps are jb_ic's R2b and jb_cl's C2.")
    print("     Named so the reader knows they were not re-measured here.")
    print()
    print("  ALL CHECKS PASSED." if not bad else f"  {bad} CHECK(S) FAILED.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
