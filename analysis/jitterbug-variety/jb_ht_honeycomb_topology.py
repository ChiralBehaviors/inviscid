"""jb_ht -- Topology rewired onto the honeycomb, and what the array's mobility is.

THE REWIRING BEAD (inviscid-ia5). `jb_x_array_linkage.build_topologies` places
neighbours at 2*v[g], where two cuboctahedra touch at ONE VERTEX. That is the
wrong packing and everything computed on it is superseded. This file is the
replacement, gated. Design of record: T2
`inviscid/design-honeycomb-topology-rewiring`.

WHAT CHANGED, AND WHAT DELIBERATELY DID NOT

  THE CONTACT FORMAT DID NOT CHANGE. `Topology.contacts` is a tuple of
  (i, k, j, l) -- unit i's vertex label k is the same point as unit j's vertex
  label l -- which is already exactly one vertex identification. It was never
  the wrong shape. `build_topologies` simply emitted ONE per neighbour where
  the real packing has three across a triangular face and four (then two, then
  one, as it folds) across a square. `assemble_free`'s three-rows-per-contact
  loop is untouched.

  THE LIST IS READ AT a = -30, NOT AT a = 0. This is the decision that bites,
  and it INVERTS the idiom `jb_x._fcc13_contacts` documents ("Read at a = 0 and
  never re-read"). The identification count is phase dependent because the
  squares fold: 48 at a = 0, 36 at a = -30, 30 at a = -60. Only the a = -30
  list is valid everywhere. Row T3 holds the a = 0 list instead and measures
  what happens: twelve of the forty-eight break, separating to 2.000000 by
  a = -60. Holding them welds the folding squares shut and forbids the
  exchange, and the array would then report as rigid -- a plausible-looking
  wrong answer of exactly the kind this epic has already published twice.

  UNITS NOW CARRY A PHASE. `topo.phases[i]` is the unit's own offset: 0 on the
  VE sublattice, 60 on the hole cells, which is the b = a + 60 of the exchange
  carried as topology data. Every array configuration in jb_x and jb_z is now
  built through `unit_corners`/`unit_verts` so a two-sublattice topology cannot
  be silently driven as though it were one. The offsets default to zero, so
  every pre-existing topology is bit-identical -- jb_x, jb_y and jb_z all
  reproduce their baseline gate rows exactly, which is the check that says this
  was a no-op where it had to be.

  `dsites` NOW TAKES v AS WELL AS dv. `dsites(dv) = sites(dv)` was correct only
  because `sites` is LINEAR in v for every other kind. The honeycomb's origins
  are lam(v)*site with lam = 0.5*|v0 - v3|, which is not linear. Row T8 is the
  control, and it carries a trap worth reading: the old shortcut returns
  |dlam/da| where the truth is signed, so it is wrong for a > -30 and
  ACCIDENTALLY RIGHT for a < -30. A test placed at a = -50 would have passed
  and hidden it.

WHAT IT MEASURES: THE MOBILITY QUESTION (ia5 scope 2)

  The bead asks the decisive question -- "if the honeycomb is globally 1-DOF
  then every cell moves in lockstep and the medium breathes rather than
  propagates, which is fatal for the wave programme". It is not. A two-hole
  patch of fourteen cells has nullity 28 against six rigid-body modes, and the
  six are MEASURED to lie in the null space rather than subtracted off (memo R2:
  naive DOF counting has failed here once already).

  READ THIS WITH THE BOUNDARY IN MIND. These are finite patches with free
  surfaces, and the record already warns that a finite-patch mechanism count
  grows linearly with patch radius against cubic volume, i.e. it is a SURFACE
  count. The bulk statement is jb_hc's Bloch result -- Maxwell -6, rigid at
  generic k, exactly floppy on the six <110> lines -- and these numbers do not
  contradict it or replace it. What they establish is the negative: the array
  is not a single-DOF linkage, so the fatal case the bead named is ruled out.

WHAT IS STILL OPEN: global consistency across cells (ia5 scope 3). Everything
here is built from one cell's neighbour list applied at each site; that the
local closure conditions are mutually satisfiable is verified by T2 on real
patches, but no proof is offered that an unbounded array closes.

Canonical prose state: T2 `inviscid/qvf-epic-consolidated-state`.
"""

from __future__ import annotations

import pickle
import sys

import numpy as np

import jb_hc_honeycomb as HC
import jb_w_honeycomb as W
from jb_x_array_linkage import (PAIRS, assemble_free, build_topologies,
                                dverts_exact, lattice_from_verts, rank_of,
                                unit_corners, verts)

# ---------------------------------------------------------------------------
# Thresholds, re-declared locally (house rule), each priced from a number this
# run prints and named in the comment beside it.
# ---------------------------------------------------------------------------

TOL_GEOM = 1e-12        # held-contact residual, measured <= 1.02e-15
TOL_LAM = 1e-12         # three-route lambda agreement, measured 4.44e-16
TOL_DSITES = 1e-8       # analytic vs finite difference, measured <= 2.64e-10
TOL_RIGID = 1e-12       # ||A @ rigid||, measured <= 7.5e-16
BREAK_MIN = 1.0         # the a=0 list's worst break, measured 2.000000
SHORTCUT_MIN = 1e-2     # linear-shortcut error above -30, measured >= 0.0140

#: Phases the held list is checked at, and a SECOND ABSOLUTE incommensurate arm.
#: The second grid is written as absolute numbers, never derived from the first.
PHASES = (0.0, -7.5, -15.0, -22.5, -30.0, -37.5, -45.0, -52.5, -60.0)
PHASES_ALT = (-3.7, -11.3, -19.1, -26.9, -33.1, -41.7, -48.3, -56.9)

#: Where the dsites shortcut MUST fail (a > -30, dlam/da < 0) and where it
#: silently agrees (a < -30). Both are gated, because the second is the trap.
SHORTCUT_FAIL_AT = (-2.0, -5.0, -10.0, -20.0)
SHORTCUT_AGREE_AT = (-40.0, -50.0, -55.0)

REF = W.HONEYCOMB_REF_PHASE


# ==========================================================================
# T0-T1: the primitives, and the census
# ==========================================================================

def t0_lambda():
    """Three independent routes to the lattice parameter. jb_w has the closed
    form, jb_hc measures it from verts and gates it against the project's own
    fold table, and jb_x needs a v-parameterised form because `sites` receives
    v. Three implementations exist for real reasons; this row is what stops
    them drifting."""
    worst_hc = worst_x = 0.0
    for a in np.linspace(0.0, -60.0, 61):
        lam = W.lam(a)
        worst_hc = max(worst_hc, abs(lam - HC.lattice(a)))
        worst_x = max(worst_x, abs(lam - lattice_from_verts(verts(a))))
    return dict(hc=worst_hc, x=worst_x, n=61)


def t1_census():
    counts = []
    for site in ((0, 0, 0), (1, 1, 1), (2, 0, 0), (-1, 1, -1)):
        tri, sq = W.neighbours(site)
        counts.append((len(tri), len(sq)))
    return dict(counts=counts, ok=all(c == (8, 6) for c in counts),
                n=len(counts))


# ==========================================================================
# T2-T4: the held list, and the control that says why it is read at -30
# ==========================================================================

def _held_residual(topo, contacts, phases):
    """Worst separation of a HELD identification list, over `phases`."""
    worst = 0.0
    for a in phases:
        origins = topo.sites(verts(a))
        xs = [unit_corners(a, topo, i) + origins[i] for i in range(topo.n)]
        for (i, k, j, l) in contacts:
            (fa, ca), (fb, cb) = PAIRS[k][0], PAIRS[l][0]
            worst = max(worst, float(np.linalg.norm(
                xs[i][fa][ca] - xs[j][fb][cb])))
    return worst


def t2_held(topos):
    rows = []
    for t in topos:
        rows.append((t.name,
                     _held_residual(t, t.contacts, PHASES),
                     _held_residual(t, t.contacts, PHASES_ALT)))
    return dict(rows=rows, worst=max(max(r[1], r[2]) for r in rows),
                n=len(rows))


def t3_read_at_zero(topos):
    """CONTROL. Build the same topologies with the list read at a = 0 and hold
    it. It must BREAK, and the row prints by how much."""
    rows = []
    for t in topos:
        if t.n < 2:
            continue
        at0 = W.honeycomb_identifications(t.lattice_sites, a=0.0)
        rows.append((t.name, len(at0), len(t.contacts),
                     _held_residual(t, at0, PHASES)))
    return dict(rows=rows, worst=max(r[3] for r in rows) if rows else 0.0,
                n=len(rows))


def t4_counts():
    """The identification count for ONE interior cell, at three phases."""
    out = {}
    for a in (0.0, REF, -60.0):
        cs = W.honeycomb_contacts((0, 0, 0), a)
        out[a] = (len(cs), sum(len(c[2]) for c in cs))
    invariant = all(
        sum(len(c[2]) for c in W.honeycomb_contacts((0, 0, 0), a))
        == out[REF][1] for a in PHASES_ALT)
    return dict(counts=out, invariant=invariant, n=len(PHASES_ALT))


def t5_old_topology():
    """CONTROL, kept from jb_hc's H1: what the topology being replaced builds."""
    return dict(old=W.wrong_packing_shared(0.0, 0),
                tri=len(W.honeycomb_contacts((0, 0, 0), REF)[0][2]))


# ==========================================================================
# T6-T8: phases default to zero; dsites needs v
# ==========================================================================

def t6_phases(topos):
    legacy = build_topologies()
    legacy_all_zero = all(set(t.phases) == {0.0} for t in legacy)
    hc_two = all(set(t.phases) <= {0.0, 60.0} for t in topos)
    hc_uses_both = any(set(t.phases) == {0.0, 60.0} for t in topos)

    # A Topology in the pre-change PICKLE FORMAT, i.e. without `phases`. This
    # is not hypothetical: jb_cache's argument trace replays last run's call
    # list into fresh workers and outlives a source change by design, so adding
    # an attribute to a persisted class breaks every prefetch task until the
    # class can read the old format. It cost one full jb_z run to find.
    revived, revived_ok = None, False
    for t in legacy:
        state = dict(t.__dict__)
        state.pop("phases", None)
        obj = object.__new__(type(t))
        obj.__setstate__(state)
        revived = obj
        revived_ok = (set(obj.phases) == {0.0} and len(obj.phases) == obj.n)
        if not revived_ok:
            break
    round_trip = all(
        set(pickle.loads(pickle.dumps(t)).phases) == set(t.phases)
        for t in list(legacy) + list(topos))
    return dict(legacy_zero=legacy_all_zero, n_legacy=len(legacy),
                hc_two=hc_two, hc_both=hc_uses_both, n_hc=len(topos),
                revived=revived_ok, revived_n=revived.n if revived else 0,
                round_trip=round_trip)


def t7_dsites(topos):
    worst = 0.0
    h = 1e-6
    for t in topos:
        for a in (-10.0, -30.0, -50.0, -19.1):
            fd = (t.sites(verts(a + h)) - t.sites(verts(a - h))) / (2 * h)
            an = t.dsites(dverts_exact(a), verts(a))
            worst = max(worst, float(np.abs(an - fd).max()))
    return dict(worst=worst, n=len(topos) * 4)


def t8_shortcut(topos):
    """The linear shortcut `sites(dv)`, which `dsites` used to be. It returns
    |dlam/da| where the truth is signed, so it is wrong above a = -30 and
    accidentally right below it."""
    t = topos[-1]
    h = 1e-6
    fail, agree = [], []
    for a in SHORTCUT_FAIL_AT:
        fd = (t.sites(verts(a + h)) - t.sites(verts(a - h))) / (2 * h)
        fail.append(float(np.abs(t.sites(dverts_exact(a)) - fd).max()))
    for a in SHORTCUT_AGREE_AT:
        fd = (t.sites(verts(a + h)) - t.sites(verts(a - h))) / (2 * h)
        agree.append(float(np.abs(t.sites(dverts_exact(a)) - fd).max()))
    return dict(fail_min=min(fail), agree_max=max(agree),
                n_fail=len(fail), n_agree=len(agree))


# ==========================================================================
# T9-T10: rigid modes, measured. Then the mobility.
# ==========================================================================

def _rigid_generators(topo, a):
    n = topo.n
    origins = topo.sites(verts(a))
    xs = [unit_corners(a, topo, i) + origins[i] for i in range(n)]
    G = np.zeros((6, 48 * n))
    for i in range(n):
        for p in range(8):
            for d in range(3):
                G[d, 48 * i + 24 + 3 * p + d] = 1.0
            c = xs[i][p].mean(axis=0)
            for d in range(3):
                e = np.zeros(3)
                e[d] = 1.0
                G[3 + d, 48 * i + 3 * p:48 * i + 3 * p + 3] += e
                G[3 + d, 48 * i + 24 + 3 * p:48 * i + 24 + 3 * p + 3] += np.cross(e, c)
    return G


def t9_mobility(topos):
    rows = []
    worst_rigid = 0.0
    for t in topos:
        A = assemble_free(REF, t)
        r, _s = rank_of(A)
        nullity = 48 * t.n - r
        G = _rigid_generators(t, REF)
        Q, _ = np.linalg.qr(G.T)
        worst_rigid = max(worst_rigid, float(np.linalg.norm(A @ Q)))
        rows.append(dict(name=t.name, n=t.n, contacts=len(t.contacts),
                         rows=A.shape[0], cols=A.shape[1], rank=r,
                         nullity=nullity, internal=nullity - 6,
                         rigid_rank=int(np.linalg.matrix_rank(G))))
    return dict(rows=rows, worst_rigid=worst_rigid, n=len(rows))


# ==========================================================================
# THE GATE
# ==========================================================================

def gate(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9):
    checks = []
    R = checks.append

    R(("T0  three independent routes to lambda agree to machine zero",
       t0["n"] > 0 and max(t0["hc"], t0["x"]) < TOL_LAM,
       f"jb_hc {t0['hc']:.2e}, jb_x {t0['x']:.2e} over {t0['n']} phases",
       f"< {TOL_LAM:.0e}"))
    R(("T1  every cell has 8 triangular + 6 square neighbours",
       t1["n"] > 0 and t1["ok"], f"{t1['counts'][0]} at {t1['n']} sites",
       "(8, 6)"))

    R(("T2  the list read at a = -30 is EXACT at every phase of the exchange",
       t2["n"] > 0 and t2["worst"] < TOL_GEOM,
       f"worst {t2['worst']:.2e} over {t2['n']} topologies, two grids",
       f"< {TOL_GEOM:.0e}"))
    R(("T3  CONTROL: the list read at a = 0 BREAKS -- this is why -30 -- "
       "CAN FAIL", t3["n"] > 0 and t3["worst"] > BREAK_MIN,
       f"worst {t3['worst']:.6f} over {t3['n']} topologies", f"> {BREAK_MIN}"))
    R(("T4  one interior cell: 36 identifications at -30, phase invariant; "
       "48 at a = 0 and 30 at a = -60 are what T3 is about",
       t4["invariant"] and t4["counts"][REF][1] == 36
       and t4["counts"][0.0][1] == 48 and t4["counts"][-60.0][1] == 30,
       f"{t4['counts'][0.0][1]} / {t4['counts'][REF][1]} / "
       f"{t4['counts'][-60.0][1]}", "48 / 36 / 30"))
    R(("T5  CONTROL: the topology being replaced shares ONE vertex per "
       "neighbour, against a shared face's three",
       t5["old"] == 1 and t5["tri"] == 3,
       f"old {t5['old']}, honeycomb triangular {t5['tri']}", "1 vs 3"))

    R(("T6  every LEGACY topology has phase offset zero -- this is what makes "
       "the rewiring a no-op for them",
       t6["n_legacy"] > 0 and t6["legacy_zero"],
       f"all {t6['n_legacy']} legacy topologies", "all 0.0"))
    R(("T6  honeycomb topologies carry exactly the two sublattice phases, and "
       "actually use both", t6["hc_two"] and t6["hc_both"],
       f"{t6['n_hc']} topologies", "{0.0, 60.0}"))

    R(("T6  a Topology in the PRE-CHANGE pickle format revives with zero "
       "offsets instead of raising -- jb_cache's trace outlives a source "
       "change, and this cost a full jb_z run to find",
       t6["revived"] and t6["round_trip"],
       f"revived n={t6['revived_n']} all 0.0, round trip ok", "no AttributeError"))

    R(("T7  dsites(dv, v) matches the finite difference",
       t7["n"] > 0 and t7["worst"] < TOL_DSITES,
       f"worst {t7['worst']:.2e} over {t7['n']} cases",
       f"< {TOL_DSITES:.0e}"))
    R(("T8  CONTROL: the old linear shortcut sites(dv) is WRONG above a = -30 "
       "-- CAN FAIL", t8["n_fail"] > 0 and t8["fail_min"] > SHORTCUT_MIN,
       f"min error {t8['fail_min']:.6f} over {t8['n_fail']} phases",
       f"> {SHORTCUT_MIN}"))
    R(("T8  THE TRAP, gated so nobody moves the test: the same shortcut is "
       "ACCIDENTALLY RIGHT below a = -30",
       t8["n_agree"] > 0 and t8["agree_max"] < TOL_DSITES,
       f"max error {t8['agree_max']:.2e} over {t8['n_agree']} phases",
       f"< {TOL_DSITES:.0e}"))

    R(("T9  the six rigid-body modes are MEASURED in the null space, never "
       "subtracted off (memo R2)",
       t9["n"] > 0 and t9["worst_rigid"] < TOL_RIGID
       and all(r["rigid_rank"] == 6 for r in t9["rows"]),
       f"||A @ rigid|| {t9['worst_rigid']:.2e}, rank 6 everywhere",
       f"< {TOL_RIGID:.0e}"))
    two = [r for r in t9["rows"] if "2HOLE" in r["name"]]
    R(("T9  THE MOBILITY ANSWER: the multi-hole patch is NOT globally 1-DOF "
       "-- the case ia5 named as fatal is ruled out",
       len(two) == 1 and two[0]["internal"] > 1,
       f"{two[0]['n']} cells, nullity {two[0]['nullity']}, "
       f"{two[0]['internal']} internal DOF" if two else "no patch",
       "> 1"))
    R(("T9  NON-VACUITY, TWO-SIDED: internal DOF is neither collapsed nor "
       "uncoupled (6 per free unit would mean no contact does anything)",
       len(two) == 1 and 1 < two[0]["internal"] < 6 * two[0]["n"],
       f"{two[0]['internal']} vs 6n = {6 * two[0]['n']}" if two else "n/a",
       "1 < dof < 6n"))
    one = [r for r in t9["rows"] if r["n"] == 1]
    R(("T9  CONTROL: the single-unit topology reproduces jb_x's own documented "
       "control -- rank 36, 6 internal DOF",
       len(one) == 1 and one[0]["rank"] == 36 and one[0]["internal"] == 6,
       f"rank {one[0]['rank']}, {one[0]['internal']} internal DOF"
       if one else "n/a", "36, 6"))

    print()
    print("=" * 78)
    print(f"GATE  {len(checks)} rows")
    print("=" * 78)
    for name, ok, val, crit in checks:
        print(f"  {'PASS' if ok else 'FAIL':4s}  {name:66s} {str(val):>26s} {str(crit):>18s}")

    print()
    print("  ROWS THAT EXIST ONLY TO STOP ANOTHER ROW BEING UNFALSIFIABLE:")
    print("   * T3, the a = 0 list. Without it 'read the list at -30' is a")
    print("     convention with no stated cost, and the next person to copy")
    print("     _fcc13_contacts' read-at-zero idiom silently welds the")
    print("     folding squares shut and gets a rigid array.")
    print("   * T5, the one-shared-vertex control. It is the whole reason")
    print("     this file exists and it stays visible in the gate that")
    print("     replaces the topology it measures.")
    print("   * T6's legacy row. 'The rewiring is a no-op for existing")
    print("     topologies' is an empirical claim, and the phase offsets")
    print("     being zero is what makes it true; jb_x, jb_y and jb_z")
    print("     reproducing their baseline rows exactly is the other half.")
    print("   * T6's pickle row. Adding an attribute to a class that jb_cache")
    print("     persists in its argument trace breaks every prefetch worker,")
    print("     silently degrading to a serial recompute rather than failing")
    print("     loudly. The row is here because the failure mode is invisible")
    print("     in the gate output it damages.")
    print("   * T8's SECOND row. The shortcut returns |dlam/da| where the")
    print("     truth is signed, so it is wrong above a = -30 and right")
    print("     below it. Gating only the failure would let someone move the")
    print("     test to a = -50, watch it pass, and delete the fix.")
    print("   * T9's rigid-mode row. Internal DOF is a subtraction, and memo")
    print("     R2 records that naive DOF counting has already failed here")
    print("     once; the six modes are measured into the null space first.")
    print("   * T9's single-unit row, which reproduces a control jb_x already")
    print("     documented (rank 36, 6 internal DOF) through the new code")
    print("     path.")
    print()
    print("  ROWS DELETED RATHER THAN FIXED: a row asserting the honeycomb")
    print("  array is over-constrained because Maxwell counting says so. It")
    print("  is over-constrained in the BULK (jb_hc: 18 dof against 24")
    print("  constraints per unit cell) and these finite patches are not, and")
    print("  the difference is entirely the free surface. Counting rows")
    print("  against columns on a patch and calling the answer the medium's")
    print("  is the mistake memo R2 is about.")
    print()
    print("  A ROW DELIBERATELY NOT BUILT: global consistency of the contact")
    print("  list over an unbounded array (ia5 scope 3). Every topology here")
    print("  is built by applying one cell's neighbour list at each site, and")
    print("  T2 verifies the result closes on real patches up to fifteen")
    print("  cells. That is evidence, not proof, and no row here claims more.")
    print()
    print("  WHAT THIS FILE DOES NOT DO: it does not rewire jb_z's CRANK path")
    print("  onto the honeycomb. jb_z's gate still drives the legacy")
    print("  topologies, unchanged and still green; what this rewiring buys")
    print("  it is that `unit_corners` is now the only way a configuration")
    print("  gets built, so a honeycomb topology handed to it will be driven")
    print("  correctly rather than silently as a single sublattice.")

    failed = [n for n, ok, _v, _c in checks if not ok]
    print()
    if failed:
        print(f"  !! {len(failed)} CHECK(S) FAILED -- this is a bug report, not a")
        print("     measurement. Nothing above may enter the record.")
        for n in failed:
            print(f"       - {n}")
        return 1
    print("  ALL CHECKS PASSED.")
    return 0


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_ht -- Topology on the honeycomb, and the array's mobility")
    print("=" * 78)
    print("  Neighbours share FACES, so a contact is three vertex")
    print("  identifications across a triangle or four-then-two-then-one")
    print("  across a folding square -- not the single vertex the topology")
    print("  being replaced emitted. The list is read at a = -30 because that")
    print("  is the only phase whose list is valid at every other one.")
    topos = W.build_honeycomb_topologies()
    t0 = t0_lambda()
    t1 = t1_census()
    t2 = t2_held(topos)
    t3 = t3_read_at_zero(topos)
    t4 = t4_counts()
    t5 = t5_old_topology()
    t6 = t6_phases(topos)
    t7 = t7_dsites(topos)
    t8 = t8_shortcut(topos)
    t9 = t9_mobility(topos)

    print()
    print("-" * 78)
    print("  MOBILITY")
    print("-" * 78)
    print(f"  {'topology':40s} {'n':>3s} {'cont':>5s} {'rank':>5s} "
          f"{'null':>5s} {'internal':>9s}")
    for r in t9["rows"]:
        print(f"  {r['name']:40s} {r['n']:3d} {r['contacts']:5d} "
              f"{r['rank']:5d} {r['nullity']:5d} {r['internal']:9d}")
    print("  internal = nullity - 6, and the 6 are measured into the null")
    print("  space by T9, not assumed. Free-surface patches: see the module")
    print("  docstring before reading these as the medium's numbers.")
    return gate(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9)


if __name__ == "__main__":
    with np.errstate(all="ignore"):
        sys.exit(main())
