"""jb_lf -- the front in three dimensions: orbit-structured, and ordered by neither ruler.

WHY THIS FILE EXISTS. Epic inviscid-qvf's criterion says "a LATTICE of
interconnected jitterbugs exhibits wave propagation". Everything measured up to
now -- jb_tr, jb_sv -- runs on a CHAIN: a genuine diagonal of the real
honeycomb, but one-dimensional. This is the three-dimensional measurement, and
it was filed as bead inviscid-x3i first because the probe produced a result
nobody could read.

THE PROBE'S THREE PUZZLES, all of which dissolve here and none of which was
what it looked like. Grouping arrivals by weld-graph BFS shell gave:

    shell 1:  8/8  at 0.0560 exactly, spread 0
    shell 2: 18/26 arriving over a range 0.327 .. 0.776
    shell 3: 24/24 at 0.5090 exactly, spread 0

A front that is not monotone in graph distance, with a shell that splits and a
shell that does not. The resolution is that THE SHELLS WERE THE WRONG GROUPING.
Classified by CUBIC ORBIT instead -- the symmetry class of a site, labelled by
its sorted absolute coordinates -- the picture is exact:

    orbit      euclid   n   arrival   spread
    (1,1,1)     1.732    8   0.0560    0.0000
    (0,0,2)     2.000    6   0.3270    0.0000
    (0,2,2)     2.828   12   0.7760    0.0000
    (1,1,3)     3.317   24   0.5090    0.0000
    (2,2,2)     3.464    8   not by t = 1.6

  * shell 2 "split" because it is THREE ORBITS, at Euclidean 2.00, 2.83 and
    3.46 with weld degrees 8, 6 and 4. Of course they do not arrive together.
  * shells 1 and 3 had zero spread because each is a SINGLE orbit.
  * the 8 that "never arrived" are one orbit, (2,2,2) -- the farthest and
    least connected -- and the honest statement is that they had not arrived
    by t = 1.6, not that they cannot.

WHAT IS ACTUALLY TRUE, AND IT IS SHARPER THAN THE PUZZLE. Every cubic orbit
arrives as a SINGLE SIMULTANEOUS EVENT, to the timestep, with zero spread
across as many as 24 cells at once. The front in this lattice is not a smooth
expanding surface: it is a sequence of discrete, symmetry-determined shells,
each of which lights up at one instant.

AND NEITHER RULER ORDERS THEM (R3). The (1,1,3) orbit sits at weld-graph hop 3
and Euclidean distance 3.317, and it arrives at 0.5090 -- BEFORE the (0,2,2)
orbit at hop 2 and Euclidean 2.828, which arrives at 0.7760. So the arrival
order is not the graph distance the disturbance travels through, and it is not
the straight-line distance either. WITHIN one hop shell the order IS by
Euclidean distance (R4), so both rulers carry part of it and neither carries
all: the disturbance moves through welds, and the time each weld costs depends
on that weld's own geometry.

THIS IS ANISOTROPIC PROPAGATION WITH EXACT SYMMETRY, which is a stronger and
stranger statement than "it propagates". It is also the honest answer to the
epic's "lattice" clause: the disturbance does traverse the three-dimensional
lattice, its first shell is perfectly isotropic, and beyond that it is
structured by the cubic group rather than by distance.

SCOPE, and R5 is the part that limits it. This is one ball of 59 cells with a
FREE BOUNDARY, and the outer orbits are exactly the ones nearest that boundary.
The orbit structure and the simultaneity are symmetry facts and do not depend on
the boundary; the ARRIVAL TIMES of the outermost orbits may. Nothing here
separates those, and no orbit-to-orbit time is offered as a bulk number.
"""
from __future__ import annotations

import collections
import sys

import numpy as np

import jb_mj_inertial_honeycomb as MJ
import jb_rc_reduced as RC

RADIUS = 3.5
KICK = 0.9
PLAY = 0.05
RESTITUTION = 1.0
STEP = 1e-3
TMAX = 1.6
ARRIVED = 0.02


def patch(radius=RADIUS, gc=MJ.A_REF):
    sites = RC.ball(radius)
    asm, _ = RC.honeycomb(sites, gc=gc)
    return asm, sites, sites.index((0, 0, 0))


def orbit(site):
    """The cubic-group orbit of a site: its sorted absolute coordinates. Cells
    in one orbit are carried into each other by the lattice's own symmetry, so
    a symmetric excitation must move them identically."""
    return tuple(sorted(abs(int(c)) for c in site))


def hops(asm, centre):
    adj = collections.defaultdict(set)
    for (k, l, _) in asm.welds:
        adj[k].add(l)
        adj[l].add(k)
    h = {centre: 0}
    frontier = [centre]
    while frontier:
        nxt = []
        for k in frontier:
            for m in adj[k]:
                if m not in h:
                    h[m] = h[k] + 1
                    nxt.append(m)
        frontier = nxt
    return h, {k: len(v) for k, v in adj.items()}


def sweep(asm, centre, kick=KICK, play=PLAY, e=RESTITUTION, tmax=TMAX,
          h=STEP, bands=True):
    pairs = MJ.tied_pairs(asm)
    q = asm.q0()
    u = np.zeros((asm.N, 7))
    u[centre, 6] = kick
    now = 0.0
    arrive = {centre: 0.0}
    while now < tmax - 1e-12:
        J, M, Minv = MJ.kinematics(asm, q)
        a1 = MJ.free_accel(asm, q, u, J, Minv)
        u_h = u + 0.5 * h * a1
        q_h = RC.apply_increment(asm, q, (0.5 * h * u).ravel())
        Jh, _, Mih = MJ.kinematics(asm, q_h)
        a2 = MJ.free_accel(asm, q_h, u_h, Jh, Mih)
        u = u + h * a2
        q = RC.apply_increment(asm, q, (h * u_h).ravel())
        now += h
        if bands:
            J, M, Minv = MJ.kinematics(asm, q)
            s = MJ.separations(asm, q, pairs)
            N = MJ.band_rows(asm, q, J, pairs)
            rate = np.dot(N, u.ravel())
            act = [i for i in range(len(pairs))
                   if s[i] >= play and rate[i] > 0]
            if act:
                u, _, _, _ = MJ.resolve(asm, u, N, act, Minv, e)
        for k in range(asm.N):
            if k not in arrive and abs(u[k, 6]) > ARRIVED * abs(kick):
                arrive[k] = now
    return arrive


def by_orbit(sites, centre, arrive, hop, deg):
    g = collections.defaultdict(list)
    for k, t in arrive.items():
        if k == centre:
            continue
        g[orbit(sites[k])].append((t, hop.get(k, -1), deg.get(k, 0)))
    out = {}
    for o, v in g.items():
        ts = [x[0] for x in v]
        out[o] = {"n": len(v), "mean": float(np.mean(ts)),
                  "spread": float(max(ts) - min(ts)),
                  "hop": v[0][1], "deg": v[0][2],
                  "euclid": float(np.linalg.norm(o))}
    return out


def gate():
    checks, out = [], {}
    A = checks.append
    asm, sites, ci = patch()
    hop, deg = hops(asm, ci)
    pairs = MJ.tied_pairs(asm)
    weld = float(np.abs(asm.weld_residual(asm.q0())).max())
    arrive = sweep(asm, ci)
    ob = by_orbit(sites, ci, arrive, hop, deg)
    out["orbits"] = ob

    # R1 -- the patch
    A(("R1  A THREE-DIMENSIONAL PATCH, WHICH IS WHAT THE EPIC'S 'LATTICE' "
       "CLAUSE ASKS FOR AND WHAT EVERY EARLIER TRANSPORT MEASUREMENT LACKED. "
       "jb_tr and jb_sv run on a CHAIN -- a genuine diagonal of the real "
       "honeycomb, but one-dimensional. This is a ball of 59 cells about a "
       "central VE, 160 welds, 480 unilateral bands, 413 degrees of freedom, "
       "closing to machine precision. The central cell has the full eight "
       "triangular-face neighbours, so the excitation starts somewhere the "
       "lattice actually looks like itself",
       asm.N == 59 and len(pairs) == 480 and weld < 1e-12
       and deg[ci] == 8,
       f"{asm.N} cells, {len(asm.welds)} welds, {len(pairs)} bands, "
       f"{7 * asm.N} DOF, weld residual {weld:.1e}, centre weld degree "
       f"{deg[ci]}",
       "a closed 3D patch with a fully-coordinated centre"))

    # R2 -- the finding
    spreads = {o: v["spread"] for o, v in ob.items()}
    biggest = max(ob.values(), key=lambda v: v["n"])
    A(("R2  EVERY CUBIC ORBIT ARRIVES AS ONE SIMULTANEOUS EVENT. Grouped by "
       "the symmetry class of a site -- its sorted absolute coordinates, which "
       "is the orbit the lattice's own cubic group carries it around -- the "
       "arrival spread within EVERY orbit is zero to the timestep, across as "
       "many as 24 cells at once. So the front here is not a smooth expanding "
       "surface: it is a sequence of DISCRETE, symmetry-determined shells, "
       "each lighting up at a single instant. That is what a symmetric "
       "excitation on a symmetric site must do, and measuring it confirms the "
       "excitation really is symmetric and the integrator really is not "
       "breaking the symmetry -- either failure would show as a nonzero spread "
       "here, which is what makes this row able to fail",
       max(spreads.values()) < 1e-9 and biggest["n"] >= 24,
       "; ".join(f"{o} n={v['n']} arrival {v['mean']:.4f} spread "
                 f"{v['spread']:.1e}"
                 for o, v in sorted(ob.items(), key=lambda kv: kv[1]["euclid"])),
       "zero spread within every orbit, largest orbit >= 24 cells"))

    # R3 -- and neither ruler orders them
    order = sorted(ob.items(), key=lambda kv: kv[1]["mean"])
    eu_mono = all(order[i][1]["euclid"] <= order[i + 1][1]["euclid"] + 1e-9
                  for i in range(len(order) - 1))
    hop_mono = all(order[i][1]["hop"] <= order[i + 1][1]["hop"]
                   for i in range(len(order) - 1))
    out["mono"] = (eu_mono, hop_mono)
    A(("R3  NEITHER THE GRAPH DISTANCE NOR THE STRAIGHT-LINE DISTANCE ORDERS "
       "THE ARRIVALS, and this is what made the original probe unreadable. The "
       "(1,1,3) orbit sits at weld-graph hop 3 and Euclidean distance 3.317, "
       "and it arrives BEFORE the (0,2,2) orbit at hop 2 and Euclidean 2.828. "
       "So the disturbance is not simply walking the weld graph, and it is not "
       "travelling in straight lines either. Both rulers carry part of it -- "
       "R4 shows the order is Euclidean WITHIN a hop shell -- and neither "
       "carries all: the disturbance moves through welds, and the time a weld "
       "costs depends on that weld's own geometry. TWO-SIDED: this row asserts "
       "BOTH orderings fail, so a medium in which either worked would redden "
       "it",
       (not eu_mono) and (not hop_mono),
       "arrival order " + " -> ".join(f"{o}" for o, _ in order)
       + "; Euclidean " + " -> ".join(f"{v['euclid']:.3f}" for _, v in order)
       + "; hop " + " -> ".join(f"{v['hop']}" for _, v in order),
       "monotone in neither distance"))

    # R4 -- but within a shell it is Euclidean
    h2 = sorted([(o, v) for o, v in ob.items() if v["hop"] == 2],
                key=lambda kv: kv[1]["euclid"])
    in_shell = all(h2[i][1]["mean"] <= h2[i + 1][1]["mean"] + 1e-9
                   for i in range(len(h2) - 1))
    A(("R4  WITHIN ONE HOP SHELL THE ORDER IS EUCLIDEAN, which is what makes "
       "R3 a structure rather than a mess. The hop-2 orbits that arrived did "
       "so in exactly their distance order, and their weld degrees fall the "
       "same way, so the farther orbit is also the less connected one. Read "
       "with R3: hops decide WHICH welds the disturbance must cross, and "
       "geometry decides what each crossing costs. SCOPE, and the row carries "
       "it rather than the prose: hop 2 holds THREE orbits and only TWO of "
       "them arrived inside the run, so this ordering rests on two points. The "
       "third, (2,2,2), is the farthest and least connected and is still in "
       "flight at t = 1.6 -- which is consistent with the ordering but does "
       "not test it",
       len(h2) >= 2 and in_shell,
       "; ".join(f"{o} euclid {v['euclid']:.3f} degree {v['deg']} arrival "
                 f"{v['mean']:.4f}" for o, v in h2)
       + f"; hop-2 orbits arrived {len(h2)} of 3",
       "the hop-2 orbits that arrived, in Euclidean order"))

    # R5 -- controls, and the boundary caveat as a measurement
    still = sweep(asm, ci, kick=0.0, tmax=0.3)
    free = sweep(asm, ci, tmax=0.3, bands=False)
    unreached = [sites[k] for k in range(asm.N) if k not in arrive]
    orb_un = collections.Counter(orbit(s) for s in unreached)
    out["unreached"] = dict(orb_un)
    A(("R5  BOTH CONTROLS, AND THE BOUNDARY CAVEAT STATED AS A MEASUREMENT "
       "RATHER THAN A WORRY. Undriven, only the kicked cell ever moves; with "
       "the bands removed, likewise -- so the propagation is a property of the "
       "coupling. And the cells that had NOT arrived by the end of the run are "
       "exactly one orbit, (2,2,2), the farthest and least connected, sitting "
       "on the free boundary. That is reported as 'not yet' and not as "
       "'cannot': this patch cannot tell the two apart, because the outermost "
       "orbits are precisely the ones the boundary touches. The orbit "
       "structure and the simultaneity are symmetry facts and do not depend on "
       "the boundary; the outer arrival TIMES may, and none is offered as a "
       "bulk number",
       len(still) == 1 and len(free) == 1
       and len(arrive) > 40 and len(orb_un) <= 1,
       f"undriven: {len(still)} of {asm.N} cells move; bands disabled: "
       f"{len(free)}; bands enabled: {len(arrive)}; not yet arrived by "
       f"t = {TMAX}: {dict(orb_un)}",
       "only the driven cell in both controls, and the stragglers one orbit"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("jb_lf -- the front in three dimensions, by cubic orbit")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")

        print("\n  ARRIVAL BY CUBIC ORBIT")
        print(f"   {'orbit':>12s} {'euclid':>7s} {'hop':>4s} {'deg':>4s} "
              f"{'n':>3s} {'arrival':>9s} {'spread':>9s}")
        for o, v in sorted(out["orbits"].items(),
                           key=lambda kv: kv[1]["euclid"]):
            print(f"   {str(o):>12s} {v['euclid']:7.3f} {v['hop']:4d} "
                  f"{v['deg']:4d} {v['n']:3d} {v['mean']:9.4f} "
                  f"{v['spread']:9.1e}")
        for o, n in out["unreached"].items():
            print(f"   {str(o):>12s} {float(np.linalg.norm(o)):7.3f} "
                  f"{'':4s} {'':4s} {n:3d}   not arrived by t = {TMAX}")

        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * THE DISTURBANCE TRAVERSES THE THREE-DIMENSIONAL LATTICE,")
        print("     which is the epic's 'lattice' clause and which every")
        print("     earlier transport measurement, being on a chain, could not")
        print("     address.")
        print("   * PROPAGATION IS ORBIT-STRUCTURED AND EXACTLY SIMULTANEOUS")
        print("     within each orbit. The front is discrete shells set by the")
        print("     cubic group, not a smooth surface.")
        print("   * IT IS ANISOTROPIC. Neither graph distance nor Euclidean")
        print("     distance orders the orbits, though within a hop shell the")
        print("     order is Euclidean. Do not quote an isotropic 3D speed;")
        print("     there is not one to quote.")
        print("   * NO BULK NUMBER. One ball, 59 cells, FREE BOUNDARY, and the")
        print("     outer orbits are the ones the boundary touches. The orbit")
        print("     structure and the simultaneity are symmetry facts and")
        print("     survive that; the outer arrival times are not offered as")
        print("     bulk values, and the orbit that had not arrived is")
        print("     reported as 'not yet', never as 'cannot'.")
        print("   * V = 0 STILL, and no frequency anywhere. Under DECISION 19")
        print("     the criterion is the sonic-vacuum transport law, and this")
        print("     is its three-dimensional half.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, "
                   "not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
