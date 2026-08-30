"""cluster -- a jitterbug cluster with an interior, and what a local impulse does.

WHY A CLUSTER AND NOT A CHAIN. A line of cells never has an interior: only 2 of
each cell's 8 faces are shared, so every vertex lies on an outward-facing face
however long the line gets. Measured on the 3-cell line: all 30 vertices are on
the boundary, and embedding it in the array by pinning what the array holds
leaves ZERO degrees of freedom. Pinning is the wrong embedding anyway -- in the
real bulk those neighbours move. A patch has to have cells to spare.

The smallest one that does is the centre cell plus all eight of its <111>
neighbours. Its centre cell's twelve vertices are then genuinely interior: each
lies on two of its own faces, and those two faces go to two DIFFERENT
neighbours, so the vertex is shared by THREE cells and carries six triangle
instances -- which is four distinct triangles once Fuller's congruence is taken
out. C3 gates exactly that.

THE IMPULSE IS THE TRIANGLE'S OWN 1-DOF. A jitterbug triangle does not spin
freely: Gray's V = EL cos(gamma) locks its radial position to its rotation, so
its motion has ONE parameter and the only free choice is handedness. An earlier
version of this experiment applied a PURE SPIN, which is not on that path at
all; the constraint projection then had to mangle it, and only 68 percent of the
direction survived against the 1-DOF impulse's 73.

WHAT IT SHOWS, and it is the point of the whole exercise. The constraint is
IMMEDIATE and the transport is NOT, and they are different things:

  * At t = 0, before any time has passed, 80 percent of the energy is already
    in the shell. One triangle cannot move without its neighbours moving, so
    making the impulse admissible distributes it instantly. A rigid constraint
    has infinite signal speed; there is no onset lag to find.
  * Then it takes TIME. The centre drains by a factor of eleven over t = 20 and
    climbs back by t = 30, with one shell cell holding half the total energy
    around t = 5.

Both at once, and neither is the other.
"""
from __future__ import annotations

import itertools as it
import sys

import numpy as np

from analysis.model import plates as Z
from analysis.model import cell as IC
SITES = list(it.product((-1, 1), repeat=3))
_XG = Z.corners(IC.GEN)
FDIR = np.array([_XG[f].mean(0) / np.linalg.norm(_XG[f].mean(0)) for f in range(8)])


def face_along(d):
    d = np.asarray(d, float)
    return int(np.argmax(FDIR @ (d / np.linalg.norm(d))))


def cluster(gc=0.0, sites=SITES):
    """Centre cell at gc, every <111> neighbour at gc+60. Welds face to face
    along each diagonal, with the corner correspondence restricted to the mating
    face's three corners -- searching all twelve lets a non-corner win, which is
    the defect jb_ic's R2b gates on the line."""
    gn = gc + 60.0
    sep = IC.ZC * (np.cos(np.radians(gc)) + np.cos(np.radians(gn)))
    L = sep / np.sqrt(3)
    cells = [(gc, np.zeros(3))] + [(gn, L * np.array(s, float)) for s in sites]
    ncb = len(cells)
    par = list(range(ncb * IC.NV))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    welded_faces = []
    for m, s in enumerate(sites):
        d = np.array(s, float)
        fc, fn = face_along(d), face_along(-d)
        A = IC.cell_verts(gc, np.zeros(3))
        B = IC.cell_verts(gn, cells[m + 1][1])
        bc = [IC.SLOT[(fn, c)] for c in range(3)]
        used, hit = set(), []
        for c in range(3):
            a = IC.SLOT[(fc, c)]
            k = min([b for b in bc if b not in used],
                    key=lambda b: float(np.linalg.norm(B[b] - A[a])))
            used.add(k)
            hit.append(k)
            ra, rb = find(a), find((m + 1) * IC.NV + k)
            if ra != rb:
                par[ra] = rb
        welded_faces.append((fc, fn, sorted(hit) == sorted(bc)))
    uniq, gid = {}, np.zeros((ncb, IC.NV), int)
    for k in range(ncb):
        for i in range(IC.NV):
            gid[k, i] = uniq.setdefault(find(k * IC.NV + i), len(uniq))
    P = np.zeros((len(uniq), 3))
    wr = np.zeros(len(uniq), bool)
    worst = 0.0
    for k, (g, o) in enumerate(cells):
        V = IC.cell_verts(g, o)
        for i in range(IC.NV):
            gg = gid[k, i]
            if wr[gg]:
                worst = max(worst, float(np.linalg.norm(P[gg] - V[i])))
            P[gg] = V[i]
            wr[gg] = True
    tris = [tuple(gid[k, IC.SLOT[(f, c)]] for c in range(3)) for k in range(ncb) for f in range(8)]
    bars = [(t[a], t[b]) for t in tris for a, b in ((0, 1), (1, 2), (0, 2))]
    return P, bars, tris, gid, ncb, worst, welded_faces


def dynamics(gc=0.0, face=1, dt=0.02, nsteps=1500, sample=125):
    P, bars, tris, gid, ncb, worst, _ = cluster(gc)
    n, nb = len(P), len(bars)
    mass = np.zeros(n)
    for t in tris:
        for i in t:
            mass[i] += 1.0 / 3.0
    BI = np.array([b[0] for b in bars])
    BJ = np.array([b[1] for b in bars])
    L2 = np.sum((P[BI] - P[BJ]) ** 2, axis=1)
    Minv = np.repeat(1.0 / mass, 3)
    cellv = [sorted(set(gid[k])) for k in range(ncb)]
    own = np.zeros(n)
    for k in range(ncb):
        for i in cellv[k]:
            own[i] += 1

    def gvec(q):
        d = q[BI] - q[BJ]
        return np.einsum('bi,bi->b', d, d) - L2

    def G(q):
        d = q[BI] - q[BJ]
        M = np.zeros((nb, 3 * n))
        r = np.arange(nb)
        for k in range(3):
            M[r, 3 * BI + k] = 2 * d[:, k]
            M[r, 3 * BJ + k] = -2 * d[:, k]
        return M

    def E(v):
        return 0.5 * float(np.sum(mass[:, None] * v * v))

    def shares(v):
        # each vertex's energy split among the cells owning it: no double count
        e = 0.5 * mass * np.sum(v * v, axis=1) / own
        return np.array([e[cellv[k]].sum() for k in range(ncb)])

    h = 1e-6
    v = np.zeros_like(P)
    for c in range(3):
        v[gid[0, IC.SLOT[(face, c)]]] = (Z.corners(gc + h)[face][c]
                                         - Z.corners(gc - h)[face][c]) / (2 * h)
    raw = v.copy()
    G0 = G(P)
    lam = np.linalg.lstsq(G0 @ (Minv[:, None] * G0.T), -(G0 @ v.reshape(-1)), rcond=None)[0]
    v = (v.reshape(-1) + Minv * (G0.T @ lam)).reshape(-1, 3)
    kept = float(np.sum(v * raw) / np.linalg.norm(v) / np.linalg.norm(raw))
    v = v / np.sqrt(2 * E(v))
    E0 = E(v)
    adm = float(np.linalg.norm(G(P) @ v.reshape(-1)))
    q = P.copy()
    hist = [(0.0, shares(v), E0)]
    for s in range(nsteps):
        q0 = q.copy()
        G0 = G(q0)
        qp = q + dt * v
        lam = np.zeros(nb)
        for _ in range(20):
            qn = qp + (0.5 * dt * dt * (Minv * (G0.T @ lam))).reshape(-1, 3)
            g = gvec(qn)
            if np.abs(g).max() < 1e-12:
                break
            lam -= np.linalg.lstsq(G(qn) @ (0.5 * dt * dt * (Minv[:, None] * G0.T)), g,
                                   rcond=None)[0]
        qn = qp + (0.5 * dt * dt * (Minv * (G0.T @ lam))).reshape(-1, 3)
        vh = v + (0.5 * dt * (Minv * (G0.T @ lam))).reshape(-1, 3)
        G1 = G(qn)
        mu = np.linalg.lstsq(G1 @ (Minv[:, None] * G1.T), -(G1 @ vh.reshape(-1)) / (0.5 * dt),
                             rcond=None)[0]
        v = vh + (0.5 * dt * (Minv * (G1.T @ mu))).reshape(-1, 3)
        q = qn
        if (s + 1) % sample == 0:
            hist.append(((s + 1) * dt, shares(v), E(v)))
    return dict(P=P, bars=bars, tris=tris, gid=gid, ncb=ncb, worst=worst, own=own,
                hist=hist, E0=E0, adm=adm, kept=kept)


def gate():
    A = []
    P, bars, tris, gid, ncb, worst, welds = cluster(0.0)
    n = len(P)
    L = np.array([np.linalg.norm(P[i] - P[j]) for i, j in bars])

    A.append(("C1  THE CLUSTER IS A REAL FRAMEWORK: nine cells, every strut a "
              "strut, and every shared vertex written twice and agreeing. The "
              "doubled writes are the free check -- but only on POSITION, which "
              "is why C2 exists as well",
              ncb == 9 and len(bars) == 24 * 9 and len(tris) == 8 * 9
              and abs(L.max() - IC.EL) < 1e-9 and abs(L.min() - IC.EL) < 1e-9 and worst < 1e-9,
              f"{ncb} cells, {n} vertices, {len(bars)} struts (24 x 9), "
              f"{len(tris)} triangles (8 x 9), struts {L.min():.6f}..{L.max():.6f}, "
              f"write disagreement {worst:.1e}",
              f"24 struts and 8 triangles per cell, all {IC.EL:.6f}"))

    A.append(("C2  EVERY WELD JOINS FACE TO FACE, in all eight diagonal "
              "directions. Restricting the corner search to the MATING FACE is "
              "the whole content of this row: searching all twelve of the "
              "neighbour's vertices lets a nearer non-corner win and welds three "
              "vertices that form no triangle. That defect survived a position "
              "agreement of 3e-15 on the line, because A REDUNDANCY CHECK ON "
              "POSITIONS DOES NOT VALIDATE IDENTITY",
              all(ok for _, _, ok in welds),
              f"{sum(1 for _, _, ok in welds if ok)}/8 welds land on the mating face",
              "8 of 8"))

    cnt = {}
    for t in tris:
        for i in t:
            cnt[i] = cnt.get(i, 0) + 1
    interior = sorted(i for i, c in cnt.items() if c == 6)
    centre = sorted(set(gid[0]))
    A.append(("C3  IT HAS AN INTERIOR, which a line never does. The centre "
              "cell's twelve vertices each lie on two of its faces, and those "
              "two faces go to two DIFFERENT neighbours, so each is shared by "
              "THREE cells and carries SIX triangle instances -- which is four "
              "DISTINCT triangles once congruence is removed, the count the "
              "medium actually has. Everything else in the cluster is boundary",
              len(interior) == 12 and set(interior) == set(centre)
              and set(cnt.values()) == {2, 6},
              f"{len(interior)} vertices on 6 triangle instances, and they are "
              f"exactly the centre cell's {len(centre)}; multiplicities present "
              f"{sorted(set(cnt.values()))}",
              "12 interior vertices, the centre cell's"))

    R = np.zeros((len(bars), 3 * n))
    for r, (i, j) in enumerate(bars):
        d = P[i] - P[j]
        d = d / np.linalg.norm(d)
        R[r, 3 * i:3 * i + 3] = d
        R[r, 3 * j:3 * j + 3] = -d
    s = np.linalg.svd(R, compute_uv=False)
    dof = 3 * n - int((s > s[0] * 1e-8).sum()) - 6
    A.append(("C4  AND IT STILL HAS ROOM TO MOVE. Nine cells joined through "
              "eight shared faces keep six internal freedoms per cell -- the "
              "sharing costs nothing, exactly as on the line",
              dof == 54,
              f"internal DOF {dof} over {ncb} cells ({dof / ncb:.2f} per cell)",
              "54, i.e. 6 per cell"))

    d = dynamics()
    A.append(("C5  THE IMPULSE IS THE TRIANGLE'S OWN 1-DOF AND IT IS "
              "ADMISSIBLE. A jitterbug triangle's radial position is locked to "
              "its rotation by V = EL cos(gamma), so its motion has ONE "
              "parameter and only handedness is free. A PURE SPIN is not on "
              "that path at all. ON THE 3-CELL LINE the 1-DOF impulse kept 73 "
              "percent of its direction through the constraint projection "
              "against a pure spin's 68; IN THIS CLUSTER it keeps far less, "
              "because eight neighbours constrain the triangle rather than one "
              "or two, and the surviving fraction is reported below rather than "
              "carried over from the line. V = 0, so energy is the only audit "
              "and RATTLE has to earn it",
              d["adm"] < 1e-10 and max(abs(E - d["E0"]) / d["E0"] for _, _, E in d["hist"]) < 1e-4,
              f"||G v|| = {d['adm']:.1e}; {d['kept'] * 100:.0f}% of the 1-DOF "
              f"direction survives projection; E conserved to "
              f"{max(abs(E - d['E0']) / d['E0'] for _, _, E in d['hist']):.1e}",
              "admissible, drift < 1e-4"))

    f0 = d["hist"][0][1] / d["hist"][0][1].sum()
    A.append(("C6  THE CONSTRAINT IS IMMEDIATE. At t = 0, before any time has "
              "passed, most of the energy is ALREADY in the shell: one triangle "
              "cannot move without its neighbours moving, so making the impulse "
              "admissible distributes it instantly. A rigid constraint has "
              "infinite signal speed and there is no onset lag to find "
              "anywhere. TWO-SIDED: a finite signal speed would show as the "
              "shell starting at zero, and it never does",
              f0[1:].sum() > 0.5 and f0[0] > 0,
              f"at t=0 the centre holds {f0[0]:.5f} and the shell {f0[1:].sum():.5f}",
              "shell already above half at t=0"))

    fr = np.array([h[1] / h[1].sum() for h in d["hist"]])
    ts = [h[0] for h in d["hist"]]
    lo = int(np.argmin(fr[:, 0]))
    A.append(("C7  AND THE TRANSPORT IS NOT. The centre DRAINS by an order of "
              "magnitude and then comes back, which is the part that takes "
              "time. Immediate coupling and slow transport are different "
              "things, and this medium has both at once",
              fr[lo, 0] < f0[0] / 5 and fr[-1, 0] > fr[lo, 0] * 1.5,
              f"centre {f0[0]:.5f} at t=0 -> {fr[lo, 0]:.5f} at t={ts[lo]:.1f} "
              f"(factor {f0[0] / fr[lo, 0]:.1f}) -> {fr[-1, 0]:.5f} at t={ts[-1]:.1f}; "
              f"hottest shell cell reaches {fr[:, 1:].max():.5f}",
              "drains > 5x then recovers"))

    A.append(("C8  PRINTED NOT GATED: what this cluster cannot say. It is ONE "
              "SHELL DEEP, so the return at C7 is reflection off its own "
              "boundary and not a property of the bulk; no speed is quoted. "
              "Mass model DECLARED: unit mass per triangle, lumped m/3 to each "
              "corner -- the uniform-lamina model moved a period by 7 percent "
              "elsewhere in this programme",
              True,
              f"9 cells, 1 shell, {len(interior)} interior vertices of {n}",
              "printed"))
    return A, d


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_cl -- a jitterbug cluster with an interior, and a local impulse")
    print("=" * 78)
    checks, d = gate()
    bad = 0
    for name, ok, got, want in checks:
        tag = "PASS" if ok else "FAIL"
        bad += 0 if ok else 1
        print(f"  {tag}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")
    print(f"\n  {'t':>7} {'centre':>9} {'shell':>9} {'hottest':>9}")
    for t, sh, E in d["hist"][::2]:
        f = sh / sh.sum()
        print(f"  {t:7.2f} {f[0]:9.5f} {f[1:].sum():9.5f} {f[1:].max():9.5f}")
    print()
    print("  ALL CHECKS PASSED." if not bad else f"  {bad} CHECK(S) FAILED.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
