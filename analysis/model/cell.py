"""cell -- a jitterbug chain with mass, integrated, and the congruence kept.

WHAT THIS IS FOR. Everything before this file that "moved" a chain drove it
through a single coherent angle: both end cells locked to theta, the middle to
theta-60, spacing from the shared-face law. That parameterisation cannot
represent a state in which one cell differs from another, so it reported one
degree of freedom and zero lag -- properties of the parameterisation, not of the
linkage. Measured here instead: the chain has SIX INTERNAL DOF PER CELL, growing
without bound (6, 12, 18, 24, 30, 36 for 1..6 cells).

CONGRUENCE IS NOT IDENTITY, which is the whole reason this builder exists.
Fuller, "Deceptiveness of Topology -- Quanta Lost by Congruence" (1977): the
jitterbug has 24 edges and 12 vertices at EVERY configuration; at the octahedron
they are "24 EDGES CONGRUENT AS 12, 12 VERTICES CONGRUENT AS 6", and at the
tetrahedron "congruent as 6" and "congruent as 4". His complaint is that
"Eulerian topological accounting as presently practiced ... accounts each of
these multicongruent topological aspects as consisting of only one of such
aspects."

That is exactly what deduplicating by POSITION does, and doing it costs two real
things: an octahedral cell loses half its struts and so half its inertia, and
its congruent vertices get welded into a rigid solid that can never come apart
again. So this builder carries all 24 struts and all 8 triangles per cell as
IDENTITY objects, welds only genuine shared-face vertices, and lets congruent
quanta coincide without merging. R1 checks the count against Fuller's table.

TWO CONSTRUCTION TRAPS, both of which bit before being caught:
  * The within-cell slot->vertex map must come from a GENERIC angle. Derived at
    gamma = 60 it would see 6 vertices instead of 12.
  * The shared-face corner correspondence FLIPS with the direction of the 60
    degree offset -- (30,-30) and (60,0) give one pairing, (0,60) the other.
    A chain alternates, so every second weld needs the other one. Using a single
    correspondence throughout produced malformed triangles whose struts were not
    sqrt(2), which is what R2 now gates.
"""
from __future__ import annotations

import itertools as it
import sys

import numpy as np

from analysis.model import plates as Z
EL = float(np.linalg.norm(Z.corners(0.0)[0][0] - Z.corners(0.0)[0][1]))
ZC = EL * np.sqrt(2.0 / 3.0)
NH = Z.corners(0.0)[0].mean(0)
NH = NH / np.linalg.norm(NH)
GEN = 30.0                     # generic angle: the topology is read off here
_XG = Z.corners(GEN)

SLOT, _P = {}, {}
for _f in range(8):
    for _c in range(3):
        SLOT[(_f, _c)] = _P.setdefault(tuple(np.round(_XG[_f][_c], 6)), len(_P))
NV = len(_P)
FP = max(range(8), key=lambda f: np.dot(_XG[f].mean(0) / np.linalg.norm(_XG[f].mean(0)), NH))
FM = min(range(8), key=lambda f: np.dot(_XG[f].mean(0) / np.linalg.norm(_XG[f].mean(0)), NH))


def cell_verts(g, o):
    X = Z.corners(g) + o
    V = np.zeros((NV, 3))
    seen = np.zeros(NV, bool)
    for f in range(8):
        for c in range(3):
            i = SLOT[(f, c)]
            if not seen[i]:
                V[i] = X[f][c]
                seen[i] = True
    return V


#: The two mating faces are fixed: cell k presents FP (+n), cell k+1 presents
#: FM (-n). Only the CORNER PERMUTATION between them varies with the direction
#: of the 60 degree offset.
_FM_CORNERS = [SLOT[(FM, c)] for c in range(3)]


def weld_for(ga, gb):
    """Corner correspondence across a shared face. DIRECTION-DEPENDENT.

    The search is restricted to the MATING FACE's three corners. Searching all
    twelve of the neighbour's vertices -- which is what this did first -- lets a
    nearer non-corner win: at the (0, 60) pair, where the second cell is the
    open VE and the first is closed, it returned {4, 5, 8}, which are not the
    corners of any face. The weld then joined three vertices that do not form a
    triangle, and nothing caught it, because the doubly-written shared vertices
    still AGREED IN POSITION to 3e-15. A redundancy check on positions does not
    validate identity. Geometry was the wrong tool for a bookkeeping problem --
    the same sentence jb_aa's node mapping already carries."""
    sep = ZC * (np.cos(np.radians(ga)) + np.cos(np.radians(gb)))
    A, B = cell_verts(ga, np.zeros(3)), cell_verts(gb, sep * NH)
    out = []
    for c in range(3):
        a = SLOT[(FP, c)]
        d = [np.linalg.norm(B[k] - A[a]) for k in _FM_CORNERS]
        out.append((a, _FM_CORNERS[int(np.argmin(d))]))
    tgt = {b for _, b in out}
    if tgt != set(_FM_CORNERS):
        raise AssertionError(f"weld for ({ga},{gb}) hit {sorted(tgt)}, "
                             f"not the mating face {sorted(_FM_CORNERS)}")
    return out


def build(gammas):
    """All 24 struts and 8 triangles per cell, by identity. Vertices welded only
    across genuine shared faces. Returns positions, struts, triangles, the
    owning cell of each triangle, the cell->vertex map, and the worst
    disagreement between the two writes of every shared vertex."""
    N = len(gammas)
    seps = [ZC * (np.cos(np.radians(gammas[k])) + np.cos(np.radians(gammas[k + 1])))
            for k in range(N - 1)]
    orig = [np.zeros(3)]
    for s in seps:
        orig.append(orig[-1] + s * NH)
    par = list(range(N * NV))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for k in range(N - 1):
        for a, b in weld_for(gammas[k], gammas[k + 1]):
            ra, rb = find(k * NV + a), find((k + 1) * NV + b)
            if ra != rb:
                par[ra] = rb
    uniq, gid = {}, np.zeros((N, NV), int)
    for k in range(N):
        for i in range(NV):
            gid[k, i] = uniq.setdefault(find(k * NV + i), len(uniq))
    P = np.zeros((len(uniq), 3))
    wr = np.zeros(len(uniq), bool)
    worst = 0.0
    for k in range(N):
        V = cell_verts(gammas[k], orig[k])
        for i in range(NV):
            g = gid[k, i]
            if wr[g]:
                worst = max(worst, float(np.linalg.norm(P[g] - V[i])))
            P[g] = V[i]
            wr[g] = True
    tris = [tuple(gid[k, SLOT[(f, c)]] for c in range(3)) for k in range(N) for f in range(8)]
    tcell = [k for k in range(N) for _ in range(8)]
    bars = [(t[a], t[b]) for t in tris for a, b in ((0, 1), (1, 2), (0, 2))]
    return P, bars, tris, tcell, gid, worst


def rigidity(P, bars):
    n = len(P)
    R = np.zeros((len(bars), 3 * n))
    for r, (i, j) in enumerate(bars):
        d = P[i] - P[j]
        d = d / np.linalg.norm(d)
        R[r, 3 * i:3 * i + 3] = d
        R[r, 3 * j:3 * j + 3] = -d
    s = np.linalg.svd(R, compute_uv=False)
    return int((s > s[0] * 1e-8).sum())


def integrate(P, bars, tris, gid, ncell, spin_cell, spin_face, h=0.005, nsteps=600):
    """RATTLE, V = 0. Redundant constraints are handled by least squares."""
    n, nb = len(P), len(bars)
    mass = np.zeros(n)
    for t in tris:
        for i in t:
            mass[i] += 1.0 / 3.0
    BI = np.array([b[0] for b in bars])
    BJ = np.array([b[1] for b in bars])
    L2 = np.sum((P[BI] - P[BJ]) ** 2, axis=1)
    Minv = np.repeat(1.0 / mass, 3)
    cellv = [sorted(set(gid[k])) for k in range(ncell)]

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

    tv = list(tris[spin_face])
    ctr = P[tv].mean(0)
    ax = ctr - P[cellv[spin_cell]].mean(0)
    ax = ax / np.linalg.norm(ax)
    v = np.zeros_like(P)
    for i in tv:
        v[i] = np.cross(ax, P[i] - ctr)      # SPIN about the triangle's own axis
    G0 = G(P)
    lam = np.linalg.lstsq(G0 @ (Minv[:, None] * G0.T), -(G0 @ v.reshape(-1)), rcond=None)[0]
    v = (v.reshape(-1) + Minv * (G0.T @ lam)).reshape(-1, 3)
    adm = float(np.linalg.norm(G(P) @ v.reshape(-1)))
    E0 = 0.5 * float(np.sum(mass[:, None] * v * v))

    def cellKE(vv):
        return np.array([0.5 * np.sum(mass[c, None] * vv[c] ** 2) for c in cellv])

    q = P.copy()
    hist = [(0.0, cellKE(v), E0)]
    for s in range(nsteps):
        q0 = q.copy()
        G0 = G(q0)
        qp = q + h * v
        lam = np.zeros(nb)
        for _ in range(20):
            qn = qp + (0.5 * h * h * (Minv * (G0.T @ lam))).reshape(-1, 3)
            g = gvec(qn)
            if np.abs(g).max() < 1e-12:
                break
            lam -= np.linalg.lstsq(G(qn) @ (0.5 * h * h * (Minv[:, None] * G0.T)), g, rcond=None)[0]
        qn = qp + (0.5 * h * h * (Minv * (G0.T @ lam))).reshape(-1, 3)
        vh = v + (0.5 * h * (Minv * (G0.T @ lam))).reshape(-1, 3)
        G1 = G(qn)
        mu = np.linalg.lstsq(G1 @ (Minv[:, None] * G1.T), -(G1 @ vh.reshape(-1)) / (0.5 * h),
                             rcond=None)[0]
        v = vh + (0.5 * h * (Minv * (G1.T @ mu))).reshape(-1, 3)
        q = qn
        if (s + 1) % 20 == 0:
            hist.append(((s + 1) * h, cellKE(v), 0.5 * float(np.sum(mass[:, None] * v * v))))
    return hist, adm, E0


def gate():
    checks = []
    A = checks.append

    # R1 -- Fuller's congruence table, from the geometry
    ident_e = 24
    rows = []
    for g, name in ((0.0, "VE"), (22.238756, "icosahedron"), (60.0, "octahedron")):
        X = Z.corners(g)
        vp = {tuple(np.round(X[f][c], 6)) for f in range(8) for c in range(3)}
        ep = {tuple(sorted((tuple(np.round(X[f][a], 6)), tuple(np.round(X[f][b], 6)))))
              for f in range(8) for a, b in ((0, 1), (1, 2), (0, 2))}
        rows.append((name, ident_e, len(ep), NV, len(vp)))
    A(("R1  FULLER'S CONGRUENCE TABLE, reproduced from the geometry rather than "
       "quoted. 'Deceptiveness of Topology -- Quanta Lost by Congruence' (1977): "
       "the VE has 24 edges and 12 vertices; the octahedron has 24 EDGES "
       "CONGRUENT AS 12 and 12 VERTICES CONGRUENT AS 6. The identity counts are "
       "24 and 12 at EVERY angle -- that invariance IS the fixed-incidence "
       "constraint, in Fuller's own terms",
       all(e_id == 24 and v_id == 12 for _, e_id, _, v_id, _ in rows)
       and rows[2][2] == 12 and rows[2][4] == 6
       and rows[0][2] == 24 and rows[0][4] == 12,
       "; ".join(f"{n}: {ei} edges congruent as {ep}, {vi} vertices congruent as {vp}"
                 for n, ei, ep, vi, vp in rows),
       "octahedron 24->12 edges, 12->6 vertices; identity always 24 and 12"))

    # R2 -- the build is a real framework at every configuration
    cfgs = {"octa-VE-octa": [60., 0., 60.], "VE-octa-VE": [0., 60., 0.],
            "generic": [30., -30., 30.], "6-cell": [60., 0., 60., 0., 60., 0.]}
    built = {k: build(v) for k, v in cfgs.items()}
    bad = []
    for k, (P, bars, tris, tcell, gid, worst) in built.items():
        L = np.array([np.linalg.norm(P[i] - P[j]) for i, j in bars])
        if abs(L.max() - EL) > 1e-9 or abs(L.min() - EL) > 1e-9 or worst > 1e-9:
            bad.append(k)
    A(("R2  EVERY STRUT IS A STRUT AND EVERY SHARED VERTEX AGREES. The corner "
       "correspondence across a shared face FLIPS with the direction of the 60 "
       "degree offset, and a chain alternates, so a single correspondence used "
       "throughout silently produces malformed triangles. The doubly-written "
       "shared vertices are the free check on it -- the same redundancy that "
       "corrupted jb_aa's initial condition proves this one right",
       not bad,
       "; ".join(f"{k}: struts {min(np.linalg.norm(built[k][0][i] - built[k][0][j]) for i, j in built[k][1]):.6f}"
                 f"..{max(np.linalg.norm(built[k][0][i] - built[k][0][j]) for i, j in built[k][1]):.6f}, "
                 f"write disagreement {built[k][5]:.1e}" for k in cfgs),
       f"all struts {EL:.6f}, disagreement < 1e-9"))

    # R2b -- the weld is FACE TO FACE, in every direction
    pairs = [(30., -30.), (-30., 30.), (60., 0.), (0., 60.), (45., -15.), (-15., 45.)]
    welds = {}
    okw = True
    for ga, gb in pairs:
        try:
            w = weld_for(ga, gb)
            welds[(ga, gb)] = sorted(b for _, b in w)
        except AssertionError:
            okw = False
            welds[(ga, gb)] = "NOT A FACE"
    A(("R2b THE WELD JOINS FACE TO FACE, in both offset directions. The corner "
       "PERMUTATION flips with the direction; the mating FACE does not. "
       "Searching all twelve of the neighbour's vertices for the nearest match "
       "-- rather than the mating face's three -- let a non-corner win at the "
       "(0, 60) pair and welded three vertices that form no triangle. Nothing "
       "caught it: the doubly-written shared vertices still agreed in POSITION "
       "to 3e-15, and a redundancy check on positions does not validate "
       "IDENTITY",
       okw and all(v == sorted(_FM_CORNERS) for v in welds.values()),
       f"mating face corners {sorted(_FM_CORNERS)}; targets "
       + ", ".join(f"{k}->{v}" for k, v in welds.items()),
       "every direction lands on the mating face"))

    # R3 -- all 24 struts and 8 triangles per cell are CARRIED, not deduplicated
    P, bars, tris, tcell, gid, worst = built["octa-VE-octa"]
    dup_b = len(bars) - len({(min(i, j), max(i, j)) for i, j in bars})
    A(("R3  THE CONGRUENT QUANTA ARE CARRIED. 24 struts and 8 triangles per "
       "cell survive into the framework even where they coincide. Deduplicating "
       "them by position is precisely the accounting Fuller calls losing the "
       "quanta, and it costs an octahedral cell HALF ITS INERTIA and welds its "
       "congruent vertices into a solid that can never come apart",
       len(bars) == 24 * 3 and len(tris) == 8 * 3 and dup_b > 0,
       f"3 cells -> {len(bars)} struts (24 x 3), {len(tris)} triangles (8 x 3), "
       f"of which {dup_b} struts are congruent-but-distinct",
       "24 struts and 8 triangles per cell, congruent ones kept"))

    # R4 -- internal DOF grows with length; the chain is NOT one degree of freedom
    dofs = [3 * len(build([30. if k % 2 == 0 else -30. for k in range(N)])[0])
            - rigidity(*build([30. if k % 2 == 0 else -30. for k in range(N)])[:2]) - 6
            for N in range(1, 7)]
    A(("R4  THE CHAIN HAS SIX INTERNAL DOF PER CELL, GROWING WITHOUT BOUND. "
       "Every earlier chain in this programme was driven through a single "
       "coherent angle and so reported ONE degree of freedom and zero lag; "
       "those were properties of the parameterisation. A disturbance needs "
       "somewhere to be, and this is where",
       dofs == [6, 12, 18, 24, 30, 36],
       f"internal DOF for 1..6 cells: {dofs}", "[6, 12, 18, 24, 30, 36]"))

    # R5 -- closing a cell to an octahedron costs the chain freedom
    d_oct = 3 * len(built["octa-VE-octa"][0]) - rigidity(*built["octa-VE-octa"][:2]) - 6
    d_gen = 3 * len(built["generic"][0]) - rigidity(*built["generic"][:2]) - 6
    A(("R5  CLOSING THE END CELLS COSTS NO FREEDOM. octa-VE-octa and the "
       "generic chain both carry 18 internal DOF -- 6 per cell either way. The "
       "congruence at the octahedron changes WHICH quanta coincide, not how "
       "many freedoms there are. An earlier version of this file reported 15 "
       "against 18 and read it as the octahedral configuration being "
       "structurally different; that was the weld of R2b searching all twelve "
       "of the neighbour's vertices and joining three that form no triangle. "
       "The number was an artifact of a bug, not a property of the medium",
       d_oct == d_gen == 18,
       f"octa-VE-octa {d_oct} internal DOF, generic {d_gen}", "both 18"))

    # R6 -- THE MEASUREMENT: impulse one triangle, integrate, energy audited
    P, bars, tris, tcell, gid, worst = built["octa-VE-octa"]
    blue = [t for t, k in enumerate(tcell) if k == 2][1]
    hist, adm, E0 = integrate(P, bars, tris, gid, 3, 2, blue)
    drift = max(abs(E - E0) / E0 for _, _, E in hist)
    A(("R6  THE IMPULSE IS ADMISSIBLE AND THE INTEGRATION CONSERVES ENERGY. A "
       "spin on one triangle is not a legal initial condition until it is "
       "projected onto the constraint tangent space; V = 0, so energy is the "
       "only audit there is and RATTLE has to earn it",
       adm < 1e-12 and drift < 1e-4,
       f"||G v|| after projection {adm:.1e}; E = {E0:.6f} conserved to "
       f"{drift:.1e} over {len(hist)} samples", "admissible, drift < 1e-4"))

    # R7 -- the onset is INSTANTANEOUS; only the transfer takes time
    ke0 = hist[0][1] / hist[0][1].sum()
    kel = hist[-1][1] / hist[-1][1].sum()
    A(("R7  THE ONSET IS INSTANTANEOUS. Every cell has kinetic energy at "
       "t = 0: a rigid constraint has infinite signal speed, so the projection "
       "that makes the impulse admissible reaches the whole chain at once and "
       "there is no onset lag to measure anywhere. This medium therefore has "
       "NO WAVEFRONT; any apparent front is mode superposition. CORRECTED "
       "2026-08-28: this row used to add 'and a finite signal speed would need "
       "compliant constraints -- the fork this programme has already "
       "rejected', which is FALSE and was the sentence that kept the wave "
       "programme looking in the wrong place. It needs CLEARANCE, which is "
       "compatible with perfectly rigid struts, and jb_ct measures it. The "
       "measurement in this row is untouched: it is the PLAY-FREE LIMIT of "
       "jb_ct's finite speed, which diverges as the joint tightens. "
       "TWO-SIDED: a lag would show as a zero here, and it never does",
       ke0[0] > 1e-6 and kel[0] > ke0[0],
       f"far cell's share of the kinetic energy: {ke0[0]:.2e} at t=0 -> "
       f"{kel[0]:.2e} at t={hist[-1][0]:.1f}, a factor {kel[0] / ke0[0]:.2f}",
       "nonzero at t=0 (no lag), and rising"))

    A(("R7b PRINTED NOT GATED: the SIZE of the transfer, which is small over "
       "this window and is NOT quoted as a rate. An earlier version of this "
       "file reported a factor of 11.8 here; that came from the R2b weld bug "
       "and is withdrawn. Three cells with free ends is also the wrong object "
       "to read transport from -- every one of its 30 vertices lies on an "
       "outward-facing face, so it has no interior at all, and embedding it in "
       "the array by pinning what the array holds leaves ZERO degrees of "
       "freedom. A transport number needs a patch with cells to spare",
       True,
       f"far cell {ke0[0]:.2e} -> {kel[0]:.2e} over t={hist[-1][0]:.1f} "
       f"(factor {kel[0] / ke0[0]:.2f})", "printed"))

    return checks, hist


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_ic -- a jitterbug chain with mass, and the congruence kept")
    print("=" * 78)
    checks, hist = gate()
    bad = 0
    for name, ok, got, want in checks:
        tag = "PASS" if ok else "FAIL"
        bad += 0 if ok else 1
        print(f"  {tag}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")
    print("\n  per-cell kinetic energy after a spin on one triangle of cell 2:")
    print(f"  {'t':>7} {'cell0':>10} {'cell1':>10} {'cell2':>10}")
    for t, ke, E in hist[::max(1, len(hist) // 12)]:
        f = ke / ke.sum()
        print(f"  {t:7.2f} {f[0]:10.5f} {f[1]:10.5f} {f[2]:10.5f}")
    print()
    print("  WHAT THIS DOES AND DOES NOT LICENSE.")
    print("   * The chain transports: a spin on one triangle redistributes")
    print("     energy toward the far end over time, and a longer chain sends")
    print("     it out, re-concentrates it at the far wall and returns it.")
    print("   * It licenses NO signal speed. The onset is simultaneous")
    print("     everywhere (R7), so there is no wavefront to time.")
    print("   * Mass model DECLARED: unit mass per triangle, lumped m/3 to each")
    print("     corner. The uniform-lamina model moved a period by 7% elsewhere")
    print("     in this programme, so no number here is model-independent.")
    print()
    print("  ALL CHECKS PASSED." if not bad else f"  {bad} CHECK(S) FAILED.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
