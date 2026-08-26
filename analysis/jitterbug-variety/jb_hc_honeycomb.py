"""jb_hc -- the VE/octa rectified cubic honeycomb, and the two waves it carries.

THE PACKING THIS EPIC IS ACTUALLY ABOUT. `jb_x_array_linkage.build_topologies`
places neighbours at 2*v[g], |offset| = 2*sqrt(2), where two cuboctahedra touch
at ONE VERTEX. That is not how these cells pack. The real structure is the
rectified cubic honeycomb -- cuboctahedra and octahedra alternating, space
filling -- and neighbours share FACES:

    VE  <-> VE    4 shared vertices (a square face)      at distance 2
    VE  <-> OCTA  3 shared vertices (a triangular face)  at distance sqrt(3)
    the codebase's own topologies:  1 shared vertex

6 square + 8 triangular = 14 neighbours, which is every face of the VE, nothing
left over. VE cells sit at even lattice points and octa cells at all-odd, and
the octahedron IS the collapsed VE: `verts(60)` has exactly 6 distinct vertices
at radius 1, edge sqrt(2), the same edge as the VE's struts.

Consequence, recorded because it retires a lot: SC7, CHAIN5, SQUARE4, CUBE8-M/R,
CUBE27-M and FCC13 are all single-vertex-contact structures, so every number
computed on them is superseded -- a "30 degree wall", a 6.7 degree dephasing, a
0.94-decades-per-unit localisation, and the whole qvf.19 lock surface. It also
retires two diagnosis rounds: the array "came apart" not because a ball joint was
missing and not only because of the wire defect (inviscid-l1d, still real and
still unfixed), but because the cells were never touching the way they touch.

WHY THIS FILE EXISTS AT ALL. Every result below was first obtained in throwaway
one-liners and written up in T2 with no code behind it. The consolidation pass of
2026-08-25 found the epic's entire current model unreproducible: results in
memory, three HTML figures in a session-scoped temp directory, and nothing in the
repo. This file is that model, with every headline number as a gate row that can
fail. Canonical prose state: T2 `inviscid/qvf-epic-consolidated-state`.

WHAT IT ESTABLISHES

  THE EXCHANGE (H2, H3). Solving the shared-face constraint gives ONE motion:
  the hole cell's phase runs exactly 60 degrees ahead, no cell rotates, and the
  lattice is the fold half-diagonal. Driven DOWN from a = 0 to a = -60 the cells
  trade places -- VE becomes octahedron, octahedron becomes VE, lattice back to
  where it started -- with the cells separated at every step. Driven UP they
  interpenetrate immediately, which is why the direction matters.

  TWO WAVES (H5-H8). As a bar-and-joint framework (nodes at the shared vertices,
  bars along the rigid triangle edges) the medium carries two orthogonal soft
  families:
    POSITIONAL   anisotropic, c[100] = c[110] = sqrt(8/27), c[111] = sqrt(4/27),
                 with an exactly FLAT ZERO BAND along the six <110> closest-
                 packing directions and nothing soft anywhere else.
    PHASE        isotropic, c = 2/3 exactly, gapless, and the Goldstone mode of
                 the exchange itself.
  They do not mix: at M the positional band is 0 while both phase bands sit at
  2/sqrt(3).

  SCOPE OF THE PHASE RESULT, added after jb_bt_band_touching.py: everything in
  H7 and H8 is computed at a = -30, and the gaplessness is a property of THAT
  PHASE. At any other reference phase the same projection gaps the phase band
  by exactly sqrt(6)*|d lambda/d a|, which vanishes at a = -30 because lambda
  is stationary there. The fixed-lattice assumption is free at the midpoint and
  nowhere else. c = 2/3 is real, and it is the midpoint's number.

METHOD NOTE THAT COST THE MOST TO LEARN. The corner correspondence across a
shared face MUST be resolved by permutation search. Pairing corner c to corner c
leaves a seed residual of 1.0 at a = 0 and the solver then converges happily to a
spurious closure at b = 110.237 whose hole cell is not an octahedron at all. Two
published results died that way. The rule the gate enforces: ASSERT THE SEED IS
MACHINE ZERO BEFORE CONTINUING ANY BRANCH. A residual that closes to 1e-16 tells
you the constraints are satisfied, never that they are the right constraints --
and a system that closes for EVERY parameter value is underdetermined, which is
how "b = 60 + a/2" survived long enough to be written down.
"""

from __future__ import annotations

import itertools as it
import sys

import numpy as np

import jb_gp_plate_geometry as Z
from jb_x_array_linkage import STRUT_LEN, verts

#: Reference phase for the linearised models. Midpoint of the exchange, where
#: the two cells are congruent and the lattice is widest.
A_REF = -30.0

#: Phase offset of the hole sublattice. THE central relation; H2 measures it
#: rather than assuming it.
PHASE_OFFSET = 60.0

#: The eight cube-diagonal directions -- a VE's triangular-face neighbours.
DIRS = [np.array(t, dtype=float) for t in it.product((1, -1), repeat=3)]

TOL_EXACT = 1e-12      # machine-zero geometric coincidence
TOL_SEED = 1e-12       # the seed assertion above
TOL_MODE = 1e-8        # relative singular-value cut for rank/nullity


def face_toward(d):
    """The cell's triangular face whose outward plate normal points along `d`.
    Plate normals are fixed along the cube diagonals for every phase, so this
    is phase independent."""
    u = np.asarray(d, dtype=float)
    u = u / np.linalg.norm(u)
    return max(range(8), key=lambda k: float(Z.plate_normal(k) @ u))


def lattice(a):
    """The honeycomb's lattice parameter at phase `a`: the FOLD HALF-DIAGONAL.

    Measured, not fitted -- H2 checks it against this project's own
    `FOLD_TABLE_TARGET` values."""
    v = verts(a)
    return 0.5 * float(np.linalg.norm(v[0] - v[3]))


def cell(phase, origin):
    """One cell's 8 rigid triangles, placed. Cells never rotate in this
    honeycomb (H2 gates that), so a placement is a translation."""
    return Z.corners(phase) + np.asarray(origin, dtype=float)


def shared_vertices(pa, oa, pb, ob, tol=1e-8):
    """How many vertices two placed cells hold in common."""
    A = {tuple(np.round(p, 8)) for p in cell(pa, oa).reshape(-1, 3)}
    B = {tuple(np.round(p, 8)) for p in cell(pb, ob).reshape(-1, 3)}
    return len(A & B)


def face_pairing(a=0.0):
    """For each of the 8 directions: (hole face, VE face, corner permutation).

    The permutation is the load-bearing part -- see the module docstring."""
    b = a + PHASE_OFFSET
    L = lattice(a)
    out = []
    worst = 0.0
    for d in DIRS:
        fv = face_toward(d)          # VE's face toward the hole at +L*d
        fo = face_toward(-d)         # hole's face back at the VE
        V = Z.corners(a)[fv]
        O = Z.corners(b)[fo] + L * d
        perm = min(it.permutations(range(3)),
                   key=lambda pm: sum(np.linalg.norm(V[pm[c]] - O[c]) for c in range(3)))
        worst = max(worst, max(np.linalg.norm(V[perm[c]] - O[c]) for c in range(3)))
        out.append((d, fo, fv, perm))
    return out, worst


def separation(a):
    """Overlap of a VE cell and its hole neighbour, measured on the shared
    face's plane. Positive means the cells have passed through each other.

    A separating-plane test rather than `Z._is_piercing`: that predicate
    false-positives on edge-adjacent triangles, which every structural contact
    in this honeycomb is."""
    L = lattice(a)
    n = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    V = cell(a, np.zeros(3)).reshape(-1, 3)
    O = cell(a + PHASE_OFFSET, L * np.ones(3)).reshape(-1, 3)
    return float((V @ n).max() - (O @ n).min())


def unit_cell(a):
    """The periodic bar-and-joint framework: nodes at the shared vertices, bars
    along the rigid triangle edges, reduced modulo the cubic lattice.

    This is the formulation that made the mobility question answerable. Cells as
    rigid bodies need a 48-DOF-per-cell Jacobian whose kernel is easy to get
    wrong; a bar-joint rigidity matrix has no such ambiguity, and a triangle is
    exactly rigid under its three edges so nothing is lost."""
    L = lattice(a)
    A = 2.0 * L
    nodes, reps, slots = {}, [], []
    bars = set()

    def reduce(p):
        n = np.floor(p / A + 1e-9)
        q = p - A * n
        key = tuple(np.round(q, 6))
        if key not in nodes:
            nodes[key] = len(reps)
            reps.append(q)
        return nodes[key], n.astype(int)

    for ci, (ph, off) in enumerate(((a, np.zeros(3)),
                                    (a + PHASE_OFFSET, L * np.ones(3)))):
        X = cell(ph, off)
        for f in range(8):
            ids, offs = [], []
            for c in range(3):
                i, nn = reduce(X[f][c])
                ids.append(i)
                offs.append(nn)
                slots.append((ci, ph, off, f, c, i, tuple(nn)))
            for u, v in it.combinations(range(3), 2):
                i, j = ids[u], ids[v]
                R = tuple(offs[v] - offs[u])
                mR = tuple(-np.array(R))
                bars.add((i, j, R) if (i, R) <= (j, mR) else (j, i, mR))
    return np.array(reps), sorted(bars), A, slots


def bloch(P, bars, A, kvec, unit=True):
    """Rigidity matrix at wavevector `kvec`. `unit=True` normalises the bar
    vectors, which is what the DYNAMICS needs (spring force acts along the
    bar); `unit=False` is equivalent for rank and nullity."""
    n = len(P)
    M = np.zeros((len(bars), 3 * n), dtype=complex)
    for r, (i, j, Rv) in enumerate(bars):
        Rw = np.array(Rv, dtype=float) * A
        d = P[i] - (P[j] + Rw)
        if unit:
            d = d / np.linalg.norm(d)
        M[r, 3 * i:3 * i + 3] += d
        M[r, 3 * j:3 * j + 3] -= d * np.exp(1j * float(np.dot(kvec, Rw)))
    return M


def zero_modes(P, bars, A, kvec, rtol=TOL_MODE):
    s = np.linalg.svd(bloch(P, bars, A, kvec, unit=False), compute_uv=False)
    return 3 * len(P) - int((s > s[0] * rtol).sum()), s


def phase_basis(P, slots, A, kvec, h=1e-5):
    """The two phase directions in node space: d(node position)/d(cell phase),
    per sublattice, with Bloch phases. Orthonormalised.

    The cell phase is a DERIVED quantity of the node positions, not an extra
    degree of freedom -- so the phase field is a SUBSPACE of the node dynamics,
    and this is how you look at it."""
    n = len(P)
    B = np.zeros((3 * n, 2), dtype=complex)
    seen = {}
    for (ci, ph, off, f, c, i, nn) in slots:
        key = (ci, i, nn)
        if key in seen:
            continue
        dp = ((Z.corners(ph + h) + off)[f][c] - (Z.corners(ph - h) + off)[f][c]) / (2 * h)
        seen[key] = True
        R = np.array(nn, dtype=float) * A
        B[3 * i:3 * i + 3, ci] += dp * np.exp(-1j * float(np.dot(kvec, R)))
    q, _ = np.linalg.qr(B)
    return q


def speed(bands_at, direction, G, eps=(0.02, 0.01, 0.005)):
    """Branch slope near Gamma, refined; returns the finest estimate."""
    u = np.asarray(direction, dtype=float)
    u = u / np.linalg.norm(u)
    out = [bands_at(e * G * u) / (e * G) for e in eps]
    return out[-1], out


# ==========================================================================
# H1: THE PACKING. Neighbours share FACES, not a vertex.
# ==========================================================================

def h1_packing():
    v0 = verts(0.0)
    octa = np.array([[1., 0, 0], [-1, 0, 0], [0, 1, 0],
                     [0, -1, 0], [0, 0, 1], [0, 0, -1]])
    A = {tuple(np.round(p, 8)) for p in v0}
    sq = len(A & {tuple(np.round(p, 8)) for p in v0 + np.array([2., 0, 0])})
    tri = len(A & {tuple(np.round(p, 8)) for p in octa + np.array([1., 1, 1])})
    # what the codebase builds, for contrast
    from jb_x_array_linkage import ANTI
    off = v0[0] - v0[ANTI[0]]
    old = len(A & {tuple(np.round(p, 8)) for p in v0 + off})
    n_sq = sum(1 for d in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
               if len(A & {tuple(np.round(p, 8)) for p in v0 + 2 * np.array(d, dtype=float)}) == 4)
    n_tri = sum(1 for d in DIRS
                if len(A & {tuple(np.round(p, 8)) for p in octa + d}) == 3)
    u60 = np.unique(np.round(verts(60.0), 6), axis=0)
    r60 = np.linalg.norm(u60, axis=1)
    return dict(sq=sq, tri=tri, old=old, n_sq=n_sq, n_tri=n_tri,
                n60=len(u60), r60max=float(r60.max()), r60min=float(r60.min()))


# ==========================================================================
# H2: THE EXCHANGE. One motion; b = a + 60; lambda = fold half-diagonal.
# ==========================================================================

def _rod(w):
    th = float(np.linalg.norm(w))
    if th < 1e-14:
        return np.eye(3)
    k = w / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


def _exchange_residual(p, a, pairs):
    """8 unknowns (VE rotation, hole rotation, hole phase, lambda) against
    8 faces x 3 vertices x 3 coords = 72 constraints. PROPERLY over-determined
    -- the earlier version that let each of 8 VEs move independently was not,
    and closed for every (a, b), which is how a false relation got published."""
    Rv, Rh, b, L = _rod(p[0:3]), _rod(p[3:6]), p[6], p[7]
    xv, xo = Z.corners(a), Z.corners(b)
    out = []
    for (d, fo, fv, perm) in pairs:
        for c in range(3):
            out.append(Rv @ xv[fv][perm[c]] - (Rh @ xo[fo][c] + L * d))
    return np.array(out).ravel()


def h2_exchange(angles=(0.0, -5.0, -15.0, -30.0, -45.0, -60.0)):
    from scipy.optimize import least_squares
    pairs, seed_err = face_pairing(0.0)
    p = np.array([0, 0, 0, 0, 0, 0, PHASE_OFFSET, 1.0], dtype=float)
    seed_res = float(np.abs(_exchange_residual(p, 0.0, pairs)).max())
    rows = []
    for a in angles:
        s = least_squares(_exchange_residual, p, args=(float(a), pairs),
                          xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=8000)
        m = float(np.abs(s.fun).max())
        if m < 1e-9:
            p = s.x.copy()
        rows.append(dict(a=float(a), resid=m, b=float(s.x[6]), lam=float(s.x[7]),
                         twist_v=float(np.degrees(np.linalg.norm(s.x[0:3]))),
                         twist_h=float(np.degrees(np.linalg.norm(s.x[3:6]))),
                         sep=separation(float(a)),
                         r_ve=float(np.linalg.norm(verts(float(a))[0])),
                         r_hole=float(np.linalg.norm(verts(float(a) + PHASE_OFFSET)[0]))))
    # uniqueness: is any OTHER hole phase admissible?
    def resid7(q, a, b):
        Rv, Rh, L = _rod(q[0:3]), _rod(q[3:6]), q[6]
        xv, xo = Z.corners(a), Z.corners(b)
        o = []
        for (d, fo, fv, perm) in pairs:
            for c in range(3):
                o.append(Rv @ xv[fv][perm[c]] - (Rh @ xo[fo][c] + L * d))
        return np.array(o).ravel()
    admissible = []
    for b in np.arange(-60.0, 181.0, 2.5):
        best = 9.0
        for seed in (np.array([0, 0, 0, 0, 0, 0, 1.0]), np.array([0, 0, 0, 0, 0, 0, -1.0])):
            s = least_squares(resid7, seed, args=(-15.0, float(b)),
                              xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=800)
            best = min(best, float(np.abs(s.fun).max()))
        if best < 1e-9:
            admissible.append(round(float(b), 1))
    # lambda against this project's own fold table
    # `lattice` IS the fold half-diagonal d1, checked against this project's
    # own FOLD_TABLE_TARGET. Quoted at POSITIVE a, where the table lives; the
    # exchange runs at negative a, where the same function rises above 1. The
    # identity is the claim, not the sign.
    fold = {5.0: 0.94588, 10.0: 0.88455, 22.238756093: 0.70711,
            30.0: 0.57735, 45.0: 0.29886}
    fold_err = max(abs(lattice(k) - v) for k, v in fold.items())
    return dict(rows=rows, seed_err=seed_err, seed_res=seed_res,
                admissible=admissible, fold_err=fold_err)


# ==========================================================================
# H3: EIGHT AROUND ONE. The reciprocal condition, and the square contact.
# ==========================================================================

def h3_cluster(angles=(0.0, -10.0, -30.0, -45.0, -60.0)):
    rows = []
    for a in angles:
        L = lattice(a)
        b = a + PHASE_OFFSET
        tri = [shared_vertices(b, np.zeros(3), a, L * d) for d in DIRS]
        sq = []
        for i, j in it.combinations(range(8), 2):
            if int(np.abs(DIRS[i] - DIRS[j]).sum()) == 2:
                sq.append(shared_vertices(a, L * DIRS[i], a, L * DIRS[j]))
        rows.append(dict(a=float(a), tri=sorted(set(tri)), sq=sorted(set(sq))))
    return rows


# ==========================================================================
# H4-H8: THE FRAMEWORK, AND THE TWO WAVES.
# ==========================================================================

def h4_framework(a=A_REF):
    P, bars, A, slots = unit_cell(a)
    deg = {}
    for (i, j, _R) in bars:
        deg[i] = deg.get(i, 0) + 1
        deg[j] = deg.get(j, 0) + 1
    return dict(n=len(P), nb=len(bars), dof=3 * len(P),
                degrees=sorted(set(deg.values())), P=P, bars=bars, A=A, slots=slots)


def h5_zero_modes(fw):
    P, bars, A = fw["P"], fw["bars"], fw["A"]
    G = np.pi / A
    out = {}
    for lbl, d in [("<110>", [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]),
                   ("<100>", [(1, 0, 0), (0, 1, 0), (0, 0, 1)]),
                   ("<111>", [(1, 1, 1), (1, 1, -1)])]:
        vals = [float(zero_modes(P, bars, A, 0.6 * G * np.array(x, dtype=float))[1][-1])
                for x in d]
        out[lbl] = (max(vals), min(vals))
    gam = zero_modes(P, bars, A, np.zeros(3))[0]
    gen = zero_modes(P, bars, A, np.array([0.41, 0.67, 0.29]) * G)[0]
    line = [zero_modes(P, bars, A, np.array([t, t, 0.0]) * G)[0]
            for t in np.linspace(0.1, 1.0, 10)]
    return dict(bydir=out, gamma=gam, generic=gen, line=line)


def h6_elastic(fw):
    P, bars, A = fw["P"], fw["bars"], fw["A"]
    G = np.pi / A

    def lowest(kv):
        s = np.linalg.svd(bloch(P, bars, A, kv), compute_uv=False)
        nz = [x for x in np.sort(s) if x > 1e-9]
        return nz[0] if nz else 0.0
    return {lbl: speed(lowest, d, G)[0]
            for lbl, d in (("[100]", (1, 0, 0)), ("[110]", (1, 1, 0)), ("[111]", (1, 1, 1)))}


def h7_phase_is_a_zero_mode(fw):
    P, bars, A, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
    n = len(P)
    M0 = bloch(P, bars, A, np.zeros(3), unit=False).real
    U, S, Vt = np.linalg.svd(M0)
    ns = Vt[np.sum(S > S[0] * TOL_MODE):]
    T = np.zeros((3, 3 * n))
    for d in range(3):
        for i in range(n):
            T[d, 3 * i + d] = 1.0
    Tq, _ = np.linalg.qr(T.T)
    proj = ns - (ns @ Tq) @ Tq.T
    q, rr = np.linalg.qr(proj.T)
    keep = [i for i in range(q.shape[1]) if abs(rr[i, i]) > 1e-8]
    dph = np.zeros(3 * n)
    seen = set()
    for (ci, ph, off, f, c, i, nn) in slots:
        if (ci, i, nn) in seen:
            continue
        seen.add((ci, i, nn))
        h = 1e-5
        dph[3 * i:3 * i + 3] += ((Z.corners(ph + h) + off)[f][c]
                                 - (Z.corners(ph - h) + off)[f][c]) / (2 * h)
    dph /= np.linalg.norm(dph)
    dil = np.zeros(3 * n)
    for i in range(n):
        dil[3 * i:3 * i + 3] = P[i] - P.mean(axis=0)
    dil /= np.linalg.norm(dil)
    overlap = abs(float(q[:, keep[0]] @ dph)) if keep else 0.0
    return dict(nullity=ns.shape[0], nontrans=len(keep), overlap=overlap,
                cost_phase=float(np.linalg.norm(M0 @ dph)),
                cost_dilation=float(np.linalg.norm(M0 @ dil)))


def h8_phase_wave(fw):
    P, bars, A, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
    G = np.pi / A

    def bands(kv):
        q = phase_basis(P, slots, A, kv)
        M = bloch(P, bars, A, kv)
        D = M.conj().T @ M
        return np.sqrt(np.maximum(np.linalg.eigvalsh(q.conj().T @ D @ q), 0))
    speeds = {lbl: speed(lambda kv: bands(kv)[0], d, G)[0]
              for lbl, d in (("[100]", (1, 0, 0)), ("[110]", (1, 1, 0)), ("[111]", (1, 1, 1)))}
    gap = float(bands(np.zeros(3))[1])
    at_M = bands(np.array([1.0, 1.0, 0.0]) * G)
    pos_at_M = zero_modes(P, bars, A, np.array([1.0, 1.0, 0.0]) * G)[0]
    return dict(speeds=speeds, gap=gap, at_M=[float(x) for x in at_M], pos_at_M=pos_at_M)


# ==========================================================================
# THE GATE
# ==========================================================================

def gate(h1, h2, h3, fw, h5, h6, h7, h8):
    checks = []
    R = checks.append

    R(("H1  VE-VE neighbours share a SQUARE FACE (4 vertices)",
       h1["sq"] == 4, f"{h1['sq']} vertices", "4"))
    R(("H1  VE-OCTA neighbours share a TRIANGULAR FACE (3 vertices)",
       h1["tri"] == 3, f"{h1['tri']} vertices", "3"))
    R(("H1  CONTROL: the codebase's own offset shares only ONE vertex "
       "(this is the defect)", h1["old"] == 1, f"{h1['old']} vertex", "1"))
    R(("H1  neighbour census closes: 6 square + 8 triangular = every VE face",
       h1["n_sq"] == 6 and h1["n_tri"] == 8,
       f"{h1['n_sq']} square + {h1['n_tri']} triangular", "6 + 8 = 14"))
    R(("H1  the octahedron IS the collapsed VE: verts(60) has 6 distinct "
       "vertices at radius 1",
       h1["n60"] == 6 and abs(h1["r60max"] - 1) < TOL_EXACT
       and abs(h1["r60min"] - 1) < TOL_EXACT,
       f"{h1['n60']} vertices, radius {h1['r60max']:.9f}", "6 at 1.0"))

    R(("H2  SEED IS EXACT before any continuation (the rule two dead results "
       "broke)", h2["seed_err"] < TOL_SEED and h2["seed_res"] < TOL_SEED,
       f"pairing {h2['seed_err']:.2e}, residual {h2['seed_res']:.2e}",
       f"< {TOL_SEED:.0e}"))
    rows = h2["rows"]
    R(("H2  the face constraint closes at EVERY phase of the exchange",
       len(rows) > 0 and all(r["resid"] < 1e-9 for r in rows),
       f"worst {max(r['resid'] for r in rows):.2e} over {len(rows)} phases", "< 1e-9"))
    R(("H2  b = a + 60 EXACTLY",
       all(abs(r["b"] - (r["a"] + PHASE_OFFSET)) < 1e-6 for r in rows),
       f"max dev {max(abs(r['b'] - (r['a'] + PHASE_OFFSET)) for r in rows):.2e}", "< 1e-6"))
    R(("H2  NO OTHER hole phase admits closure -- one motion, not a family "
       "-- CAN FAIL",
       h2["admissible"] == [45.0], f"admissible b = {h2['admissible']}", "[45.0] only"))
    R(("H2  lambda IS the fold half-diagonal (this project's own fold table)",
       h2["fold_err"] < 1e-5, f"max dev {h2['fold_err']:.2e}", "< 1e-5"))
    R(("H2  NO CELL ROTATES anywhere in the exchange",
       all(max(r["twist_v"], r["twist_h"]) < 1e-6 for r in rows),
       f"max twist {max(max(r['twist_v'], r['twist_h']) for r in rows):.2e} deg", "< 1e-6"))
    R(("H2  the cells stay SEPARATED at every step (a valid packing throughout) "
       "-- CAN FAIL",
       all(r["sep"] < 1e-9 for r in rows),
       f"worst overlap {max(r['sep'] for r in rows):.2e}", "< 1e-9"))
    first, last = rows[0], rows[-1]
    R(("H2  the ROLES SWAP: VE radius sqrt(2)->1 and hole radius 1->sqrt(2)",
       abs(first["r_ve"] - np.sqrt(2)) < 1e-6 and abs(first["r_hole"] - 1) < 1e-6
       and abs(last["r_ve"] - 1) < 1e-6 and abs(last["r_hole"] - np.sqrt(2)) < 1e-6,
       f"{first['r_ve']:.5f}/{first['r_hole']:.5f} -> "
       f"{last['r_ve']:.5f}/{last['r_hole']:.5f}", "swap"))
    R(("H2  the LATTICE RETURNS: lambda 1 -> 1 across the full exchange",
       abs(first["lam"] - 1) < 1e-6 and abs(last["lam"] - 1) < 1e-6,
       f"{first['lam']:.6f} -> {last['lam']:.6f}", "1 -> 1"))
    mid = [r for r in rows if abs(r["a"] - A_REF) < 1e-9]
    R(("H2  midpoint a=-30: cells CONGRUENT and lambda peaks at 2/sqrt(3)",
       len(mid) == 1 and abs(mid[0]["r_ve"] - mid[0]["r_hole"]) < 1e-6
       and abs(mid[0]["lam"] - 2 / np.sqrt(3)) < 1e-6,
       f"radii {mid[0]['r_ve']:.5f}/{mid[0]['r_hole']:.5f}, lambda {mid[0]['lam']:.6f}"
       if mid else "n/a", "equal, 1.154701"))

    R(("H3  RECIPROCAL: all 8 hole-VE contacts hold 3 vertices at every phase",
       all(r["tri"] == [3] for r in h3), f"{sorted({t for r in h3 for t in r['tri']})}", "[3]"))
    R(("H3  the SQUARE contact decays 4 -> 2 -> 1 as the folding squares open "
       "-- CAN FAIL",
       h3[0]["sq"] == [4] and h3[-1]["sq"] == [1] and h3[len(h3) // 2]["sq"] == [2],
       " -> ".join(str(r["sq"]) for r in h3), "[4] .. [2] .. [1]"))

    R(("H4  unit cell: 6 nodes, 24 bars, uniform node degree 8",
       fw["n"] == 6 and fw["nb"] == 24 and fw["degrees"] == [8],
       f"{fw['n']} nodes, {fw['nb']} bars, degree {fw['degrees']}", "6, 24, [8]"))
    R(("H4  Maxwell count is OVER-CONSTRAINED, so counting cannot decide "
       "mobility here", fw["dof"] - fw["nb"] < 0,
       f"{fw['dof']} - {fw['nb']} = {fw['dof'] - fw['nb']}", "< 0"))

    R(("H5  EXACT zero modes along ALL SIX <110> directions",
       h5["bydir"]["<110>"][0] < 1e-12,
       f"worst {h5['bydir']['<110>'][0]:.2e}", "< 1e-12"))
    R(("H5  and NONE along <100> or <111> -- CAN FAIL",
       h5["bydir"]["<100>"][1] > 0.1 and h5["bydir"]["<111>"][1] > 0.1,
       f"<100> {h5['bydir']['<100>'][1]:.3f}, <111> {h5['bydir']['<111>'][1]:.3f}", "> 0.1"))
    R(("H5  the <110> zero mode is a LINE, not a point (10 samples)",
       all(v >= 1 for v in h5["line"]), f"counts {h5['line']}", "all >= 1"))
    R(("H5  generic k has NO zero mode: the bulk is rigid off those lines",
       h5["generic"] == 0, f"{h5['generic']}", "0"))
    R(("H5  Gamma carries 4: three translations plus one more",
       h5["gamma"] == 4, f"{h5['gamma']}", "4"))

    R(("H6  ELASTIC speeds: c[100] = c[110] = sqrt(8/27)",
       abs(h6["[100]"] - np.sqrt(8 / 27)) < 1e-4 and abs(h6["[110]"] - np.sqrt(8 / 27)) < 1e-4,
       f"{h6['[100]']:.5f}, {h6['[110]']:.5f}", f"{np.sqrt(8/27):.5f}"))
    R(("H6  ELASTIC c[111] = sqrt(4/27), slower by exactly sqrt(2)",
       abs(h6["[111]"] - np.sqrt(4 / 27)) < 1e-4
       and abs(h6["[100]"] / h6["[111]"] - np.sqrt(2)) < 1e-3,
       f"{h6['[111]']:.5f}, ratio {h6['[100]'] / h6['[111]']:.5f}",
       f"{np.sqrt(4/27):.5f}, sqrt(2)"))

    R(("H7  the non-translational Gamma zero mode IS the phase mode",
       h7["nontrans"] == 1 and h7["overlap"] > 1 - 1e-6,
       f"overlap {h7['overlap']:.6f}", "> 0.999999"))
    R(("H7  phase motion costs NO bar energy",
       h7["cost_phase"] < 1e-8, f"{h7['cost_phase']:.2e}", "< 1e-8"))
    R(("H7  CONTROL: a uniform dilation DOES cost energy (so H7 is not vacuous)",
       h7["cost_dilation"] > 0.1, f"{h7['cost_dilation']:.4f}", "> 0.1"))

    R(("H8  PHASE wave speed is 2/3 EXACTLY",
       all(abs(v - 2 / 3) < 1e-4 for v in h8["speeds"].values()),
       ", ".join(f"{v:.6f}" for v in h8["speeds"].values()), "0.666667"))
    R(("H8  PHASE wave is ISOTROPIC -- CAN FAIL",
       max(h8["speeds"].values()) - min(h8["speeds"].values()) < 1e-5,
       f"spread {max(h8['speeds'].values()) - min(h8['speeds'].values()):.2e}", "< 1e-5"))
    R(("H8  PHASE wave is FASTER than either elastic branch",
       min(h8["speeds"].values()) > max(h6.values()),
       f"{min(h8['speeds'].values()):.5f} vs {max(h6.values()):.5f}", "faster"))
    R(("H8  the two soft families are ORTHOGONAL: at M the positional band is "
       "0 while both phase bands sit at 2/sqrt(3)",
       h8["pos_at_M"] >= 1 and all(abs(x - 2 / np.sqrt(3)) < 1e-6 for x in h8["at_M"]),
       f"positional {h8['pos_at_M']} zero mode(s), phase {h8['at_M'][0]:.6f}/"
       f"{h8['at_M'][1]:.6f}", "0 vs 1.154701"))

    print()
    print("=" * 78)
    print(f"GATE  {len(checks)} rows")
    print("=" * 78)
    for name, ok, val, crit in checks:
        print(f"  {'PASS' if ok else 'FAIL':4s}  {name:66s} {str(val):>26s} {str(crit):>18s}")

    print()
    print("  ROWS THAT EXIST ONLY TO STOP ANOTHER ROW BEING UNFALSIFIABLE:")
    print("   * H1's one-shared-vertex CONTROL -- without it 'neighbours share")
    print("     faces' is a claim with nothing to contrast against, and the")
    print("     whole reason this file exists (the codebase builds the wrong")
    print("     packing) goes unstated.")
    print("   * H2's NO-OTHER-CLOSURE row. Without it 'b = a+60' is satisfied by")
    print("     any relation the optimiser happens to walk to -- which is")
    print("     exactly how 'b = 60 + a/2' was published and retracted.")
    print("   * H2's SEED-IS-EXACT row. The corner permutation search is")
    print("     load-bearing; pairing corner c to corner c leaves residual 1.0")
    print("     and the solver converges to a spurious b = 110.237 whose hole")
    print("     cell is not an octahedron. The row fails loudly instead.")
    print("   * H2's SEPARATION row. Driven the other way (a increasing) this")
    print("     same closure drives the cells through each other; without the")
    print("     row the exchange would look valid in both directions.")
    print("   * H5's <100>/<111> row -- without it 'zero modes along <110>' is")
    print("     satisfiable by a framework that is floppy everywhere.")
    print("   * H7's dilation CONTROL -- without it 'the phase mode costs no")
    print("     energy' could be satisfied by a rigidity matrix that is simply")
    print("     small, and the identification would mean nothing.")
    print()
    print("  A ROW THIS FILE DID NOT BUILD, AND NO LONGER OWES: the")
    print("  band-touching PLANE family in the faithful model. H8 checks M")
    print("  only. It is now built, in jb_bt_band_touching.py, and the answer")
    print("  is not the one the reduced model implied: the LOCUS survives")
    print("  exactly (it is the simple-cubic zone boundary) but the degenerate")
    print("  value DISPERSES across each plane, from sqrt(2/3) at the face")
    print("  centre to 2/sqrt(3) at its edges, the bands cross linearly, and")
    print("  no group velocity vanishes. There is no localisation mechanism")
    print("  there. jb_bt also finds that the whole structure -- and H8's own")
    print("  gaplessness -- belongs to a = -30 alone. See below.")
    print()
    print("  WHAT THIS FILE DOES NOT MODEL: anharmonic terms (everything is")
    print("  harmonic about a = -30, so the phase/positional orthogonality is a")
    print("  statement about a LINEAR subspace and need not survive finite")
    print("  amplitude); a lattice that responds dynamically (lambda is held at")
    print("  its reference value while the true coherent motion has lambda")
    print("  tracking the phase); and inviscid-l1d, the PARALLEL_TOL branch flip")
    print("  in jb_z's signed_gap, which is a real and still-unfixed code defect")
    print("  that is independent of everything here.")
    print()
    print("  AND THE FIRST OF THOSE IS NOT A CAVEAT, IT IS THE LEADING TERM.")
    print("  jb_bt_band_touching.py measures the Gamma gap of the phase band at")
    print("  other reference phases and finds gap = sqrt(6)*|d lambda/d a|, to")
    print("  10 figures, isotropic. lambda is the fold half-diagonal and a =")
    print("  -30 is its MAXIMUM, so d lambda/d a is zero there and ONLY there.")
    print("  H8's c = 2/3, its gaplessness and the Goldstone identification are")
    print("  therefore statements about the exchange MIDPOINT under a fixed")
    print("  lattice, not about the medium at a general phase.")

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
    print("jb_hc -- the VE/octa rectified cubic honeycomb, and its two waves")
    print("=" * 78)
    print("  The packing this epic is about: cuboctahedra and octahedra")
    print("  alternating, neighbours sharing FACES. The octahedra ARE collapsed")
    print("  VEs. Driving the exchange DOWNWARD swaps them, with the lattice")
    print("  returning to where it started. As a bar-and-joint framework the")
    print("  medium carries two orthogonal waves: positional (anisotropic, with")
    print("  an exactly flat zero band along the six closest-packing")
    print("  directions) and phase (isotropic, c = 2/3, the Goldstone mode of")
    print("  the exchange).")
    h1 = h1_packing()
    h2 = h2_exchange()
    h3 = h3_cluster()
    fw = h4_framework()
    h5 = h5_zero_modes(fw)
    h6 = h6_elastic(fw)
    h7 = h7_phase_is_a_zero_mode(fw)
    h8 = h8_phase_wave(fw)
    return gate(h1, h2, h3, fw, h5, h6, h7, h8)


if __name__ == "__main__":
    with np.errstate(all="ignore"):
        sys.exit(main())
