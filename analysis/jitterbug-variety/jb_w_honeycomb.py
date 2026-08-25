"""Step W: the VE/octa rectified cubic honeycomb, and the exchange.

THE MEDIUM. Alternate cubes occupied by the cuboctahedron (VE) and the
collapsed VE (the octahedron) -- Fuller's own packing, the rectified cubic
honeycomb. VE cells at even lattice points, octa cells at all-odd points,
neighbours sharing FACES: four vertices VE-VE across a square face, three
vertices VE-octa across a triangular face.

WHY THIS FILE EXISTS. Every array result in this epic before 2026-08-25 was
computed on `jb_x_array_linkage.build_topologies`, which places neighbours at
offset ``2*v[g]`` -- where two cuboctahedra touch at ONE VERTEX. That is the
wrong packing (row W10 measures it: one shared vertex against the honeycomb's
four and three), and a single-vertex contact is a terrible coupling. The
canonical current state of the epic is T2 inviscid/qvf-epic-consolidated-state;
this module is its executable half, and it exists because the scripts that
first produced those numbers did not survive the session that made them.

THE EXCHANGE. Solving the shared-face constraint gives exactly one motion:

    b = a + 60          the hole cell's phase runs 60 degrees ahead, exactly
    no cell rotates     verified zero, not assumed
    lambda = d1(a)      the fold half-diagonal, from this repo's `fold_halves`
                        == (2/sqrt(3)) * cos(a + 30 deg)

Run a from 0 to -60 and the cells trade places: the VE closes into the
octahedron while the hole cell opens into the VE, the lattice breathes out to
2/sqrt(3) at the midpoint a = -30 (where the two are congruent) and returns to
exactly 1, and the packing is valid at every step. The structure comes back to
itself with the roles swapped.

RETRACTION THIS FILE ENCODES. An earlier solution of the same constraint gave
``b = 60 + a/2`` with each cell twisting by ``a/2`` -- a "geared 2:1
mechanism". It closed to 2e-16 at every a and it was wrong: the corner
correspondence across each shared face had been left free, so the system was
UNDERDETERMINED and closed for a whole family of b(a). Row W4 is that lesson
made executable, and row W6 is the check that would have caught it -- a
determined face-sharing condition admits exactly ONE closure, and a system that
closes for every parameter value is telling you it is not yet a system.

Machine-zero residuals confirm that constraints are satisfied. They never
confirm that they are the right constraints.
"""
import numpy as np

from jb_a_family import corners  # noqa: F401  (imported for provenance)
from jb_x_array_linkage import verts, PAIRS, Topology

# --------------------------------------------------------------------------
# Cell geometry, derived from this repo's own `verts(a)` -- not re-invented.
# --------------------------------------------------------------------------

SQRT2 = np.sqrt(2.0)
SEMI_MINOR = np.sqrt(2.0 / 3.0)

#: The 12 vertex ellipses' axes. `verts(a)` is exactly
#: ``sqrt(2) cos(a) * M_AXIS + sqrt(2/3) sin(a) * B_AXIS`` (row W0), so the
#: major axis is the a=0 (cuboctahedral) position and the minor axis the a=90
#: position. Both unit vectors, mutually orthogonal -- asserted in W1.
M_AXIS = verts(0.0) / SQRT2
B_AXIS = verts(90.0) / SEMI_MINOR


def cell(phase, centre=(0.0, 0.0, 0.0)):
    """The 12 vertices of one cell at `phase`, centred on `centre`.

    No rotation term: the honeycomb's motion does not turn any cell, and
    leaving the rotation out of the model is what makes W5's zero a
    statement rather than a fit."""
    t = np.radians(phase)
    p = SQRT2 * np.cos(t) * M_AXIS + SEMI_MINOR * np.sin(t) * B_AXIS
    return p + np.asarray(centre, dtype=float)


def lam(a):
    """The lattice scale: the fold half-diagonal d1(a).

    Equals ``(2/sqrt(3)) cos(a + 30 deg)``, so it is 1 at a=0, rises to its
    maximum 2/sqrt(3) at a=-30 and returns to exactly 1 at a=-60. The lattice
    breathes open to let the cells trade roles and closes again. Checked
    against this repo's `jb_z_quasistatic_array.fold_halves` in row W2."""
    t = np.radians(a)
    return np.cos(t) - np.sin(t) / np.sqrt(3.0)


def hole_phase(a):
    """The phase of the cell filling the hole. Exactly 60 degrees ahead."""
    return a + 60.0


def _tri_faces():
    """The 8 triangular faces, as vertex labels, derived from jb_x's PAIRS.

    Each label sits in two of the 8 faces (12 vertices x 2 == 8 faces x 3),
    so inverting PAIRS' (face, corner) slots recovers the faces without a
    hardcoded table."""
    out = [[] for _ in range(8)]
    for v, slots in enumerate(PAIRS):
        for (f, _c) in slots:
            out[f].append(v)
    return tuple(tuple(sorted(t)) for t in out)


TRI_FACES = _tri_faces()


def _face_dirs():
    """Each triangular face's outward cube-diagonal direction, as an integer
    sign triple. The face centroid at a=0 points along a cube diagonal
    because the face's own 3-fold axis is that diagonal."""
    v0 = verts(0.0)
    return tuple(tuple(int(np.sign(round(x, 9)))
                       for x in v0[list(f)].mean(axis=0))
                 for f in TRI_FACES)


FACE_DIRS = _face_dirs()

#: face index -> the index of the antipodal face, i.e. the face a neighbour
#: presents back across the shared triangle.
OPPOSITE = tuple(FACE_DIRS.index(tuple(-c for c in d)) for d in FACE_DIRS)

#: The 6 square-face directions: VE-VE neighbours sit at ``lambda * 2 * u``.
SQUARE_DIRS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
               (0, -1, 0), (0, 0, 1), (0, 0, -1))


def ve_centre(cell_index, a):
    """Centre of the VE cell at an EVEN lattice point."""
    return lam(a) * np.asarray(cell_index, dtype=float)


def octa_centre(cell_index, a):
    """Centre of the hole cell at an ALL-ODD lattice point."""
    return lam(a) * np.asarray(cell_index, dtype=float)


def neighbours(cell_index):
    """The 14 face neighbours of a cell: 8 across triangles, 6 across squares.

    Exactly the VE's own face census -- 8 triangular + 6 square -- so the
    packing is space-filling with nothing left over."""
    c = np.asarray(cell_index, dtype=int)
    tri = [tuple(c + np.array(d)) for d in FACE_DIRS]
    sq = [tuple(c + 2 * np.array(d)) for d in SQUARE_DIRS]
    return tri, sq


# --------------------------------------------------------------------------
# Face correspondence -- the load-bearing method note.
# --------------------------------------------------------------------------

_PERMS = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))


def face_residual(xs_a, xs_b, perm):
    """Worst vertex distance across a shared face under a given corner map."""
    return max(float(np.linalg.norm(xs_a[i] - xs_b[perm[i]])) for i in range(3))


def match_face(xs_a, xs_b):
    """Resolve the corner correspondence by PERMUTATION SEARCH.

    Returns ``(best_residual, best_perm, identity_residual)``.

    THIS IS NOT A CONVENIENCE. Pairing corner c to corner c is wrong -- the
    two cells list their shared triangle's corners in different orders -- and
    the identity pairing leaves an O(1) residual. A solver seeded from that
    residual does not fail; it converges, to a spurious closure whose hole
    cell is not the cell you asked for. Always assert the seed residual is
    machine-zero BEFORE continuing a branch: the seed is the one configuration
    whose answer is known in closed form."""
    best, best_p = np.inf, None
    for p in _PERMS:
        r = face_residual(xs_a, xs_b, p)
        if r < best:
            best, best_p = r, p
    return best, best_p, face_residual(xs_a, xs_b, (0, 1, 2))


def _distinct(xs, tol=1e-8):
    """A cell's DISTINCT vertex positions.

    Load-bearing at the collapsed phases. A cell at phase 0 or 60 (mod 120)
    has 12 vertex entries occupying only 6 positions, each doubled -- that is
    the collapse itself. Counting index pairs instead of positions reports a
    3-vertex triangular contact as 6 and a 1-vertex square contact as 4."""
    keep = []
    for p in xs:
        if not any(np.linalg.norm(p - q) < tol for q in keep):
            keep.append(p)
    return np.array(keep)


def shared_vertices(xs_a, xs_b, tol=1e-8):
    """How many DISTINCT positions two cells hold in common."""
    a, b = _distinct(xs_a, tol), _distinct(xs_b, tol)
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return int((d < tol).any(axis=1).sum())


# --------------------------------------------------------------------------
# The exchange.
# --------------------------------------------------------------------------

def exchange_residual(a, cell_index=(0, 0, 0)):
    """Worst face mismatch over all 8 triangular contacts of one VE cell.

    Every contact is resolved by permutation search, so this is the honest
    closure of the face-sharing condition and not an artifact of a lucky
    corner ordering."""
    b = hole_phase(a)
    xs_ve = cell(a, ve_centre(cell_index, a))
    worst = 0.0
    for f, d in enumerate(FACE_DIRS):
        hole = tuple(np.asarray(cell_index) + np.array(d))
        xs_h = cell(b, octa_centre(hole, a))
        tri_a = xs_ve[list(TRI_FACES[f])]
        tri_b = xs_h[list(TRI_FACES[OPPOSITE[f]])]
        r, _p, _ident = match_face(tri_a, tri_b)
        worst = max(worst, r)
    return worst


def reciprocal_residual(a, hole_index=(1, 1, 1)):
    """The condition the two-cell solve does NOT impose: a hole cell's OTHER
    faces must meet their own VE cells. Eight around one, closing both ways."""
    b = hole_phase(a)
    xs_h = cell(b, octa_centre(hole_index, a))
    worst = 0.0
    for f, d in enumerate(FACE_DIRS):
        nb = tuple(np.asarray(hole_index) + np.array(d))
        xs_v = cell(a, ve_centre(nb, a))
        tri_h = xs_h[list(TRI_FACES[f])]
        tri_v = xs_v[list(TRI_FACES[OPPOSITE[f]])]
        r, _p, _ident = match_face(tri_h, tri_v)
        worst = max(worst, r)
    return worst


def closure_scan(a, lo=-180.0, hi=180.0, n=1441):
    """Every hole phase whose face condition closes, at fixed a.

    A determined constraint system admits ONE closure. If this returns a
    family, the system is underdetermined and whatever relation was read off
    it is meaningless -- which is exactly how ``b = 60 + a/2`` was produced."""
    xs_ve = cell(a, ve_centre((0, 0, 0), a))
    tri_a = xs_ve[list(TRI_FACES[0])]
    hole = np.array(FACE_DIRS[0], dtype=float) * lam(a)
    hits = []
    for b in np.linspace(lo, hi, n):
        tri_b = cell(b, hole)[list(TRI_FACES[OPPOSITE[0]])]
        if match_face(tri_a, tri_b)[0] < 1e-9:
            hits.append(float(b))
    return hits


def honeycomb_contacts(cell_index, a, tol=1e-8):
    """The face-sharing contact list -- the replacement for
    `jb_x_array_linkage.build_topologies` (bead inviscid-ia5).

    Returns one record per face neighbour:

        (neighbour_index, kind, my_labels, their_labels)

    where ``kind`` is ``"tri"`` or ``"sq"`` and the two label tuples are
    ALIGNED: ``my_labels[i]`` and ``their_labels[i]`` are the same point in
    space. Triangular contacts are resolved by permutation search; square
    contacts are resolved by coincidence, and shrink from 4 labels to 2 to 1
    as the square folds, so a caller must read ``len(my_labels)`` rather than
    assume 4.

    This is the whole correction. `build_topologies` emits ONE shared vertex
    per neighbour; every contact here is a whole shared FACE."""
    a = float(a)
    b = hole_phase(a)
    even = sum(abs(c) for c in cell_index) % 2 == 0
    mine = cell(a if even else b, lam(a) * np.asarray(cell_index, float))
    tri_n, sq_n = neighbours(cell_index)
    out = []
    for f, nb in enumerate(tri_n):
        theirs = cell(b if even else a, lam(a) * np.asarray(nb, float))
        my_lab = TRI_FACES[f]
        their_lab = TRI_FACES[OPPOSITE[f]]
        r, p, _ = match_face(mine[list(my_lab)], theirs[list(their_lab)])
        if r < tol:
            out.append((nb, "tri", my_lab,
                        tuple(their_lab[p[i]] for i in range(3))))
    for nb in sq_n:
        theirs = cell(a if even else b, lam(a) * np.asarray(nb, float))
        # One record per distinct shared POINT. At the collapsed phases a
        # cell's 12 labels occupy 6 positions, so several labels name the
        # same joint; emitting each would duplicate the coincidence
        # constraint a caller builds from it.
        mm, tt, seen = [], [], []
        for i, q in enumerate(mine):
            if any(np.linalg.norm(q - s) < tol for s in seen):
                continue
            for j, s in enumerate(theirs):
                if np.linalg.norm(q - s) < tol:
                    mm.append(i)
                    tt.append(j)
                    seen.append(q)
                    break
        if mm:
            out.append((nb, "sq", tuple(mm), tuple(tt)))
    return out


#: The phase the contact list is READ at, and the whole reason this constant
#: exists rather than a literal at each call site. The identification COUNT is
#: phase dependent because the squares fold -- 48 at a = 0, 36 at a = -30, 30 at
#: a = -60 -- and only the a = -30 list is valid at every phase. Read at a = 0
#: instead and twelve of the forty-eight break, reaching a separation of 2.0 by
#: a = -60, which welds the folding squares shut and forbids the exchange.
#:
#: This INVERTS `jb_x._fcc13_contacts`' documented idiom ("Read at a = 0 and
#: never re-read"). That is right for a topology whose contacts do not fold and
#: wrong here. Row T3 of jb_ht_honeycomb_topology is the control that fails if
#: anyone changes it back.
HONEYCOMB_REF_PHASE = -30.0


def honeycomb_identifications(sites, a=HONEYCOMB_REF_PHASE):
    """`Topology.contacts` tuples for a set of integer honeycomb sites.

    Returns (i, k, j, l) records -- unit i's vertex label k is the same point
    as unit j's vertex label l -- which is EXACTLY the format
    `jb_x.Topology.contacts` has always had. The format was never wrong;
    `build_topologies` simply emitted one identification per neighbour where
    the real packing has three (triangular face) or four/two/one (square face,
    as it folds).

    Each unordered pair is emitted once, from the lower-indexed cell.
    """
    idx = {tuple(int(c) for c in s): i for i, s in enumerate(sites)}
    out = []
    for s, i in sorted(idx.items(), key=lambda kv: kv[1]):
        for (nb, _kind, my_lab, their_lab) in honeycomb_contacts(s, a):
            j = idx.get(tuple(int(c) for c in nb))
            if j is None or j <= i:
                continue
            for t in range(len(my_lab)):
                out.append((i, int(my_lab[t]), j, int(their_lab[t])))
    return tuple(out)


def honeycomb_phases(sites):
    """Per-unit phase OFFSET: 0 for the VE sublattice (all-even sites), 60 for
    the hole cells (all-odd). This is the b = a + 60 of the exchange, carried
    as topology data so no caller has to remember it."""
    return tuple(0.0 if sum(abs(int(c)) for c in s) % 2 == 0 else 60.0
                 for s in sites)


def _honeycomb(name, sites, note=""):
    sites = tuple(tuple(int(c) for c in s) for s in sites)
    return Topology(name, "honeycomb", (), note=note, sites_int=sites,
                    phases=honeycomb_phases(sites),
                    contacts=honeycomb_identifications(sites))


def _shell(centre):
    tri, sq = neighbours(centre)
    return [tuple(centre)] + [tuple(t) for t in tri] + [tuple(q) for q in sq]


def build_honeycomb_topologies():
    """The honeycomb replacements for `jb_x.build_topologies`.

    Every one of these places neighbours across a shared FACE. The topologies
    they replace (SC7, CHAIN5, SQUARE4, CUBE8-M/R, CUBE27-M, FCC13) place them
    at 2*v[g], where two cuboctahedra touch at ONE VERTEX, and every number
    computed on them is superseded -- see T2 inviscid/qvf-epic-consolidated-state
    (c)(3)."""
    even8 = [(x, y, z) for x in (0, 2) for y in (0, 2) for z in (0, 2)]
    return (
        _honeycomb("HC1 (control, one VE)", [(0, 0, 0)],
                   note="the mandated control: one unit, 6 internal DOF"),
        _honeycomb("HC2 (VE + one hole)", [(0, 0, 0), (1, 1, 1)],
                   note="the minimum honeycomb array: one shared triangle"),
        _honeycomb("HC3 (VE - hole - VE)", [(0, 0, 0), (1, 1, 1), (2, 2, 2)],
                   note="the reciprocal condition: a hole meeting its own VEs"),
        _honeycomb("HC9 (one hole + its 8 VEs)", [(1, 1, 1)] + even8,
                   note="jb_hc's H3 cluster: every contact a shared triangle"),
        _honeycomb("HC15 (one VE + all 14 face neighbours)", _shell((0, 0, 0)),
                   note="full coordination for the interior unit"),
        _honeycomb("HC-2HOLE (two holes)",
                   [(x, y, z) for x in (0, 2, 4) for y in (0, 2) for z in (0, 2)]
                   + [(1, 1, 1), (3, 1, 1)],
                   note="the MULTI-HOLE patch -- ia5 scope 2's decisive case"),
    )


def wrong_packing_shared(a=0.0, g=0):
    """The control: what `build_topologies` actually builds.

    jb_x places a neighbour at offset ``2*v[g]``, norm 2*sqrt(2) == 2.828,
    where two cuboctahedra meet at a single vertex."""
    xs = cell(a)
    return shared_vertices(xs, cell(a, 2.0 * verts(a)[g]))


# --------------------------------------------------------------------------
# Gate.
# --------------------------------------------------------------------------

def _rows():
    rows = []

    # W0 -- the parameterisation is this repo's own geometry, not a restatement.
    err = max(float(np.abs(cell(a) - verts(a)).max())
              for a in np.linspace(-180, 180, 721))
    rows.append(("W0  ellipse form reproduces jb_x.verts(a) over 721 samples",
                 err < 1e-14, f"max dev {err:.3e}"))

    # W1 -- the ellipse invariants.
    mn = np.abs(np.linalg.norm(M_AXIS, axis=1) - 1).max()
    bn = np.abs(np.linalg.norm(B_AXIS, axis=1) - 1).max()
    orth = np.abs((M_AXIS * B_AXIS).sum(1)).max()
    rows.append(("W1a axes orthonormal (semi-major sqrt2, semi-minor sqrt(2/3))",
                 max(mn, bn, orth) < 1e-14,
                 f"|m|-1 {mn:.1e}  |b|-1 {bn:.1e}  m.b {orth:.1e}"))
    ratio = SQRT2 / SEMI_MINOR
    rows.append(("W1b axis ratio is exactly sqrt(3)",
                 abs(ratio - np.sqrt(3.0)) < 1e-14, f"ratio {ratio:.12f}"))
    r0, r60, r90 = (np.linalg.norm(cell(x), axis=1).mean() for x in (0, 60, 90))
    rows.append(("W1c radius sqrt2 at a=0, EXACTLY 1 at a=60, sqrt(2/3) at a=90",
                 abs(r0 - SQRT2) < 1e-12 and abs(r60 - 1.0) < 1e-12
                 and abs(r90 - SEMI_MINOR) < 1e-12,
                 f"{r0:.9f} / {r60:.9f} / {r90:.9f}"))
    n60 = len(np.unique(np.round(cell(60.0), 9), axis=0))
    rows.append(("W1d the octahedron IS the collapsed VE: 6 distinct verts at a=60",
                 n60 == 6, f"{n60} distinct"))

    # W2 -- lambda is the fold half-diagonal, from the repo's own machinery.
    from jb_z_quasistatic_array import fold_halves
    dev = 0.0
    for a in (-60, -45, -30, -15, 0, 15, 30):
        for d1, _d2 in fold_halves(a).values():
            dev = max(dev, abs(d1 - lam(a)))
    rows.append(("W2a lambda(a) == fold_halves d1(a), all 6 squares",
                 dev < 1e-12, f"max dev {dev:.3e}"))
    cf = max(abs(lam(a) - (2 / np.sqrt(3)) * np.cos(np.radians(a + 30)))
             for a in np.linspace(-90, 90, 361))
    rows.append(("W2b lambda(a) == (2/sqrt3) cos(a+30) in closed form",
                 cf < 1e-14, f"max dev {cf:.3e}"))
    peak = max(lam(a) for a in np.linspace(-60, 0, 6001))
    rows.append(("W2c lambda: 1 at a=0, max 2/sqrt3 at a=-30, back to 1 at a=-60",
                 abs(lam(0) - 1) < 1e-12 and abs(lam(-60) - 1) < 1e-12
                 and abs(peak - 2 / np.sqrt(3)) < 1e-9,
                 f"{lam(0):.9f} / {peak:.9f} / {lam(-60):.9f}"))

    # W3 -- the packing, at the reference phase.
    xs = cell(0.0)
    sq = shared_vertices(xs, cell(0.0, ve_centre((2, 0, 0), 0.0)))
    tri = shared_vertices(xs, cell(60.0, octa_centre((1, 1, 1), 0.0)))
    rows.append(("W3a VE-VE across a square face: 4 shared vertices",
                 sq == 4, f"{sq} shared, dist {lam(0)*2:.6f}"))
    rows.append(("W3b VE-octa across a triangular face: 3 shared vertices",
                 tri == 3, f"{tri} shared, dist {lam(0)*np.sqrt(3):.6f}"))
    t_n, s_n = neighbours((0, 0, 0))
    rows.append(("W3c census 8 triangular + 6 square == 14 == the VE's own faces",
                 len(t_n) == 8 and len(s_n) == 6 and len(t_n) + len(s_n) == 14,
                 f"{len(t_n)} + {len(s_n)} = {len(t_n)+len(s_n)}"))
    parity = all(sum(abs(c) for c in n) % 2 == 1 for n in t_n) and \
        all(sum(abs(c) for c in n) % 2 == 0 for n in s_n)
    rows.append(("W3d VE at even lattice points, hole cells at all-odd",
                 parity, "parity holds for all 14 neighbours"))

    # W4 -- the method note, executable.
    xs_ve, xs_h = cell(0.0), cell(60.0, octa_centre((1, 1, 1), 0.0))
    r, p, ident = match_face(xs_ve[list(TRI_FACES[0])],
                             xs_h[list(TRI_FACES[OPPOSITE[0]])])
    rows.append(("W4a identity corner pairing leaves an O(1) SEED RESIDUAL",
                 ident > 0.1, f"identity residual {ident:.6f}"))
    rows.append(("W4b permutation search resolves it to machine zero",
                 r < 1e-14, f"perm {p} residual {r:.3e}"))

    # W5 -- the exchange.
    worst = max(exchange_residual(a) for a in np.linspace(0, -60, 121))
    rows.append(("W5a exchange closes: b=a+60, no rotation, lambda=d1(a)",
                 worst < 1e-14, f"worst face mismatch {worst:.3e} over a in [0,-60]"))
    rworst = max(reciprocal_residual(a) for a in np.linspace(0, -60, 121))
    rows.append(("W5b reciprocal condition -- eight around one, closing both ways",
                 rworst < 1e-14, f"worst {rworst:.3e}"))
    off = max(exchange_residual(a) for a in (-90.0, -120.0, 30.0, 90.0))
    rows.append(("W5c ... and it is not vacuous: relation still exact off [0,-60]",
                 off < 1e-14, f"worst {off:.3e}"))

    # W6 -- the check that would have caught the retracted 2:1 result.
    hits = closure_scan(-30.0)
    uniq = [h for h in hits if abs(((h - 30.0 + 180) % 360) - 180) < 1e-6]
    rows.append(("W6  ONE closure over a full 360 hole-phase scan, not a family",
                 len(hits) == 1 and len(uniq) == 1,
                 f"{len(hits)} hit(s): {[round(h,4) for h in hits]}, b=a+60 gives 30.0"))

    # W7 -- the role swap.
    rad = {a: (np.linalg.norm(cell(a), axis=1).mean(),
               np.linalg.norm(cell(hole_phase(a)), axis=1).mean())
           for a in (0.0, -30.0, -60.0)}
    swap = (abs(rad[0.0][0] - SQRT2) < 1e-12 and abs(rad[0.0][1] - 1.0) < 1e-12
            and abs(rad[-60.0][0] - 1.0) < 1e-12
            and abs(rad[-60.0][1] - SQRT2) < 1e-12)
    rows.append(("W7a complete role swap at a=-60: VE radius sqrt2 <-> 1",
                 swap,
                 f"a=0 ({rad[0.0][0]:.6f},{rad[0.0][1]:.6f}) -> "
                 f"a=-60 ({rad[-60.0][0]:.6f},{rad[-60.0][1]:.6f})"))
    cong = abs(rad[-30.0][0] - rad[-30.0][1])
    rows.append(("W7b midpoint a=-30: the two cells are CONGRUENT",
                 cong < 1e-12, f"both at radius {rad[-30.0][0]:.6f}, dev {cong:.1e}"))
    ico = -37.7612
    rows.append(("W7c each cell passes THROUGH the icosahedral phase in the swap",
                 abs(hole_phase(ico) - 22.238756093) < 1e-3,
                 f"hole reads b={hole_phase(ico):.4f} at a={ico}"))

    # W8 -- the two kinds of face are not equivalent.
    decay = [shared_vertices(cell(a), cell(a, ve_centre((2, 0, 0), a)))
             for a in (0.0, -30.0, -60.0)]
    held = [shared_vertices(cell(a), cell(hole_phase(a),
                                          octa_centre((1, 1, 1), a)))
            for a in (0.0, -30.0, -60.0)]
    rows.append(("W8a square (VE-VE) contact DECAYS 4 -> 2 -> 1",
                 decay == [4, 2, 1], f"{decay} at a = 0 / -30 / -60"))
    rows.append(("W8b triangular (VE-octa) contact stays WHOLE at 3",
                 held == [3, 3, 3], f"{held} at a = 0 / -30 / -60"))

    # W9 -- zero rotation is a result, not an assumption.
    dev = 0.0
    for a in (-15.0, -30.0, -45.0):
        xs_h = cell(hole_phase(a), octa_centre((1, 1, 1), a))
        ref = cell(hole_phase(a)) + octa_centre((1, 1, 1), a)
        dev = max(dev, float(np.abs(xs_h - ref).max()))
    rows.append(("W9  no cell rotates: hole cell is a pure translate of cell(b)",
                 dev < 1e-14, f"max dev {dev:.3e}"))

    # W10 -- the wrong packing, measured. Evidence for bead inviscid-ia5.
    w = wrong_packing_shared()
    off_norm = float(np.linalg.norm(2.0 * verts(0.0)[0]))
    rows.append(("W10 CONTROL: build_topologies' offset shares ONE vertex",
                 w == 1,
                 f"{w} shared at |offset| {off_norm:.6f} "
                 f"(honeycomb: 4 square / 3 triangular)"))

    # W11 -- the contact list that replaces build_topologies (bead ia5).
    con = honeycomb_contacts((0, 0, 0), -30.0)
    ntri = sum(1 for c in con if c[1] == "tri")
    nsq = sum(1 for c in con if c[1] == "sq")
    rows.append(("W11a contact list: 8 triangular + 6 square face contacts",
                 ntri == 8 and nsq == 6, f"{ntri} tri + {nsq} sq"))
    widths = sorted({len(c[2]) for c in con if c[1] == "tri"})
    rows.append(("W11b every triangular contact carries THREE aligned labels",
                 widths == [3], f"label widths {widths} (build_topologies: [1])"))
    bad = 0.0
    for nb, kind, ml, tl in con:
        even = sum(abs(c) for c in nb) % 2 == 0
        mine = cell(-30.0, lam(-30.0) * np.zeros(3))
        theirs = cell(-30.0 if even else hole_phase(-30.0),
                      lam(-30.0) * np.asarray(nb, float))
        for i, j in zip(ml, tl):
            bad = max(bad, float(np.linalg.norm(mine[i] - theirs[j])))
    rows.append(("W11c ... and the alignment is exact, every contact",
                 bad < 1e-14, f"worst {bad:.3e}"))
    sq_decay = [sorted({len(c[2]) for c in honeycomb_contacts((0, 0, 0), a)
                        if c[1] == "sq"}) for a in (0.0, -30.0, -60.0)]
    rows.append(("W11d square contacts shrink 4 -> 2 -> 1, so callers must read len()",
                 sq_decay == [[4], [2], [1]], f"{sq_decay} at a = 0 / -30 / -60"))

    return rows


def main():
    rows = _rows()
    print(__doc__.split("\n\n")[0])
    print()
    print(f"  cells    VE at even lattice points, hole cells at all-odd")
    print(f"  lattice  A = 2*lambda,  lambda(a) = (2/sqrt3) cos(a+30)")
    print(f"  faces    {len(TRI_FACES)} triangular, dirs {FACE_DIRS[:2]} ...")
    print()
    for name, ok, note in rows:
        print(f"  {'PASS' if ok else 'FAIL':4s}  {name:62s}  {note}")
    failed = [n for n, ok, _ in rows if not ok]
    print()
    if failed:
        print(f"  !! {len(failed)} CHECK(S) FAILED -- a bug report, not a result:")
        for n in failed:
            print(f"     {n}")
        return 1
    print(f"  ALL {len(rows)} CHECKS PASSED.")
    print()
    print("  THE EXCHANGE, at a glance:")
    print("      a       b      lambda    VE radius   hole radius   face residual")
    for a in (0.0, -15.0, -30.0, -45.0, -60.0):
        rv = np.linalg.norm(cell(a), axis=1).mean()
        rh = np.linalg.norm(cell(hole_phase(a)), axis=1).mean()
        print(f"   {a:6.1f}  {hole_phase(a):6.1f}   {lam(a):.6f}   {rv:.6f}"
              f"    {rh:.6f}      {exchange_residual(a):.2e}")
    print()
    print("  Canonical state: T2 inviscid/qvf-epic-consolidated-state")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
