"""Step H (M0): control -- reproduce C5's hull-volume landmarks, and check the
two claims the V = -k*Vol_hull hypothesis silently rests on.

C5 in T2 `inviscid/jitterbug-geometry-and-derived-dynamics.md` records:
    V/V_tet = 20.000000 at a=0 (a LOCAL MINIMUM), maxima 20.4430915 at
    a = +/-7.2189519, 18.512296 at a_ico = 22.238756093, 4 at a=60.
That is inherited, so it is re-measured here from scratch. Nothing downstream is
worth anything if these do not come back.

THE KNOWN HULL TRAP, guarded explicitly below: a brute-force "accept every
supporting triple" hull OVERCOUNTS non-simplicial faces. The cuboctahedron's 6
SQUARE faces each admit 4 coplanar supporting triples where only 2 triangles
belong, which inflates V(0) from 20 to 32. `_brute_hull_volume` reproduces that
failure on purpose so the guard is a check that CAN fail; the measurements use
scipy's Qhull.

Two structural claims are also tested here because M1/M2/M3 lean on them:

  COSPHERICITY. C5 asserts all 12 shared vertices are cospherical at every a
  (r: 1.4142136 -> 1.0). If true, all 12 are extreme points of the hull at every
  a and none ever leaves it -- which is the whole argument that Vol_hull could be
  smooth. Measured here independently rather than inherited.

  SQUARE-FACE PLANARITY. The 6 square faces of the cuboctahedron are what makes
  a=0 non-simplicial. Whether their 4 vertices STAY coplanar off a=0, and at what
  order in a they leave, decides whether a=0 is a special point for the hull's
  combinatorics -- and therefore whether Vol_hull can be C2 there. This is the
  measurement that sets up M1.
"""
import numpy as np
from scipy.spatial import ConvexHull

from jb_a_family import corners, cluster, L_EDGE
from jb_b_variety import PAIRS

V_TET = L_EDGE ** 3 / (6 * np.sqrt(2.0))
A_ICO = 22.238756093


def vertices(a_deg, X=None):
    """The 12 shared vertices, in a labelling that is STABLE across a.

    `cluster` relabels by first encounter and collapses 12 -> 6 at the merge
    angles, so it cannot be used to track a vertex. PAIRS is combinatorial and
    fixed, so vertex v is simply corner PAIRS[v][0], which exists at every a.
    """
    if X is None:
        X = corners(a_deg)
    return np.array([X[i, j] for (i, j), _ in PAIRS])


def hull_volume(a_deg, X=None):
    """Convex-hull volume of the 12 shared vertices via Qhull."""
    return ConvexHull(vertices(a_deg, X)).volume


def _brute_hull_volume(a_deg):
    """DELIBERATELY WRONG hull, kept as a tripwire for the overcount trap.

    Accepts every supporting triple as a facet and sums the tetrahedra to the
    centroid. On a simplicial hull it agrees with Qhull; on the cuboctahedron
    each square face contributes 4 coplanar triples instead of 2 triangles.
    """
    P = vertices(a_deg)
    c = P.mean(axis=0)
    n = len(P)
    total = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                nrm = np.cross(P[j] - P[i], P[k] - P[i])
                s = np.linalg.norm(nrm)
                if s < 1e-12:
                    continue
                nrm = nrm / s
                d = P @ nrm - P[i] @ nrm
                if d.max() < 1e-9 or d.min() > -1e-9:      # supporting plane
                    total += abs((P[i] - c) @ np.cross(P[j] - c, P[k] - c)) / 6.0
    return total


def radii(a_deg):
    return np.linalg.norm(vertices(a_deg), axis=1)


def square_faces():
    """The 6 square faces of the cuboctahedron, as 4-cycles of vertex indices.

    Read off at a=0 from the hull itself (Qhull merges the coplanar triangles
    into one 4-vertex facet), then held fixed -- the labelling is combinatorial.
    """
    P = vertices(0.0)
    h = ConvexHull(P, qhull_options="Qt")
    out = []
    for eq, simplices in _facet_groups(P, h):
        if len(simplices) == 4:
            out.append(sorted(simplices))
    return out


def _facet_groups(P, h, tol=1e-9):
    """Group hull simplices by their supporting plane."""
    groups = {}
    for s, eq in zip(h.simplices, h.equations):
        key = None
        for k in groups:
            if np.allclose(np.array(k), eq, atol=tol):
                key = k
                break
        if key is None:
            key = tuple(eq)
            groups[key] = set()
        groups[key].update(s.tolist())
    return [(np.array(k), sorted(v)) for k, v in groups.items()]


def square_out_of_plane(a_deg, quad):
    """Signed volume of the tetrahedron on a square face's 4 vertices.

    Zero iff the 4 vertices are coplanar. This is the quantity whose ABSOLUTE
    VALUE the convex hull adds to the volume when the quad goes skew, so its
    order of vanishing in `a` is exactly what decides smoothness at a=0.
    """
    P = vertices(a_deg)
    A, B, C, D = (P[i] for i in quad)
    return np.dot(B - A, np.cross(C - A, D - A)) / 6.0


if __name__ == "__main__":
    print("=== M0: hull-volume landmarks (normalised V/V_tet) ===")
    print(f"    V_tet = {V_TET:.12f}")
    targets = [("VE            a=0", 0.0, 20.000000),
               ("peak          a=+7.2189519", 7.2189519, 20.4430915),
               ("peak          a=-7.2189519", -7.2189519, 20.4430915),
               ("icosahedron   a=22.238756093", A_ICO, 18.512296),
               ("octahedron    a=60", 60.0, 4.0)]
    ok = True
    for name, a, want in targets:
        got = hull_volume(a) / V_TET
        d = abs(got - want)
        ok &= d < 5e-6
        print(f"  {name:32s} {got:.7f}   memo {want:.7f}   |diff| {d:.2e}"
              f"   {'OK' if d < 5e-6 else 'MISMATCH'}")
    print(f"  --> M0 {'REPRODUCED' if ok else 'FAILED'}")

    print("\n=== the overcount trap, reproduced so the guard can fail ===")
    for a in (0.0, 7.2189519, A_ICO):
        q, b = hull_volume(a) / V_TET, _brute_hull_volume(a) / V_TET
        tag = "OVERCOUNTS" if abs(q - b) > 1e-9 else "agrees (hull is simplicial)"
        print(f"  a={a:12.7f}  Qhull {q:10.6f}   brute {b:10.6f}   {tag}")

    print("\n=== locate the maximum independently (no inherited 7.2189519) ===")
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda a: -hull_volume(a) / V_TET,
                        bracket=(2.0, 8.0, 15.0), method="brent",
                        options={"xtol": 1e-12})
    print(f"  argmax a = {r.x:.9f} deg   V/V_tet = {-r.fun:.9f}")
    print(f"  memo     a = 7.2189519 deg  V/V_tet = 20.4430915")

    print("\n=== is V even about a=0? (C5 says to ten digits) ===")
    for a in (0.5, 2.0, 7.2189519, 20.0, 45.0):
        vp, vm = hull_volume(a) / V_TET, hull_volume(-a) / V_TET
        print(f"  a={a:11.7f}  V(+a)-V(-a) = {vp - vm:+.3e}")

    print("\n=== COSPHERICITY of the 12 shared vertices (C5's smoothness argument) ===")
    print("     a        r_min          r_max          spread")
    worst = 0.0
    for a in np.concatenate([np.linspace(0, 180, 19), [A_ICO, 7.2189519, 59, 61, 90]]):
        r = radii(a)
        sp = r.max() - r.min()
        worst = max(worst, sp)
        print(f"  {a:7.3f}   {r.min():.9f}    {r.max():.9f}    {sp:.2e}")
    print(f"  --> worst spread over the sweep: {worst:.2e}"
          f"   ({'COSPHERICAL' if worst < 1e-12 else 'NOT cospherical'})")
    print(f"  r(0) = {radii(0.0)[0]:.9f} (memo 1.4142136), "
          f"r(60) = {radii(60.0)[0]:.9f} (memo 1.0)")
    print("  NOTE: cospherical + 12 points => every vertex is extreme => the hull")
    print("  never loses a vertex. That is necessary for smoothness, NOT sufficient:")
    print("  the FACE combinatorics can still change. See the square-face test.")

    SQ = square_faces()
    print(f"\n=== SQUARE FACES of the cuboctahedron: found {len(SQ)} 4-vertex facets ===")
    for q in SQ:
        print(f"  {q}")
    print("\n  out-of-plane tetrahedron volume of one square's 4 vertices vs a:")
    print("     a           tet vol         |tet vol|/a       |tet vol|/a^2")
    for a in (1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1.0, 2.0, 5.0, 7.2189519):
        t = square_out_of_plane(a, SQ[0])
        print(f"  {a:11.6f}   {t:+.9e}   {abs(t)/a:.9e}   {abs(t)/a**2:.9e}")
    print("\n  and on the other side (evenness of |tet vol| in a):")
    for a in (1e-3, 1e-2, 1e-1, 1.0):
        print(f"  a=+-{a:<8.4f} tet(+a) = {square_out_of_plane(a, SQ[0]):+.6e}   "
              f"tet(-a) = {square_out_of_plane(-a, SQ[0]):+.6e}")
