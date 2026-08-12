"""Step A: reproduce the symmetric jitterbug family headlessly and pin its facts.

Parameterisation derived from Jitterbug.java (construct() bakes vertices at
-sigma*60 about the face axis with pivot at the centroid, then translates by
-centroid and RESETS the live transforms; rotateTo(a) then applies sigma*a about
the same axis with the same pivot, and translate() sets u*Z*cos(a)).

Because the pivot c is parallel to the axis u, R(u,theta)*c == c, so the
pivot contribution (c - R c) vanishes identically and the whole stack collapses to

    x(a) = R(u, sigma*(a - 60deg)) @ (v - c)  +  u * Z * cos(a)

with v an octahedron face vertex, c the face centroid, u = c/|c|.

SIGN CONVENTION: these scripts use sigma = +(sx*sy*sz). The Java
(PhiCoordinates.Octahedrons[4] together with JITTERBUG_INVERSES) uses
sigma = -(sx*sy*sz). Verified immaterial: the two families are related by
a -> -a exactly (set-Hausdorff distance 4e-16 between them), with identical
distance spectra and the octahedron at a = +/-60. No measured quantity changes.
Do NOT "fix" the sign -- it would only mirror every result.
"""
import itertools
import numpy as np

R_CIRC = 1.0                      # octahedron circumradius
L_EDGE = R_CIRC * np.sqrt(2.0)    # octahedron edge length
Z = L_EDGE * np.sqrt(2.0) / np.sqrt(3.0)   # == 2R/sqrt3, matches Jitterbug.java


def rot(axis, deg):
    """Rodrigues rotation matrix, axis assumed unit."""
    t = np.radians(deg)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * (K @ K)


SIGNS = [np.array(s, dtype=float) for s in itertools.product((1, -1), repeat=3)]


def faces(uniform_sigma=False):
    """The 8 octahedron faces: (vertices, centroid, axis, sigma)."""
    out = []
    for s in SIGNS:
        v = np.array([[s[0], 0, 0], [0, s[1], 0], [0, 0, s[2]]]) * R_CIRC
        c = v.mean(axis=0)
        u = c / np.linalg.norm(c)
        sigma = 1.0 if uniform_sigma else s[0] * s[1] * s[2]
        out.append((v, c, u, sigma))
    return out


def corners(a_deg, uniform_sigma=False):
    """24 corner positions, shape (8, 3, 3): [face][corner][xyz]."""
    F = faces(uniform_sigma)
    out = np.empty((8, 3, 3))
    for i, (v, c, u, sigma) in enumerate(F):
        Rm = rot(u, sigma * (a_deg - 60.0))
        out[i] = (Rm @ (v - c).T).T + u * Z * np.cos(np.radians(a_deg))
    return out


def cluster(points, tol=1e-9):
    """Group near-identical points. Returns (reps, multiplicities, labels)."""
    pts = points.reshape(-1, 3)
    reps, labels = [], []
    for p in pts:
        for k, r in enumerate(reps):
            if np.linalg.norm(p - r) < tol:
                labels.append(k)
                break
        else:
            labels.append(len(reps))
            reps.append(p)
    mult = np.bincount(labels, minlength=len(reps))
    return np.array(reps), mult, np.array(labels)


def edge_lengths(a_deg, uniform_sigma=False):
    """All 24 strut lengths (3 per triangle)."""
    X = corners(a_deg, uniform_sigma)
    return np.array([np.linalg.norm(X[i, j] - X[i, (j + 1) % 3])
                     for i in range(8) for j in range(3)])


def convex_volume(a_deg):
    from scipy.spatial import ConvexHull
    reps, _, _ = cluster(corners(a_deg), tol=1e-7)
    return ConvexHull(reps).volume


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    print(f"R={R_CIRC}  L_edge={L_EDGE:.6f}  Z={Z:.6f}  (Z should be 2/sqrt3={2/np.sqrt(3):.6f})")

    print("\n--- struts are rigid and equal at every a (they must be: rigid triangles) ---")
    for a in (0, 22.24, 30, 60, 90, 180):
        e = edge_lengths(a)
        print(f"  a={a:7.2f}  strut len min/max = {e.min():.12f} / {e.max():.12f}")

    print("\n--- vertex sharing, sigma = sx*sy*sz (FORCED) ---")
    for a in (0, 15, 22.24, 30, 45, 59, 60, 61, 90, 180):
        reps, mult, _ = cluster(corners(a), tol=1e-7)
        print(f"  a={a:7.2f}  distinct corners = {len(reps):3d}  multiplicities = {sorted(set(mult.tolist()))}")

    print("\n--- vertex sharing, sigma = +1 UNIFORM (the wrong choice) ---")
    for a in (0, 15, 30, 60):
        reps, mult, _ = cluster(corners(a, uniform_sigma=True), tol=1e-7)
        print(f"  a={a:7.2f}  distinct corners = {len(reps):3d}  multiplicities = {sorted(set(mult.tolist()))}")

    print("\n--- pinned forms ---")
    reps0, _, _ = cluster(corners(0.0), tol=1e-7)
    d0 = [np.linalg.norm(p) for p in reps0]
    print(f"  a=0    : {len(reps0)} vertices, circumradius {min(d0):.9f}..{max(d0):.9f}"
          f"  (VE: radius == edge {L_EDGE:.9f}? {np.allclose(d0, L_EDGE)})")

    # icosahedron: all 30 short edges equal. The 24 struts are fixed; the other 6
    # are the "open square" diagonals. Find a where they match the strut length.
    from scipy.optimize import brentq

    def min_nonstrut(a):
        """Shortest distance between two shared vertices NOT joined by a strut.

        Every edge of the cuboctahedron is a strut, so the nearest non-strut pair
        is the short diagonal of a folding square. The icosahedron is where that
        diagonal reaches strut length: all 30 short edges equal."""
        reps, _, labels = cluster(corners(a), tol=1e-7)
        lab = labels.reshape(8, 3)
        strut = {tuple(sorted((lab[i, j], lab[i, (j + 1) % 3])))
                 for i in range(8) for j in range(3)}
        best = np.inf
        for p in range(len(reps)):
            for q in range(p + 1, len(reps)):
                if (p, q) in strut:
                    continue
                best = min(best, np.linalg.norm(reps[p] - reps[q]))
        return best

    a_ico = brentq(lambda a: min_nonstrut(a) - L_EDGE, 1.0, 40.0, xtol=1e-12)
    print(f"  a_ico  : {a_ico:.9f} deg   (memo says 22.24)")
    reps_i, _, _ = cluster(corners(a_ico), tol=1e-7)
    Di = np.linalg.norm(reps_i[:, None, :] - reps_i[None, :, :], axis=-1)
    short = np.sort(Di[Di > 1e-9])[:60]
    print(f"           30 shortest edges span {short.min():.12f}..{short.max():.12f}"
          f"  (count within 1e-9 of min: {np.sum(short < short.min()+1e-9)//2})")

    print("\n--- volumes, normalised so the tetrahedron of edge L is 1 ---")
    V_tet = L_EDGE ** 3 / (6 * np.sqrt(2))
    for name, a in (("VE (a=0)", 0.0), ("icosa", a_ico), ("octa (a=60)", 60.0)):
        print(f"  {name:14s} V/V_tet = {convex_volume(a)/V_tet:.6f}")
    print(f"  {'tetrahedron':14s} V/V_tet = 1.000000  (by definition)")
