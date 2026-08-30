"""REVIEW CHECK 1: does the Python parameterisation match Jitterbug.java?

Java ground truth dumped from PhiCoordinates.Octahedrons[4] (see DumpJb.java):

  face 0  inverse=false  centroid octant (-,-,-)  product=-1  => sigma_java=+1
  face 1  inverse=false  (+,+,-)  product=-1  => +1
  face 2  inverse=false  (+,-,+)  product=-1  => +1
  face 3  inverse=false  (-,+,+)  product=-1  => +1
  face 4  inverse=true   (+,+,+)  product=+1  => -1
  face 5  inverse=true   (-,-,+)  product=+1  => -1
  face 6  inverse=true   (-,+,-)  product=+1  => -1
  face 7  inverse=true   (+,-,-)  product=+1  => -1

so in Java  sigma_i = -(sx*sy*sz).  The Python model uses sigma_i = +(sx*sy*sz).
"""
import itertools
import numpy as np
from jb_a_family import rot, corners, cluster, Z, R_CIRC, L_EDGE


def corners_signed(a_deg, sign=+1.0):
    """sigma_i = sign * (sx*sy*sz)."""
    out = np.empty((8, 3, 3))
    for i, s in enumerate(itertools.product((1, -1), repeat=3)):
        s = np.array(s, dtype=float)
        v = np.array([[s[0], 0, 0], [0, s[1], 0], [0, 0, s[2]]]) * R_CIRC
        c = v.mean(axis=0)
        u = c / np.linalg.norm(c)
        sigma = sign * s[0] * s[1] * s[2]
        out[i] = (rot(u, sigma * (a_deg - 60.0)) @ (v - c).T).T + u * Z * np.cos(np.radians(a_deg))
    return out


def dist_spectrum(X):
    p = X.reshape(-1, 3)
    D = np.linalg.norm(p[:, None] - p[None, :], axis=-1)
    return np.sort(D[np.triu_indices(24, 1)])


print("=== A. does the JAVA sign convention (sigma = -product) share vertices? ===")
for a in (0, 15, 30, 45, 60, 90):
    rp, mp, _ = cluster(corners_signed(a, +1), tol=1e-7)
    rj, mj, _ = cluster(corners_signed(a, -1), tol=1e-7)
    print(f"  a={a:6.2f}   python(+prod): {len(rp):2d} corners mult{sorted(set(mp.tolist()))}"
          f"     java(-prod): {len(rj):2d} corners mult{sorted(set(mj.tolist()))}")

print("\n=== B. are the two families CONGRUENT at each a? (sorted 276 pairwise distances) ===")
for a in (0, 10, 22.238756093, 30, 45, 59, 60, 61, 75, 90, 120, 180):
    dp, dj = dist_spectrum(corners_signed(a, +1)), dist_spectrum(corners_signed(a, -1))
    print(f"  a={a:10.5f}   max |dist-spectrum difference| = {np.abs(dp - dj).max():.3e}")

print("\n=== C. strut lengths in the JAVA convention ===")
for a in (0, 30, 60, 90):
    X = corners_signed(a, -1)
    e = np.array([np.linalg.norm(X[i, j] - X[i, (j + 1) % 3]) for i in range(8) for j in range(3)])
    print(f"  a={a:6.2f}  strut {e.min():.12f}..{e.max():.12f}")

print("\n=== D. explicit reflection test: is java-family = mirror of python-family? ===")
for M, name in ((np.diag([-1.0, -1, -1]), "point inversion"),
                (np.diag([1.0, 1, -1]), "reflect z"),
                (np.diag([1.0, -1, 1]), "reflect y")):
    worst = 0.0
    for a in (0, 30, 45, 75):
        Xp = corners_signed(a, +1).reshape(-1, 3) @ M.T
        Xj = corners_signed(a, -1).reshape(-1, 3)
        # compare as point SETS
        d = max(min(np.linalg.norm(p - q) for q in Xj) for p in Xp)
        worst = max(worst, d)
    print(f"  {name:18s}: worst set-Hausdorff over a in (0,30,45,75) = {worst:.3e}")

print("\n=== E. does a -> -a or a -> 120-a relate them? ===")
for shift, name in ((lambda a: -a, "a -> -a"), (lambda a: 120 - a, "a -> 120-a")):
    worst = 0.0
    for a in (10, 30, 45, 75):
        Xp = corners_signed(shift(a), +1).reshape(-1, 3)
        Xj = corners_signed(a, -1).reshape(-1, 3)
        d = max(min(np.linalg.norm(p - q) for q in Xj) for p in Xp)
        worst = max(worst, d)
    print(f"  {name:12s}: worst set-Hausdorff = {worst:.3e}")
