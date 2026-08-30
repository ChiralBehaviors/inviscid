"""REVIEW CHECK 10: walk() scores success as index-wise RMS to the target after a
rigid alignment.  That correspondence is FIXED, so a walk that lands on a
tetrahedron with a different corner->vertex labelling is scored as a failure.
Ask instead the labelling-free question: IS the endpoint a regular tetrahedron?
"""
import numpy as np
from jb_a_family import corners, L_EDGE
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import build_tet_config
from jb_e_tighten import is_tet, walk

tets = [s for s in tet_assignments() if is_tet(s)]
X0 = corners(0.0)

print(f"{'t':>3s} {'RMS to target':>14s} {'verdict of walk()':>18s} | "
      f"{'#clusters@1e-9':>15s} {'@1e-6':>7s} {'cluster spread':>15s} "
      f"{'edge min..max':>30s} {'hinge resid':>12s}")
n_tet = 0
for t in range(8):
    Xt = build_tet_config(tets[t])
    X, dist, wa, ws = walk(X0, Xt, n=80)
    pts = X.reshape(-1, 3)

    def nclust(tol):
        u = []
        for p in pts:
            for q in u:
                if np.linalg.norm(p - q) < tol:
                    break
            else:
                u.append(p)
        return u

    u9, u6 = nclust(1e-9), nclust(1e-6)
    # spread inside each of the 1e-6 clusters
    spread = 0.0
    for c in u6:
        mem = [p for p in pts if np.linalg.norm(p - c) < 1e-6]
        spread = max(spread, max(np.linalg.norm(p - q) for p in mem for q in mem))
    D = [np.linalg.norm(u6[i] - u6[j]) for i in range(len(u6)) for j in range(i + 1, len(u6))]
    ok = len(u6) == 4 and max(D) - min(D) < 1e-6 and abs(min(D) - L_EDGE) < 1e-6
    n_tet += ok
    print(f"{t:3d} {dist:14.3e} {('REACHED' if dist < 1e-6 else 'not reached'):>18s} | "
          f"{len(u9):15d} {len(u6):7d} {spread:15.2e} "
          f"{f'{min(D):.9f}..{max(D):.9f}':>30s} "
          f"{np.linalg.norm(residual(X)):12.2e}   {'TETRAHEDRON' if ok else ''}")
print(f"\n  walk()'s own score: 5/8.   labelling-free score: {n_tet}/8 land on a regular tetrahedron.")
