"""REVIEW CHECK 6: how big are walk()'s steps relative to the whole journey, and
does the reachability verdict survive refinement of the discretisation?
"""
import sys
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners, L_EDGE
from jb_b_variety import jacobian
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import TET, build_tet_config, kabsch
from jb_e_tighten import is_tet, project
from rv_4_continuation import shape_dist
from rv_5_walk import describe


def walk_rec(X0, Xt, n=80, pull=1.0, w=1e4):
    z = np.zeros(48)
    steps, configs = [], []
    for s in np.linspace(1.0 / n, 1.0, n):
        Xprev = place(X0, z)
        Xa = kabsch(Xt.reshape(-1, 3), Xprev.reshape(-1, 3)).reshape(8, 3, 3)
        way = (1 - s) * Xprev + s * Xa
        sol = least_squares(lambda zz: np.concatenate(
            [w * residual(place(X0, zz)), pull * (place(X0, zz) - way).reshape(-1)]),
            z, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        zp, _ = project(X0, sol.x)
        Xnew = place(X0, zp)
        steps.append(shape_dist(Xnew, Xprev))
        configs.append(Xnew)
        z = zp
    X = place(X0, z)
    Xa = kabsch(Xt.reshape(-1, 3), X.reshape(-1, 3)).reshape(8, 3, 3)
    return X, np.linalg.norm(X - Xa) / np.sqrt(24), np.array(steps), configs


def n_tet_points(X):
    pts, uniq = X.reshape(-1, 3), []
    for p in pts:
        if not any(np.linalg.norm(p - q) < 1e-6 for q in uniq):
            uniq.append(p)
    return len(uniq)


if __name__ == "__main__":
    tets = [s for s in tet_assignments() if is_tet(s)]
    X0 = corners(0.0)

    print("=== A. scale: how far is the tetrahedron from the VE? ===")
    for t in (0, 1, 5):
        Xt = build_tet_config(tets[t])
        print(f"  target {t}: shape-quotient distance VE -> tetrahedron = "
              f"{shape_dist(Xt, X0):.4f}   (config coords, strut length {L_EDGE:.4f})")

    print("\n=== B. does the verdict survive refining the discretisation? ===")
    print(f"  {'target':>7s} {'n':>5s} {'final RMS':>12s} {'max step':>10s} "
          f"{'max/journey':>12s} {'#steps>0.1':>11s} {'reached tet?':>28s}")
    for t in [int(x) for x in sys.argv[1:]] or [1, 5]:
        Xt = build_tet_config(tets[t])
        journey = shape_dist(Xt, X0)
        for n in (80, 200, 400):
            X, d, st, cfg = walk_rec(X0, Xt, n=n)
            k = n_tet_points(X)
            print(f"  {t:7d} {n:5d} {d:12.3e} {st.max():10.3e} "
                  f"{st.max()/journey:12.3f} {int(np.sum(st > 0.1)):11d} "
                  f"{('YES (' + str(k) + ' pts)') if k == 4 else ('no (' + str(k) + ' pts)'):>28s}")
            if n == 80:
                first = next((i for i, c in enumerate(cfg) if n_tet_points(c) == 4), None)
                print(f"          -> first step at which the config IS a 4-point tetrahedron: {first}"
                      f" of {n}   (afterwards the walk is pinned and adds nothing)")
