"""Does MY demonstrated VE->tetrahedron path stay collision-free?

The reachability claim is mine; the qualification has to be measured on my own
path, not inherited from the critic's.
"""
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import build_tet_config, kabsch
from jb_e_tighten import is_tet, project
from verify_critic import interpenetrations


def walk_recording(X0, Xt, n=80, pull=1.0, w=1e4):
    z, path = np.zeros(48), [X0.copy()]
    for s in np.linspace(1.0 / n, 1.0, n):
        Xa = kabsch(Xt.reshape(-1, 3), place(X0, z).reshape(-1, 3)).reshape(8, 3, 3)
        way = (1 - s) * place(X0, z) + s * Xa
        sol = least_squares(lambda zz: np.concatenate(
            [w * residual(place(X0, zz)), pull * (place(X0, zz) - way).reshape(-1)]),
            z, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        z, _ = project(X0, sol.x)
        path.append(place(X0, z))
    return path


if __name__ == "__main__":
    tets = [s for s in tet_assignments() if is_tet(s)]
    X0 = corners(0.0)
    for t in (1, 2):
        path = walk_recording(X0, build_tet_config(tets[t]), n=80)
        counts = [interpenetrations(X) for X in path]
        bad = sum(c > 0 for c in counts)
        print(f"target {t}: {bad}/{len(path)} waypoints interpenetrate  "
              f"(max simultaneous {max(counts)}, endpoint {counts[-1]})")
        first = next((i for i, c in enumerate(counts) if c > 0), None)
        print(f"          first violation at waypoint {first}"
              f"{'' if first is None else f' of {len(path)}'}")
