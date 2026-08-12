"""REVIEW CHECK 9: is the variety actually SINGULAR at a=90 / at the tetrahedron,
or merely a point where the defining map loses transversality?

Estimate the LOCAL DIMENSION of the solution set: perturb by eps in random
directions, project onto {residual = 0}, quotient out the global rigid motions,
and take the SVD of the resulting displacement cloud.  A smooth 6-manifold gives
6 significant directions at every eps.  A branch point / higher stratum gives 7.
"""
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners
from jb_b_variety import jacobian
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import build_tet_config
from jb_e_tighten import is_tet
from rv_4_continuation import kabsch_align


def local_cloud(Xc, eps, K=80, seed=0):
    rng = np.random.default_rng(seed)
    disp = []
    for _ in range(K):
        z0 = rng.normal(size=48)
        z0 *= eps / np.linalg.norm(z0)
        sol = least_squares(lambda zz: residual(place(Xc, zz)), z0,
                            xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if np.linalg.norm(residual(place(Xc, sol.x))) > 1e-11:
            continue
        X = place(Xc, sol.x)
        Xal = kabsch_align(X, Xc).reshape(Xc.shape)      # quotient out SE(3)
        disp.append((Xal - Xc).reshape(-1))
    return np.array(disp)


cases = [("a=30  (generic, rank 36)", corners(30.0)),
         ("a=60  (octahedron, rank 36)", corners(60.0)),
         ("a=90  (rank 35)", corners(90.0)),
         ("tetrahedron (rank 35)", build_tet_config(
             [s for s in tet_assignments() if is_tet(s)][0]))]

for name, Xc in cases:
    s = np.linalg.svd(jacobian(Xc), compute_uv=False)
    print(f"\n=== {name}   J-rank = {int(np.sum(s > 1e-9*max(1,s[0])))} ===")
    for eps in (1e-2, 1e-3):
        D = local_cloud(Xc, eps)
        sv = np.linalg.svd(D, compute_uv=False)
        sv = sv / sv[0]
        n6 = int(np.sum(sv > 1e-2))
        print(f"  eps={eps:.0e}  {len(D)} samples   normalised spectrum "
              f"{np.array2string(sv[:9], precision=4, suppress_small=True)}")
        print(f"            significant directions (>1% of the leading one): {n6}")
