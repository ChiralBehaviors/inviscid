"""REVIEW CHECK 8: 'rank 35 => genuine singularity' is asserted at a=90 and at
the tetrahedron.  A Jacobian rank drop makes the constraint map fail to be a
submersion; it does NOT by itself make the variety singular.  The discriminating
test is whether the EXTRA null direction integrates to a finite motion:

  extends  -> extra branch / higher-dimensional stratum (a real branch point)
  dies     -> infinitesimal flex only; the extra nullity is first-order noise
"""
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners
from jb_b_variety import jacobian
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import build_tet_config
from jb_e_tighten import is_tet
from rv_4_continuation import shape_dist


def globals_at(X):
    cen = X.mean(axis=1)
    g = np.zeros((6, 48))
    for d in range(3):
        for i in range(8):
            g[d, 24 + 3*i + d] = 1.0
    for d in range(3):
        e = np.eye(3)[d]
        for i in range(8):
            g[3 + d, 3*i:3*i+3] = e
            g[3 + d, 24 + 3*i:24 + 3*i+3] = np.cross(e, cen[i])
    return g


def internal_modes(X):
    J = jacobian(X)
    s = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(s > 1e-9 * max(1.0, s[0])))
    _, _, Vt = np.linalg.svd(J)
    null = Vt[rank:]
    Q, _ = np.linalg.qr(globals_at(X).T)
    internal = null - (null @ Q) @ Q.T
    Ui, Si, _ = np.linalg.svd(internal.T, full_matrices=False)
    return rank, Ui[:, Si > 1e-8].T


def march(X0, xi, step=0.01, n=40):
    """Fixed-direction predictor + exact corrector; quotient arclength."""
    z = np.zeros(48)
    d = xi / np.linalg.norm(xi)
    arc, done = 0.0, 0
    for _ in range(n):
        sol = least_squares(lambda zz: residual(place(X0, zz)), z + step * d,
                            xtol=1e-14, ftol=1e-14, gtol=1e-14)
        if np.linalg.norm(residual(place(X0, sol.x))) > 1e-10:
            break
        arc += shape_dist(place(X0, sol.x), place(X0, z))
        z = sol.x
        done += 1
    return arc, done, place(X0, z)


for name, X in (("a=90 (faces through centre)", corners(90.0)),
                ("a=270", corners(270.0)),
                ("tetrahedron", build_tet_config(
                    [s for s in tet_assignments() if is_tet(s)][0])),
                ("a=30 (control, generic)", corners(30.0))):
    rank, basis = internal_modes(X)
    print(f"\n=== {name}: rank {rank}, {len(basis)} internal modes ===")
    arcs = []
    for m, xi in enumerate(basis):
        arc, steps, Xend = march(X, xi, 0.01, 40)
        arcs.append(arc)
        print(f"  mode {m}: quotient arclength = {arc:.4f}  steps {steps}/40"
              f"   {'FINITE' if arc > 0.1 else 'infinitesimal only'}")
    print(f"  -> modes that integrate to finite motion: "
          f"{sum(1 for a in arcs if a > 0.1)} of {len(arcs)}")
