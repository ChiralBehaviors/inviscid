"""Step B: the jitterbug as a LINKAGE, not a 1-parameter family.

Eight rigid triangles, free in space: 8 x 6 = 48 DOF.
Twelve shared vertices, each identifying one corner of one triangle with one
corner of another: 12 x 3 = 36 scalar constraints.
Global rigid motions: 6, always in the null space and always discarded.

The symmetric jitterbug path is a 1-D curve inside whatever variety that leaves.
This computes the constraint-Jacobian rank ALONG the path, which is what locates
the singularity and says how many branches can leave it.

Local parameterisation at a configuration: body i gets an infinitesimal rotation
dw_i about its own current centroid plus a translation dt_i, so for corner j

    d p_ij = dw_i x r_ij + dt_i,   r_ij = X[i][j] - centroid_i
"""
import numpy as np
from analysis.model.jitterbug import corners, cluster, L_EDGE


def pairing(a_probe=30.0, tol=1e-7):
    """The linkage's fixed combinatorics: 12 pairs of (face, corner).

    Read off at a generic a, then held FIXED — at a=60 the 12 shared vertices
    merge to 6 but the hinge incidences do not change.
    """
    _, mult, labels = cluster(corners(a_probe), tol=tol)
    assert len(mult) == 12 and set(mult.tolist()) == {2}, "not a 12x2 configuration"
    lab = labels.reshape(8, 3)
    pairs = []
    for v in range(12):
        where = [(i, j) for i in range(8) for j in range(3) if lab[i, j] == v]
        assert len(where) == 2
        pairs.append(tuple(where))
    return pairs


PAIRS = pairing()


def jacobian(X):
    """36 x 48 constraint Jacobian at corner configuration X (8,3,3)."""
    cen = X.mean(axis=1)                      # (8,3) triangle centroids
    J = np.zeros((36, 48))
    for c, ((i, j), (k, l)) in enumerate(PAIRS):
        rij, rkl = X[i, j] - cen[i], X[k, l] - cen[k]
        for d in range(3):
            row = 3 * c + d
            # d p_d / d w  where  d p = dw x r  =>  cross_row(r, d) . dw
            J[row, 3 * i:3 * i + 3] += cross_row(rij, d)
            J[row, 3 * k:3 * k + 3] -= cross_row(rkl, d)
            J[row, 24 + 3 * i + d] += 1.0
            J[row, 24 + 3 * k + d] -= 1.0
    return J


def cross_row(r, d):
    """Row d of the matrix M with M @ w = w x r  (i.e. -[r]_x)."""
    M = np.array([[0, r[2], -r[1]],
                  [-r[2], 0, r[0]],
                  [r[1], -r[0], 0]])
    return M[d]


def analyse(a, tol=1e-9):
    X = corners(a)
    J = jacobian(X)
    s = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(s > tol * max(1.0, s[0])))
    return rank, 48 - rank, s


def path_tangent(a, h=1e-6):
    """d/da of the configuration, expressed in the (dw, dt) coordinates."""
    Xp, Xm = corners(a + h), corners(a - h)
    dX = (Xp - Xm) / (2 * h)
    cen = corners(a).mean(axis=1)
    X = corners(a)
    xi = np.zeros(48)
    for i in range(8):
        # least-squares fit of dX[i] = dw x r + dt over the 3 corners
        A = np.zeros((9, 6))
        b = dX[i].reshape(9)
        for j in range(3):
            r = X[i, j] - cen[i]
            for d in range(3):
                A[3 * j + d, 0:3] = cross_row(r, d)
                A[3 * j + d, 3 + d] = 1.0
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
        resid = np.linalg.norm(A @ sol - b)
        assert resid < 1e-6, f"face {i} motion is not rigid: residual {resid:.2e}"
        xi[3 * i:3 * i + 3] = sol[0:3]
        xi[24 + 3 * i:24 + 3 * i + 3] = sol[3:6]
    return xi


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    print("hinge incidences (12 shared vertices, held fixed for all a):")
    for v, ((i, j), (k, l)) in enumerate(PAIRS):
        print(f"   v{v:2d}: face {i} corner {j}   ==   face {k} corner {l}")

    print("\n--- constraint-Jacobian rank along the symmetric path ---")
    print("  (48 unknowns, 36 constraints; internal DOF = nullity - 6 global)")
    for a in (0, 5, 15, 22.238756093, 30, 45, 55, 59, 59.9, 59.99, 60,
              60.01, 60.1, 61, 75, 90, 120, 180):
        rank, nullity, s = analyse(a)
        print(f"  a={a:10.5f}  rank={rank:3d}  nullity={nullity:3d}  "
              f"internal DOF={nullity - 6:3d}   sigma_min={s[-1]:.3e} "
              f"sigma_35={s[34]:.3e} sigma_36={s[35]:.3e}")

    print("\n--- is the symmetric path tangent actually in the null space? ---")
    for a in (0, 22.238756093, 30, 59, 61, 90):
        xi = path_tangent(a)
        J = jacobian(corners(a))
        print(f"  a={a:10.5f}  ||J xi|| = {np.linalg.norm(J @ xi):.3e}   "
              f"||xi|| = {np.linalg.norm(xi):.6f}")
