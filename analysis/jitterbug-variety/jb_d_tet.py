"""Step D: construct the tetrahedron configuration explicitly, verify it against
the hinge constraints, and test whether it is REACHABLE from the vector
equilibrium inside the linkage's configuration variety.

Existence is not reachability. Step C showed 1344 hinge-consistent ways to put
2 triangles on each face of a tetrahedron. This checks (a) that such a placement
really satisfies the constraints as rigid bodies, (b) what the Jacobian does
there, and (c) whether a continuous constraint-respecting path joins it to a=0.
"""
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners, L_EDGE
from jb_b_variety import PAIRS, jacobian
from jb_c_branches import tet_assignments, place, residual

TET = np.array([[0.5, 0.5, 0.5], [0.5, -0.5, -0.5],
                [-0.5, 0.5, -0.5], [-0.5, -0.5, 0.5]])   # regular, edge sqrt(2)


def build_tet_config(assign):
    return np.array([[TET[assign[i][j]] for j in range(3)] for i in range(8)])


def kabsch(P, Q):
    """Rigid transform aligning P onto Q (both (N,3)); returns rotated P."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (R @ Pc.T).T + Q.mean(0)


def homotopy_to_target(X0, Xt, n=60, pull=1.0):
    """March from X0 to Xt staying on the constraint set.

    At each step the target is the current config moved a fraction of the way
    toward Xt; the solve keeps the hinge residual at zero while getting as close
    to that waypoint as it can. Records the worst constraint violation seen."""
    z = np.zeros(48)
    worst_r = 0.0
    for s in np.linspace(1.0 / n, 1.0, n):
        Xa = kabsch(Xt.reshape(-1, 3), place(X0, z).reshape(-1, 3)).reshape(8, 3, 3)
        way = (1 - s) * place(X0, z) + s * Xa

        def f(zz):
            X = place(X0, zz)
            return np.concatenate([1e3 * residual(X),
                                   pull * (X - way).reshape(-1)])

        sol = least_squares(f, z, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        z = sol.x
        worst_r = max(worst_r, np.linalg.norm(residual(place(X0, z))))
    X = place(X0, z)
    Xa = kabsch(Xt.reshape(-1, 3), X.reshape(-1, 3)).reshape(8, 3, 3)
    return X, np.linalg.norm(X - Xa) / np.sqrt(24), worst_r


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    def is_tet(s):
        faces = [tuple(sorted(d.values())) for d in s]
        return len(set(faces)) == 4 and all(faces.count(f) == 2 for f in faces)

    tets = [s for s in tet_assignments() if is_tet(s)]
    print(f"tetrahedron assignments available: {len(tets)}")

    print("\n=== (a) DOES THE PLACEMENT SATISFY THE HINGE CONSTRAINTS? ===")
    for t in range(3):
        Xt = build_tet_config(tets[t])
        e = np.array([np.linalg.norm(Xt[i, j] - Xt[i, (j + 1) % 3])
                      for i in range(8) for j in range(3)])
        print(f"  candidate {t}: hinge residual = {np.linalg.norm(residual(Xt)):.3e}   "
              f"strut len {e.min():.12f}..{e.max():.12f}  (need {L_EDGE:.12f})")

    Xt = build_tet_config(tets[0])
    print("\n=== (b) THE JACOBIAN AT THE TETRAHEDRON ===")
    s = np.linalg.svd(jacobian(Xt), compute_uv=False)
    rank = int(np.sum(s > 1e-9 * s[0]))
    print(f"  rank = {rank}   nullity = {48-rank}   internal DOF = {48-rank-6}"
          f"   (generic point of the variety: rank 36, internal DOF 6)")
    print(f"  smallest 6 singular values: {s[-6:]}")

    print("\n=== (c) IS IT REACHABLE FROM THE VECTOR EQUILIBRIUM? ===")
    X0 = corners(0.0)
    print(f"  start: a=0 vector equilibrium, hinge residual {np.linalg.norm(residual(X0)):.2e}")
    for t in range(3):
        Xt = build_tet_config(tets[t])
        Xf, dist, worst = homotopy_to_target(X0, Xt, n=60)
        print(f"  target {t}: final RMS distance to tetrahedron = {dist:.6e}   "
              f"worst hinge residual along the path = {worst:.2e}   "
              f"{'REACHED' if dist < 1e-6 else 'NOT reached'}")
