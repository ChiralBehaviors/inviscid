"""Step E: harden the reachability result.

The step-D homotopy held the hinge residual near 1e-6, which is small but not
zero -- so it showed a path NEAR the variety, not ON it. Here every waypoint is
Newton-projected exactly back onto the constraint set, and the projection
displacement is recorded, so the path provably lies on the variety.
"""
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners
from jb_b_variety import jacobian
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import TET, build_tet_config, kabsch


def classify_endpoint(X, tol=1e-7):
    """LABELLING-FREE endpoint test.

    The RMS-to-target metric used below scores a walk against a FIXED corner
    correspondence, so a walk that lands on a genuine tetrahedron carrying a
    different labelling is scored a failure (its RMS is 1/sqrt(3) = 0.5774,
    the relabelling signature). This asks the only question that matters:
    are the 24 corners at 4 points with all pairwise distances equal?
    """
    pts = X.reshape(-1, 3)
    uniq = []
    for p in pts:
        if not any(np.linalg.norm(p - q) < tol for q in uniq):
            uniq.append(p)
    if len(uniq) != 4:
        return f"{len(uniq)} distinct corners - NOT a tetrahedron"
    D = [np.linalg.norm(uniq[i] - uniq[j]) for i in range(4) for j in range(i + 1, 4)]
    ok = max(D) - min(D) < 1e-7
    return (f"TETRAHEDRON (4 corners, edges {min(D):.9f}..{max(D):.9f})" if ok
            else f"4 corners but irregular ({min(D):.6f}..{max(D):.6f})")


def is_tet(s):
    faces = [tuple(sorted(d.values())) for d in s]
    return len(set(faces)) == 4 and all(faces.count(f) == 2 for f in faces)


def project(X0, z):
    """Newton-project onto {residual = 0}."""
    sol = least_squares(lambda zz: residual(place(X0, zz)), z,
                        xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    return sol.x, np.linalg.norm(residual(place(X0, sol.x)))


def walk(X0, Xt, n=80, pull=1.0, w=1e4):
    z = np.zeros(48)
    worst_after, worst_shift = 0.0, 0.0
    for s in np.linspace(1.0 / n, 1.0, n):
        Xa = kabsch(Xt.reshape(-1, 3), place(X0, z).reshape(-1, 3)).reshape(8, 3, 3)
        way = (1 - s) * place(X0, z) + s * Xa
        sol = least_squares(lambda zz: np.concatenate(
            [w * residual(place(X0, zz)), pull * (place(X0, zz) - way).reshape(-1)]),
            z, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        zp, r = project(X0, sol.x)
        worst_shift = max(worst_shift, np.linalg.norm(place(X0, zp) - place(X0, sol.x)))
        worst_after = max(worst_after, r)
        z = zp
    X = place(X0, z)
    Xa = kabsch(Xt.reshape(-1, 3), X.reshape(-1, 3)).reshape(8, 3, 3)
    return X, np.linalg.norm(X - Xa) / np.sqrt(24), worst_after, worst_shift


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    tets = [s for s in tet_assignments() if is_tet(s)]
    X0 = corners(0.0)

    print("=== reachability with EXACT projection onto the variety at every step ===")
    print("  (worst-after = hinge residual after projection; shift = how far projection moved it)")
    reached = 0
    for t in range(8):
        Xt = build_tet_config(tets[t])
        X, dist, wa, ws = walk(X0, Xt, n=80)
        hit = dist < 1e-6
        reached += hit
        print(f"  target {t}: RMS to tet = {dist:.3e}  worst-after = {wa:.2e}  "
              f"shift = {ws:.2e}   {'REACHED' if hit else 'not reached'}")
    print(f"\n  reached {reached}/8 sampled tetrahedron targets")

    print("\n=== the reached configuration really is a tetrahedron ===")
    Xt = build_tet_config(tets[1])
    X, dist, wa, ws = walk(X0, Xt, n=80)
    pts = X.reshape(-1, 3)
    uniq = []
    for p in pts:
        if not any(np.linalg.norm(p - q) < 1e-8 for q in uniq):
            uniq.append(p)
    D = [np.linalg.norm(uniq[i] - uniq[j]) for i in range(len(uniq)) for j in range(i + 1, len(uniq))]
    print(f"  distinct corner positions: {len(uniq)}  (tetrahedron => 4)")
    print(f"  pairwise distances: {np.min(D):.12f}..{np.max(D):.12f}  (regular => all equal)")
    print(f"  final hinge residual: {np.linalg.norm(residual(X)):.2e}")

    print("\n=== rank along the reached path vs the symmetric path ===")
    for name, XX in (("VE a=0", corners(0.0)), ("octa a=60", corners(60.0)),
                     ("faces-through-centre a=90", corners(90.0)),
                     ("tetrahedron (reached)", X)):
        s = np.linalg.svd(jacobian(XX), compute_uv=False)
        rank = int(np.sum(s > 1e-9 * s[0]))
        print(f"  {name:28s} rank={rank:3d}  internal DOF={48-rank-6:3d}  sigma_36={s[35]:.3e}")
