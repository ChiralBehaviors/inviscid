"""Step C: (1) is a=60 really a branch point? (2) do the 6 DOF extend to finite
motions? (3) does the TETRAHEDRON configuration exist in the linkage?
"""
import itertools
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners, L_EDGE
from jb_b_variety import PAIRS, jacobian, analyse, path_tangent


# ---------------------------------------------------------------- 1. rank scan
def rank_scan(lo, hi, n=401):
    """Smallest singular value of J across a dense sweep."""
    out = []
    for a in np.linspace(lo, hi, n):
        s = np.linalg.svd(jacobian(corners(a)), compute_uv=False)
        out.append((a, s[35], s[-1]))
    return np.array(out)


# ------------------------------------------------- 2. finite-motion continuation
def unpack(z):
    """z (48,) -> rotations (8,3,3) and translations (8,3) via exp of rot vectors."""
    Rs, ts = [], []
    for i in range(8):
        w = z[3 * i:3 * i + 3]
        th = np.linalg.norm(w)
        if th < 1e-14:
            Rs.append(np.eye(3))
        else:
            k = w / th
            K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
            Rs.append(np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K))
        ts.append(z[24 + 3 * i:24 + 3 * i + 3])
    return np.array(Rs), np.array(ts)


def place(X0, z):
    """Apply body motions z to base configuration X0, each about its centroid."""
    Rs, ts = unpack(z)
    cen = X0.mean(axis=1)
    return np.array([(Rs[i] @ (X0[i] - cen[i]).T).T + cen[i] + ts[i] for i in range(8)])


def residual(X):
    return np.array([X[i, j] - X[k, l] for (i, j), (k, l) in PAIRS]).reshape(-1)


def continue_along(a0, xi, step, n_steps=40):
    """March along null direction xi from the config at a0, Newton-correcting
    back onto the constraint set at each step. Returns the achieved arclength."""
    X0 = corners(a0)
    z = np.zeros(48)
    total = 0.0
    for _ in range(n_steps):
        z_pred = z + step * xi / np.linalg.norm(xi)
        sol = least_squares(lambda zz: residual(place(X0, zz)), z_pred,
                            xtol=1e-14, ftol=1e-14, gtol=1e-14)
        if np.linalg.norm(residual(place(X0, sol.x))) > 1e-10:
            break
        total += np.linalg.norm(place(X0, sol.x) - place(X0, z))
        z = sol.x
    return total, place(X0, z)


# ------------------------------------------------------- 3. the tetrahedron test
def tet_assignments():
    """Every way to send each triangle's 3 corners to 3 distinct tetrahedron
    vertices such that all 12 hinge pairings agree. An equilateral triangle of
    edge L maps rigidly onto a tet face of edge L under ANY of the 6 bijections
    (flipping over is a rigid motion in 3D), so combinatorial consistency is
    exactly realisability here."""
    options = [dict(zip(range(3), p))
               for p in itertools.permutations(range(4), 3)]   # 24 per face
    sols, assign = [], [None] * 8

    def ok(i, opt):
        for (p, j), (q, l) in PAIRS:
            if p == i and assign[q] is not None and opt[j] != assign[q][l]:
                return False
            if q == i and assign[p] is not None and opt[l] != assign[p][j]:
                return False
        return True

    def rec(i):
        if i == 8:
            sols.append([dict(d) for d in assign])
            return
        for opt in options:
            if ok(i, opt):
                assign[i] = opt
                rec(i + 1)
                assign[i] = None

    rec(0)
    return sols


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)

    print("=== 1. IS a=60 A SINGULARITY? dense scan of sigma_36 (the rank-deciding one) ===")
    for lo, hi in ((50.0, 70.0), (85.0, 95.0)):
        sc = rank_scan(lo, hi, 401)
        k = np.argmin(sc[:, 1])
        print(f"  sweep [{lo},{hi}]: min sigma_36 = {sc[k,1]:.6e} at a = {sc[k,0]:.4f}")
    for a in (59.999999, 60.0, 60.000001):
        r, n, s = analyse(a)
        print(f"  a={a:.6f}: rank={r} sigma_36={s[35]:.6e}")
    for a in (89.999999, 90.0, 90.000001):
        r, n, s = analyse(a)
        print(f"  a={a:.6f}: rank={r} sigma_36={s[35]:.6e}")

    print("\n=== 2. DO THE 6 INTERNAL DOF EXTEND TO FINITE MOTIONS? ===")
    a0 = 30.0
    J = jacobian(corners(a0))
    U, S, Vt = np.linalg.svd(J)
    null = Vt[36:]                      # 12 null directions (6 global + 6 internal)
    # project out the global rigid motions
    X0 = corners(a0)
    cen = X0.mean(axis=1)
    glob = np.zeros((6, 48))
    for d in range(3):
        for i in range(8):
            glob[d, 24 + 3 * i + d] = 1.0                    # translations
    for d in range(3):
        e = np.eye(3)[d]
        for i in range(8):
            glob[3 + d, 3 * i:3 * i + 3] = e                 # rotations
            glob[3 + d, 24 + 3 * i:24 + 3 * i + 3] = np.cross(e, cen[i])
    gerr = max(np.linalg.norm(J @ glob[k]) for k in range(6))
    print(f"  global motions satisfy J: max ||J g|| = {gerr:.2e}")
    assert gerr < 1e-9, "global rigid motions are not in the null space - convention bug"
    Q, _ = np.linalg.qr(glob.T)
    internal = null - (null @ Q) @ Q.T
    Ui, Si, _ = np.linalg.svd(internal.T, full_matrices=False)
    internal_basis = Ui[:, Si > 1e-8].T
    print(f"  internal null directions after removing globals: {len(internal_basis)}")
    tangent = path_tangent(a0)
    tproj = tangent - Q @ (Q.T @ tangent)
    ov = internal_basis @ tproj / np.linalg.norm(tproj)
    print(f"  symmetric-path tangent lies in that space: ||proj|| = {np.linalg.norm(ov):.9f}")
    for m, xi in enumerate(internal_basis):
        arc, _ = continue_along(a0, xi, 0.02, 40)
        print(f"  internal mode {m}: continued arclength = {arc:.4f} "
              f"({'FINITE MOTION' if arc > 0.5 else 'first-order flex only'})")

    print("\n=== 3. DOES THE TETRAHEDRON CONFIGURATION EXIST? ===")
    sols = tet_assignments()
    print(f"  hinge-consistent corner->tet-vertex assignments (no occupancy rule): {len(sols)}")
    profile = {}
    for s in sols:
        load = {}
        for d in s:
            load[tuple(sorted(d.values()))] = load.get(tuple(sorted(d.values())), 0) + 1
        profile.setdefault(tuple(sorted(load.values(), reverse=True)), []).append(s)
    print("  occupancy profiles (triangles per distinct tet face used):")
    for k in sorted(profile, key=lambda t: (-len(t), t)):
        print(f"    {str(k):18s} -> {len(profile[k]):6d} assignments"
              f"   {'<== THE TETRAHEDRON (2 on each of 4 faces)' if k == (2,2,2,2) else ''}")
    tets = profile.get((2, 2, 2, 2), [])
    print(f"\n  TETRAHEDRON configurations: {len(tets)}")
    if tets:
        s = tets[0]
        print(f"  example assignment (face -> corner:tetvertex):")
        for i, d in enumerate(s):
            print(f"    face {i}: {dict(sorted(d.items()))}  on tet face {tuple(sorted(d.values()))}")
