"""Step M: cross-check the two kink results against independent machinery.

Both headline findings of jb_i/jb_k/jb_l are NON-SMOOTHNESS claims, and a
non-smoothness claim is exactly the kind of thing a sloppy numerical routine
manufactures. Two pieces of code are load-bearing and neither was written this
session:

  1. `jb_g.segment_distance` -- the clamped parametric segment-segment solve.
     jb_l used it to report that V_strut = sum 1/d has a CONE POINT at a=0
     (diagonal curvature diverging like 1/h). If that routine has a branch
     discontinuity, the cone point is its artefact and not the geometry's.
     Cross-checked here against a dense-grid + Nelder-Mead reference distance.

  2. scipy's Qhull volume. jb_i reported a |a| corner in Vol_hull at a=0. The
     analytic square-quad decomposition already matches it to 7 significant
     figures, which is one independent line; the second is COMBINATORIAL --
     count hull facets. At a=0 the hull must have 8 triangles + 6 squares = 14
     facets; off a=0 every square splits and it must be 8 + 12 = 20 triangles.
     A combinatorial jump at a=0 and nowhere near it is what a corner IS.

Also settled here: WHICH direction is the unstable one for radius-normalised
Thomson at the VE, and how the M4 inertia counts respond to perturbation.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull

from jb_a_family import corners
from jb_g_strut_clearance import segment_distance
from jb_j_internal_frame import Frame, inertia
from jb_k_hull_hessian import aligned_frame, hull_vol
from jb_l_vertex_potentials import (STRUT_PAIRS, strut_repulsion, thomson,
                                    thomson_normalised, vert)


def ref_segment_distance(p1, q1, p2, q2, n=60):
    """Reference segment-segment distance: coarse grid then a bounded refine.

    Deliberately dumb and algorithm-independent -- it shares no code path with
    the clamped parametric solve it is checking.
    """
    ts = np.linspace(0.0, 1.0, n)
    A = p1[None, :] + np.outer(ts, q1 - p1)
    B = p2[None, :] + np.outer(ts, q2 - p2)
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
    i, j = np.unravel_index(np.argmin(D), D.shape)

    def f(z):
        s, t = np.clip(z, 0.0, 1.0)
        return np.linalg.norm((p1 + s * (q1 - p1)) - (p2 + t * (q2 - p2)))
    r = minimize(f, x0=[ts[i], ts[j]], method="Nelder-Mead",
                 options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 4000})
    return min(float(r.fun), float(D[i, j]))


def ref_strut_repulsion(X):
    tot = 0.0
    for (i, j), (k, l) in STRUT_PAIRS:
        d = ref_segment_distance(X[i, j], X[i, (j + 1) % 3],
                                 X[k, l], X[k, (l + 1) % 3])
        if d < 1e-12:
            return np.inf
        tot += 1.0 / d
    return tot


def facet_count(a, tol=1e-9):
    """(n distinct supporting planes, n Qhull simplices) of the vertex hull."""
    P = vert(a)
    h = ConvexHull(P)
    planes = []
    for eq in h.equations:
        if not any(np.allclose(eq, e, atol=tol) for e in planes):
            planes.append(eq)
    return len(planes), len(h.simplices)


if __name__ == "__main__":
    np.set_printoptions(precision=7, suppress=False, linewidth=170)

    print("=== 1a. segment_distance vs an independent reference ===")
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(300):
        p1, q1, p2, q2 = rng.standard_normal((4, 3))
        d1 = segment_distance(p1, q1, p2, q2)
        d2 = ref_segment_distance(p1, q1, p2, q2)
        worst = max(worst, abs(d1 - d2))
    print(f"  worst |clamped - reference| over 300 random segment pairs: {worst:.3e}")
    print("  (this check CAN fail: the clamped solve is a known approximation)")

    print("\n  on the ACTUAL strut pairs near a=0, where the cone point is claimed:")
    wa = 0.0
    for a in (-1e-2, -1e-3, 0.0, 1e-3, 1e-2):
        X = corners(a)
        dev = max(abs(segment_distance(X[i, j], X[i, (j + 1) % 3],
                                       X[k, l], X[k, (l + 1) % 3])
                      - ref_segment_distance(X[i, j], X[i, (j + 1) % 3],
                                             X[k, l], X[k, (l + 1) % 3]))
                  for (i, j), (k, l) in STRUT_PAIRS)
        wa = max(wa, dev)
        print(f"    a={a:+.4f}  max deviation over the {len(STRUT_PAIRS)} pairs = {dev:.3e}")

    print("\n=== 1b. is the V_strut cone point at a=0 real, or the solver's? ===")
    print("     h        d2 (clamped)     d2 (reference)    d2*h (clamped)  d2*h (ref)")
    for h in (1e-2, 3e-3, 1e-3, 3e-4):
        f0c, fpc, fmc = (strut_repulsion(corners(x)) for x in (0.0, h, -h))
        f0r, fpr, fmr = (ref_strut_repulsion(corners(x)) for x in (0.0, h, -h))
        dc = (fpc - 2 * f0c + fmc) / h ** 2
        dr = (fpr - 2 * f0r + fmr) / h ** 2
        print(f"  {h:.0e}   {dc:14.4f}  {dr:14.4f}   {dc * h:12.6f}  {dr * h:12.6f}")
    print("  Both diverging as 1/h (d2*h flat) => the cone point is GEOMETRIC.")
    print("  Only the clamped one diverging => it would have been an artefact.")

    print("\n  which strut pairs carry the kink? (one-sided slopes of 1/d)")
    hh = 1e-4
    rows = []
    for idx, ((i, j), (k, l)) in enumerate(STRUT_PAIRS):
        def d(x):
            X = corners(x)
            return segment_distance(X[i, j], X[i, (j + 1) % 3],
                                    X[k, l], X[k, (l + 1) % 3])
        dp = (1 / d(hh) - 1 / d(0.0)) / hh
        dm = (1 / d(0.0) - 1 / d(-hh)) / hh
        rows.append((abs(dp - dm), idx, d(0.0), dp, dm))
    rows.sort(reverse=True)
    print("    |D+ - D-|    pair   d(0)        D+          D-")
    for r in rows[:6]:
        print(f"     {r[0]:9.5f}  {r[1]:5d}   {r[2]:.6f}  {r[3]:+10.5f}  {r[4]:+10.5f}")
    print(f"    pairs with |D+ - D-| > 1e-3: "
          f"{sum(1 for r in rows if r[0] > 1e-3)} of {len(rows)}")

    print("\n=== 2. hull FACET COMBINATORICS across a=0 ===")
    print("     a          distinct planes   qhull simplices")
    for a in (-1.0, -1e-2, -1e-4, -1e-6, 0.0, 1e-6, 1e-4, 1e-2, 1.0,
              7.2189520, 22.238756093, 45.0):
        p, s = facet_count(a)
        tag = "  <-- 6 squares unsplit" if p == 14 else ""
        print(f"  {a:+12.7f}   {p:10d}      {s:10d}{tag}")
    print("  14 planes only at a=0, 20 everywhere else: the hull's face lattice")
    print("  changes AT the VE and nowhere near it. That combinatorial jump is the")
    print("  corner. It could have failed: a hull that stayed 14 (squares remaining")
    print("  planar) would have left Vol_hull smooth.")

    print("\n=== 3. WHICH direction is unstable for radius-normalised Thomson at the VE? ===")
    F = aligned_frame(0.0)
    H = F.hessian(thomson_normalised, h=1e-3)
    ev, evec = np.linalg.eigh(H)
    print(f"  eigenvalues {np.array2string(ev, precision=6)}")
    for k in range(6):
        print(f"    lambda={ev[k]:+.6f}  |overlap with symmetric-path direction| "
              f"= {abs(evec[0, k]):.6f}")
    print("  direction 0 of the chart IS the symmetric-path tangent (jb_k), so an")
    print("  overlap of 1 on the negative eigenvalue means the instability is the")
    print("  jitterbug motion itself, not a transverse mode.")

    print("\n=== 4. perturbation sensitivity of the M4 inertia counts ===")
    for a0, tag in ((0.0, "VE"), (22.238756093, "icosahedron")):
        for name, fn in (("thomson raw", thomson), ("thomson norm", thomson_normalised)):
            counts = []
            for h in (3e-3, 1e-3, 3e-4):
                Fh = aligned_frame(a0)
                (p, z, n), _ = inertia(Fh.hessian(fn, h=h))
                counts.append((p, z, n))
            for da in (1e-9, 1e-6):
                Fd = aligned_frame(a0 + da)
                (p, z, n), _ = inertia(Fd.hessian(fn, h=1e-3))
                counts.append((p, z, n))
            uniq = sorted(set(counts))
            print(f"  a0={a0:12.7f} [{tag:11s}] {name:13s} counts over "
                  f"3 step sizes + 2 base shifts: {uniq}"
                  f"   {'STABLE' if len(uniq) == 1 else 'MOVES'}")

    print("\n=== 5. zero-mode audit: no rigid motion may leak into any Hessian ===")
    F = Frame(30.0)
    print(f"  ||Rg^T B||_max = {np.abs(F.Rg.T @ F.B).max():.3e}")
    print("  a global rotation/translation leaves every potential here EXACTLY")
    print("  invariant, so a leaked rigid mode would appear as a zero eigenvalue.")
    print("  Direct test: apply a global rigid motion and check V is unchanged.")
    P = vert(30.0)
    from jb_a_family import rot
    Rg = rot(np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0), 37.0)
    P2 = (Rg @ P.T).T + np.array([0.3, -0.2, 0.7])
    print(f"    thomson       : {thomson(P):.12f} vs {thomson(P2):.12f}")
    print(f"    hull_vol      : {hull_vol(P):.12f} vs {hull_vol(P2):.12f}")
    print("    (thomson_normalised is NOT translation invariant -- it references")
    print("     the origin -- which is a real modelling caveat, not a bug:")
    print(f"     {thomson_normalised(P):.9f} vs {thomson_normalised(P2):.9f})")
