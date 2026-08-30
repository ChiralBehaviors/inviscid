"""Step L (M4, M5): the competing V candidates, measured in the same harness.

M4  V_thomson = sum_{i<j} 1/r_ij over the 12 shared vertices.
    The Thomson minimum for 12 points on a sphere IS the icosahedron, which
    sits at a = 22.238756093 in this family, so this potential looks like it
    hands us the icosahedron as ground state for free. THE CAVEAT the memo
    records as unmeasured: the sphere RADIUS also shrinks 1.4142136 -> 1.0 along
    the path, and 1/r_ij scales as 1/R, so shape gain competes against radius
    loss exactly as it did for hull volume. Measured here both ways:
      - raw          : sum 1/r_ij on the actual vertices (shape AND scale)
      - radius-normalised : project every vertex to the unit sphere first, which
                       isolates SHAPE. On the symmetric path this is exactly
                       R * V_raw, but it is defined off the path too, where
                       cosphericity is not guaranteed.
    Independent non-trivial check available: the radius-normalised value at
    a_ico must equal the known N=12 Thomson minimum, 49.165253058.

M5  V_strut = sum over non-hinged strut pairs 1/d_ij, d = segment-segment
    distance from jb_g's clearance machinery.
    Struts CROSS (d = 0 exactly) throughout 60 < a < 120 and 240 < a < 300, and
    COINCIDE in pairs at 60/120/240/300. So 1/d is +infinity on a whole INTERVAL,
    not merely at four angles. That is reported raw. A regularised variant
    1/(d+eps) is reported alongside with eps STATED, and it is a modelling
    choice, not a fix -- the raw divergence is the honest content.
"""
import numpy as np

from jb_a_family import corners, L_EDGE
from jb_b_variety import PAIRS
from jb_g_strut_clearance import segment_distance, _hinged
from jb_j_internal_frame import Frame, inertia
from jb_k_hull_hessian import aligned_frame, hull_vol, V_TET

A_ICO = 22.238756093
THOMSON_12 = 49.165253058          # published N=12 minimum (icosahedron)


def thomson(P):
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    iu = np.triu_indices(len(P), 1)
    d = D[iu]
    if d.min() < 1e-12:
        return np.inf
    return float(np.sum(1.0 / d))


def thomson_normalised(P):
    """Radius-normalised about the ORIGIN -- the variant as the memo states it.

    CAVEAT MEASURED IN jb_m: this is NOT translation invariant. A rigid global
    translation changes it (49.313952 -> 57.216668 for a 0.76-length shift), so
    as written it is not a function on the linkage variety at all -- it would
    generate spurious forces on the centre of mass. It is retained because it is
    the stated candidate and its symmetric-path profile is what the memo asks
    about; use `thomson_normalised_centroid` for anything dynamical.
    """
    r = np.linalg.norm(P, axis=1)
    if r.min() < 1e-12:
        return np.inf
    return thomson(P / r[:, None])


def thomson_normalised_centroid(P):
    """Radius-normalised about the vertex CENTROID. Translation invariant.

    Identical to `thomson_normalised` on the symmetric path (the centroid sits
    at the origin there), so it changes no profile in this file -- but it is a
    different function OFF the path, which is exactly where the 6-D Hessian
    lives. Whether the inertia counts survive the swap is a check that can fail.
    """
    Q = P - P.mean(axis=0)
    r = np.linalg.norm(Q, axis=1)
    if r.min() < 1e-12:
        return np.inf
    return thomson(Q / r[:, None])


def _strut_pairs():
    out = []
    for m in range(24):
        i, j = divmod(m, 3)
        for n in range(m + 1, 24):
            k, l = divmod(n, 3)
            ea = {(i, j), (i, (j + 1) % 3)}
            eb = {(k, l), (k, (l + 1) % 3)}
            if not _hinged(ea, eb):
                out.append(((i, j), (k, l)))
    return out


STRUT_PAIRS = _strut_pairs()


def strut_gaps(X):
    g = []
    for (i, j), (k, l) in STRUT_PAIRS:
        g.append(segment_distance(X[i, j], X[i, (j + 1) % 3],
                                  X[k, l], X[k, (l + 1) % 3]))
    return np.array(g)


def strut_repulsion(X, eps=0.0):
    g = strut_gaps(X)
    if eps == 0.0 and g.min() < 1e-12:
        return np.inf
    return float(np.sum(1.0 / (g + eps)))


def vert(a):
    X = corners(a)
    return np.array([X[i, j] for (i, j), _ in PAIRS])


if __name__ == "__main__":
    np.set_printoptions(precision=7, suppress=False, linewidth=170)

    print("=== check the harness against a number it did not choose ===")
    print(f"  radius-normalised Thomson at a_ico = "
          f"{thomson_normalised(vert(A_ICO)):.9f}")
    print(f"  published N=12 Thomson minimum     = {THOMSON_12:.9f}")
    print(f"  |diff| = {abs(thomson_normalised(vert(A_ICO)) - THOMSON_12):.2e}"
          "   (this check CAN fail: a wrong vertex set misses it by ~0.4)")
    print(f"  cuboctahedron (a=0) normalised     = "
          f"{thomson_normalised(vert(0.0)):.9f}   (must be LARGER)")

    print("\n=== M4/M5 profiles on the symmetric path (DECISION 16: no region excluded) ===")
    print("      a       Vol_hull/V_tet   V_thomson(raw)  V_thomson(norm)   R"
          "        V_strut(raw)   V_strut(eps=0.01)")
    grid = sorted(set(list(np.arange(0, 181, 10.0)) +
                      [7.218951982, A_ICO, 59.0, 60.0, 61.0, 90.0, 119.0,
                       120.0, 121.0, 145.0, 175.0]))
    for a in grid:
        X = corners(a)
        P = vert(a)
        R = np.linalg.norm(P, axis=1).mean()
        tr = thomson(P)
        tn = thomson_normalised(P)
        sr = strut_repulsion(X)
        se = strut_repulsion(X, eps=0.01)
        try:
            hv = hull_vol(P)
        except Exception:
            hv = float("nan")
        print(f"  {a:8.4f}   {hv:12.6f}   {tr:14.6f}  {tn:14.6f}  {R:8.6f}  "
              f"{sr:14.4f}  {se:14.4f}")

    print("\n=== M4: where is each Thomson variant minimised on the path? ===")
    from scipy.optimize import minimize_scalar
    for name, fn in (("raw", lambda a: thomson(vert(a))),
                     ("radius-normalised", lambda a: thomson_normalised(vert(a)))):
        best = min(((fn(a), a) for a in np.linspace(-59.5, 59.5, 2387)),
                   key=lambda t: t[0])
        r = minimize_scalar(fn, bracket=(best[1] - 1.0, best[1], best[1] + 1.0),
                            method="brent", options={"xtol": 1e-12})
        print(f"  {name:18s} coarse min at a={best[1]:.4f}  refined a={r.x:.9f}"
              f"   V={r.fun:.9f}")
        print(f"                     value at a=0     {fn(0.0):.9f}")
        print(f"                     value at a_ico   {fn(A_ICO):.9f}")

    print("\n  centroid-referenced variant agrees with the origin-referenced one")
    print("  ON the symmetric path (centroid at the origin) -- check, then use it:")
    for a in (0.0, A_ICO, 30.0, 90.0):
        print(f"    a={a:9.5f}  origin-ref {thomson_normalised(vert(a)):.9f}   "
              f"centroid-ref {thomson_normalised_centroid(vert(a)):.9f}   "
              f"|centroid| = {np.linalg.norm(vert(a).mean(axis=0)):.2e}")

    print("\n=== M4: gradient and Hessian inertia in all 6 internal directions ===")
    for a0, tag in ((0.0, "VE"), (A_ICO, "icosahedron"),
                    (7.218951982, "hull max")):
        F = aligned_frame(a0)
        for name, fn in (("thomson raw", thomson),
                         ("thomson norm(O)", thomson_normalised),
                         ("thomson norm(c)", thomson_normalised_centroid)):
            g = F.grad(fn, h=1e-3)
            H = F.hessian(fn, h=1e-3)
            (p, z, n), ev = inertia(H)
            print(f"  a0={a0:11.7f} [{tag:12s}] {name:13s} "
                  f"|grad|={np.linalg.norm(g):.3e} transv={np.linalg.norm(g[1:]):.3e}")
            print(f"      eig(Hess V) = {np.array2string(ev, precision=6)}")
            print(f"      inertia (pos,zero,neg) = ({p},{z},{n})   "
                  f"{'MINIMUM in all 6' if n == 0 and z == 0 else 'not a minimum'}")

    print("\n=== M4: are the Thomson candidates SMOOTH where the hull is not? ===")
    print("  same discriminator as jb_k: |D_+ - D_-| falling with h = smooth,")
    print("  flat in h = corner. Run at a=0, where Vol_hull has its cone point.")
    F0 = aligned_frame(0.0)
    for name, fn in (("thomson raw", thomson),
                     ("thomson norm(c)", thomson_normalised_centroid),
                     ("Vol_hull (control)", hull_vol)):
        line = []
        for h in (1e-2, 1e-3, 1e-4):
            g, gp, gm = F0.grad(fn, h=h, onesided=True)
            line.append(np.abs(gp - gm).max())
        print(f"  {name:20s} max|D+ - D-| at h=1e-2,1e-3,1e-4: "
              f"{line[0]:.6f} {line[1]:.6f} {line[2]:.6f}   "
              f"{'SMOOTH' if line[2] < 0.2 * line[1] else 'CORNER'}")

    print("\n=== M5: the strut-repulsion divergence, stated rather than smoothed ===")
    print("  minimum non-hinged strut gap, and whether 1/d is finite:")
    for a in (0, 30, 55, 59, 59.9, 60, 65, 90, 110, 119.9, 120, 121, 130, 180):
        X = corners(a)
        g = strut_gaps(X)
        print(f"  a={a:7.2f}  min gap {g.min():.9f}  "
              f"V_raw = {strut_repulsion(X)!r:>10}   "
              f"V_eps0.01 = {strut_repulsion(X, 0.01):12.4f}")
    print("  V_strut is +inf on the WHOLE interval (60,120), not only at 60/120.")
    print("  Consequence: with 1/d and no regularisation, strut repulsion re-imposes")
    print("  as an infinite barrier exactly the region DECISION 16 declared legal.")

    print("\n  where is V_strut minimised on the finite part [-59.9, 59.9]?")
    best = min(((strut_repulsion(corners(a)), a)
                for a in np.linspace(-59.9, 59.9, 1199)), key=lambda t: t[0])
    print(f"    coarse min at a={best[1]:.4f}, V={best[0]:.6f}")
    r = minimize_scalar(lambda a: strut_repulsion(corners(a)),
                        bracket=(-5.0, 0.0, 5.0), method="brent",
                        options={"xtol": 1e-10})
    print(f"    refined near 0: a={r.x:.9f}  V={r.fun:.9f}")
    print(f"    V_strut(0) = {strut_repulsion(corners(0.0)):.9f}")

    print("\n  M5 in the full variety at a=0 (is the VE stationary for V_strut?)")
    F = aligned_frame(0.0)
    g = F.grad(lambda P: 0.0, h=1e-3)   # placeholder shape; strut V needs X not P
    F0 = Frame(0.0)

    def strut_of_s(s):
        y, res = F0.solve(s)
        assert res < 1e-9
        return strut_repulsion(F0.config(y))

    n = 6
    f0 = strut_of_s(np.zeros(n))
    for h in (1e-3, 3e-4):
        gr = np.zeros(n)
        H = np.zeros(n)
        for i in range(n):
            e = np.zeros(n)
            e[i] = h
            fp, fm = strut_of_s(e), strut_of_s(-e)
            gr[i] = (fp - fm) / (2 * h)
            H[i] = (fp - 2 * f0 + fm) / h ** 2
        print(f"    h={h:.0e}  V_strut(0)={f0:.6f}  grad={np.array2string(gr, precision=5)}")
        print(f"             diagonal curvature = {np.array2string(H, precision=5)}")
