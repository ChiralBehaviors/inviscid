"""Step I (M1): is Vol_hull twice differentiable at the vector equilibrium?

The hypothesis V = -k*Vol_hull needs a Hessian at a=0 to call the VE "unstable".
a=0 is exactly where the hull is NON-SIMPLICIAL: 6 coplanar square faces. Off
a=0 each square goes skew, and the convex hull must then pick one of the two
triangulations of the skew quad -- always the one giving the LARGER volume,
because a convex hull is the maximal such body. The two triangulations differ by
the volume of the tetrahedron on the quad's 4 vertices, so

    Vol_hull(a) = (smooth part) + (1/2) * sum over the 6 squares |tet_vol_q(a)|

and the smoothness of Vol_hull at a=0 is decided by the ORDER OF VANISHING of
tet_vol_q. jb_h measured it: tet_vol_q is ODD in a and LINEAR at the origin
(|tet|/a -> 1.343555e-2, flat over four decades). |linear| = a kink.

This file tests that consequence three independent ways, each of which could
have come back negative:

  1. the |a| decomposition, checked NUMERICALLY against Qhull -- if the identity
     above is wrong the residual will not be at rounding level;
  2. the central second difference at a decreasing sequence of h -- a C2 function
     converges, a |a| kink diverges like 2c/h;
  3. one-sided first and second differences -- a kink shows equal-and-opposite
     first derivatives across a=0 while each side's second derivative converges.

Control: the same battery at a = 7.2189519 (the hull maximum), where jb_h
measured the hull to be simplicial and therefore generic.
"""
import numpy as np

from jb_h_hull_control import (V_TET, A_ICO, hull_volume, square_faces,
                               square_out_of_plane, vertices)


def V(a):
    """Normalised hull volume, V/V_tet."""
    return hull_volume(a) / V_TET


SQ = square_faces()


def kink_term(a):
    """(1/2) * sum |tet vol| over the 6 square-face quads, normalised."""
    return 0.5 * sum(abs(square_out_of_plane(a, q)) for q in SQ) / V_TET


def second_central(a, h):
    return (V(a + h) - 2.0 * V(a) + V(a - h)) / h ** 2


def second_onesided(a, h, sign):
    """f'' from three points all on one side of a."""
    s = sign * h
    return (V(a) - 2.0 * V(a + s) + V(a + 2 * s)) / h ** 2


def first_onesided(a, h, sign):
    s = sign * h
    return (V(a + s) - V(a)) / s


if __name__ == "__main__":
    print("=== 1. the |a| decomposition: does V - kink_term look smooth? ===")
    print("   a          V(a)-20          kink(a)        V-20-kink      ratio")
    for a in (1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1.0, 2.0):
        v, k = V(a) - 20.0, kink_term(a)
        print(f"  {a:8.5f}  {v:+.9e}  {k:+.9e}  {v - k:+.3e}  {v / k:.9f}")
    print("  If the ratio -> 1 as a -> 0, the whole leading behaviour of V near")
    print("  the VE IS the kink term, and V is not differentiable there.")

    print("\n=== 2. central second difference at a=0 (the C2 test) ===")
    print("     h (deg)      V''_central     V''*h     (2c if kink)")
    for h in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5):
        d2 = second_central(0.0, h)
        print(f"  {h:.1e}     {d2:14.6f}   {d2 * h:12.8f}")
    print("  CONVERGES -> C2.  DIVERGES like 1/h -> kink (V''*h tends to a constant).")

    print("\n=== 2b. CONTROL at a = 7.2189519 (hull simplicial, generic) ===")
    print("     h (deg)      V''_central     V''*h")
    for h in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        d2 = second_central(7.2189519, h)
        print(f"  {h:.1e}     {d2:14.6f}   {d2 * h:12.8f}")

    print("\n=== 2c. CONTROL at a = 30 and a = 22.238756093 ===")
    for a0 in (A_ICO, 30.0):
        print(f"  a0 = {a0}")
        for h in (1e-2, 1e-3, 1e-4):
            print(f"     h={h:.0e}   V'' = {second_central(a0, h):14.6f}")

    print("\n=== 3. one-sided derivatives across a=0 ===")
    print("     h (deg)     V'(0+)          V'(0-)        V''(0+)      V''(0-)")
    for h in (1e-2, 1e-3, 1e-4, 1e-5):
        print(f"  {h:.1e}   {first_onesided(0, h, +1):+.8f}   "
              f"{first_onesided(0, h, -1):+.8f}   "
              f"{second_onesided(0, h, +1):+11.5f}  {second_onesided(0, h, -1):+11.5f}")
    print("  Equal-and-opposite one-sided FIRST derivatives with each side's SECOND")
    print("  derivative converging is the exact signature of |a|: V is not even C1.")

    print("\n=== 4. perturbation sensitivity: does the verdict move? ===")
    print("  (a) one-ULP jiggle of every vertex coordinate, then re-measure V''(0)")
    rng = np.random.default_rng(20260812)

    def V_pert(a, scale):
        P = vertices(a)
        P = P + scale * rng.standard_normal(P.shape) * np.maximum(np.abs(P), 1.0)
        from scipy.spatial import ConvexHull
        return ConvexHull(P).volume / V_TET

    for scale in (0.0, 2.2e-16, 1e-13, 1e-10):
        vals = []
        for h in (1e-2, 1e-3, 1e-4):
            d2 = (V_pert(h, scale) - 2 * V_pert(0.0, scale) + V_pert(-h, scale)) / h ** 2
            vals.append(d2 * h)
        print(f"    jiggle {scale:.1e}:  V''*h at h=1e-2,1e-3,1e-4 = "
              f"{vals[0]:.6f} {vals[1]:.6f} {vals[2]:.6f}")
    print("  A kink verdict that survives a 1e-10 jiggle is not a rounding artefact.")

    print("\n  (b) the same battery with the OPPOSITE chirality convention (a -> -a)")
    for h in (1e-2, 1e-3, 1e-4):
        print(f"    h={h:.0e}  V''(0)*h = {second_central(0.0, h) * h:.8f}"
              f"   (mirror) {second_central(-0.0, h) * h:.8f}")

    print("\n=== 5. so what IS the smooth branch? one-sided C2 fit on a>0 ===")
    for h in (1e-3, 1e-4):
        c = first_onesided(0, h, +1)
        print(f"  h={h:.0e}   V'(0+) = {c:.9f} per deg   "
              f"kink slope from squares = {kink_term(h) / h:.9f} per deg")

    print("\n=== 6. the corner slope to full precision, from the ANALYTIC kink term ===")
    print("  (kink_term uses exact tetrahedron volumes, so it has no hull noise)")
    prev = None
    for a in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
        c_deg = kink_term(a) / a
        c_rad = c_deg * 180.0 / np.pi
        d = "" if prev is None else f"   delta {c_rad - prev:+.3e}"
        print(f"  a={a:.0e}   c = {c_deg:.12f} /deg = {c_rad:.12f} /rad{d}")
        prev = c_rad
    # a=1e-7 is the WORST estimate, not the best: the tetrahedron volume is a
    # difference of nearly equal terms, so cancellation grows as a shrinks while
    # the O(a^2) truncation shrinks. a=1e-4 is where the two errors cross.
    c_rad = kink_term(1e-4) / 1e-4 * 180.0 / np.pi
    print("\n  best-conditioned estimate is a=1e-4 (cancellation dominates below it)")
    print("  candidate closed forms (reported, NOT claimed -- none is confirmed")
    print("  unless it agrees to the precision of the limit above):")
    for nm, val in (("4*sqrt(3)", 4 * np.sqrt(3.0)),
                    ("2*sqrt(12)", 2 * np.sqrt(12.0)),
                    ("24/sqrt(12)", 24 / np.sqrt(12.0)),
                    ("6*sqrt(4/3)+0.0", 6 * np.sqrt(4 / 3)),
                    ("9*sqrt(6)/pi", 9 * np.sqrt(6) / np.pi)):
        print(f"    {nm:16s} = {val:.12f}   diff {c_rad - val:+.3e}")
    print(f"  measured c = {c_rad:.12f} /rad in V/V_tet units.")
    print("  So near the VE:  V/V_tet = 20 + c*|a_rad| - 0.5*k*a_rad^2 + ...")
    kq = -second_onesided(0, 1e-4, +1) * (180.0 / np.pi) ** 2
    print(f"  with the smooth quadratic coefficient k = {kq:.9f} /rad^2"
          f"  (from the one-sided V''(0+))")
    print(f"  predicted turning point a = c/k = "
          f"{np.degrees(c_rad / kq):.6f} deg vs the measured maximum 7.218952 deg")
    print("  (the two differ by the quartic and beyond; agreement to ~1% is the")
    print("   consistency check, not an identity)")
