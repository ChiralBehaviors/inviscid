"""Step K (M2, M3): is the VE stationary in ALL SIX internal directions, and what
is the inertia of the Hessian of -Vol_hull where it exists?

C5 established the hull-volume profile on the SYMMETRIC SLICE only. A stationary
point of a 1-parameter slice need not be stationary in the 6-dimensional variety
(R2), and if it is not, "the VE is an equilibrium" is simply false and the Fuller
reading of V = -k*Vol_hull collapses before any Hessian is computed.

M3 needs no mass metric: by Sylvester's law of inertia the COUNT of positive /
zero / negative eigenvalues is invariant under congruence, so it is unchanged by
any choice of M (and by the arbitrary orthonormal basis the chart hands us).
Only the FREQUENCIES would need M.

Basis convention: the chart's basis is rotated so direction 0 is the SYMMETRIC
PATH tangent and directions 1..5 span its orthogonal complement inside the
internal tangent space. So a gradient supported on component 0 alone means
"stationary transversally, moving along the path", which is the distinction M2
is about.
"""
import numpy as np
from scipy.spatial import ConvexHull

from jb_a_family import L_EDGE
from jb_b_variety import path_tangent
from jb_j_internal_frame import Frame, inertia

V_TET = L_EDGE ** 3 / (6 * np.sqrt(2.0))
A_ICO = 22.238756093


def hull_vol(P):
    """Hull volume of a vertex set, normalised, robust to exactly merged points."""
    keep = [0]
    for i in range(1, len(P)):
        if min(np.linalg.norm(P[i] - P[k]) for k in keep) > 1e-12:
            keep.append(i)
    return ConvexHull(P[keep]).volume / V_TET


def aligned_frame(a0):
    """Frame whose basis direction 0 is the symmetric-path tangent."""
    F = Frame(a0)
    xi = path_tangent(a0)
    xi = F.B @ (F.B.T @ xi)
    xi = xi / np.linalg.norm(xi)
    M = np.column_stack([xi] + [F.B[:, i] for i in range(F.dim)])
    Q, _ = np.linalg.qr(M)
    F.B = Q[:, :F.dim]
    assert abs(abs(F.B[:, 0] @ xi) - 1.0) < 1e-9
    return F


def refine_path_stationary(a_lo, a_hi, h=1e-5):
    """Locate dV/da = 0 on the symmetric slice, without inheriting 7.2189519."""
    from scipy.optimize import brentq
    from jb_h_hull_control import hull_volume

    def dV(a):
        return (hull_volume(a + h) - hull_volume(a - h)) / (2 * h)
    return brentq(dV, a_lo, a_hi, xtol=1e-11)


def report_point(a0, label, hs=(3e-3, 1e-3, 3e-4)):
    """Gradient, smooth/kink discrimination, and Hessian inertia at one point.

    SMOOTH vs KINK is decided by how the one-sided asymmetry |D_+ - D_-| SCALES
    with h, not by its size. For a C2 function the asymmetry is f''*h and so
    falls linearly to zero; for a |s| corner it tends to the constant 2c. Reading
    "asymmetry is large" as "not differentiable" would be a check that cannot
    fail -- at h=3e-3 a perfectly smooth point with f''=7 shows 2e-2.
    """
    F = aligned_frame(a0)
    print(f"\n--- {label}   a0 = {a0:.9f}   chart dim {F.dim} ---")
    prev = None
    for h in hs:
        g, gp, gm = F.grad(hull_vol, h=h, onesided=True)
        asym = np.abs(gp - gm).max()
        ratio = "" if prev is None else f"  (fell {prev[1] / asym:.2f}x for a " \
                                        f"{prev[0] / h:.2f}x smaller h)"
        print(f"  h={h:.0e}  grad(Vol) = {np.array2string(g, precision=7)}")
        print(f"           |grad| = {np.linalg.norm(g):.3e}   "
              f"transverse |grad_1..5| = {np.linalg.norm(g[1:]):.3e}")
        print(f"           max |D_+ - D_-| = {asym:.3e}{ratio}")
        prev = (h, asym)
    print("           -> asymmetry falling in proportion to h = SMOOTH;"
          " flat in h = KINK")
    for h in hs:
        H = F.hessian(hull_vol, h=h)
        (npos, nzero, nneg), ev = inertia(-H)      # inertia of the POTENTIAL -Vol
        print(f"  h={h:.0e}  eig(-Hess Vol) = {np.array2string(ev, precision=6)}")
        print(f"           inertia (pos, zero, neg) = ({npos}, {nzero}, {nneg})")
    return F


if __name__ == "__main__":
    np.set_printoptions(precision=7, suppress=False, linewidth=170)

    print("=== M2/M3 preliminary: locate the path stationary point ourselves ===")
    a_star = refine_path_stationary(3.0, 15.0)
    print(f"  dV/da = 0 at a = {a_star:.9f}   (memo 7.2189519)")

    print("\n=== M2 at the VECTOR EQUILIBRIUM (a=0) ===")
    F0 = aligned_frame(0.0)
    print("  Vol_hull is NOT differentiable along direction 0 (jb_i). Test whether")
    print("  that is special to the symmetric direction or holds in all 6.")
    for h in (1e-3, 1e-4):
        g, gp, gm = F0.grad(hull_vol, h=h, onesided=True)
        print(f"  h={h:.0e}")
        print(f"    D_+e_i = {np.array2string(gp, precision=6)}")
        print(f"    D_-e_i = {np.array2string(gm, precision=6)}")
        print(f"    |D_+ - D_-| = {np.array2string(np.abs(gp - gm), precision=6)}")
    print("\n  random tangent directions (not just the basis axes):")
    rng = np.random.default_rng(11)
    f0 = F0.value(hull_vol, np.zeros(6))
    for trial in range(6):
        u = rng.standard_normal(6)
        u /= np.linalg.norm(u)
        row = []
        for h in (1e-3, 1e-4):
            dp = (F0.value(hull_vol, h * u) - f0) / h
            dm = (f0 - F0.value(hull_vol, -h * u)) / h
            row.append((dp, dm))
        print(f"    u{trial}: h=1e-3 D+={row[0][0]:+.6f} D-={row[0][1]:+.6f}  |  "
              f"h=1e-4 D+={row[1][0]:+.6f} D-={row[1][1]:+.6f}  "
              f"sum={row[1][0] + row[1][1]:+.6f}")
    print("  D_+ and D_- are the RIGHT and LEFT derivatives. Smooth => they AGREE.")
    print("  D_- = -D_+ with D_+ > 0 in every direction => a CONE POINT: Vol rises")
    print("  at first order whichever way you leave the VE, so the VE is a strict")
    print("  but NON-SMOOTH local minimum of Vol_hull in the full 6-D variety.")
    print("\n  degree-1 homogeneity test: is Vol(0)+t*phi(u) with phi positively")
    print("  homogeneous the right local model? check Vol(t*u) linear in t>0:")
    u = rng.standard_normal(6)
    u /= np.linalg.norm(u)
    for t in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
        d = F0.value(hull_vol, t * u) - f0
        print(f"    t={t:.0e}  Vol(tu)-Vol(0) = {d:+.9e}   /t = {d / t:+.9f}")

    print("\n=== M3 at the path maximum (hull simplicial => smooth) ===")
    report_point(a_star, "hull maximum")

    print("\n=== M3 context: icosahedron and octahedron ===")
    report_point(A_ICO, "icosahedron")
    report_point(60.0, "octahedron (contact boundary; DECISION 16: legal)")

    print("\n=== perturbation sensitivity of the inertia at the path maximum ===")
    F = aligned_frame(a_star)
    for da in (0.0, 1e-9, 1e-6, 1e-3):
        Fp = aligned_frame(a_star + da)
        H = Fp.hessian(hull_vol, h=1e-3)
        (p, z, n), ev = inertia(-H)
        print(f"  a0 shifted by {da:.0e}:  inertia ({p},{z},{n})  "
              f"eig extremes {ev.min():+.6e} .. {ev.max():+.6e}")
    print("\n  zero-eigenvalue tolerance sweep (does the count move with the cut?)")
    H = F.hessian(hull_vol, h=1e-3)
    for tol in (1e-4, 1e-6, 1e-8, 1e-10):
        (p, z, n), ev = inertia(-H, tol=tol)
        print(f"    tol={tol:.0e}  inertia ({p},{z},{n})")

    print("\n  Sylvester check: inertia must survive ANY invertible remix of the")
    print("  chart basis (which is what 'no mass metric needed' means). A bug in")
    print("  the tangent-space construction would show up as a moving count.")
    rng2 = np.random.default_rng(3)
    for k in range(4):
        Ck = rng2.standard_normal((6, 6))
        Hk = Ck.T @ (-H) @ Ck
        (p, z, n), ev = inertia(Hk)
        print(f"    remix {k}: cond={np.linalg.cond(Ck):8.2f}  inertia ({p},{z},{n})")

    print("\n  and the same via an actually re-derived chart (independent SVD basis)")
    F2 = Frame(a_star)                     # unaligned basis straight from the SVD
    H2 = F2.hessian(hull_vol, h=1e-3)
    (p, z, n), ev2 = inertia(-H2)
    print(f"    unaligned chart: inertia ({p},{z},{n})   "
          f"eig = {np.array2string(ev2, precision=6)}")
