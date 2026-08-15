"""Step R: the INTERNAL MASS METRIC in Python, and its cross-check against Java.

Every spectrum in jb_s / jb_t is a generalised eigenproblem `H v = omega^2 M v`.
`H` came from jb_k/jb_o. `M` did not exist on the Python side. This builds it,
mirroring `src/test/java/.../jitterbug/InternalMassMetric.java` WITHOUT importing
from src -- two independent implementations of the same object is how this
project has caught errors before (memo: "the 6-dimensional internal tangent
space is now constructed TWICE INDEPENDENTLY ... use both").

    M(x) = P^T diag(m) P  restricted to the internal tangent space,

with `P` the 72x48 position Jacobian (corner velocities from body motions) and
the six global rigid motions already projected out by `jb_j_internal_frame`.

TWO MASS MODELS, PEERS, NEITHER A DEFAULT (memo C2 -- a frequency reported
without naming its mass model is not a reported frequency):

  point masses   24 equal masses at the corners, spin coefficient k = 1/3,
                 M_eff(a) = (2/3) sin^2 a + 1/3      (3:1 swing)
  uniform laminae 8 uniform triangular plates,       k = 1/12,
                 M_eff(a) = (2/3) sin^2 a + 1/12     (9:1 swing)

A lamina is realised, exactly and not approximately, as three corner masses
`k*faceMass` plus a virtual centroid mass `faceMass*(1-3k)`: the centroid
velocity of a rigid triangle is exactly the mean of its corner velocities, so
the four-point body reproduces the plate's total mass AND polar second moment
with no residual. Total mass is 1/2 under BOTH models, by construction.

THE VALIDATION ANCHOR. Restricted to the symmetric-path direction the metric
built here from general 6-DOF machinery must reproduce the two closed forms
above, which were pinned by a completely different route (central differences of
the geometry, `FreeDynamicsTest`) and reproduced in Java. This is the check that
makes every frequency downstream believable, and it can fail: a wrong row layout
in P, a wrong pivot, a leaked rigid mode, or wrong lamina weights all break it.
"""
import numpy as np

from jb_a_family import corners
from jb_b_variety import cross_row, jacobian, path_tangent
from jb_j_internal_frame import Frame

# Total mass is 1/2 under both models. The value is a convention; what is NOT a
# convention is that the two models must AGREE on it, which is what makes their
# frequencies comparable at all.
CORNER_MASS = 1.0 / 48.0          # 24 corners  -> total 1/2
FACE_MASS = 1.0 / 16.0            # 8 faces     -> total 1/2
K_LAMINA = 1.0 / 12.0             # uniform triangle, polar 2nd moment L^2/12
K_POINT = 1.0 / 3.0               # 3 equal corner masses, polar 2nd moment L^2/3


def position_jacobian(X):
    """72x48 corner-velocity map. Row = 9*face + 3*corner + xyz.

    Column layout matches jb_j's chart coordinate y = (w_0..w_7, t_0..t_7):
    an infinitesimal body motion (dw_i, dt_i) moves corner j of face i by
    dw_i x r_ij + dt_i with r_ij measured from the face's own centroid --
    exactly the pivot `Frame.config` rotates about.
    """
    cen = X.mean(axis=1)
    P = np.zeros((72, 48))
    for i in range(8):
        for j in range(3):
            r = X[i, j] - cen[i]
            for d in range(3):
                row = 9 * i + 3 * j + d
                P[row, 3 * i:3 * i + 3] = cross_row(r, d)
                P[row, 24 + 3 * i + d] = 1.0
    return P


def weights(model):
    """(24 corner weights, 8 centroid weights) for a named mass model."""
    if model == "point":
        return np.full(24, CORNER_MASS), np.zeros(8)
    if model == "lamina":
        return (np.full(24, K_LAMINA * FACE_MASS),
                np.full(8, FACE_MASS * (1.0 - 3.0 * K_LAMINA)))
    raise ValueError(f"unknown mass model {model!r}")


MODELS = ("point", "lamina")


def ambient_quadratic(X, y, model):
    """Kinetic energy 2T for body-motion vector y in R^48, before restriction.

    Kept separate from `metric` so it can be probed on directions the chart
    EXCLUDES -- global translations and rotations -- which is how the weight
    bookkeeping gets a check that does not go through the chart at all.
    """
    wc, wcen = weights(model)
    v = position_jacobian(X) @ y                     # 72 corner velocities
    acc = float(np.sum(np.repeat(wc, 3) * v * v))
    vc = v.reshape(8, 3, 3).mean(axis=1)             # 8 centroid velocities
    acc += float(np.sum(np.repeat(wcen, 3) * vc.reshape(-1) ** 2))
    return acc


def metric(F, model):
    """The n x n internal mass metric in the chart `F`'s coordinates.

    `F.B` is 48 x n with orthonormal columns and dy/ds(0) = B, so a chart
    velocity sdot induces body motion B sdot and corner velocities P B sdot.
    """
    wc, wcen = weights(model)
    Vc = position_jacobian(F.X0) @ F.B                # 72 x n
    Vcen = Vc.reshape(8, 3, 3, -1).mean(axis=1).reshape(24, -1)   # 24 x n
    M = Vc.T @ (np.repeat(wc, 3)[:, None] * Vc)
    M += Vcen.T @ (np.repeat(wcen, 3)[:, None] * Vcen)
    return 0.5 * (M + M.T)


# --------------------------------------------------------------------------
# the symmetric-path generator, in closed form
# --------------------------------------------------------------------------

from jb_a_family import Z, faces                                   # noqa: E402


def symmetric_generator(a_deg):
    """d(body motion)/d(a in RADIANS) along the symmetric jitterbug path.

    Closed form, not a finite difference: face i spins about its own fixed axis
    u_i at rate sigma_i (jb_a bakes rotation angle sigma_i*(a-60deg)) and its
    centroid travels along u_i as Z*cos(a_rad).
    """
    a = np.radians(a_deg)
    z = np.zeros(48)
    for i, (_, _, u, sigma) in enumerate(faces()):
        z[3 * i:3 * i + 3] = sigma * u
        z[24 + 3 * i:24 + 3 * i + 3] = -Z * np.sin(a) * u
    return z


def internal_basis(a_deg, tol=1e-9):
    """Orthonormal basis of the internal tangent space, straight from the SVD.

    Same object `Frame` builds, but without the Newton chart -- usable at the
    rank-drop angles 90/270 where the chart's Newton solve fails and where the
    memo's C2 note says the M_eff swing must NOT be skipped.
    """
    X = corners(a_deg)
    F = Frame.__new__(Frame)
    F.X0 = X
    F.cen0 = X.mean(axis=1)
    Rg = F._global_modes()
    A = np.vstack([jacobian(X), Rg.T])
    U, S, Vt = np.linalg.svd(A)
    B = Vt[np.sum(S > tol * max(1.0, S[0])):].T
    F.B = B
    F.Rg = Rg
    F.dim = B.shape[1]
    return F


def effective_mass(a_deg, model):
    """M_eff(a): the metric contracted with the symmetric-path generator."""
    F = internal_basis(a_deg)
    z = symmetric_generator(a_deg)
    c = F.B.T @ z
    leak = np.linalg.norm(z - F.B @ c)
    return float(c @ metric(F, model) @ c), leak, F.dim


def M_eff_closed(a_deg, model):
    k = K_POINT if model == "point" else K_LAMINA
    return (2.0 / 3.0) * np.sin(np.radians(a_deg)) ** 2 + k


if __name__ == "__main__":
    np.set_printoptions(precision=9, suppress=True, linewidth=170)

    print("=" * 78)
    print("R1  WEIGHT BOOKKEEPING, probed OUTSIDE the chart")
    print("=" * 78)
    print("  A global translation moves every corner (and every centroid) at unit")
    print("  speed, so 2T must be the TOTAL MASS -- 1/2 by construction, and the")
    print("  same for both models. This is where wrong lamina weights show up:")
    print("  3*k*faceMass + faceMass*(1-3k) = faceMass identically, so a typo in")
    print("  either weight breaks it while leaving the corner weights plausible.")
    X0 = corners(0.0)
    tr = np.zeros(48)
    tr[24::3] = 1.0                      # every t_i = e_x
    for m in MODELS:
        print(f"    {m:7s} 2T(global x-translation) = {ambient_quadratic(X0, tr, m):.15f}"
              f"   (want 0.500000000000000)")
    print("  and a direction with NO corner motion must carry NO mass:")
    print(f"    2T(zero motion) = {ambient_quadratic(X0, np.zeros(48), 'lamina'):.3e}")

    print()
    print("=" * 78)
    print("R2  THE VALIDATION ANCHOR: does the 6-DOF metric reproduce M_eff(a)?")
    print("=" * 78)
    print("  point masses  ->  (2/3) sin^2 a + 1/3      (3:1)")
    print("  laminae       ->  (2/3) sin^2 a + 1/12     (9:1)")
    print("  Independently pinned in Java by a different route. Reproducing them")
    print("  from the chart machinery is a cross-check, not a restatement.")
    angles = [0, 5, 17, 22.238756093, 30, 45, 60, 75, 85, 90, 95, 120, 150,
              179, 200, 230, 260, 270, 300, 330, 355]
    worst = {m: 0.0 for m in MODELS}
    print(f"  {'a':>12s} {'dim':>4s} {'leak':>10s} "
          f"{'M_eff point':>15s} {'err':>10s} {'M_eff lamina':>15s} {'err':>10s}")
    for a in angles:
        row = [f"  {a:12.6f}"]
        dim = leak = None
        for m in MODELS:
            v, leak, dim = effective_mass(a, m)
            e = abs(v - M_eff_closed(a, m))
            worst[m] = max(worst[m], e)
            row.append(f"{v:15.12f} {e:10.2e}")
        print(f"  {a:12.6f} {dim:4d} {leak:10.2e} " + " ".join(row[1:]))
    for m in MODELS:
        print(f"  worst |measured - closed form|, {m:7s}: {worst[m]:.3e}")

    print()
    print("  SWING RATIO (memo C2: must include a=90, whose chart rank drops --")
    print("  skipping it undershoots the lamina ratio to 8.9976 by missing the peak)")
    for m in MODELS:
        vals = [effective_mass(a, m)[0] for a in np.arange(0.0, 180.5, 0.5)]
        vals.append(effective_mass(90.0, m)[0])
        print(f"    {m:7s} min {min(vals):.12f}  max {max(vals):.12f}  "
              f"ratio {max(vals) / min(vals):.9f}")

    print()
    print("=" * 78)
    print("R3  IS M POSITIVE DEFINITE? (a zero mode of M would be a zero mass)")
    print("=" * 78)
    print("  P is injective on body motions (a rigid triangle with three corners")
    print("  at rest is at rest), so M should be PD wherever the weights are")
    print("  positive. Measured, not assumed -- if it were false, every omega^2")
    print("  below would be meaningless.")
    for a in (0.0, 22.238756093, 45.0, 60.0, 90.0, 137.0):
        F = internal_basis(a)
        line = f"  a={a:11.6f} dim {F.dim}"
        for m in MODELS:
            ev = np.linalg.eigvalsh(metric(F, m))
            line += f"   {m}: eig {ev.min():.6e}..{ev.max():.6e}"
        print(line)

    print()
    print("=" * 78)
    print("R4  CLOSED-FORM GENERATOR vs FINITE-DIFFERENCE PATH TANGENT")
    print("=" * 78)
    print("  jb_b's path_tangent differentiates the geometry numerically in")
    print("  DEGREES; symmetric_generator is exact in RADIANS. They must agree")
    print("  after the 180/pi conversion. Two different derivations of the same")
    print("  vector; a sign-convention error in either shows up here.")
    for a in (0.0, 17.0, 22.238756093, 45.0, 77.0, 133.0):
        z = symmetric_generator(a)
        xi = path_tangent(a) * (180.0 / np.pi)
        F = internal_basis(a)
        print(f"  a={a:11.6f}  ||z - xi|| = {np.linalg.norm(z - xi):.3e}"
              f"   ||z|| = {np.linalg.norm(z):.9f}"
              f"   leak out of tangent space = "
              f"{np.linalg.norm(z - F.B @ (F.B.T @ z)):.3e}")

    print()
    print("=" * 78)
    print("R5  CHART COVARIANCE: M must transform as C^T M C under a basis remix")
    print("=" * 78)
    print("  The generalised eigenvalues in jb_s are chart-independent ONLY if")
    print("  H and M are both tensors. H is one at a critical point (Sylvester")
    print("  argument, jb_k); this checks M.")
    F = Frame(30.0)
    M0 = metric(F, "lamina")
    rng = np.random.default_rng(5)
    for k in range(4):
        C = rng.standard_normal((F.dim, F.dim))
        F2 = Frame(30.0)
        F2.B = F.B @ C
        M2 = metric(F2, "lamina")
        err = np.abs(M2 - C.T @ M0 @ C).max() / np.abs(M0).max()
        print(f"    remix {k}: cond {np.linalg.cond(C):8.2f}   "
              f"max rel dev from C^T M C = {err:.3e}")

    print()
    print("=" * 78)
    print("R6  GAUGE: no global rigid motion may leak into the internal metric")
    print("=" * 78)
    print("  A leaked rigid mode would be a spurious LOW-frequency (large-mass,")
    print("  zero-stiffness) direction -- the exact bug the zero-mode count in")
    print("  qvf.2 has to rule out.")
    for a in (0.0, 22.238756093, 90.0):
        F = internal_basis(a)
        print(f"  a={a:11.6f}  max |Rg^T B| = {np.abs(F.Rg.T @ F.B).max():.3e}"
              f"   max |J B| = "
              f"{np.abs(jacobian(F.X0) @ F.B).max():.3e}")
