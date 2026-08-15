"""Step O (K0, K1, K2): a FAMILY of smooth all-pairs kernels on the 12 shared vertices.

Two V candidates died by measurement (`-k*Vol_hull`, raw Thomson-as-icosahedron-
producer) and both of the non-smooth ones -- hull volume and strut clearance --
are WITNESS-SELECTION functionals: hull volume depends on which triples support a
face, segment clearance on which pair of points realises the minimum. Raw Thomson
selects nothing and came back smooth. So the surviving family is *smooth all-pairs
kernels*, and this file sweeps it:

    V_f(P) = sum_{i<j} f(|r_i - r_j|^2)

over the 12 shared vertices, for f in {1/r^p, gaussian, quadratic spread}.

STRUCTURAL FACT worth knowing before reading any number here, and CHECKED below
because it could be false: 24 of the 66 vertex pairs ARE STRUTS, whose lengths are
frozen by the rigid-triangle parameterisation. Those 24 terms are a CONSTANT for
every kernel at every configuration. Only the remaining 42 pairs carry any
a-dependence at all. A kernel sweep that did not know this would keep
re-discovering the same 42-pair functional under different names.

Each kernel is reported in TWO variants (K2), and where they disagree the
disagreement is the finding:
  raw        : f applied to the actual vertex separations (shape AND scale)
  normalised : vertices re-referenced to their CENTROID and projected to the unit
               sphere first (shape only). Centroid- not origin-referenced: jb_l /
               jb_m measured that the origin form is not translation invariant
               (49.314 -> 57.217 under a rigid shift) and so is not a function on
               the linkage variety at all.

Smoothness is decided by the SAME discriminator jb_k uses -- how the one-sided
asymmetry |D_+ - D_-| SCALES with h, not how big it is -- and is run against two
CONTROLS whose answers are already known and opposite: Vol_hull at a=0 must come
back CORNER, raw Thomson at a=0 must come back SMOOTH. If the discriminator gets
either control wrong, every smoothness verdict in this file is void.

Inertia needs no mass model: by Sylvester's law the (n_pos, n_zero, n_neg) count
is invariant under congruence. Only frequencies would need M.
"""
import numpy as np
from scipy.optimize import minimize_scalar

from jb_a_family import corners, L_EDGE
from jb_b_variety import PAIRS
from jb_j_internal_frame import Frame, inertia
from jb_k_hull_hessian import aligned_frame, hull_vol

A_ICO = 22.238756093
A_HULLMAX = 7.218951982
THOMSON_12 = 49.165253058           # published N=12 minimum (icosahedron)


# --------------------------------------------------------------------------
# vertices
# --------------------------------------------------------------------------

def vert(a):
    X = corners(a)
    return np.array([X[i, j] for (i, j), _ in PAIRS])


def sq_dists(P):
    """The 66 squared pairwise separations of a 12-point set."""
    D = P[:, None, :] - P[None, :, :]
    s = np.einsum('ijk,ijk->ij', D, D)
    iu = np.triu_indices(len(P), 1)
    return s[iu]


def to_shape(P):
    """Centroid-reference, then project to the unit sphere. Translation invariant."""
    Q = P - P.mean(axis=0)
    r = np.linalg.norm(Q, axis=1)
    if r.min() < 1e-12:
        return None
    return Q / r[:, None]


# --------------------------------------------------------------------------
# the kernel family
# --------------------------------------------------------------------------

def k_invpow(p):
    """f(s) = s^(-p/2), i.e. 1/r^p. p=1 is Thomson."""
    def f(s):
        if s.min() < 1e-24:
            return np.inf
        return float(np.sum(s ** (-0.5 * p)))
    return f


def k_gauss(sigma):
    """f(s) = exp(-s / (2 sigma^2)). Finite even where vertices merge."""
    c = 1.0 / (2.0 * sigma * sigma)
    return lambda s: float(np.sum(np.exp(-c * s)))


def k_spread(s):
    """Quadratic edge regularity: the spread of the squared pairwise distances.

    sum_(i<j) (s_ij - mean s)^2. Zero would mean all 66 separations equal (not
    attainable for 12 points), so this is a pure SHAPE-REGULARITY functional and
    the only member of the family that is not monotone in separation. It is the
    one member that could plausibly prefer the icosahedron in the RAW variant.
    """
    return float(np.sum((s - s.mean()) ** 2))


KERNELS = (
    ("1/r^1  (Thomson)", k_invpow(1.0)),
    ("1/r^2", k_invpow(2.0)),
    ("1/r^3", k_invpow(3.0)),
    ("1/r^6", k_invpow(6.0)),
    ("1/r^12", k_invpow(12.0)),
    ("gauss s=0.5", k_gauss(0.5)),
    ("gauss s=1.0", k_gauss(1.0)),
    ("gauss s=1.5", k_gauss(1.5)),
    ("gauss s=2.5", k_gauss(2.5)),
    ("spread(quadratic)", k_spread),
)


def make_V(kern, normalised):
    """Turn a kernel on squared separations into a potential on a 12-point set."""
    if not normalised:
        return lambda P: kern(sq_dists(P))

    def V(P):
        U = to_shape(P)
        if U is None:
            return np.inf
        return kern(sq_dists(U))
    return V


# --------------------------------------------------------------------------
# a shared Newton cache: the chart solve does not depend on the kernel
# --------------------------------------------------------------------------

class Probe:
    """Finite differences in a chart, sharing Newton solves across kernels.

    Frame.solve() is the expensive part and depends only on s, so solving once
    per stencil point and evaluating every kernel on the result is exact, not an
    approximation -- the stencil is identical to the one jb_k uses per kernel.
    """

    def __init__(self, F, tol=1e-9):
        self.F = F
        self.n = F.dim
        self.tol = tol
        self._c = {}

    def pts(self, s):
        key = tuple(np.round(s, 15))
        if key not in self._c:
            y, r = self.F.solve(s)
            assert r < self.tol, f"Newton failed at s={s}: {r:.2e}"
            self._c[key] = (self.F.vertices(y), self.F.config(y))
        return self._c[key]

    def val(self, f, s, on_config=False):
        P, X = self.pts(s)
        return f(X) if on_config else f(P)

    def grad(self, f, h, on_config=False):
        n = self.n
        f0 = self.val(f, np.zeros(n), on_config)
        g, gp, gm = np.zeros(n), np.zeros(n), np.zeros(n)
        for i in range(n):
            e = np.zeros(n)
            e[i] = h
            fp = self.val(f, e, on_config)
            fm = self.val(f, -e, on_config)
            gp[i] = (fp - f0) / h
            gm[i] = (f0 - fm) / h
            g[i] = (fp - fm) / (2 * h)
        return g, gp, gm

    def hess(self, f, h, on_config=False):
        n = self.n
        f0 = self.val(f, np.zeros(n), on_config)
        H = np.zeros((n, n))
        for i in range(n):
            e = np.zeros(n)
            e[i] = h
            H[i, i] = (self.val(f, e, on_config) - 2 * f0
                       + self.val(f, -e, on_config)) / h ** 2
        for i in range(n):
            for j in range(i + 1, n):
                ei, ej = np.zeros(n), np.zeros(n)
                ei[i] = h
                ej[j] = h
                H[i, j] = H[j, i] = (
                    self.val(f, ei + ej, on_config)
                    - self.val(f, ei - ej, on_config)
                    - self.val(f, -ei + ej, on_config)
                    + self.val(f, -ei - ej, on_config)) / (4 * h ** 2)
        return 0.5 * (H + H.T)


def smoothness(probe, f, hs=(1e-2, 1e-3, 1e-4), on_config=False):
    """CORNER vs SMOOTH from how |D_+ - D_-| scales with h.

    A C2 function has asymmetry f''*h -> falls linearly. A |s| corner has
    asymmetry 2c -> flat. Reading 'asymmetry is large' as 'not differentiable'
    would be a check that cannot fail, so the verdict is the RATIO, and a
    borderline ratio is reported as such rather than forced to a side.
    """
    a = []
    for h in hs:
        _, gp, gm = probe.grad(f, h, on_config)
        a.append(float(np.abs(gp - gm).max()))
    r = a[-1] / a[-2] if a[-2] > 0 else 0.0
    if a[0] < 1e-11:
        v = "SMOOTH(flat)"
    elif r < 0.25:
        v = "SMOOTH"
    elif r > 0.6:
        v = "CORNER"
    else:
        v = "AMBIGUOUS"
    return v, a, r


# --------------------------------------------------------------------------
# path scan
# --------------------------------------------------------------------------

def minimise_on_path(f_of_a, grid=None):
    """Global minimum of a scalar function of the path angle, over the WHOLE
    circle (DECISION 16: no region cut).

    Coarse grid, then a BOUNDED refinement on the bracketing grid cell. Bounded
    rather than Brent-with-a-3-point-bracket: when the coarse minimum sits next
    to the true one the 3-point bracket is invalid, scipy declines, and the
    routine silently returns the GRID value -- which then looks like an unstable
    argmin under a grid change. That false alarm was observed and is what this
    form removes.

    Merge angles 60/120/240/300 are poles for the inverse powers; the grid steps
    over them rather than excluding any interval, so a minimum inside (60,120)
    would still be found.
    """
    if grid is None:
        grid = np.arange(-179.9, 180.0, 0.1)
    vals = np.array([f_of_a(a) for a in grid])
    vals = np.where(np.isfinite(vals), vals, np.inf)
    k = int(np.argmin(vals))
    a0 = grid[k]
    step = float(grid[1] - grid[0])
    try:
        r = minimize_scalar(f_of_a, bounds=(a0 - step, a0 + step),
                            method="bounded", options={"xatol": 1e-11})
        if np.isfinite(r.fun) and r.fun <= vals[k] + 1e-12:
            return float(r.x), float(r.fun), a0, vals
    except Exception:
        pass
    return float(a0), float(vals[k]), a0, vals


def path_min(V, grid=None):
    """minimise_on_path for a potential defined on the 12-vertex set."""
    return minimise_on_path(lambda a: V(vert(a)), grid)


def fund(a):
    """Fold an angle to the fundamental domain [0,90] of the MEASURED symmetries
    a -> -a and a -> 180-a. Two argmins are the same state iff they fold to the
    same representative; comparing raw angles calls a mirror image a
    disagreement."""
    b = abs(a) % 360.0
    if b > 180.0:
        b = 360.0 - b
    return 180.0 - b if b > 90.0 else b


def name_angle(a, tol=1e-4):
    """Report the LOCATED angle; naming is a comment, never a snap.

    The mirror images under the a -> 180-a ISOMETRY measured below are named as
    such: 157.761... is not a new configuration, it is the icosahedron again.
    """
    b = abs(a)
    for nm, ref in (("VE", 0.0), ("hull-max", A_HULLMAX), ("icosahedron", A_ICO),
                    ("octahedron", 60.0), ("a=90", 90.0)):
        if abs(b - ref) < tol:
            return f"= {nm}"
        if abs(b - (180.0 - ref)) < tol:
            return f"= {nm}'" if ref else "= VE'"
    return "(unnamed)"


def gram_spectrum(P):
    """Eigenvalues of the centred Gram matrix -- an isometry invariant."""
    Q = P - P.mean(axis=0)
    return np.sort(np.linalg.eigvalsh(Q @ Q.T))


def sorted_distance_matrix(P):
    """Doubly sorted pairwise-distance matrix. Finer than the pooled spectrum."""
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    return np.sort(np.sort(D, axis=1), axis=0)


# --------------------------------------------------------------------------

def report_kernel(label, V, a_star, hs=(1e-2, 1e-3, 1e-4), hess_hs=(1e-3, 3e-4)):
    F = aligned_frame(a_star)
    pr = Probe(F)
    verd, asym, ratio = smoothness(pr, V, hs)
    g, _, _ = pr.grad(V, 1e-3)
    rows = []
    for h in hess_hs:
        H = pr.hess(V, h)
        (p, z, n), ev = inertia(H)
        rows.append((h, (p, z, n), ev))
    return dict(label=label, dim=F.dim, verdict=verd, asym=asym, ratio=ratio,
                grad=g, gnorm=float(np.linalg.norm(g)),
                gtrans=float(np.linalg.norm(g[1:])), hess=rows)


if __name__ == "__main__":
    np.set_printoptions(precision=7, suppress=False, linewidth=180)

    print("=" * 78)
    print("STRUCTURAL PRE-CHECK: are 24 of the 66 vertex pairs frozen struts?")
    print("=" * 78)
    s0 = np.sort(sq_dists(vert(0.0)))
    frozen = []
    for a in (0.0, 13.7, A_ICO, 45.0, 90.0, 137.0, 180.0):
        s = np.sort(sq_dists(vert(a)))
        frozen.append(np.sum(np.abs(s - L_EDGE ** 2) < 1e-12))
    print(f"  pairs at exactly L_edge^2 = {L_EDGE**2:.12f}, over a = "
          "0/13.7/a_ico/45/90/137/180:")
    print(f"    counts = {frozen}")
    print("  A kernel therefore sees a CONSTANT from those pairs and varies only")
    print("  over the other 42. (Check can fail: a non-rigid parameterisation, or")
    print("  the wrong shared-vertex pairing, would not give a constant 24.)")

    print("\n" + "=" * 78)
    print("STRUCTURAL PRE-CHECK 2: is  a <-> 180-a  an ISOMETRY of the 12 vertices?")
    print("=" * 78)
    print("  This decides how every minimum below must be READ. If the two are")
    print("  isometric then no pairwise-distance functional can tell them apart, the")
    print("  fundamental domain of the symmetric path is [0,90], and a kernel that")
    print("  'minimises at 157.761' has minimised at the icosahedron.")
    print("  Three invariants, each of which could disagree:")
    print(f"  {'a':>10s} {'180-a':>10s} {'sorted dists':>14s} {'Gram spectrum':>15s}"
          f" {'per-vertex dists':>17s} {'mean radius ratio':>18s}")
    for a in (0.0, 7.5, A_ICO, 30.0, 45.0, 59.0, 61.0, 88.0):
        Pa, Pb = vert(a), vert(180.0 - a)
        d1 = np.abs(np.sort(sq_dists(Pa)) - np.sort(sq_dists(Pb))).max()
        d2 = np.abs(gram_spectrum(Pa) - gram_spectrum(Pb)).max()
        d3 = np.abs(sorted_distance_matrix(Pa) - sorted_distance_matrix(Pb)).max()
        ra = np.linalg.norm(Pa - Pa.mean(0), axis=1).mean()
        rb = np.linalg.norm(Pb - Pb.mean(0), axis=1).mean()
        print(f"  {a:10.4f} {180.0 - a:10.4f} {d1:14.3e} {d2:15.3e} {d3:17.3e}"
              f" {rb / ra:18.12f}")
    print("  (could have failed: the two arms of the path are built from DIFFERENT")
    print("   rotation angles sigma*(a-60) and different translations Z*cos a, so")
    print("   there is no parameterisation-level reason for them to coincide)")
    print("  CONSEQUENCE: the symmetric path has TWO copies of the VE (a=0, a=180)")
    print("  and FOUR of the icosahedron (+/-22.2388, +/-157.7612). Every 'ground")
    print("  state' below is an ORBIT, not a point, and ties between mirror images")
    print("  are exact, not numerical.")

    print("\n" + "=" * 78)
    print("K0 CONTROL: reproduce raw Thomson's prior numbers with this session's code")
    print("=" * 78)
    Vraw1 = make_V(k_invpow(1.0), False)
    Vnrm1 = make_V(k_invpow(1.0), True)
    print(f"  raw Thomson at a=0    : {Vraw1(vert(0.0)):.9f}   (prior 34.889842)")
    print(f"  raw Thomson at a_ico  : {Vraw1(vert(A_ICO)):.9f}   (prior 36.554172)")
    a_r, v_r, a_c, _ = path_min(Vraw1)
    print(f"  raw Thomson argmin    : a = {a_r:.9f}  V = {v_r:.9f}  "
          f"{name_angle(a_r, 1e-3)}")
    print(f"  normalised at a_ico   : {Vnrm1(vert(A_ICO)):.9f}   "
          f"published N=12 Thomson {THOMSON_12:.9f}   "
          f"|diff| {abs(Vnrm1(vert(A_ICO)) - THOMSON_12):.2e}")
    print("     (this check CAN fail: a wrong vertex set misses by ~0.4;")
    print(f"      the cuboctahedron gives {Vnrm1(vert(0.0)):.9f}, which must be larger)")
    F0 = aligned_frame(0.0)
    pr0 = Probe(F0)
    g0, _, _ = pr0.grad(Vraw1, 1e-3)
    H0 = pr0.hess(Vraw1, 1e-3)
    (p0, z0, n0), ev0 = inertia(H0)
    print(f"  raw Thomson at VE: |grad| = {np.linalg.norm(g0):.3e}  "
          f"transverse {np.linalg.norm(g0[1:]):.3e}")
    print(f"     eig = {np.array2string(ev0, precision=6)}")
    print(f"     inertia = ({p0},{z0},{n0})   (prior: (6,0,0), spectrum "
          "[1.348861 x2, 1.537655 x3, 2.515439])")
    v, asym, ratio = smoothness(pr0, Vraw1)
    print(f"  raw Thomson smoothness at VE: {v}  asym {asym}  ratio {ratio:.3f}")
    vh, ah, rh = smoothness(pr0, hull_vol)
    print(f"  CONTROL Vol_hull at VE      : {vh}  asym "
          f"{[f'{x:.6f}' for x in ah]}  ratio {rh:.3f}")
    print("  The two controls must come back OPPOSITE. If Vol_hull reads SMOOTH or")
    print("  Thomson reads CORNER, every verdict below is void.")

    print("\n" + "=" * 78)
    print("K1/K2: the kernel sweep -- where is the minimum, and is V smooth there?")
    print("=" * 78)
    print("  scan: a in [-179.9, 180) step 0.1, then Brent. No region excluded.")
    results = []
    for norm in (False, True):
        tag = "normalised(centroid)" if norm else "raw"
        print(f"\n  ---------- variant: {tag} ----------")
        print(f"  {'kernel':20s} {'argmin a':>14s} {'V(min)':>16s}  {'named':14s}"
              f" {'V(VE)':>16s} {'V(ico)':>16s}")
        for kn, kf in KERNELS:
            V = make_V(kf, norm)
            a_s, v_s, a_c, vals = path_min(V)
            vve, vic = V(vert(0.0)), V(vert(A_ICO))
            print(f"  {kn:20s} {a_s:14.8f} {v_s:16.8f}  {name_angle(a_s, 1e-3):14s}"
                  f" {vve:16.8f} {vic:16.8f}")
            results.append((tag, kn, V, a_s, v_s, vve, vic))

    print("\n" + "=" * 78)
    print("K1: smoothness, 6-D gradient, and Hessian INERTIA at each located minimum")
    print("=" * 78)
    print("  (chart basis direction 0 IS the symmetric-path tangent, so 'transverse'")
    print("   is the genuinely 5-dimensional off-slice gradient)")
    for tag, kn, V, a_s, v_s, vve, vic in results:
        try:
            R = report_kernel(kn, V, a_s)
        except AssertionError as e:
            print(f"  {tag:20s} {kn:20s} CHART FAILURE: {e}")
            continue
        print(f"\n  [{tag}] {kn}   at a = {a_s:.9f}   chart dim {R['dim']}")
        print(f"     |grad| = {R['gnorm']:.3e}   transverse |grad_1..5| = "
              f"{R['gtrans']:.3e}")
        print(f"     smoothness: {R['verdict']}   |D+-D-| at h=1e-2,1e-3,1e-4 = "
              f"{R['asym'][0]:.3e} {R['asym'][1]:.3e} {R['asym'][2]:.3e}"
              f"   (ratio {R['ratio']:.3f})")
        for h, (p, z, n), ev in R['hess']:
            print(f"     h={h:.0e}  eig(Hess V) = {np.array2string(ev, precision=6)}")
            print(f"              inertia (pos,zero,neg) = ({p},{z},{n})"
                  f"   {'MINIMUM in all 6' if n == 0 and z == 0 else 'NOT a minimum'}")

    print("\n" + "=" * 78)
    print("THE SPREAD KERNEL IS NOT A SHAPE FUNCTIONAL. Both variants collapse.")
    print("=" * 78)
    Vsr, Vsn = make_V(k_spread, False), make_V(k_spread, True)
    print("  (a) raw spread == 768/11 * R^4 EXACTLY on the symmetric path, where the")
    print("      12 vertices are cospherical. It therefore carries ZERO shape")
    print("      information there and its 'ground state' a=+/-90 is just the")
    print("      SMALLEST configuration (R = 0.8165, the total-collision state).")
    print(f"  {'a':>10s} {'R':>16s} {'spread_raw':>16s} {'768/11 * R^4':>16s} {'diff':>11s}")
    for a in (0.0, 5.0, A_ICO, 45.0, 60.0, 75.0, 90.0, 137.0, 180.0):
        P = vert(a)
        R = np.linalg.norm(P - P.mean(0), axis=1).mean()
        pred = 768.0 / 11.0 * R ** 4
        print(f"  {a:10.4f} {R:16.12f} {Vsr(P):16.9f} {pred:16.9f} {Vsr(P) - pred:11.2e}")
    print("      (an identity, not a fit: 768/11 was READ OFF the normalised value")
    print("       and R^4 is forced by the quartic homogeneity of the functional;")
    print("       the check can fail off-cosphericity, and does -- see (c))")
    print("\n  (b) normalised spread is CONSTANT along the whole path:")
    vals = [Vsn(vert(a)) for a in (0.0, 5.0, A_ICO, 45.0, 60.0, 75.0, 90.0, 120.0,
                                   180.0, -33.3)]
    print(f"      values at a = 0/5/a_ico/45/60/75/90/120/180/-33.3 span "
          f"{min(vals):.12f} .. {max(vals):.12f}")
    print(f"      768/11 = {768/11:.12f}")
    print("  (c) so it has an exact ZERO MODE along the symmetric-path tangent,")
    print("      which is the SAME failure mode as soft joints: no restoring force")
    print("      in the direction the jitterbug actually moves. Off the path it is")
    print("      not constant, which is why the other 5 eigenvalues are positive.")
    for a0 in (30.0, A_ICO, 50.0):
        F = aligned_frame(a0)
        pr = Probe(F)
        H = pr.hess(Vsn, 1e-3)
        ev, evec = np.linalg.eigh(H)
        g, _, _ = pr.grad(Vsn, 1e-3)
        print(f"      a0={a0:9.5f}  |grad|={np.linalg.norm(g):.2e}  "
              f"eig = {np.array2string(ev, precision=8)}")
        print(f"                    |<smallest eigenvector, path tangent>| = "
              f"{abs(evec[0, 0]):.12f}  (1 => the flat direction IS the jitterbug)")

    print("\n" + "=" * 78)
    print("RAW SPREAD'S MINIMUM SITS ON A BRANCH POINT. The 6-D chart degrades there.")
    print("=" * 78)
    print("  R3 records local dimension 7 at a=90/270. Brent stopped at 90.000000084,")
    print("  where the chart still reports 6 -- so the inertia (6,0,0) printed above")
    print("  was taken NEXT TO the minimum, not AT it. Redone honestly:")
    for a0 in (89.99, 89.999, 90.0):
        F = Frame(a0)
        try:
            pr = Probe(F, tol=1e-7)
            g, _, _ = pr.grad(Vsr, 1e-3)
            H = pr.hess(Vsr, 1e-3)
            (p, z, n), ev = inertia(H)
            print(f"    a0={a0:10.5f}  chart dim {F.dim}  sv[41]={F.sv[41]:.2e}  "
                  f"|grad|={np.linalg.norm(g):.2e}  inertia=({p},{z},{n})")
            print(f"        eig = {np.array2string(ev, precision=6)}")
        except AssertionError as e:
            print(f"    a0={a0:10.5f}  chart dim {F.dim}  NEWTON FAILS: {e}")
    print("    A ground state at a point where the chart's own Newton solve fails is")
    print("    not a measured spectrum. Reported as a chart limitation, not a result.")

    print("\n" + "=" * 78)
    print("PERTURBATION: do the located minima and the inertias MOVE?")
    print("=" * 78)
    print("  (a) re-locate the argmin from a coarser and a finer grid.")
    print("      Where the answer jumps between an angle and its 180-a mirror that is")
    print("      an EXACT TIE (checked above to 1e-14), not an unstable argmin --")
    print("      the grid simply meets a different member of the same orbit first.")
    for kn, kf in (KERNELS[0], KERNELS[4], KERNELS[6], KERNELS[9]):
        for norm in (False, True):
            V = make_V(kf, norm)
            outs = []
            for step in (0.5, 0.1, 0.037):
                a_s, v_s, _, _ = path_min(V, np.arange(-179.9, 180.0, step))
                outs.append(a_s)
            print(f"    {kn:20s} {'norm' if norm else 'raw ':4s}  argmin at grid "
                  f"0.5/0.1/0.037 = {outs[0]:.7f} {outs[1]:.7f} {outs[2]:.7f}")

    print("\n  (b) shift the chart anchor off the located minimum and re-take inertia")
    for kn, kf in (KERNELS[0], KERNELS[9]):
        for norm in (False, True):
            V = make_V(kf, norm)
            a_s, _, _, _ = path_min(V)
            line = []
            for da in (0.0, 1e-9, 1e-6, 1e-3):
                try:
                    F = aligned_frame(a_s + da)
                    H = Probe(F).hess(V, 1e-3)
                    (p, z, n), _ = inertia(H)
                    line.append(f"{da:.0e}:({p},{z},{n})")
                except AssertionError:
                    line.append(f"{da:.0e}:chart-fail")
            print(f"    {kn:20s} {'norm' if norm else 'raw ':4s}  {'  '.join(line)}")

    print("\n  (c) Sylvester: inertia must survive any invertible remix of the basis")
    rng = np.random.default_rng(5)
    for kn, kf in (KERNELS[0], KERNELS[4], KERNELS[9]):
        V = make_V(kf, False)
        a_s, _, _, _ = path_min(V)
        try:
            H = Probe(aligned_frame(a_s)).hess(V, 1e-3)
        except AssertionError as e:
            print(f"    {kn:20s} raw  NO HESSIAN AT ITS OWN MINIMUM: {e}")
            print("      (this is the a=90 branch point again -- the refined argmin now")
            print("       lands exactly on it, where the 6-D chart's Newton solve fails)")
            continue
        out = []
        for k in range(4):
            C = rng.standard_normal((6, 6))
            (p, z, n), _ = inertia(C.T @ H @ C)
            out.append(f"({p},{z},{n})")
        print(f"    {kn:20s} raw  remixes: {' '.join(out)}")

    print("\n  (d) one-ULP-scale jiggle of the anchor vertex coordinates")
    rngj = np.random.default_rng(20260812)
    for kn, kf in (KERNELS[0], KERNELS[4], KERNELS[9]):
        V = make_V(kf, False)
        line = []
        for sc in (0.0, 2.2e-16, 1e-13, 1e-10):
            X = corners(0.0)
            X = X + sc * rngj.standard_normal(X.shape) * np.maximum(np.abs(X), 1.0)
            try:
                F = Frame(0.0, X0=X)
                H = Probe(F).hess(V, 1e-3)
                (p, z, n), _ = inertia(H)
                line.append(f"{sc:.1e}:({p},{z},{n})")
            except AssertionError:
                line.append(f"{sc:.1e}:chart-fail")
        print(f"    {kn:20s} at VE: {'  '.join(line)}")

    print("\n" + "=" * 78)
    print("K1 headline question: does ANY kernel put the ground state somewhere")
    print("other than the vector equilibrium, while staying SMOOTH there?")
    print("=" * 78)
    print("  'the VE' means the ORBIT {0, 180} under the measured a -> 180-a isometry;")
    print("  'the icosahedron' means {+/-22.2388, +/-157.7612}. Margin is against the")
    print("  VE value, so a positive margin is a genuinely lower ground state.")
    print(f"  {'variant':22s} {'kernel':20s} {'argmin':>13s} {'where':14s} "
          f"{'V(VE) - V(min)':>16s}")
    for tag, kn, V, a_s, v_s, vve, vic in results:
        at_ve = min(abs(abs(a_s) - 0.0), abs(abs(a_s) - 180.0)) < 1e-3
        print(f"  {tag:22s} {kn:20s} {a_s:13.7f} "
              f"{('AT the VE' if at_ve else 'OFF the VE'):14s} {vve - v_s:16.9f}")
