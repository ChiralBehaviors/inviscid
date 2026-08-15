"""Step Q (K4): kernels on the STRUTS, not on the vertices.

Everything in jb_l and jb_o is a kernel on the 12 shared VERTICES. That is an
unstated modelling choice of exactly the kind memo C2 caught for the mass model:
Fuller's jitterbug (Synergetics 460.011) is 24 STRUTS cohered at 12 flexible
joints -- there are no faces, and the vertices are where struts meet, not the
primitives. If a strut-based potential and a vertex-based potential disagree
about the ground state, the choice of primitive is load-bearing and must be
declared rather than inherited.

Strut representatives measured here, all SMOOTH BY CONSTRUCTION (unlike
segment-to-segment clearance, which is a min over the segments and is exactly
the witness-selection object jb_p shows kinks):

  midpoints  : the 24 strut midpoints, all-pairs kernel over 276 pairs. The
               construction the brief names.
  directions : struts also carry an ORIENTATION the midpoints throw away. A
               kernel on the 24 unsigned axes (via |u_i . u_j|, which is
               sign-free and smooth) tests whether the discarded information
               changes the answer. Included because "midpoints only" is itself
               a modelling choice smuggled in one level down.
  m+d combo  : midpoint separation weighted by axis alignment, the cheapest
               functional that sees both.

Same two variants as jb_o: raw and centroid-normalised (shape only). Same
smoothness discriminator, same controls, same Sylvester-invariant inertia.

A degeneracy to know before reading anything: at a = 60/120/240/300 the struts
COINCIDE IN PAIRS, so 12 midpoint pairs are at distance 0 and the inverse
powers have a pole -- the same failure the vertex kernels have at the same
angles, for a different reason.
"""
import numpy as np

from jb_a_family import corners
from jb_j_internal_frame import Frame, inertia
from jb_k_hull_hessian import aligned_frame
from jb_o_kernel_family import (KERNELS, Probe, fund, gram_spectrum, k_gauss,
                                k_invpow, k_spread, make_V, minimise_on_path,
                                name_angle, path_min, smoothness, sq_dists,
                                vert, A_ICO)


def strut_ends(X):
    """The 24 struts as (p0, p1) arrays, straight off the rigid triangles."""
    p0 = np.array([X[i, j] for i in range(8) for j in range(3)])
    p1 = np.array([X[i, (j + 1) % 3] for i in range(8) for j in range(3)])
    return p0, p1


def midpoints(X):
    p0, p1 = strut_ends(X)
    return 0.5 * (p0 + p1)


def axes(X):
    p0, p1 = strut_ends(X)
    u = p1 - p0
    return u / np.linalg.norm(u, axis=1)[:, None]


def to_shape(P):
    Q = P - P.mean(axis=0)
    r = np.linalg.norm(Q, axis=1)
    if r.min() < 1e-12:
        return None
    return Q / r[:, None]


def mid_V(kern, normalised):
    """All-pairs kernel on the 24 strut midpoints. Takes a CONFIG X, not P."""
    def V(X):
        M = midpoints(X)
        if normalised:
            M = to_shape(M)
            if M is None:
                return np.inf
        return kern(sq_dists(M))
    return V


def axis_V(kern):
    """All-pairs kernel on the 24 unsigned strut AXES.

    Uses s_ij = 2 - 2*|u_i . u_j| in [0, 2]: zero for parallel axes, 2 for
    perpendicular. Smooth EXCEPT where two axes are exactly perpendicular
    (|.| has a corner at 0) -- which does happen on this family, so the
    smoothness verdict for this member is a genuine question, not a formality.
    """
    def V(X):
        U = axes(X)
        C = np.abs(U @ U.T)
        iu = np.triu_indices(len(U), 1)
        return kern(2.0 - 2.0 * C[iu])
    return V


def combo_V(kern, lam=1.0):
    """Midpoint separation modulated by axis alignment. Sees both.

    s_ij = |m_i - m_j|^2 + lam*(1 - (u_i . u_j)^2). The squared dot product is
    smooth and sign-free, so unlike axis_V this member has no |.| anywhere.
    """
    def V(X):
        M = midpoints(X)
        U = axes(X)
        D = M[:, None, :] - M[None, :, :]
        s = np.einsum('ijk,ijk->ij', D, D)
        C = U @ U.T
        s = s + lam * (1.0 - C ** 2)
        iu = np.triu_indices(len(M), 1)
        return kern(s[iu])
    return V


def path_min_cfg(V, grid=None):
    """Minimum over the whole circle for a potential defined on the CONFIG."""
    a, v, _, _ = minimise_on_path(lambda t: V(corners(t)), grid)
    return a, v


if __name__ == "__main__":
    np.set_printoptions(precision=7, suppress=False, linewidth=180)

    print("=" * 78)
    print("PRE-CHECK: the strut set, and where it degenerates")
    print("=" * 78)
    for a in (0.0, A_ICO, 45.0, 59.0, 60.0, 90.0, 120.0, 180.0):
        X = corners(a)
        M = midpoints(X)
        d = np.sqrt(sq_dists(M))
        U = axes(X)
        C = np.abs(U @ U.T)
        iu = np.triu_indices(24, 1)
        rM = np.linalg.norm(M - M.mean(0), axis=1)
        print(f"  a={a:8.3f}  24 midpoints: min sep {d.min():.9f}  "
              f"radius {rM.min():.6f}..{rM.max():.6f}   "
              f"parallel axis pairs {int(np.sum(C[iu] > 1 - 1e-9)):3d}  "
              f"perpendicular {int(np.sum(C[iu] < 1e-9)):3d}")
    print("  Midpoints merge at 60/120 (struts coincide) -> inverse powers pole there,")
    print("  exactly as the vertex kernels do, and for an independent reason.")

    print("\n  is a <-> 180-a an isometry of the MIDPOINT set too? (jb_o measured it")
    print("  for the vertices; it need not carry over, and the answer decides how")
    print("  the strut table is read)")
    for a in (0.0, 7.5, A_ICO, 45.0, 88.0):
        Ma, Mb = midpoints(corners(a)), midpoints(corners(180.0 - a))
        d1 = np.abs(np.sort(sq_dists(Ma)) - np.sort(sq_dists(Mb))).max()
        d2 = np.abs(gram_spectrum(Ma) - gram_spectrum(Mb)).max()
        print(f"    a={a:8.3f} vs {180.0-a:8.3f}: sorted dists {d1:.3e}   "
              f"Gram spectrum {d2:.3e}")
    print("\n  WARNING, and it is a check-that-cannot-fail caught in THIS file's own")
    print("  machinery: the Gram-spectrum column above is VACUOUS for the NORMALISED")
    print("  midpoints. Once projected to the unit sphere the 24 midpoints are")
    print("  ISOTROPIC at every angle -- spectrum (0 x21, 8, 8, 8) -- so it matches")
    print("  everything against everything. It discriminates for the RAW midpoints")
    print("  only. Demonstrated rather than asserted:")
    for a in (0.0, 17.0, 43.0, 79.1066058):
        w = gram_spectrum(to_shape(midpoints(corners(a))))
        print(f"    a={a:11.6f}  normalised-midpoint Gram spectrum, top 3 = "
              f"{np.array2string(w[-3:], precision=12)}   rest max |.| = "
              f"{np.abs(w[:-3]).max():.2e}")

    print("\n" + "=" * 78)
    print("K4 TABLE: strut-MIDPOINT kernels vs the vertex kernels of jb_o")
    print("=" * 78)
    print(f"  {'kernel':20s} {'variant':10s} {'argmin a':>14s} {'V(min)':>16s}"
          f"  {'named':16s} {'vertex-kernel argmin':>22s}")
    rows = []
    for kn, kf in KERNELS:
        for norm in (False, True):
            Vs = mid_V(kf, norm)
            a_s, v_s = path_min_cfg(Vs)
            Vv = make_V(kf, norm)
            a_v, v_v, _, _ = path_min(Vv)
            print(f"  {kn:20s} {'norm' if norm else 'raw':10s} {a_s:14.8f} "
                  f"{v_s:16.8f}  {name_angle(a_s, 1e-3):16s} "
                  f"{a_v:14.8f} {name_angle(a_v, 1e-3):>8s}")
            rows.append((kn, norm, Vs, a_s, v_s, a_v))

    print("\n" + "=" * 78)
    print("K4: smoothness, 6-D gradient and inertia at the strut-kernel minima")
    print("=" * 78)
    for kn, norm, Vs, a_s, v_s, a_v in rows:
        if kn == "spread(quadratic)":
            continue                     # handled separately below
        try:
            F = aligned_frame(a_s)
            pr = Probe(F)
            verd, asym, ratio = smoothness(pr, Vs, on_config=True)
            g, _, _ = pr.grad(Vs, 1e-3, on_config=True)
            H = pr.hess(Vs, 1e-3, on_config=True)
            (p, z, n), ev = inertia(H)
        except AssertionError as e:
            print(f"  {kn:20s} {'norm' if norm else 'raw':5s} CHART FAILURE: {e}")
            continue
        print(f"  {kn:20s} {'norm' if norm else 'raw':5s} a={a_s:13.8f}  "
              f"{verd:12s} ratio {ratio:.3f}  |grad|={np.linalg.norm(g):.2e} "
              f"transv={np.linalg.norm(g[1:]):.2e}  inertia=({p},{z},{n})")
        print(f"       eig = {np.array2string(ev, precision=6)}")

    print("\n" + "=" * 78)
    print("K4: does the SPREAD kernel collapse on struts the way it did on vertices?")
    print("=" * 78)
    print("  On the vertices, raw spread was exactly 768/11 * R^4 (pure scale) and")
    print("  normalised spread was constant. MEASURED ABOVE, and it decides this:")
    print("  the 24 midpoints are ALSO cospherical at every a (radius min == max to")
    print("  the printed precision), so the same collapse is available to them. The")
    print("  question is whether it actually happens -- a different point set can be")
    print("  cospherical and still have a varying normalised spread.")
    Vsr, Vsn = mid_V(k_spread, False), mid_V(k_spread, True)
    print(f"  {'a':>10s} {'midpoint R spread':>20s} {'spread_raw':>18s} "
          f"{'spread_norm':>18s}")
    for a in (0.0, 5.0, A_ICO, 45.0, 60.0, 75.0, 90.0, 180.0):
        X = corners(a)
        M = midpoints(X)
        r = np.linalg.norm(M - M.mean(0), axis=1)
        print(f"  {a:10.4f} {r.max()-r.min():20.12f} {Vsr(X):18.9f} "
              f"{Vsn(X):18.9f}")
    print("  A non-constant normalised column means the strut version carries SHAPE")
    print("  information the vertex version did not. A constant one means the")
    print("  collapse is structural, not an artefact of choosing vertices.")
    print("\n  explicit cosphericity check on the midpoints (this CAN fail: the 24")
    print("  midpoints have no a-priori reason to lie on a sphere, and if they did")
    print("  not, normalising them would not be a pure rescale):")
    for a in (0.0, 3.7, A_ICO, 41.0, 88.0, 133.0):
        r = np.linalg.norm(midpoints(corners(a)) - midpoints(corners(a)).mean(0),
                           axis=1)
        print(f"    a={a:8.3f}  radius spread over the 24 midpoints = "
              f"{r.max()-r.min():.3e}   R = {r.mean():.12f}")

    print("\n" + "=" * 78)
    print("K4b: kernels that use the strut AXES, which midpoints discard")
    print("=" * 78)
    print(f"  {'functional':28s} {'argmin a':>14s} {'V(min)':>16s}  {'named':16s}"
          f" {'smooth at min?':>16s}")
    print("  NOT INCLUDED, and the reason is a measurement: an inverse-power kernel")
    print("  on the axes is ILL-POSED on this family. s_ij = 2-2|u.u'| is ZERO for")
    print("  parallel struts, and the pre-check counts 12-36 parallel pairs at EVERY")
    print("  angle tested, so 1/s^(p/2) is +inf everywhere. Not a near-miss to")
    print("  regularise -- there is no angle at which it is finite.")
    for a in (0.0, 30.0, 77.0, 133.0):
        U = axes(corners(a))
        C = np.abs(U @ U.T)
        iu = np.triu_indices(24, 1)
        print(f"    a={a:7.2f}  min over pairs of 2-2|u.u'| = "
              f"{(2 - 2 * C[iu]).min():.3e}   parallel pairs "
              f"{int(np.sum(C[iu] > 1 - 1e-9))}")
    EXTRA = [("axis gauss s=1.0", axis_V(k_gauss(1.0)))]
    EXTRA += [("axis gauss s=0.5", axis_V(k_gauss(0.5)))]
    EXTRA += [(f"mid+axis combo 1/r^{p:g}", combo_V(k_invpow(p))) for p in (1, 6)]
    EXTRA += [("mid+axis combo gauss1.0", combo_V(k_gauss(1.0)))]
    for nm, Vx in EXTRA:
        a_s, v_s = path_min_cfg(Vx)
        try:
            F = aligned_frame(a_s)
            pr = Probe(F)
            verd, asym, ratio = smoothness(pr, Vx, on_config=True)
            g, _, _ = pr.grad(Vx, 1e-3, on_config=True)
            H = pr.hess(Vx, 1e-3, on_config=True)
            (p, z, n), ev = inertia(H)
            tail = f"{verd} r={ratio:.3f}"
        except AssertionError:
            tail, p, z, n, g, ev = "chart-fail", 0, 0, 0, np.zeros(6), np.zeros(6)
        print(f"  {nm:28s} {a_s:14.8f} {v_s:16.8f}  {name_angle(a_s, 1e-3):16s}"
              f" {tail:>16s}")
        print(f"       |grad|={np.linalg.norm(g):.2e}  inertia=({p},{z},{n})  "
              f"eig = {np.array2string(ev, precision=6)}")
    print("  The axis kernels use |u_i . u_j|, which has a corner where two struts")
    print("  are exactly PERPENDICULAR -- and the pre-check shows that happens on")
    print("  this family. So a CORNER verdict here is expected and is NOT evidence")
    print("  about witness selection; it is the absolute value in the definition.")

    print("\n" + "=" * 78)
    print("PERTURBATION: do the strut-kernel minima move?")
    print("=" * 78)
    print(f"  {'kernel':20s} {'var':5s} {'folded argmin at grid 0.5 / 0.1 / 0.037':>44s}"
          f"  {'stable?':>8s}")
    for kn, kf in KERNELS:
        for norm in (False, True):
            Vs = mid_V(kf, norm)
            outs = [fund(path_min_cfg(Vs, np.arange(-179.9, 180.0, st))[0])
                    for st in (0.5, 0.1, 0.037)]
            spread = max(outs) - min(outs)
            print(f"  {kn:20s} {'norm' if norm else 'raw':5s} "
                  f"{outs[0]:14.7f} {outs[1]:14.7f} {outs[2]:14.7f}  "
                  f"{'YES' if spread < 1e-3 else 'NO -- moves %.3f deg' % spread:>8s}")
    print("  Folding first is what makes this test meaningful: mirror-image jumps are")
    print("  the measured a -> 180-a orbit and are NOT instability, but a folded")
    print("  spread that survives IS. Any 'NO' row's minimum is not a located number.")

    print("\n  WHY the surviving 'NO' rows move: an EXACT tie at two angles that the")
    print("  a->-a and a->180-a symmetries do NOT relate. Measured, unexplained:")
    print(f"  {'kernel':16s} {'a1':>12s} {'V(a1)':>22s} {'a2':>12s} {'V(a2)':>22s}"
          f" {'rel diff':>10s}")
    for nm, kf, a1, a2 in (("gauss s=0.5", k_gauss(0.5), 0.0, 79.1066058),
                           ("gauss s=1.0", k_gauss(1.0), 0.0, 79.1066280),
                           ("1/r^3", k_invpow(3.0), 18.7687268, 81.7385673)):
        Vx = mid_V(kf, True)
        v1, v2 = Vx(corners(a1)), Vx(corners(a2))
        print(f"  {nm:16s} {a1:12.7f} {v1:22.15f} {a2:12.7f} {v2:22.15f} "
              f"{abs(v1-v2)/abs(v1):10.2e}")
    print("  Agreement at 1-2 ULP is a tie, not a coincidence of two separate minima.")
    print("  But the sorted-distance signatures of the two normalised midpoint sets")
    print("  differ by ~5e-8, so an EXACT isometry between them was NOT demonstrated.")
    print("  Recorded as an open structural question. The consequence for K4 stands")
    print("  either way: the normalised strut ground state is at least TWO-FOLD")
    print("  degenerate, so 'the strut kernel minimises at 17.16' is an incomplete")
    print("  statement of where its ground state is.")

    print("\n" + "=" * 78)
    print("K4 VERDICT INPUT: vertex vs strut ground state, side by side")
    print("=" * 78)
    print("  Both columns FOLDED to the fundamental domain [0,90], so a difference")
    print("  here is a real difference of state and not a mirror image.")
    print(f"  {'kernel':20s} {'variant':8s} {'VERTEX (folded)':>16s} "
          f"{'STRUT (folded)':>16s}  {'gap':>9s}  same state?")
    for kn, norm, Vs, a_s, v_s, a_v in rows:
        fv, fs = fund(a_v), fund(a_s)
        print(f"  {kn:20s} {'norm' if norm else 'raw':8s} {fv:16.7f} "
              f"{fs:16.7f}  {abs(fv - fs):9.4f}  "
              f"{'YES' if abs(fv - fs) < 1e-3 else 'NO'}")
    print("\n  a_ico folded = %.7f, VE folded = 0, a=90 folded = 90" % fund(A_ICO))
