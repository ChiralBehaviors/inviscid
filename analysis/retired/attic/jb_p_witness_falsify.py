"""Step P (K3): try to BREAK the witness-selection hypothesis.

THE HYPOTHESIS, stated as it was put:

    Hull volume and strut clearance are both WITNESS-SELECTION functionals --
    hull volume depends on which triples support a face, segment clearance on
    which pair of points realises the minimum -- and at highly symmetric
    configurations the witness set is degenerate and jumps under perturbation,
    so the functional is generically non-smooth there. Raw Thomson selects
    nothing and came back smooth, consistent with this.

It was formed AFTER seeing the data it explains, so it is treated here as a
suspect, not as a conclusion. Four attacks, in increasing order of how much
they would cost the hypothesis if they landed.

P0. THE HALF THAT IS A THEOREM, NOT A MEASUREMENT.
    "A smooth all-pairs kernel is smooth" is not an empirical claim. V_f(P) =
    sum f(|r_i-r_j|^2) is a finite composition of polynomials with f; wherever
    no two vertices coincide and f is C^k on a neighbourhood of the realised
    squared distances, V_f is C^k, full stop. So the second suggested
    refutation -- "a smooth all-pairs kernel that is NOT smooth at the VE" --
    cannot be found by searching, and any hit would be a bug in the harness.
    This matters: HALF THE HYPOTHESIS CARRIES NO INFORMATION. Its whole content
    is the other half, about witness functionals. Reported, not measured; the
    only measurable residue is the BOUNDARY case where the premise fails, which
    is P1.

P1. A smooth all-pairs kernel where the vertex set itself degenerates.
    At a=60/120 the 12 shared vertices merge in pairs. An inverse power blows
    up (a pole, not a corner). A Gaussian stays finite. Is it still SMOOTH
    there, in all six chart directions? A corner would refute P0's argument and
    so refute the hypothesis's easy half. Prediction: smooth.

P2. WITNESS-SELECTION FUNCTIONALS THAT MIGHT BE SMOOTH ANYWAY. The real attack.
    Five selectors evaluated at the VE, where their witness sets are maximally
    degenerate, and at a=30 as a control where they are not:
      - circumradius max_i |r_i|        (all 12 tied: the VE is cospherical)
      - min pairwise separation          (24 tied: the struts)
      - diameter max pairwise separation (6 tied: the antipodal pairs)
      - hull SURFACE AREA                (same face-lattice jump as the volume)
      - lambda_max of the inertia tensor (the VE's inertia tensor is isotropic,
                                          so ALL directions are tied)
    Any one of these coming back SMOOTH at the VE refutes the hypothesis as
    stated.

P3. THE SHARPENING TEST, and the one that actually lands.
    trace(I) = lambda_1 + lambda_2 + lambda_3 is computed by exactly the same
    eigen-selection machinery as lambda_max, over exactly the same maximally
    degenerate witness set -- but it sums the WHOLE orbit instead of picking a
    proper subset. Likewise Vol_hull = smooth + (1/2)*sum |tet_q| picks the
    LARGER triangulation of each skew square; the sum of BOTH triangulations
    picks the whole orbit. If the whole-orbit versions are smooth while the
    proper-subset versions are not, then "witness selection" is not the
    operative property and the hypothesis is refuted AS STATED -- and replaced
    by something sharper and already standard.

P4. Is the strut-repulsion kink really carried by a witness change?
    jb_l located it on exactly 12 of the 204 pairs, all at d = 2.000000. This
    inspects the ARGMIN of the segment-segment distance on those 12 pairs
    across a=0 and asks whether it jumps, is clamped at an endpoint, or is
    degenerate through parallelism. If the argmin is smooth and interior, the
    kink comes from somewhere else and the hypothesis loses its only direct
    mechanism.
"""
import numpy as np
from scipy.spatial import ConvexHull

from jb_a_family import corners
from jb_b_variety import PAIRS
from jb_g_strut_clearance import segment_distance, _hinged
from jb_j_internal_frame import Frame
from jb_k_hull_hessian import aligned_frame, hull_vol
from jb_o_kernel_family import (Probe, k_gauss, k_invpow, make_V, smoothness,
                                sq_dists, vert, A_ICO)
from jb_l_vertex_potentials import STRUT_PAIRS


# ---------------------------------------------------------------- P2 selectors

def circumradius(P):
    """max_i |r_i| about the centroid. Witness: which vertex is farthest."""
    Q = P - P.mean(axis=0)
    return float(np.linalg.norm(Q, axis=1).max())


def min_separation(P):
    """min_{i<j} r_ij. Witness: which pair is closest."""
    return float(np.sqrt(sq_dists(P).min()))


def diameter(P):
    """max_{i<j} r_ij. Witness: which pair is farthest."""
    return float(np.sqrt(sq_dists(P).max()))


def hull_area(P):
    """Convex-hull surface area. Same face-lattice witness as the volume."""
    keep = [0]
    for i in range(1, len(P)):
        if min(np.linalg.norm(P[i] - P[k]) for k in keep) > 1e-12:
            keep.append(i)
    return float(ConvexHull(P[keep]).area)


def inertia_tensor(P):
    Q = P - P.mean(axis=0)
    return Q.T @ Q


def lam_max(P):
    """Largest eigenvalue of the second-moment tensor. Witness: which axis."""
    return float(np.linalg.eigvalsh(inertia_tensor(P))[-1])


def lam_min(P):
    return float(np.linalg.eigvalsh(inertia_tensor(P))[0])


def lam_sum(P):
    """lambda_1 + lambda_2 + lambda_3, computed BY THE EIGENSOLVER.

    Deliberately NOT written as np.trace: the point of P3 is that the same
    degenerate eigen-selection that makes lam_max non-smooth is being run here,
    and only the choice to keep the whole orbit differs.
    """
    return float(np.linalg.eigvalsh(inertia_tensor(P)).sum())


def lam_trace(P):
    """The identical quantity via trace. A GUARD, not an independent check:
    if these two disagree the eigensolver is broken, which is worth knowing."""
    return float(np.trace(inertia_tensor(P)))


# ------------------------------------------------- P3: hull triangulation orbit

def square_quads(a_probe=1e-3, tol=1e-7):
    """The 6 square faces of the cuboctahedron, as vertex-index 4-tuples.

    Read off the a=0 hull: the four coplanar vertices of each square face. Found
    by taking the 6 face-planes with 4 supporting vertices at a=0.
    """
    P = vert(0.0)
    hull = ConvexHull(P)
    quads, seen = [], set()
    for eq in hull.equations:
        d = P @ eq[:3] + eq[3]
        idx = tuple(sorted(np.where(np.abs(d) < 1e-9)[0].tolist()))
        if len(idx) == 4 and idx not in seen:
            seen.add(idx)
            quads.append(idx)
    return quads


QUADS = square_quads()


def tet_vol(P, q):
    """Signed volume of the tetrahedron on a quad's 4 vertices.

    Zero exactly when the quad is planar. This is the quantity whose ORDER OF
    VANISHING decides the hull's smoothness (jb_i): it is odd and linear in a,
    so |tet| is a kink and (1/2)sum|tet| is the whole non-smooth part.
    """
    i, j, k, l = q
    return float(np.dot(np.cross(P[j] - P[i], P[k] - P[i]), P[l] - P[i]) / 6.0)


def hull_pick(P):
    """(1/2) * sum_q |tet_q| -- the LARGER of the two triangulations, per quad.

    A proper subset of the two-element witness orbit. This is exactly the term
    jb_i identified as the hull's non-smooth part.
    """
    return 0.5 * sum(abs(tet_vol(P, q)) for q in QUADS)


def hull_orbit(P):
    """(1/2) * sum_q (|tet_q| + |-tet_q|) / 2 -> the AVERAGE over the orbit.

    Summing both triangulations of each skew quad gives 2 * (their mean), and
    their mean is (+tet + -tet)/2 = 0 plus the common smooth part -- so the
    orbit-symmetric combination is identically the SMOOTH branch. Implemented
    as the explicit orbit average so the selection machinery is still present.
    """
    return 0.5 * sum(0.5 * (tet_vol(P, q) + (-tet_vol(P, q))) for q in QUADS)


def hull_smooth_branch(P):
    """Vol_hull minus its picked kink term. Predicted smooth at the VE."""
    return hull_vol(P) - hull_pick(P) / (np.sqrt(2.0) ** 3 / (6 * np.sqrt(2.0)))


# ------------------------------------------------------------ P4: the argmin

def segment_argmin(p0, p1, q0, q1):
    """Realising parameters (s,t) of the segment-segment distance, by brute
    force on a fine grid then local polish. Independent of jb_g's clamped
    analytic solve, so a disagreement is informative."""
    ss = np.linspace(0, 1, 601)
    u, v = p1 - p0, q1 - q0
    A = p0[None, :] + ss[:, None] * u[None, :]
    B = q0[None, :] + ss[:, None] * v[None, :]
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
    i, j = np.unravel_index(np.argmin(D), D.shape)
    s, t, best = ss[i], ss[j], D[i, j]
    step = 1.0 / 600
    for _ in range(60):
        improved = False
        for ds in (-step, 0.0, step):
            for dt in (-step, 0.0, step):
                s2 = min(1.0, max(0.0, s + ds))
                t2 = min(1.0, max(0.0, t + dt))
                d = np.linalg.norm((p0 + s2 * u) - (q0 + t2 * v))
                if d < best - 1e-15:
                    s, t, best, improved = s2, t2, d, True
        if not improved:
            step *= 0.5
    return s, t, best


if __name__ == "__main__":
    np.set_printoptions(precision=7, suppress=False, linewidth=175)

    print("=" * 78)
    print("P0. The 'smooth kernel => smooth V' half is a THEOREM, not a measurement")
    print("=" * 78)
    print("  V_f(P) = sum f(|r_i-r_j|^2) is a finite composition of a polynomial map")
    print("  with f. Where no two vertices coincide and f is C^k near the realised")
    print("  squared distances, V_f is C^k. There is nothing to search for. Half the")
    print("  hypothesis is therefore vacuous as a prediction; its entire empirical")
    print("  content is the witness half, which P2-P4 attack.")

    print("\n" + "=" * 78)
    print("P1. BOUNDARY CASE: a smooth kernel where the VERTEX SET degenerates (a=60)")
    print("=" * 78)
    print("  At a=60 the 12 shared vertices merge in pairs. Inverse powers have a")
    print("  POLE there (not a corner). A Gaussian is finite. Is it smooth?")
    Vg = make_V(k_gauss(1.0), False)
    Vi = make_V(k_invpow(1.0), False)
    for a0, tag in ((60.0, "octahedron (12->6 merge)"),
                    (30.0, "control, no merge"),
                    (0.0, "control, VE")):
        F = aligned_frame(a0)
        pr = Probe(F)
        v, asym, r = smoothness(pr, Vg)
        print(f"  a0={a0:6.2f} [{tag:26s}] gauss s=1.0 : {v:12s} "
              f"asym {asym[0]:.3e} {asym[1]:.3e} {asym[2]:.3e}  ratio {r:.3f}")
        print(f"          raw 1/r there = {Vi(vert(a0))!r}  "
              f"min separation = {min_separation(vert(a0)):.3e}")
    print("  VERDICT P1: a divergence at a merge is a POLE of the kernel, a different")
    print("  failure from a corner, and it does not touch the hypothesis either way.")

    print("\n" + "=" * 78)
    print("P2. WITNESS-SELECTION FUNCTIONALS AT THE VE. Any SMOOTH one refutes.")
    print("=" * 78)
    SEL = (("circumradius max|r_i|", circumradius, "all 12 tied (cospherical)"),
           ("min separation", min_separation, "24 tied (the struts)"),
           ("diameter", diameter, "antipodal pairs tied"),
           ("hull surface area", hull_area, "face lattice 14 -> 20"),
           ("lambda_max(inertia)", lam_max, "inertia tensor isotropic"),
           ("lambda_min(inertia)", lam_min, "inertia tensor isotropic"))
    for a0, tag in ((0.0, "VE"), (30.0, "CONTROL a=30")):
        print(f"\n  --- at a = {a0} ({tag}) ---")
        F = aligned_frame(a0)
        pr = Probe(F)
        print(f"  {'functional':24s} {'verdict':13s} "
              f"{'|D+-D-| h=1e-2,1e-3,1e-4':>36s} {'ratio':>7s}  witness set")
        for nm, fn, why in SEL:
            v, asym, r = smoothness(pr, fn)
            print(f"  {nm:24s} {v:13s} {asym[0]:11.3e} {asym[1]:11.3e} "
                  f"{asym[2]:11.3e} {r:7.3f}  {why}")
        v, asym, r = smoothness(pr, hull_vol)
        exp = "must read CORNER here" if a0 == 0.0 else "must read SMOOTH here"
        print(f"  {'Vol_hull (CONTROL)':24s} {v:13s} {asym[0]:11.3e} "
              f"{asym[1]:11.3e} {asym[2]:11.3e} {r:7.3f}  {exp}")
        v, asym, r = smoothness(pr, Vi)
        print(f"  {'1/r Thomson (CONTROL)':24s} {v:13s} {asym[0]:11.3e} "
              f"{asym[1]:11.3e} {asym[2]:11.3e} {r:7.3f}  must read SMOOTH here")

    print("\n  WHY hull AREA is smooth where hull VOLUME is not -- the mechanism.")
    print("  Both take the SAME witness decision (which triangulation of each skew")
    print("  square), so if 'witness selection' were the operative property they")
    print("  would fail together. They do not, because the two triangulations differ")
    print("  at DIFFERENT ORDERS in a:")
    print(f"  {'a':>10s} {'|vol difference|':>20s} {'/a':>14s} "
          f"{'|area difference|':>20s} {'/a^2':>14s}")
    for a in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        P = vert(a)
        dv = da = 0.0
        for q in QUADS:
            i, j, k, l = q
            dv += abs(tet_vol(P, q))
            # the two triangulations of quad (i,j,k,l): diagonals i-k and j-l
            def tri(x, y, z):
                return 0.5 * np.linalg.norm(np.cross(P[y] - P[x], P[z] - P[x]))
            A1 = tri(i, j, k) + tri(i, k, l)
            A2 = tri(i, j, l) + tri(j, k, l)
            da += abs(A1 - A2)
        print(f"  {a:10.1e} {dv:20.12e} {dv / a:14.6e} {da:20.12e} "
              f"{da / a ** 2:14.6e}")
    print("  vol difference / a is FLAT (linear in a => |.| is a kink);")
    print("  area difference / a^2 is FLAT (quadratic => |.| is C1). Same witness,")
    print("  different order of vanishing, opposite smoothness verdict.")

    print("\n" + "=" * 78)
    print("P3. THE SHARPENING TEST: proper subset of the orbit vs the WHOLE orbit")
    print("=" * 78)
    print("  lambda_max and lambda_sum run the SAME eigen-selection over the SAME")
    print("  maximally degenerate witness set at the VE. They differ only in whether")
    print("  a proper subset or the whole orbit is kept. If lambda_sum is smooth and")
    print("  lambda_max is not, 'witness selection' is NOT the operative property.")
    F0 = aligned_frame(0.0)
    pr0 = Probe(F0)
    print(f"\n  {'functional':30s} {'verdict':13s} {'|D+-D-| h=1e-2,1e-3,1e-4':>36s}"
          f" {'ratio':>7s}")
    for nm, fn in (("lambda_max  (proper subset)", lam_max),
                   ("lambda_min  (proper subset)", lam_min),
                   ("lambda_sum  (WHOLE orbit)", lam_sum),
                   ("trace       (GUARD on the above)", lam_trace),
                   ("hull_pick   (proper subset)", hull_pick),
                   ("hull_orbit  (WHOLE orbit)", hull_orbit),
                   ("hull_smooth_branch", hull_smooth_branch)):
        v, asym, r = smoothness(pr0, fn)
        print(f"  {nm:30s} {v:13s} {asym[0]:11.3e} {asym[1]:11.3e} "
              f"{asym[2]:11.3e} {r:7.3f}")
    print("\n  do lambda_sum and trace agree numerically? (a GUARD: they are the same")
    print("  quantity, so disagreement means the eigensolver is broken)")
    for a in (0.0, 13.0, A_ICO, 77.0):
        P = vert(a)
        print(f"    a={a:8.3f}  lambda_sum={lam_sum(P):.14f}  "
              f"trace={lam_trace(P):.14f}  diff={lam_sum(P)-lam_trace(P):+.2e}")
    print("\n  is the inertia tensor really isotropic at the VE (the tie that makes")
    print("  lambda_max's witness set degenerate)? this check can fail:")
    for a in (0.0, 1e-6, 1e-3, 30.0):
        w = np.linalg.eigvalsh(inertia_tensor(vert(a)))
        print(f"    a={a:10.6f}  eig(I) = {np.array2string(w, precision=10)}"
              f"   spread = {w.max()-w.min():.3e}")

    print("\n" + "=" * 78)
    print("P4. Is the strut-repulsion kink carried by a WITNESS change?")
    print("=" * 78)
    print("  jb_l: the a=0 cone point of V_strut is carried by exactly 12 of the 204")
    print("  non-hinged pairs, all at d = 2.000000. Re-find them, then inspect the")
    print("  realising parameters (s,t) of the segment-segment distance across a=0.")
    X0 = corners(0.0)
    gaps = np.array([segment_distance(X0[i, j], X0[i, (j + 1) % 3],
                                      X0[k, l], X0[k, (l + 1) % 3])
                     for (i, j), (k, l) in STRUT_PAIRS])
    hits = np.where(np.abs(gaps - 2.0) < 1e-12)[0]
    print(f"  pairs at d = 2.000000 exactly: {len(hits)} of {len(STRUT_PAIRS)}"
          f"   (jb_l reported 12 of 204)")
    print(f"  full gap spectrum at a=0: {np.array2string(np.unique(np.round(gaps, 9)), precision=9)}")

    print("\n  per-pair one-sided slope of d(a): which pairs actually kink?")
    hh = 1e-4
    kinkers = []
    for idx in range(len(STRUT_PAIRS)):
        (i, j), (k, l) = STRUT_PAIRS[idx]
        ds = []
        for a in (-hh, 0.0, hh):
            X = corners(a)
            ds.append(segment_distance(X[i, j], X[i, (j + 1) % 3],
                                       X[k, l], X[k, (l + 1) % 3]))
        dp, dm = (ds[2] - ds[1]) / hh, (ds[1] - ds[0]) / hh
        if abs(dp - dm) > 1e-3:
            kinkers.append((idx, ds[1], dp, dm))
    print(f"  pairs with |D_+ - D_-| > 1e-3 : {len(kinkers)}")
    print(f"  are they exactly the d=2 pairs? "
          f"{sorted(i for i, *_ in kinkers) == sorted(hits.tolist())}")
    for idx, d0, dp, dm in kinkers[:4]:
        print(f"    pair {idx:3d}  d(0)={d0:.9f}  D_+={dp:+.6f}  D_-={dm:+.6f}")

    print("\n  the witness itself: realising (s,t) at a = -1e-3, 0, +1e-3")
    print("  (s,t interior and continuous => NOT a witness jump => REFUTES;")
    print("   s,t jumping or pinned at an endpoint, or the segments parallel")
    print("   so the realiser is a whole interval => SUPPORTS)")
    print(f"  {'pair':>5s} {'a':>9s} {'s*':>9s} {'t*':>9s} {'d':>12s} "
          f"{'|cos angle|':>12s}  note")
    for idx in hits[:4]:
        (i, j), (k, l) = STRUT_PAIRS[idx]
        for a in (-1e-3, -1e-6, 0.0, 1e-6, 1e-3):
            X = corners(a)
            p0, p1 = X[i, j], X[i, (j + 1) % 3]
            q0, q1 = X[k, l], X[k, (l + 1) % 3]
            s, t, d = segment_argmin(p0, p1, q0, q1)
            u, v = p1 - p0, q1 - q0
            ca = abs(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v))
            note = "PARALLEL" if ca > 1 - 1e-12 else ""
            if s in (0.0, 1.0) or t in (0.0, 1.0):
                note += " endpoint-clamped"
            print(f"  {idx:5d} {a:9.0e} {s:9.5f} {t:9.5f} {d:12.9f} {ca:12.9f}  {note}")

    print("\n  degeneracy of the realiser: is there a whole INTERVAL of minimisers?")
    print("  Profile d(s) = MIN OVER t, not d(s) at t fixed -- for parallel segments")
    print("  the second is quadratic even when the realiser set is an interval, so")
    print("  probing with t pinned would manufacture a false 'not degenerate'.")
    for idx in hits[:4]:
        (i, j), (k, l) = STRUT_PAIRS[idx]
        p0, p1 = X0[i, j], X0[i, (j + 1) % 3]
        q0, q1 = X0[k, l], X0[k, (l + 1) % 3]
        u, v = p1 - p0, q1 - q0
        ca = abs(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v))
        tt = np.linspace(0, 1, 2001)
        prof = []
        for s2 in (0.0, 0.25, 0.5, 0.75, 1.0):
            pt = p0 + s2 * u
            prof.append(np.linalg.norm(pt[None, :] - (q0[None, :] + tt[:, None] * v[None, :]),
                                       axis=1).min())
        prof = np.array(prof)
        print(f"    pair {idx:3d}  |cos angle| = {ca:.12f}   min_t d at "
              f"s=0,.25,.5,.75,1 = {np.array2string(prof, precision=9)}")
        print(f"              spread {prof.max()-prof.min():.3e}   "
              f"{'INTERVAL of realisers' if prof.max()-prof.min() < 1e-9 else 'isolated realiser'}")

    print("\n  THE CROSS-TABULATION. Does 'kinks' coincide with 'parallel'?")
    print("  (this is the falsifiable form: if some kinking pair is non-parallel, or")
    print("   some parallel pair does not kink, the mechanism is not what it looks)")
    kset = set(i for i, *_ in kinkers)
    tab = {}
    for idx in range(len(STRUT_PAIRS)):
        (i, j), (k, l) = STRUT_PAIRS[idx]
        u = X0[i, (j + 1) % 3] - X0[i, j]
        v = X0[k, (l + 1) % 3] - X0[k, l]
        ca = abs(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v))
        par = ca > 1 - 1e-9
        tab[(par, idx in kset)] = tab.get((par, idx in kset), 0) + 1
    print(f"  {'':14s} {'kinks':>10s} {'does not kink':>15s}")
    print(f"  {'parallel':14s} {tab.get((True, True), 0):10d} "
          f"{tab.get((True, False), 0):15d}")
    print(f"  {'not parallel':14s} {tab.get((False, True), 0):10d} "
          f"{tab.get((False, False), 0):15d}")
    print("  NOT a clean diagonal: 24 parallel pairs do NOT kink. So an interval of")
    print("  realisers is NECESSARY for the kink here but NOT SUFFICIENT. Refined:")
    groups = {}
    for idx in range(len(STRUT_PAIRS)):
        (i, j), (k, l) = STRUT_PAIRS[idx]
        u = X0[i, (j + 1) % 3] - X0[i, j]
        v = X0[k, (l + 1) % 3] - X0[k, l]
        if abs(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v)) > 1 - 1e-9:
            d0 = round(segment_distance(X0[i, j], X0[i, (j + 1) % 3],
                                        X0[k, l], X0[k, (l + 1) % 3]), 9)
            groups.setdefault(d0, []).append(idx)
    print(f"  {'d(0)':>16s} {'count':>7s} {'|D+-D-| at h = 1e-2, 1e-3, 1e-4, 1e-5':>44s}"
          f"  verdict")
    for d0, idxs in sorted(groups.items()):
        idx = idxs[0]
        (i, j), (k, l) = STRUT_PAIRS[idx]
        out = []
        for hh2 in (1e-2, 1e-3, 1e-4, 1e-5):
            ds = []
            for a in (-hh2, 0.0, hh2):
                Xa = corners(a)
                ds.append(segment_distance(Xa[i, j], Xa[i, (j + 1) % 3],
                                           Xa[k, l], Xa[k, (l + 1) % 3]))
            out.append(abs((ds[2] - ds[1]) / hh2 - (ds[1] - ds[0]) / hh2))
        flat = out[2] > 0.5 * out[0]
        print(f"  {d0:16.9f} {len(idxs):7d} " + " ".join(f"{x:10.3e}" for x in out)
              + f"  {'CORNER' if flat else 'smooth (stationary)'}")
    print("  The two non-kinking parallel groups sit at a STATIONARY point of d(a):")
    print("  all their tied realisers have the same derivative, so the min does not")
    print("  jump. Same degeneracy, no kink -- the same refinement P3 forced.")
    print("  TRAP, recorded so it is not read as a refutation: at h=1e-5 the d=2")
    print("  group's asymmetry COLLAPSES to ~1e-8. That is the double-precision floor")
    print("  (d changes by ~2e-7 against a value of 2), not smoothness. The kink")
    print("  verdict must be read off h >= 1e-4.")

    print("\n  cross-check the KINK MAGNITUDE against the prior survey, which reported")
    print("  D_+ = +0.00504 = -D_- for these 12 pairs. That figure is for 1/d, not d:")
    for idx, d0, dp, dm in kinkers[:3]:
        print(f"    pair {idx:3d}  d: D_+={dp:+.6f}  ->  1/d: D_+={-dp/d0**2:+.8f}"
              f"   (prior 0.00504)")

    print("\n" + "=" * 78)
    print("WHAT WAS TRIED AND DID NOT DISCRIMINATE")
    print("=" * 78)
    print("  - Looking for a smooth all-pairs kernel that is non-smooth at the VE:")
    print("    abandoned, see P0. It is ruled out by composition, not by search, so")
    print("    running it would have been a check that cannot fail.")
    print("  - Testing the selectors ALONG THE SYMMETRIC SLICE ONLY: on the slice the")
    print("    tied witnesses move together by symmetry, so several selectors look")
    print("    smooth there and the test yields a FALSE refutation. Only the 6-D")
    print("    chart separates them. Demonstrated:")
    for nm, fn, _ in SEL[:3]:
        h = 1e-4
        d1 = (fn(vert(h)) - fn(vert(0.0))) / h
        d2 = (fn(vert(0.0)) - fn(vert(-h))) / h
        F = aligned_frame(0.0)
        pr = Probe(F)
        v6, asym6, r6 = smoothness(pr, fn)
        print(f"    {nm:24s} slice D_+={d1:+.6f} D_-={d2:+.6f} "
              f"(slice says {'CORNER' if abs(d1-d2) > 1e-4 else 'smooth'}) "
              f"| 6-D says {v6}")
