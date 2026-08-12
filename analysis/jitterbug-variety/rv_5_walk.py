"""REVIEW CHECK 5: does walk() establish that a PATH lies on the variety, or only
that 80 sampled waypoints do?

Instrument the walk to record, per step:
  - the shape-quotient displacement between consecutive on-variety configs
    (a discontinuous jump shows up here and nowhere else)
  - sigma_36 of the Jacobian at each waypoint (if the manifold stays rank-36
    between waypoints, small steps can be joined by short manifold arcs; a dip
    toward 0 means the step may have crossed a singular locus)
  - the pull-solve residual BEFORE projection (how far off-variety the penalty
    solve actually was)

Also: are the three "not reached" targets genuinely unreachable, or did the walk
reach a tetrahedron with a DIFFERENT corner labelling than the fixed index-wise
comparison in walk() can see?
"""
import sys
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners, L_EDGE
from jb_b_variety import jacobian
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import TET, build_tet_config, kabsch
from jb_e_tighten import is_tet, project
from rv_4_continuation import shape_dist


def walk_instrumented(X0, Xt, n=80, pull=1.0, w=1e4):
    z = np.zeros(48)
    rec = []
    for s in np.linspace(1.0 / n, 1.0, n):
        Xprev = place(X0, z)
        Xa = kabsch(Xt.reshape(-1, 3), Xprev.reshape(-1, 3)).reshape(8, 3, 3)
        way = (1 - s) * Xprev + s * Xa
        sol = least_squares(lambda zz: np.concatenate(
            [w * residual(place(X0, zz)), pull * (place(X0, zz) - way).reshape(-1)]),
            z, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        r_before = np.linalg.norm(residual(place(X0, sol.x)))
        zp, r_after = project(X0, sol.x)
        Xnew = place(X0, zp)
        rec.append(dict(step_raw=np.linalg.norm(Xnew - Xprev),
                        step_quot=shape_dist(Xnew, Xprev),
                        r_before=r_before, r_after=r_after,
                        s36=np.linalg.svd(jacobian(Xnew), compute_uv=False)[35]))
        z = zp
    X = place(X0, z)
    Xa = kabsch(Xt.reshape(-1, 3), X.reshape(-1, 3)).reshape(8, 3, 3)
    return X, np.linalg.norm(X - Xa) / np.sqrt(24), rec


def describe(X):
    """Is this configuration a regular tetrahedron, whatever the labelling?"""
    pts = X.reshape(-1, 3)
    uniq = []
    for p in pts:
        if not any(np.linalg.norm(p - q) < 1e-6 for q in uniq):
            uniq.append(p)
    if len(uniq) < 2:
        return f"{len(uniq)} distinct points (degenerate)"
    D = [np.linalg.norm(uniq[i] - uniq[j])
         for i in range(len(uniq)) for j in range(i + 1, len(uniq))]
    return (f"{len(uniq)} distinct corner positions, pairwise distances "
            f"{min(D):.9f}..{max(D):.9f}"
            + ("  == REGULAR TETRAHEDRON" if len(uniq) == 4 and
               max(D) - min(D) < 1e-6 and abs(min(D) - L_EDGE) < 1e-6 else ""))


if __name__ == "__main__":
    which = [int(x) for x in sys.argv[1:]] or [1, 0]
    tets = [s for s in tet_assignments() if is_tet(s)]
    X0 = corners(0.0)
    for t in which:
        Xt = build_tet_config(tets[t])
        X, dist, rec = walk_instrumented(X0, Xt, n=80)
        sq = np.array([r['step_quot'] for r in rec])
        sr = np.array([r['step_raw'] for r in rec])
        s36 = np.array([r['s36'] for r in rec])
        rb = np.array([r['r_before'] for r in rec])
        ra = np.array([r['r_after'] for r in rec])
        print(f"\n===== target {t}   final index-wise RMS to target = {dist:.3e}"
              f"   ({'REACHED' if dist < 1e-6 else 'NOT reached'}) =====")
        print(f"  per-step quotient displacement: median {np.median(sq):.3e}  "
              f"max {sq.max():.3e} at step {int(np.argmax(sq))}  total {sq.sum():.4f}")
        print(f"  per-step raw displacement:      max {sr.max():.3e}   total {sr.sum():.4f}")
        print(f"  jump ratio max/median = {sq.max()/max(np.median(sq),1e-30):.2f}"
              f"   (a discontinuous jump would be orders of magnitude)")
        print(f"  hinge residual BEFORE projection: max {rb.max():.3e}")
        print(f"  hinge residual AFTER  projection: max {ra.max():.3e}")
        print(f"  sigma_36 along the path: min {s36.min():.3e} at step {int(np.argmin(s36))}"
              f"   (rank 36 needs it >> 0)")
        print(f"  last 6 sigma_36: {np.array2string(s36[-6:], precision=3)}")
        print(f"  largest 5 quotient steps: "
              f"{np.array2string(np.sort(sq)[-5:], precision=4)}")
        print(f"  WHAT DID IT ACTUALLY REACH: {describe(X)}")
