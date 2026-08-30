"""Step F: are the 3 unreached tetrahedron targets disconnected, or was the
homotopy just stuck?

A straight-line homotopy is a greedy path and can park on an obstruction. This
retries each target via a RANDOM INTERMEDIATE configuration on the variety, so
the path is no longer forced to be monotone toward the target. Reaching a target
on any restart proves connectivity. Failing every restart proves nothing --
recorded as such, not as disconnection.
"""
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import build_tet_config, kabsch
from jb_e_tighten import is_tet, project, walk

RNG = np.random.default_rng(20260811)


def random_on_variety(X0, scale=0.35):
    """A random configuration reached from X0 by projecting a random perturbation."""
    z = RNG.normal(scale=scale, size=48)
    z, r = project(X0, z)
    return z, r


def walk_from(Xstart, Xt, n=80, pull=1.0, w=1e4):
    """Same as jb_e walk() but from an arbitrary on-variety start configuration."""
    z = np.zeros(48)
    worst = 0.0
    for s in np.linspace(1.0 / n, 1.0, n):
        Xa = kabsch(Xt.reshape(-1, 3), place(Xstart, z).reshape(-1, 3)).reshape(8, 3, 3)
        way = (1 - s) * place(Xstart, z) + s * Xa
        sol = least_squares(lambda zz: np.concatenate(
            [w * residual(place(Xstart, zz)), pull * (place(Xstart, zz) - way).reshape(-1)]),
            z, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        z, r = project(Xstart, sol.x)
        worst = max(worst, r)
    X = place(Xstart, z)
    Xa = kabsch(Xt.reshape(-1, 3), X.reshape(-1, 3)).reshape(8, 3, 3)
    return X, np.linalg.norm(X - Xa) / np.sqrt(24), worst


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    tets = [s for s in tet_assignments() if is_tet(s)]
    X0 = corners(0.0)
    STUCK = [0, 5, 7]
    N_RESTART = 6

    print("=== retrying the 3 unreached targets via random intermediate configurations ===")
    for t in STUCK:
        Xt = build_tet_config(tets[t])
        best = np.inf
        hit_at = None
        for k in range(N_RESTART):
            zi, ri = random_on_variety(X0)
            if ri > 1e-10:
                continue
            Xmid = place(X0, zi)
            X, d, worst = walk_from(Xmid, Xt, n=60)
            best = min(best, d)
            if d < 1e-6:
                hit_at = k
                break
        verdict = (f"REACHED via restart {hit_at}" if hit_at is not None
                   else f"still unreached after {N_RESTART} restarts (best RMS {best:.3e})")
        print(f"  target {t}: {verdict}")

    print("\n=== control: do the same restarts still reach a known-reachable target? ===")
    Xt = build_tet_config(tets[1])
    zi, ri = random_on_variety(X0)
    X, d, worst = walk_from(place(X0, zi), Xt, n=60)
    print(f"  target 1 via random intermediate: RMS {d:.3e}  worst residual {worst:.2e}"
          f"   {'REACHED' if d < 1e-6 else 'NOT reached -- restart machinery itself is suspect'}")
