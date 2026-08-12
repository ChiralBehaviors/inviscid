"""Step N: is the claimed ground state a GLOBAL minimum, or only a slice minimum?

M3 showed a = +/-7.218951982 is a strict local minimum of -Vol_hull in all six
internal directions. "Ground state" claims more than that: it claims nothing
else on the 6-dimensional variety does better. C5 and M0-M3 only ever looked at
a ONE-dimensional slice of a SIX-dimensional space, so the claim is untested by
construction.

Method: random-restart projected ascent. Start from a configuration on the
variety, build the local 6-D chart (jb_j), step along the chart gradient,
re-anchor the chart at the new configuration, repeat. Every iterate is Newton-
projected back onto the variety, so the walk never leaves it. Restarts are
seeded both on and OFF the symmetric path (random chart displacements from
several base angles), so the search is not confined to the sector the answer is
expected in.

DECISION 16: interference is permitted, so the walk is NOT stopped when struts
cross. Excluding those configurations would have pre-decided the answer.

This is a falsification attempt, not a proof. Finding nothing better than
20.4430915 is weak evidence; finding something better would be decisive.
"""
import numpy as np

from jb_a_family import corners, L_EDGE
from jb_b_variety import PAIRS
from jb_j_internal_frame import Frame
from jb_k_hull_hessian import hull_vol
from jb_l_vertex_potentials import thomson, thomson_normalised_centroid

V_TET = L_EDGE ** 3 / (6 * np.sqrt(2.0))


def verts_of(X):
    return np.array([X[i, j] for (i, j), _ in PAIRS])


def walk(X0, f, sense=+1, steps=140, step0=6e-2, seed=0, verbose=False):
    """Projected gradient walk on the variety. sense=+1 ascends f, -1 descends."""
    X = np.array(X0)
    best = f(verts_of(X))
    step = step0
    for it in range(steps):
        F = Frame(0.0, X0=X)
        if F.dim < 6:
            break
        g = np.zeros(F.dim)
        h = 1e-4
        f0 = f(verts_of(F.config(F.solve(np.zeros(F.dim))[0])))
        for i in range(F.dim):
            e = np.zeros(F.dim)
            e[i] = h
            yp, rp = F.solve(e)
            ym, rm = F.solve(-e)
            if rp > 1e-9 or rm > 1e-9:
                return X, best, "newton failure"
            g[i] = (f(verts_of(F.config(yp)))
                    - f(verts_of(F.config(ym)))) / (2 * h)
        n = np.linalg.norm(g)
        if not np.isfinite(n) or n < 1e-10:
            break
        d = sense * g / n
        improved = False
        for _ in range(14):
            y, r = F.solve(step * d)
            if r < 1e-9:
                cand = f(verts_of(F.config(y)))
                if np.isfinite(cand) and sense * (cand - best) > 0:
                    X, best, improved = F.config(y), cand, True
                    break
            step *= 0.5
        if not improved:
            break
        step = min(step * 1.6, step0)
    return X, best, f"stopped at iter {it}"


def seeds(n_off=10, rng=None):
    rng = rng or np.random.default_rng(0)
    out = [("path a=%g" % a, corners(a)) for a in
           (0.0, 3.0, 7.218951982, 15.0, 22.238756093, 40.0, 55.0, 75.0, 90.0,
            110.0, 150.0, 180.0)]
    for k in range(n_off):
        a0 = rng.uniform(0.0, 180.0)
        F = Frame(a0)
        s = rng.uniform(0.05, 0.6) * rng.standard_normal(F.dim)
        y, r = F.solve(s)
        if r < 1e-9:
            out.append((f"off-path from a={a0:.1f} |s|={np.linalg.norm(s):.2f}",
                        F.config(y)))
    return out


if __name__ == "__main__":
    np.set_printoptions(precision=7, suppress=False, linewidth=160)
    rng = np.random.default_rng(2026)
    S = seeds(12, rng)
    print(f"{len(S)} seeds ({sum(1 for n, _ in S if n.startswith('off'))} off the "
          "symmetric path)\n")

    print("=== A. MAXIMISE Vol_hull  (ground state of V = -k*Vol) ===")
    print("  reference: symmetric-path maximum is 20.4430915")
    bestall = (-np.inf, None, None)
    for name, X0 in S:
        X, v, why = walk(X0, hull_vol, sense=+1, seed=1)
        P = verts_of(X)
        r = np.linalg.norm(P, axis=1)
        flag = "  <<< BEATS THE SLICE MAXIMUM" if v > 20.4430915 + 1e-6 else ""
        print(f"  {name:38s} start {hull_vol(verts_of(X0)):9.6f} -> {v:11.7f}"
              f"   r spread {r.max() - r.min():.2e}{flag}")
        if v > bestall[0]:
            bestall = (v, X, name)
    print(f"\n  best found: {bestall[0]:.9f}  from seed '{bestall[2]}'")
    print(f"  slice maximum: 20.443091521   difference: "
          f"{bestall[0] - 20.443091521:+.3e}")
    Pb = verts_of(bestall[1])
    rb = np.linalg.norm(Pb, axis=1)
    print(f"  is the winner ON the symmetric path? radius spread "
          f"{rb.max() - rb.min():.3e} (0 => cospherical, consistent with the path)")

    print("\n=== B. MINIMISE raw Thomson (its claimed ground state is the VE) ===")
    best = (np.inf, None, None)
    for name, X0 in S:
        X, v, why = walk(X0, thomson, sense=-1)
        flag = "  <<< BEATS THE VE" if v < 34.889842063 - 1e-6 else ""
        print(f"  {name:38s} -> {v:12.7f}{flag}")
        if v < best[0]:
            best = (v, X, name)
    print(f"  best found {best[0]:.9f}   VE value 34.889842063   "
          f"difference {best[0] - 34.889842063:+.3e}")

    print("\n=== C. MINIMISE centroid-normalised Thomson (claimed: icosahedron) ===")
    best = (np.inf, None, None)
    for name, X0 in S:
        X, v, why = walk(X0, thomson_normalised_centroid, sense=-1)
        flag = "  <<< BEATS THE ICOSAHEDRON" if v < 49.165253058 - 1e-6 else ""
        print(f"  {name:38s} -> {v:12.7f}{flag}")
        if v < best[0]:
            best = (v, X, name)
    print(f"  best found {best[0]:.9f}   icosahedron 49.165253058   "
          f"difference {best[0] - 49.165253058:+.3e}")
    print("  NOTE: 49.165253058 is the Thomson minimum over ALL 12-point sphere")
    print("  configurations, so it cannot be beaten by anything -- this leg of the")
    print("  search is a consistency check on the walk, not a test of the potential.")

    print("\n=== D. does the walk actually leave the symmetric path? ===")
    print("  (if every walk collapsed onto the path, A-C would prove nothing)")
    for name, X0 in S[:6]:
        X, v, _ = walk(X0, hull_vol, sense=+1, steps=40)
        d = np.linalg.norm(verts_of(X) - verts_of(X0))
        r = np.linalg.norm(verts_of(X), axis=1)
        print(f"  {name:38s} moved {d:.4f} in vertex space, "
              f"radius spread {r.max() - r.min():.2e}")
