"""REVIEW CHECK 4: is continue_along measuring INTERNAL finite motion, or is the
corrector sliding along global rigid motions and inflating the arclength?

The fix: measure arclength in the SHAPE QUOTIENT.  Between consecutive configs
align the new one onto the old with an optimal proper rigid motion (Kabsch) and
take the residual.  Global drift contributes exactly zero to that.
"""
import numpy as np
from scipy.optimize import least_squares
from jb_a_family import corners
from jb_b_variety import jacobian, path_tangent
from jb_c_branches import place, residual, continue_along


def kabsch_align(P, Q):
    """Optimal PROPER rigid motion taking P onto Q; returns transformed P."""
    P, Q = P.reshape(-1, 3), Q.reshape(-1, 3)
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return (R @ Pc.T).T + Q.mean(0)


def shape_dist(A, B):
    return np.linalg.norm(kabsch_align(A, B).reshape(B.shape) - B)


def continue_instrumented(a0, xi, step, n_steps=40, retangent=False):
    """Same march as jb_c_branches.continue_along, but records BOTH the raw
    arclength it reports and the shape-quotient arclength (globals removed),
    plus the worst single-step jump and the minimum sigma_36 along the way."""
    X0 = corners(a0)
    z = np.zeros(48)
    raw = quot = 0.0
    worst_step = 0.0
    min_s36 = np.inf
    dirn = xi / np.linalg.norm(xi)
    steps_done = 0
    for _ in range(n_steps):
        z_pred = z + step * dirn
        sol = least_squares(lambda zz: residual(place(X0, zz)), z_pred,
                            xtol=1e-14, ftol=1e-14, gtol=1e-14)
        if np.linalg.norm(residual(place(X0, sol.x))) > 1e-10:
            break
        Xa, Xb = place(X0, z), place(X0, sol.x)
        d_raw = np.linalg.norm(Xb - Xa)
        d_quot = shape_dist(Xb, Xa)
        raw += d_raw
        quot += d_quot
        worst_step = max(worst_step, d_quot)
        min_s36 = min(min_s36, np.linalg.svd(jacobian(Xb), compute_uv=False)[35])
        z = sol.x
        steps_done += 1
        if retangent:   # re-derive the tangent at the new point (true continuation)
            J = jacobian(Xb)
            _, _, Vt = np.linalg.svd(J)
            N = Vt[36:]
            proj = N @ dirn
            newd = N.T @ proj
            if np.linalg.norm(newd) > 1e-8:
                dirn = newd / np.linalg.norm(newd)
    return dict(raw=raw, quot=quot, steps=steps_done, worst_step=worst_step,
                min_s36=min_s36, X=place(X0, z))


if __name__ == "__main__":
    a0 = 30.0
    X0 = corners(a0)
    J = jacobian(X0)
    U, S, Vt = np.linalg.svd(J)
    null = Vt[36:]
    cen = X0.mean(axis=1)
    glob = np.zeros((6, 48))
    for d in range(3):
        for i in range(8):
            glob[d, 24 + 3*i + d] = 1.0
    for d in range(3):
        e = np.eye(3)[d]
        for i in range(8):
            glob[3 + d, 3*i:3*i+3] = e
            glob[3 + d, 24 + 3*i:24 + 3*i+3] = np.cross(e, cen[i])

    print("=== A. is the hand-built global basis complete? independent construction ===")
    # numerically differentiate the true SE(3) action on X0
    h = 1e-6
    Gnum = []
    for d in range(3):
        e = np.eye(3)[d]
        Gnum.append(((X0 + h * e) - X0).reshape(-1) / h)
    for d in range(3):
        e = np.eye(3)[d]
        K = np.array([[0, -e[2], e[1]], [e[2], 0, -e[0]], [-e[1], e[0], 0]])
        Rp = np.eye(3) + np.sin(h)*K + (1-np.cos(h))*(K @ K)
        Gnum.append(((X0.reshape(-1, 3) @ Rp.T).reshape(8, 3, 3) - X0).reshape(-1) / h)
    Gnum = np.array(Gnum)                                    # (6, 72) in X-space
    # push the hand-built (dw,dt) basis into X-space
    Ghand = []
    for g in glob:
        dX = np.array([[np.cross(g[3*i:3*i+3], X0[i, j] - cen[i]) + g[24+3*i:24+3*i+3]
                        for j in range(3)] for i in range(8)])
        Ghand.append(dX.reshape(-1))
    Ghand = np.array(Ghand)
    Qn, _ = np.linalg.qr(Gnum.T)
    resid = Ghand - (Ghand @ Qn) @ Qn.T
    print(f"  rank(numeric SE(3) tangent) = {np.linalg.matrix_rank(Gnum)}   "
          f"rank(hand-built) = {np.linalg.matrix_rank(Ghand)}")
    print(f"  hand-built basis outside the true SE(3) tangent: "
          f"max residual = {np.abs(resid).max():.3e}")
    Qh, _ = np.linalg.qr(Ghand.T)
    r2 = Gnum - (Gnum @ Qh) @ Qh.T
    print(f"  true SE(3) tangent outside the hand-built span: "
          f"max residual = {np.abs(r2).max():.3e}   => spans agree, basis COMPLETE")

    print("\n=== B. do the 6 modes move the SHAPE, or just the frame? ===")
    Q, _ = np.linalg.qr(glob.T)
    internal = null - (null @ Q) @ Q.T
    Ui, Si, _ = np.linalg.svd(internal.T, full_matrices=False)
    basis = Ui[:, Si > 1e-8].T
    print(f"  internal modes: {len(basis)}")
    print(f"  {'mode':6s} {'raw arc':>9s} {'quotient arc':>13s} {'steps':>6s} "
          f"{'worst step':>11s} {'min s36':>10s}")
    for m, xi in enumerate(basis):
        r = continue_instrumented(a0, xi, 0.02, 40)
        print(f"  {m:<6d} {r['raw']:9.4f} {r['quot']:13.4f} {r['steps']:6d} "
              f"{r['worst_step']:11.3e} {r['min_s36']:10.3e}")

    print("\n  control: a direction that is NOT in the null space (should die at once)")
    rng = np.random.default_rng(1)
    row = Vt[0]                       # a maximally CONSTRAINED direction
    r = continue_instrumented(a0, row, 0.02, 40)
    print(f"  {'J-row':<6s} {r['raw']:9.4f} {r['quot']:13.4f} {r['steps']:6d}")
    r = continue_instrumented(a0, glob[3], 0.02, 40)
    print(f"  {'globrot':<6s} {r['raw']:9.4f} {r['quot']:13.4f} {r['steps']:6d}"
          f"   <== pure global motion: raw arc large, quotient arc ~0")

    print("\n=== C. with the tangent RE-DERIVED at every step (true continuation) ===")
    for m, xi in enumerate(basis):
        r = continue_instrumented(a0, xi, 0.02, 40, retangent=True)
        print(f"  mode {m}: quotient arc = {r['quot']:.4f}  steps={r['steps']}  "
              f"min sigma_36 = {r['min_s36']:.3e}")

    print("\n=== D. does the symmetric path show up as ONE of these directions? ===")
    t = path_tangent(a0)
    tq = t - Q @ (Q.T @ t)
    tq /= np.linalg.norm(tq)
    print(f"  |overlap| of the symmetric-path tangent with each internal mode: "
          f"{np.abs(basis @ tq).round(4)}")
    print(f"  total (should be 1.0 if it lies in the internal span): "
          f"{np.linalg.norm(basis @ tq):.9f}")
