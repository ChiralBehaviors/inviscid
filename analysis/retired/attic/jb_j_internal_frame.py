"""Step J: a 6-dimensional local chart on the linkage variety.

Everything in M2-M5 needs the same thing: a way to move off a configuration in
the SIX internal directions (R2: 48 body DOF - 36 hinge constraints - 6 global
rigid motions) while staying exactly on the variety. jb_b measured the tangent
space; this builds a chart.

Parameterisation. Body i carries a rotation vector w_i about its own BASE
centroid and a translation t_i, so y = (w_0..w_7, t_0..t_7) in R^48 and

    X_i(y)_j = R(w_i) @ (X0[i,j] - cen0[i]) + cen0[i] + t_i

is exactly rigid for every y, by construction -- no rigidity residual to chase.

Gauge. The 6 global rigid motions are removed by an explicit LINEAR side
condition Rg^T y = 0, with Rg the 48x6 matrix of global modes at the base
configuration (global rotation omega: w_i = omega, t_i = omega x cen_i; global
translation tau: w_i = 0, t_i = tau). This matters: a "zero eigenvalue" that is
really a leaked rigid mode would be a BUG masquerading as physics, so the gauge
is fixed by construction and then re-checked numerically.

Chart. Let B (48x6) be an orthonormal basis of ker[DC(0); Rg^T]. For s in R^6,
Newton-solve

    C(y) = 0,   Rg^T y = 0,   B^T y = s

for y(s). y(0) = 0 and dy/ds(0) = B, so s is a genuine local coordinate whose
tangent at the origin is the internal tangent space. Any function of the
configuration becomes a function of s, and at a CRITICAL POINT its Hessian in s
has retraction-independent inertia (Sylvester).
"""
import numpy as np

from jb_a_family import corners, rot
from jb_b_variety import PAIRS, jacobian


def rodrigues(w):
    """Exact rotation matrix from a rotation vector (no small-angle assumption)."""
    th = np.linalg.norm(w)
    if th < 1e-14:
        K = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
        return np.eye(3) + K + 0.5 * (K @ K)
    u = w / th
    return rot(u, np.degrees(th))


class Frame:
    """A 6-D chart on the variety, anchored at the configuration `corners(a0)`."""

    def __init__(self, a0, X0=None, tol=1e-9):
        self.a0 = a0
        self.X0 = corners(a0) if X0 is None else np.array(X0, dtype=float)
        self.cen0 = self.X0.mean(axis=1)
        self.Rg = self._global_modes()
        J = jacobian(self.X0)
        A = np.vstack([J, self.Rg.T])                    # 42 x 48
        U, S, Vt = np.linalg.svd(A)
        self.sv = S
        ns = Vt[np.sum(S > tol * max(1.0, S[0])):]       # rows spanning the kernel
        self.B = ns.T                                    # 48 x k
        self.dim = self.B.shape[1]
        self.A_lin = A

    def _global_modes(self):
        Rg = np.zeros((48, 6))
        for k in range(3):
            om = np.eye(3)[k]
            for i in range(8):
                Rg[3 * i:3 * i + 3, k] = om
                Rg[24 + 3 * i:24 + 3 * i + 3, k] = np.cross(om, self.cen0[i])
            Rg[24 + 0::3, 3 + k] = 0.0
        for k in range(3):
            for i in range(8):
                Rg[24 + 3 * i + k, 3 + k] = 1.0
        # orthonormalise so the gauge is well conditioned
        Q, _ = np.linalg.qr(Rg)
        return Q

    def config(self, y):
        X = np.empty((8, 3, 3))
        for i in range(8):
            R = rodrigues(y[3 * i:3 * i + 3])
            X[i] = (R @ (self.X0[i] - self.cen0[i]).T).T + self.cen0[i] \
                + y[24 + 3 * i:24 + 3 * i + 3]
        return X

    def residual(self, y):
        X = self.config(y)
        return np.concatenate([X[i, j] - X[k, l] for (i, j), (k, l) in PAIRS])

    def solve(self, s, y0=None, iters=60, tol=1e-13):
        """Newton-solve for the configuration at chart coordinate s."""
        s = np.atleast_1d(np.asarray(s, dtype=float))
        y = np.zeros(48) if y0 is None else y0.copy()
        for _ in range(iters):
            F = np.concatenate([self.residual(y), self.Rg.T @ y, self.B.T @ y - s])
            if np.linalg.norm(F) < tol:
                break
            X = self.config(y)
            DC = jacobian(X)      # exact only at y=0; a good Newton matrix elsewhere
            Jf = np.vstack([DC, self.Rg.T, self.B.T])
            dy, *_ = np.linalg.lstsq(Jf, -F, rcond=None)
            step = 1.0
            for _ in range(30):
                yn = y + step * dy
                Fn = np.concatenate([self.residual(yn), self.Rg.T @ yn,
                                     self.B.T @ yn - s])
                if np.linalg.norm(Fn) < np.linalg.norm(F):
                    break
                step *= 0.5
            y = y + step * dy
        return y, np.linalg.norm(np.concatenate(
            [self.residual(y), self.Rg.T @ y, self.B.T @ y - s]))

    def vertices(self, y):
        X = self.config(y)
        return np.array([X[i, j] for (i, j), _ in PAIRS])

    # ---- differentiation of a scalar f(vertices) in the chart ----

    def value(self, f, s, cache=None):
        key = tuple(np.round(s, 15))
        if cache is not None and key in cache:
            return cache[key]
        y, res = self.solve(s)
        assert res < 1e-9, f"Newton failed at s={s}: residual {res:.2e}"
        v = f(self.vertices(y))
        if cache is not None:
            cache[key] = v
        return v

    def grad(self, f, h=1e-4, onesided=False):
        n = self.dim
        g = np.zeros(n)
        gp, gm = np.zeros(n), np.zeros(n)
        f0 = self.value(f, np.zeros(n))
        for i in range(n):
            e = np.zeros(n)
            e[i] = h
            fp, fm = self.value(f, e), self.value(f, -e)
            gp[i] = (fp - f0) / h
            gm[i] = (f0 - fm) / h
            g[i] = (fp - fm) / (2 * h)
        return (g, gp, gm) if onesided else g

    def hessian(self, f, h=1e-4):
        n = self.dim
        cache = {}
        f0 = self.value(f, np.zeros(n), cache)
        H = np.zeros((n, n))
        for i in range(n):
            e = np.zeros(n)
            e[i] = h
            H[i, i] = (self.value(f, e, cache) - 2 * f0
                       + self.value(f, -e, cache)) / h ** 2
        for i in range(n):
            for j in range(i + 1, n):
                ei, ej = np.zeros(n), np.zeros(n)
                ei[i] = h
                ej[j] = h
                v = (self.value(f, ei + ej, cache) - self.value(f, ei - ej, cache)
                     - self.value(f, -ei + ej, cache)
                     + self.value(f, -ei - ej, cache)) / (4 * h ** 2)
                H[i, j] = H[j, i] = v
        return 0.5 * (H + H.T)


def inertia(H, tol=None):
    """(n_pos, n_zero, n_neg) of a symmetric matrix, with an explicit tolerance."""
    ev = np.linalg.eigvalsh(H)
    if tol is None:
        tol = 1e-6 * max(1.0, np.abs(ev).max())
    return (int(np.sum(ev > tol)), int(np.sum(np.abs(ev) <= tol)),
            int(np.sum(ev < -tol))), ev


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True, linewidth=150)
    print("=== chart sanity: dimension, gauge, and Newton ===")
    for a0 in (0.0, 7.218951362, 22.238756093, 30.0, 60.0, 90.0):
        F = Frame(a0)
        print(f"  a0={a0:12.7f}  chart dim = {F.dim}   "
              f"sv[41]={F.sv[41]:.3e}  sv[35]={F.sv[35]:.3e}")

    F = Frame(30.0)
    print("\n=== does the chart actually stay ON the variety? ===")
    rng = np.random.default_rng(7)
    for mag in (1e-4, 1e-3, 1e-2, 1e-1):
        s = mag * rng.standard_normal(F.dim)
        y, res = F.solve(s)
        X = F.config(y)
        hinge = max(np.linalg.norm(X[i, j] - X[k, l]) for (i, j), (k, l) in PAIRS)
        # rigidity is exact by construction; check it anyway (could fail if the
        # rodrigues/rotation bookkeeping were wrong)
        from jb_a_family import L_EDGE
        el = np.array([np.linalg.norm(X[i, j] - X[i, (j + 1) % 3])
                       for i in range(8) for j in range(3)])
        print(f"  |s|={mag:.0e}  newton residual {res:.2e}   max hinge gap "
              f"{hinge:.2e}   strut len dev {np.abs(el - L_EDGE).max():.2e}")

    print("\n=== is the symmetric-path direction inside the chart's tangent space? ===")
    from jb_b_variety import path_tangent
    for a0 in (0.0, 7.218951362, 22.238756093, 30.0):
        F = Frame(a0)
        xi = path_tangent(a0)
        xi = xi / np.linalg.norm(xi)
        proj = F.B @ (F.B.T @ xi)
        print(f"  a0={a0:12.7f}  ||xi - P_B xi|| = {np.linalg.norm(xi - proj):.3e}"
              f"   (0 => the symmetric tangent is one of the 6)")

    print("\n=== gauge check: are the 6 global modes really excluded? ===")
    F = Frame(30.0)
    print(f"  max |Rg^T B| = {np.abs(F.Rg.T @ F.B).max():.3e}  (0 => no rigid leak)")
    print("  a rigid motion moves no vertex RELATIVE to the others; if one leaked")
    print("  into B it would show up as a spurious zero eigenvalue of every Hessian.")
