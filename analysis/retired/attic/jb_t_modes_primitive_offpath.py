"""Step T: what the six modes ARE (S3), whether the primitive matters for the
SPECTRUM (S4), and what happens off the equilibrium (S5). Bead inviscid-qvf.2.

S3. Six numbers are not a result until you can say what moves. The 2+3+1
degeneracy means individual eigenvectors inside the doublet and the triplet are
NOT determined -- any orthonormal rotation inside a block is an equally valid
eigenbasis -- so every geometric descriptor here is a TRACE over a block, which
is invariant under that rotation. Reporting a single eigenvector's shape would
be reporting an artefact of the eigensolver's tie-breaking.

S4. DECISION 17 calls the vertex-vs-strut primitive question MOOT because raw
kernels on both share a ground state. Sharing an equilibrium does not imply
sharing a spectrum. This measures whether they do.

S5. Off the equilibrium the Hessian is NOT a tensor: under a reparameterisation
s = psi(t) it picks up (dV)·D^2 psi, which vanishes only where dV = 0. A
"spectrum at the icosahedron" is therefore chart-dependent, and the honest thing
is to MEASURE the dependence rather than quote one chart's answer. Two genuinely
different charts are built for that: the existing one, which pivots each body's
rotation at its own CENTROID, and `OriginFrame`, which pivots at the ORIGIN.
They agree to first order (same tangent space, linearly related coordinates) and
differ at second order, so:

    at the VE (dV = 0) they MUST agree            <- the control that can fail
    at the icosahedron (dV =/= 0) they need not   <- the measurement
"""
import numpy as np
import scipy.linalg as sla

from jb_a_family import corners
from jb_b_variety import PAIRS, cross_row, jacobian, path_tangent
from jb_j_internal_frame import Frame, inertia, rodrigues
from jb_k_hull_hessian import aligned_frame
from jb_o_kernel_family import (KERNELS, Probe, make_V, minimise_on_path, vert,
                                A_ICO)
from jb_q_strut_kernels import mid_V
from jb_r_mass_metric import MODELS, position_jacobian, weights
from jb_s_frequency_spectrum import (LABEL, SHORT, degeneracy_blocks,
                                     fmt_blocks, gen_eig, irrep_blocks,
                                     block_stiffness, RAW_KERNELS)

np.seterr(all="ignore")      # spurious numpy/BLAS matmul warnings; see jb_r


# --------------------------------------------------------------------------
# a SECOND, genuinely different chart: rotation pivoted at the ORIGIN
# --------------------------------------------------------------------------

class OriginFrame(Frame):
    """Same variety, same tangent space, DIFFERENT second-order structure.

    `Frame` rotates body i about its own centroid; this rotates it about the
    origin. The two body parameterisations have linearly related tangent maps
    (a rotation about the origin is a rotation about the centroid plus a
    translation) but the finite maps differ, so the induced charts differ by a
    nonlinear reparameterisation. That is exactly the freedom a Hessian is
    invariant under AT a critical point and not invariant under elsewhere.
    """

    def _global_modes(self):
        Rg = np.zeros((48, 6))
        for k in range(3):
            Rg[3 * np.arange(8) + k, k] = 1.0            # global rotation
            Rg[24 + 3 * np.arange(8) + k, 3 + k] = 1.0   # global translation
        Q, _ = np.linalg.qr(Rg)
        return Q

    def config(self, y):
        X = np.empty((8, 3, 3))
        for i in range(8):
            R = rodrigues(y[3 * i:3 * i + 3])
            X[i] = (R @ self.X0[i].T).T + y[24 + 3 * i:24 + 3 * i + 3]
        return X

    def _jac(self, X):
        """36x48 constraint Jacobian in the ORIGIN-pivot coordinates."""
        J = np.zeros((36, 48))
        for c, ((i, j), (k, l)) in enumerate(PAIRS):
            for d in range(3):
                row = 3 * c + d
                J[row, 3 * i:3 * i + 3] += cross_row(X[i, j], d)
                J[row, 3 * k:3 * k + 3] -= cross_row(X[k, l], d)
                J[row, 24 + 3 * i + d] += 1.0
                J[row, 24 + 3 * k + d] -= 1.0
        return J

    def __init__(self, a0, X0=None, tol=1e-9):
        self.a0 = a0
        self.X0 = corners(a0) if X0 is None else np.array(X0, dtype=float)
        self.cen0 = self.X0.mean(axis=1)
        self.Rg = self._global_modes()
        A = np.vstack([self._jac(self.X0), self.Rg.T])
        U, S, Vt = np.linalg.svd(A)
        self.sv = S
        self.B = Vt[np.sum(S > tol * max(1.0, S[0])):].T
        self.dim = self.B.shape[1]

    def solve(self, s, y0=None, iters=60, tol=1e-13):
        s = np.atleast_1d(np.asarray(s, dtype=float))
        y = np.zeros(48) if y0 is None else y0.copy()
        for _ in range(iters):
            F = np.concatenate([self.residual(y), self.Rg.T @ y, self.B.T @ y - s])
            if np.linalg.norm(F) < tol:
                break
            Jf = np.vstack([self._jac(self.config(y)), self.Rg.T, self.B.T])
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


def origin_position_jacobian(X):
    """72x48 corner-velocity map for the ORIGIN-pivot parameterisation."""
    P = np.zeros((72, 48))
    for i in range(8):
        for j in range(3):
            for d in range(3):
                row = 9 * i + 3 * j + d
                P[row, 3 * i:3 * i + 3] = cross_row(X[i, j], d)
                P[row, 24 + 3 * i + d] = 1.0
    return P


def metric_for(F, model, origin_pivot=False):
    wc, wcen = weights(model)
    P = origin_position_jacobian(F.X0) if origin_pivot else position_jacobian(F.X0)
    Vc = P @ F.B
    Vcen = Vc.reshape(8, 3, 3, -1).mean(axis=1).reshape(24, -1)
    M = Vc.T @ (np.repeat(wc, 3)[:, None] * Vc)
    M += Vcen.T @ (np.repeat(wcen, 3)[:, None] * Vcen)
    return 0.5 * (M + M.T)


# --------------------------------------------------------------------------
# S3 -- mode geometry
# --------------------------------------------------------------------------

def vertex_velocity(F, v):
    """Displacement field of the 12 SHARED vertices for chart direction v.

    Also returns the hinge discrepancy: the two corners identified at a shared
    vertex must move identically, which they do only if B really lies in the
    constraint Jacobian's null space. A non-zero value here would mean the
    'internal' directions tear the linkage apart.
    """
    vel = (position_jacobian(F.X0) @ (F.B @ v)).reshape(8, 3, 3)
    u = np.array([vel[i, j] for (i, j), _ in PAIRS])
    u2 = np.array([vel[k, l] for _, (k, l) in PAIRS])
    return u, float(np.abs(u - u2).max())


def s3_modes(kernel_name="1/r^1  (Thomson)"):
    print("=" * 78)
    print(f"S3  WHAT THE SIX MODES DO, GEOMETRICALLY  ({kernel_name})")
    print("=" * 78)
    print("  Every number below is a TRACE over a degenerate block, hence")
    print("  invariant under the eigensolver's arbitrary choice of basis inside")
    print("  the doublet and the triplet. A per-eigenvector figure would not be.")
    F = aligned_frame(0.0)          # basis direction 0 IS the symmetric-path tangent
    pr = Probe(F)
    kern = dict(RAW_KERNELS)[kernel_name]
    H = pr.hess(make_V(kern, normalised=False), 1e-3)
    P12 = vert(0.0)
    rhat = P12 / np.linalg.norm(P12, axis=1)[:, None]

    for model in MODELS:
        M = metric_for(F, model)
        print(f"\n  --- mass model: {model} ---")
        for mult, mval, Q in irrep_blocks(M):
            h, dev = block_stiffness(H, Q)
            omega = np.sqrt(h / mval)
            # trace invariants over the block
            e0 = np.zeros(F.dim)
            e0[0] = 1.0
            overlap2 = float(np.sum((Q.T @ e0) ** 2))
            tot = rad = 0.0
            breath = 0.0
            breath_abs = 0.0
            part = np.zeros(12)
            spin = 0.0
            hinge = 0.0
            for k in range(Q.shape[1]):
                u, hg = vertex_velocity(F, Q[:, k])
                hinge = max(hinge, hg)
                tot += float(np.sum(u * u))
                rad += float(np.sum(np.sum(u * rhat, axis=1) ** 2))
                per_vertex = np.sum(P12 * u, axis=1)
                breath += float(per_vertex.sum()) ** 2
                breath_abs += float(np.abs(per_vertex).sum())
                part += np.sum(u * u, axis=1)
                yv = F.B @ Q[:, k]
                spin += float(yv[:24] @ yv[:24]) / float(yv @ yv)
            print(f"   {LABEL[mult]:11s} omega={omega:10.6f}")
            print(f"      overlap^2 with the symmetric-path tangent = {overlap2:.9f}")
            print(f"      radial fraction of vertex motion           = {rad / tot:.9f}")
            print(f"      breathing  (sum_i r.u)^2 over block        = {breath:.3e}")
            print(f"      sum_i |r.u| over block (is it a real       = "
                  f"{breath_abs:.6f}")
            print(f"         cancellation, or is r.u zero per vertex?)")
            print(f"      body-rotation share of the chart vector    = "
                  f"{spin / Q.shape[1]:.6f}")
            print(f"      per-vertex participation min/max           = "
                  f"{part.min() / part.max():.9f}   (1 = every vertex equally)")
            print(f"      hinge tear (must be 0)                     = {hinge:.2e}")
    print()
    print("  READ: overlap^2 = 1 means the block IS the jitterbug's own motion;")
    print("  0 means the block is entirely transverse to the symmetric path.")
    print("  radial fraction 1 = pure breathing, 0 = pure tangential shear.")
    print("  (sum r.u) is the rate of change of sum |r_i|^2 -- the configuration's")
    print("  size. It must VANISH on any non-trivial irrep, because it is a")
    print("  linear functional invariant under the symmetry group; a non-zero")
    print("  value on the doublet or triplet would refute the block picture.")
    return F, pr, H


# --------------------------------------------------------------------------
# S4 -- vertices vs strut midpoints, RAW only (DECISION 17)
# --------------------------------------------------------------------------

def s4_primitive():
    print()
    print("=" * 78)
    print("S4  RAW VERTEX vs RAW STRUT-MIDPOINT: same ground state, same spectrum?")
    print("=" * 78)
    F = aligned_frame(0.0)
    pr = Probe(F)
    Ms = {m: metric_for(F, m) for m in MODELS}
    print("  First: is the VE critical for the strut kernels too? (re-measured)")
    for name, kern in RAW_KERNELS[:2]:
        g, _, _ = pr.grad(mid_V(kern, False), 1e-4, on_config=True)
        print(f"    {name:20s} |grad V_strut| = {np.linalg.norm(g):.3e}"
              f"   transverse {np.linalg.norm(g[1:]):.3e}")
    print()
    print(f"  {'kernel':20s} {'primitive':10s} {'omega_D':>11s} {'omega_T':>11s} "
          f"{'omega_S':>11s} {'T/D':>9s} {'S/D':>9s} {'S/T':>9s}")
    diffs = []
    for name, kern in RAW_KERNELS:
        rows = {}
        for prim, V, oc in (("vertex", make_V(kern, False), False),
                            ("strut", mid_V(kern, False), True)):
            H = pr.hess(V, 1e-3, on_config=oc)
            d = {}
            for mult, mval, Q in irrep_blocks(Ms["lamina"]):
                h, dev = block_stiffness(H, Q)
                d[mult] = np.sqrt(h / mval)
            rows[prim] = d
            print(f"  {name if prim == 'vertex' else '':20s} {prim:10s} "
                  f"{d[2]:11.6f} {d[3]:11.6f} {d[1]:11.6f} "
                  f"{d[3] / d[2]:9.6f} {d[1] / d[2]:9.6f} {d[1] / d[3]:9.6f}")
        rv, rs = rows["vertex"], rows["strut"]
        diffs.append((name,
                      abs(rs[3] / rs[2] - rv[3] / rv[2]),
                      abs(rs[1] / rs[2] - rv[1] / rv[2]),
                      (rs[2] / rv[2], rs[3] / rv[3], rs[1] / rv[1])))
    print()
    print("  A pure change of overall stiffness would leave the RATIO columns")
    print("  identical between the two primitives and give ONE common scale")
    print("  factor for all three blocks. Measured:")
    print(f"  {'kernel':20s} {'|dT/D|':>10s} {'|dS/D|':>10s}   "
          f"{'scale_D':>9s} {'scale_T':>9s} {'scale_S':>9s}  common?")
    for name, dtd, dsd, sc in diffs:
        common = max(sc) / min(sc)
        print(f"  {name:20s} {dtd:10.6f} {dsd:10.6f}   "
              f"{sc[0]:9.4f} {sc[1]:9.4f} {sc[2]:9.4f}  spread {common:.4f}x")


# --------------------------------------------------------------------------
# S5 -- off the equilibrium
# --------------------------------------------------------------------------

def s5_chart_control(kernel_name="1/r^1  (Thomson)"):
    print()
    print("=" * 78)
    print("S5a  CONTROL: two genuinely different charts must agree AT the VE")
    print("=" * 78)
    kern = dict(RAW_KERNELS)[kernel_name]
    V = make_V(kern, normalised=False)
    out = {}
    for tag, F, orig in (("centroid pivot", aligned_frame(0.0), False),
                         ("origin pivot", OriginFrame(0.0), True)):
        pr = Probe(F)
        g, _, _ = pr.grad(V, 1e-4)
        H = pr.hess(V, 1e-3)
        M = metric_for(F, "lamina", origin_pivot=orig)
        w, _ = gen_eig(H, M)
        out[tag] = w
        print(f"  {tag:16s} dim {F.dim}  |grad| {np.linalg.norm(g):.2e}   "
              f"omega^2 {fmt_blocks(degeneracy_blocks(w))}")
    d = np.abs(np.sort(out["centroid pivot"]) - np.sort(out["origin pivot"]))
    print(f"  max |domega^2| between the two charts = {d.max():.3e}"
          f"   (relative {d.max() / np.abs(out['origin pivot']).max():.2e})")
    print("  This is the control. If it failed, S5b would be measuring a bug.")
    return out


def s5_icosahedron(kernel_name="1/r^1  (Thomson)"):
    print()
    print("=" * 78)
    print("S5b  AT THE ICOSAHEDRON (a = 22.238756093): is it even an equilibrium?")
    print("=" * 78)
    kern = dict(RAW_KERNELS)[kernel_name]
    V = make_V(kern, normalised=False)
    res = {}
    for tag, F, orig in (("centroid pivot", aligned_frame(A_ICO), False),
                         ("origin pivot", OriginFrame(A_ICO), True)):
        pr = Probe(F)
        g, _, _ = pr.grad(V, 1e-4)
        H = pr.hess(V, 1e-3)
        print(f"\n  {tag}: dim {F.dim}   |grad V| = {np.linalg.norm(g):.6e}")
        if tag == "centroid pivot":
            print(f"    gradient in the aligned basis = "
                  f"{np.array2string(g, precision=6)}")
            print(f"    transverse part |g_1..5| = {np.linalg.norm(g[1:]):.3e}"
                  f"   -> the gradient is along the SYMMETRIC PATH, as the")
            print("       symmetry requires (the path is a fixed-point set)")
        (p, z, n), ev = inertia(H)
        print(f"    Hessian inertia ({p},{z},{n})  eig {fmt_blocks(degeneracy_blocks(ev))}")
        for m in MODELS:
            M = metric_for(F, m, origin_pivot=orig)
            w, _ = gen_eig(H, M)
            res[(tag, m)] = w
            om = np.sqrt(np.clip(w, 0, None))
            print(f"    {m:7s} omega^2 {fmt_blocks(degeneracy_blocks(w))}")
            print(f"    {m:7s} omega   {fmt_blocks(degeneracy_blocks(om))}"
                  f"   all positive: {bool(np.all(w > 0))}")
    print()
    print("  CHART DEPENDENCE, which is the point of running two charts:")
    for m in MODELS:
        a = np.sort(res[("centroid pivot", m)])
        b = np.sort(res[("origin pivot", m)])
        print(f"    {m:7s} max |domega^2| = {np.abs(a - b).max():.6f}"
              f"   relative to spectrum span = "
              f"{np.abs(a - b).max() / (a.max() - a.min()):.4f}")
    print("  Compare with the VE control above (agreement at ~1e-5 absolute).")


def s5_path_critical_points():
    print()
    print("=" * 78)
    print("S5c  THE WELL-POSED QUESTION: WHERE ELSE ON THE PATH IS V STATIONARY?")
    print("=" * 78)
    print("  'Is the VE the only stable point on the path' only has an answer at")
    print("  points that ARE equilibria. Scan dV/da over the whole circle")
    print("  (DECISION 16: no region excluded) and locate the sign changes.")
    print("  Note the a -> 180-a isometry: the fundamental domain is [0,90], so")
    print("  every stationary point found there has mirror copies.")
    for name in ("1/r^1  (Thomson)", "gauss s=1.0"):
        kern = dict(RAW_KERNELS)[name]
        V = make_V(kern, normalised=False)

        def f(a):
            return V(vert(a))

        h = 1e-5
        grid = np.arange(0.0, 360.0, 0.05)
        vals = np.array([f(a) for a in grid])
        finite = np.isfinite(vals)
        d = np.array([np.nan if not finite[i] else
                      (f(a + h) - f(a - h)) / (2 * h) for i, a in enumerate(grid)])
        roots = []
        for i in range(len(grid) - 1):
            if np.isfinite(d[i]) and np.isfinite(d[i + 1]) and d[i] == 0:
                roots.append(grid[i])
            elif np.isfinite(d[i]) and np.isfinite(d[i + 1]) and d[i] * d[i + 1] < 0:
                lo, hi = grid[i], grid[i + 1]
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    dm = (f(mid + h) - f(mid - h)) / (2 * h)
                    if not np.isfinite(dm):
                        break
                    if d[i] * dm < 0:
                        hi = mid
                    else:
                        lo = mid
                roots.append(0.5 * (lo + hi))
        print(f"\n  {name}: dV/da sign changes at "
              + ", ".join(f"{r:.6f}" for r in roots))
        print(f"    non-finite V on the grid at {int(np.sum(~finite))} of "
              f"{len(grid)} angles (poles where vertices merge)")
        for r in roots:
            if min(abs(r - x) for x in (90.0, 270.0)) < 1e-6:
                print(f"    a={r:11.6f}: BRANCH POINT (local dim 7). The chart's")
                print("               Newton solve fails there; no spectrum claimed,")
                print("               matching the prior survey's refusal.")
                continue
            F = aligned_frame(r)
            pr = Probe(F)
            g, _, _ = pr.grad(V, 1e-4)
            H = pr.hess(V, 1e-3)
            (p, z, n), ev = inertia(H)
            line = (f"    a={r:11.6f}: dim {F.dim}  |grad_6D| "
                    f"{np.linalg.norm(g):.2e}  inertia ({p},{z},{n})")
            for m in MODELS:
                w, _ = gen_eig(H, metric_for(F, m))
                line += f"   {m}: omega {np.sqrt(np.abs(w)).min():.4f}.." \
                        f"{np.sqrt(np.abs(w)).max():.4f}"
            print(line)
            print(f"               all omega^2 > 0: "
                  f"{bool(np.all(np.linalg.eigvalsh(H) > 0))}")


def path_critical(V, lo, hi, step=0.05, h=1e-5, verbose=False):
    """Verified sign changes of dV/da on [lo, hi], refined by bisection.

    TRAP, caught in this session and kept as a guard rather than quietly fixed:
    a POLE manufactures a false stationary point. An inverse power is +inf at
    a=60 (the 12 vertices merge in pairs), so dV/da runs +inf just below and
    -inf just above -- a textbook sign change that bisection happily refines to
    a=60, where V is not even finite and the Hessian is NaN. Every root is
    therefore VERIFIED: V finite, and |dV/da| actually small at the refined
    root. Rejected brackets are reported, not silently dropped.
    """
    def f(a):
        return V(vert(a))

    grid = np.arange(lo, hi, step)
    d = []
    for a in grid:
        vp, vm = f(a + h), f(a - h)
        d.append((vp - vm) / (2 * h) if np.isfinite(vp) and np.isfinite(vm)
                 else np.nan)
    d = np.array(d)
    roots, rejected = [], []
    for i in range(len(grid) - 1):
        if not (np.isfinite(d[i]) and np.isfinite(d[i + 1])):
            continue
        if d[i] == 0.0:
            cand = float(grid[i])
        elif d[i] * d[i + 1] < 0:
            a0, a1 = grid[i], grid[i + 1]
            for _ in range(80):
                mid = 0.5 * (a0 + a1)
                dm = (f(mid + h) - f(mid - h)) / (2 * h)
                if not np.isfinite(dm):
                    break
                if d[i] * dm < 0:
                    a1 = mid
                else:
                    a0 = mid
            cand = 0.5 * (a0 + a1)
        else:
            continue
        vc = f(cand)
        dc = (f(cand + h) - f(cand - h)) / (2 * h)
        if np.isfinite(vc) and np.isfinite(dc) and abs(dc) < 1e-3 * (1 + abs(vc)):
            roots.append(cand)
        else:
            rejected.append((cand, vc, dc))
    if verbose and rejected:
        for c, vc, dc in rejected:
            print(f"      REJECTED false root at a={c:.6f}: V={vc}, "
                  f"dV/da={dc}  (pole, not a stationary point)")
    return roots


def s5d_second_minimum():
    print()
    print("=" * 78)
    print("S5d  IS THE VE THE ONLY STABLE POINT ON THE PATH? -- all nine kernels")
    print("=" * 78)
    print("  Scanned on the fundamental domain [0, 90] (the a -> 180-a isometry")
    print("  mirrors everything found here; DECISION 16 excludes no region, so")
    print("  the interference band (60,120) IS scanned).")
    print("  MECHANISM TO WATCH FOR: an inverse power has a POLE at a=60 and")
    print("  a=120 where the 12 vertices merge in pairs. A function that is +inf")
    print("  at both ends of (60,120) must have an interior minimum. A Gaussian")
    print("  has no pole and so need not. If the two families split exactly that")
    print("  way, the extra equilibrium is a property of the pole, not of shape.")
    print()
    print(f"  {'kernel':20s} {'stationary a in (0,90)':>24s} {'V there':>14s} "
          f"{'V(VE)':>14s} {'inertia':>9s} {'|grad_6D|':>10s}")
    found = {}
    for name, kern in RAW_KERNELS:
        V = make_V(kern, normalised=False)
        roots = [r for r in path_critical(V, 0.0, 90.0, verbose=True)
                 if 1e-3 < r < 90.0 - 1e-3]
        v_ve = V(vert(0.0))
        if not roots:
            print(f"  {name:20s} {'none':>24s} {'--':>14s} {v_ve:14.8f} "
                  f"{'--':>9s} {'--':>10s}")
            found[name] = None
            continue
        for r in roots:
            F = aligned_frame(r)
            pr = Probe(F)
            g, _, _ = pr.grad(V, 1e-4)
            H = pr.hess(V, 1e-3)
            (p, z, n), ev = inertia(H)
            print(f"  {name:20s} {r:24.6f} {V(vert(r)):14.8f} {v_ve:14.8f} "
                  f"({p},{z},{n})   {np.linalg.norm(g):10.2e}")
            found[name] = (r, V(vert(r)), (p, z, n))
    print()
    print("  Spectra at the second minimum, both mass models "
          "(1/r^1 and 1/r^12):")
    for name in ("1/r^1  (Thomson)", "1/r^12"):
        if found.get(name) is None:
            continue
        r = found[name][0]
        F = aligned_frame(r)
        H = Probe(F).hess(make_V(dict(RAW_KERNELS)[name], False), 1e-3)
        print(f"    {name} at a = {r:.6f}")
        for m in MODELS:
            M = metric_for(F, m)
            for mult, mval, Q in irrep_blocks(M):
                h, dev = block_stiffness(H, Q)
                print(f"      {m:7s} {LABEL[mult]:11s} omega = "
                      f"{np.sqrt(h / mval):11.6f}   |H-scalar| = {dev:.2e}")
    print()
    print("  ROBUSTNESS of the extra root (do not fit a number to one run):")
    V = make_V(dict(RAW_KERNELS)["1/r^1  (Thomson)"], False)
    for step in (0.2, 0.05, 0.01):
        for h in (1e-4, 1e-5, 1e-6):
            rr = [r for r in path_critical(V, 0.0, 90.0, step, h)
                  if 1e-3 < r < 90.0 - 1e-3]
            print(f"    grid {step:5.2f} deg, h {h:.0e} -> "
                  + (", ".join(f"{r:.7f}" for r in rr) if rr else "none"))
    print("  and the pole structure that would explain it:")
    for a in (59.9, 59.99, 60.0, 60.01, 60.1, 75.262042, 89.9):
        print(f"    V_thomson({a:10.6f}) = {V(vert(a)):.6f}")


def path_critical_f(f, lo, hi, step=0.05, h=1e-5):
    """`path_critical` for a bare f(a). Same pole-verification guard."""
    grid = np.arange(lo, hi, step)
    d = []
    for a in grid:
        vp, vm = f(a + h), f(a - h)
        d.append((vp - vm) / (2 * h) if np.isfinite(vp) and np.isfinite(vm)
                 else np.nan)
    d = np.array(d)
    roots = []
    for i in range(len(grid) - 1):
        if not (np.isfinite(d[i]) and np.isfinite(d[i + 1])):
            continue
        if d[i] * d[i + 1] >= 0:
            continue
        a0, a1 = grid[i], grid[i + 1]
        for _ in range(80):
            mid = 0.5 * (a0 + a1)
            dm = (f(mid + h) - f(mid - h)) / (2 * h)
            if not np.isfinite(dm):
                break
            if d[i] * dm < 0:
                a1 = mid
            else:
                a0 = mid
        cand = 0.5 * (a0 + a1)
        vc, dc = f(cand), (f(cand + h) - f(cand - h)) / (2 * h)
        if np.isfinite(vc) and np.isfinite(dc) and abs(dc) < 1e-3 * (1 + abs(vc)):
            roots.append((cand, float(vc), "min" if d[i] < 0 else "max"))
    return roots


def s5e_energy_and_barrier():
    """Is the second equilibrium ABOVE or BELOW the VE, and is it separated?

    The stability verdict (inertia (6,0,0)) says only that the second point is a
    LOCAL minimum. Whether DECISION 17 survives turns on the ENERGY: above the
    VE means metastable and the decision stands; below means the VE was never
    the ground state and the earlier surveys missed it by not scanning inside
    the interference band.

    HYPOTHESIS BEING TESTED, recorded so it can be broken rather than
    confirmed: a monotone-repulsive kernel should prefer the configuration with
    the largest separations, and the VE has the largest circumradius on the
    path, so the VE should win. NOTE the heuristic is NOT a proof and has
    already failed once in this project -- the centroid-NORMALISED kernels
    prefer the icosahedron, so 'repulsion picks the biggest thing' does not
    survive contact with shape. The energies below are what decides it.
    """
    print()
    print("=" * 78)
    print("S5e  ENERGY OF THE SECOND EQUILIBRIUM, AND THE BARRIER")
    print("=" * 78)
    prims = (("vertex", lambda kern: (lambda a: make_V(kern, False)(vert(a)))),
             ("strut", lambda kern: (lambda a: mid_V(kern, False)(corners(a)))))
    for pname, build in prims:
        print(f"\n  --- primitive: {pname} ---")
        print(f"  {'kernel':20s} {'a2':>11s} {'V(a2)':>15s} {'V(VE)':>15s} "
              f"{'V(a2)-V(VE)':>15s}  verdict")
        for name, kern in RAW_KERNELS:
            f = build(kern)
            v0 = f(0.0)
            roots = [r for r in path_critical_f(f, 0.0, 90.0)
                     if r[2] == "min" and r[0] > 1e-3]
            if not roots:
                print(f"  {name:20s} {'none':>11s} {'--':>15s} {v0:15.8f} "
                      f"{'--':>15s}  no second minimum")
                continue
            for a2, v2, _ in roots:
                verdict = ("ABOVE -> metastable, DECISION 17 stands" if v2 > v0
                           else "BELOW -> THE VE IS NOT THE GROUND STATE")
                print(f"  {name:20s} {a2:11.6f} {v2:15.8f} {v0:15.8f} "
                      f"{v2 - v0:+15.8f}  {verdict}")

    print()
    print("  --- the barrier between the two minima, along the path ---")
    print("  DERIVED, no run needed, for every inverse power and BOTH primitives:")
    print("  at a=60 the 12 shared vertices merge in pairs (and the 24 struts")
    print("  coincide in pairs, so 12 midpoint pairs also merge), so some r_ij")
    print("  is exactly 0 and 1/r^p is +inf. The barrier is INFINITE. Numeric")
    print("  corroboration, not proof of the limit:")
    fv = lambda a: make_V(dict(RAW_KERNELS)["1/r^1  (Thomson)"], False)(vert(a))
    for a in (59.0, 59.9, 59.99, 59.999, 60.0):
        print(f"    V_thomson(vertex, {a:9.5f}) = {fv(a):.6f}")
    print()
    print("  MEASURED for the members with no pole (Gaussians):")
    for name, kern in RAW_KERNELS:
        if not name.startswith("gauss"):
            continue
        for pname, build in prims:
            f = build(kern)
            v0 = f(0.0)
            roots = path_critical_f(f, 0.0, 90.0)
            mins = [r for r in roots if r[2] == "min" and r[0] > 1e-3]
            maxs = [r for r in roots if r[2] == "max"]
            if not mins:
                continue
            a2, v2, _ = mins[0]
            if not maxs:
                print(f"    {name} / {pname}: second min at {a2:.6f} but NO "
                      f"interior maximum found -- barrier not located")
                continue
            ab, vb, _ = max(maxs, key=lambda r: r[1])
            print(f"    {name} / {pname}: min {a2:.6f} V={v2:.8f} | "
                  f"barrier top {ab:.6f} V={vb:.8f} | VE V={v0:.8f}")
            print(f"        barrier above the second well = {vb - v2:.8f} "
                  f"({(vb - v2) / (v2 - v0) * 100:.2f}% of the well's excitation "
                  f"above the VE)")
    print()
    print("  SCOPE LIMIT, and it is load-bearing: every barrier here is measured")
    print("  ALONG THE SYMMETRIC PATH ONLY. The variety is 6-dimensional, so a")
    print("  lower route around the a=60 wall through the five transverse")
    print("  directions is NOT excluded by anything measured. 'Infinite barrier'")
    print("  means infinite ON THE PATH, not dynamically disconnected.")


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=False, linewidth=170)
    s3_modes()
    s4_primitive()
    s5_chart_control()
    s5_icosahedron()
    s5_path_critical_points()
    s5d_second_minimum()
