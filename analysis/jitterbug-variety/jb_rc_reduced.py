"""jb_rc -- the medium in reduced coordinates, where the ellipses cannot be left.

WHAT THIS REPLACES. Every integration in this programme so far ran in the
BAR-ONLY constraint space: strut lengths held, nothing else. That space is
strictly weaker than the jitterbug. It admits SIX internal degrees of freedom
per cell where the linkage has ONE, and the five extra are not physics -- the
jumbled cluster render was those modes made visible. The fix is not a stiffer
constraint set. It is a different coordinate system.

THE COORDINATES. Per cell k: a centre c_k (3), an orientation R_k (3), and a
fold angle g_k (1). Seven generalized coordinates. A cell's vertex instance i
sits at

    x_{k,i} = c_k + R_k . v_i(g_k)

with v_i the body-frame jitterbug vertex -- jb_ic's `cell_verts(g, 0)`. Two
things then hold BY CONSTRUCTION rather than by constraint, and they are the
whole reason for the change:

  * All 24 strut lengths. The struts are the edges of the eight rigid triangles,
    and v(g) turns each triangle rigidly about its own axis. R2 gates that the
    fold tangent has d(strut^2)/dg identically zero, so bars are not constraints.
  * The ellipses. v_i is a function of g ALONE. There is no coordinate in the
    state that could move a vertex off its ellipse, so the ellipse residual is
    not small, it is absent: R1 gates that the constraint set contains ZERO
    ellipse rows. The invariants themselves -- perpendicular distance EL/sqrt(3)
    to the axis of each of the two faces the vertex belongs to, and one Cartesian
    coordinate identically zero -- are properties of the fixed body function,
    checked once over a sweep in g, not re-checked every step of every run.

So the ONLY constraints left are the welds. Nine scalar rows per shared face, of
which exactly SIX are independent: both mating triangles are equilateral of edge
EL at every g, so their shapes always agree and three rows are dependent. R5
measures that rank rather than assuming it.

WHAT THE COUNT THEN SAYS. For a tree of welds, 7N - 6(N-1) = N + 6, and taking
out the six global rigid motions leaves N. One fold angle per cell survives, at
every size. That is the field claim, and it is now a measurement in the
coordinates the model is written in. Freeze the orientations and the same
Jacobian -- restricted to its centre and fold columns -- has rank 4 per weld
instead of 6, leaving 4N - 4(N-1) - 3 = 1 at every size: a single coordinate, and
nowhere for a disturbance to be.

THE CARRIED-FORWARD DEFECT IS RESOLVED HERE. An earlier fold-map-rank run
projected out a six-dimensional "global motion" space built with an identity
orientation block and no omega x c term on the centres, and reported a DOF column
two too high. The correct basis is explicit below -- translations are
(dc_k = t, w_k = 0), rotations are (dc_k = w x c_k, w_k = w) -- and R5a gates
that C annihilates it, printing the residual of the defective basis alongside so
the row cannot pass vacuously.

V = 0, SO ENERGY IS THE ONLY AUDIT. The motion is a geodesic of the kinetic
energy metric on the weld manifold. This integrates the acceleration-level
equations with DOP853 at tight tolerance, deliberately NOT with a symplectic
scheme: a symplectic integrator conserves a shadow Hamiltonian and would make
the energy row nearly vacuous, whereas a high-order non-symplectic one has no
such protection and will drift visibly if the equations of motion are wrong.

MASS MODEL, DECLARED: unit mass per triangle, lumped m/3 to each corner, which
is jb_ic's model so the two are comparable. Each cell carries 24 corner
instances at 1/3 -- 2/3 on each of its twelve vertex IDENTITIES, eight units per
cell at every configuration. R3 gates that this survives the octahedron, where
the twelve identities occupy six positions and a position-dedup build would
carry four.
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.integrate import solve_ivp

import jb_gp_plate_geometry as Z
import jb_ic_inertial_chain as IC
import jb_cl_cluster as CL
from jb_a_family import Z as ZAX, rot

EL = IC.EL
NV = IC.NV                 # 12 vertex identities per cell, at every angle
SLOT = IC.SLOT
FACES = Z.faces()
NCORNER = 24               # (face, corner) slots per cell

#: the two faces each vertex identity belongs to -- the pair of axes whose
#: cylinders intersect in that vertex's ellipse
VFACES = {i: [] for i in range(NV)}
for (_f, _c), _i in SLOT.items():
    VFACES[_i].append(_f)

#: unit mass per triangle, lumped m/3 to each corner. Counted from the slot
#: incidence rather than hardcoded, so congruence cannot be lost by a constant.
VMASS = np.zeros(NV)
for (_f, _c), _i in SLOT.items():
    VMASS[_i] += 1.0 / 3.0

TRIS = [tuple(SLOT[(f, c)] for c in range(3)) for f in range(8)]
BARS = [(t[a], t[b]) for t in TRIS for a, b in ((0, 1), (1, 2), (0, 2))]


# --------------------------------------------------------------------------
# the body function and its two derivatives, analytic in the fold angle
# --------------------------------------------------------------------------

def body(a_deg, nder=0):
    """Body-frame vertices v(g) and, on request, dv/dg and d2v/dg2 in RADIANS.

    x(a) = R(u, sigma*(a-60)) @ (v - c) + u * Z * cos(a), so with K the
    cross-product matrix of the face axis u and R(theta) = exp(theta K),

        dx/da  = sigma * K R (v-c)     - u Z sin(a)
        d2x/da2 = sigma^2 K^2 R (v-c)  - u Z cos(a)

    R3b checks both against central differences, because a sign on sigma here
    would be invisible in the positions and wrong in every velocity.
    """
    a = np.radians(a_deg)
    out = [np.zeros((NV, 3)) for _ in range(nder + 1)]
    seen = np.zeros(NV, bool)
    for f, (v, c, u, sg) in enumerate(FACES):
        Rm = rot(u, sg * (a_deg - 60.0))
        K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
        base = (Rm @ (v - c).T).T
        for cc in range(3):
            i = SLOT[(f, cc)]
            if seen[i]:
                continue
            seen[i] = True
            out[0][i] = base[cc] + u * ZAX * np.cos(a)
            if nder >= 1:
                out[1][i] = sg * (K @ base[cc]) - u * ZAX * np.sin(a)
            if nder >= 2:
                out[2][i] = sg * sg * (K @ K @ base[cc]) - u * ZAX * np.cos(a)
    return out[0] if nder == 0 else tuple(out)


def hat(w):
    return np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])


def _solve(C, MiCt, rhs):
    """Least-squares multiplier solve. Do NOT rewrite the np.dot as `@`.

    On numpy 1.26 against Apple's Accelerate BLAS the `@` operator's small-matrix
    SIMD path raises spurious divide-by-zero / overflow / invalid RuntimeWarnings
    on these operands, while np.dot returns a bit-identical, finite result. The
    warnings are a status-flag artifact, not a numerical problem -- verified by
    np.array_equal on the two products -- but a gate that prints floating point
    warnings on every run trains its reader to ignore them.
    """
    return np.linalg.lstsq(np.dot(C, MiCt), rhs, rcond=None)[0]


def quat_to_R(q):
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
                     [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
                     [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def quat_dot(q, w):
    """Spatial angular velocity: qdot = 0.5 * (0, w) x q."""
    qw, qv = q[0], q[1:]
    return 0.5 * np.concatenate(([-np.dot(w, qv)], qw * w + np.cross(w, qv)))


# --------------------------------------------------------------------------
# a welded assembly of cells, in reduced coordinates
# --------------------------------------------------------------------------

class Assembly:
    """Cells and the welds between them. The welds are the whole constraint set.

    `welds` is a list of (k, l, [(a, b) x 3]) -- cell k's vertex a coincides with
    cell l's vertex b. The corner correspondence comes from jb_ic and jb_cl, whose
    R2b row gates that it lands on the mating FACE in both offset directions.
    """

    def __init__(self, gam, ctr, welds):
        self.N = len(gam)
        self.gam0 = np.asarray(gam, float)
        self.ctr0 = np.asarray(ctr, float)
        self.welds = welds
        self.nc = 9 * len(welds)

    # -- state packing: centre(3) + quaternion(4) + gamma(1) per cell ------
    def q0(self):
        q = np.zeros((self.N, 8))
        q[:, 0:3] = self.ctr0
        q[:, 3] = 1.0
        q[:, 7] = self.gam0
        return q.ravel()

    @staticmethod
    def unpack(q):
        Q = q.reshape(-1, 8)
        return Q[:, 0:3], Q[:, 3:7], Q[:, 7]

    def frames(self, q):
        ctr, quat, gam = self.unpack(q)
        R = np.array([quat_to_R(qq) for qq in quat])
        B = [body(g, 2) for g in gam]
        return ctr, R, gam, B

    # -- kinematics -------------------------------------------------------
    def positions(self, q):
        ctr, R, gam, B = self.frames(q)
        return np.array([ctr[k] + np.dot(R[k], B[k][0].T).T for k in range(self.N)])

    def cell_jacobians(self, ctr, R, B):
        """J_k : (36, 7), the map from (cdot, omega, gdot) to vertex velocities."""
        J = np.zeros((self.N, 3 * NV, 7))
        for k in range(self.N):
            RB0 = np.dot(R[k], B[k][0].T).T
            RB1 = np.dot(R[k], B[k][1].T).T
            for i in range(NV):
                r = 3 * i
                J[k, r:r + 3, 0:3] = np.eye(3)
                J[k, r:r + 3, 3:6] = -hat(RB0[i])
                J[k, r:r + 3, 6] = RB1[i]
        return J

    def cell_gyro(self, R, B, u):
        """The velocity-only part of the ambient acceleration, Jdot . u.

        d2/dt2 (c + R v) = cddot + wdot x Rv + R v' gddot
                           + [ w x (w x Rv) + 2 w x (R v' gdot) + R v'' gdot^2 ]
        and the bracket is what this returns, per vertex.
        """
        A = np.zeros((self.N, NV, 3))
        for k in range(self.N):
            w, gd = u[k, 3:6], u[k, 6]
            RB0 = np.dot(R[k], B[k][0].T).T
            RB1 = np.dot(R[k], B[k][1].T).T
            RB2 = np.dot(R[k], B[k][2].T).T
            A[k] = (np.cross(w, np.cross(w, RB0))
                    + 2.0 * np.cross(w, RB1 * gd)
                    + RB2 * gd * gd)
        return A

    def mass_blocks(self, J):
        """M is block diagonal over cells: a cell's vertices depend on its own
        seven coordinates and nothing else."""
        m3 = np.repeat(VMASS, 3)
        return np.array([np.dot(J[k].T, m3[:, None] * J[k]) for k in range(self.N)])

    # -- constraints ------------------------------------------------------
    def weld_residual(self, q):
        X = self.positions(q)
        g = np.zeros(self.nc)
        for r, (k, l, pairs) in enumerate(self.welds):
            for m, (a, b) in enumerate(pairs):
                g[9 * r + 3 * m:9 * r + 3 * m + 3] = X[k][a] - X[l][b]
        return g

    def constraint_jacobian(self, J):
        C = np.zeros((self.nc, 7 * self.N))
        for r, (k, l, pairs) in enumerate(self.welds):
            for m, (a, b) in enumerate(pairs):
                row = 9 * r + 3 * m
                C[row:row + 3, 7 * k:7 * k + 7] = J[k][3 * a:3 * a + 3]
                C[row:row + 3, 7 * l:7 * l + 7] = -J[l][3 * b:3 * b + 3]
        return C

    def constraint_gyro(self, A):
        d = np.zeros(self.nc)
        for r, (k, l, pairs) in enumerate(self.welds):
            for m, (a, b) in enumerate(pairs):
                d[9 * r + 3 * m:9 * r + 3 * m + 3] = A[k][a] - A[l][b]
        return d

    def globals(self, ctr, defective=False):
        """The six global rigid motions, in the reduced coordinates.

        translations  dc_k = t,        w_k = 0,  gdot_k = 0
        rotations     dc_k = w x c_k,  w_k = w,  gdot_k = 0

        `defective=True` reproduces the basis that produced the off-by-2: an
        identity orientation block and no omega x c term on the centres. It is
        here so R5a can print what a wrong basis does, rather than assert it.
        """
        G = np.zeros((6, 7 * self.N))
        for d in range(3):
            for k in range(self.N):
                G[d, 7 * k + d] = 1.0
        for d in range(3):
            w = np.eye(3)[d]
            for k in range(self.N):
                if not defective:
                    G[3 + d, 7 * k:7 * k + 3] = np.cross(w, ctr[k])
                G[3 + d, 7 * k + 3 + d] = 1.0
        return G

    def momentum_rows(self, ctr, R, B, J):
        """Total linear and angular momentum as linear functionals of u."""
        m3 = np.repeat(VMASS, 3)
        P = np.zeros((6, 7 * self.N))
        for k in range(self.N):
            X = ctr[k] + np.dot(R[k], B[k][0].T).T
            S = np.zeros((6, 3 * NV))
            for i in range(NV):
                S[0:3, 3 * i:3 * i + 3] = VMASS[i] * np.eye(3)
                S[3:6, 3 * i:3 * i + 3] = VMASS[i] * hat(X[i])
            P[:, 7 * k:7 * k + 7] = np.dot(S, J[k])
        return P

    # -- dynamics ---------------------------------------------------------
    def accel(self, q, u):
        """u-dot from the acceleration-level equations. V = 0, so the only
        forces are the weld reactions. C is rank deficient by three rows per
        weld, so the multiplier solve is least squares."""
        ctr, R, gam, B = self.frames(q)
        J = self.cell_jacobians(ctr, R, B)
        A = self.cell_gyro(R, B, u)
        M = self.mass_blocks(J)
        m3 = np.repeat(VMASS, 3)
        f = np.array([-np.dot(J[k].T, m3 * A[k].ravel()) for k in range(self.N)])
        Minv = np.array([np.linalg.inv(M[k]) for k in range(self.N)])
        a_free = np.einsum('kij,kj->ki', Minv, f)
        if not self.welds:
            return a_free
        C = self.constraint_jacobian(J)
        d = -self.constraint_gyro(A)
        MiCt = np.zeros((7 * self.N, self.nc))
        for k in range(self.N):
            MiCt[7 * k:7 * k + 7] = np.dot(Minv[k], C[:, 7 * k:7 * k + 7].T)
        lam = _solve(C, MiCt, d - np.dot(C, a_free.ravel()))
        return a_free + np.dot(MiCt, lam).reshape(self.N, 7)

    def energy(self, q, u):
        ctr, R, gam, B = self.frames(q)
        M = self.mass_blocks(self.cell_jacobians(ctr, R, B))
        return np.array([0.5 * float(np.dot(u[k], np.dot(M[k], u[k])))
                         for k in range(self.N)])

    def project_velocity(self, q, u0, momentum=True):
        """Nearest admissible velocity in the kinetic-energy metric.

        Minimises (u-u0)' M (u-u0) subject to C u = 0, and -- when `momentum` --
        to zero total linear and angular momentum. The momentum rows matter: a
        fold impulse on one cell does induce net momentum, and without them part
        of what follows would be the whole patch drifting rather than transport.
        """
        ctr, R, gam, B = self.frames(q)
        J = self.cell_jacobians(ctr, R, B)
        M = self.mass_blocks(J)
        Minv = np.array([np.linalg.inv(M[k]) for k in range(self.N)])
        rows = [self.constraint_jacobian(J)] if self.welds else []
        if momentum:
            rows.append(self.momentum_rows(ctr, R, B, J))
        C = np.vstack(rows)
        MiCt = np.zeros((7 * self.N, C.shape[0]))
        for k in range(self.N):
            MiCt[7 * k:7 * k + 7] = np.dot(Minv[k], C[:, 7 * k:7 * k + 7].T)
        lam = _solve(C, MiCt, -np.dot(C, u0.ravel()))
        return u0 + np.dot(MiCt, lam).reshape(self.N, 7), C

    def run(self, u0, tmax, nsample, rtol=1e-11, atol=1e-12):
        q = self.q0()

        def rhs(_t, y):
            qq, uu = y[:8 * self.N], y[8 * self.N:].reshape(self.N, 7)
            ctr, quat, gam = self.unpack(qq)
            dq = np.zeros((self.N, 8))
            dq[:, 0:3] = uu[:, 0:3]
            for k in range(self.N):
                dq[k, 3:7] = quat_dot(quat[k] / np.linalg.norm(quat[k]), uu[k, 3:6])
            dq[:, 7] = np.degrees(uu[:, 6])
            return np.concatenate([dq.ravel(), self.accel(qq, uu).ravel()])

        ts = np.linspace(0.0, tmax, nsample + 1)
        sol = solve_ivp(rhs, (0.0, tmax), np.concatenate([q, u0.ravel()]),
                        method="DOP853", t_eval=ts, rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        rec = []
        for j, t in enumerate(sol.t):
            qq = sol.y[:8 * self.N, j]
            uu = sol.y[8 * self.N:, j].reshape(self.N, 7)
            rec.append((t, self.energy(qq, uu), self.unpack(qq)[2].copy(),
                        float(np.abs(self.weld_residual(qq)).max())))
        return rec


#: gamma is carried in DEGREES in the state and gdot in RADIANS/s, because the
#: body function is written in degrees and its derivatives in radians. `run`
#: converts once, here, and nowhere else.


def chain(N, g0=30.0, flip=True, drop=()):
    """A line of N cells. Alternating +-g0 satisfies the shared-face law
    b = a + 60 at 30/-30, and the corner correspondence flips with the direction
    of the offset -- jb_ic's second construction trap, handled by asking
    `weld_for` per link rather than reusing one pairing. `drop` omits links,
    for the mutation probe."""
    gam = [g0 if (k % 2 == 0 or not flip) else -g0 for k in range(N)]
    seps = [IC.ZC * (np.cos(np.radians(gam[k])) + np.cos(np.radians(gam[k + 1])))
            for k in range(N - 1)]
    ctr = [np.zeros(3)]
    for s in seps:
        ctr.append(ctr[-1] + s * IC.NH)
    welds = [(k, k + 1, IC.weld_for(gam[k], gam[k + 1]))
             for k in range(N - 1) if k not in drop]
    return Assembly(gam, ctr, welds)


def cluster(gc=0.0, sites=None, drop=()):
    """The centre cell plus its <111> neighbours -- the smallest patch with an
    interior. `drop` omits welds, for the mutation probe."""
    sites = CL.SITES if sites is None else sites
    gn = gc + 60.0
    sep = IC.ZC * (np.cos(np.radians(gc)) + np.cos(np.radians(gn)))
    L = sep / np.sqrt(3)
    ctr = [np.zeros(3)] + [L * np.array(s, float) for s in sites]
    gam = [gc] + [gn] * len(sites)
    A = IC.cell_verts(gc, np.zeros(3))
    welds = []
    for m, s in enumerate(sites):
        if m in drop:
            continue
        d = np.array(s, float)
        fc, fn = CL.face_along(d), CL.face_along(-d)
        Bv = IC.cell_verts(gn, ctr[m + 1])
        bc = [SLOT[(fn, c)] for c in range(3)]
        used, pairs = set(), []
        for c in range(3):
            a = SLOT[(fc, c)]
            k = min([b for b in bc if b not in used],
                    key=lambda b: float(np.linalg.norm(Bv[b] - A[a])))
            used.add(k)
            pairs.append((a, k))
        welds.append((0, m + 1, pairs))
    return Assembly(gam, ctr, welds)


def rank_of(M, tol=1e-8):
    if M.size == 0:
        return 0, 0.0
    s = np.linalg.svd(M, compute_uv=False)
    r = int((s > s[0] * tol).sum())
    gap = (s[r - 1] / s[r]) if 0 < r < len(s) else np.inf
    return r, gap


def fold_impulse(asm, cell, momentum=True):
    """The local disturbance the field framing asks for: a phase-rate kick on
    ONE cell, gdot = 1, everything else zero, then projected onto the admissible
    set. A cell has no other internal freedom -- the triangle spin used before
    this coordinate change was never on the jitterbug path, which is why its
    projection had to mangle it."""
    u0 = np.zeros((asm.N, 7))
    u0[cell, 6] = 1.0
    return asm.project_velocity(asm.q0(), u0, momentum=momentum)


# --------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------

def gate():
    checks, out = [], {}
    A = checks.append

    # R1 -- the ellipse is not in the constraint set, and never could be
    ch6 = chain(6)
    cl9 = cluster()
    only_welds = (ch6.nc == 9 * 5 and cl9.nc == 9 * 8)
    worst_r, worst_plane = 0.0, 0.0
    for g in np.linspace(-95.0, 95.0, 61):
        V = body(g)
        for i in range(NV):
            worst_plane = max(worst_plane, float(np.abs(V[i]).min()))
            for f in VFACES[i]:
                u = FACES[f][2]
                worst_r = max(worst_r, abs(float(np.linalg.norm(V[i] - np.dot(V[i], u) * u))
                                           - EL / np.sqrt(3.0)))
    A(("R1  THE ELLIPSE RESIDUAL IS ABSENT, NOT SMALL. The constraint set is "
       "welds and nothing else -- nine rows per shared face, zero ellipse rows "
       "-- because the body-frame vertex is a function of the fold angle ALONE "
       "and no coordinate in the state could move it off its ellipse. The "
       "invariants are checked ONCE here, on the fixed body function: each "
       "vertex stays at perpendicular distance EL/sqrt(3) from the axis of "
       "BOTH faces it belongs to, and one Cartesian coordinate is identically "
       "zero. A run does not re-verify them because a run cannot violate them",
       only_welds and worst_r < 1e-14 and worst_plane < 1e-14,
       f"6-chain {ch6.nc} constraint rows = 9 x {len(ch6.welds)} welds, "
       f"9-cluster {cl9.nc} = 9 x {len(cl9.welds)}, ellipse rows 0; over 61 "
       f"angles in [-95, 95]: max |axis distance - {EL / np.sqrt(3.0):.6f}| = "
       f"{worst_r:.1e}, max |in-plane coordinate| = {worst_plane:.1e}",
       "welds only, and the body function's invariants at machine precision"))

    # R2 -- struts are structural too
    Lmin, Lmax, tang = np.inf, 0.0, 0.0
    for g in np.linspace(-95.0, 95.0, 41):
        V0, V1, _ = body(g, 2)
        for i, j in BARS:
            d, dv = V0[i] - V0[j], V1[i] - V1[j]
            Lmin, Lmax = min(Lmin, np.linalg.norm(d)), max(Lmax, np.linalg.norm(d))
            tang = max(tang, abs(2.0 * float(np.dot(d, dv))))
    A(("R2  THE BARS ARE NOT CONSTRAINTS EITHER. All 24 struts hold length at "
       "every fold angle, and the fold tangent has d(strut^2)/dg identically "
       "zero, so the entire bar-only constraint set of every previous "
       "integration is satisfied by construction here. What that set did NOT "
       "do is forbid the other five motions per cell, which is R6",
       abs(Lmax - EL) < 1e-13 and abs(Lmin - EL) < 1e-13 and tang < 1e-13,
       f"strut length over the sweep {Lmin:.12f}..{Lmax:.12f} against "
       f"EL = {EL:.12f}; max |d(strut^2)/dg| along the fold {tang:.1e}",
       "length EL everywhere, tangent derivative at machine zero"))

    # R3 -- congruence survives into the mass, and the derivatives are right
    m_tot = float(VMASS.sum())
    npos60 = len({tuple(np.round(p, 9)) for p in body(60.0)})
    A(("R3  THE CONGRUENT QUANTA ARE STILL CARRIED, in the coordinates as well "
       "as in the geometry. A cell has 24 corner instances at mass 1/3, so "
       "2/3 on each of its TWELVE vertex identities and 8 units per cell at "
       "EVERY configuration. At the octahedron those twelve identities occupy "
       "six positions; a build that deduplicated by position would carry four, "
       "which is the half-inertia Fuller's accounting objects to losing",
       abs(m_tot - 8.0) < 1e-12 and len(VMASS) == 12 and npos60 == 6
       and all(len(v) == 2 for v in VFACES.values()),
       f"{NCORNER} corner instances -> {NV} identities at mass "
       f"{VMASS[0]:.6f}, total {m_tot:.6f}; at g=60 they occupy {npos60} "
       f"distinct positions (a position build would carry "
       f"{npos60 * VMASS[0]:.1f})",
       "12 identities, 8 units of mass, 6 positions at the octahedron"))

    g0, h = 23.7, 1e-5
    V0, V1, V2 = body(g0, 2)
    fd1 = (body(g0 + np.degrees(h)) - body(g0 - np.degrees(h))) / (2 * h)
    fd2 = (body(g0 + np.degrees(h)) - 2 * V0 + body(g0 - np.degrees(h))) / h ** 2
    e1, e2 = float(np.abs(V1 - fd1).max()), float(np.abs(V2 - fd2).max())
    A(("R3b THE ANALYTIC DERIVATIVES ARE THE RIGHT ONES. dv/dg and d2v/dg2 "
       "carry the per-face sign sigma, which is invisible in the POSITIONS and "
       "wrong in every velocity if it is dropped. Checked against central "
       "differences; the second difference is roundoff limited at ~eps/h^2, so "
       "1e-6 there is agreement, not error",
       e1 < 1e-9 and e2 < 1e-4,
       f"max |dv/dg - central difference| {e1:.1e}; "
       f"max |d2v/dg2 - central difference| {e2:.1e} (floor ~{2e-16 / h ** 2:.0e})",
       "first derivative at 1e-9, second at its roundoff floor"))

    # R4 -- the assemblies close
    res = {"6-chain": float(np.abs(ch6.weld_residual(ch6.q0())).max()),
           "9-cluster": float(np.abs(cl9.weld_residual(cl9.q0())).max())}
    A(("R4  THE WELDS CLOSE AT THE INITIAL CONFIGURATION. The corner "
       "correspondence flips with the direction of the 60 degree offset and a "
       "chain alternates, so this is asked per link rather than reused. Note "
       "what this row does NOT prove: jb_ic's R2b established that agreement "
       "in POSITION does not validate IDENTITY, which is why the "
       "correspondence comes from the gated builders rather than being "
       "recomputed here",
       max(res.values()) < 1e-14,
       "; ".join(f"{k}: max |weld residual| {v:.1e}" for k, v in res.items()),
       "closed to machine precision"))

    # R5a -- the global-motion basis, and what the defective one does
    ctr, R, gam, B = ch6.frames(ch6.q0())
    J6 = ch6.cell_jacobians(ctr, R, B)
    C6 = ch6.constraint_jacobian(J6)
    good = float(np.abs(np.dot(C6, ch6.globals(ctr).T)).max())
    bad = float(np.abs(np.dot(C6, ch6.globals(ctr, defective=True).T)).max())
    A(("R5a THE GLOBAL-MOTION BASIS IS THE CORRECT ONE. Translations are "
       "(dc = t, w = 0); rotations are (dc = w x c, w = w). The omega x c term "
       "on the centres is what an earlier fold-map-rank run omitted, together "
       "with using an identity orientation block, and it reported a DOF column "
       "two too high. TWO-SIDED: the defective basis is NOT annihilated by the "
       "constraint Jacobian, and its residual is printed here, so this row "
       "cannot pass by measuring nothing",
       good < 1e-12 and bad > 1e-3,
       f"||C G|| with the correct basis {good:.1e}; with the defective basis "
       f"{bad:.3f}", "correct basis in the nullspace, defective basis not"))

    # R5b -- THE DOF TABLE, from these coordinates
    free, fixed, ranks = [], [], []
    for N in range(1, 8):
        a = chain(N)
        ct, Rm, gm, Bd = a.frames(a.q0())
        Jn = a.cell_jacobians(ct, Rm, Bd)
        Cn = a.constraint_jacobian(Jn) if a.welds else np.zeros((0, 7 * N))
        r, _ = rank_of(Cn)
        keep = [i for i in range(7 * N) if i % 7 not in (3, 4, 5)]
        rf, _ = rank_of(Cn[:, keep])
        free.append(7 * N - r - 6)
        fixed.append(4 * N - rf - 3)
        if a.welds:
            ranks.append((r / len(a.welds), rf / len(a.welds)))
    ctr9, R9, g9, B9 = cl9.frames(cl9.q0())
    J9 = cl9.cell_jacobians(ctr9, R9, B9)
    C9 = cl9.constraint_jacobian(J9)
    r9, gap9 = rank_of(C9)
    keep9 = [i for i in range(63) if i % 7 not in (3, 4, 5)]
    rf9, _ = rank_of(C9[:, keep9])
    free9, fixed9 = 63 - r9 - 6, 36 - rf9 - 3
    out["dof"] = (free, fixed, free9, fixed9)
    A(("R5b THE DOF TABLE, REPRODUCED FROM THE REDUCED COORDINATES. With cells "
       "free to turn the array has ONE internal degree of freedom PER CELL at "
       "every size -- a fold angle each, which is what makes the medium a "
       "FIELD. Freeze the orientations and the SAME Jacobian, restricted to "
       "its centre and fold columns, has rank 4 per weld instead of 6 and "
       "leaves a SINGLE coordinate at every size, so a disturbance has nowhere "
       "to be. Both rows are column subsets of one matrix, which is a stronger "
       "statement than two separate measurements",
       free == [1, 2, 3, 4, 5, 6, 7] and fixed == [1] * 7
       and free9 == 9 and fixed9 == 1
       and all(abs(a - 6.0) < 1e-9 and abs(b - 4.0) < 1e-9 for a, b in ranks),
       f"chains of 1..7 cells, may rotate: {free}; orientation fixed: {fixed}; "
       f"9-cluster {free9} and {fixed9}; weld rank {ranks[0][0]:.0f} of 9 "
       f"(orientation fixed {ranks[0][1]:.0f}), rank gap at the cluster "
       f"{gap9:.1e}",
       "[1..7] and 9 free, 1 everywhere with orientation fixed"))

    # R6 -- the five spurious DOF per cell, measured rather than asserted
    Pb, bb, _, _, _, _ = IC.build(list(ch6.gam0))
    bar_dof = 3 * len(Pb) - IC.rigidity(Pb, bb) - 6
    A(("R6  THE SPURIOUS MODES ARE GONE BY CONSTRUCTION, and the difference is "
       "measured on the SAME six-cell chain rather than quoted. The bar-only "
       "space -- strut lengths held, nothing else, which is what every earlier "
       "integration in this programme used -- carries six internal degrees of "
       "freedom per cell. Five of them are not jitterbug motions at all, and "
       "they are what the jumbled cluster render was showing",
       bar_dof == 6 * ch6.N and free[5] == ch6.N,
       f"same 6-cell chain: bar-only {bar_dof} internal DOF "
       f"({bar_dof // ch6.N} per cell), reduced coordinates {free[5]} "
       f"({free[5] // ch6.N} per cell)", "36 against 6, i.e. 6 per cell against 1"))

    # R7 -- the impulse is admissible and carries no net momentum
    u6, _ = fold_impulse(ch6, 0)
    adm = float(np.abs(np.dot(C6, u6.ravel())).max())
    mom = float(np.abs(np.dot(ch6.momentum_rows(ctr, R, B, J6), u6.ravel())).max())
    E6 = ch6.energy(ch6.q0(), u6)
    A(("R7  THE PHASE IMPULSE IS ADMISSIBLE AND INTERNAL. gdot = 1 on one cell "
       "and nothing else is not a legal state until it is projected onto the "
       "constraint tangent space; the projection also removes net linear and "
       "angular momentum, without which part of what follows would be the "
       "whole patch drifting rather than transport",
       adm < 1e-12 and mom < 1e-12 and E6.sum() > 0,
       f"||C u|| after projection {adm:.1e}, ||momentum|| {mom:.1e}, "
       f"E = {E6.sum():.6f}", "admissible, momentum free, energy nonzero"))

    # R8 -- energy audit over a long run
    rec6 = ch6.run(u6, 40.0, 200)
    Etot = np.array([r[1].sum() for r in rec6])
    drift = float(np.abs(Etot - Etot[0]).max() / Etot[0])
    wmax = max(r[3] for r in rec6)
    gspan = (min(r[2].min() for r in rec6), max(r[2].max() for r in rec6))
    out["chain"] = rec6
    A(("R8  ENERGY IS CONSERVED WITH V = 0, AND THE WELDS STAY CLOSED. The "
       "integrator is DELIBERATELY NOT symplectic: a symplectic scheme "
       "conserves a shadow Hamiltonian and would make this row nearly vacuous, "
       "whereas DOP853 has no such protection and would drift visibly if the "
       "equations of motion were wrong. The weld residual is reported "
       "UNPROJECTED -- nothing is put back onto the manifold between samples, "
       "so it is a second independent check on the same equations",
       drift < 1e-8 and wmax < 1e-8,
       f"E = {Etot[0]:.9f} held to {drift:.1e} over t = 0..40 (200 samples); "
       f"max weld residual {wmax:.1e}, unprojected; fold angles stayed in "
       f"[{gspan[0]:.1f}, {gspan[1]:.1f}] degrees",
       "relative drift < 1e-8, weld residual < 1e-8"))

    # R9 -- transport along the chain
    sh = np.array([r[1] / r[1].sum() for r in rec6])
    ts = np.array([r[0] for r in rec6])
    i0 = int(sh[:, 0].argmin())
    back = i0 + int(sh[i0:, 0].argmax())
    j1 = int(sh[:, 1].argmin())
    j1b = j1 + int(sh[j1:, 1].argmax())
    inner = float(sh[:, 1:-1].max())
    A(("R9  THE CHAIN TRANSPORTS, AND IT COMES BACK. The driven cell gives up "
       "most of its share, an interior cell ends up holding more than the "
       "driven one has left, and cell 1 empties almost completely and refills. "
       "This is the six-cell propagation-and-reflection run re-measured in the "
       "reduced coordinates; the bar-space numbers it replaces are withdrawn, "
       "not carried over. NOT A FRONT: R10 shows the onset is simultaneous, so "
       "this is mode superposition and no arrival time here means a signal "
       "speed",
       sh[i0, 0] < 0.45 * sh[0, 0] and inner > 0.30
       and sh[j1, 1] < 0.005 and sh[j1b, 1] > 0.20
       and sh[:, -1].max() > 1.5 * sh[0, -1],
       f"driven cell {sh[0, 0]:.4f} at t=0, down to {sh[i0, 0]:.4f} at "
       f"t={ts[i0]:.1f}, back to {sh[back, 0]:.4f} by t={ts[back]:.1f}; the "
       f"hottest interior cell reaches {inner:.4f}; cell 1 empties to "
       f"{sh[j1, 1]:.5f} at t={ts[j1]:.1f} and refills to {sh[j1b, 1]:.4f} at "
       f"t={ts[j1b]:.1f}; the far end goes {sh[0, -1]:.4f} -> "
       f"{sh[:, -1].max():.4f}",
       "driven cell below 45% of its share, an interior cell above 0.30, "
       "cell 1 empties below 0.005 and refills above 0.20"))

    # R10 -- the cluster, and what the constraint does before any time passes
    u9, _ = fold_impulse(cl9, 0)
    rec9 = cl9.run(u9, 40.0, 200)
    sh9 = np.array([r[1] / r[1].sum() for r in rec9])
    E9 = np.array([r[1].sum() for r in rec9])
    d9 = float(np.abs(E9 - E9[0]).max() / E9[0])
    w9 = max(r[3] for r in rec9)
    out["cluster"] = rec9
    A(("R10 THE CONSTRAINT IS IMMEDIATE, AND ITS SPLIT IS EXACT. At t = 0, "
       "before any time has passed, the projection that makes a phase kick on "
       "the centre admissible has already put 32/37 of the energy in the shell "
       "-- one cell cannot fold without its eight neighbours folding, and a "
       "rigid constraint has infinite signal speed, so there is no onset lag "
       "to find and no wavefront to time. The split is exactly rational, "
       "5/37 to the centre and 4/37 to each neighbour, which is a check on the "
       "projection as well as a result. The bar-space 80/20 is withdrawn, not "
       "rescaled into this",
       abs(sh9[0, 0] - 5.0 / 37.0) < 1e-12
       and float(np.abs(sh9[0, 1:] - 4.0 / 37.0).max()) < 1e-12
       and abs(float(u9[0, 6]) - 5.0 / 37.0) < 1e-12
       and d9 < 1e-8 and w9 < 1e-8,
       f"centre {sh9[0, 0]:.9f} against 5/37 = {5 / 37:.9f}; each of the eight "
       f"shell cells {sh9[0, 1]:.9f} against 4/37 = {4 / 37:.9f}; projected "
       f"fold rates 5/37 and 1/37; E = {E9[0]:.9f} held to {d9:.1e} over "
       f"t = 0..40, max weld residual {w9:.1e}",
       "exactly 5/37 : 32/37 at t=0, energy conserved"))

    spread = max(float(s[1:].max() - s[1:].min()) for s in sh9)
    rise = float(np.diff(sh9[:, 0]).max())
    A(("R10b THE DRAIN-AND-RETURN DOES NOT REPRODUCE UNDER A SYMMETRIC KICK, "
       "and this corrects the framing of the claim it replaces rather than "
       "restating it with new numbers. A phase kick on the centre is invariant "
       "under the cube group, so it excites ONLY the symmetric mode: all eight "
       "shell cells stay identical to machine precision for the whole run, and "
       "the centre drains monotonically with no return anywhere in t = 0..40. "
       "The earlier drain-and-return used a spin on ONE triangle, which "
       "selects a face and breaks the symmetry -- so it was measuring a "
       "different initial condition, not just a different constraint space",
       spread < 1e-10 and int(sh9[:, 0].argmax()) == 0 and rise < 1e-9
       and sh9[-1, 0] < 0.2 * sh9[0, 0],
       f"spread across the eight shell cells, worst over the run "
       f"{spread:.1e}; centre {sh9[0, 0]:.5f} -> {sh9[-1, 0]:.5f} (factor "
       f"{sh9[0, 0] / sh9[-1, 0]:.1f}); largest sample-to-sample CHANGE "
       f"{rise:.1e}, negative throughout, so strictly decreasing",
       "shell cells identical, centre monotone down by more than 5x"))

    u9s, _ = fold_impulse(cl9, 1)
    rec9s = cl9.run(u9s, 40.0, 200)
    sh9s = np.array([r[1] / r[1].sum() for r in rec9s])
    ts9s = np.array([r[0] for r in rec9s])
    ks = int(sh9s[:, 1].argmin())
    ksb = ks + int(sh9s[ks:, 1].argmax())
    out["shell"] = rec9s
    A(("R10c BREAK THE SYMMETRY AND THE RETURN IS THERE. The same cluster, the "
       "same kind of kick, applied to a SHELL cell instead of the centre: the "
       "driven cell drains and refills. So the return is a property of "
       "asymmetric initial conditions, not of the patch, and R10b is a "
       "statement about the symmetric mode rather than about the medium",
       sh9s[ks, 1] < 0.9 * sh9s[0, 1] and sh9s[ksb, 1] > 1.1 * sh9s[ks, 1],
       f"driven shell cell {sh9s[0, 1]:.5f} at t=0, down to {sh9s[ks, 1]:.5f} "
       f"at t={ts9s[ks]:.1f}, back to {sh9s[ksb, 1]:.5f} at "
       f"t={ts9s[ksb]:.1f}; the centre takes at most {sh9s[:, 0].max():.5f}",
       "the driven cell drains and refills"))

    # R11 -- mutation probes, one per row that could otherwise pass vacuously
    probes = []

    cl_drop = cluster(drop=(0,))
    ctd, Rd, gd, Bd = cl_drop.frames(cl_drop.q0())
    Cd = cl_drop.constraint_jacobian(cl_drop.cell_jacobians(ctd, Rd, Bd))
    rd, _ = rank_of(Cd)
    probes.append(("drop a weld -> R5b", 63 - rd - 6, free9, (63 - rd - 6) != free9))

    def _cycle(pairs):
        """Send each corner to its NEIGHBOUR's partner. Note that reversing
        the list of pairs is not a mutation at all -- the same three
        correspondences in another order -- and gives 4e-16; the corner map
        itself has to change."""
        tgt = [b for _, b in pairs]
        return [(a, tgt[(i + 1) % 3]) for i, (a, _) in enumerate(pairs)]

    cl_bad = Assembly(cl9.gam0, cl9.ctr0,
                      [(k, l, _cycle(pairs)) for k, l, pairs in cl9.welds])
    rbad = float(np.abs(cl_bad.weld_residual(cl_bad.q0())).max())
    probes.append(("permute a weld's corners -> R4", f"{rbad:.2e}", "< 1e-14",
                   rbad > 1e-9))

    ctf, Rf, gf, Bf = ch6.frames(ch6.q0())
    Jf = ch6.cell_jacobians(ctf, Rf, Bf)
    Cf = np.vstack([ch6.constraint_jacobian(Jf),
                    ch6.momentum_rows(ctf, Rf, Bf, Jf),
                    np.eye(7 * ch6.N)[[7 * k + 6 for k in range(ch6.N)]]])
    Mf = ch6.mass_blocks(Jf)
    MiCf = np.zeros((7 * ch6.N, Cf.shape[0]))
    for k in range(ch6.N):
        MiCf[7 * k:7 * k + 7] = np.linalg.solve(Mf[k], Cf[:, 7 * k:7 * k + 7].T)
    u_frozen = np.zeros((ch6.N, 7))
    u_frozen[0, 6] = 1.0
    lamf = _solve(Cf, MiCf, -np.dot(Cf, u_frozen.ravel()))
    uf = u_frozen + np.dot(MiCf, lamf).reshape(ch6.N, 7)
    Ef = float(ch6.energy(ch6.q0(), uf).sum())
    probes.append(("delete the fold coordinate -> R7", f"E = {Ef:.1e}",
                   f"E = {E6.sum():.4f}", Ef < 1e-20))

    ch_cut = chain(6, drop=(2,))
    u_cut, _ = fold_impulse(ch_cut, 0)
    rec_cut = ch_cut.run(u_cut, 40.0, 100)
    sh_cut = np.array([r[1] / r[1].sum() for r in rec_cut])
    far_cut = sh_cut[:, 3:].sum(axis=1)
    far_int = sh[:, 3:].sum(axis=1)
    probes.append(("sever the chain between cells 2 and 3 -> R9",
                   f"far fragment varies {far_cut.max() - far_cut.min():.1e}",
                   f"intact varies {far_int.max() - far_int.min():.4f}",
                   (far_cut.max() - far_cut.min()) < 1e-8))

    A(("R11 MUTATION PROBES. Four, because one probe only proves one row is "
       "non-vacuous. Dropping a weld changes the DOF count. Permuting a weld's "
       "corner correspondence opens the residual R4 gates -- and note it does "
       "so only because the permutation is a genuine mismatch, whereas "
       "jb_ic's R2b defect agreed in POSITION to 3e-15 while the identities "
       "were wrong, which is why R4 leans on the gated builders rather than on "
       "its own arithmetic. Deleting the fold coordinate leaves the phase kick "
       "with nowhere to go at all: the admissible projection of it is "
       "identically zero, so there is no disturbance to transport, which is "
       "the same statement as the orientation-fixed column of R5b arrived at "
       "from the other side. Severing the chain stops transport dead -- "
       "disconnected fragments cannot exchange energy, so the far fragment's "
       "share is frozen where the intact chain's moves",
       all(p[3] for p in probes),
       "; ".join(f"{n}: {g} vs {w} -> {'reddens' if ok else 'DID NOT REDDEN'}"
                 for n, g, w, ok in probes),
       "every probe reddens its row"))

    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_rc -- the medium in reduced coordinates, ellipses by construction")
    print("=" * 78)
    checks, out = gate()
    bad = 0
    for name, ok, got, want in checks:
        tag = "PASS" if ok else "FAIL"
        bad += 0 if ok else 1
        print(f"  {tag}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")

    free, fixed, free9, fixed9 = out["dof"]
    print("\n  INTERNAL DEGREES OF FREEDOM, global rigid motions removed with a")
    print("  correct (omega x c) basis:")
    print(f"    {'cells':<22}" + "".join(f"{n:>4}" for n in range(1, 8)) + "   9-cluster")
    print(f"    {'orientation FIXED':<22}" + "".join(f"{n:>4}" for n in fixed) + f"{fixed9:>12}")
    print(f"    {'cells MAY ROTATE':<22}" + "".join(f"{n:>4}" for n in free) + f"{free9:>12}")

    print("\n  6-cell chain, phase-rate kick on cell 0 -- share of kinetic energy:")
    rec = out["chain"]
    print("  " + f"{'t':>7}" + "".join(f"{'cell' + str(k):>9}" for k in range(6)))
    for t, ke, _g, _w in rec[::max(1, len(rec) // 12)]:
        f = ke / ke.sum()
        print(f"  {t:7.2f}" + "".join(f"{x:9.5f}" for x in f))

    print("\n  9-cell cluster, phase-rate kick on the CENTRE. Symmetric, so the")
    print("  eight shell cells stay identical and the centre only drains:")
    rec = out["cluster"]
    print(f"  {'t':>7} {'centre':>10} {'shell':>10} {'per shell cell':>16}")
    for t, ke, _g, _w in rec[::max(1, len(rec) // 10)]:
        f = ke / ke.sum()
        print(f"  {t:7.2f} {f[0]:10.5f} {1 - f[0]:10.5f} {f[1:].max():16.5f}")

    print("\n  Same cluster, same kick on a SHELL cell. The symmetry is broken,")
    print("  and the driven cell drains and refills:")
    rec = out["shell"]
    print(f"  {'t':>7} {'driven':>10} {'centre':>10} {'hottest other':>15}")
    for t, ke, _g, _w in rec[::max(1, len(rec) // 10)]:
        f = ke / ke.sum()
        print(f"  {t:7.2f} {f[1]:10.5f} {f[0]:10.5f} "
              f"{np.delete(f, [0, 1]).max():15.5f}")

    print()
    print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
    print("   * The ellipses are enforced by the COORDINATES, so the five")
    print("     spurious modes per cell the bar-only space admitted are absent")
    print("     rather than suppressed (R6).")
    print("   * The medium is a FIELD: one internal degree of freedom per cell")
    print("     at every size, and exactly one in total if cells may not turn.")
    print("   * It licenses NO signal speed. The onset is simultaneous")
    print("     everywhere (R10), so there is no wavefront to time. A finite")
    print("     signal speed would need compliant constraints, which is the")
    print("     fork this programme has rejected (T2 inviscid 23562).")
    print("   * V = 0 still. Nothing here finds a potential energy, and no")
    print("     frequency, speed, band or gap is quoted or implied. It is")
    print("     also why the fold angles WIND without bound in these runs")
    print("     -- hundreds of degrees over t = 40 -- instead of oscillating:")
    print("     nothing restores a cell's phase. What redistributes is")
    print("     energy, not a bounded displacement, and the open problem is")
    print("     still finding V for the rigid linkage.")
    print("   * Mass model DECLARED: unit mass per triangle, lumped m/3 to each")
    print("     corner, matching jb_ic so the two are comparable. The")
    print("     uniform-lamina alternative moved a period by 7% elsewhere in")
    print("     this programme, so no number here is model-independent.")
    print()
    print("  ALL CHECKS PASSED." if not bad else f"  {bad} CHECK(S) FAILED.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
