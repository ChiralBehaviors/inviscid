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

WHAT THE COUNT SAYS, AND THE SIZE OF WHAT IT DOES NOT. For a TREE of welds,
7N - 6(N-1) = N + 6, and taking out the six global rigid motions leaves N: one
fold angle per cell, at every size. Freeze the orientations and the same
Jacobian -- restricted to its centre and fold columns -- has rank 4 per weld
instead of 6, leaving 4N - 4(N-1) - 3 = 1: a single coordinate, and nowhere for
a disturbance to be.

    THE FIRST OF THOSE IS A TREE RESULT AND IS NOT TRUE OF THE HONEYCOMB.
    Read the qualifier in the first sentence as load bearing, because for a
    while it was not read at all. Chains have no cycles and the 9-cluster is a
    star, so BOTH patches the field claim was measured on are trees. The
    rectified cubic honeycomb is not: brick 5x5x5 is 35 cells, 64 welds and 30
    independent cycles. With cycles the count does not grow with N. It is
    1 + (number of degree-1 cells), so a COMPACT patch of 113 cells has ONE
    internal degree of freedom -- coherent breathing -- and a disturbance has
    nowhere to be in the interior after all. R5e measures it. Bead
    inviscid-qvf.26 and T2 inviscid [23595] carry the finding and the audit.

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
        # ROW OFFSET PER WELD. `honeycomb` welds a shared TRIANGULAR face and
        # always carries three pairs, so this was `9 * len(welds)` and the
        # offsets were `9 * r`. `honeycomb_single` welds AXIS neighbours and
        # carries TWO (DECISION 21), so the layout is a running sum. For any
        # all-three-pair assembly `_woff[r] == 9 * r` and nothing moves.
        self._woff, off = [], 0
        for (_k, _l, _pairs) in welds:
            self._woff.append(off)
            off += 3 * len(_pairs)
        self.nc = off

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
                row = self._woff[r] + 3 * m
                g[row:row + 3] = X[k][a] - X[l][b]
        return g

    def constraint_jacobian(self, J):
        C = np.zeros((self.nc, 7 * self.N))
        for r, (k, l, pairs) in enumerate(self.welds):
            for m, (a, b) in enumerate(pairs):
                row = self._woff[r] + 3 * m
                C[row:row + 3, 7 * k:7 * k + 7] = J[k][3 * a:3 * a + 3]
                C[row:row + 3, 7 * l:7 * l + 7] = -J[l][3 * b:3 * b + 3]
        return C

    def constraint_gyro(self, A):
        d = np.zeros(self.nc)
        for r, (k, l, pairs) in enumerate(self.welds):
            for m, (a, b) in enumerate(pairs):
                row = self._woff[r] + 3 * m
                d[row:row + 3] = A[k][a] - A[l][b]
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


def honeycomb(sites, gc=-30.0):
    """A patch of the rectified cubic honeycomb, with its CYCLES.

    Integer `sites` scaled by the breathing lattice constant L: all-EVEN sites
    carry fold angle `gc`, all-ODD sites carry `gc + 60`, and two cells are
    welded exactly when their integer difference is (+-1, +-1, +-1) -- the eight
    triangular-face neighbours. (The other six of the 6+8 census are SQUARE-face
    contacts, which are open at every angle but 0 and are not welds.)

    `cluster` builds a star and `chain` builds a line; both are TREES, and the
    difference is not cosmetic. This builder is the one that closes cycles, and
    closing them is what takes the internal freedom away.
    """
    gn = gc + 60.0
    sep = IC.ZC * (np.cos(np.radians(gc)) + np.cos(np.radians(gn)))
    L = sep / np.sqrt(3)
    gam = [gc if all(c % 2 == 0 for c in s) else gn for s in sites]
    ctr = [L * np.array(s, float) for s in sites]
    welds, deg = [], [0] * len(sites)
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            d = tuple(sites[j][t] - sites[i][t] for t in range(3))
            if sorted(map(abs, d)) != [1, 1, 1]:
                continue
            deg[i] += 1
            deg[j] += 1
            A = IC.cell_verts(gam[i], ctr[i])
            B = IC.cell_verts(gam[j], ctr[j])
            dv = np.array(d, float)
            fa, fb = CL.face_along(dv), CL.face_along(-dv)
            bc = [SLOT[(fb, c)] for c in range(3)]
            used, pairs = set(), []
            for c in range(3):
                a = SLOT[(fa, c)]
                k = min([b for b in bc if b not in used],
                        key=lambda b: float(np.linalg.norm(B[b] - A[a])))
                used.add(k)
                pairs.append((a, k))
            welds.append((i, j, pairs))
    return Assembly(gam, ctr, welds), deg


#: The phase the SINGLE covering's axis welds are read at, and held from.
#: The axis contact's arity is phase dependent -- four coincident vertex pairs
#: at a = 0, TWO through the interior, one at a = -60 -- because the square is
#: an OPENING that closes as the exchange passes through it, not a face. Only
#: the interior pair set is valid at every phase; jb_ht T2 gates exactly that
#: and it is measured here in R2. Same constant, same reason, as
#: `jb_w_honeycomb.HONEYCOMB_REF_PHASE`.
SINGLE_REF_PHASE = -30.0


def honeycomb_single(sites, gc=-30.0, ref=SINGLE_REF_PHASE):
    """A patch of the honeycomb with ONE COVERING: the voids stay EMPTY.

    OWNER DECISION 21, 2026-08-29 (T2 [23727]). The octahedral spaces between
    the cells are EMPTY. Through the jitterbug exchange the SOLID cells run
    VE -> octahedron while the VOIDS run octahedron -> VE: the shapes swap and
    the OCCUPANCY does not, so one sublattice is solid at every phase and the
    other is void at every phase.

    `honeycomb` puts a cell at EVERY site, which fills the voids. Because each
    triangular face is shared by one even cell and one odd cell and both carry
    a plate there, that builder draws every interior triangle TWICE (measured:
    HC9 72 plates / 64 distinct, box r=3 728 / 512, ratio -> 2 in the bulk).
    This builder keeps only the ALL-EVEN sublattice, so every triangle is
    drawn once.

    THE WELDS ARE THEREFORE DIFFERENT, and that is the whole structural
    consequence. With the voids empty a cell's nearest sites are voids, so it
    welds to its six AXIS neighbours -- the second-nearest cells, at
    (+-2, 0, 0) and its permutations -- and each of those carries TWO
    coincident vertex pairs rather than a triangular face's three. The pairs
    are read at `ref` and HELD, because their count is phase dependent while
    that pair set is not.

    Returns (Assembly, degree) exactly as `honeycomb` does. `honeycomb` is
    left intact: the two are meant to be compared module by module, not
    swapped blind.
    """
    sites = [tuple(int(c) for c in s) for s in sites]
    solid = [s for s in sites if all(c % 2 == 0 for c in s)]
    if not solid:
        raise ValueError("honeycomb_single: no all-even sites in `sites`")

    def _ctr(g, ss):
        sep = IC.ZC * (np.cos(np.radians(g)) + np.cos(np.radians(g + 60.0)))
        L = sep / np.sqrt(3)
        return [L * np.array(t, float) for t in ss]

    # Correspondence read ONCE, at `ref`, and held -- never re-read at gc.
    probe = Assembly([ref] * len(solid), _ctr(ref, solid), [])
    Xr = probe.positions(probe.q0())
    welds, deg = [], [0] * len(solid)
    for i in range(len(solid)):
        for j in range(i + 1, len(solid)):
            d = tuple(solid[j][t] - solid[i][t] for t in range(3))
            if sorted(map(abs, d)) != [0, 0, 2]:
                continue
            # ONE RECORD PER DISTINCT SHARED POINT. At the collapsed phases a
            # cell's twelve labels occupy six positions, so several label
            # pairs name the SAME joint -- at a = -60 the axis bond has four
            # label pairs and ONE distinct point -- and emitting each would
            # duplicate that joint's constraint fourfold. Same dedup, same
            # reason, as `jb_w_honeycomb.honeycomb_contacts`' square branch;
            # without it this builder is only correct at a reference phase
            # where no vertices have merged.
            pairs, seen = [], []
            for a in range(NV):
                v = Xr[i][a]
                if any(np.linalg.norm(v - u) < 1e-9 for u in seen):
                    continue
                for b in range(NV):
                    if np.linalg.norm(v - Xr[j][b]) < 1e-9:
                        pairs.append((a, b))
                        seen.append(v)
                        break
            if not pairs:
                continue
            deg[i] += 1
            deg[j] += 1
            welds.append((i, j, pairs))
    return Assembly([gc] * len(solid), _ctr(gc, solid), welds), deg


def brick(nx, ny, nz):
    return [(x, y, z) for x in range(nx) for y in range(ny) for z in range(nz)
            if all(c % 2 == 0 for c in (x, y, z))
            or all(c % 2 == 1 for c in (x, y, z))]


def ball(radius):
    n = int(radius) + 2
    return [(x, y, z) for x in range(-n, n + 1) for y in range(-n, n + 1)
            for z in range(-n, n + 1)
            if (all(c % 2 == 0 for c in (x, y, z))
                or all(c % 2 == 1 for c in (x, y, z)))
            and x * x + y * y + z * z <= radius * radius]


def shared_vertices(asm, q=None):
    """Vertex INSTANCES that coincide, grouped into classes.

    The honeycomb shares vertices between cells that do not share a FACE: in the
    bulk a vertex is met by four cells, and a face-weld identifies only the three
    corners of one triangle. What is left over is not a numerical coincidence --
    R5f checks the class structure is identical at five fold angles -- but it is
    only ever imposed here, never by `honeycomb`'s welds.
    """
    X = asm.positions(asm.q0() if q is None else q)
    buck = {}
    for k in range(asm.N):
        for i in range(NV):
            buck.setdefault(tuple(np.round(X[k, i], 6)), []).append((k, i))
    return sorted(tuple(sorted(v)) for v in buck.values() if len(v) > 1)


def identification_rows(asm, classes):
    """Constraint rows forcing every instance in a class to stay coincident:
    3(m-1) rows for a class of m, spanning it rather than over-specifying it."""
    ctr, R, gam, B = asm.frames(asm.q0())
    J = asm.cell_jacobians(ctr, R, B)
    rows = []
    for cls in classes:
        ka, ia = cls[0]
        for kb, ib in cls[1:]:
            for t in range(3):
                r = np.zeros(7 * asm.N)
                r[7 * ka:7 * ka + 7] = J[ka][3 * ia + t]
                r[7 * kb:7 * kb + 7] -= J[kb][3 * ib + t]
                rows.append(r)
    return np.array(rows)


def apply_increment(asm, q, d):
    """Move a configuration by a generalized displacement, rotating each
    orientation by the exponential map so it stays on SO(3)."""
    q = q.copy()
    for k in range(asm.N):
        q[8 * k:8 * k + 3] += d[7 * k:7 * k + 3]
        w = d[7 * k + 3:7 * k + 6]
        th = float(np.linalg.norm(w))
        if th > 0:
            dq = np.concatenate(([np.cos(th / 2)], np.sin(th / 2) * w / th))
            a, b = dq, q[8 * k + 3:8 * k + 7]
            q[8 * k + 3:8 * k + 7] = np.array([
                    a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
                    a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
                    a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
                    a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]])
        q[8 * k + 7] += np.degrees(d[7 * k + 6])
        q[8 * k + 3:8 * k + 7] /= np.linalg.norm(q[8 * k + 3:8 * k + 7])
    return q


def walk(asm, q0, u, eps):
    """A FINITE step along a generalized velocity, Gauss-Newton'd back onto the
    weld manifold. Finite rather than infinitesimal on purpose: whether a contact
    opens or closes is second order, so a tangent alone cannot answer it."""
    q = apply_increment(asm, q0, eps * u)
    for _ in range(30):
        g = asm.weld_residual(q)
        if np.abs(g).max() < 1e-13:
            break
        ctr, R, gam, B = asm.frames(q)
        C = asm.constraint_jacobian(asm.cell_jacobians(ctr, R, B))
        # np.dot, not `@` -- see _solve on the Accelerate BLAS status flags.
        q = apply_increment(asm, q, -np.dot(
                C.T, np.linalg.lstsq(np.dot(C, C.T), g, rcond=None)[0]))
    return q


#: A cell folded past its own octahedron self-intersects, so |gamma| <= 60 is a
#: HARD STOP, not a convention. The shared-face law b = a + 60 forces BOTH
#: sublattices into that window at once, so the medium's admissible range is
#: a in [-60, 0] -- sixty degrees, over which the sublattices exchange roles
#: completely and the lattice constant breathes 1 -> 2/sqrt(3) -> 1.
FOLD_LIMIT = 60.0


def lattice_constant(a_deg):
    """L(a): the centre-to-centre spacing per unit site, from the shared-face
    law. Stationary at a = -30, which is the only angle where a coherently
    breathing patch's centres do not move."""
    return (IC.ZC * (np.cos(np.radians(a_deg))
                     + np.cos(np.radians(a_deg + 60.0))) / np.sqrt(3.0))


def coherent(sites, a_deg, rate=1.0):
    """THE MEDIUM'S ONE MOTION, in closed form.

    A compact patch has exactly one internal degree of freedom, and R5h measures
    what it is: every cell folds at the SAME rate, NO cell rotates, and every
    centre rides the lattice constant. So the whole configuration and its whole
    velocity are functions of one angle:

        gamma_even = a,  gamma_odd = a + 60,  R = I,  centre = L(a) * site
        cdot = dL/da * adot * site,  omega = 0,  gammadot = adot

    Returned as (q, u) in the Assembly's own packing, so it can be checked
    against the constrained model rather than trusted -- which is what R5h does.
    `rate` is adot in RADIANS per unit time.
    """
    n = len(sites)
    L, dL, _ = lattice_derivatives(a_deg)
    q = np.zeros((n, 8))
    u = np.zeros((n, 7))
    for k, site in enumerate(sites):
        even = all(c % 2 == 0 for c in site)
        q[k, 0:3] = L * np.array(site, float)
        q[k, 3] = 1.0
        q[k, 7] = a_deg if even else a_deg + 60.0
        u[k, 0:3] = dL * rate * np.array(site, float)
        u[k, 6] = rate
    return q.ravel(), u


def lattice_derivatives(a_deg):
    """L, dL/da and d2L/da2 in RADIANS. L is a sum of two cosines, so
    d2L/da2 = -L exactly -- worth writing down rather than differencing, because
    a finite-differenced dM/da makes the equation of motion noisy at ~1e-8 and an
    adaptive integrator asked for 1e-11 will thrash against it forever. That is
    not hypothetical; it is what the first version of `swing` did."""
    r = np.radians(a_deg)
    z = IC.ZC / np.sqrt(3.0)
    L = z * (np.cos(r) + np.cos(r + np.pi / 3.0))
    return L, -z * (np.sin(r) + np.sin(r + np.pi / 3.0)), -L


def effective_mass(sites, a_deg, derivative=False):
    """M_eff(a) for a coherently breathing patch: twice the kinetic energy at
    unit fold rate, with dM/da on request.

    Summed straight over the cells rather than through an `Assembly`, because the
    coherent velocity is closed form and needs no weld structure at all -- a
    vertex moves at dL/da * site + v'(gamma). Going through the assembly per call
    instead was two orders of magnitude slower. R5h checks this against the
    constrained model's own energy, so the shortcut is verified rather than
    assumed.

    M_eff is NOT constant in a, which is why the breathe is not uniform in time
    and why an animation of it is dynamics rather than a sweep.
    """
    _, dL, ddL = lattice_derivatives(a_deg)
    ev = body(a_deg, 2)
    od = body(a_deg + 60.0, 2)
    m, dm = 0.0, 0.0
    for site in sites:
        s = np.array(site, float)
        d = ev if all(c % 2 == 0 for c in site) else od
        vel = dL * s + d[1]
        acc = ddL * s + d[2]
        m += float(np.dot(VMASS, np.einsum('ij,ij->i', vel, vel)))
        if derivative:
            dm += 2.0 * float(np.dot(VMASS, np.einsum('ij,ij->i', vel, acc)))
    return (m, dm) if derivative else m


def swing(sites, a0, adot0, tmax, nsample):
    """The 1-DOF motion, bounced elastically off both fold limits.

    V = 0, so the only force is the configuration-dependent inertia:
    addot = -(1/2)(M'/M) adot^2, which is FreeDynamics' equation with this
    patch's effective mass. The impact is velocity reversal, and that is exact
    rather than convenient: with one degree of freedom every admissible velocity
    is a multiple of one mode, so reversing it conserves energy identically and
    cannot leave the weld tangent space.
    """
    from scipy.integrate import solve_ivp

    def rhs(_t, y):
        a, ad = y
        m, dm = effective_mass(sites, a, derivative=True)
        return [np.degrees(ad), -0.5 * (dm / m) * ad * ad]

    def hi(_t, y):
        return -y[0]

    def lo(_t, y):
        return y[0] + FOLD_LIMIT

    hi.terminal = lo.terminal = True
    #: DIRECTION MATTERS, and omitting it hangs. After a bounce the state sits
    #: exactly ON the limit, so an undirected event fires again at once, at zero
    #: elapsed time, forever. Requiring a downward crossing means the outgoing
    #: leg -- which moves away from the limit -- does not retrigger it.
    hi.direction = lo.direction = -1.0
    y = [a0, adot0]
    ts = np.linspace(0.0, tmax, nsample + 1)
    rec, bounces, t0 = [], 0, 0.0
    while t0 < tmax - 1e-12:
        want = ts[(ts > t0 - 1e-12) & (ts <= tmax)]
        sol = solve_ivp(rhs, (t0, tmax), y, method="DOP853", t_eval=want,
                        rtol=1e-11, atol=1e-13, events=(hi, lo))
        for j, t in enumerate(sol.t):
            a, ad = sol.y[0, j], sol.y[1, j]
            rec.append((t, a, ad,
                        0.5 * effective_mass(sites, a) * ad * ad, bounces))
        if sol.status != 1:
            break
        t0 = float(sol.t_events[0][0] if len(sol.t_events[0])
                   else sol.t_events[1][0])
        ye = (sol.y_events[0][0] if len(sol.y_events[0])
              else sol.y_events[1][0])
        y = [float(ye[0]), -float(ye[1])]      # the elastic impact, in full
        bounces += 1
    return rec, bounces


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
    A(("R5b THE DOF TABLE FOR A TREE OF WELDS -- WHICH IS WHAT A CHAIN AND A "
       "SINGLE-SHELL CLUSTER ARE, AND WHICH THE HONEYCOMB IS NOT. Read with "
       "R5e or not at all: with cells free to turn a TREE has one internal "
       "degree of freedom per cell at every size, and this was recorded as the "
       "field claim before anyone counted the cycles in the real lattice. "
       "Freeze the orientations and the SAME Jacobian, restricted to "
       "its centre and fold columns, has rank 4 per weld instead of 6 and "
       "leaves a SINGLE coordinate at every size. Both rows are column subsets "
       "of one matrix, which is a stronger statement than two separate "
       "measurements -- and R5e reaches the frozen column's answer WITHOUT "
       "freezing anything, just by closing the cycles",
       free == [1, 2, 3, 4, 5, 6, 7] and fixed == [1] * 7
       and free9 == 9 and fixed9 == 1
       and all(abs(a - 6.0) < 1e-9 and abs(b - 4.0) < 1e-9 for a, b in ranks),
       f"chains of 1..7 cells, may rotate: {free}; orientation fixed: {fixed}; "
       f"9-cluster {free9} and {fixed9}; weld rank {ranks[0][0]:.0f} of 9 "
       f"(orientation fixed {ranks[0][1]:.0f}), rank gap at the cluster "
       f"{gap9:.1e}",
       "[1..7] and 9 free, 1 everywhere with orientation fixed"))

    # R5c -- the three-cell V, which is what the animation builds
    vee = cluster(sites=[(1, 1, 1), (1, 1, -1)])
    ctrv, Rv, gv, Bv = vee.frames(vee.q0())
    Jv = vee.cell_jacobians(ctrv, Rv, Bv)
    Cv = vee.constraint_jacobian(Jv)
    rv, _ = rank_of(Cv)
    keepv = [i for i in range(21) if i % 7 not in (3, 4, 5)]
    rvf, _ = rank_of(Cv[:, keepv])
    n1 = ctrv[1] / np.linalg.norm(ctrv[1])
    n2 = ctrv[2] / np.linalg.norm(ctrv[2])
    uv, _ = fold_impulse(vee, 0)
    Ev = vee.energy(vee.q0(), uv)
    uvo, _ = fold_impulse(vee, 1)
    recv = vee.run(uvo, 30.0, 300)
    shv = np.array([r[1] / r[1].sum() for r in recv])
    tsv = np.array([r[0] for r in recv])
    mid, far = int(shv[:, 0].argmax()), int(shv[:, 2].argmax())
    out["vee"] = (float(uv[0, 6]), float(uv[1, 6]), float(Ev.sum()),
                  float(shv[mid, 0]), float(tsv[mid]),
                  float(shv[far, 2]), float(tsv[far]))
    A(("R5c THE THREE-CELL V, and its numbers are the reference the JAVA port is "
       "checked against. It is the configuration ThreeCellAnimation builds -- a "
       "centre cell and two neighbours down ADJACENT diagonals, so the outer "
       "cells sit at n1.n2 = 1/3 -- and in these coordinates it has THREE "
       "internal degrees of freedom where that animation's single coherent "
       "angle had one. A phase kick on the centre projects to EXACTLY RATIONAL "
       "fold rates, 57/137 and 12/137 at E = 152/137, and an exact rational is "
       "worth more than a decimal as a cross-implementation check: a wrong port "
       "does not land near a different rational, it lands on one",
       abs(float(np.dot(n1, n2)) - 1 / 3) < 1e-14
       and 21 - rv - 6 == 3 and 12 - rvf - 3 == 1
       and abs(float(uv[0, 6]) - 57 / 137) < 1e-12
       and abs(float(uv[1, 6]) - 12 / 137) < 1e-12
       and abs(float(uv[2, 6]) - 12 / 137) < 1e-12
       and abs(float(Ev.sum()) - 152 / 137) < 1e-12
       and float(np.abs(vee.weld_residual(vee.q0())).max()) < 1e-14,
       f"n1.n2 = {float(np.dot(n1, n2)):.15f}; internal DOF {21 - rv - 6} free, "
       f"{12 - rvf - 3} with orientation fixed; centre-kick fold rates "
       f"{uv[0, 6]:.9f} and {uv[1, 6]:.9f} against 57/137 = {57 / 137:.9f} and "
       f"12/137 = {12 / 137:.9f}; E = {Ev.sum():.9f} against 152/137 = "
       f"{152 / 137:.9f}; shares {Ev[0] / Ev.sum():.6f} / {Ev[1] / Ev.sum():.6f} "
       f"/ {Ev[2] / Ev.sum():.6f}",
       "n1.n2 = 1/3, 3 free and 1 fixed, rates 57/137 and 12/137 at 152/137"))

    A(("R5d TRANSPORT ACROSS THE V, which is what the animation shows and what "
       "the Java test pins. Kick ONE OUTER cell instead of the centre and the "
       "disturbance reaches the far one THROUGH the middle. Two-sided on the far "
       "cell: it has to end up with several times the share the projection gave "
       "it at t = 0, or there is no transport to look at. Still NOT a front -- "
       "R10's onset is simultaneous, so no time on this trajectory is an "
       "arrival time",
       shv[mid, 0] > 3 * shv[0, 0] and shv[far, 2] > 3 * shv[0, 2],
       f"driven outer cell fold rate {uvo[1, 6]:.9f}; middle cell "
       f"{shv[0, 0]:.5f} at t=0, peaks {shv[mid, 0]:.5f} at t={tsv[mid]:.2f}; "
       f"far cell {shv[0, 2]:.5f} at t=0, peaks {shv[far, 2]:.5f} at "
       f"t={tsv[far]:.2f}",
       "both the middle and the far cell more than treble their t=0 share"))

    # R5e -- and what happens when the welds close a cycle
    hc = {}
    for nm, sites in (("4-cycle", [(0, 0, 0), (1, 1, 1), (2, 0, 0), (1, 1, -1)]),
                      ("brick 5x5x5", brick(5, 5, 5)),
                      ("brick 6x6x6", brick(6, 6, 6)),
                      ("ball r=3.5", ball(3.5)),
                      ("ball r=4.5", ball(4.5))):
        asm, deg = honeycomb(sites)
        cq = asm.q0()
        cc, cR, cg, cB = asm.frames(cq)
        rhc, gaphc = rank_of(asm.constraint_jacobian(
                asm.cell_jacobians(cc, cR, cB)))
        hc[nm] = (asm.N, len(asm.welds), 7 * asm.N - rhc - 6,
                  sum(1 for d in deg if d == 1), gaphc,
                  float(np.abs(asm.weld_residual(cq)).max()))
    trimmed = brick(5, 5, 5)
    for _ in range(5):
        asm, deg = honeycomb(trimmed)
        keep = [i for i, d in enumerate(deg) if d > 1]
        if len(keep) == len(trimmed):
            break
        trimmed = [trimmed[i] for i in keep]
    asm, deg = honeycomb(trimmed)
    tq = asm.q0()
    tc, tR, tg, tB = asm.frames(tq)
    rt, gapt = rank_of(asm.constraint_jacobian(asm.cell_jacobians(tc, tR, tB)))
    hc["5x5x5 trimmed"] = (asm.N, len(asm.welds), 7 * asm.N - rt - 6, 0, gapt,
                           float(np.abs(asm.weld_residual(tq)).max()))

    #: the 4-cycle's single surviving mode, read as a fold-rate pattern
    a4, _ = honeycomb([(0, 0, 0), (1, 1, 1), (2, 0, 0), (1, 1, -1)])
    q4 = a4.q0()
    c4, R4, g4, B4 = a4.frames(q4)
    M4 = np.vstack([a4.constraint_jacobian(a4.cell_jacobians(c4, R4, B4)),
                    a4.globals(c4)])
    r4, _ = rank_of(M4)
    fold4 = np.linalg.svd(M4)[2][r4:][0][[7 * k + 6 for k in range(4)]]
    spread = float(fold4.max() - fold4.min()) / float(np.abs(fold4).max())
    out["honeycomb"] = hc
    A(("R5e CLOSE THE CYCLES AND THE FIELD GOES AWAY. R5b's per-cell count is a "
       "TREE result -- a chain has no cycles and the 9-cluster is a star, so "
       "BOTH patches it was measured on are trees, and the honeycomb is not: "
       "brick 5x5x5 is 35 cells, 64 welds, 30 independent cycles. With cycles "
       "the count is 1 + (number of DEGREE-1 cells), which is why odd-sided "
       "bricks all report 9 whatever their size -- they carry eight dangling "
       "corner cells at every size, and that constancy is a property of the "
       "SHAPE, not of the lattice. A COMPACT patch has ONE internal degree of "
       "freedom, coherent breathing, at 59 cells and at 113 alike, and the "
       "minimal 4-cycle locks all four cells to the same fold rate. So a "
       "disturbance has nowhere to be in the interior after all -- which is "
       "R5b's orientation-FIXED answer, reached without freezing anything",
       all(v[2] == 1 + v[3] for v in hc.values())
       and hc["ball r=4.5"][2] == 1 and hc["ball r=3.5"][2] == 1
       and hc["5x5x5 trimmed"][2] == 1 and hc["4-cycle"][2] == 1
       and hc["brick 5x5x5"][2] == 9 and hc["brick 6x6x6"][2] == 3
       and spread < 1e-9
       and max(v[5] for v in hc.values()) < 1e-14,
       "; ".join(f"{k}: N={v[0]} welds={v[1]} dangling={v[3]} -> DOF {v[2]}"
                 for k, v in hc.items())
       + f"; the 4-cycle's one mode locks all four fold rates together, spread "
         f"{spread:.1e}; worst weld residual "
         f"{max(v[5] for v in hc.values()):.1e}; smallest rank gap "
         f"{min(v[4] for v in hc.values()):.1e}",
       "DOF = 1 + dangling cells; compact patches give 1 at any size"))

    # R5f -- and the sharing the face-welds never imposed
    sv = {}
    for nm, sites in (("9-cluster", brick(3, 3, 3)),
                      ("brick 5x5x5", brick(5, 5, 5)),
                      ("ball r=3.5", ball(3.5))):
        struct = True
        ref = None
        for ang in (-45.0, -30.0, -10.0, -52.3, -7.7):
            asm_a, _ = honeycomb(sites, ang)
            cls_a = shared_vertices(asm_a)
            if ref is None:
                ref = cls_a
            elif cls_a != ref:
                struct = False
        asm_s, _ = honeycomb(sites)
        cls = shared_vertices(asm_s)
        hist = {}
        for c in cls:
            hist[len(c)] = hist.get(len(c), 0) + 1
        weld_ids = 3 * len(asm_s.welds)
        span_ids = sum(len(c) - 1 for c in cls)
        Cf = identification_rows(asm_s, cls)
        rf2, gapf = rank_of(Cf)
        ctr_s = asm_s.frames(asm_s.q0())[0]
        Mf = np.vstack([Cf, asm_s.globals(ctr_s)])
        rM, _ = rank_of(Mf)
        null = np.linalg.svd(Mf)[2][rM:]
        fold = null[0][[7 * k + 6 for k in range(asm_s.N)]]
        sv[nm] = (asm_s.N, dict(sorted(hist.items())), struct, weld_ids, span_ids,
                  7 * asm_s.N - rf2 - 6,
                  float(fold.max() - fold.min()) / float(np.abs(fold).max()),
                  float(np.abs(np.dot(Cf, asm_s.globals(ctr_s).T)).max()), gapf)
    out["shared"] = sv
    A(("R5f THE HONEYCOMB SHARES VERTICES THE FACE-WELDS NEVER IMPOSE, and "
       "imposing them leaves ONE degree of freedom EVERYWHERE -- the tree-shaped "
       "9-cluster included, 9 -> 1. In the bulk a vertex is met by FOUR cells, "
       "while a face-weld identifies only the three corners of one triangle, so "
       "the 9-cluster leaves 12 coincidences un-imposed. This is not numerical "
       "luck: the class structure is BIT-IDENTICAL at five fold angles, so it is "
       "structural. The surviving mode is exactly uniform -- every cell folds at "
       "one rate, both sublattices together. CAVEAT, and it is the owner's call: "
       "a vertex contact is UNILATERAL in a real build, so imposing it as a "
       "bilateral joint is a modelling choice. The BULK answer does not depend "
       "on that choice -- R5e already gets 1 for a compact patch from welds "
       "alone -- only the dangling-corner modes at a cut boundary do",
       all(v[2] for v in sv.values())
       and all(v[5] == 1 for v in sv.values())
       and all(v[6] < 1e-9 for v in sv.values())
       and all(v[7] < 1e-12 for v in sv.values())
       and sv["9-cluster"][4] - sv["9-cluster"][3] == 12
       and max(sv["brick 5x5x5"][1]) == 4,
       "; ".join(f"{k}: N={v[0]} classes {v[1]} structural={v[2]} "
                 f"weld-imposed {v[3]} of {v[4]} spanning -> DOF {v[5]} "
                 f"(fold spread {v[6]:.0e})" for k, v in sv.items()),
       "structural at 5 angles, DOF 1 everywhere, one uniform fold rate"))

    # R5g -- and whether that sharing is a constraint at all
    asm_t, _ = honeycomb(brick(3, 3, 3))
    q0t = asm_t.q0()
    weldset = {frozenset([(k, a), (l, b)])
               for k, l, prs in asm_t.welds for a, b in prs}
    touch = [c for c in shared_vertices(asm_t)
             if len(c) == 2 and frozenset(c) not in weldset]
    ctr_t, R_t, g_t, B_t = asm_t.frames(q0t)
    Mt = np.vstack([asm_t.constraint_jacobian(
            asm_t.cell_jacobians(ctr_t, R_t, B_t)), asm_t.globals(ctr_t)])
    rt2, _ = rank_of(Mt)
    modes = np.linalg.svd(Mt)[2][rt2:]
    held, worst_w, smallest = 0, 0.0, np.inf
    for m in range(modes.shape[0]):
        for eps in (0.08, -0.08):
            qw = walk(asm_t, q0t, modes[m], eps)
            worst_w = max(worst_w, float(np.abs(asm_t.weld_residual(qw)).max()))
            Xw = asm_t.positions(qw)
            gaps = np.array([np.linalg.norm(Xw[a[0]][a[1]] - Xw[b[0]][b[1]])
                             for a, b in touch])
            held += int((gaps < 1e-9).sum())
            smallest = min(smallest, float(gaps.min()))
    out["contact"] = (len(touch), modes.shape[0], held, smallest, worst_w)
    A(("R5g THE VERTEX SHARING IS A TANGENCY, NOT A CONSTRAINT, so R5f's caveat "
       "resolves itself for the MEDIUM. The 12 identifications the face-welds "
       "miss are all between EVEN cells that share no face -- the square-face "
       "neighbours, which is ThreeCellAnimation's 'the outer cells touch each "
       "other without being asked to'. Walk a FINITE distance along every "
       "weld-only internal mode, in BOTH directions, and every one of those "
       "contacts OPENS: never once held closed, and the gap grows like eps^2, "
       "so the modes graze the contact rather than cross it. A real vertex "
       "contact is UNILATERAL -- two cells touching at a point can separate -- "
       "so it never binds and imposes nothing. TWO-SIDED: a mode that closed a "
       "contact would be counted here. What R5f measures is therefore the RIG's "
       "answer, where vertices are joined; the MEDIUM's is R5e's",
       len(touch) == 12 and out["contact"][2] == 0
       and out["contact"][3] > 1e-4 and worst_w < 1e-12,
       f"{len(touch)} even-even contacts, {modes.shape[0]} modes x 2 directions: "
       f"{held} contacts held closed, smallest gap opened "
       f"{out['contact'][3]:.1e} at |eps| = 0.08 (quadratic: eps^2 = 6.4e-3); "
       f"worst weld residual along the walks {worst_w:.1e}",
       "12 contacts, none ever held closed, all opening quadratically"))

    # R5h -- the medium's one motion, in closed form and bounded by its geometry
    bs = ball(2.0)
    basm, bdeg = honeycomb(bs, -30.0)
    bc, bR, bg, bB = basm.frames(basm.q0())
    brk, _ = rank_of(basm.constraint_jacobian(basm.cell_jacobians(bc, bR, bB)))
    bdof = 7 * basm.N - brk - 6
    bu, _ = fold_impulse(basm, 0)
    #: the closed form IS the constrained motion -- checked, not assumed
    worst_rot = worst_ctr = worst_weld = worst_mass = 0.0
    for ang in (-7.0, -20.0, -30.0, -45.0, -55.0):
        am, _ = honeycomb(bs, ang)
        cq, cu = coherent(bs, ang)
        ku, _ = fold_impulse(am, 0)
        ku = ku / ku[0, 6]                       # normalise to unit fold rate
        worst_rot = max(worst_rot, float(np.abs(ku[:, 3:6]).max()))
        worst_ctr = max(worst_ctr, float(np.abs(ku[:, 0:3] - cu[:, 0:3]).max()))
        worst_weld = max(worst_weld, float(np.abs(am.weld_residual(cq)).max()))
        worst_mass = max(worst_mass, abs(effective_mass(bs, ang)
                                         - 2.0 * float(am.energy(cq, cu).sum())))
    hfd = 1e-5
    dm_err = 0.0
    for ang in (-7.0, -30.0, -52.0):
        _, dm = effective_mass(bs, ang, derivative=True)
        fd = (effective_mass(bs, ang + np.degrees(hfd))
              - effective_mass(bs, ang - np.degrees(hfd))) / (2 * hfd)
        dm_err = max(dm_err, abs(dm - fd) / abs(fd))
    srec, bounces = swing(bs, -30.0, 0.30, 40.0, 200)
    sE = np.array([r[3] for r in srec])
    sa = np.array([r[1] for r in srec])
    lam = [lattice_constant(x) for x in (0.0, -30.0, -60.0)]
    out["breathe"] = (basm.N, len(basm.welds), bdof, bounces, srec, bs)
    A(("R5h THE MEDIUM'S ONE MOTION, IN CLOSED FORM, BOUNDED BY ITS OWN "
       "GEOMETRY. A compact patch -- 15 cells, no dangling ones, ONE internal "
       "degree of freedom -- answers a phase kick on a SINGLE cell with a "
       "perfectly coherent breathe, because there is no other motion available "
       "to it. And that motion is closed form: every cell folds at the same "
       "rate, NO cell rotates, every centre just rides the lattice constant "
       "L(a). Checked against the constrained model rather than assumed, which "
       "is what lets the animation skip the constraint solve entirely. Past "
       "|gamma| = 60 a cell folds through its own octahedron and self-"
       "intersects, so with b = a + 60 forcing both sublattices into the window "
       "at once the range is a in [-60, 0] -- the limits are the geometry's, "
       "not a convention. The impact is VELOCITY REVERSAL, exact rather than "
       "convenient: with one degree of freedom every admissible velocity is a "
       "multiple of one mode, so energy survives the bounce identically",
       bdof == 1 and sum(1 for d in bdeg if d == 1) == 0
       and float(bu[:, 6].max() - bu[:, 6].min()) < 1e-12
       and worst_rot < 1e-11 and worst_ctr < 1e-10
       and worst_weld < 1e-14 and worst_mass < 1e-9 and dm_err < 1e-6
       and bounces >= 6
       and sa.min() > -FOLD_LIMIT and sa.max() < 0.0
       and float(np.abs(sE - sE[0]).max() / sE[0]) < 1e-9,
       f"N={basm.N} welds={len(basm.welds)} dangling=0 -> DOF {bdof}; one-cell "
       f"kick gives fold rates equal to "
       f"{float(bu[:, 6].max() - bu[:, 6].min()):.0e}; over five angles the "
       f"closed form has cell rotation {worst_rot:.0e}, centre velocity "
       f"matching dL/da*site to {worst_ctr:.0e}, welds closed to "
       f"{worst_weld:.0e}, and M_eff agreeing with the assembly's own energy to "
       f"{worst_mass:.0e}; analytic dM/da matches central differences to "
       f"{dm_err:.0e}; {bounces} elastic bounces over t=0..40 with a held in "
       f"[{sa.min():.2f}, {sa.max():.2f}] and E = {sE[0]:.9f} conserved to "
       f"{float(np.abs(sE - sE[0]).max() / sE[0]):.0e} ACROSS them; M_eff "
       f"{effective_mass(bs, 0.0):.1f} -> {effective_mass(bs, -30.0):.1f} -> "
       f"{effective_mass(bs, -60.0):.1f}, so the breathe is fastest at "
       f"mid-swing; lattice constant ratio {lam[1] / lam[0]:.9f} against "
       f"2/sqrt(3) = {2 / np.sqrt(3):.9f}",
       "1 DOF, closed form verified, bounded in [-60,0], E exact across bounces"))

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

    print("\n  THE HONEYCOMB, WITH ITS CYCLES -- read this before the tree table:")
    print(f"    {'patch':<18}{'cells':>7}{'welds':>7}{'dangling':>10}{'internal DOF':>14}")
    for k, v in out["honeycomb"].items():
        print(f"    {k:<18}{v[0]:>7}{v[1]:>7}{v[3]:>10}{v[2]:>14}")
    print("    A COMPACT patch has ONE internal degree of freedom at any size.")
    print("    The 9 that odd-sided bricks report is eight dangling corner cells.")
    print("\n  ...and imposing the vertex sharing the face-welds never did:")
    for k, v in out["shared"].items():
        print(f"    {k:<14} N={v[0]:<4} vertex classes {str(v[1]):<24} "
              f"-> DOF {v[5]}  (one uniform fold rate)")
    print("    ONE degree of freedom EVERYWHERE, the tree-shaped cluster included.")

    gdc, gdo, ev, mp, mt, fp, ft = out["vee"]
    print("\n  THREE-CELL V -- the configuration ThreeCellAnimation builds, and")
    print("  the reference for the Java port in "
          "com.chiralbehaviors.inviscid.jitterbug:")
    print(f"    centre-kick fold rates {gdc:.9f} = 57/137 and {gdo:.9f} = 12/137,")
    print(f"    total energy {ev:.9f} = 152/137")
    print(f"    outer kick: middle peaks {mp:.5f} at t={mt:.2f}, "
          f"far peaks {fp:.5f} at t={ft:.2f}")

    bn, bw, bdf, bb, brec, bsites = out["breathe"]
    print(f"\n  A COMPACT PATCH BREATHING -- {bn} cells, {bw} welds, {bdf} internal")
    print("  degree of freedom, bounced elastically off both fold limits:")
    print(f"  {'t':>7} {'a':>9} {'a+60':>8} {'adot':>9} {'L(a)':>9} {'E':>13} {'bounces':>8}")
    for t, a, ad, e, nb in brec[::max(1, len(brec) // 12)]:
        print(f"  {t:7.2f} {a:9.3f} {a + 60:8.3f} {ad:+9.4f} "
              f"{lattice_constant(a):9.6f} {e:13.9f} {nb:8d}")

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
    print("     everywhere (R10), so there is no wavefront to time -- and")
    print("     that is the PLAY-FREE LIMIT, not a property of the medium.")
    print("     CORRECTED 2026-08-28: this line used to say a finite signal")
    print("     speed would need compliant constraints, which is FALSE. It")
    print("     needs CLEARANCE, which perfectly rigid struts can have.")
    print("     jb_ct measures the finite speed and jb_pr the bounded phase")
    print("     gradient that goes with it; both diverge to this file's")
    print("     answer as the play goes to zero.")
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
