"""Kinematics of the reduced-coordinate assembly: tied vertex pairs, mass
matrices under the declared mass model, the free (gyroscopic) acceleration,
joint separations and band rows, kinetic energy, and the standard kicks.

Split out of the hard-wall impact-law module (once jb_mj) on 2026-08-30 so
that nothing in the model of record imports a retired file. The functions are
unchanged; their docstrings still say what they were written for.
"""
from __future__ import annotations

import itertools as it

import numpy as np

from analysis.model import assembly as RC

NV = RC.NV
SLOT = RC.SLOT

#: The reference phase. Midpoint of the exchange, widest lattice, and the one
#: phase at which the front is uniform (jb_ct R4).
A_REF = -30.0

#: Joint clearance, in the same units as the edge (EL = sqrt(2)). EXAGGERATED
#: relative to a real build's ~1% for the same reason jb_ct's is: at 1e-2 the
#: contacts are too rare to exercise the LCP in a short run. Every scaling in
#: jb_ct is 1/play, so a larger band only slows things down.
PLAY = 0.05

#: Lamina spin coefficient. 1/12 for a uniform triangular plate against 1/3 for
#: three corner point masses -- the 9:1 / 3:1 second-moment split C2 records.
K_LAMINA = 1.0 / 12.0


def hc15_sites():
    """One VE and all fourteen of its face neighbours: 8 triangular at
    (+-1,+-1,+-1) and 6 square at (+-2,0,0). 6 + 8 = the VE's whole
    coordination, nothing left over (jb_hc's census)."""
    s = [(0, 0, 0)]
    s += [tuple(t) for t in it.product((1, -1), repeat=3)]
    for d in range(3):
        for sg in (1, -1):
            s.append(tuple(int(2 * sg * v) for v in np.eye(3, dtype=int)[d]))
    return [tuple(int(x) for x in t) for t in s]


def build(gc=A_REF):
    asm, deg = RC.honeycomb(hc15_sites(), gc=gc)
    return asm


def tied_pairs(asm):
    """One UNILATERAL band per tied vertex pair -- DECISION 18's wire, which on
    this packing is three per shared face rather than one per neighbour."""
    return [(k, a, l, b) for (k, l, ps) in asm.welds for (a, b) in ps]


def kinematics(asm, q, lamina=False):
    """(J, M, Minv) at `q`. J maps a cell's seven coordinates to its vertex
    velocities; M is the mass matrix under the declared model."""
    ctr, R, gam, B = asm.frames(q)
    J = asm.cell_jacobians(ctr, R, B)
    if not lamina:
        M = asm.mass_blocks(J)
    else:
        M = np.zeros((asm.N, 7, 7))
        for k in range(asm.N):
            for f in range(8):
                rows = [J[k][3 * SLOT[(f, c)]:3 * SLOT[(f, c)] + 3]
                        for c in range(3)]
                cen = (rows[0] + rows[1] + rows[2]) / 3.0
                for r in rows:
                    M[k] += K_LAMINA * np.dot(r.T, r)
                M[k] += (1.0 - 3.0 * K_LAMINA) * np.dot(cen.T, cen)
    Minv = np.array([np.linalg.inv(M[k]) for k in range(asm.N)])
    return J, M, Minv


def free_accel(asm, q, u, J, Minv, lamina=False):
    """V = 0 and the bands are inactive between impacts, so the only term is
    the gyroscopic one. This is `Assembly.accel`'s `a_free`, without the
    bilateral weld reaction that file adds.

    THE MASS MODEL ENTERS HERE TOO. `Assembly.accel` weights the gyroscopic
    term with `VMASS`, the point-mass distribution; under laminae that would
    pair a lamina METRIC with a point-mass FORCE and the run would conserve
    nothing. The centroid's gyroscopic acceleration is the mean of its three
    corners' by the same linearity that makes its velocity the mean of theirs.
    """
    ctr, R, gam, B = asm.frames(q)
    A = asm.cell_gyro(R, B, u)
    if not lamina:
        m3 = np.repeat(RC.VMASS, 3)
        f = np.array([-np.dot(J[k].T, m3 * A[k].ravel())
                      for k in range(asm.N)])
    else:
        f = np.zeros((asm.N, 7))
        for k in range(asm.N):
            for fc in range(8):
                idx = [SLOT[(fc, c)] for c in range(3)]
                rows = [J[k][3 * i:3 * i + 3] for i in idx]
                accs = [A[k][i] for i in idx]
                cen_r = (rows[0] + rows[1] + rows[2]) / 3.0
                cen_a = (accs[0] + accs[1] + accs[2]) / 3.0
                for r, ac in zip(rows, accs):
                    f[k] -= K_LAMINA * np.dot(r.T, ac)
                f[k] -= (1.0 - 3.0 * K_LAMINA) * np.dot(cen_r.T, cen_a)
    return np.einsum('kij,kj->ki', Minv, f)


def separations(asm, q, pairs):
    X = asm.positions(q)
    return np.array([float(np.linalg.norm(X[k][a] - X[l][b]))
                     for (k, a, l, b) in pairs])


def band_rows(asm, q, J, pairs):
    """N: one row per tied pair, giving d/dt ||va - vb||.

    The band is ||va - vb|| <= t, so a row is the unit separation direction
    contracted with the difference of the two vertices' Jacobians. Where the
    pair is exactly coincident the direction is undefined and the row is zero --
    which is correct: a joint sitting dead centre in its clearance is not
    approaching a stop in any direction.
    """
    X = asm.positions(q)
    N = np.zeros((len(pairs), 7 * asm.N))
    for r, (k, a, l, b) in enumerate(pairs):
        d = X[k][a] - X[l][b]
        n = float(np.linalg.norm(d))
        if n < 1e-14:
            continue
        u = d / n
        N[r, 7 * k:7 * k + 7] = np.dot(u, J[k][3 * a:3 * a + 3])
        N[r, 7 * l:7 * l + 7] -= np.dot(u, J[l][3 * b:3 * b + 3])
    return N


def kinetic(M, u):
    """Kinetic energy under the metric supplied, which is the only correct way
    to ask for it here -- see the note in `run`."""
    return float(np.sum([0.5 * u[k] @ M[k] @ u[k] for k in range(len(u))]))


def staggered_kick(asm, amp=0.6):
    """The two sublattices folding in OPPOSITE senses -- the theta = pi plane
    wave on this lattice, and the initial condition that loads every joint at
    once. VE cells sit at all-even sites and hole cells at all-odd."""
    even = np.array([all(c % 2 == 0 for c in s) for s in hc15_sites()])
    u = np.zeros((asm.N, 7))
    u[:, 6] = np.where(even, amp, -amp)
    return u


def coherent_kick(asm, amp=0.6):
    """Every cell's fold rate driven together -- the medium's one motion, and
    the initial condition that loads every joint at once."""
    u = np.zeros((asm.N, 7))
    u[:, 6] = amp
    return u


def single_kick(asm, amp=0.6, cell=0):
    u = np.zeros((asm.N, 7))
    u[cell, 6] = amp
    return u
