"""geometry -- the few measuring tools the first-principles modules share.

Coincidence classes of triangle corners, the twelve joints of a cell read at
a generic angle, segment-segment distance (scalar and vectorised), a rigid
fit, and the corner triples of a cell's eight plates in the assembly's
vertex-label space. Nothing here models anything; it only measures.
"""
from __future__ import annotations

import itertools as it

import numpy as np

from analysis.model import cell as IC
from analysis.model import plates as Z

#: corner triples of the eight plates, as assembly vertex labels
TRI = [[IC.SLOT[(f, c)] for c in range(3)] for f in range(8)]


def classes(C, tol=1e-9):
    """Group the 24 corners of one cell (8, 3, 3) into coincident classes."""
    cls = []
    for fc in [(f, c) for f in range(8) for c in range(3)]:
        p = C[fc[0]][fc[1]]
        for g in cls:
            if np.linalg.norm(C[g[0][0]][g[0][1]] - p) < tol:
                g.append(fc)
                break
        else:
            cls.append([fc])
    return cls


def joints(a_ref=-30.0):
    """The twelve joints of one cell as corner pairs, read at a generic angle."""
    cls = classes(Z.corners(a_ref))
    assert len(cls) == 12 and all(len(g) == 2 for g in cls)
    return [tuple(sorted(g)) for g in cls]


def joint_gap(C, J):
    """Largest separation of any joint's two corners."""
    return max(float(np.linalg.norm(C[f1][c1] - C[f2][c2])) for ((f1, c1), (f2, c2)) in J)


def segdist(p1, q1, p2, q2):
    """Minimum distance between segments p1q1 and p2q2, with the parameters."""
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    b, c = d1 @ d2, d1 @ r
    den = a * e - b * b
    s = np.clip((b * f - c * e) / den, 0, 1) if den > 1e-14 else 0.0
    t = (b * s + f) / e
    if t < 0:
        t = 0.0
        s = np.clip(-c / a, 0, 1)
    elif t > 1:
        t = 1.0
        s = np.clip((b - c) / a, 0, 1)
    return float(np.linalg.norm((p1 + d1 * s) - (p2 + d2 * t))), float(s), float(t)


def cell_struts(C):
    """The 24 struts of one cell as (face, p, q)."""
    return [(f, C[f][i], C[f][j]) for f in range(8) for i, j in ((0, 1), (1, 2), (2, 0))]


def self_clearance(C):
    """Nearest two struts of one cell that belong to different plates and share no joint."""
    S = cell_struts(C)
    best = np.inf
    for (f1, p1, q1), (f2, p2, q2) in it.combinations(S, 2):
        if f1 == f2:
            continue
        if min(np.linalg.norm(x - y) for x in (p1, q1) for y in (p2, q2)) < 1e-9:
            continue
        best = min(best, segdist(p1, q1, p2, q2)[0])
    return float(best)


def crossings(P, Q, owner, tol=1e-6):
    """Count crossing strut pairs among segments P[i]Q[i] (vectorised).

    A crossing is two struts of DIFFERENT plates (owner differs), sharing no
    endpoint, at distance < tol with both parameters strictly interior.
    """
    P, Q, owner = np.asarray(P), np.asarray(Q), np.asarray(owner)
    n = len(P)
    i, j = np.triu_indices(n, 1)
    p1, q1, p2, q2 = P[i], Q[i], P[j], Q[j]
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a = np.einsum("ij,ij->i", d1, d1)
    e = np.einsum("ij,ij->i", d2, d2)
    f = np.einsum("ij,ij->i", d2, r)
    b = np.einsum("ij,ij->i", d1, d2)
    c = np.einsum("ij,ij->i", d1, r)
    den = a * e - b * b
    safe = np.where(den > 1e-14, den, 1.0)
    s = np.where(den > 1e-14, np.clip((b * f - c * e) / safe, 0, 1), 0.0)
    t = (b * s + f) / e
    lo, hi = t < 0, t > 1
    t = np.clip(t, 0, 1)
    s = np.where(lo, np.clip(-c / a, 0, 1), s)
    s = np.where(hi, np.clip((b - c) / a, 0, 1), s)
    dist = np.linalg.norm((p1 + d1 * s[:, None]) - (p2 + d2 * t[:, None]), axis=1)
    ends = np.min(np.stack([np.linalg.norm(p1 - p2, axis=1), np.linalg.norm(p1 - q2, axis=1),
                            np.linalg.norm(q1 - p2, axis=1), np.linalg.norm(q1 - q2, axis=1)]), axis=0)
    ok = ((owner[i] != owner[j]) & (ends > tol) & (dist < tol)
          & (s > 1e-4) & (s < 1 - 1e-4) & (t > 1e-4) & (t < 1 - 1e-4))
    return int(ok.sum())


def assembly_struts(X):
    """All struts of an assembly's cells from its (N, 12, 3) vertex positions."""
    P, Q, owner = [], [], []
    for k in range(X.shape[0]):
        for f in range(8):
            for i, j in ((0, 1), (1, 2), (2, 0)):
                P.append(X[k][TRI[f][i]])
                Q.append(X[k][TRI[f][j]])
                owner.append(k * 8 + f)
    return np.array(P), np.array(Q), np.array(owner)


def kabsch(P, Q):
    """Proper rotation R and translation t with R P + t = Q for matched rows."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def rot_angle_about(R, axis):
    """Signed rotation angle of R about a unit axis, in degrees."""
    ax = np.asarray(axis, float)
    return float(np.degrees(np.arctan2(
        (R[2, 1] - R[1, 2]) * ax[0] + (R[0, 2] - R[2, 0]) * ax[1] + (R[1, 0] - R[0, 1]) * ax[2],
        np.trace(R) - 1.0)))


def front_signs(X0, ctr0):
    """Per plate, the sign that makes (v1-v0)x(v2-v0) point OUT of its body at the reference."""
    N = X0.shape[0]
    S = np.zeros((N, 8))
    for k in range(N):
        for f in range(8):
            P = X0[k][TRI[f]]
            n = np.cross(P[1] - P[0], P[2] - P[0])
            S[k, f] = np.sign(n @ (P.mean(0) - ctr0[k]))
    return S


def fronts_out(X, ctr, S):
    """Per body, how many of its eight plates have their front pointing out of it."""
    N = X.shape[0]
    out = np.zeros(N, int)
    for k in range(N):
        for f in range(8):
            P = X[k][TRI[f]]
            n = np.cross(P[1] - P[0], P[2] - P[0]) * S[k, f]
            out[k] += int(n @ (P.mean(0) - ctr[k]) > 0)
    return out
