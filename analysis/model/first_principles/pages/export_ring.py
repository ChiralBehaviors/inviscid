"""Frames of the four-body ring along the exchange, from the model's
honeycomb builder at each fold: two VEs (O, D) and the two voids (U, L)
between them, carried as jitterbugs whose fold is a + 60. a in [-60, 0]
is the ring's range; +2.5..+10 are exported only to show the void folded
past its own octahedron (its struts crossing)."""
import itertools as it, json
import numpy as np

from analysis.model.first_principles.pages import common
from analysis.model import assembly as RC
from analysis.model import cell as IC

SITES = [(0, 0, 0), (1, 1, 1), (2, 2, 0), (1, 1, -1)]
NAMES = ["O", "U", "D", "L"]
TRI = [[IC.SLOT[(f, c)] for c in range(3)] for f in range(8)]

def segdist(p1, q1, p2, q2):
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    b, c = d1 @ d2, d1 @ r
    den = a * e - b * b
    s = np.clip((b * f - c * e) / den, 0, 1) if den > 1e-14 else 0.0
    t = (b * s + f) / e
    if t < 0: t = 0.0; s = np.clip(-c / a, 0, 1)
    elif t > 1: t = 1.0; s = np.clip((b - c) / a, 0, 1)
    return float(np.linalg.norm((p1 + d1 * s) - (p2 + d2 * t))), float(s), float(t)

def self_cross(X):  # X (12,3) vertex labels of one cell: nearest struts of different triangles not sharing a joint
    S = [(f, X[TRI[f][i]], X[TRI[f][j]]) for f in range(8) for i, j in ((0,1),(1,2),(2,0))]
    best = (np.inf, None)
    for (f1, p1, q1), (f2, p2, q2) in it.combinations(S, 2):
        if f1 == f2: continue
        if min(np.linalg.norm(x - y) for x in (p1, q1) for y in (p2, q2)) < 1e-9: continue
        d, s, t = segdist(p1, q1, p2, q2)
        if d < best[0]: best = (d, (f1, f2, s, t))
    d, (f1, f2, s, t) = best
    return d, bool(d < 1e-9 and 1e-6 < s < 1-1e-6 and 1e-6 < t < 1-1e-6), [f1, f2]

angles = list(np.arange(-60.0, 0.0 + 1e-9, 2.5)) + [2.5, 5.0, 7.5, 10.0]
frames = []
for a in angles:
    asm, _ = RC.honeycomb(SITES, gc=a)
    q = asm.q0(); X = asm.positions(q)
    res = float(np.abs(asm.weld_residual(q)).max())
    gam = list(asm.gam0)
    cells = [[[np.round(X[k][i], 6).tolist() for i in TRI[f]] for f in range(8)] for k in range(4)]
    # corner census at P = O's vertex shared with U and L (the vertex at (1,1,0)*L)
    L = RC.lattice_constant(a)
    dU, crossU, pairU = self_cross(X[1])
    frames.append({"a": float(a), "b": float(gam[1]), "L": float(L), "res": res, "cells": cells,
                   "void_clear": dU, "void_cross": crossU, "void_pair": pairU})
    if a in (-60.0, -30.0, 0.0, 5.0, 10.0):
        print(f"a={a:+6.1f}: void fold {gam[1]:+6.1f}, spacing {L:.4f}, weld {res:.0e}, void self-clearance {dU:.4f} crossing {crossU}")
json.dump({"names": NAMES, "sites": SITES, "strut": IC.EL, "frames": frames, "ive": angles.index(0.0)},
          open(str(common.out("ring")), "w"))
print("frames", len(frames))
