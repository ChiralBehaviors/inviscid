"""Corrected export: the model's own lattice at every frame (cells at a,
spacing lattice_constant(a)); ties are identities and never change. Range
-60..0 is the tied array's; +2.5..+10 is exported ONLY to show the strut
crossing that stops it there. Per frame: coincidence groups, the gap between
the untied corner pair at this vertex, and the nearest foreign strut pair
among the four drawn cells (highlighted when it is a crossing)."""
import itertools as it, json
import numpy as np

from analysis.model.first_principles.pages import common
from analysis.model import assembly as RC
from analysis.model import plates as Z
from analysis.model import cell as IC

SITES4 = [(0, 0, 0), (2, 0, 0), (0, 2, 0), (2, 2, 0)]
NAMES = ["O", "A", "B", "D"]
L0 = RC.lattice_constant(0.0)
VTX = L0 * np.array([1.0, 1.0, 0.0])

def segpts(p1, q1, p2, q2):
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    b, c = d1 @ d2, d1 @ r
    den = a * e - b * b
    s = np.clip((b * f - c * e) / den, 0, 1) if den > 1e-14 else 0.0
    t = (b * s + f) / e
    if t < 0: t = 0.0; s = np.clip(-c / a, 0, 1)
    elif t > 1: t = 1.0; s = np.clip((b - c) / a, 0, 1)
    return float(np.linalg.norm((p1 + d1 * s) - (p2 + d2 * t))), float(s), float(t)

C0 = Z.corners(0.0)
at = [(k, f, c) for k, s in enumerate(SITES4) for f in range(8) for c in range(3)
      if np.linalg.norm(C0[f][c] + L0 * np.array(s, float) - VTX) < 1e-9]
assert len(at) == 8
angles = list(np.arange(-60.0, 0.0 + 1e-9, 2.5)) + [2.5, 5.0, 7.5, 10.0]
frames = []
for a in angles:
    L = RC.lattice_constant(a); C = Z.corners(a)
    cells = [C + L * np.array(s, float) for s in SITES4]
    pts = [cells[k][f][c] for (k, f, c) in at]
    cls = []
    for i, p in enumerate(pts):
        for g in cls:
            if np.linalg.norm(pts[g[0]] - p) < 1e-9:
                g.append(i); break
        else:
            cls.append([i])
    groups = [sorted({NAMES[at[i][0]] for i in g}) for g in cls]
    # gap between the untied corners at this vertex: O's vs B's
    iO = next(i for i, (k, f, c) in enumerate(at) if k == 0)
    iB = next(i for i, (k, f, c) in enumerate(at) if k == 2)
    gap = float(np.linalg.norm(pts[iO] - pts[iB]))
    # nearest foreign strut pair among the drawn cells (no shared endpoint)
    best = (np.inf, None)
    for k1, k2 in it.combinations(range(4), 2):
        for f1 in range(8):
            for i1, j1 in ((0, 1), (1, 2), (2, 0)):
                p1, q1 = cells[k1][f1][i1], cells[k1][f1][j1]
                for f2 in range(8):
                    for i2, j2 in ((0, 1), (1, 2), (2, 0)):
                        p2, q2 = cells[k2][f2][i2], cells[k2][f2][j2]
                        if min(np.linalg.norm(x - y) for x in (p1, q1) for y in (p2, q2)) < 1e-9:
                            continue
                        d, sp, tp = segpts(p1, q1, p2, q2)
                        if d < best[0]:
                            best = (d, [k1, f1, i1, j1, k2, f2, i2, j2, sp, tp])
    d, pair = best
    crossing = d < 1e-9 and 1e-6 < pair[8] < 1 - 1e-6 and 1e-6 < pair[9] < 1 - 1e-6
    frames.append({"a": float(a), "L": float(L), "cells": [np.round(P, 6).tolist() for P in cells],
                   "groups": groups, "npts": len(cls), "gap": gap,
                   "near": float(d), "pair": pair[:8], "crossing": bool(crossing)})
    if a in (-60.0, -30.0, -5.0, 0.0, 2.5, 10.0):
        print(f"a={a:+6.1f} L={L:.4f} points {len(cls)} {groups} gap O-B {gap:.4f} nearest {d:.4f} crossing {crossing}")
voids = []
for s in ((1, 1, 1), (1, 1, -1)):
    c = L0 * np.array(s, float)
    voids.append(np.round([c + L0 * e for e in np.eye(3)] + [c - L0 * e for e in np.eye(3)], 6).tolist())
out = {"strut": IC.EL, "L0": L0, "vertex": VTX.tolist(), "at": at, "names": NAMES,
       "sites": SITES4, "voids": voids, "frames": frames, "ive": angles.index(0.0)}
json.dump(out, open(str(common.out("vertex")), "w"))
print("frames", len(frames), "VE index", angles.index(0.0))
