"""Fifteen bodies -- one VE, its eight voids, its six axis VEs (the
model's hc15 census) -- driven through 360 deg of the fold. Joints held by
identity from a = -30; residual on every frame; fronts out per body;
strut crossings among all 360 struts (vectorised)."""
import itertools as it, json
import numpy as np
from analysis.model import assembly as RC
from analysis.model import cell as IC
from analysis.model import kinematics as MJ
from analysis.model.first_principles.pages import common

import sys
PATCH = sys.argv[1] if len(sys.argv) > 1 else "ring"
PATCHES = {"ring": [(0, 0, 0), (1, 1, 1), (2, 2, 0), (1, 1, -1)],
           "hc15": [tuple(int(c) for c in s) for s in MJ.hc15_sites()],
           "block": [tuple(int(c) for c in s) for s in RC.brick(5, 5, 5)]}
SITES = PATCHES[PATCH]
STEP = {"ring": 2.0, "hc15": 3.0, "block": 3.0}[PATCH]
NB = len(SITES)
KIND = ["cell" if all(c % 2 == 0 for c in s) else "void" for s in SITES]
TRI = [[IC.SLOT[(f, c)] for c in range(3)] for f in range(8)]
ref, deg = RC.honeycomb(SITES, gc=-30.0)
WELDS = ref.welds
print(f"bodies {NB} ({KIND.count('cell')} cells, {KIND.count('void')} voids), face welds {len(WELDS)}, degrees {deg}")

def crossings(P, Q, owner):
    """P,Q (n,3) strut endpoints; owner (n,) body id. Count crossing pairs:
    different triangles, no shared endpoint, distance < 1e-6 at interior params."""
    n = len(P); i, j = np.triu_indices(n, 1)
    p1, q1, p2, q2 = P[i], Q[i], P[j], Q[j]
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a = np.einsum('ij,ij->i', d1, d1); e = np.einsum('ij,ij->i', d2, d2); f = np.einsum('ij,ij->i', d2, r)
    b = np.einsum('ij,ij->i', d1, d2); c = np.einsum('ij,ij->i', d1, r)
    den = a * e - b * b
    s = np.where(den > 1e-14, np.clip((b * f - c * e) / np.where(den > 1e-14, den, 1), 0, 1), 0.0)
    t = (b * s + f) / e
    lo, hi = t < 0, t > 1
    t = np.clip(t, 0, 1)
    s = np.where(lo, np.clip(-c / a, 0, 1), s); s = np.where(hi, np.clip((b - c) / a, 0, 1), s)
    dist = np.linalg.norm((p1 + d1 * s[:, None]) - (p2 + d2 * t[:, None]), axis=1)
    ends = np.min(np.stack([np.linalg.norm(p1 - p2, axis=1), np.linalg.norm(p1 - q2, axis=1),
                            np.linalg.norm(q1 - p2, axis=1), np.linalg.norm(q1 - q2, axis=1)]), axis=0)
    same_tri = owner[i] == owner[j]
    ok = (~same_tri) & (ends > 1e-6) & (dist < 1e-6) & (s > 1e-4) & (s < 1 - 1e-4) & (t > 1e-4) & (t < 1 - 1e-4)
    return int(ok.sum())

asm0, _ = RC.honeycomb(SITES, gc=0.0); X0 = asm0.positions(asm0.q0()); C0 = asm0.ctr0
SIGN = np.zeros((NB, 8))
for k in range(NB):
    for f in range(8):
        Pt = X0[k][TRI[f]]; n = np.cross(Pt[1]-Pt[0], Pt[2]-Pt[0]); SIGN[k, f] = np.sign(n @ (Pt.mean(0) - C0[k]))
frames = []
for a in np.arange(-60.0, 300.0 + 1e-9, STEP):
    asm, _ = RC.honeycomb(SITES, gc=float(a))
    held = RC.Assembly(asm.gam0, asm.ctr0, WELDS)
    q = held.q0(); X = held.positions(q)
    res = float(np.abs(held.weld_residual(q)).max())
    cells, outs = [], []
    P, Q, owner = [], [], []
    for k in range(NB):
        tris = []; n_out = 0
        for f in range(8):
            Pt = X[k][TRI[f]]
            n = np.cross(Pt[1]-Pt[0], Pt[2]-Pt[0]) * SIGN[k, f]
            out = bool(n @ (Pt.mean(0) - asm.ctr0[k]) > 0); n_out += out
            tris.append({"p": np.round(Pt, 4).tolist(), "s": int(SIGN[k, f]), "out": out})
            for i2, j2 in ((0, 1), (1, 2), (2, 0)):
                P.append(Pt[i2]); Q.append(Pt[j2]); owner.append(k * 8 + f)
        cells.append(tris); outs.append(n_out)
    ncross = crossings(np.array(P), np.array(Q), np.array(owner))
    ctr = asm.ctr0
    spread = float(np.max(np.linalg.norm(ctr - ctr.mean(0), axis=1)))
    frames.append({"a": float(a), "b": float(a + 60.0), "L": float(RC.lattice_constant(a)), "res": res,
                   "cells": cells, "out": outs, "out_cells": int(sum(o for o, k in zip(outs, KIND) if k == "cell")),
                   "out_voids": int(sum(o for o, k in zip(outs, KIND) if k == "void")),
                   "cross": ncross, "ctr": np.round(ctr, 4).tolist(), "spread": spread})
    if a % 30 == 0:
        print(f"a={a:+6.0f}  L={RC.lattice_constant(a):+.3f}  weld {res:.0e}  fronts out cells {frames[-1]["out_cells"]}/{8*KIND.count("cell")} voids {frames[-1]["out_voids"]}/{8*KIND.count("void")}  crossings {ncross:5d}  centres spread {spread:.3f}")
json.dump({"frames": frames, "kind": KIND, "strut": IC.EL, "nb": NB, "patch": PATCH},
          open(str(common.out(f"overdrive_{PATCH}")), "w"))
print("frames", len(frames))
