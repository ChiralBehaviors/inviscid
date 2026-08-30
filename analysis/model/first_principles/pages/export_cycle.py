"""The four-body ring at every whole degree of its range, from the model's
builder, for a smooth full cycle VE -> octahedron -> VE."""
import json
import numpy as np

from analysis.model.first_principles.pages import common
from analysis.model import assembly as RC
from analysis.model import cell as IC
SITES = [(0, 0, 0), (1, 1, 1), (2, 2, 0), (1, 1, -1)]
TRI = [[IC.SLOT[(f, c)] for c in range(3)] for f in range(8)]
frames = []
for a in np.arange(-60.0, 0.0 + 1e-9, 1.0):
    asm, _ = RC.honeycomb(SITES, gc=float(a))
    q = asm.q0(); X = asm.positions(q)
    res = float(np.abs(asm.weld_residual(q)).max())
    frames.append({"a": float(a), "b": float(asm.gam0[1]), "L": float(RC.lattice_constant(a)), "res": res,
                   "cells": [[[np.round(X[k][i], 5).tolist() for i in TRI[f]] for f in range(8)] for k in range(4)]})
print("frames", len(frames), "max weld residual", max(f["res"] for f in frames))
json.dump({"frames": frames, "strut": IC.EL}, open(str(common.out("cycle")), "w"))
