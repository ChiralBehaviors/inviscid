"""One vertex point. The O+A joint and the B+D joint that meet there at the
VE; where each goes as the fold runs -15..+15 (model spacing throughout);
the eight triangles at the point; the crossing strut pair on the far side."""
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
C0 = Z.corners(0.0)
at = [(k, f, c) for k, s in enumerate(SITES4) for f in range(8) for c in range(3)
      if np.linalg.norm(C0[f][c] + L0 * np.array(s, float) - VTX) < 1e-9]
tris = sorted({(k, f) for k, f, c in at})           # the 8 triangles at the point
assert len(at) == 8 and len(tris) == 8

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

frames = []
for a in np.arange(-15.0, 15.0 + 1e-9, 0.5):
    L = RC.lattice_constant(a); C = Z.corners(a)
    cells = [C + L * np.array(s, float) for s in SITES4]
    P = {k: cells[k][f][c] for (k, f, c) in at}      # one corner per cell at this point
    jOA = P[0]; jBD = P[2]
    assert np.linalg.norm(P[0] - P[1]) < 1e-9 and np.linalg.norm(P[2] - P[3]) < 1e-9
    d = jBD - jOA
    # crossing among the eight triangles' struts (no shared endpoint)
    best = (np.inf, None)
    S = [((k, f, i, j), cells[k][f][i], cells[k][f][j]) for (k, f) in tris for i, j in ((0,1),(1,2),(2,0))]
    for (l1, p1, q1), (l2, p2, q2) in it.combinations(S, 2):
        if l1[0] == l2[0]: continue
        if min(np.linalg.norm(x - y) for x in (p1, q1) for y in (p2, q2)) < 1e-9: continue
        dd, sp, tp = segpts(p1, q1, p2, q2)
        if dd < best[0]: best = (dd, (l1, l2, sp, tp))
    dd, (l1, l2, sp, tp) = best
    crossing = dd < 1e-9 and 1e-6 < sp < 1 - 1e-6 and 1e-6 < tp < 1 - 1e-6
    frames.append({"a": float(a), "L": float(L),
                   "tris": [np.round(cells[k][f], 6).tolist() for (k, f) in tris],
                   "jOA": np.round(jOA, 6).tolist(), "jBD": np.round(jBD, 6).tolist(),
                   "sep": float(np.linalg.norm(d)), "crossing": bool(crossing),
                   "pair": [list(l1), list(l2)] if crossing else None})
    if abs(a) in (15.0, 7.5, 2.5, 0.0) :
        print(f"a={a:+6.1f}: joints {np.linalg.norm(d):.4f} apart, direction {np.round(d/ (np.linalg.norm(d) or 1), 3)}, crossing {crossing}")
# is the pass head-on? compare direction at -a and +a
def dirat(a):
    f = next(fr for fr in frames if abs(fr["a"] - a) < 1e-9)
    d = np.array(f["jBD"]) - np.array(f["jOA"]); return d / np.linalg.norm(d)
print("cos(angle) between separation directions at -7.5 and +7.5:", round(float(dirat(-7.5) @ dirat(7.5)), 4), "(-1 = straight through)")
out = {"strut": IC.EL, "L0": L0, "vertex": VTX.tolist(), "tris": [[k, f] for (k, f) in tris],
       "names": NAMES, "frames": frames}
json.dump(out, open(str(common.out("joints")), "w"))
print("frames", len(frames))
