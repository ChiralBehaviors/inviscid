"""A single jitterbug driven -60 .. +60 (octahedron, VE, octahedron).
Per frame: the 8 triangles (24 corners); which corners coincide; whether
the 12 joints are the SAME corner pairs on every frame; distinct vertex
count (Fuller: 12, congruent as 6 at the octahedron); strut length; the
nearest pair of struts that do not share a joint (self-clearance); the
driven triangle's rotation about its own axis."""
import itertools as it, json
import numpy as np

from analysis.model.first_principles.pages import common
from analysis.model import plates as Z
from analysis.model import cell as IC

def segdist(p1, q1, p2, q2):
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    b, c = d1 @ d2, d1 @ r
    den = a * e - b * b
    s = np.clip((b * f - c * e) / den, 0, 1) if den > 1e-14 else 0.0
    t = (b * s + f) / e
    if t < 0: t = 0.0; s = np.clip(-c / a, 0, 1)
    elif t > 1: t = 1.0; s = np.clip((b - c) / a, 0, 1)
    return float(np.linalg.norm((p1 + d1 * s) - (p2 + d2 * t)))

def classes(C, tol=1e-9):
    pts = [(f, c) for f in range(8) for c in range(3)]
    cls = []
    for fc in pts:
        p = C[fc[0]][fc[1]]
        for g in cls:
            if np.linalg.norm(C[g[0][0]][g[0][1]] - p) < tol:
                g.append(fc); break
        else:
            cls.append([fc])
    return cls

# the joints: read at a generic angle and held as identities
ref = classes(Z.corners(-30.0))
assert len(ref) == 12 and all(len(g) == 2 for g in ref)
joints = [tuple(sorted(g)) for g in ref]
faces = Z.faces()
axes = [f[2] for f in faces]; sig = [f[3] for f in faces]
DRIVEN = 0
frames = []
for a in np.arange(-60.0, 60.0 + 1e-9, 2.0):
    C = Z.corners(a)
    # every joint still coincident?
    jd = max(float(np.linalg.norm(C[f1][c1] - C[f2][c2])) for ((f1, c1), (f2, c2)) in joints)
    cl = classes(C)
    nvert = len(cl)
    # strut lengths
    sl = [float(np.linalg.norm(C[f][i] - C[f][j])) for f in range(8) for i, j in ((0,1),(1,2),(2,0))]
    # self clearance: struts of different triangles not sharing a joint
    S = [((f, i, j), C[f][i], C[f][j]) for f in range(8) for i, j in ((0,1),(1,2),(2,0))]
    dmin = np.inf
    for (l1, p1, q1), (l2, p2, q2) in it.combinations(S, 2):
        if l1[0] == l2[0]: continue
        if min(np.linalg.norm(x - y) for x in (p1, q1) for y in (p2, q2)) < 1e-9: continue
        dmin = min(dmin, segdist(p1, q1, p2, q2))
    # driven triangle: rotation about its own axis relative to the VE
    u = axes[DRIVEN]
    v0 = Z.corners(0.0)[DRIVEN][0] - Z.corners(0.0)[DRIVEN].mean(0)
    v1 = C[DRIVEN][0] - C[DRIVEN].mean(0)
    v0p = v0 - u * (v0 @ u); v1p = v1 - u * (v1 @ u)
    ang = float(np.degrees(np.arctan2(np.cross(v0p, v1p) @ u, v0p @ v1p)))
    h = float(C[DRIVEN].mean(0) @ u)
    frames.append({"a": float(a), "tris": np.round(C, 6).tolist(), "joint_max_gap": jd,
                   "nvert": nvert, "strut": [min(sl), max(sl)], "clear": float(dmin),
                   "spin": ang, "height": h,
                   "verts": [np.round(C[g[0][0]][g[0][1]], 6).tolist() for g in cl]})
    if a in (-60.0, -30.0, 0.0, 30.0, 60.0):
        print(f"a={a:+5.0f}: joints max gap {jd:.1e}  distinct vertices {nvert:2d}  strut {min(sl):.6f}..{max(sl):.6f}  "
              f"self-clearance {dmin:.4f}  driven triangle spin {ang:+7.2f} deg, height {h:.4f}")
out = {"joints": [[list(p), list(q)] for p, q in joints], "driven": DRIVEN,
       "axis": np.round(axes[DRIVEN], 6).tolist(), "sigma": [float(s) for s in sig],
       "axes": np.round(axes, 6).tolist(), "strut": IC.EL, "frames": frames}
json.dump(out, open(str(common.out("one_cell")), "w"))
print("frames", len(frames), "joints permanent:", all(f["joint_max_gap"] < 1e-9 for f in frames))
