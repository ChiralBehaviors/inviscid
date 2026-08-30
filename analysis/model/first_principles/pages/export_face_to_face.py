"""VE + adjacent octahedron cell, welded on one triangular face. The pair's
configuration space is (a, b) with the octa's pose fixed by the shared
triangle. Export: cell-1 corners per a; cell-2 body corners per b; the pose
(R, t) per (a, b) from the three welded corners (checked against the
model's cluster builder on the b = a + 60 line); per (a, b): the octa's
turn about the shared axis, centre distance, corners per shared joint,
nearest strut pair between the cells."""
import itertools as it, json
import numpy as np

from analysis.model.first_principles.pages import common
from analysis.model import assembly as RC
from analysis.model import plates as Z
from analysis.model import cell as IC

asm0 = RC.cluster(gc=0.0, sites=[(1, 1, 1)])
(_, _, PAIRS) = asm0.welds[0]                      # [(VE label, octa label)] x3
AX = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
STEP = 3.0
ANG = [int(x) for x in np.arange(-60.0, 60.0 + 1e-9, STEP)]

def kabsch(P, Q):  # R, t with R P + t = Q  (P, Q: (3,3) rows)
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Dm = np.diag([1, 1, d])
    R = Vt.T @ Dm @ U.T
    return R, qc - R @ pc

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

def struts(P):
    return [(f, P[f][i], P[f][j]) for f in range(8) for i, j in ((0,1),(1,2),(2,0))]

cell1 = {a: Z.corners(a) for a in ANG}
verts1 = {a: IC.cell_verts(a, np.zeros(3)) for a in ANG}
body2 = {b: Z.corners(b) for b in ANG}
vbody2 = {b: IC.cell_verts(b, np.zeros(3)) for b in ANG}
pose, turn, dist, corners, near = {}, {}, {}, {}, {}
maxres = 0.0; checks = 0
for a in ANG:
    P1 = np.array([verts1[a][i] for i, _ in PAIRS])
    for b in ANG:
        P2 = np.array([vbody2[b][j] for _, j in PAIRS])
        R, t = kabsch(P2, P1)
        res = np.abs((R @ P2.T).T + t - P1).max(); maxres = max(maxres, res)
        pose[(a, b)] = (R, t)
        # octa's turn about the shared axis: angle of R about AX
        ang = np.degrees(np.arctan2(((R[2,1]-R[1,2])*AX[0] + (R[0,2]-R[2,0])*AX[1] + (R[1,0]-R[0,1])*AX[2]), np.trace(R) - 1))
        turn[(a, b)] = float(ang)
        c2 = t + R @ np.zeros(3)
        dist[(a, b)] = float(np.linalg.norm(c2))
        # corners on each shared joint
        W2 = (R @ body2[b].reshape(-1, 3).T).T + t
        cnt = []
        for i, _ in PAIRS:
            p = verts1[a][i]
            n1 = sum(1 for x in cell1[a].reshape(-1, 3) if np.linalg.norm(x - p) < 1e-9)
            n2 = sum(1 for x in W2 if np.linalg.norm(x - p) < 1e-9)
            cnt.append([n1, n2])
        corners[(a, b)] = cnt
        # nearest struts between the cells, excluding pairs sharing a joint
        S1 = struts(cell1[a]); S2 = struts(W2.reshape(8, 3, 3))
        dmin = np.inf
        for (f1, p1, q1) in S1:
            for (f2, p2, q2) in S2:
                if min(np.linalg.norm(x - y) for x in (p1, q1) for y in (p2, q2)) < 1e-9: continue
                dmin = min(dmin, segdist(p1, q1, p2, q2))
        near[(a, b)] = float(dmin)
        # check against the model's builder on the exchange line
        if abs(b - (a + 60.0)) < 1e-9 and -60 <= a <= 0:
            asm = RC.cluster(gc=a, sites=[(1, 1, 1)])
            cm = asm.q0().reshape(-1, 8)[1, 0:3]
            assert np.allclose(cm, t, atol=1e-9) and np.allclose(R, np.eye(3), atol=1e-9), (a, b, cm, t)
            checks += 1
print(f"poses {len(pose)}, max weld residual {maxres:.1e}, builder checks on b=a+60: {checks} passed")
for (a, b) in ((0, 60), (0, 30), (-30, 30), (-60, 0), (30, 60), (0, -60), (60, 60), (-60, -60)):
    print(f"  a={a:+4d} b={b:+4d}: octa turned {turn[(a,b)]:+7.2f} deg, centres {dist[(a,b)]:.4f}, corners per shared joint {corners[(a,b)]}, nearest struts {near[(a,b)]:.4f}")
out = {"angles": ANG, "pairs": [list(p) for p in PAIRS], "axis": AX.tolist(), "strut": IC.EL,
       "cell1": {str(a): np.round(cell1[a], 6).tolist() for a in ANG},
       "body2": {str(b): np.round(body2[b], 6).tolist() for b in ANG},
       "shared1": {str(a): np.round([verts1[a][i] for i, _ in PAIRS], 6).tolist() for a in ANG},
       "pose": {f"{a}|{b}": [np.round(pose[(a,b)][0], 8).tolist(), np.round(pose[(a,b)][1], 8).tolist()] for a in ANG for b in ANG},
       "turn": {f"{a}|{b}": round(turn[(a,b)], 4) for a in ANG for b in ANG},
       "dist": {f"{a}|{b}": round(dist[(a,b)], 5) for a in ANG for b in ANG},
       "corners": {f"{a}|{b}": corners[(a,b)] for a in ANG for b in ANG},
       "near": {f"{a}|{b}": round(near[(a,b)], 5) for a in ANG for b in ANG}}
json.dump(out, open(str(common.out("face_to_face")), "w"))
print("written")
