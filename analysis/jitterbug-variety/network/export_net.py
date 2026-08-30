import sys, json, base64, time, collections
sys.path.insert(0, "/Users/hal.hildebrand/git/inviscid/analysis/jitterbug-variety")
sys.path.insert(0, "/private/tmp/claude-501/-Users-hal-hildebrand-git-inviscid/e58634db-a827-44db-9d7f-54d7fe219a1b/scratchpad")
import numpy as np, scipy.sparse as sp
import jb_rc_reduced as RC, jb_mj_inertial_honeycomb as MJ
from tri_network import Net, hat, A

S = sys.argv[1]
d = json.load(open(S + "/frames.json"))
sites = [tuple(s) for s in d["sites"]]; centre = d["centre"]
net = Net(sites); T = net.T
print(f"network: {net.N} cells, {T} triangles, {6*T} DOF, {len(net.junctions)} junctions {dict(sorted(net.valence.items()))}")
# sparse blob rows and mass
rows, cols, vals, r = [], [], [], 0
for g in net.junctions:
    m = len(g); Js = [(c[0], net.corner_rows(c)) for c in g]
    for (tc, Jc) in Js:
        blk = collections.defaultdict(lambda: np.zeros((3, 6)))
        blk[tc] += Jc
        for (t2, J2) in Js: blk[t2] -= J2 / m
        for t2, Bk in blk.items():
            for i in range(3):
                for j in range(6):
                    if Bk[i, j] != 0: rows.append(r + i); cols.append(6 * t2 + j); vals.append(Bk[i, j])
        r += 3
B = sp.csr_matrix((vals, (rows, cols)), shape=(r, 6 * T))
K = 1.0
H = (K * (B.T @ B)).tocsr()
M = net.mass(); Minv = np.array([np.linalg.inv(M[t]) for t in range(T)])
def accel(u):
    f = -(H @ u)
    return np.concatenate([Minv[t] @ f[6*t:6*t+6] for t in range(T)])
# the same kick: the centre cell's fold rate 0.2, mapped onto its eight triangles
asm, _ = RC.honeycomb_single(sites, gc=A)
ctr0, R0, gam0, B0 = asm.frames(asm.q0()); Jcell = asm.cell_jacobians(ctr0, R0, B0)
uc = np.zeros(7); uc[6] = 0.2
V = (Jcell[centre] @ uc).reshape(12, 3)
v = np.zeros(6 * T)
for t, (k, f, c, P) in enumerate(net.tri):
    if k != centre: continue
    vs = V[list(RC.TRIS[f])]; rs = P - c
    Amat = np.vstack([np.hstack([np.eye(3), -hat(rr)]) for rr in rs])
    v[6*t:6*t+6] = np.linalg.lstsq(Amat, vs.ravel(), rcond=None)[0]
# normals and face senses for the fold-equivalent colour
normals = np.zeros((T, 3)); sense = np.zeros(T)
for t, (k, f, c, P) in enumerate(net.tri):
    n = c - ctr0[k]; normals[t] = n / np.linalg.norm(n); sense[t] = RC.FACES[f][3]
tri_cell = [net.tri[t][0] for t in range(T)]
corner_index = {c: i for i, c in enumerate(net.corners)}
junction_groups = [[corner_index[c] for c in g] for g in net.junctions]
P0 = np.array([net.tri[t][3] for t in range(T)])          # (T,3,3)
C0 = np.array([net.tri[t][2] for t in range(T)])
def positions(x):
    out = np.empty_like(P0)
    for t in range(T):
        dc = x[6*t:6*t+3]; th = x[6*t+3:6*t+6]; a = np.linalg.norm(th)
        Rm = np.eye(3) if a < 1e-15 else np.eye(3) + np.sin(a) / a * hat(th) + (1 - np.cos(a)) / a**2 * (hat(th) @ hat(th))
        out[t] = C0[t] + dc + (Rm @ (P0[t] - C0[t]).T).T
    return out
def kinetic_split(u):
    ks, kt, kp = 0.0, 0.0, 0.0
    for t in range(T):
        w = u[6*t+3:6*t+6]; n = normals[t]
        wn = np.dot(w, n) * n; wp = w - wn
        kt += 0.5 * u[6*t:6*t+3] @ M[t][0:3, 0:3] @ u[6*t:6*t+3]
        ks += 0.5 * wn @ M[t][3:6, 3:6] @ wn
        kp += 0.5 * wp @ M[t][3:6, 3:6] @ wp
    return ks, kt, kp
x = np.zeros(6 * T); u = v.copy(); h, steps, every = 0.01, 1200, 12
frames, foldeq, energy, times = [positions(x)], [np.zeros(T)], [], [0.0]
ks, kt, kp = kinetic_split(u); pe = 0.5 * K * float(u @ (H @ u)) * 0  # PE at x=0 is 0
E0 = ks + kt + kp; energy.append(dict(spin=ks, trans=kt, tilt=kp, pe=0.0))
a = np.zeros(6 * T); t0 = time.time()   # x = 0 at t = 0, so the force is zero
for s in range(1, steps + 1):
    u_h = u + 0.5 * h * a; x = x + h * u_h
    a = -(H @ x); a = np.concatenate([Minv[t] @ a[6*t:6*t+6] for t in range(T)])
    u = u_h + 0.5 * h * a
    if s % every == 0:
        frames.append(positions(x)); times.append(s * h)
        th = x.reshape(T, 6)[:, 3:6]
        foldeq.append(np.degrees(np.einsum('ij,ij->i', th, normals)) * sense)   # rotation about own normal, in the fold's sign
        ks, kt, kp = kinetic_split(u); pe = 0.5 * K * float(x @ (H @ x))
        energy.append(dict(spin=ks, trans=kt, tilt=kp, pe=pe))
print(f"{steps} steps in {time.time()-t0:.0f}s; E0 {E0:.6f}, E end {sum(energy[-1].values()):.6f}")
F = np.array(frames); FE = np.array(foldeq)
disp = np.abs(F - F[0]).max() / d["strut"]
gaps = max(max(np.linalg.norm(F[:, c1 // 3, c1 % 3] - F[:, c2 // 3, c2 % 3], axis=-1).max() for c1 in g for c2 in g) for g in junction_groups if len(g) > 1)
cent = [t for t in range(T) if tri_cell[t] == centre]
print(f"max vertex displacement {100*disp:.1f}% of strut; max junction opening {gaps:.4f}; centre-cell fold-equivalent max {np.abs(FE[:, cent]).max():.2f} deg; overall max {np.abs(FE).max():.2f} deg")
print("energy split by frame (spin about normal / translation / tilt / joints), fractions of E0:")
for i in range(0, len(energy), 10):
    e = energy[i]; print(f"   t={times[i] if i < len(times) else '?':>5}  {e['spin']/E0:.3f} / {e['trans']/E0:.3f} / {e['tilt']/E0:.3f} / {e['pe']/E0:.3f}")
def q16(arr):
    arr = np.asarray(arr, np.float64); sc = float(np.abs(arr).max()) / 32000.0
    return base64.b64encode(np.round(arr / sc).astype(np.int16).tobytes()).decode("ascii"), sc
pos_b64, sc = q16(F)
d["net"] = dict(T=T, tri_cell=tri_cell, junctions=junction_groups, times=times, pos=pos_b64, scale=sc,
                foldeq=[[float(v) for v in row] for row in FE], energy=energy[1:], kick=0.2, k=K, h=h,
                maxdisp=[float(np.abs(F[i] - F[0]).max() / d["strut"]) for i in range(F.shape[0])], maxgap=float(gaps))
json.dump(d, open(S + "/frames.json", "w"))
print("wrote", len(open(S + "/frames.json").read()))
