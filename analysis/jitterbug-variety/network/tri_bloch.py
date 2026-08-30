import sys, time, collections
sys.path.insert(0, "/Users/hal.hildebrand/git/inviscid/analysis/jitterbug-variety")
sys.path.insert(0, "/private/tmp/claude-501/-Users-hal-hildebrand-git-inviscid/e58634db-a827-44db-9d7f-54d7fe219a1b/scratchpad")
import numpy as np
import jb_rc_reduced as RC, jb_mj_inertial_honeycomb as MJ
from tri_network import Net, hat, A, block

# ---- periodic triangle network: 8 triangles per primitive cell ----
four, _ = RC.honeycomb_single([(0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2)], gc=A)
X = four.positions(four.q0())
home = X[0]                                   # (12,3) home cell vertices (centre at origin)
tri_c = [home[list(RC.TRIS[f])].mean(0) for f in range(8)]
corners_of = collections.defaultdict(list)    # vertex -> [(f, i)]
for f in range(8):
    for i, v in enumerate(RC.TRIS[f]): corners_of[v].append((f, i))
def Jc(f, i):
    r = home[RC.TRIS[f][i]] - tri_c[f]
    J = np.zeros((3, 6)); J[:, 0:3] = np.eye(3); J[:, 3:6] = -hat(r); return J
M = np.zeros((8, 6, 6))
for f in range(8):
    for i in range(3): M[f] += (1.0 / 3.0) * Jc(f, i).T @ Jc(f, i)
Mf = np.zeros((48, 48))
for f in range(8): Mf[6*f:6*f+6, 6*f:6*f+6] = M[f]
Lm = np.linalg.cholesky(Mf)
bonds = [(np.array(e), four.welds[n][2]) for n, e in enumerate(((1,0,0),(0,1,0),(0,0,1)))]  # lattice vectors, tied (a,b)
assert all(w[0] == 0 for w in four.welds) and len(four.welds) == 3
def bands(kvec, k=1.0):
    rows = []
    for e, prs in bonds:
        ph = np.exp(1j * float(np.dot(kvec, e)))
        for (a, b) in prs:
            cs = [(f, i, 1.0) for (f, i) in corners_of[a]] + [(f, i, ph) for (f, i) in corners_of[b]]
            m = len(cs)
            for (f, i, p) in cs:
                row = np.zeros((3, 48), complex)
                row[:, 6*f:6*f+6] += p * Jc(f, i)
                for (f2, i2, p2) in cs: row[:, 6*f2:6*f2+6] -= p2 * Jc(f2, i2) / m
                rows.append(row)
    B = np.vstack(rows)
    H = k * (B.conj().T @ B)
    Aq = np.linalg.solve(Lm, np.linalg.solve(Lm, H).conj().T).conj().T
    ev = np.linalg.eigvalsh((Aq + Aq.conj().T) / 2)
    return np.sqrt(np.clip(ev, 0, None)), ev
w0, ev0 = bands(np.zeros(3))
cut = 1e-9 * ev0.max()
print(f"PERIODIC triangle network: {len(w0)} bands; zero at Gamma: {int((ev0 < cut).sum())}; Gamma spectrum (omega^2, distinct): {sorted({round(float(x), 6) for x in ev0})}")
wg, evg = bands(np.array([0.7, 1.1, 1.9]))
print(f"   zero at a generic k: {int((evg < 1e-9 * evg.max()).sum())}  (flat zero bands = local mechanisms)")
for kk in (1e-3, 1e-4):
    w, _ = bands(np.array([kk, 0.0, 0.0]))
    print(f"   [100] k={kk:g}: lowest 8 omega/k = {np.round(np.sort(w)[:8] / kk, 4)}")
w, _ = bands(np.array([1e-4, 1e-4, 0.0]) / np.sqrt(2)); print(f"   [110] k=1e-4: lowest 8 omega/k = {np.round(np.sort(w)[:8] / 1e-4, 4)}")
w, _ = bands(np.array([1e-4]*3) / np.sqrt(3)); print(f"   [111] k=1e-4: lowest 8 omega/k = {np.round(np.sort(w)[:8] / 1e-4, 4)}")
ks = np.linspace(0, np.pi, 13); mx = 0; arg = None
for i, kx in enumerate(ks):
    for j, ky in enumerate(ks[:i+1]):
        for kz in ks[:j+1]:
            w, _ = bands(np.array([kx, ky, kz]))
            if w.max() > mx: mx, arg = w.max(), (kx, ky, kz)
print(f"   band max {mx:.6f} at k={np.round(arg, 3)}  (cell model: sqrt3 = 1.732051 at R)")
wR, _ = bands(np.array([np.pi]*3)); print(f"   at R: distinct omega^2 {sorted({round(float(x**2), 5) for x in wR})}")

# ---- where do a finite block's mechanisms live? side 4: 8 interior cells ----
t0 = time.time()
net = Net(block(4))
C = net.rigid_C()
U, s, Vt = np.linalg.svd(C, full_matrices=True)
rank = int((s > 1e-9).sum()); Z = Vt[rank:].T
interior = [k for k, st in enumerate(net.sites) if all(2 <= c <= 4 for c in st)]
tri_int = [t for t in range(net.T) if net.tri[t][0] in interior]
wts = np.array([np.linalg.norm(Z[[6*t+d for t in tri_int for d in range(6)], j]) ** 2 for j in range(Z.shape[1])])
print(f"block 4^3: {net.N} cells ({len(interior)} interior), nullity {Z.shape[1]} = 6 + {Z.shape[1]-6} mechanisms; "
      f"null-space weight on interior triangles: total {wts.sum():.3f} of {Z.shape[1]} (interior fraction of DOF {len(tri_int)/net.T:.3f}); "
      f"vectors with >50% interior weight: {int((wts > 0.5).sum())}   ({time.time()-t0:.0f}s)")
# which motion is the interior one? test the coherent breathe: every cell's fold rate 1 with centres on dL
q, u = RC.coherent(net.sites, A, rate=1.0)
asm, _ = RC.honeycomb_single(net.sites, gc=A)
ctr, R, gam, B = asm.frames(asm.q0()); J = asm.cell_jacobians(ctr, R, B)
V = np.array([J[k] @ u[k] for k in range(asm.N)])       # vertex velocities (N,36)->(N,12,3)
V = V.reshape(asm.N, 12, 3)
# map to triangle rigid velocities: fit (cdot, omega) per triangle from its 3 corner velocities
vt = np.zeros(6 * net.T)
for t, (k, f, c, P) in enumerate(net.tri):
    vs = V[k][list(RC.TRIS[f])]; rs = P - c
    Amat = np.vstack([np.hstack([np.eye(3), -hat(r)]) for r in rs]); vt[6*t:6*t+6] = np.linalg.lstsq(Amat, vs.ravel(), rcond=None)[0]
res = np.linalg.norm(C @ vt) / np.linalg.norm(vt)
print(f"   coherent breathe as triangle motion: constraint residual {res:.1e} (0 = it is one of the mechanisms)")
