"""Triangle network: every triangle a rigid body, every junction the same blob. No cells."""
import sys, time, collections, itertools as it
sys.path.insert(0, "/Users/hal.hildebrand/git/inviscid/analysis/jitterbug-variety")
import numpy as np
import jb_rc_reduced as RC
import jb_mj_inertial_honeycomb as MJ

A = MJ.A_REF
def hat(w): return np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])

class Net:
    """Triangles from a single-covering patch at phase A; junctions by coincidence."""
    def __init__(self, sites, periodic=None):
        asm, _ = RC.honeycomb_single(sites, gc=A)
        X = asm.positions(asm.q0())              # (N, 12, 3)
        self.sites = [tuple(s) for s in sites]; self.N = asm.N
        self.tri = []                            # list of (cell, face, centroid, corners(3,3))
        for k in range(asm.N):
            for f in range(8):
                P = X[k][list(RC.TRIS[f])]
                self.tri.append((k, f, P.mean(0), P))
        self.T = len(self.tri)
        # corner instances -> junctions by position
        self.corners = [(t, i) for t in range(self.T) for i in range(3)]
        pos = {c: self.tri[c[0]][3][c[1]] for c in self.corners}
        groups = collections.defaultdict(list)
        for c, p in pos.items():
            groups[tuple(np.round(p, 6))].append(c)
        self.junctions = list(groups.values())
        self.valence = collections.Counter(len(g) for g in self.junctions)
        self.pos = pos
    def corner_rows(self, c):
        """3x6 Jacobian of corner c's velocity w.r.t. its triangle's (cdot, omega)."""
        t, i = c
        r = self.pos[c] - self.tri[t][2]
        J = np.zeros((3, 6)); J[:, 0:3] = np.eye(3); J[:, 3:6] = -hat(r)
        return J
    def mass(self):
        M = np.zeros((self.T, 6, 6))
        for t in range(self.T):
            for i in range(3):
                J = self.corner_rows((t, i)); M[t] += (1.0 / 3.0) * J.T @ J
        return M
    def rigid_C(self):
        rows = []
        for g in self.junctions:
            J0 = self.corner_rows(g[0]); t0 = g[0][0]
            for c in g[1:]:
                row = np.zeros((3, 6 * self.T))
                row[:, 6 * c[0]:6 * c[0] + 6] += self.corner_rows(c)
                row[:, 6 * t0:6 * t0 + 6] -= J0
                rows.append(row)
        return np.vstack(rows)
    def blob_B(self):
        """rows of x_c - mean over the junction, for V = k/2 sum |.|^2"""
        rows = []
        for g in self.junctions:
            m = len(g)
            Js = [(c[0], self.corner_rows(c)) for c in g]
            for (tc, Jc) in Js:
                row = np.zeros((3, 6 * self.T))
                row[:, 6 * tc:6 * tc + 6] += Jc
                for (t2, J2) in Js:
                    row[:, 6 * t2:6 * t2 + 6] -= J2 / m
                rows.append(row)
        return np.vstack(rows)

def spectrum(B, M, k=1.0):
    n = M.shape[0] * 6
    H = k * (B.T @ B)
    Mf = np.zeros((n, n))
    for t in range(M.shape[0]): Mf[6 * t:6 * t + 6, 6 * t:6 * t + 6] = M[t]
    L = np.linalg.cholesky(Mf)
    Aq = np.linalg.solve(L, np.linalg.solve(L, H).T).T
    ev = np.linalg.eigvalsh((Aq + Aq.T) / 2)
    return np.sqrt(np.clip(ev, 0, None)), ev

def block(side):
    return [(2 * x, 2 * y, 2 * z) for x in range(side) for y in range(side) for z in range(side)]

if __name__ == "__main__":
    for side in (2, 3):
        t0 = time.time()
        net = Net(block(side))
        C = net.rigid_C()
        rank = np.linalg.matrix_rank(C, tol=1e-9)
        nul = 6 * net.T - rank
        w, ev = spectrum(net.blob_B(), net.mass())
        cut = 1e-9 * ev.max()
        nz = int((ev < cut).sum())
        print(f"block {side}^3: cells {net.N}, triangles {net.T}, DOF {6*net.T}, junctions {len(net.junctions)} (valence {dict(sorted(net.valence.items()))})")
        print(f"   RIGID junctions: constraint rank {rank}, nullity {nul} = 6 rigid + {nul-6} internal mechanisms   (cell model: 1)")
        print(f"   SOFT junctions (k=1): zero modes {nz} of {6*net.T}; lowest nonzero {w[nz]:.5f}, highest {w[-1]:.5f}   ({time.time()-t0:.0f}s)")
