import sys, json, base64, time
sys.path.insert(0, "/Users/hal.hildebrand/git/inviscid/analysis/jitterbug-variety")
import numpy as np
import jb_rc_reduced as RC
import jb_mj_inertial_honeycomb as MJ
import jb_je_joint_exponent as JE

OUT = sys.argv[1]
R = 4  # all-even sites -4..4 -> 5x5x5 = 125 cells
sites = [(x, y, z) for x in range(-R, R + 1, 2) for y in range(-R, R + 1, 2) for z in range(-R, R + 1, 2)]
centre = sites.index((0, 0, 0))

def q16(arr):
    a = np.asarray(arr, np.float64)
    scale = float(np.abs(a).max()) / 32000.0
    q = np.round(a / scale).astype(np.int16)
    return base64.b64encode(q.tobytes()).decode("ascii"), scale

# ---- breathe: every cell at the same phase, centres on L(a) ----
phases = np.linspace(0.0, -60.0, 61)
pos_b, L_b = [], []
for a in phases:
    asm, _ = RC.honeycomb_single(sites, gc=float(a))
    pos_b.append(asm.positions(asm.q0()))
    L_b.append(float(RC.lattice_constant(float(a))))
pos_b = np.array(pos_b)  # (F, N, 12, 3)
b64_b, sc_b = q16(pos_b)

# ---- kick: soft joint k=1, corner point masses, centre fold rate 0.5 rad/time ----
asm, _ = RC.honeycomb_single(sites, gc=MJ.A_REF)
pairs = MJ.tied_pairs(asm)
KICK, K, H, STEPS, EVERY = 0.2, 1.0, 5e-3, 2400, 24
q = asm.q0(); u = np.zeros((asm.N, 7)); u[centre, 6] = KICK
pos_k, gam_k, times = [asm.positions(q)], [RC.Assembly.unpack(q)[2].copy()], [0.0]
e0 = None; t0 = time.time()
for s in range(1, STEPS + 1):
    a1, sep, M = JE.state(asm, q, u, pairs, 2.0, K)
    if e0 is None:
        e0 = MJ.kinetic(M, u) + 0.5 * K * float(np.sum(sep ** 2))
    u_h = u + 0.5 * H * a1
    q_h = RC.apply_increment(asm, q, (0.5 * H * u).ravel())
    a2, _, _ = JE.state(asm, q_h, u_h, pairs, 2.0, K)
    u = u + H * a2
    q = RC.apply_increment(asm, q, (H * u_h).ravel())
    if s % EVERY == 0:
        pos_k.append(asm.positions(q)); gam_k.append(RC.Assembly.unpack(q)[2].copy()); times.append(s * H)
a1, sep, M = JE.state(asm, q, u, pairs, 2.0, K)
e1 = MJ.kinetic(M, u) + 0.5 * K * float(np.sum(sep ** 2))
pos_k = np.array(pos_k); gam_k = np.array(gam_k)
b64_k, sc_k = q16(pos_k)
dev = gam_k - MJ.A_REF
print(f"kick run: {STEPS} steps in {time.time()-t0:.0f}s; energy {e0:.6f} -> {e1:.6f}; centre fold excursion max {np.abs(dev[:, centre]).max():.2f} deg; far-corner max {np.abs(dev[:, sites.index((R,R,R))]).max():.3f} deg; overall max dev {np.abs(dev).max():.2f}")
# joint separations at the end, for the record
print(f"max joint separation at end {sep.max():.4f} (strut {RC.EL:.4f})")

data = dict(
    N=asm.N, sites=sites, centre=centre, tris=[list(map(int, t)) for t in RC.TRIS],
    bars=[list(map(int, b)) for b in RC.BARS], strut=float(RC.EL),
    bonds=[[int(i), int(j), [[int(a), int(b)] for (a, b) in ps]] for (i, j, ps) in asm.welds],
    breathe=dict(phases=[float(a) for a in phases], L=L_b, pos=b64_b, scale=sc_b),
    kick=dict(times=times, gamma=[[float(g) for g in row] for row in gam_k], pos=b64_k, scale=sc_k,
              kick=KICK, k=K, h=H, energy=[e0, e1]),
)
with open(OUT, "w") as f:
    json.dump(data, f)
print("wrote", OUT, "bytes", len(open(OUT).read()))
