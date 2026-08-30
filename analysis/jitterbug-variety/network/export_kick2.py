import sys, json, base64, time
sys.path.insert(0, "/Users/hal.hildebrand/git/inviscid/analysis/jitterbug-variety")
import numpy as np
import jb_rc_reduced as RC
import jb_mj_inertial_honeycomb as MJ
import jb_je_joint_exponent as JE

S = sys.argv[1]
d = json.load(open(S + "/frames.json"))
sites = [tuple(s) for s in d["sites"]]; centre = d["centre"]; N = d["N"]
asm, _ = RC.honeycomb_single(sites, gc=MJ.A_REF)
pairs = MJ.tied_pairs(asm)
KICK, K, H, STEPS, EVERY = 0.2, 1.0, 5e-3, 2400, 24
q = asm.q0(); u = np.zeros((N, 7)); u[centre, 6] = KICK
shell = np.array([max(abs(c) for c in s) // 2 for s in sites])  # 0..2, Chebyshev shell
B_ref = RC.body(MJ.A_REF)

def split(q, u):
    ctr, R, gam, B = asm.frames(q)
    J = asm.cell_jacobians(ctr, R, B)
    M = asm.mass_blocks(J)
    kf = np.array([0.5 * M[k][6, 6] * u[k, 6] ** 2 for k in range(N)])
    kt = np.array([0.5 * u[k, 0:3] @ M[k][0:3, 0:3] @ u[k, 0:3] for k in range(N)])
    kr = np.array([0.5 * u[k, 3:6] @ M[k][3:6, 3:6] @ u[k, 3:6] for k in range(N)])
    kk = np.array([0.5 * u[k] @ M[k] @ u[k] for k in range(N)])
    return kf, kt, kr, kk, ctr, R, gam

def decomp(ctr, R, gam):
    fold_only = np.array([asm.ctr0[k] + RC.body(float(gam[k])) for k in range(N)])
    rigid_only = np.array([ctr[k] + (R[k] @ B_ref.T).T for k in range(N)])
    return fold_only, rigid_only

def q16(arr):
    a = np.asarray(arr, np.float64); scale = float(np.abs(a).max()) / 32000.0
    return base64.b64encode(np.round(a / scale).astype(np.int16).tobytes()).decode("ascii"), scale

rows, fo, ro = [], [], []
def record(t, q, u, sep):
    kf, kt, kr, kk, ctr, R, gam = split(q, u)
    pe = 0.5 * K * float(np.sum(sep ** 2))
    per_shell = [[float(kf[shell == s].sum()), float(kt[shell == s].sum()), float(kr[shell == s].sum())] for s in range(3)]
    rows.append(dict(t=t, fold=float(kf.sum()), trans=float(kt.sum()), rot=float(kr.sum()), cross=float(kk.sum() - kf.sum() - kt.sum() - kr.sum()), pe=pe, shells=per_shell))
    f_, r_ = decomp(ctr, R, gam); fo.append(f_); ro.append(r_)

a1, sep, M = JE.state(asm, q, u, pairs, 2.0, K); record(0.0, q, u, sep)
t0 = time.time()
for s in range(1, STEPS + 1):
    a1, sep, M = JE.state(asm, q, u, pairs, 2.0, K)
    u_h = u + 0.5 * H * a1
    q_h = RC.apply_increment(asm, q, (0.5 * H * u).ravel())
    a2, _, _ = JE.state(asm, q_h, u_h, pairs, 2.0, K)
    u = u + H * a2
    q = RC.apply_increment(asm, q, (H * u_h).ravel())
    if s % EVERY == 0:
        _, sep, _ = JE.state(asm, q, u, pairs, 2.0, K)
        record(s * H, q, u, sep)
print(f"{STEPS} steps in {time.time()-t0:.0f}s")
E = rows[0]["fold"] + rows[0]["pe"]
print(f"{'t':>5s} {'fold':>7s} {'trans':>7s} {'rot':>7s} {'cross':>7s} {'PE':>7s}   shell0 f/t/r      shell1 f/t/r      shell2 f/t/r   (fractions of E0 = {E:.4f})")
for r in rows[::10]:
    sh = "  ".join(f"{a/E:5.3f}/{b/E:5.3f}/{c/E:5.3f}" for a, b, c in r["shells"])
    print(f"{r['t']:5.2f} {r['fold']/E:7.3f} {r['trans']/E:7.3f} {r['rot']/E:7.3f} {r['cross']/E:7.3f} {r['pe']/E:7.3f}   {sh}")
d["kick"]["energy_split"] = rows
d["kick"]["pos_fold"], d["kick"]["scale_fold"] = q16(np.array(fo))
d["kick"]["pos_rigid"], d["kick"]["scale_rigid"] = q16(np.array(ro))
json.dump(d, open(S + "/frames.json", "w"))
print("wrote frames.json", len(open(S + "/frames.json").read()))
