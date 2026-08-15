"""Step S: THE FIRST FREQUENCY SPECTRUM. USER DECISION 17, bead inviscid-qvf.2.

Everything this project has computed about stability so far was INERTIA -- the
count of positive/zero/negative Hessian eigenvalues, which is metric-independent
by Sylvester's law and therefore says nothing about how fast anything moves.
This file combines the potential Hessian (jb_k / jb_o machinery) with the
internal mass metric (jb_r) and solves

    H v = omega^2 M v

on the six-dimensional internal tangent space at the vector equilibrium.

DECISION 17 fixes the ground state (the VE) and the variant (RAW, on real
distances -- the struts are rigid, so there is no free scale DOF and normalising
to unit radius evaluates the energy at a configuration the linkage cannot
occupy). It does NOT fix the kernel, and the kernel is what sets the spectrum.
So the question this file exists to answer is S2: are the frequency RATIOS
kernel-invariant?

WHY omega^2 IS CHART-INDEPENDENT, which is what makes any of this a measurement
of the linkage rather than of our parameterisation. At a critical point dV = 0,
so under a reparameterisation s = psi(t) the Hessian transforms as
H -> Dpsi^T H Dpsi (the D^2 psi term is contracted with dV = 0). The mass metric
is a tensor everywhere: M -> Dpsi^T M Dpsi. A congruence applied to BOTH sides of
H v = omega^2 M v leaves the generalised eigenvalues fixed. Measured below under
random basis remixes, not merely argued.

Two things are NOT free choices here and both are reported with every number:
the MASS MODEL (memo C2: point masses k=1/3 vs uniform laminae k=1/12; a
frequency without its mass model is not a reported frequency) and the ZERO-MODE
COUNT, which qvf.2's acceptance criterion demands be reported AND accounted for.
"""
import numpy as np
import scipy.linalg as sla

from jb_a_family import corners
from jb_j_internal_frame import Frame, inertia
from jb_k_hull_hessian import aligned_frame, hull_vol
from jb_o_kernel_family import (KERNELS, Probe, make_V, minimise_on_path,
                                smoothness, vert, A_ICO)
from jb_r_mass_metric import MODELS, metric

np.seterr(all="ignore")   # numpy 1.26 + this BLAS emits spurious matmul
                          # RuntimeWarnings even for np.zeros @ np.zeros;
                          # verified spurious, results are exact.

RAW_KERNELS = [(n, k) for n, k in KERNELS if not n.startswith("spread")]

# Prior values, quoted ONLY so the re-measurement has something to disagree
# with. Nothing downstream reads them.
PRIOR = {
    "thomson_VE": 34.889842063,
    "thomson_ico": 36.554172376,
    "hess_VE": [1.348861, 1.348861, 1.537655, 1.537655, 1.537655, 2.515439],
}


def gen_eig(H, M):
    """(omega^2 ascending, M-orthonormal modes as columns). Symmetric-definite."""
    w, V = sla.eigh(0.5 * (H + H.T), 0.5 * (M + M.T))
    return w, V


def zero_modes(w, rel_tol=1e-6):
    """Count of omega^2 indistinguishable from zero, at an explicit tolerance.

    Relative to the largest |omega^2| in the same spectrum, so the count means
    the same thing for a 1/r^12 kernel (spectrum ~1e-2) and a Gaussian (~1).
    A count reported at one tolerance is not a measurement; the caller sweeps.
    """
    return int(np.sum(np.abs(w) <= rel_tol * max(1.0, np.abs(w).max())))


def degeneracy_blocks(w, rel_tol=1e-6):
    """Group a spectrum into degenerate blocks; returns [(value, multiplicity)]."""
    out = []
    scale = max(1.0, np.abs(w).max())
    for x in w:
        if out and abs(x - out[-1][0]) <= rel_tol * scale:
            out[-1] = (out[-1][0], out[-1][1] + 1)
        else:
            out.append((float(x), 1))
    return out


def fmt_blocks(blocks, prec=6):
    return "[" + ", ".join(f"{v:.{prec}f}" + (f" x{m}" if m > 1 else "")
                           for v, m in blocks) + "]"


# --------------------------------------------------------------------------
# S0 -- the control. If this does not reproduce, everything below is void.
# --------------------------------------------------------------------------

def s0_control():
    print("=" * 78)
    print("S0  CONTROL -- re-measured this session, NOT inherited")
    print("=" * 78)
    V_thom = make_V(RAW_KERNELS[0][1], normalised=False)

    v_ve = V_thom(vert(0.0))
    v_ico = V_thom(vert(A_ICO))
    print(f"  raw Thomson at the VE            {v_ve:.9f}   "
          f"(prior {PRIOR['thomson_VE']:.9f}, diff {abs(v_ve - PRIOR['thomson_VE']):.2e})")
    print(f"  raw Thomson at the icosahedron   {v_ico:.9f}   "
          f"(prior {PRIOR['thomson_ico']:.9f}, diff {abs(v_ico - PRIOR['thomson_ico']):.2e})")
    print(f"  VE is LOWER by {v_ico - v_ve:.9f}  -> the VE is the ground state "
          f"(DECISION 17)")

    a_min, v_min, _, _ = minimise_on_path(lambda a: V_thom(vert(a)))
    print(f"  argmin over the WHOLE circle at 0.1 deg + refinement: a = {a_min:.7f}"
          f"  V = {v_min:.9f}")

    F = aligned_frame(0.0)
    pr = Probe(F)
    print(f"  chart dim at the VE = {F.dim}")

    # TWO OPPOSITE smoothness controls in the same run (jb_o's discipline): if
    # either flips, the discriminator is broken and no verdict here is valid.
    v_th, a_th, r_th = smoothness(pr, V_thom)
    v_hu, a_hu, r_hu = smoothness(pr, hull_vol)
    print(f"  discriminator control A: raw Thomson at the VE -> {v_th} "
          f"(asym {a_th[0]:.3e}/{a_th[1]:.3e}/{a_th[2]:.3e}, ratio {r_th:.3f})")
    print(f"  discriminator control B: Vol_hull   at the VE -> {v_hu} "
          f"(asym {a_hu[0]:.3e}/{a_hu[1]:.3e}/{a_hu[2]:.3e}, ratio {r_hu:.3f})")

    g, gp, gm = pr.grad(V_thom, 1e-4)
    print(f"  |grad V| in all six internal directions = {np.linalg.norm(g):.3e}"
          f"   (transverse {np.linalg.norm(g[1:]):.3e})")

    H = pr.hess(V_thom, 1e-3)
    (p, z, n), ev = inertia(H)
    print(f"  Hessian spectrum  {fmt_blocks(degeneracy_blocks(ev))}")
    print(f"  prior             {fmt_blocks(degeneracy_blocks(np.array(PRIOR['hess_VE'])))}")
    print(f"  max |diff| = {np.abs(np.sort(ev) - np.array(PRIOR['hess_VE'])).max():.2e}")
    print(f"  inertia (pos, zero, neg) = ({p}, {z}, {n})")
    ok = (abs(v_ve - PRIOR["thomson_VE"]) < 1e-8
          and abs(v_ico - PRIOR["thomson_ico"]) < 1e-8
          and (p, z, n) == (6, 0, 0)
          and np.abs(np.sort(ev) - np.array(PRIOR["hess_VE"])).max() < 1e-5)
    print(f"\n  S0 REPRODUCED: {ok}")
    return ok, F, pr


# --------------------------------------------------------------------------
# S1 -- the headline
# --------------------------------------------------------------------------

def s1_spectrum(F, pr, h=1e-3):
    print()
    print("=" * 78)
    print("S1  THE FREQUENCY SPECTRUM AT THE VECTOR EQUILIBRIUM")
    print("=" * 78)
    Ms = {m: metric(F, m) for m in MODELS}
    for m in MODELS:
        ev = np.linalg.eigvalsh(Ms[m])
        print(f"  mass metric M, {m:7s}: {fmt_blocks(degeneracy_blocks(ev), 9)}")
    print("  (M is kernel-INDEPENDENT: it is the kinetic energy, not the potential)")

    results = {}
    for name, kern in RAW_KERNELS:
        V = make_V(kern, normalised=False)
        H = pr.hess(V, h)
        hev = np.linalg.eigvalsh(H)
        print()
        print(f"--- {name} ---")
        print(f"  Hessian eigenvalues      {fmt_blocks(degeneracy_blocks(hev))}")
        for m in MODELS:
            w, Vec = gen_eig(H, Ms[m])
            om = np.sqrt(np.clip(w, 0, None))
            nz = zero_modes(w)
            results[(name, m)] = (w, om, Vec, H, hev)
            lab = " ".join(
                f"{SHORT[mult]}({mult}) omega={np.sqrt(h / mv):.6f}"
                for mult, mv, b in irrep_blocks(Ms[m])
                for h, _ in [block_stiffness(H, b)])
            print(f"  {m:7s} omega^2 {fmt_blocks(degeneracy_blocks(w))}")
            print(f"  {m:7s} omega   {fmt_blocks(degeneracy_blocks(om))}")
            print(f"  {m:7s} by irrep: {lab}")
            print(f"  {m:7s} zero modes at rel tol 1e-6: {nz}   "
                  f"min omega^2 = {w.min():.6e}  "
                  f"(min/max = {w.min() / w.max():.4f})")
    return results, Ms


def s1_zero_mode_accounting(results):
    print()
    print("-" * 78)
    print("ZERO-MODE ACCOUNTING (qvf.2 acceptance criterion)")
    print("-" * 78)
    print("  Tolerance sweep -- a count that moves with the cut is not a count.")
    tols = (1e-4, 1e-6, 1e-8, 1e-10, 1e-12)
    print(f"  {'kernel / model':28s} " + " ".join(f"{t:>7.0e}" for t in tols))
    for (name, m), (w, om, Vec, H, hev) in results.items():
        counts = [zero_modes(w, t) for t in tols]
        print(f"  {name + ' / ' + m:28s} " + " ".join(f"{c:7d}" for c in counts))
    print()
    print("  WHY zero is the expected answer, and where a NON-zero one would come")
    print("  from: the six global rigid motions are removed from the chart by an")
    print("  explicit linear gauge (jb_j), and jb_r R6 measures the leak at 1e-16.")
    print("  A leaked rigid motion has zero stiffness and finite mass and would")
    print("  appear here as omega^2 = 0 -- a BUG in the tangent-space construction,")
    print("  not a physical mode. M is positive definite (jb_r R3), so M cannot")
    print("  manufacture one either. Any zero mode found is therefore a genuine")
    print("  flat direction of V on the variety -- which is what killed soft")
    print("  joints and normalised spread.")


def s1_robustness(F, pr):
    print()
    print("-" * 78)
    print("S1 ROBUSTNESS -- does any of it move? (do not fit a number to one run)")
    print("-" * 78)
    V = make_V(RAW_KERNELS[0][1], normalised=False)
    M = metric(F, "lamina")

    print("  (a) finite-difference step h for the Hessian")
    for h in (3e-3, 1e-3, 3e-4, 1e-4):
        H = pr.hess(V, h)
        w, _ = gen_eig(H, M)
        print(f"    h={h:.0e}  omega = {np.array2string(np.sqrt(w), precision=6)}")

    print("  (b) random remix of the chart basis -- omega^2 must be INVARIANT")
    print("      (H and M are both congruence tensors at a critical point)")
    H = pr.hess(V, 1e-3)
    w0, _ = gen_eig(H, M)
    rng = np.random.default_rng(17)
    for k in range(4):
        C = rng.standard_normal((F.dim, F.dim))
        w, _ = gen_eig(C.T @ H @ C, C.T @ M @ C)
        print(f"    remix {k}: cond {np.linalg.cond(C):8.2f}  "
              f"max |domega^2| = {np.abs(w - w0).max():.3e}")

    print("  (c) an INDEPENDENTLY re-derived chart (unaligned SVD basis)")
    F2 = Frame(0.0)
    pr2 = Probe(F2)
    w2, _ = gen_eig(pr2.hess(V, 1e-3), metric(F2, "lamina"))
    print(f"    max |domega^2| vs aligned chart = {np.abs(w2 - w0).max():.3e}")

    print("  (d) one-ULP-scale jiggle of the anchor angle "
          "(a pinned count has moved on this before)")
    for da in (1e-12, 1e-9, 1e-6, 1e-3):
        Fp = aligned_frame(0.0 + da)
        wp, _ = gen_eig(Probe(Fp).hess(V, 1e-3), metric(Fp, "lamina"))
        print(f"    a0 + {da:.0e}:  max |domega^2| = {np.abs(wp - w0).max():.3e}"
              f"   min omega^2 = {wp.min():.6e}")

    print("  (e) is the POINT-MASS doublet/triplet ORDERING FLIP robust?")
    print("      The flip is a ~1-2% gap, so it is the one qualitative claim in")
    print("      S2 that finite-difference noise could plausibly manufacture.")
    Mp = metric(F, "point")
    bp = {k: v for k, v, _ in irrep_blocks(Mp)}
    for name in ("1/r^2", "1/r^3", "1/r^6", "gauss s=1.0"):
        kern = dict(RAW_KERNELS)[name]
        Vk = make_V(kern, normalised=False)
        row = []
        for h in (3e-3, 1e-3, 3e-4, 1e-4):
            Hk = pr.hess(Vk, h)
            d = {}
            for mult, mval, Q in irrep_blocks(Mp):
                hh, _ = block_stiffness(Hk, Q)
                d[mult] = np.sqrt(hh / mval)
            row.append(d[3] / d[2])
        print(f"    {name:14s} omega_T/omega_D at h=3e-3/1e-3/3e-4/1e-4: "
              + " ".join(f"{r:.6f}" for r in row)
              + f"   spread {max(row) - min(row):.2e}")


# --------------------------------------------------------------------------
# S2 -- what is kernel-invariant
# --------------------------------------------------------------------------

LABEL = {1: "S(singlet)", 2: "D(doublet)", 3: "T(triplet)"}
SHORT = {1: "S", 2: "D", 3: "T"}


def irrep_blocks(M, rel_tol=1e-6):
    """The three symmetry blocks, read off the KERNEL-INDEPENDENT mass metric.

    Labelling blocks by their position in a sorted spectrum is WRONG here --
    the measured ordering flips between mass models and between kernels, and an
    index-based label silently relabels the modes when it does. The mass metric
    is the same object for every kernel, so its eigenspaces are a stable
    labelling: multiplicity 2 = doublet, 3 = triplet, 1 = singlet.

    Returns [(multiplicity, m_value, 6xk basis)].
    """
    ev, EV = np.linalg.eigh(M)
    out = []
    i = 0
    scale = max(1.0, abs(ev).max())
    while i < len(ev):
        j = i + 1
        while j < len(ev) and abs(ev[j] - ev[i]) <= rel_tol * scale:
            j += 1
        out.append((j - i, float(ev[i:j].mean()), EV[:, i:j]))
        i = j
    return out


def block_stiffness(H, basis):
    """(scalar h on the block, max deviation from scalar).

    H restricted to a symmetry block must be a MULTIPLE OF THE IDENTITY if the
    block really is an irreducible representation. The deviation is reported
    because it is exactly what would be non-zero if the 2+3+1 pattern were an
    accident of the kernel rather than a symmetry fact.
    """
    Hb = basis.T @ H @ basis
    h = float(np.trace(Hb) / Hb.shape[0])
    dev = float(np.abs(Hb - h * np.eye(Hb.shape[0])).max())
    return h, dev


def block_table(results, Ms):
    """{(kernel, model): [(mult, omega, h, m, scalar-deviation)]} by irrep."""
    out = {}
    for (name, m), (w, om, Vec, H, hev) in results.items():
        rows = []
        for mult, mval, basis in irrep_blocks(Ms[m]):
            h, dev = block_stiffness(H, basis)
            rows.append((mult, np.sqrt(h / mval), h, mval, dev))
        out[(name, m)] = rows
    return out


def s2_block_structure(results, Ms):
    print()
    print("=" * 78)
    print("S2a  BLOCK STRUCTURE: why the spectrum has only THREE free numbers")
    print("=" * 78)
    print("  The 2+3+1 degeneracy is the chiral tetrahedral symmetry of the VE")
    print("  acting on the 6-D internal tangent space (6 = 1 + 2 + 3). If that is")
    print("  what it is, then H and M are BOTH scalar on each block, and")
    print("  omega^2_block = h_block / m_block exactly.")
    print("  FALSIFIABLE THREE WAYS: (i) H restricted to an M-eigenspace must be")
    print("  a multiple of the identity -- reported as 'dev' below, and it is the")
    print("  H eigenvalues that would have to conspire; (ii) sqrt(h/m) must")
    print("  reproduce the generalised eigenvalue actually solved for; (iii) the")
    print("  multiplicities must match between H and M, which they need not.")
    for m in MODELS:
        blocks = irrep_blocks(Ms[m])
        print(f"\n  mass metric M, {m}: "
              + ", ".join(f"{LABEL[k]} m={v:.9f}" for k, v, _ in blocks))
    tbl = block_table(results, Ms)
    for name, m in (("1/r^1  (Thomson)", "point"), ("1/r^12", "lamina"),
                    ("gauss s=1.0", "point"), ("gauss s=2.5", "lamina")):
        w = np.sort(results[(name, m)][0])
        print(f"\n  {name} / {m}")
        for mult, omega, h, mval, dev in tbl[(name, m)]:
            closest = np.abs(w - omega ** 2).min()
            print(f"    {LABEL[mult]:11s} h={h:13.9f}  m={mval:.9f}  "
                  f"omega^2=h/m={omega ** 2:13.9f}  omega={omega:11.6f}   "
                  f"|H-scalar| on block = {dev:.2e}   "
                  f"dist to solved spectrum = {closest:.2e}")


def s2_ratio_table(results, Ms):
    print()
    print("=" * 78)
    print("S2b  ARE THE FREQUENCY RATIOS KERNEL-INVARIANT?")
    print("=" * 78)
    print("  Modes labelled by IRREP (stable), not by sort position (not stable).")
    tbl = block_table(results, Ms)
    for m in MODELS:
        print(f"\n  MASS MODEL: {m}")
        print(f"  {'kernel':20s} {'omega_D(x2)':>12s} {'omega_T(x3)':>12s} "
              f"{'omega_S':>12s} {'T/D':>9s} {'S/D':>9s} {'S/T':>9s}  ascending")
        for name, _ in RAW_KERNELS:
            d = {mult: om for mult, om, _, _, _ in tbl[(name, m)]}
            order = "<".join(SHORT[k] for k, _ in
                             sorted(((k, v) for k, v in d.items()),
                                    key=lambda kv: kv[1]))
            print(f"  {name:20s} {d[2]:12.6f} {d[3]:12.6f} {d[1]:12.6f} "
                  f"{d[3] / d[2]:9.6f} {d[1] / d[2]:9.6f} {d[1] / d[3]:9.6f}"
                  f"  {order}")
    print()
    print("  A constant ratio column would mean the kernel is a pure stiffness")
    print("  SCALE and the spectrum is a structural fact about the linkage.")
    print("  A moving one means the kernel is a real physical parameter that the")
    print("  owner has to choose.")
    for m in MODELS:
        tb = block_table(results, Ms)
        r = np.array([tb[(n, m)][0][1] for n, _ in RAW_KERNELS])   # placeholder
        sd = {}
        for lab, idx in (("T/D", (3, 2)), ("S/D", (1, 2)), ("S/T", (1, 3))):
            vals = []
            for n, _ in RAW_KERNELS:
                d = {mult: om for mult, om, _, _, _ in tb[(n, m)]}
                vals.append(d[idx[0]] / d[idx[1]])
            sd[lab] = (min(vals), max(vals), max(vals) / min(vals))
        print(f"    {m:7s} spread over the nine kernels: " + "  ".join(
            f"{k} {v[0]:.6f}..{v[1]:.6f} ({v[2]:.3f}x)" for k, v in sd.items()))


def s2_mass_model_effect(results, Ms):
    print()
    print("=" * 78)
    print("S2c  DOES THE MASS MODEL CHANGE THE ORDERING, OR ONLY THE SCALE?")
    print("=" * 78)
    print("  Because H and M are both block-scalar, omega_b = sqrt(h_b/m_b) and")
    print("  the mass model enters ONLY through m_b. So switching models must")
    print("  multiply block b by sqrt(m_b^point / m_b^lamina) -- a factor that is")
    print("  the SAME for every kernel. If the columns below were not constant,")
    print("  the block picture would be wrong.")
    tbl = block_table(results, Ms)
    print(f"  {'kernel':20s} {'f_D':>11s} {'f_T':>11s} {'f_S':>11s}"
          f"   {'point':>9s} {'lamina':>9s}")
    for name, _ in RAW_KERNELS:
        dp = {k: om for k, om, _, _, _ in tbl[(name, "point")]}
        dl = {k: om for k, om, _, _, _ in tbl[(name, "lamina")]}
        op = "<".join(SHORT[k] for k, _ in sorted(dp.items(), key=lambda kv: kv[1]))
        ol = "<".join(SHORT[k] for k, _ in sorted(dl.items(), key=lambda kv: kv[1]))
        print(f"  {name:20s} {dl[2] / dp[2]:11.6f} {dl[3] / dp[3]:11.6f} "
              f"{dl[1] / dp[1]:11.6f}   {op:>9s} {ol:>9s}")
    print("\n  exact factors predicted from the mass metrics alone:")
    bp = {k: v for k, v, _ in irrep_blocks(Ms["point"])}
    bl = {k: v for k, v, _ in irrep_blocks(Ms["lamina"])}
    for k in (2, 3, 1):
        print(f"    {LABEL[k]:11s} sqrt(m_point/m_lamina) = "
              f"sqrt({bp[k]:.9f}/{bl[k]:.9f}) = {np.sqrt(bp[k] / bl[k]):.9f}")


def s2d_crossing(F, pr):
    """Where in the inverse-power family does the point-mass ordering flip?

    The Newton solves are cached per chart coordinate, so sweeping the exponent
    costs almost nothing beyond the first kernel: the stencil is identical and
    only the kernel evaluated on it changes.
    """
    from jb_o_kernel_family import k_invpow
    print()
    print("=" * 78)
    print("S2d  WHERE THE POINT-MASS ORDERING FLIPS, as a function of the exponent")
    print("=" * 78)
    print("  omega_T/omega_D under POINT masses crosses 1 between p=2 and p=3.")
    print("  Locating it turns 'the kernel is a real parameter' from a table")
    print("  observation into a number. Under LAMINAE the same ratio never")
    print("  crosses (reported alongside, as the contrast).")
    Mp, Ml = metric(F, "point"), metric(F, "lamina")
    bp = irrep_blocks(Mp)
    bl = irrep_blocks(Ml)

    def ratio(p, blocks):
        H = pr.hess(make_V(k_invpow(p), normalised=False), 1e-3)
        d = {}
        for mult, mval, Q in blocks:
            h, _ = block_stiffness(H, Q)
            d[mult] = np.sqrt(h / mval)
        return d[3] / d[2]

    print(f"  {'p':>7s} {'point T/D':>12s} {'lamina T/D':>12s}")
    for p in (0.5, 1.0, 1.5, 2.0, 2.2, 2.4, 2.6, 3.0, 4.0, 6.0, 12.0):
        print(f"  {p:7.2f} {ratio(p, bp):12.6f} {ratio(p, bl):12.6f}")
    lo, hi = 2.0, 3.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if (ratio(mid, bp) - 1.0) * (ratio(lo, bp) - 1.0) < 0:
            hi = mid
        else:
            lo = mid
    print(f"  point-mass crossing omega_T = omega_D at p = {0.5 * (lo + hi):.6f}")
    print("  Below it the doublet is the soft mode; above it the triplet is.")


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=False, linewidth=170)
    ok, F, pr = s0_control()
    if not ok:
        raise SystemExit("S0 DID NOT REPRODUCE -- stopping, per the brief.")
    results, Ms = s1_spectrum(F, pr)
    s1_zero_mode_accounting(results)
    s2_block_structure(results, Ms)
    s2_ratio_table(results, Ms)
    s2_mass_model_effect(results, Ms)
    s2d_crossing(F, pr)
    s1_robustness(F, pr)
