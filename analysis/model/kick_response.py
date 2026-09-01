"""kick_response -- drive ONE cell and watch the medium react.

THE QUESTION. Everything measured so far drives the coherent coordinate --
every cell together (the breathe) -- or reads the linearised operator's
spectrum. The owner's ask (2026-09-01): a significantly large patch, ONE
cell driven, the reaction watched. The 2026-08-30 kick (T2 [23764]) did the
unwatched version on a free 5^3 box -- fold impulse, soft joints k = 1,
fold rate 0.2 -- and found the impulse STAYS a fold: 99.7% of the energy in
fold + joint springs, cells translating < 0.3% of a strut and rotating
< 0.1 deg. Its script was never committed. This module is that measurement
committed, gated, scaled up, and exported as frames a page can show.

WHAT THE OPERATOR TURNED OUT TO BE. Linearised at a = -30, the soft-joint
model's fold sector DECOUPLES EXACTLY: the stiffness k C'^T C' and the
mass matrix both have zero fold-to-translation and fold-to-rotation blocks
to machine precision (R5), and the bulk fold row is integer-exact --
m_gamma = 8, onsite 12k, six axis couplings of -2k, row sum zero (the
breathe is the gapless point). So a fold kick is EXACTLY a scalar lattice
wave,

    8 gddot_i = -k (12 g_i - 2 sum_nbr g_j),   i.e.
    gddot = (k/4) (sum_nbr g - 6 g),

with dispersion omega^2 = (k/2)(3 - sum cos) -- gapless at Gamma, sound
speed EXACTLY 1/2 at k = 1, and R5c measures that this scalar band IS one
of dispersion.py's seven to machine precision (it is one member of jb_1c's
doubly degenerate speed-1/2 pair; its partner is translational). The
recorded 99.7% is the finite-amplitude shadow of this exact linear fact:
at linear order the fraction is 100% identically, and the recorded 0.3%
translation / 0.1 deg rotation are the nonlinear corrections at rate 0.2,
which this linear module does not model. "A fold impulse stays a fold" is
not approximately true of the linearised medium; it is its structure.

THE CLOSURE IS THE DOUBLE (OWNER DECISION 22, T2 [23932]): the patch of
record is doubled_block.double(side) -- every boundary half-weld closed
onto the twin sheet. The kick lands mid-sheet-1; when the front reaches
the wall of the box it does not reflect off a free surface and does not
wrap a torus -- it crosses into the twin sheet. The frames export carries
BOTH sheets so the page can show exactly that.

THE FRONT. Arrivals are read at 20% of each cell's own peak (an absolute
threshold reads the exponential tail and "outruns" the band -- measured
before this row was written). The acoustic front disperses (Airy), so the
threshold speed approaches 1/2 FROM BELOW as the patch grows; R7 fits the
far cells of double(16) and gates the window [0.85, 1.02] x 1/2, and R8
doubles the front by quadrupling k (the sqrt k law, statement 2's family).

UNITS. k = 1, plate mass 1 per triangle (m_gamma = 8 is measured, not
assumed). Distance for speeds is counted in primitive lattice steps
(axis-neighbour count), dispersion.bands' own convention, so no physical
length enters. Time as in every spectrum module.

SCOPE. Harmonic, at a = -30, linearised kinematics. Amplitude never
matters in a linear run; the joint-law amplitude family is
joint_exponent.py's line, and the finite-amplitude kick corrections are
open work, not modelled here.
"""
from __future__ import annotations

import base64
import json
import pathlib
import sys

import numpy as np
from scipy import sparse

from analysis.model import assembly as RC
from analysis.model import dispersion as OC
from analysis.model import kinematics as MJ
from analysis.model.double_covering import soft_joint_spectrum as SJ
from analysis.model.su2 import doubled_block as DB

A_REF = MJ.A_REF          # -30, where every published number lives
K_JOINT = SJ.K_JOINT      # 1.0, the convention
RATE = 0.2                # the recorded 2026-08-30 fold-rate kick
DT = 0.01                 # verlet step; omega_max ~ 1.7 so ~370 steps/period
S_FOLD = 0.5              # the fold band's sound speed, gated exact in R5

STRUT = float(np.sqrt(2.0))


# --------------------------------------------------------------------------
# the sparse operator, full seven-coordinate sector
# --------------------------------------------------------------------------

def operator(asm, k=K_JOINT):
    """(Minv blocks (N,7,7), M blocks, C sparse (nc x 7N)) at q0.

    C is constraint_jacobian assembled sparse: row block (weld r, pair m)
    touches only cells k and l. Identical entries to the dense one -- R1
    gates that at machine precision.
    """
    q = asm.q0()
    ctr, R, gam, B = asm.frames(q)
    J = asm.cell_jacobians(ctr, R, B)
    M = asm.mass_blocks(J)
    Minv = np.array([np.linalg.inv(m) for m in M])
    rows, cols, vals = [], [], []
    for r, (kk, ll, pairs) in enumerate(asm.welds):
        for m, (a, b) in enumerate(pairs):
            row = asm._woff[r] + 3 * m
            for d in range(3):
                for c in range(7):
                    rows.append(row + d)
                    cols.append(7 * kk + c)
                    vals.append(J[kk][3 * a + d, c])
                    rows.append(row + d)
                    cols.append(7 * ll + c)
                    vals.append(-J[ll][3 * b + d, c])
    C = sparse.csr_matrix((vals, (rows, cols)), shape=(asm.nc, 7 * asm.N))
    return Minv, M, C


def _minv_apply(Minv, f):
    return np.einsum("nij,nj->ni", Minv, f.reshape(-1, 7)).ravel()


# --------------------------------------------------------------------------
# the fold sector: what the kick actually excites
# --------------------------------------------------------------------------

def fold_operator(M, C, k=K_JOINT):
    """(m_gamma (N,), Hff sparse (N,N)) -- the scalar fold operator, sliced
    out of the full one. R5 gates that the slicing loses nothing: the
    cross blocks it discards are zero to machine precision."""
    n = M.shape[0]
    H = k * (C.T @ C).tocsr()
    fold = np.arange(6, 7 * n, 7)
    return M[:, 6, 6].copy(), H[fold][:, fold].tocsr()


def fold_cross_norms(M, C):
    """(stiffness cross max, mass cross max): the decoupling, measured."""
    n = M.shape[0]
    H = (C.T @ C).tocsr()
    fold = np.arange(6, 7 * n, 7)
    other = np.array([i for i in range(7 * n) if i % 7 != 6])
    hx = H[fold][:, other]
    hmax = float(np.max(np.abs(hx.toarray()))) if hx.nnz else 0.0
    mmax = float(max(np.max(np.abs(M[i][6, :6])) for i in range(n)))
    return hmax, mmax


def verlet_fold(mg, Hff, cell, rate, tmax, dt=DT, sample=0.25):
    """Velocity Verlet on the scalar fold field. `cell` is a single index
    (kicked at `rate`) or a list of (index, rate) pairs. Returns
    (times, fold frames (T, N), energy rows (T, 2): KE, PE)."""
    n = len(mg)
    g = np.zeros(n)
    v = np.zeros(n)
    if isinstance(cell, (list, tuple)):
        for (ci, ri) in cell:
            v[ci] = ri
    else:
        v[cell] = rate
    minv = 1.0 / mg

    def force(gg):
        return -(Hff @ gg)

    a = minv * force(g)
    steps = int(round(tmax / dt))
    every = max(1, int(round(sample / dt)))
    times, frames, erows = [], [], []
    for s in range(steps + 1):
        if s % every == 0:
            times.append(s * dt)
            frames.append(g.copy())
            erows.append((0.5 * float(np.sum(mg * v * v)),
                          0.5 * float(g @ (Hff @ g))))
        if s == steps:
            break
        v += 0.5 * dt * a
        g += dt * v
        a = minv * force(g)
        v += 0.5 * dt * a
    return np.array(times), np.array(frames), np.array(erows)


# --------------------------------------------------------------------------
# full-sector integration (the instrument the reduction is checked against)
# --------------------------------------------------------------------------

def kick_u0(n, cell, rate=RATE):
    u = np.zeros(7 * n)
    u[7 * cell + 6] = rate
    return u


def energies(M, C, k, d, v):
    """(fold KE, translation KE, rotation KE, cross KE, spring PE)."""
    V = v.reshape(-1, 7)
    fold = 0.5 * float(np.sum(M[:, 6, 6] * V[:, 6] ** 2))
    trans = 0.5 * float(np.einsum("ni,nij,nj->", V[:, 0:3], M[:, 0:3, 0:3],
                                  V[:, 0:3]))
    rot = 0.5 * float(np.einsum("ni,nij,nj->", V[:, 3:6], M[:, 3:6, 3:6],
                                V[:, 3:6]))
    total = 0.5 * float(np.einsum("ni,nij,nj->", V, M, V))
    cd = C @ d
    return fold, trans, rot, total - fold - trans - rot, \
        0.5 * k * float(cd @ cd)


def verlet(Minv, M, C, k, u0, tmax, dt=DT, sample=0.25):
    """Velocity Verlet on the full 7N sector from rest, kicked velocity.
    Returns (times, fold frames (T, N), energy rows (T, 5), final (d, v))."""
    n = M.shape[0]
    d = np.zeros(7 * n)
    v = u0.copy()

    def force(dd):
        return -k * (C.T @ (C @ dd))

    a = _minv_apply(Minv, force(d))
    steps = int(round(tmax / dt))
    every = max(1, int(round(sample / dt)))
    times, folds, erows = [], [], []
    for s in range(steps + 1):
        if s % every == 0:
            times.append(s * dt)
            folds.append(d.reshape(-1, 7)[:, 6].copy())
            erows.append(energies(M, C, k, d, v))
        if s == steps:
            break
        v += 0.5 * dt * a
        d += dt * v
        a = _minv_apply(Minv, force(d))
        v += 0.5 * dt * a
    return (np.array(times), np.array(folds), np.array(erows), (d, v))


def modal_folds(asm, k, u0, times):
    """Exact evolution via the dense generalised eigenproblem: what the
    verlet instrument is checked against. Zero modes drift linearly."""
    q = asm.q0()
    ctr, R, gam, B = asm.frames(q)
    J = asm.cell_jacobians(ctr, R, B)
    M = asm.mass_blocks(J)
    n = asm.N
    Mf = np.zeros((7 * n, 7 * n))
    for i in range(n):
        Mf[7 * i:7 * i + 7, 7 * i:7 * i + 7] = M[i]
    L = np.linalg.cholesky(Mf)
    Cd = asm.constraint_jacobian(J)
    H = k * np.dot(Cd.T, Cd)
    A = np.linalg.solve(L, np.linalg.solve(L, H).T).T
    w2, V = np.linalg.eigh((A + A.T) / 2.0)
    w2 = np.clip(w2, 0.0, None)
    y0 = np.dot(L.T, u0)
    c = np.dot(V.T, y0)
    w = np.sqrt(w2)
    nzm = w > 1e-8
    out = []
    for t in times:
        yt = np.empty_like(c)
        yt[nzm] = c[nzm] * np.sin(w[nzm] * t) / w[nzm]
        yt[~nzm] = c[~nzm] * t
        d_full = np.linalg.solve(L.T, np.dot(V, yt))
        out.append(d_full.reshape(-1, 7)[:, 6])
    return np.array(out)


# --------------------------------------------------------------------------
# the doubled block at scale: the fold-sector graph, built directly
# --------------------------------------------------------------------------

def double_graph(side, gc=A_REF):
    """(centres (2N, 3), m_gamma (2N,), Hff sparse (2N, 2N)) for the doubled
    block, built from the lattice graph and the R5b-gated integers -- onsite
    12k, -2k per weld, k = 1 here (scale Hff for other k).

    WHY THIS EXISTS. honeycomb_single discovers welds by an O(N^2) scan over
    sites, which is the model's own careful path and is kept; at side 28 the
    ghost-padded build behind double() would take tens of minutes. The fold
    sector needs none of that: R5 gates that it decouples exactly and R5b
    that its couplings are integers, so the operator IS the graph -- interior
    axis bonds within each sheet, and one -2k tie to the twin cell per
    missing boundary bond (a corner cell ties to its twin three times).
    R9 gates this builder equal to the model's own constructor -- operator
    entries and centres, exactly -- on the sizes both can build, which is
    the only reason it may be trusted at the sizes only this one can.
    """
    sites = DB.boxsites(side)
    idx = {s: i for i, s in enumerate(sites)}
    n = len(sites)
    rows, cols, vals = [], [], []

    def couple(i, j, w):
        rows.append(i)
        cols.append(j)
        vals.append(-2.0 * w)
        rows.append(j)
        cols.append(i)
        vals.append(-2.0 * w)

    axes = [(2, 0, 0), (0, 2, 0), (0, 0, 2)]
    for i, s in enumerate(sites):
        for e in axes:
            t = (s[0] + e[0], s[1] + e[1], s[2] + e[2])
            j = idx.get(t)
            if j is not None:
                couple(i, j, 1.0)          # sheet 1 interior bond
                couple(i + n, j + n, 1.0)  # sheet 2 mirror
        # missing bonds in BOTH directions tie the cell to its twin
        miss = 0
        for e in axes:
            for sgn in (1, -1):
                t = (s[0] + sgn * e[0], s[1] + sgn * e[1], s[2] + sgn * e[2])
                if t not in idx:
                    miss += 1
        if miss:
            couple(i, i + n, float(miss))
    for i in range(2 * n):
        rows.append(i)
        cols.append(i)
        vals.append(12.0)
    Hff = sparse.csr_matrix((vals, (rows, cols)), shape=(2 * n, 2 * n))
    # per-unit-site spacing from the model's own two-cell build, never
    # re-derived: ctr0[1] sits at (2, 0, 0) in site units
    two, _ = RC.honeycomb_single([(0, 0, 0), (2, 0, 0)], gc=gc)
    lu = float(two.ctr0[1][0]) / 2.0
    ctr = np.array(sites, float) * lu
    return np.vstack([ctr, ctr]), np.full(2 * n, 8.0), Hff


# --------------------------------------------------------------------------
# patches, axes, fronts
# --------------------------------------------------------------------------

def centre_cell(side):
    c = side // 2
    return (c * side + c) * side + c


def site_index(side, ix, iy, iz):
    return (ix * side + iy) * side + iz


def mirror_index(side, i):
    """The x-mirror (ix -> side-1-ix) of a cell index, sheet-preserving."""
    n = side * side * side
    sheet, j = divmod(i, n)
    ix, rest = divmod(j, side * side)
    return sheet * n + (side - 1 - ix) * side * side + rest


def axis_cells(side):
    """(distance in primitive steps, sheet-1 cell index) along +x from the
    kicked cell, out to the sheet boundary (the medium continues into the
    twin sheet there, so the boundary cell's arrival is still bulk-like)."""
    c = side // 2
    return [(j, ((c + j) * side + c) * side + c)
            for j in range(1, side - c)]


def arrival_times(times, frames, cells, frac=0.2):
    """First crossing of `frac` of the cell's own peak |fold|. An absolute
    threshold reads the exponential tail and outruns the band; the
    relative one reads the wavefront."""
    out = []
    for (r, i) in cells:
        a = np.abs(frames[:, i])
        peak = float(a.max())
        if peak <= 0.0:
            continue
        hit = int(np.argmax(a > frac * peak))
        if a[hit] > frac * peak:
            out.append((r, float(times[hit])))
    return out


def front_speed(times, frames, side, frac=0.2, tail=4):
    """Slope of r(t) over the LAST `tail` axis arrivals -- the far field,
    where the threshold speed has settled toward the band speed."""
    arr = arrival_times(times, frames, axis_cells(side), frac)
    if len(arr) < tail:
        return float("nan"), arr
    r = np.array([a[0] for a in arr[-tail:]], float)
    t = np.array([a[1] for a in arr[-tail:]], float)
    return float(np.polyfit(t, r, 1)[0]), arr


# --------------------------------------------------------------------------
# frames export for the page
# --------------------------------------------------------------------------

def export(side, tmax, path, sample=0.3):
    """Run the doubled-block kick on the fold sector and write the page's
    frame data: both sheets, int8-quantised per frame, energy split,
    arrivals. Uses the graph-built operator (R9 gates it equal to the
    model's constructor), which is what makes large sides affordable; the
    decoupling (R5) is what licenses the fold-only run."""
    ctr2, mg, Hff = double_graph(side)
    cell = centre_cell(side)
    times, frames, erows = verlet_fold(mg, Hff, cell, RATE, tmax,
                                       sample=sample)
    speed, arr = front_speed(times, frames, side)
    scales = np.max(np.abs(frames), axis=1)
    scales[scales == 0.0] = 1.0
    quant = np.clip(np.round(frames / scales[:, None] * 127.0),
                    -127, 127).astype(np.int8)
    data = {
        "side": side,
        "kicked": cell,
        "rate": RATE,
        "k": K_JOINT,
        "sound": S_FOLD,
        "times": [round(float(t), 4) for t in times],
        "scales": [float(s) for s in scales],
        "centres": base64.b64encode(
            np.asarray(ctr2, np.float32).tobytes()).decode(),
        "folds_i8": base64.b64encode(quant.tobytes()).decode(),
        "energy": [[round(float(x), 9) for x in row] for row in erows],
        "front_speed": None if np.isnan(speed) else round(speed, 4),
        "arrivals": [[int(r), round(t, 3)] for (r, t) in arr],
    }
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data))
    return p, times, frames, erows


def export_pair(side, tmax, path, sample=0.6, x1=None):
    """Frame data for the TWO-KICK page: ONE single-source run, kicked
    off-centre at x-index `x1` (mirror partner at side-1-x1). The page
    superposes g1 +/- mirror(g1) itself -- exactly the two-source field by
    the linearity and mirror symmetry that R10 gates, so both phase
    choices ride one dataset."""
    if x1 is None:
        x1 = max(1, side // 3)
    c = side // 2
    ctr2, mg, Hff = double_graph(side)
    cell = site_index(side, x1, c, c)
    times, frames, erows = verlet_fold(mg, Hff, cell, RATE, tmax,
                                       sample=sample)
    scales = np.max(np.abs(frames), axis=1)
    scales[scales == 0.0] = 1.0
    quant = np.clip(np.round(frames / scales[:, None] * 127.0),
                    -127, 127).astype(np.int8)
    data = {
        "side": side,
        "kicked": cell,
        "kicked2": mirror_index(side, cell),
        "sep": side - 1 - 2 * x1,
        "rate": RATE,
        "k": K_JOINT,
        "sound": S_FOLD,
        "times": [round(float(t), 4) for t in times],
        "scales": [float(s) for s in scales],
        "centres": base64.b64encode(
            np.asarray(ctr2, np.float32).tobytes()).decode(),
        "folds_i8": base64.b64encode(quant.tobytes()).decode(),
    }
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data))
    return p, times, frames


# --------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------

def gate():
    checks, out = [], {}
    A = checks.append

    # R1: the sparse operator IS the dense one
    dbl2, _ = DB.double(2)
    Minv2, M2, C2 = operator(dbl2)
    q = dbl2.q0()
    ctr, R, gam, B = dbl2.frames(q)
    Cd = dbl2.constraint_jacobian(dbl2.cell_jacobians(ctr, R, B))
    r1 = float(np.max(np.abs(C2.toarray() - Cd)))
    out["r1"] = r1
    A(("R1 sparse constraint jacobian == dense, double(2)", r1 < 1e-14,
       f"max abs diff {r1:.2e}"))

    # R2: verlet against exact modal evolution, at second order in dt
    u0 = kick_u0(dbl2.N, centre_cell(2), RATE)
    times, folds, _e2, _ = verlet(Minv2, M2, C2, K_JOINT, u0, 20.0, dt=DT)
    ref = modal_folds(dbl2, K_JOINT, u0, times)
    err1 = float(np.max(np.abs(folds - ref)))
    t2h, f2h, _e2h, _ = verlet(Minv2, M2, C2, K_JOINT, u0, 20.0, dt=DT / 2)
    ref2 = modal_folds(dbl2, K_JOINT, u0, t2h)
    err2 = float(np.max(np.abs(f2h - ref2)))
    order = err1 / err2
    scale = float(np.max(np.abs(ref)))
    out["r2"] = (err1, err2, order, scale)
    A(("R2 verlet == modal on double(2), t <= 20, converging at dt^2",
       err1 < 1e-4 * max(scale, 1e-12) and 3.0 < order < 5.0,
       f"err(dt) {err1:.2e} on amplitude {scale:.2e}, "
       f"err(dt)/err(dt/2) = {order:.2f}"))

    # R3: energy conservation -- verlet's bounded oscillation, no secular
    # drift (the bound is O((omega dt)^2) ~ 3e-4 at this dt; measured 3e-5)
    dbl3, _ = DB.double(3)
    Minv3, M3, C3 = operator(dbl3)
    t3, f3, e3, _ = verlet(Minv3, M3, C3, K_JOINT,
                           kick_u0(dbl3.N, centre_cell(3), RATE), 40.0)
    tot = e3.sum(axis=1)
    r3 = float(np.max(np.abs(tot - tot[0])) / tot[0])
    out["r3"] = r3
    A(("R3 energy conserved on double(3), t <= 40", r3 < 1e-4,
       f"max relative deviation {r3:.2e}"))

    # R4: the double keeps exactly the seven zeros
    ev, _Z, _ = SJ.spectrum(dbl3, k=K_JOINT)
    cut = 1e-8 * float(ev.max())
    nz = int(np.sum(ev < cut))
    out["r4"] = nz
    A(("R4 zero space of double(3) is seven (6 rigid + breathe)", nz == 7,
       f"{nz} zeros"))

    # R5: THE DECOUPLING -- the fold sector separates exactly, and its
    # operator is integer-exact: the kick is a scalar lattice wave
    hmax, mmax = fold_cross_norms(M3, C3)
    mg, Hff = fold_operator(M3, C3)
    i = centre_cell(3)
    row = Hff[i].toarray().ravel()
    nbrs = np.sort(row[np.nonzero(np.abs(row) > 1e-12)])
    ok_row = (abs(row[i] - 12.0) < 1e-12 and len(nbrs) == 7
              and np.all(np.abs(nbrs[:6] + 2.0) < 1e-12)
              and abs(float(row.sum())) < 1e-12
              and np.all(np.abs(mg - 8.0) < 1e-12))
    out["r5"] = (hmax, mmax)
    A(("R5 fold sector decouples EXACTLY: a fold impulse stays a fold at "
       "linear order (T2 [23764]'s 99.7% is its finite-amplitude shadow)",
       hmax < 1e-12 and mmax < 1e-12,
       f"stiffness cross {hmax:.2e}, mass cross {mmax:.2e}"))
    A(("R5b the scalar law is integer-exact: m = 8, onsite 12k, six axis "
       "couplings -2k, row sum 0 (breathe gapless, sound speed 1/2)",
       ok_row, f"onsite {row[i]:.12f}, couplings {nbrs[:6]}"))

    # R5c: the scalar band is one of dispersion.py's seven bands
    cell7 = OC.periodic_cell()
    worst = 0.0
    for kv in ([0.7, 0.0, 0.0], [0.9, 0.5, 0.3], [2.1, 1.1, 0.4]):
        wf = float(np.sqrt((3.0 - np.sum(np.cos(kv))) / 2.0))
        wall = OC.bands(np.array(kv), cell7)
        worst = max(worst, float(np.min(np.abs(wall - wf))))
    out["r5c"] = worst
    A(("R5c omega^2 = (3 - sum cos k)/2 sits inside the seven bands",
       worst < 1e-10, f"worst membership residual {worst:.2e}"))

    # R6: the fold-only integrator IS the full one on the fold sector
    tf, ff, ef = verlet_fold(*fold_operator(M3, C3), centre_cell(3), RATE,
                             40.0)
    r6 = float(np.max(np.abs(ff - f3)))
    out["r6"] = r6
    A(("R6 fold-only verlet == full-sector verlet on double(3)", r6 < 1e-10,
       f"max abs diff {r6:.2e}"))

    # R7: the front approaches the fold sound speed from below, double(16)
    side = 16
    dbl16, _ = DB.double(side)
    Minv16, M16, C16 = operator(dbl16)
    hx, mx = fold_cross_norms(M16, C16)
    mg16, Hff16 = fold_operator(M16, C16)
    t16, g16, _e16 = verlet_fold(mg16, Hff16, centre_cell(side), RATE, 45.0,
                                 sample=0.1)
    vf, arr = front_speed(t16, g16, side)
    out["r7"] = (vf, arr, hx, mx)
    A(("R7 far-field front speed on double(16) vs the exact 1/2",
       (not np.isnan(vf)) and 0.85 * S_FOLD <= vf <= 1.02 * S_FOLD,
       f"front {vf:.4f} steps/time over {len(arr)} arrivals "
       f"(window [{0.85 * S_FOLD:.3f}, {1.02 * S_FOLD:.3f}])"))

    # R8: 4x stiffness doubles the front -- the sqrt(k) law, and the
    # detector detects
    t16b, g16b, _e = verlet_fold(mg16, (4.0 * Hff16).tocsr(),
                                 centre_cell(side), RATE, 23.0,
                                 dt=DT / 2, sample=0.05)
    vf4, _ = front_speed(t16b, g16b, side)
    ratio = vf4 / vf
    out["r8"] = (vf4, ratio)
    A(("R8 4x stiffness doubles the front (sqrt k)", abs(ratio - 2.0) < 0.1,
       f"front(4k)/front(k) = {ratio:.3f}"))

    # R9: the graph builder IS the model's own constructor, on the sizes
    # both can build -- the only reason double_graph may be trusted at the
    # sizes only it can build
    worst_h, worst_v = 0.0, 0.0
    for (dblg, Mg, Cg) in ((dbl2, M2, C2), (dbl3, M3, C3)):
        s9 = round((dblg.N / 2) ** (1 / 3))
        mg9, H9 = fold_operator(Mg, Cg)
        ctrg, mgg, Hg = double_graph(s9)
        worst_h = max(worst_h, float(np.max(np.abs((Hg - H9).toarray()))))
        worst_v = max(worst_v, float(np.max(np.abs(ctrg - dblg.ctr0))))
        worst_v = max(worst_v, float(np.max(np.abs(mgg - mg9))))
    out["r9"] = (worst_h, worst_v)
    A(("R9 graph-built double == the model's own constructor (sides 2, 3)",
       worst_h < 1e-12 and worst_v < 1e-12,
       f"operator diff {worst_h:.2e}, centres/mass diff {worst_v:.2e}"))

    # R10: two kicks -- superposition and the x-mirror, measured against a
    # direct two-source integration. What licenses the two-kick page to
    # superpose g1 +/- mirror(g1) from one exported run.
    sideT = 8
    ctrT, mgT, HffT = double_graph(sideT)
    cT = sideT // 2
    i1 = site_index(sideT, 2, cT, cT)
    i2 = site_index(sideT, sideT - 3, cT, cT)
    _t, gA, _e = verlet_fold(mgT, HffT, i1, RATE, 16.0, sample=0.2)
    _t, gB, _e = verlet_fold(mgT, HffT, i2, RATE, 16.0, sample=0.2)
    _t, gD, _e = verlet_fold(mgT, HffT, [(i1, RATE), (i2, -RATE)], 0.0,
                             16.0, sample=0.2)
    sup = float(np.max(np.abs(gD - (gA - gB))))
    mmap = np.array([mirror_index(sideT, i) for i in range(2 * sideT ** 3)])
    mir = float(np.max(np.abs(gB - gA[:, mmap])))
    out["r10"] = (sup, mir)
    A(("R10 two kicks: direct two-source run == g1 - g2 (superposition), "
       "and g2 == mirror(g1) (lattice symmetry), double(8)",
       sup < 1e-10 and mir < 1e-10,
       f"superposition {sup:.2e}, mirror {mir:.2e}"))

    return checks, out


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "export":
        side = int(sys.argv[2]) if len(sys.argv) > 2 else 28
        tmax = float(sys.argv[3]) if len(sys.argv) > 3 else 72.0
        sample = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
        p, times, frames, erows = export(
            side, tmax, "analysis/.pages/data/kick.json", sample=sample)
        print(f"exported double({side}) kick, {len(times)} frames -> {p}")
        return 0
    if len(sys.argv) > 1 and sys.argv[1] == "export-pair":
        side = int(sys.argv[2]) if len(sys.argv) > 2 else 28
        tmax = float(sys.argv[3]) if len(sys.argv) > 3 else 72.0
        sample = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
        x1 = int(sys.argv[5]) if len(sys.argv) > 5 else 9
        p, times, frames = export_pair(
            side, tmax, "analysis/.pages/data/kick_pair.json",
            sample=sample, x1=x1)
        print(f"exported double({side}) pair source (x1={x1}), "
              f"{len(times)} frames -> {p}")
        return 0
    checks, _ = gate()
    fails = 0
    for name, ok, detail in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}  [{detail}]")
        fails += 0 if ok else 1
    print(f"{len(checks) - fails}/{len(checks)} rows pass")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
