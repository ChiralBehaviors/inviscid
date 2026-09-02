"""shear_response -- shove ONE cell sideways, and stir one: polarized waves.

THE QUESTION (owner, 2026-09-01): can this medium carry POLARIZED waves?
kick_response.py drives the fold and finds a scalar wave at speed 1/2 --
and nothing else, because the fold sector decouples exactly (its R5). The
other acoustic branches are TRANSLATIONAL, and a translational wave has a
direction of displacement: that is what polarization is. This module
drives them.

WHAT THE ACOUSTIC SECTOR IS, BY EIGENVECTOR (R2, read from dispersion.py's
seven bands at small k, along [100], [110], [111] and random directions):

    two lowest    speed 1/(2 sqrt 2)   pure translation, displacement
                                        PERPENDICULAR to k: TRANSVERSE,
                                        doubly degenerate
    next two      speed 1/2            one pure fold (kick_response's
                                        wave) and one pure translation
                                        ALONG k: LONGITUDINAL
    top three     gapped at 12/5       pure rotation

and the speeds are the SAME in every direction to 1e-4 at |k| = 0.02: at
long wavelength the medium is an ISOTROPIC elastic solid with
c_L / c_T = sqrt 2 exactly -- Poisson ratio zero. Not what six axis bonds
on a cubic lattice usually give; the welds carry the whole vertex
displacement, so shear is as stiff as it is. So: YES. Any transverse
polarization propagates unchanged (no birefringence at long wavelength),
and the two degenerate shear modes in quadrature are a CIRCULARLY
polarized wave that keeps its handedness. The geometry is chiral; at
linear order, long wavelength, nothing splits left from right.

AT FINITE WAVELENGTH THE CHIRALITY SHOWS -- AS BIREFRINGENCE, NOT ACTIVITY
(R2b, R2c). The cell at a = -30 has only the tetrahedral rotations, so a
lattice axis is two-fold, not four-fold, and nothing forces the two
transverse modes along it to stay degenerate: they split, 1.3% at k = 1
(wavelength 6.3 cells), 3.4% at k = 2. Their eigenmodes are LINEAR
(ellipticity 1.000 at every k tried), so the split is linear
birefringence -- a circular wave along an axis beats between them with a
length pi c_T / (dw) -- and NOT acoustic activity (which would make the
eigenmodes circular and rotate a linear polarization). At the stir's
wavelength of four cells the beat length is 75 cells, six times the
distance from the stirred cell to the wall, so on the exported box the
circle is clean (R6: rotation rate 1.6% under the drive's). One branch's group velocity also
peaks 2% above c_T at k ~ 1.25 (mild anomalous dispersion; a transverse
front may run up to it, R5's ceiling).

TWO DRIVES, TWO PAGES' WORTH OF FRAMES.
  SHOVE: a translational velocity impulse on the centre cell along one
    axis. Linear polarization. The field is a dipole: longitudinal lobes
    along the shove axis at speed 1/2, transverse lobes across it at
    1/(2 sqrt 2), displacement along the shove axis everywhere. Both
    speeds in one picture (R5; the export shoves along x so the page's
    top-down cutaway holds both).
  STIR: a rotating force on the centre cell, (cos w t, sin w t, 0) under
    a raised-cosine envelope for a few periods, w chosen on the
    transverse band at wavelength four cells. Circular polarization. Along
    +-z the wave is transverse and its displacement vector ROTATES at w in
    the drive's sense -- the same lab sense on both sides, so the two
    outgoing waves have OPPOSITE helicity about their own propagation
    directions (R6).

THE OPERATOR AT SCALE. kick_response.double_graph builds the scalar fold
operator from the lattice graph. This module's double_graph7 builds the
FULL seven-coordinate operator the same way: every cell is the same body
in the same orientation, so ONE cell's 36x7 vertex Jacobian and its
diagonal 7x7 mass (8, 8, 8, 80/9 x3, 8 -- measured, not assumed) serve
every cell; an axis bond contributes k (J_a - J_b)^T (J_a - J_b) per
coincident vertex pair; a missing boundary bond ties the cell's own bond
vertices to the same vertices of its twin, exactly as doubled_block.double
closes the half-weld. R1 gates this builder equal to the model's own
constructor -- stiffness and mass, exactly -- on the sizes both can build.

UNITS as kick_response: k = 1, plate mass 1 per triangle, distance in
primitive lattice steps, time as in every spectrum module.

SCOPE. Harmonic, at a = -30, linearised. The fold stays EXACTLY zero under
both drives (R3: the converse of kick_response R5). Rotation IS excited
(R4): translation and rotation couple at k != 0 -- a shear gradient is a
torque -- so the integrated sector is all seven coordinates, and the
frames carry the translation.
"""
from __future__ import annotations

import base64
import json
import pathlib
import sys

import numpy as np
from scipy import sparse

from analysis.model import dispersion as OC
from analysis.model import kick_response as KR
from analysis.model.su2 import doubled_block as DB

A_REF = KR.A_REF
K_JOINT = KR.K_JOINT
DT = KR.DT
RATE = KR.RATE            # shove: translational velocity impulse, one cell
S_LONG = 0.5
S_TRANS = 1.0 / (2.0 * np.sqrt(2.0))
STIR_WAVELENGTH = 4.0     # cells, sets the stir frequency on the band
STIR_PERIODS = 2.0


# --------------------------------------------------------------------------
# the cell, once; the operator, from the graph
# --------------------------------------------------------------------------

def cell_blocks():
    """(J (36, 7), M (7, 7), bonds) of ONE cell from dispersion.periodic_cell:
    the vertex Jacobian, the mass block, and the three +axis bonds with
    their coincident vertex pairs (a on this cell, b on the neighbour)."""
    J, M, bonds = OC.periodic_cell()
    return J, M, bonds


def double_graph7(side, k=K_JOINT):
    """(centres (2N, 3), M (7, 7), H sparse (14N, 14N)) for the doubled
    block: the full seven-coordinate stiffness k C^T C, assembled from the
    lattice graph and one cell's blocks. R1 gates it equal to
    kick_response.operator(doubled_block.double(side)) exactly."""
    J, M, bonds = cell_blocks()
    sites = DB.boxsites(side)
    idx = {s: i for i, s in enumerate(sites)}
    n = len(sites)
    blocks = {}   # (i, j) -> 7x7 accumulated

    def add(i, j, blk):
        key = (i, j)
        if key in blocks:
            blocks[key] = blocks[key] + blk
        else:
            blocks[key] = blk.copy()

    def bond(i, j, Ja, Jb):
        # k (Ja d_i - Jb d_j)^2
        add(i, i, Ja.T @ Ja)
        add(j, j, Jb.T @ Jb)
        add(i, j, -(Ja.T @ Jb))
        add(j, i, -(Jb.T @ Ja))

    for i, s in enumerate(sites):
        for (e, pairs) in bonds:
            step = (2 * e[0], 2 * e[1], 2 * e[2])
            for sgn in (1, -1):
                t = (s[0] + sgn * step[0], s[1] + sgn * step[1],
                     s[2] + sgn * step[2])
                j = idx.get(t)
                if sgn == 1 and j is not None:
                    for (a, b) in pairs:
                        Ja, Jb = J[3 * a:3 * a + 3], J[3 * b:3 * b + 3]
                        bond(i, j, Ja, Jb)
                        bond(i + n, j + n, Ja, Jb)
                elif j is None:
                    # missing bond: this cell's own vertices on that face
                    # (a's on the +e face, b's on the -e face) tie to the
                    # same vertices of the twin
                    for (a, b) in pairs:
                        v = a if sgn == 1 else b
                        Jv = J[3 * v:3 * v + 3]
                        bond(i, i + n, Jv, Jv)
    rows, cols, vals = [], [], []
    for (i, j), blk in blocks.items():
        for r in range(7):
            for c in range(7):
                if blk[r, c] != 0.0:
                    rows.append(7 * i + r)
                    cols.append(7 * j + c)
                    vals.append(k * blk[r, c])
    H = sparse.csr_matrix((vals, (rows, cols)), shape=(14 * n, 14 * n))
    H.sum_duplicates()
    ctr2, _mg, _Hff = KR.double_graph(side)
    return ctr2, M, H


# --------------------------------------------------------------------------
# the integrator: seven coordinates, sparse, with an optional drive
# --------------------------------------------------------------------------

def stir_omega(k=K_JOINT, wavelength=STIR_WAVELENGTH):
    """The transverse band's frequency at |k| = 2 pi / wavelength along an
    axis, read from dispersion.bands (the lowest branch), scaled sqrt k."""
    cell = OC.periodic_cell()
    kk = 2.0 * np.pi / wavelength
    w = np.sort(OC.bands(np.array([0.0, 0.0, kk]), cell))
    return float(w[0]) * float(np.sqrt(k))


def stir_force(cell, omega, amp=1.0, periods=STIR_PERIODS):
    """drive(t) -> (cell, 3-force): rotating in the xy-plane, raised-cosine
    envelope over `periods` periods, zero after."""
    tdrive = periods * 2.0 * np.pi / abs(omega)

    def drive(t):
        if t >= tdrive:
            return None
        env = np.sin(np.pi * t / tdrive) ** 2
        return cell, amp * env * np.array([np.cos(omega * t),
                                           np.sin(omega * t), 0.0])
    drive.tdrive = tdrive
    return drive


def energies7(M, H, d, v):
    """(trans KE, rot KE, fold KE, spring PE)."""
    V = v.reshape(-1, 7)
    m = np.diag(M)
    ke = 0.5 * (V * V) * m[None, :]
    return (float(ke[:, 0:3].sum()), float(ke[:, 3:6].sum()),
            float(ke[:, 6].sum()), 0.5 * float(np.sum(d * (H @ d))))


def verlet7(M, H, u0, tmax, dt=DT, sample=0.25, drive=None):
    """Velocity Verlet on the full seven-coordinate sector from rest with
    velocity u0 and optional drive(t) -> (cell, 3-force) or None.
    Returns (times, translation frames (T, Ncells, 3), fold frames (T, N),
    energy rows (T, 4), final (d, v))."""
    ncell = H.shape[0] // 7
    minv = np.tile(1.0 / np.diag(M), ncell)
    d = np.zeros(7 * ncell)
    v = u0.copy()

    def force(dd, t):
        f = -(H @ dd)
        if drive is not None:
            dr = drive(t)
            if dr is not None:
                c, fv = dr
                f[7 * c:7 * c + 3] += fv
        return f

    a = minv * force(d, 0.0)
    steps = int(round(tmax / dt))
    every = max(1, int(round(sample / dt)))
    times, trans, folds, erows = [], [], [], []
    for s in range(steps + 1):
        if s % every == 0:
            times.append(s * dt)
            D7 = d.reshape(-1, 7)
            trans.append(D7[:, 0:3].copy())
            folds.append(D7[:, 6].copy())
            erows.append(energies7(M, H, d, v))
        if s == steps:
            break
        v += 0.5 * dt * a
        d += dt * v
        a = minv * force(d, (s + 1) * dt)
        v += 0.5 * dt * a
    return (np.array(times), np.array(trans), np.array(folds),
            np.array(erows), (d, v))


def shove_u0(ncell, cell, axis=1, rate=RATE):
    u = np.zeros(7 * ncell)
    u[7 * cell + axis] = rate
    return u


# --------------------------------------------------------------------------
# reading the field: polarization and fronts
# --------------------------------------------------------------------------

def polarization_at_k(kvec, cell=None):
    """(speeds sorted (4,), sector weights (7, 3), longitudinal fraction (7,))
    of the seven bands at kvec: which are translation / rotation / fold,
    and how much of each translational eigenvector lies along k."""
    if cell is None:
        cell = OC.periodic_cell()
    J, M, bonds = cell
    rows = []
    for e, prs in bonds:
        ph = np.exp(1j * float(np.dot(kvec, e)))
        for (a, b) in prs:
            rows.append(J[3 * a:3 * a + 3] - ph * J[3 * b:3 * b + 3])
    C = np.vstack(rows)
    Hk = K_JOINT * (C.conj().T @ C)
    Hk = (Hk + Hk.conj().T) / 2.0
    L = np.linalg.cholesky(M)
    A = np.linalg.solve(L, np.linalg.solve(L, Hk).conj().T).conj().T
    w2, Y = np.linalg.eigh((A + A.conj().T) / 2.0)
    V = np.linalg.solve(L.T, Y)
    w = np.sqrt(np.clip(w2, 0.0, None))
    kk = float(np.linalg.norm(kvec))
    khat = np.asarray(kvec, float) / kk
    m = np.diag(M)
    weights = np.zeros((7, 3))
    longf = np.zeros(7)
    for i in range(7):
        v = V[:, i]
        e = np.abs(v) ** 2 * m
        tot = float(e.sum())
        weights[i] = [e[0:3].sum() / tot, e[3:6].sum() / tot, e[6] / tot]
        tt = v[0:3]
        nt = float(np.vdot(tt, tt).real)
        longf[i] = abs(np.vdot(khat, tt)) ** 2 / nt if nt > 1e-20 else 0.0
    return w / kk, weights, longf


def group_velocity_max(band, axis=(1, 0, 0), npts=400):
    """max_k d omega/dk of acoustic branch `band` (0, 1 transverse; 2, 3
    longitudinal + fold) along `axis`, numerically from dispersion.bands."""
    cell = OC.periodic_cell()
    d = np.array(axis, float)
    d /= np.linalg.norm(d)
    ks = np.linspace(0.01, np.pi, npts)
    W = np.array([np.sort(OC.bands(k * d, cell))[:4] for k in ks])
    vg = np.gradient(W[:, band], ks)
    i = int(np.argmax(vg))
    return float(vg[i]), float(ks[i])


def axis_split(kmag, cell=None):
    """(w_lo, w_hi, ellipticity_lo, ellipticity_hi) of the transverse pair
    at k = kmag along [100]: the split, and |t.t|/|t|^2 of each transverse
    eigenvector (1 = linear, 0 = circular)."""
    if cell is None:
        cell = OC.periodic_cell()
    J, M, bonds = cell
    kvec = np.array([kmag, 0.0, 0.0])
    rows = []
    for e, prs in bonds:
        ph = np.exp(1j * float(np.dot(kvec, e)))
        for (a, b) in prs:
            rows.append(J[3 * a:3 * a + 3] - ph * J[3 * b:3 * b + 3])
    C = np.vstack(rows)
    Hk = K_JOINT * (C.conj().T @ C)
    Hk = (Hk + Hk.conj().T) / 2.0
    L = np.linalg.cholesky(M)
    A = np.linalg.solve(L, np.linalg.solve(L, Hk).conj().T).conj().T
    w2, Y = np.linalg.eigh((A + A.conj().T) / 2.0)
    V = np.linalg.solve(L.T, Y)
    w = np.sqrt(np.clip(w2, 0.0, None))
    o = np.argsort(w)
    ell = []
    for i in o[:2]:
        t = V[1:3, i]
        ell.append(float(abs(t @ t) / np.vdot(t, t).real))
    return (float(w[o[0]]), float(w[o[1]]), ell[0], ell[1],
            float(abs(w[o[2]] - w[o[3]])))


def line_cells(side, axis, sgn, skip_wall=True):
    """(distance, sheet-1 index) from the centre cell along +-axis. The
    wall cell -- tied to its twin, where the field is the sum of both
    sheets' -- is left out of front fits: measured to cross the 20%
    threshold 0.1 early at double(24)."""
    c = side // 2
    jmax = (side - 1 - c) if sgn > 0 else c
    if skip_wall:
        jmax -= 1
    out = []
    for j in range(1, jmax + 1):
        p = [c, c, c]
        p[axis] = c + sgn * j
        out.append((j, KR.site_index(side, *p)))
    return out


def line_front(times, amp, side, axis, sgn=1, frac=0.2, tail=3, tcut=None):
    """Far-field front speed along one axis from the scalar field `amp`
    (T, N): slope of r(t) over the last `tail` arrivals, read before
    `tcut` (default 1.5 side): the double is closed, and what crosses into
    the twin sheet comes back -- a cell's own peak read over a long run
    is that return, and the 20%-of-peak threshold then reads late. The
    threshold on a dispersing pulse also reads a precursor on the nearer
    cells (measured: r = 7 of double(24) crosses 0.6 early on the
    longitudinal axis), so the fit is over the three farthest cells."""
    if tcut is None:
        tcut = 1.5 * side
    sel = times <= tcut
    times, amp = times[sel], amp[sel]
    arr = KR.arrival_times(times, amp, line_cells(side, axis, sgn), frac)
    if len(arr) < tail:
        return float("nan"), arr
    r = np.array([a[0] for a in arr[-tail:]], float)
    t = np.array([a[1] for a in arr[-tail:]], float)
    return float(np.polyfit(t, r, 1)[0]), arr


def rotation_rate(times, trans, cell, t0, t1):
    """(dphi/dt, transverse fraction) of the displacement vector of `cell`
    over [t0, t1]: phi = atan2(d_y, d_x), unwrapped, least-squares slope."""
    sel = (times >= t0) & (times <= t1)
    d = trans[sel, cell, :]
    phi = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
    slope = float(np.polyfit(times[sel], phi, 1)[0])
    xy = float(np.sum(d[:, 0] ** 2 + d[:, 1] ** 2))
    zz = float(np.sum(d[:, 2] ** 2))
    return slope, xy / (xy + zz)


# --------------------------------------------------------------------------
# frames export for the page
# --------------------------------------------------------------------------

def _quantise(trans):
    """Per-frame int8 quantisation of (T, N, 3) by the frame's max |d|, after
    subtracting each frame's mean displacement over all cells: the closed
    double conserves momentum, so a shove sets the whole block drifting
    (a rigid zero mode, 8 x rate shared by 2N cells) and by t ~ 40 that
    drift would dominate the per-frame normalisation of the decayed field."""
    trans = trans - trans.mean(axis=1, keepdims=True)
    scales = np.max(np.linalg.norm(trans, axis=2), axis=1)
    scales[scales == 0.0] = 1.0
    q = np.clip(np.round(trans / scales[:, None, None] * 127.0),
                -127, 127).astype(np.int8)
    return scales, q


def export(side, tmax, path, sample=1.0):
    """Both drives on double(side), one JSON: translation frames int8x3
    per cell per frame, both sheets, plus energies and fronts."""
    ctr2, M, H = double_graph7(side)
    ncell = H.shape[0] // 7
    c = KR.centre_cell(side)

    # the shove runs fine-sampled so the fronts resolve; the page's frames
    # are every `sample` of that run
    fine = 0.1
    every = max(1, int(round(sample / fine)))
    t_f, tr_f, _f, e_f, _ = verlet7(M, H, shove_u0(ncell, c, axis=0), tmax,
                                    sample=fine)
    amp = np.linalg.norm(tr_f, axis=2)
    vL, arrL = line_front(t_f, amp, side, 0)
    vT, arrT = line_front(t_f, amp, side, 2)
    t_s, tr_s, e_s = t_f[::every], tr_f[::every], e_f[::every]
    del tr_f, amp
    sc_s, q_s = _quantise(tr_s)

    w0 = stir_omega()
    drive = stir_force(c, w0)
    t_r, tr_r, _f, e_r, _ = verlet7(M, H, np.zeros(7 * ncell), tmax,
                                    sample=sample, drive=drive)
    sc_r, q_r = _quantise(tr_r)

    data = {
        "side": side, "kicked": c, "rate": RATE, "k": K_JOINT,
        "sound_long": S_LONG, "sound_trans": S_TRANS,
        "centres": base64.b64encode(
            np.asarray(ctr2, np.float32).tobytes()).decode(),
        "shove": {
            "axis": 0,
            "times": [round(float(t), 4) for t in t_s],
            "scales": [float(s) for s in sc_s],
            "disp_i8": base64.b64encode(q_s.tobytes()).decode(),
            "energy": [[round(float(x), 9) for x in row] for row in e_s],
            "front_long": None if np.isnan(vL) else round(vL, 4),
            "front_trans": None if np.isnan(vT) else round(vT, 4),
            "arrivals_long": [[int(r), round(t, 3)] for (r, t) in arrL],
            "arrivals_trans": [[int(r), round(t, 3)] for (r, t) in arrT],
        },
        "stir": {
            "omega": w0, "tdrive": drive.tdrive,
            "wavelength": STIR_WAVELENGTH, "periods": STIR_PERIODS,
            "times": [round(float(t), 4) for t in t_r],
            "scales": [float(s) for s in sc_r],
            "disp_i8": base64.b64encode(q_r.tobytes()).decode(),
            "energy": [[round(float(x), 9) for x in row] for row in e_r],
        },
    }
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data))
    return p, data


# --------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------

def gate():
    checks, out = [], {}
    A = checks.append

    # R1: the graph-built seven-coordinate operator IS the model's own
    worst_h, worst_m, worst_c = 0.0, 0.0, 0.0
    small = {}
    for s in (2, 3):
        dbl, _ = DB.double(s)
        Minv, Mb, C = KR.operator(dbl)
        Hm = (K_JOINT * (C.T @ C)).tocsr()
        ctrg, Mg, Hg = double_graph7(s)
        worst_h = max(worst_h, float(np.max(np.abs((Hg - Hm).toarray()))))
        worst_m = max(worst_m, float(np.max(np.abs(Mb - Mg[None, :, :]))))
        worst_c = max(worst_c, float(np.max(np.abs(ctrg - dbl.ctr0))))
        small[s] = (Mg, Hg)
    out["r1"] = (worst_h, worst_m, worst_c)
    A(("R1 graph-built 7-coordinate double == the model's own constructor "
       "(sides 2, 3): stiffness, mass, centres",
       worst_h < 1e-12 and worst_m < 1e-12 and worst_c < 1e-12,
       f"stiffness diff {worst_h:.2e}, mass diff {worst_m:.2e}, "
       f"centres diff {worst_c:.2e}"))

    # R2: the acoustic sector by eigenvector -- transverse pair, longitudinal
    # + fold, gapped rotation; isotropic at long wavelength
    cell7 = OC.periodic_cell()
    rng = np.random.default_rng(1)
    dirs = [np.array(v, float) for v in ((1, 0, 0), (1, 1, 0), (1, 1, 1))]
    dirs += [rng.normal(size=3) for _ in range(40)]
    tmin, tmax_, lmin, lmax = 1e9, -1e9, 1e9, -1e9
    ok_id = True
    for dv in dirs:
        dv = dv / np.linalg.norm(dv)
        sp, wt, lf = polarization_at_k(0.02 * dv, cell7)
        o = np.argsort(sp)
        sp, wt, lf = sp[o], wt[o], lf[o]
        tmin, tmax_ = min(tmin, sp[0], sp[1]), max(tmax_, sp[0], sp[1])
        lmin, lmax = min(lmin, sp[2], sp[3]), max(lmax, sp[2], sp[3])
        # two lowest: translation, transverse
        ok_id &= bool(np.all(wt[0:2, 0] > 0.999) and np.all(lf[0:2] < 1e-3))
        # next two span fold + longitudinal translation (degenerate pair:
        # any basis of the plane): summed weights 1 fold + 1 trans, and the
        # translational part is along k
        ok_id &= abs(wt[2:4, 2].sum() - 1.0) < 1e-6
        ok_id &= abs(wt[2:4, 0].sum() - 1.0) < 1e-6
        ok_id &= bool(np.all((lf[2:4] > 0.999) | (wt[2:4, 0] < 1e-6)))
        # top three: rotation, gapped
        ok_id &= bool(np.all(wt[4:7, 1] > 0.999) and sp[4] * 0.02 > 1.0)
    vg0, k0 = group_velocity_max(0)
    vg1, k1 = group_velocity_max(1)
    vgT, kT = max((vg0, k0), (vg1, k1))
    vgL, kL = group_velocity_max(2)
    out["r2"] = (tmin, tmax_, lmin, lmax, vgT, kT, vgL)
    A(("R2 the acoustic sector by eigenvector: two TRANSVERSE translational "
       "branches at 1/(2 sqrt 2), LONGITUDINAL translation degenerate with "
       "the fold at 1/2, rotation gapped -- and isotropic (43 directions)",
       ok_id and abs(tmin - S_TRANS) < 1e-4 and abs(tmax_ - S_TRANS) < 1e-4
       and abs(lmin - S_LONG) < 1e-4 and abs(lmax - S_LONG) < 1e-4,
       f"transverse [{tmin:.5f}, {tmax_:.5f}] vs {S_TRANS:.5f}; "
       f"longitudinal [{lmin:.5f}, {lmax:.5f}] vs {S_LONG:.5f}"))

    A(("R2b group velocity along an axis: the transverse branch peaks ABOVE "
       "its sound speed at finite k (mild anomalous dispersion, so a "
       "transverse front may run up to it); the longitudinal never exceeds 1/2",
       S_TRANS < vgT < 1.05 * S_TRANS and abs(vgL - S_LONG) < 1e-3,
       f"transverse max dw/dk {vgT:.4f} at k {kT:.2f} (c_T {S_TRANS:.4f}); "
       f"longitudinal max {vgL:.4f}"))

    # R2c: along an axis the transverse pair SPLITS at finite k -- the cell
    # is chiral, the axis only two-fold -- and the eigenmodes are LINEAR:
    # birefringence, not acoustic activity. At the stir's wavelength the
    # beat length dwarfs the box, which is why R6's circle survives.
    wl1, wh1, e1, e2, dg1 = axis_split(1.0)
    kS = 2.0 * np.pi / STIR_WAVELENGTH
    wlS, whS, _e3, _e4, dgS = axis_split(kS)
    dgap = max(dg1, dgS)
    split1 = (wh1 - wl1) / wl1
    beat = np.pi * S_TRANS / (whS - wlS)   # cells: pi / dk, dk ~ dw / c_T
    out["r2c"] = (split1, e1, e2, beat, dgap)
    A(("R2c along [100] the transverse pair splits at finite k with LINEAR "
       "eigenmodes (birefringent, not acoustically active); at the stir's "
       "k = pi/2 the linear-birefringence beat length dwarfs the exported "
       "box; the longitudinal branch stays degenerate with the fold there",
       0.005 < split1 < 0.05 and e1 > 0.999 and e2 > 0.999 and beat > 2 * 24
       and dgap < 1e-10,
       f"split {100 * split1:.2f}% at k = 1; ellipticity {e1:.4f}, {e2:.4f}; "
       f"beat length {beat:.0f} cells at k = pi/2; fold-longitudinal gap "
       f"{dgap:.1e}"))

    # R3: the fold stays EXACTLY zero under a shove -- the converse of
    # kick_response R5 -- and R4: rotation IS excited; energy conserved
    M3, H3 = small[3]
    n3 = H3.shape[0] // 7
    t3, tr3, fo3, e3, _ = verlet7(M3, H3, shove_u0(n3, KR.centre_cell(3)),
                                  40.0)
    fmax = float(np.max(np.abs(fo3)))
    tot = e3.sum(axis=1)
    drift = float(np.max(np.abs(tot - tot[0])) / tot[0])
    rot_frac = float(np.max(e3[:, 1]) / tot[0])
    out["r3"] = (fmax, drift, rot_frac)
    A(("R3 a shove never folds: fold coordinate zero to 1e-12 on double(3), "
       "t <= 40; energy conserved", fmax < 1e-12 and drift < 1e-4,
       f"max |fold| {fmax:.1e}, energy drift {drift:.2e}"))
    A(("R4 a shove DOES rotate cells: translation and rotation couple at "
       "k != 0 (a shear gradient is a torque), so the run is the full sector",
       rot_frac > 1e-3, f"peak rotational KE fraction {rot_frac:.3f}"))

    # R5: both sound speeds from ONE shove, double(24): longitudinal along
    # the shove axis, transverse across it (the fit's far cells are r = 7..10)
    side = 24
    _c, M16, H16 = double_graph7(side)
    n16 = H16.shape[0] // 7
    t16, tr16, _f, _e, _ = verlet7(M16, H16, shove_u0(n16, KR.centre_cell(side)),
                                   40.0, sample=0.1)
    amp = np.linalg.norm(tr16, axis=2)
    vL, arrL = line_front(t16, amp, side, 1)
    vT, arrT = line_front(t16, amp, side, 0)
    out["r5"] = (vL, vT, arrL, arrT)
    A(("R5 one shove, two fronts on double(24): longitudinal along the "
       "shove axis vs 1/2, transverse across it vs 1/(2 sqrt 2)",
       (not np.isnan(vL)) and 0.85 * S_LONG <= vL <= 1.05 * S_LONG
       and (not np.isnan(vT)) and 0.85 * S_TRANS <= vT <= 1.03 * vgT,
       f"longitudinal {vL:.4f} ({len(arrL)} arrivals, window "
       f"[{0.85 * S_LONG:.3f}, {1.05 * S_LONG:.3f}]); transverse {vT:.4f} "
       f"({len(arrT)} arrivals, window [{0.85 * S_TRANS:.3f}, "
       f"{1.03 * vgT:.3f}] -- ceiling is the branch's max group velocity "
       f"plus the threshold's precursor bias)"))

    # R6: the stir -- circular polarization. On +-z the displacement
    # rotates at the drive frequency in the drive's sense (same lab sense
    # both sides: opposite helicity about their propagation directions),
    # and is transverse.
    sideS = 12
    _c, MS, HS = double_graph7(sideS)
    nS = HS.shape[0] // 7
    cS = KR.centre_cell(sideS)
    w0 = stir_omega()
    drv = stir_force(cS, w0)
    tS, trS, _f, _e, _ = verlet7(MS, HS, np.zeros(7 * nS), 30.0, sample=0.1,
                                 drive=drv)
    cc = sideS // 2
    r = 4
    up, dn = KR.site_index(sideS, cc, cc, cc + r), KR.site_index(sideS, cc, cc, cc - r)
    # window: after arrival (r / c_T) through the end of the drive's passage
    t0 = r / S_TRANS + 0.5
    t1 = min(t0 + drv.tdrive - 1.0, 30.0)
    su, fu = rotation_rate(tS, trS, up, t0, t1)
    sd, fd = rotation_rate(tS, trS, dn, t0, t1)
    out["r6"] = (w0, su, sd, fu, fd, t0, t1)
    A(("R6 the stir is CIRCULARLY polarized: on +z and -z the displacement "
       "rotates at the drive frequency, both in the drive's lab sense "
       "(opposite helicity about propagation), and is transverse",
       0.9 < su / w0 < 1.1 and 0.9 < sd / w0 < 1.1 and fu > 0.95 and fd > 0.95,
       f"omega {w0:.4f}; dphi/dt +z {su:.4f}, -z {sd:.4f}; transverse "
       f"fraction +z {fu:.3f}, -z {fd:.3f}; window [{t0:.1f}, {t1:.1f}]"))

    # R7: the detector detects -- reverse the drive, the sense reverses
    drv_r = stir_force(cS, -w0)
    tR, trR, _f, _e, _ = verlet7(MS, HS, np.zeros(7 * nS), 30.0, sample=0.1,
                                 drive=drv_r)
    sr, _fr = rotation_rate(tR, trR, up, t0, t1)
    out["r7"] = sr
    A(("R7 reversing the stir reverses the sense (the detector detects)",
       0.9 < -sr / w0 < 1.1, f"dphi/dt +z {sr:.4f} for omega {-w0:.4f}"))

    return checks, out


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "export":
        side = int(sys.argv[2]) if len(sys.argv) > 2 else 24
        tmax = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0
        sample = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
        path = sys.argv[5] if len(sys.argv) > 5 else "analysis/.pages/data/shear.json"
        p, data = export(side, tmax, path, sample=sample)
        print(f"exported double({side}) shove+stir, "
              f"{len(data['shove']['times'])} frames -> {p}; fronts "
              f"L {data['shove']['front_long']} T {data['shove']['front_trans']}")
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
