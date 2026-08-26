"""jb_fl -- letting the lattice breathe, and what that does and does not fix.

BEAD inviscid-k36. T2 23500 measured that the phase band's Gamma gap is exactly
sqrt(6)*|d lambda/d a| and concluded that the fixed lattice was the defect and
letting lambda respond was a prerequisite for everything. That is HALF right,
and this file establishes which half.

  IT FIXES k = 0 COMPLETELY. With the lattice vectors held fixed the coherent
  exchange stops being a zero mode of the bar framework the moment lambda has a
  slope: Gamma nullity is 4 at a = -30 and 3 at every other phase, and the mode
  that goes missing leaves behind a singular value of exactly
  2*sqrt(3)*|d lambda/d a|. Give the framework the nine lattice-velocity
  columns and the exchange is an exact zero mode AT EVERY PHASE, cost 3e-10 to
  3e-9 against a fixed-lattice cost rising to 1.78. The Goldstone claim is
  therefore true at every phase, not just the midpoint, and T2 23486 (b)(6)'s
  version of it was under-stated rather than wrong.

  IT DOES NOT FIX k != 0, AND THE GAP THERE IS REAL. The lattice-velocity
  degrees of freedom are a k = 0 object -- a uniform strain -- so they cannot
  reach finite wavevector. The obvious next hypothesis was that the projection
  onto the 2-d phase subspace was discarding the accompanying relaxation, which
  a constrained calculation cannot represent.
  It is not. Identifying bands in the FULL 18-band spectrum by overlap with the
  in-phase direction finds the phase-carrying band at the projected value to
  4e-5 (0.491151 at a = -20, 0.967379 at a = -10) and FLAT as k -> 0, while at
  a = -30 the same procedure gives a branch linear in k. Independently, the
  phase SPECTRAL WEIGHT does not migrate downward: below 0.05 it saturates at
  0.668 at a = -30 and at 0.138 / 0.080 at a = -20 / -10.

  A CARE NOTE, because the first draft of this file got it wrong and the gate
  caught it. Rayleigh-Ritz bounds the projected LOWEST eigenvalue against the
  full LOWEST one. It says NOTHING about the MAX-OVERLAP band, which is a
  different selection -- and indeed the projected value sits 3.6e-5 BELOW the
  max-overlap band at a = -20. The two agreeing away from the midpoint is a
  measurement, not a theorem, and F7 gates it as one.

SO THE CONCLUSION IS SHARPER AND WORSE THAN THE ONE IT REPLACES.

    A UNIFORM phase advance is free at every phase, once the lattice may breathe.
    A MODULATED one is not, and its cost does not vanish with wavelength.
    The phase wave is a propagating GAPLESS mode ONLY at the exchange midpoint.

AND THEN F10 SOFTENS THAT, WHICH IS WHY IT IS HERE AND NOT LEFT AS A CAVEAT.
A gap is not immobility, and the gapped sector turns out to disperse.

  AT the midpoint the phase spectral peak is CONTINUOUS across the whole [100]
  line -- zero jumps in 199 steps, biggest step 0.00371 against a median of
  0.00365 -- with median |d omega/d k| = 0.5344, which is F9's sqrt(8/27) read a
  completely different way. So it is a genuine branch, not merely a
  long-wavelength mode.

  AWAY from the midpoint it is PIECEWISE continuous: one handoff at a = -20 and
  two at a = -10, where the phase character transfers between bands. Between
  the handoffs it disperses with real group velocity, median 0.381 and 0.111.
  So the gapped phase sector DOES carry packets over the smooth stretches --
  a packet whose k-content straddles a handoff will split, but one that does
  not, propagates. And the velocity falls monotonically as the reference phase
  leaves the midpoint: 0.534, 0.381, 0.111 at a = -30, -20, -10.

  The response is never smeared: the dominant band always carries at least 0.319
  of the weight, against 1/18 = 0.056 for an even spread across all eighteen.

  METHOD, and why it is not the obvious one. Eigenvector continuity-following --
  matching each step's modes to the previous step's -- CANNOT do this here. Its
  worst matched overlap is 0.21 and REFINING THE SAMPLING EIGHTFOLD DOES NOT
  MOVE IT (0.2152, 0.2124, 0.2110, 0.2103 at n = 30, 60, 120, 240), so the
  ambiguity is a genuine degeneracy and not a step-size artifact. Tracking the
  PEAK FREQUENCY instead sidesteps it: band index hopping is a labelling
  artifact and harmless, while a jump in the peak is physical. That distinction
  is the whole content of the method.

THE MECHANISM, offered as an interpretation and labelled as one. Each cell's
lambda is slaved to its phase, so a modulated phase field demands a local
dilation proportional to d lambda/d a. Accommodating it needs a displacement
field of amplitude ~ delta-lambda / k, which diverges as k -> 0 -- but the
STRAIN that displacement produces is k*u ~ delta-lambda, independent of
wavelength. Elastic energy goes as strain squared, so the cost is wavelength
independent: a gap proportional to |d lambda/d a|, which is the measured law.
At a = -30 the slope is zero, phase decouples from dilation at linear order,
nothing has to strain, and the branch is gapless with c = 2/3.

METHOD NOTE THAT COST TWO ATTEMPTS, AND IS THE SAME TRAP IN A NEW PLACE. The
exchange direction has THREE contributions in the framework's own variables:
    1. the cell's own FOLD;
    2. the translation of its CENTRE, because origins are lam(a)*site;
    3. MINUS dL/da . n, because the variable is a node's REDUCED representative
       and the reduction is modulo a lattice that is moving.
All three vanish or coincide at a = -30. A construction missing either of the
last two therefore PASSES at the midpoint and fails everywhere else -- exactly
the shape of T2 23486 method note 6, where every verification to date had been
run on the one pose where the defect could not appear. Row F6 gates both terms
AND gates that they are invisible at the midpoint, so the trap is stated rather
than merely avoided.

A ROW DELIBERATELY NOT BUILT: the gapped branch's shape across the zone. Near
Gamma, max-overlap band identification is stable and converged, which is where
the gap claim lives. Away from Gamma it HOPS between bands as phase character
redistributes -- at a = -10 it reads 0.967, 1.293, 1.198, 1.174, 0.197, 0.220
along [100], and the last two are plainly a different branch. Characterising
bandwidth or group velocity needs continuity-following, which is not built here,
so no bandwidth is claimed. A gap is not immobility and the gapped sector may
well still carry packets; this file does not say either way.

Canonical prose state: T2 `inviscid/qvf-epic-consolidated-state`.
"""

from __future__ import annotations

import sys

import numpy as np

import jb_hc_honeycomb as HC
import jb_gp_plate_geometry as Z
from jb_x_array_linkage import PAIRS

A_REF = -30.0
PHASE_OFFSET = 60.0

#: The two sublattice sites of the periodic unit cell, in integer lattice units.
SITE = {0: np.zeros(3), 1: np.ones(3)}

TOL_ZERO = 1e-7         # flexible-lattice exchange cost, measured <= 3e-9
TOL_AGREE = 1e-8        # shared-node velocity agreement, measured <= 3e-11
TOL_RATIO = 1e-3        # lost sv / (2 sqrt3 |dlam/da|), measured spread 6e-4
TOL_PROJ = 1e-4         # projected-vs-full excess away from the midpoint,
                        # measured 3.6e-5 -- small but NOT zero: the
                        # projection is an upper bound everywhere, merely
                        # near-tight where the phase band is isolated.
COST_MIN = 0.1          # controls that must cost energy, measured >= 0.43
MID_LOW_WEIGHT = 0.5    # phase weight below 0.05 at the midpoint, measured 0.668
OFF_LOW_WEIGHT = 0.25   # ... away from it, measured 0.138 and 0.080

#: Phases, and a SECOND ABSOLUTE incommensurate arm written as absolute numbers.
PHASES = (-29.0, -25.0, -20.0, -15.0, -10.0, -5.0, -40.0, -45.0, -50.0, -55.0)
PHASES_ALT = (-7.3, -13.9, -21.1, -34.7, -43.3, -51.9)

FD = 1e-5


def dlambda(a):
    """d lambda / d a per DEGREE, Richardson-extrapolated."""
    d1 = (HC.lattice(a + 1e-3) - HC.lattice(a - 1e-3)) / 2e-3
    d2 = (HC.lattice(a + 5e-4) - HC.lattice(a - 5e-4)) / 1e-3
    return (4.0 * d2 - d1) / 3.0


def flex_rigidity(P, bars, A):
    """Periodic rigidity matrix with a FLEXIBLE lattice.

    Columns are 3 per node then 9 for the lattice velocity Ldot, row-major. For
    a bar from node i to node j in cell R,
        d/dt (|d|^2 / 2) = d . (v_i - v_j - Ldot R)
    so the lattice block is -(d_a R_b). Holding the lattice fixed simply deletes
    those nine columns, which is the model this file is correcting."""
    n = len(P)
    m = np.zeros((len(bars), 3 * n + 9))
    for r, (i, j, Rv) in enumerate(bars):
        R = np.array(Rv, dtype=float)
        d = P[i] - (P[j] + R * A)
        m[r, 3 * i:3 * i + 3] += d
        m[r, 3 * j:3 * j + 3] -= d
        m[r, 3 * n:] -= np.outer(d, R).ravel()
    return m


def exchange_direction(a, centre=True, reduce_corr=True, h=FD):
    """d/da of the coherent exchange, in the framework's OWN variables.

    `centre` and `reduce_corr` switch off contributions 2 and 3 of the method
    note, so F6 can measure what each is worth and show that both are invisible
    at the midpoint."""
    fw = HC.h4_framework(a)
    P, bars, A, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
    n = len(P)
    dl = (HC.lattice(a + h) - HC.lattice(a - h)) / (2 * h)
    dA = 2.0 * dl
    per = {}
    for (ci, ph, off, f, c, i, nn) in slots:
        key = (ci, i, tuple(nn))
        if key in per:
            continue
        fold = ((Z.corners(ph + h) + off)[f][c]
                - (Z.corners(ph - h) + off)[f][c]) / (2 * h)
        val = fold + (dl * SITE[ci] if centre else 0.0)
        if reduce_corr:
            val = val - dA * np.array(nn, dtype=float)
        per[key] = val
    by = {}
    for (ci, i, nn), val in per.items():
        by.setdefault(i, []).append(val)
    spread = max(float(np.abs(np.array(v) - np.array(v).mean(axis=0)).max())
                 for v in by.values())
    v = np.zeros(3 * n + 9)
    for i, vals in by.items():
        v[3 * i:3 * i + 3] = np.mean(vals, axis=0)
    v[3 * n:] = (dA * np.eye(3)).ravel()
    return v, P, bars, A, n, dl, spread


def cost(m, v):
    return float(np.linalg.norm(m @ v) / np.linalg.norm(v))


# ==========================================================================
# F1-F2: what the fixed lattice destroys, and by exactly how much
# ==========================================================================

def f1_nullity():
    rows = []
    for a in (A_REF,) + PHASES + PHASES_ALT:
        fw = HC.h4_framework(a)
        nul, s = HC.zero_modes(fw["P"], fw["bars"], fw["A"], np.zeros(3))
        rows.append((a, nul, float(np.sort(s)[3])))
    mid = [r for r in rows if r[0] == A_REF]
    off = [r for r in rows if r[0] != A_REF]
    ratios = [r[2] / (2.0 * np.sqrt(3.0) * abs(dlambda(r[0]) * 180.0 / np.pi))
              for r in off]
    return dict(rows=rows, mid_nullity=mid[0][1] if mid else -1,
                off_nullities=sorted({r[1] for r in off}),
                ratio_dev=max(abs(x - 1.0) for x in ratios), n=len(off))


# ==========================================================================
# F3-F6: the flexible lattice, its controls, and the construction trap
# ==========================================================================

def f3_flexible():
    flexi, fixed, spreads = [], [], []
    for a in (A_REF,) + PHASES + PHASES_ALT:
        v, P, bars, A, n, _dl, sp = exchange_direction(a)
        M = flex_rigidity(P, bars, A)
        Mfix = HC.bloch(P, bars, A, np.zeros(3), unit=False).real
        flexi.append(cost(M, v))
        fixed.append(cost(Mfix, v[:3 * n]))
        spreads.append(sp)
    return dict(flex_worst=max(flexi), fixed_max=max(fixed),
                spread=max(spreads), n=len(flexi))


def f4_controls():
    nb, ws, lo = [], [], []
    for a in (-25.0, -10.0, -45.0):
        v, P, bars, A, n, dl, _ = exchange_direction(a)
        M = flex_rigidity(P, bars, A)
        w = v.copy()
        w[3 * n:] = 0.0                                   # no breathing
        u = v.copy()
        u[3 * n:] = (-2.0 * dl * np.eye(3)).ravel()       # breathing backwards
        r = np.zeros_like(v)
        r[3 * n:] = v[3 * n:]                             # lattice only
        nb.append(cost(M, w))
        ws.append(cost(M, u))
        lo.append(cost(M, r))
    return dict(no_breath=min(nb), wrong_sign=min(ws), lattice_only=min(lo),
                n=len(nb))


def f6_construction():
    """Which of the two correction terms is load-bearing, and which is gauge.

    The first draft of this file asserted both were necessary. The gate said
    otherwise and it was right: dropping the CENTRE term costs exactly nothing.
    Every node is owned by one cell of each sublattice, so the per-node average
    shifts by half the centre term uniformly -- a rigid TRANSLATION, which is a
    zero mode of anything. It is a gauge choice, not a contribution. Only the
    reduced-representative term is load-bearing."""
    off_r, mid_r, trans = [], [], []
    for a in (-20.0, -10.0, -45.0):
        vr, P, bars, A, n, _d, _s = exchange_direction(a, reduce_corr=False)
        off_r.append(cost(flex_rigidity(P, bars, A), vr))
        full, _P, _b, _A, n2, _d2, _s2 = exchange_direction(a)
        noc, _P3, _b3, _A3, _n3, _d3, _s3 = exchange_direction(a, centre=False)
        diff = (full - noc)[:3 * n2].reshape(-1, 3)
        trans.append(float(np.abs(diff - diff.mean(axis=0)).max()))
    vr, P, bars, A, n, _d, _s = exchange_direction(A_REF, reduce_corr=False)
    mid_r.append(cost(flex_rigidity(P, bars, A), vr))
    return dict(off_reduce=min(off_r), mid_reduce=max(mid_r),
                translation_resid=max(trans), n=len(off_r))


# ==========================================================================
# F7-F8: and why none of it reaches finite wavevector
# ==========================================================================

def _phase_vector(P, slots, A, kv):
    B = HC.phase_basis(P, slots, A, kv)
    u = B[:, 0] + B[:, 1]
    return u / np.linalg.norm(u)


def f9_speed():
    """The midpoint speed, projected against free. The projection CONSTRAINS
    the motion to be pure phase, which no solution of the equations of motion
    is; it is a Rayleigh-Ritz upper bound and it is 22% high."""
    fw = HC.h4_framework(A_REF)
    P, bars, A, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
    G = np.pi / A
    eps = 0.005
    proj, full, elas = [], [], []
    for d in ((1, 0, 0), (1, 1, 0), (1, 1, 1)):
        u = np.array(d, dtype=float)
        u = u / np.linalg.norm(u)
        kv = eps * G * u
        q = HC.phase_basis(P, slots, A, kv)
        M = HC.bloch(P, bars, A, kv)
        D = M.conj().T @ M
        proj.append(float(np.sqrt(max(
            np.linalg.eigvalsh(q.conj().T @ D @ q)[0], 0.0))) / (eps * G))
        w, V = np.linalg.eigh(D)
        w = np.sqrt(np.maximum(w, 0.0))
        wt = np.abs(V.conj().T @ _phase_vector(P, slots, A, kv)) ** 2
        full.append(float(w[int(np.argmax(wt))]) / (eps * G))
        sv = np.sort(np.linalg.svd(M, compute_uv=False))
        elas.append(float([x for x in sv if x > 1e-9][0]) / (eps * G))
    return dict(proj=proj, full=full, elastic=elas,
                full_iso=max(full) - min(full),
                full_err=max(abs(x - np.sqrt(8 / 27)) for x in full),
                ratio=proj[0] / full[0],
                ratio_err=abs(proj[0] / full[0] - np.sqrt(1.5)),
                elastic_aniso=max(elas) - min(elas))


def f7_projection_is_not_the_defect():
    rows = []
    for a in (A_REF, -20.0, -10.0):
        fw = HC.h4_framework(a)
        P, bars, A, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
        G = np.pi / A
        kv = 0.005 * G * np.array([1.0, 0.0, 0.0])
        q = HC.phase_basis(P, slots, A, kv)
        M = HC.bloch(P, bars, A, kv)
        D = M.conj().T @ M
        proj = float(np.sqrt(max(np.linalg.eigvalsh(q.conj().T @ D @ q)[0], 0.0)))
        w, V = np.linalg.eigh(D)
        w = np.sqrt(np.maximum(w, 0.0))
        wt = np.abs(V.conj().T @ _phase_vector(P, slots, A, kv)) ** 2
        full = float(w[int(np.argmax(wt))])
        rows.append((a, proj, full, abs(proj - full)))
    off = [r for r in rows if r[0] != A_REF]
    mid = [r for r in rows if r[0] == A_REF]
    return dict(rows=rows, off_worst=max(r[3] for r in off),
                mid_excess=(mid[0][1] / mid[0][2]) if mid else float("nan"),
                signed_min=min(r[1] - r[2] for r in rows), n=len(rows))


def f8_spectral_weight():
    out = {}
    for a in (A_REF, -20.0, -10.0):
        fw = HC.h4_framework(a)
        P, bars, A, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
        G = np.pi / A
        series = []
        for eps in (0.05, 0.02, 0.01, 0.005, 0.002):
            kv = eps * G * np.array([1.0, 0.0, 0.0])
            M = HC.bloch(P, bars, A, kv)
            D = M.conj().T @ M
            w, V = np.linalg.eigh(D)
            w = np.sqrt(np.maximum(w, 0.0))
            wt = np.abs(V.conj().T @ _phase_vector(P, slots, A, kv)) ** 2
            wt = wt / wt.sum()
            series.append(float(wt[w < 0.05].sum()))
        out[a] = series
    return dict(series=out, mid=out[A_REF][-1],
                off=max(out[-20.0][-1], out[-10.0][-1]), n=len(out))


# ==========================================================================
# F10: the branch's SHAPE across the zone, which this file previously declared
#      a row deliberately not built. Bead inviscid-46y.
# ==========================================================================

#: k-path leg used for the shape rows, and the sampling that showed the
#: continuity-following method cannot be rescued by refinement.
PEAK_N = 200
FOLLOW_N = (30, 60, 120, 240)


def _peak_track(a0, n=PEAK_N):
    """omega of the MAX-WEIGHT band along [100], densely sampled.

    Band INDEX hopping is a labelling artifact and harmless; a jump in the PEAK
    FREQUENCY is not, because a packet needs a locally smooth omega(k). So this
    tracks the peak and not the label."""
    fw = HC.h4_framework(a0)
    P, bars, A, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
    G = np.pi / A
    ts = np.linspace(0.002, 1.0, n)
    peak, weight = [], []
    for t in ts:
        kv = t * G * np.array([1.0, 0.0, 0.0])
        M = HC.bloch(P, bars, A, kv)
        D = M.conj().T @ M
        w, V = np.linalg.eigh(D)
        w = np.sqrt(np.maximum(w, 0.0))
        wt = np.abs(V.conj().T @ _phase_vector(P, slots, A, kv)) ** 2
        wt = wt / wt.sum()
        i = int(np.argmax(wt))
        peak.append(w[i])
        weight.append(float(wt[i]))
    peak = np.array(peak)
    step = (ts[1] - ts[0]) * G
    d = np.abs(np.diff(peak))
    smooth = d[d < 0.05]
    return dict(jumps=int((d > 0.05).sum()), biggest=float(d.max()),
                vel_med=float(np.median(smooth)) / step,
                lo=float(peak.min()), hi=float(peak.max()),
                min_weight=min(weight))


def _follow_confidence(a0, n):
    """Worst matched eigenvector overlap along G-X-M-G at sampling `n`."""
    from scipy.optimize import linear_sum_assignment
    fw = HC.h4_framework(a0)
    P, bars, A = fw["P"], fw["bars"], fw["A"]
    G = np.pi / A
    path = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 0, 0)]
    ks = []
    for i in range(len(path) - 1):
        s0 = np.array(path[i], dtype=float)
        e0 = np.array(path[i + 1], dtype=float)
        for t in np.linspace(0, 1, n, endpoint=(i == len(path) - 2)):
            ks.append(G * (s0 + t * (e0 - s0)))
    Vp, worst = None, 1.0
    for kv in ks:
        M = HC.bloch(P, bars, A, kv)
        w, V = np.linalg.eigh(M.conj().T @ M)
        if Vp is not None:
            ov = np.abs(Vp.conj().T @ V) ** 2
            r, c = linear_sum_assignment(-ov)
            worst = min(worst, float(ov[r, c].min()))
            V = V[:, c]
        Vp = V
    return worst


def f10_shape():
    tracks = {a: _peak_track(a) for a in (A_REF, -20.0, -10.0)}
    conf = {n: _follow_confidence(A_REF, n) for n in FOLLOW_N}
    return dict(tracks=tracks, conf=conf,
                conf_spread=max(conf.values()) - min(conf.values()),
                conf_worst=max(conf.values()))


# ==========================================================================
# C: THE COMPLIANCE DECISION, which the record has carried as open since the
#    first dispersion was computed. It is a FORK, not a missing scale factor.
# ==========================================================================

def c_compliance():
    """Where the compliance lives, and what each choice does to the phase field.

    A perfectly rigid jitterbug carries no waves: its mechanism motions cost
    nothing and everything else is infinitely stiff. So a wave requires putting
    compliance SOMEWHERE, and where it goes is a modelling decision the geometry
    cannot make.

        STRUTS compliant, hinges free -- what this epic implements. The bars are
        springs, the jitterbug motion stretches none of them at k = 0, so the
        phase field is a GOLDSTONE mode and gapless. This is the right choice
        for a strut-and-pin structure: real struts are elastic, ideal pins are
        frictionless.

        HINGES compliant, struts rigid -- the alternative, right for a folded
        sheet rather than a strut frame. The plates' outward normals are phase
        INVARIANT (Z0), so the dihedral between any two fixed plates never
        moves and the single hinge degree of freedom per cell IS the fold angle.
        A torsional spring on it is therefore exactly (kappa/2)(a - a0)^2 PER
        CELL -- on-site, a MASS term, and it gaps the phase field at ANY
        stiffness however small.

    The two differ in KIND. Under the first the exchange is free and the phase
    wave is gapless; under the second the exchange costs energy at every hinge
    and the phase field is massive everywhere, midpoint included. Every gapless
    result in this file is downstream of the first choice, and that choice has
    never been declared until now."""
    dihedral_spread = 0.0
    for a in (0.0, -15.0, -30.0, -45.0, -60.0):
        vals = []
        for _v, slots in enumerate(PAIRS):
            (fa, _ca), (fb, _cb) = slots
            na, nb = Z.plate_normal(fa), Z.plate_normal(fb)
            vals.append(abs(float(na @ nb)))
        dihedral_spread = max(dihedral_spread, max(vals) - min(vals))
    fw = HC.h4_framework(A_REF)
    P, bars, A, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
    G = np.pi / A
    q, _ = np.linalg.qr(HC.phase_basis(P, slots, A, 1e-6 * G * np.array([1.0, 0, 0])))
    M = HC.bloch(P, bars, A, 1e-6 * G * np.array([1.0, 0, 0]))
    w0 = float(np.sqrt(max(np.linalg.eigvalsh(
        q.conj().T @ (M.conj().T @ M) @ q)[0], 0.0)))
    gapped = [float(np.sqrt(k + w0 ** 2)) for k in (1e-4, 1e-2, 1.0)]
    return dict(dihedral_spread=dihedral_spread, strut_gap=w0,
                hinge_gaps=gapped, n=3)


# ==========================================================================
# A: ANHARMONIC. Do the two soft families stay orthogonal past linear order?
#    T2 23486 open question 1. They do not.
# ==========================================================================

#: Directions probed for the second-harmonic channel, and their wavevector
#: fractions. [011] is a second <110> member (the channel must not be an
#: artifact of one representative) and [210] is a deliberate NON-symmetry
#: direction, so "only <110>" is measured rather than assumed from the three
#: high-symmetry labels.
ANH_DIRS = (("[110]", (1, 1, 0)), ("[100]", (1, 0, 0)), ("[111]", (1, 1, 1)),
            ("[011]", (0, 1, 1)), ("[210]", (2, 1, 0)))
ANH_T = (0.15, 0.30, 0.45)
ANH_SOFT_TOL = 1e-6


def _cubic(P, bars, A, k1, u1, k2, u2, k3, u3):
    """Three-mode coupling from the bar energy past quadratic order.

    Expanding (|d| - L)^2 / 2 with d = d0 + delta gives an extension
    nhat.delta + (|delta|^2 - (nhat.delta)^2) / 2L, and the cross term of those
    two is the cubic vertex. Symmetrised over the three legs."""
    tot = 0.0 + 0j
    for (i, j, Rv) in bars:
        Rw = np.array(Rv, dtype=float) * A
        dv = P[i] - (P[j] + Rw)
        L = np.linalg.norm(dv)
        nh = dv / L

        def leg(k, u):
            return u[3 * i:3 * i + 3] - u[3 * j:3 * j + 3] * np.exp(1j * float(k @ Rw))
        d1, d2, d3 = leg(k1, u1), leg(k2, u2), leg(k3, u3)
        tot += ((nh @ d1) * ((d2 @ d3) - (nh @ d2) * (nh @ d3))
                + (nh @ d2) * ((d1 @ d3) - (nh @ d1) * (nh @ d3))
                + (nh @ d3) * ((d1 @ d2) - (nh @ d1) * (nh @ d2))) / (6.0 * L)
    return abs(complex(tot))


def a_anharmonic():
    """A phase wave at q forces the lattice at 2q. Where does that force land?

    The two soft families are orthogonal at LINEAR order -- jb_hc H8 measures it
    at M, where the positional band is 0 while both phase bands sit at
    2/sqrt(3). This asks whether that survives finite amplitude, and it does
    not: along <110> the second harmonic of a <110> wavevector is still <110>,
    which is exactly floppy, so the forcing lands on a mode with NO RESTORING
    FORCE and accumulates rather than oscillating.

    The distinction the control row protects: the anharmonic coupling itself is
    ISOTROPIC. What is anisotropic is the SINK."""
    fw = HC.h4_framework(A_REF)
    P, bars, A, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
    G = np.pi / A
    rows = []
    for lbl, d in ANH_DIRS:
        u = np.array(d, dtype=float)
        u = u / np.linalg.norm(u)
        for t in ANH_T:
            q = t * G * u
            k3 = -2.0 * q
            B = HC.phase_basis(P, slots, A, q)
            up = B[:, 0] + B[:, 1]
            up = up / np.linalg.norm(up)
            M3 = HC.bloch(P, bars, A, k3)
            w3, V3 = np.linalg.eigh(M3.conj().T @ M3)
            w3 = np.sqrt(np.maximum(w3, 0.0))
            cs = [_cubic(P, bars, A, q, up, q, up, k3, V3[:, m])
                  for m in range(V3.shape[1])]
            soft = [c for c, w in zip(cs, w3) if w < ANH_SOFT_TOL]
            rows.append(dict(dir=lbl, t=t, total=max(cs), n_soft=len(soft),
                             into_soft=max(soft) if soft else 0.0))
    is110 = {"[110]", "[011]"}
    return dict(rows=rows,
                soft_110=min(r["into_soft"] for r in rows if r["dir"] in is110),
                soft_other=max(r["into_soft"] for r in rows if r["dir"] not in is110),
                n_soft_110=min(r["n_soft"] for r in rows if r["dir"] in is110),
                n_soft_other=max(r["n_soft"] for r in rows if r["dir"] not in is110),
                total_lo=min(r["total"] for r in rows),
                total_hi=max(r["total"] for r in rows), n=len(rows))


# ==========================================================================
# THE GATE
# ==========================================================================

def gate(f1, f3, f4, f6, f7, f8, f9, f10, c, an):
    checks = []
    R = checks.append

    R(("F1  the FIXED lattice destroys the exchange as a zero mode: nullity 4 "
       "at a = -30, 3 at every other phase -- CAN FAIL both ways",
       f1["mid_nullity"] == 4 and f1["off_nullities"] == [3] and f1["n"] > 0,
       f"mid {f1['mid_nullity']}, off {f1['off_nullities']} "
       f"over {f1['n']} phases", "4 and [3]"))
    R(("F2  and the mode it destroys leaves EXACTLY 2*sqrt(3)*|d lambda/d a| "
       "behind", f1["ratio_dev"] < TOL_RATIO,
       f"worst |ratio-1| {f1['ratio_dev']:.2e}", f"< {TOL_RATIO:.0e}"))

    R(("F3  with the lattice FREE the exchange is an exact zero mode at EVERY "
       "phase", f3["n"] > 0 and f3["flex_worst"] < TOL_ZERO,
       f"worst {f3['flex_worst']:.2e} over {f3['n']} phases",
       f"< {TOL_ZERO:.0e}"))
    R(("F3  CONTROL: the same direction costs real energy with the lattice "
       "pinned, so F3 is not vacuous -- CAN FAIL",
       f3["fixed_max"] > 1.0, f"max fixed cost {f3['fixed_max']:.6f}", "> 1.0"))
    R(("F3  the two cells sharing a node AGREE on its velocity -- the contact "
       "is maintained, which is what validates the direction",
       f3["spread"] < TOL_AGREE, f"worst {f3['spread']:.2e}",
       f"< {TOL_AGREE:.0e}"))

    R(("F4  CONTROLS: no breathing, backwards breathing and lattice-only all "
       "cost energy -- CAN FAIL",
       f4["n"] > 0 and min(f4["no_breath"], f4["wrong_sign"],
                           f4["lattice_only"]) > COST_MIN,
       f"{f4['no_breath']:.3f} / {f4['wrong_sign']:.3f} / "
       f"{f4['lattice_only']:.3f}", f"all > {COST_MIN}"))

    R(("F6  the reduced-representative term IS load-bearing away from the "
       "midpoint", f6["n"] > 0 and f6["off_reduce"] > COST_MIN,
       f"cost without it {f6['off_reduce']:.3f}", f"> {COST_MIN}"))
    R(("F6  THE TRAP, gated so it is stated and not merely avoided: it is "
       "INVISIBLE at a = -30, so a construction missing it passes there and "
       "only there", f6["mid_reduce"] < TOL_ZERO,
       f"cost at the midpoint {f6['mid_reduce']:.2e}", f"< {TOL_ZERO:.0e}"))
    R(("F6  the CENTRE term is a pure GAUGE -- dropping it shifts the whole "
       "direction by a rigid translation, which this row measures rather "
       "than the first draft's claim that it was necessary",
       f6["translation_resid"] < TOL_AGREE,
       f"non-uniform residual {f6['translation_resid']:.2e}",
       f"< {TOL_AGREE:.0e}"))

    R(("F7  AWAY from the midpoint the projected band and the FULL spectrum's "
       "max-overlap band agree, so the gap there is not a projection artifact",
       f7["n"] > 0 and f7["off_worst"] < TOL_PROJ,
       f"worst difference {f7['off_worst']:.2e}", f"< {TOL_PROJ:.0e}"))
    R(("F7  AT the midpoint it is NOT exact -- it is a variational UPPER "
       "bound, high by exactly sqrt(3/2) -- CAN FAIL",
       abs(f7["mid_excess"] - np.sqrt(1.5)) < 1e-4,
       f"projected/full {f7['mid_excess']:.7f}",
       f"sqrt(3/2) = {np.sqrt(1.5):.7f}"))

    R(("F9  THE MIDPOINT SPEED IS sqrt(8/27), NOT 2/3, and it is ISOTROPIC",
       f9["full_err"] < 1e-4 and f9["full_iso"] < 1e-4,
       f"{f9['full'][0]:.7f} isotropic to {f9['full_iso']:.1e}",
       f"sqrt(8/27) = {np.sqrt(8 / 27):.7f}"))
    R(("F9  2/3 is the CONSTRAINED value and exceeds it by exactly sqrt(3/2)",
       f9["ratio_err"] < 1e-4, f"ratio {f9['ratio']:.7f}",
       f"{np.sqrt(1.5):.7f}"))
    R(("F9  CONTROL: the lowest ELASTIC branch is anisotropic over the same "
       "three directions, so F9's isotropy is not a property of every band "
       "-- CAN FAIL", f9["elastic_aniso"] > 0.1,
       f"elastic spread {f9['elastic_aniso']:.6f}", "> 0.1"))
    R(("F8  phase spectral weight stays UP away from the midpoint and goes "
       "DOWN at it -- the gap is real, not a labelling choice",
       f8["mid"] > MID_LOW_WEIGHT and f8["off"] < OFF_LOW_WEIGHT,
       f"midpoint {f8['mid']:.3f}, away {f8['off']:.3f}",
       f"> {MID_LOW_WEIGHT} and < {OFF_LOW_WEIGHT}"))

    mid = f10["tracks"][A_REF]
    off = [f10["tracks"][a] for a in (-20.0, -10.0)]
    R(("F10 AT the midpoint the phase peak is CONTINUOUS across the whole "
       "[100] line -- a genuine branch, not just a long-wavelength mode",
       mid["jumps"] == 0,
       f"{mid['jumps']} jumps in {PEAK_N - 1} steps, biggest "
       f"{mid['biggest']:.5f}", "0 jumps"))
    R(("F10 and its group velocity is sqrt(8/27), matching F9's speed read a "
       "different way", abs(mid["vel_med"] - np.sqrt(8 / 27)) < 2e-2,
       f"median |dw/dk| {mid['vel_med']:.4f}",
       f"{np.sqrt(8 / 27):.4f}"))
    R(("F10 AWAY from the midpoint it is only PIECEWISE continuous, so the "
       "phase sector is NOT a single branch there -- CAN FAIL",
       all(1 <= t["jumps"] <= 3 for t in off),
       f"jumps {[t['jumps'] for t in off]}, biggest "
       f"{max(t['biggest'] for t in off):.4f}", "1..3 each"))
    R(("F10 BUT IT STILL DISPERSES between the handoffs: a gap is not "
       "immobility, and this SOFTENS the k != 0 conclusion above",
       all(t["vel_med"] > 0.05 for t in off),
       f"median |dw/dk| {[round(t['vel_med'], 4) for t in off]}", "> 0.05"))
    R(("F10 the phase response is never SMEARED: one band always carries a "
       "large share, against 1/18 = 0.056 for an even spread",
       min(t["min_weight"] for t in f10["tracks"].values()) > 0.2,
       f"min max-weight {min(t['min_weight'] for t in f10['tracks'].values()):.3f}",
       "> 0.2 (even = 0.056)"))
    R(("F10 CONTROL: eigenvector continuity-following CANNOT do this, and "
       "refinement does not rescue it -- which is why the peak is tracked "
       "instead of the label -- CAN FAIL",
       f10["conf_worst"] < 0.4 and f10["conf_spread"] < 0.01,
       f"confidence {[round(v, 4) for v in f10['conf'].values()]} over "
       f"8x refinement", "< 0.4, flat"))

    R(("C   the ONE hinge degree of freedom per cell IS the fold angle: plate "
       "dihedrals are phase-invariant, so a hinge spring is an ON-SITE term",
       c["dihedral_spread"] < 1e-12,
       f"dihedral spread {c['dihedral_spread']:.1e} over 5 phases", "< 1e-12"))
    R(("C   so THE GAPLESSNESS IS A CONSEQUENCE OF PUTTING COMPLIANCE IN THE "
       "STRUTS, not of the geometry: struts-compliant is gapless",
       c["strut_gap"] < 1e-5, f"omega(k->0) {c['strut_gap']:.2e}", "< 1e-5"))
    R(("C   CONTROL: ANY hinge stiffness gaps it, however small -- the two "
       "choices differ in KIND and not by a scale factor -- CAN FAIL",
       c["n"] > 0 and all(g > 1e-3 for g in c["hinge_gaps"]),
       f"kappa/I = 1e-4, 1e-2, 1 -> {[round(g, 4) for g in c['hinge_gaps']]}",
       "all gapped"))
    R(("C   PRINTED NOT GATED, the way jb_x discloses box_forced: omega = "
       "sqrt(K/M)*sigma, so every speed in this file is a PURE NUMBER and none "
       "is a physical speed until K and M are named",
       True, f"c[100] = {np.sqrt(8 / 27):.7f} x sqrt(K/M)", "printed"))

    R(("A   the SECOND HARMONIC of a <110> phase wave is still <110>, hence "
       "exactly FLOPPY -- a channel with no restoring force",
       an["n_soft_110"] >= 1 and an["n_soft_other"] == 0,
       f"soft modes at 2q: {an['n_soft_110']} along <110>, "
       f"{an['n_soft_other']} elsewhere", ">=1 vs 0"))
    R(("A   and the cubic coupling INTO that channel is NONZERO, so a finite-"
       "amplitude phase wave forces a mode that cannot push back -- CAN FAIL",
       an["soft_110"] > 1e-3, f"min |C| into it {an['soft_110']:.3e}", "> 1e-3"))
    R(("A   CONTROL: elsewhere it is EXACTLY zero, and NOT because the coupling "
       "vanishes -- the TOTAL anharmonic strength is comparable in every "
       "direction, so the anisotropy is in the SINK and not the vertex",
       an["soft_other"] == 0.0 and an["total_lo"] > 0.1
       and an["total_hi"] / an["total_lo"] < 2.0,
       f"into-soft {an['soft_other']:.1e}, total "
       f"{an['total_lo']:.3f}..{an['total_hi']:.3f}", "0, and total flat"))
    R(("A   SO THE TWO SOFT FAMILIES ARE NOT ANHARMONICALLY ORTHOGONAL -- "
       "T2 23486 open question 1, answered NO",
       an["n"] > 0 and an["soft_110"] > 1e-3 and an["soft_other"] == 0.0,
       f"coupled along the six <110>, decoupled elsewhere", "measured"))

    print()
    print("=" * 78)
    print(f"GATE  {len(checks)} rows")
    print("=" * 78)
    for name, ok, val, crit in checks:
        print(f"  {'PASS' if ok else 'FAIL':4s}  {name:66s} {str(val):>26s} {str(crit):>18s}")

    print()
    print("  ROWS THAT EXIST ONLY TO STOP ANOTHER ROW BEING UNFALSIFIABLE:")
    print("   * F3's pinned-lattice control. 'The exchange costs nothing'")
    print("     is satisfied by a rigidity matrix that is simply small, or by")
    print("     a direction vector that is mostly zero.")
    print("   * F3's shared-node agreement row. The exchange direction is")
    print("     assembled per cell, and two cells own every node; if they")
    print("     disagreed about its velocity the direction would not describe")
    print("     a motion of the structure at all, however cheap it looked.")
    print("   * F4's three controls, and F6's two. Without them 'the lattice")
    print("     response is correct' is satisfied by ANY lattice response,")
    print("     including none and including backwards.")
    print("   * F6's SECOND row above all. Both construction terms vanish at")
    print("     a = -30, so a wrong construction passes there. That is the")
    print("     same shape as T2 23486 method note 6, where every check to")
    print("     date had been run on the one pose the defect could not reach,")
    print("     and it cost two attempts here before being noticed.")
    print("   * F8. F7 alone shows the projected and full calculations agree,")
    print("     which they could do while BOTH mislabel which band is the")
    print("     phase one. The spectral weight is label-free: it sums to 1")
    print("     over the whole spectrum by construction.")
    print()
    print("  ROWS DELETED RATHER THAN FIXED: a free-boundary CHAIN test, which")
    print("  was the first plan for deciding whether a modulated phase field")
    print("  is soft. A free chain of N jitterbugs has nullity 6N+6 -- roughly")
    print("  six floppy modes per cell -- so every mode is soft and there is")
    print("  no signal to read. The bulk medium is over-constrained (Maxwell")
    print("  -6); a chain is not the medium, and the test would have measured")
    print("  its boundary.")
    print()
    print("  A ROW DELIBERATELY NOT BUILT: the gapped branch's shape across the")
    print("  zone. Max-overlap band identification is stable near Gamma, which")
    print("  is where the gap claim lives, and HOPS away from it as phase")
    print("  character redistributes -- at a = -10 it reads 0.967, 1.293,")
    print("  1.198, 1.174, 0.197, 0.220 along [100] and the last two are")
    print("  plainly another branch. No bandwidth or group velocity is claimed.")
    print("  A gap is not immobility; whether the gapped sector still carries")
    print("  packets is open, and needs continuity-following to answer.")
    print()
    print("  WHAT REMAINS UNMODELLED: the COMPLIANCE decision (rigid bars give")
    print("  infinite wave speed; K and M enter only as sqrt(K/M), so every")
    print("  speed here is a shape ratio and none is physical until that is")
    print("  made). And the medium in MOTION is described by none of this: the")
    print("  exchange path is flat, so a real medium does not sit at any a0,")
    print("  and harmonic analysis about a fixed point is the wrong tool for")
    print("  it. That is anharmonic territory and is not entered here.")

    failed = [n for n, ok, _v, _c in checks if not ok]
    print()
    if failed:
        print(f"  !! {len(failed)} CHECK(S) FAILED -- this is a bug report, not a")
        print("     measurement. Nothing above may enter the record.")
        for n in failed:
            print(f"       - {n}")
        return 1
    print("  ALL CHECKS PASSED.")
    return 0


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_fl -- letting the lattice breathe, and the half it does not fix")
    print("=" * 78)
    print("  A UNIFORM phase advance is free at every phase once the lattice")
    print("  may breathe: the fixed lattice destroyed it, and restoring nine")
    print("  lattice-velocity columns restores it exactly. A MODULATED one is")
    print("  not, and its cost does not vanish with wavelength. So the phase")
    print("  wave is a propagating gapless mode only at the exchange midpoint,")
    print("  where lambda is stationary and phase decouples from dilation.")
    f1 = f1_nullity()
    f3 = f3_flexible()
    f4 = f4_controls()
    f6 = f6_construction()
    f7 = f7_projection_is_not_the_defect()
    f8 = f8_spectral_weight()
    f9 = f9_speed()
    f10 = f10_shape()
    c = c_compliance()
    an = a_anharmonic()

    print()
    print("-" * 78)
    print("  THE GAMMA POINT, FIXED LATTICE vs FREE")
    print("-" * 78)
    print(f"    {'a0':>7s} {'d lam/d a':>12s} {'nullity':>8s} {'4th sv':>12s}")
    for (a, nul, sv) in f1["rows"][:9]:
        print(f"    {a:7.1f} {dlambda(a) * 180 / np.pi:+12.7f} {nul:8d} "
              f"{sv:12.7f}")
    print(f"  flexible-lattice exchange cost, worst over "
          f"{f3['n']} phases: {f3['flex_worst']:.3e}")
    print()
    print("  PHASE SPECTRAL WEIGHT BELOW 0.05, as k -> 0")
    print(f"    {'a0':>7s}   eps = 0.05, 0.02, 0.01, 0.005, 0.002")
    for a, ser in f8["series"].items():
        print(f"    {a:7.1f}   " + "  ".join(f"{x:.4f}" for x in ser))
    return gate(f1, f3, f4, f6, f7, f8, f9, f10, c, an)


if __name__ == "__main__":
    with np.errstate(all="ignore"):
        sys.exit(main())
