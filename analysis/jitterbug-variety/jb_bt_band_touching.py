"""jb_bt -- the band-touching planes, re-scanned in the FAITHFUL model.

THE ROW jb_hc DELIBERATELY DID NOT BUILD. `jb_hc_honeycomb.py` closes with:
"A ROW DELIBERATELY NOT BUILT: the band-touching PLANE family in the faithful
(projected) model. H8 checks M only." This file builds it, and the answer is
not the one the record expected.

WHAT WAS INHERITED. T2 `inviscid/qvf-epic-consolidated-state` (b)(5), from the
REDUCED phase-field model in which cell centres are pinned and shared vertices
are allowed to SPLIT under a penalty alpha:

    omega_-^2 = alpha*(8 - |f|),  omega_+^2 = alpha*(8 + |f|),
    f(k) = 8 * prod_i cos(k_i * lambda)

so the two bands are degenerate wherever any k_i = +-pi/(2*lambda), and the
record reads that as "Group velocity goes to zero there. A wave packet built
near those planes does not propagate -- it STALLS", flagged as consistent-but-
not-established because only M and R had been checked faithfully.

WHAT THIS FILE FINDS

  B1  THE LOCUS SURVIVES, EXACTLY. In the faithful model the two projected
      phase bands are degenerate on the six planes k_i = +-pi/A and NOWHERE
      else. Since A = 2*lambda those planes ARE the simple-cubic Brillouin
      zone boundary, which is the more useful way to say it.

  B2  THE DEGENERATE VALUE DISPERSES, and the reduced model said it would not.
      On the plane k_1 = pi/A, with L = A/2,

          omega^2 = 4*(1 - sin^2(k_2 L) sin^2(k_3 L))
                    / (3*(cos^2(k_2 L) + cos^2(k_3 L)))

      verified against the projection to 5e-11. It runs from sqrt(2/3) at the
      face centre to 2/sqrt(3) at the face edges -- a ratio of exactly sqrt(2)
      -- where the reduced model has it constant at sqrt(8*alpha). The extra
      dispersion is the VE-VE square-face contact, which the reduced model has
      no term for: it coupled cells only through the 8 triangular faces.

  B3  THE STALL IS REFUTED. In-plane group velocity is 0.2-0.37 (zero only at
      the face centre, where the point symmetry forces it); across the plane
      the two branches meet with EQUAL AND OPPOSITE slopes ~0.41-0.48, so this
      is a linear CROSSING, not a tangential touching and not a stall. There
      is no localisation mechanism here. The recorded claim was over-stated
      even in the reduced model it came from, where the in-plane velocity is
      zero but the normal velocity is not.

  B4  THE MECHANISM, in two conditions that are separately measurable. The
      projection is a 2x2 pencil (B^H D B, B^H B). Degeneracy needs BOTH
      diagonals equal -- the two sublattices must be equivalent -- AND the
      off-diagonal to vanish. The off-diagonal is proportional to the
      inter-sublattice structure factor 8*prod cos(k_i L) to a relative spread
      of 1e-15, because the ONLY inter-sublattice coupling is the eight
      body-diagonal triangular-face neighbours. That factor vanishes on the
      zone boundary. Nothing else does.

  B5  AND IT IS A PROPERTY OF ONE PHASE, NOT OF THE MEDIUM. The sublattices
      are equivalent only at a0 = -30, where the VE and the hole cell are
      congruent: 24 of the 48 octahedral operations carry one onto the other,
      and combined with tau = (A/2)(1,1,1) they are symmetries of the whole
      honeycomb. Off a0 = -30 the degeneracy opens LINEARLY in |a0 + 30| and
      the whole plane family is gone. It is not bcc band folding -- no pure
      translation maps the sublattices onto each other (Hausdorff 0.8165), so
      the translation group really is simple cubic with two cells per cell.

  B6  THE LARGER CASUALTY: THE PHASE BAND IS GAPLESS ONLY AT a0 = -30, AND THE
      GAP HAS A CLOSED FORM.

          gap(a0) = sqrt(6) * |d lambda / d a|        (a in radians)

      to 10 significant figures at every phase tested, isotropic, and zero
      exactly where d lambda/d a = 0 -- which is a0 = -30, because lambda is
      the fold half-diagonal and -30 is its MAXIMUM (exactly 2/sqrt(3)). Both
      models hold cell centres on a fixed lattice. A uniform phase shift is
      free only when the lattice would not have had to breathe to accommodate
      it. So c = 2/3, the gapless acoustic branch and the Goldstone
      identification are all statements about the exchange midpoint under a
      fixed-lattice assumption, not about the medium in general.

      This is not a new caveat invented here: T2 open question 2 ("letting
      lambda respond dynamically") already flagged the fixed lattice as the
      obvious next refinement. What is new is that it is not optional. Every
      wave number this epic has is measured at the one phase where the defect
      it introduces happens to vanish.

METHOD NOTE. The two conditions in B4 are why the locus survived a model
change that moved every number by 7-9x: one is a structure factor, fixed by
which cells touch, and the other is a symmetry. Neither depends on the elastic
constants. That is also why the localisation did NOT survive -- it was never
carried by either condition, only by the reduced model's missing term.

Canonical prose state: T2 `inviscid/qvf-epic-consolidated-state`.
"""

from __future__ import annotations

import itertools as it
import sys

import numpy as np

import jb_hc_honeycomb as HC
import jb_z_quasistatic_array as Z

# ---------------------------------------------------------------------------
# CONSTANTS. Thresholds are re-declared locally (house rule) and every one is
# priced from a number this run prints, named in the comment beside it.
# ---------------------------------------------------------------------------

#: Reference phase. Re-declared rather than imported, per the mutation-probe rule.
A_REF = -30.0
PHASE_OFFSET = 60.0

#: Reduced-model face-energy Hessian, T2 23484. Used ONLY to evaluate that
#: model's own closed form for the B3 contrast; nothing here depends on it.
ALPHA_REDUCED = 1.827705e-03

TOL_DEGEN = 1e-8        # priced from on-plane splitting, measured <= 2.5e-10
TOL_FORM = 1e-8         # priced from closed-form residual, measured 5.2e-11
TOL_PROP = 1e-10        # priced from |S01|/|8 prod cos| spread, measured 1.0e-15
TOL_SLOPE = 1e-5        # priced from |slope+| - |slope-|, measured <= 5.9e-7
TOL_GAP = 1e-7          # priced from gap/(sqrt(6)|dlam/da|) - 1, measured <= 1.1e-9
TOL_EXACT = 1e-12

SPLIT_OFF_MIN = 1.0e-2  # off-plane splitting at matched |k|, measured >= 0.127
SPAN_BAND = (0.20, 0.50)      # in-plane band span, measured 0.338204 -- TWO SIDED
VZERO = 1.0e-3          # "zero" group velocity; measured in-plane at the face
                        # centre is 0.0, elsewhere >= 0.198
VIN_MIN = 5.0e-2        # in-plane |v| off the face centre, measured >= 0.198
VN_MIN = 1.0e-1         # normal |v|, measured >= 0.409
EXP_BAND = (0.90, 1.10)       # splitting exponent in |a0+30|, measured 1.000

#: In-plane sample grid, and a SECOND ABSOLUTE incommensurate arm. The second
#: grid is written out as absolute numbers, NOT derived from the first by any
#: ratio -- ratio-derived second arms are a recorded bug in this project.
PLANE_GRID = tuple(round(x, 6) for x in np.linspace(-0.94, 0.94, 9))
PLANE_GRID_ALT = (-0.8713, -0.6047, -0.3391, -0.1127, 0.0713,
                  0.2887, 0.5231, 0.7079, 0.9337)

#: Phase sweep for the a0-specificity rows, and its own absolute second arm.
PHASE_SWEEP = (-29.0, -27.0, -25.0, -20.0, -15.0, -10.0, -5.0)
PHASE_SWEEP_ALT = (-33.7, -36.1, -41.3, -44.9, -51.7, -55.3, -57.1)

#: Offsets from -30 for the linearity fit, absolute, spanning three decades.
DELTA_SWEEP = (0.01, 0.03, 0.1, 0.3, 1.0, 3.0)

#: Finite-difference steps for B7's invariance control. Six decades: a
#: roundoff-carried quantity moves across these, a real one does not. This is
#: the control that would have caught the recorded optical gap 1.486784.
H_SWEEP = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)

#: B8's detuning sweep, straddling the finite-difference roundoff floor
#: (basis cond ~ 3.3e-10). STRONG values dominate it; WEAK values do not.
DETUNE_STRONG = (1e-4, 1e-2, 0.5, 1.0)
DETUNE_WEAK = (0.0, 1e-12, 1e-10)

#: Interior scan size and seed. Fixed so two runs are byte-identical.
INTERIOR_N = 1500
INTERIOR_SEED = 20260825


# ---------------------------------------------------------------------------
# THE MODEL. jb_hc's `phase_basis` recomputes every cell's phase-derivative
# vector on every call; those vectors do not depend on k, only the Bloch
# factors do. Precomputing them is a pure speedup and is verified against
# jb_hc's own function in row B0.
# ---------------------------------------------------------------------------

class Model:
    """The faithful projection at one reference phase, with the k-independent
    part of the phase basis precomputed."""

    def __init__(self, a, h=1e-5, reverse=False):
        fw = HC.h4_framework(a)
        self.a = a
        self.P, self.bars, self.A, self.slots = (fw["P"], fw["bars"],
                                                 fw["A"], fw["slots"])
        self.n = len(self.P)
        self.G = np.pi / self.A
        self.L = self.A / 2.0
        stems, seen = [], set()
        for (ci, ph, off, f, c, i, nn) in (list(reversed(self.slots))
                                           if reverse else self.slots):
            if (ci, i, nn) in seen:
                continue
            seen.add((ci, i, nn))
            dp = ((Z.corners(ph + h) + off)[f][c]
                  - (Z.corners(ph - h) + off)[f][c]) / (2 * h)
            stems.append((ci, i, np.array(nn, dtype=float) * self.A, dp))
        self.stems = stems

    def basis(self, kvec):
        B = np.zeros((3 * self.n, 2), dtype=complex)
        for (ci, i, R, dp) in self.stems:
            B[3 * i:3 * i + 3, ci] += dp * np.exp(-1j * float(np.dot(kvec, R)))
        return B

    def pencil(self, kvec):
        B = self.basis(kvec)
        M = HC.bloch(self.P, self.bars, self.A, kvec)
        D = M.conj().T @ M
        return B.conj().T @ D @ B, B.conj().T @ B

    def bands(self, kvec):
        q, _ = np.linalg.qr(self.basis(kvec))
        M = HC.bloch(self.P, self.bars, self.A, kvec)
        D = M.conj().T @ M
        w = np.linalg.eigvalsh(q.conj().T @ D @ q)
        return np.sqrt(np.maximum(w, 0.0))

    def split(self, kvec):
        b = self.bands(kvec)
        return float(b[1] - b[0])

    def structure_factor(self, kvec):
        return 8.0 * float(np.prod(np.cos(np.asarray(kvec) * self.L)))

    def phase_direction(self, col=0):
        """The single phase direction at Gamma, in node space.

        At Gamma the two sublattice columns are PARALLEL -- a shared vertex is
        driven along one line by both cells that own it -- so the basis is rank
        one and the 2x2 projection there is ill-posed. This returns the one
        real direction, and B8 gates the rank deficiency rather than hiding it."""
        v = np.zeros(3 * self.n)
        for (ci, i, _R, dp) in self.stems:
            if ci == col:
                v[3 * i:3 * i + 3] += dp
        return v

    def gamma_cost(self, col=0):
        """Rayleigh quotient of that direction: the phase mode's energy at
        Gamma. No QR, no 2x2, no noise column -- this is the statistic B7 is
        gated on, and it is invariant under everything that moved 1.486784."""
        v = self.phase_direction(col)
        M = HC.bloch(self.P, self.bars, self.A, np.zeros(3), unit=True).real
        return float(np.linalg.norm(M @ v) / np.linalg.norm(v))

    def basis_cond(self, kvec):
        sv = np.linalg.svd(self.basis(kvec), compute_uv=False)
        return float(sv[-1] / sv[0])


def reduced_bands(model, kvec):
    """The REDUCED (splitting-vertex) model's closed form, T2 23484. Present
    only so B3 can contrast against a computed curve instead of a quotation."""
    f = abs(model.structure_factor(kvec))
    return (np.sqrt(ALPHA_REDUCED * (8.0 - f)),
            np.sqrt(ALPHA_REDUCED * (8.0 + f)))


def inplane_form(model, k2, k3):
    """B2's closed form for the degenerate value on the plane k_1 = pi/A."""
    s2, s3 = np.sin(k2 * model.L) ** 2, np.sin(k3 * model.L) ** 2
    c2, c3 = np.cos(k2 * model.L) ** 2, np.cos(k3 * model.L) ** 2
    den = 3.0 * (c2 + c3)
    if den < 1e-14:                      # the face corner, R: 0/0, removable
        return float("nan")
    return 4.0 * (1.0 - s2 * s3) / den


# ==========================================================================
# B0: THE PRECOMPUTED BASIS AGREES WITH jb_hc's OWN
# ==========================================================================

def b0_agreement():
    m = Model(A_REF)
    worst = 0.0
    pts = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0),
           (0.41, 0.67, 0.29), (1.0, 0.37, 0.61)]
    for v in pts:
        kv = m.G * np.array(v, dtype=float)
        mine = m.bands(kv)
        q = HC.phase_basis(m.P, m.slots, m.A, kv)
        M = HC.bloch(m.P, m.bars, m.A, kv)
        D = M.conj().T @ M
        theirs = np.sqrt(np.maximum(np.linalg.eigvalsh(q.conj().T @ D @ q), 0.0))
        worst = max(worst, float(np.max(np.abs(mine - theirs))))
    return dict(worst=worst, n=len(pts), model=m)


# ==========================================================================
# B1: THE LOCUS. Degenerate on the six zone-boundary planes and nowhere else.
# ==========================================================================

def b1_locus(m):
    on = []
    for axis in range(3):
        for sgn in (1.0, -1.0):
            for grid in (PLANE_GRID, PLANE_GRID_ALT):
                for t in grid:
                    for u in grid:
                        v = [t, u]
                        kv = np.empty(3)
                        kv[axis] = sgn
                        kv[(axis + 1) % 3] = v[0]
                        kv[(axis + 2) % 3] = v[1]
                        on.append(m.split(m.G * kv))
    matched = [m.split(m.G * np.array([t, 0.37, 0.61]))
               for t in (0.2, 0.3, 0.5, 0.7, 0.9)]
    rng = np.random.default_rng(INTERIOR_SEED)
    interior, small = [], 0
    for p in rng.uniform(-1.0, 1.0, size=(INTERIOR_N, 3)):
        if float(np.max(np.abs(np.abs(p) - 1.0))) < 1e-9:
            continue
        s = m.split(m.G * p)
        interior.append(s)
        if s < TOL_DEGEN:
            small += 1
    return dict(on_max=max(on) if on else float("nan"), on_n=len(on),
                matched_min=min(matched) if matched else float("nan"),
                interior_min=min(interior) if interior else float("nan"),
                interior_n=len(interior), interior_small=small)


# ==========================================================================
# B2: THE IN-PLANE DISPERSION, in closed form.
# ==========================================================================

def b2_inplane(m):
    resid, vals = [], []
    for grid in (PLANE_GRID, PLANE_GRID_ALT):
        for t in grid:
            for u in grid:
                kv = m.G * np.array([1.0, t, u])
                w2 = float(m.bands(kv)[0]) ** 2
                pred = inplane_form(m, m.G * t, m.G * u)
                vals.append(w2)
                if np.isfinite(pred):
                    resid.append(abs(pred - w2))
    centre = float(m.bands(m.G * np.array([1.0, 0.0, 0.0]))[0])
    edge = float(m.bands(m.G * np.array([1.0, 1.0, 0.0]))[0])
    wrong = [abs(v - 8.0 / 9.0) for v in vals]      # a deliberately wrong form
    return dict(resid=max(resid) if resid else float("nan"), n=len(resid),
                span=float(np.sqrt(max(vals)) - np.sqrt(min(vals))),
                centre=centre, edge=edge, ratio=edge / centre,
                wrong_max=max(wrong) if wrong else float("nan"))


# ==========================================================================
# B3: THE REDUCED MODEL IS FLAT ON THE PLANE. Computed, not quoted.
# ==========================================================================

def b3_reduced_contrast(m):
    red, faithful = [], []
    for t in PLANE_GRID:
        for u in PLANE_GRID:
            kv = m.G * np.array([1.0, t, u])
            red.append(reduced_bands(m, kv)[0])
            faithful.append(float(m.bands(kv)[0]))
    rs = float(np.ptp(red)) / float(np.mean(red))
    fs = float(np.ptp(faithful)) / float(np.mean(faithful))
    return dict(reduced_rel_span=rs, faithful_rel_span=fs,
                ratio=(fs / rs) if rs > 0 else float("inf"), n=len(red))


# ==========================================================================
# B4: GROUP VELOCITY, AND THE THREE-WAY VERDICT.
# ==========================================================================

def b4_velocity(m, h=1e-6):
    def band(kv, br):
        return float(m.bands(kv)[br])

    off_centre = [(0.37, 0.61), (0.5, 0.0), (0.25, 0.75), (0.9, 0.1),
                  (-0.43, 0.28)]
    vin, vn, mismatch = [], [], []
    for (t, u) in off_centre:
        kv = m.G * np.array([1.0, t, u])
        g = []
        for d in (1, 2):
            e = np.zeros(3)
            e[d] = h
            g.append((band(kv + e, 0) - band(kv - e, 0)) / (2 * h))
        vin.append(float(np.linalg.norm(g)))
        e = np.array([h, 0.0, 0.0])
        lo0, hi0 = m.bands(kv)
        slo = abs(float(lo0) - band(kv - e, 0)) / h
        shi = abs(band(kv - e, 1) - float(hi0)) / h
        vn.append(0.5 * (slo + shi))
        mismatch.append(abs(slo - shi))

    kc = m.G * np.array([1.0, 0.0, 0.0])
    gc = []
    for d in (1, 2):
        e = np.zeros(3)
        e[d] = h
        gc.append((band(kc + e, 0) - band(kc - e, 0)) / (2 * h))
    centre_v = float(np.linalg.norm(gc))

    vin_min = min(vin) if vin else float("nan")
    vn_min = min(vn) if vn else float("nan")
    computable = bool(np.isfinite(vin_min) and np.isfinite(vn_min))
    if not computable:
        verdict = "NOT-COMPUTABLE"
    elif vin_min < VZERO and vn_min < VZERO:
        verdict = "STALL"
    elif vin_min < VZERO:
        verdict = "CHANNEL"
    else:
        verdict = "CROSS"
    return dict(vin_min=vin_min, vin_max=max(vin), vn_min=vn_min,
                vn_max=max(vn), mismatch=max(mismatch), centre_v=centre_v,
                verdict=verdict, computable=computable, n=len(off_centre))


# ==========================================================================
# B5: THE MECHANISM -- equal diagonals AND a vanishing structure factor.
# ==========================================================================

def b5_mechanism(m):
    rng = np.random.default_rng(INTERIOR_SEED + 1)
    ratios, diag = [], []
    for p in rng.uniform(-1.0, 1.0, size=(24, 3)):
        kv = m.G * p
        D, S = m.pencil(kv)
        f = abs(m.structure_factor(kv))
        diag.append(abs(D[0, 0] - D[1, 1]) + abs(S[0, 0] - S[1, 1]))
        if f > 1e-3:
            ratios.append(abs(S[0, 1]) / f)
    spread = (float(np.std(ratios) / np.mean(ratios))
              if ratios else float("nan"))
    on_off = [abs(m.pencil(m.G * np.array([1.0, t, u]))[0][0, 1])
              for t in (0.37, -0.21) for u in (0.61, -0.53)]
    off_off = [abs(m.pencil(m.G * np.array([s, 0.37, 0.61]))[0][0, 1])
               for s in (0.2, 0.5, 0.8)]
    other = Model(-20.0)
    kv2 = other.G * np.array([1.0, 0.37, 0.61])
    D2, S2 = other.pencil(kv2)
    return dict(prop_spread=spread, prop_n=len(ratios),
                diag_max=max(diag) if diag else float("nan"), diag_n=len(diag),
                onplane_off=max(on_off), offplane_off=min(off_off),
                other_diag=abs(D2[0, 0] - D2[1, 1]) + abs(S2[0, 0] - S2[1, 1]),
                other_offdiag=abs(D2[0, 1]))


# ==========================================================================
# B6: a0 SPECIFICITY -- the symmetry, and how fast it dies.
# ==========================================================================

def _point_ops():
    ops = []
    for perm in it.permutations(range(3)):
        for sg in it.product((1, -1), repeat=3):
            M = np.zeros((3, 3))
            for i, p in enumerate(perm):
                M[i, p] = sg[i]
            ops.append(M)
    return ops


def _hausdorff(X, Y):
    return max(max(float(np.linalg.norm(Y - x, axis=1).min()) for x in X),
               max(float(np.linalg.norm(X - y, axis=1).min()) for y in Y))


def b6_phase_specificity():
    m0 = Model(A_REF)
    at_ref = m0.split(m0.G * np.array([1.0, 0.37, 0.61]))
    away, away_alt = [], []
    for a in PHASE_SWEEP:
        mm = Model(a)
        away.append(mm.split(mm.G * np.array([1.0, 0.37, 0.61])))
    for a in PHASE_SWEEP_ALT:
        mm = Model(a)
        away_alt.append(mm.split(mm.G * np.array([1.0, 0.37, 0.61])))
    swap = []
    for a in (-20.0, -25.0, -10.0):
        ma, mb = Model(a), Model(-60.0 - a)
        swap.append(abs(ma.split(ma.G * np.array([1.0, 0.37, 0.61]))
                        - mb.split(mb.G * np.array([1.0, 0.37, 0.61]))))
    d, s = [], []
    for da in DELTA_SWEEP:
        mm = Model(A_REF + da)
        d.append(da)
        s.append(mm.split(mm.G * np.array([1.0, 0.37, 0.61])))
    expo = float(np.polyfit(np.log(d), np.log(s), 1)[0])

    ops = _point_ops()
    VE = np.unique(np.round(Z.corners(A_REF).reshape(-1, 3), 9), axis=0)
    OC = np.unique(np.round(Z.corners(A_REF + PHASE_OFFSET).reshape(-1, 3), 9),
                   axis=0)
    maps = [M for M in ops if _hausdorff(OC, VE @ M.T) < 1e-9]
    selfm = [M for M in ops if _hausdorff(VE, VE @ M.T) < 1e-9]
    VE2 = np.unique(np.round(Z.corners(-20.0).reshape(-1, 3), 9), axis=0)
    OC2 = np.unique(np.round(Z.corners(-20.0 + PHASE_OFFSET).reshape(-1, 3), 9),
                    axis=0)
    maps2 = [M for M in ops if _hausdorff(OC2, VE2 @ M.T) < 1e-9]

    L = HC.lattice(A_REF)
    A = 2.0 * L
    tau = L * np.ones(3)
    pts = []
    for c in it.product(range(-2, 3), repeat=3):
        o = A * np.array(c, dtype=float)
        pts.append(Z.corners(A_REF).reshape(-1, 3) + o)
        pts.append(Z.corners(A_REF + PHASE_OFFSET).reshape(-1, 3) + o + tau)
    Pt = np.unique(np.round(np.vstack(pts), 8), axis=0)
    core = Pt[np.max(np.abs(Pt), axis=1) < 1.5 * A]
    worst_sym = 0.0
    for M in maps[:6]:
        img = core @ M.T + tau
        worst_sym = max(worst_sym,
                        max(float(np.linalg.norm(Pt - p, axis=1).min())
                            for p in img))
    return dict(at_ref=at_ref, away_min=min(away), away_alt_min=min(away_alt),
                swap_max=max(swap), expo=expo, n_maps=len(maps),
                n_self=len(selfm), n_ops=len(ops), n_maps_other=len(maps2),
                sym_resid=worst_sym, translation_haus=_hausdorff(OC, VE))


# ==========================================================================
# B7: THE GAMMA GAP, AND THE LATTICE THAT IS NOT ALLOWED TO BREATHE.
# ==========================================================================

def _dlambda(a):
    d1 = (HC.lattice(a + 1e-3) - HC.lattice(a - 1e-3)) / 2e-3
    d2 = (HC.lattice(a + 5e-4) - HC.lattice(a - 5e-4)) / 1e-3
    return (4.0 * d2 - d1) / 3.0 * (180.0 / np.pi)


def b7_gap():
    rows = []
    for a in PHASE_SWEEP + PHASE_SWEEP_ALT:
        mm = Model(a)
        cost = mm.gamma_cost()
        dl = _dlambda(a)
        rows.append((a, cost, dl,
                     cost / (np.sqrt(6.0) * abs(dl)) if abs(dl) > 1e-12
                     else float("nan")))
    m0 = Model(A_REF)
    # both columns must agree: at Gamma they are parallel, so this is a
    # consistency check on the rank-one claim B8 gates
    col_gap = max(abs(Model(a).gamma_cost(0) - Model(a).gamma_cost(1))
                  for a in (-20.0, -10.0, -41.3))
    # the two controls that separate a real number from the 1.486784 class
    h_spread = max(
        max(Model(a, h=hh).gamma_cost() for hh in H_SWEEP)
        - min(Model(a, h=hh).gamma_cost() for hh in H_SWEEP)
        for a in (-20.0, -10.0))
    rev_diff = max(abs(Model(a).gamma_cost() - Model(a, reverse=True).gamma_cost())
                   for a in (-20.0, -10.0, -5.0))
    iso = [float(Model(-20.0).bands(1e-6 * Model(-20.0).G
                                    * (np.array(d, dtype=float)
                                       / np.linalg.norm(d)))[0])
           for d in ((1, 0, 0), (1, 1, 0), (1, 1, 1))]
    eps = (0.02, 0.005, 0.001)
    lin = [float(m0.bands(e * m0.G * np.array([1.0, 0.0, 0.0]))[0]) for e in eps]
    agree = max(abs(Model(a).gamma_cost()
                    - float(Model(a).bands(1e-8 * Model(a).G
                                           * np.array([1.0, 0.0, 0.0]))[0]))
                for a in (-20.0, -10.0))
    return dict(rows=rows, worst=max(abs(r[3] - 1.0) for r in rows),
                gap0=m0.gamma_cost(), dl0=_dlambda(A_REF),
                lam0=HC.lattice(A_REF), iso_spread=max(iso) - min(iso),
                lin_ratio=lin[-1] / lin[0], col_gap=col_gap,
                h_spread=h_spread, rev_diff=rev_diff, agree=agree)


# ==========================================================================
# B8: THE GAMMA POINT IS RANK ONE, and the "optical gap" recorded there was
#     an accumulation-order artifact. Independent confirmation of B7.
# ==========================================================================

def b8_gamma_rank():
    conds, mm = [], None
    for a in (A_REF, -20.0, -10.0, -41.3):
        mm = Model(a)
        conds.append(mm.basis_cond(np.zeros(3)))
    m0 = Model(A_REF)
    cond_M = m0.basis_cond(m0.G * np.array([1.0, 1.0, 0.0]))
    # the upper branch: its LIMIT is well defined, its value AT Gamma is not
    lim = [float(m0.bands(e * m0.G * np.array([1.0, 0.0, 0.0]))[1])
           for e in (1e-3, 1e-4, 1e-5)]
    at0 = float(m0.bands(np.zeros(3))[1])
    target = 2.0 / np.sqrt(3.0)
    # CONTROL, TWO SIDED. Detune one node's octa entry so the basis stops
    # being rank one, and sweep the detuning ACROSS the finite-difference
    # roundoff floor (cond ~ 3.3e-10 here). Above the floor the second band is
    # a stable number; below it, it drifts back to the unperturbed value. That
    # crossover is the demonstration that the unperturbed value is carried by
    # roundoff -- without it, "the value at Gamma is noise" would be
    # indistinguishable from "the value at Gamma is a number we mis-predicted".
    def band1(d):
        q, _ = np.linalg.qr(_detuned_basis(m0, d))
        M = HC.bloch(m0.P, m0.bars, m0.A, np.zeros(3))
        D = M.conj().T @ M
        return float(np.sqrt(max(np.linalg.eigvalsh(q.conj().T @ D @ q)[1], 0.0)))
    strong = [band1(d) for d in DETUNE_STRONG]
    weak = [band1(d) for d in DETUNE_WEAK]
    return dict(cond_gamma=max(conds), cond_M=cond_M, n=len(conds),
                lim_err=max(abs(x - target) for x in lim),
                at0_err=abs(at0 - target), at0=at0,
                strong_spread=max(strong) - min(strong), strong=strong[-1],
                weak_far=min(abs(x - strong[-1]) for x in weak),
                weak_near=max(abs(x - at0) for x in weak))


def _detuned_basis(m, d):
    """m's Gamma basis with ONE node's octa entry rescaled, so the two columns
    stop being parallel and the basis is honestly rank two."""
    B = m.basis(np.zeros(3))
    B[0:3, 1] *= (1.0 + d)
    return B


# ==========================================================================
# THE GATE
# ==========================================================================

def gate(b0, b1, b2, b3, b4, b5, b6, b7, b8):
    checks = []
    R = checks.append

    R(("B0  precomputed phase basis reproduces jb_hc's own to machine zero",
       b0["n"] > 0 and b0["worst"] < TOL_EXACT,
       f"worst {b0['worst']:.2e} over {b0['n']} k", f"< {TOL_EXACT:.0e}"))

    R(("B1  bands DEGENERATE on all six zone-boundary planes",
       b1["on_n"] > 0 and b1["on_max"] < TOL_DEGEN,
       f"worst {b1['on_max']:.2e} over {b1['on_n']} k", f"< {TOL_DEGEN:.0e}"))
    R(("B1  CONTROL: SPLIT off the planes at matched |k| -- CAN FAIL",
       b1["matched_min"] > SPLIT_OFF_MIN,
       f"min {b1['matched_min']:.4f}", f"> {SPLIT_OFF_MIN}"))
    R(("B1  the locus is ONLY the boundary: no interior k is degenerate",
       b1["interior_n"] > 0 and b1["interior_small"] == 0,
       f"{b1['interior_small']}/{b1['interior_n']} below tol, "
       f"min {b1['interior_min']:.2e}", "0 of them"))

    R(("B2  in-plane CLOSED FORM matches the projection",
       b2["n"] > 0 and b2["resid"] < TOL_FORM,
       f"worst {b2['resid']:.2e} over {b2['n']} k", f"< {TOL_FORM:.0e}"))
    R(("B2  CONTROL: a deliberately wrong constant form does NOT match "
       "-- CAN FAIL", b2["wrong_max"] > 1e-2,
       f"worst {b2['wrong_max']:.4f}", "> 1e-2"))
    R(("B2  face centre = sqrt(2/3), face edge = 2/sqrt(3), ratio EXACTLY "
       "sqrt(2)",
       abs(b2["centre"] - np.sqrt(2 / 3)) < 1e-6
       and abs(b2["edge"] - 2 / np.sqrt(3)) < 1e-6
       and abs(b2["ratio"] - np.sqrt(2)) < 1e-6,
       f"{b2['centre']:.9f} / {b2['edge']:.9f} / {b2['ratio']:.9f}",
       "0.816497 / 1.154701 / 1.414214"))
    R(("B2  NON-VACUITY, TWO-SIDED: the in-plane span is neither stuck nor "
       "divergent",
       np.isfinite(b2["span"]) and SPAN_BAND[0] < b2["span"] < SPAN_BAND[1],
       f"span {b2['span']:.6f}", f"in {SPAN_BAND}"))

    R(("B3  the REDUCED model really is FLAT on the plane (computed, not "
       "quoted)",
       b3["n"] > 0 and b3["reduced_rel_span"] < 1e-12,
       f"rel span {b3['reduced_rel_span']:.2e}", "< 1e-12"))
    R(("B3  the FAITHFUL model is not -- the two disagree by orders of "
       "magnitude",
       b3["faithful_rel_span"] > 0.1 and b3["ratio"] > 1e6,
       f"faithful {b3['faithful_rel_span']:.4f}, ratio {b3['ratio']:.2e}",
       "> 0.1, > 1e6"))

    R(("B4  IN-PLANE group velocity is NONZERO off the face centre",
       b4["n"] > 0 and b4["vin_min"] > VIN_MIN,
       f"min {b4['vin_min']:.6f} max {b4['vin_max']:.6f}", f"> {VIN_MIN}"))
    R(("B4  CONTROL: at the face CENTRE it vanishes, so B4 is not vacuous "
       "-- CAN FAIL", b4["centre_v"] < VZERO,
       f"{b4['centre_v']:.2e}", f"< {VZERO:.0e}"))
    R(("B4  NORMAL group velocity is nonzero: the bands CROSS, not touch",
       b4["vn_min"] > VN_MIN, f"min {b4['vn_min']:.6f}", f"> {VN_MIN}"))
    R(("B4  CROSSING not REPULSION: the two slopes match in magnitude",
       b4["mismatch"] < TOL_SLOPE, f"worst {b4['mismatch']:.2e}",
       f"< {TOL_SLOPE:.0e}"))
    R(("B4  VERDICT is computable and is CROSS -- STALL and CHANNEL were both "
       "reachable", b4["computable"] and b4["verdict"] == "CROSS",
       b4["verdict"], "CROSS"))

    R(("B5  off-diagonal IS the inter-sublattice structure factor "
       "8*prod cos(k L)",
       b5["prop_n"] > 0 and b5["prop_spread"] < TOL_PROP,
       f"rel spread {b5['prop_spread']:.2e} over {b5['prop_n']} k",
       f"< {TOL_PROP:.0e}"))
    R(("B5  the pencil's DIAGONALS are equal at a0 = -30, at every k",
       b5["diag_n"] > 0 and b5["diag_max"] < 1e-9,
       f"worst {b5['diag_max']:.2e}", "< 1e-9"))
    R(("B5  CONTROL: off-diagonal is LARGE off the planes -- CAN FAIL",
       b5["offplane_off"] > 1e-5 and b5["onplane_off"] < 1e-9,
       f"off {b5['offplane_off']:.2e} vs on {b5['onplane_off']:.2e}",
       "> 1e-5 vs < 1e-9"))
    R(("B5  CONTROL: at a0 = -20 the DIAGONALS differ while the off-diagonal "
       "still vanishes -- so it is the diagonals that a0 = -30 buys",
       b5["other_diag"] > 1e-6 and b5["other_offdiag"] < 1e-9,
       f"diag {b5['other_diag']:.2e}, offdiag {b5['other_offdiag']:.2e}",
       "> 1e-6, < 1e-9"))

    R(("B6  the degeneracy exists at a0 = -30 and NOWHERE ELSE on the exchange",
       b6["at_ref"] < TOL_DEGEN and b6["away_min"] > 1e-3
       and b6["away_alt_min"] > 1e-3,
       f"ref {b6['at_ref']:.2e}, away >= {min(b6['away_min'], b6['away_alt_min']):.4f}",
       "< 1e-8 vs > 1e-3"))
    R(("B6  it opens LINEARLY in |a0 + 30| (TWO-SIDED band on the exponent)",
       EXP_BAND[0] < b6["expo"] < EXP_BAND[1],
       f"exponent {b6['expo']:.6f}", f"in {EXP_BAND}"))
    R(("B6  ROLE-SWAP: the spectrum at a and at -60-a is identical",
       b6["swap_max"] < 1e-9, f"worst {b6['swap_max']:.2e}", "< 1e-9"))
    R(("B6  the two sublattices ARE equivalent at a0 = -30: half the "
       "octahedral group carries one onto the other",
       b6["n_maps"] == 24 and b6["n_self"] == 24
       and b6["n_maps"] + b6["n_self"] == b6["n_ops"],
       f"{b6['n_maps']} map, {b6['n_self']} fix, of {b6['n_ops']}", "24 + 24 = 48"))
    R(("B6  CONTROL: at a0 = -20 NO octahedral operation carries one onto the "
       "other -- CAN FAIL", b6["n_maps_other"] == 0,
       f"{b6['n_maps_other']} ops", "0"))
    R(("B6  {g | (A/2)(1,1,1)} really is a symmetry of the whole honeycomb",
       b6["sym_resid"] < 1e-7, f"worst {b6['sym_resid']:.2e}", "< 1e-7"))
    R(("B6  and it is NOT bcc folding: no PURE translation does it -- CAN FAIL",
       b6["translation_haus"] > 1e-3,
       f"Hausdorff {b6['translation_haus']:.6f}", "> 1e-3"))

    R(("B7  the Gamma gap is EXACTLY sqrt(6)*|d lambda/d a| at every phase",
       len(b7["rows"]) > 0 and b7["worst"] < TOL_GAP,
       f"worst |ratio-1| {b7['worst']:.2e} over {len(b7['rows'])} phases",
       f"< {TOL_GAP:.0e}"))
    R(("B7  CONTROL: invariant across SIX DECADES of finite-difference step "
       "-- the test that exposes an accumulation artifact",
       b7["h_spread"] < TOL_GAP, f"spread {b7['h_spread']:.2e}",
       f"< {TOL_GAP:.0e}"))
    R(("B7  CONTROL: invariant under REVERSED slot order -- the exact "
       "perturbation that produced the retracted 1.486784",
       b7["rev_diff"] < 1e-15, f"diff {b7['rev_diff']:.2e}", "< 1e-15"))
    R(("B7  both sublattice columns give the same gap (they are parallel at "
       "Gamma, per B8)", b7["col_gap"] < 1e-12,
       f"worst {b7['col_gap']:.2e}", "< 1e-12"))
    R(("B7  the well-conditioned Rayleigh value AGREES with the projected "
       "band near Gamma", b7["agree"] < 1e-8, f"worst {b7['agree']:.2e}",
       "< 1e-8"))
    R(("B7  lambda is MAXIMAL at a0 = -30, exactly 2/sqrt(3), so the gap "
       "closes there",
       abs(b7["dl0"]) < 1e-9 and abs(b7["lam0"] - 2 / np.sqrt(3)) < TOL_EXACT
       and b7["gap0"] < 1e-4,
       f"dlam/da {b7['dl0']:.2e}, lambda {b7['lam0']:.9f}, gap {b7['gap0']:.2e}",
       "0, 1.154701, 0"))
    R(("B7  at a0 = -30 the branch is LINEAR in |k| (gapless), not flat "
       "-- CAN FAIL",
       b7["lin_ratio"] < 0.2, f"ratio(1e-3/2e-2) {b7['lin_ratio']:.4f}",
       "< 0.2 (0.05 = linear)"))
    R(("B7  the gap is ISOTROPIC, as a Gamma-point gap must be",
       b7["iso_spread"] < 1e-9, f"spread {b7['iso_spread']:.2e}", "< 1e-9"))

    R(("B8  the phase basis is RANK ONE at Gamma, at every phase -- the two "
       "sublattice columns are parallel",
       b8["n"] > 0 and b8["cond_gamma"] < 1e-8,
       f"worst cond {b8['cond_gamma']:.2e} over {b8['n']} phases", "< 1e-8"))
    R(("B8  CONTROL: on the zone boundary the same basis is PERFECTLY "
       "conditioned -- CAN FAIL", abs(b8["cond_M"] - 1.0) < 1e-9,
       f"cond {b8['cond_M']:.9f}", "1.0"))
    R(("B8  so the upper branch has a LIMIT (2/sqrt(3)) but no value at "
       "Gamma: they disagree",
       b8["lim_err"] < 1e-5 and b8["at0_err"] > 1e-2,
       f"lim err {b8['lim_err']:.2e}, at k=0 {b8['at0']:.6f} "
       f"(err {b8['at0_err']:.3f})", "< 1e-5 vs > 1e-2"))
    R(("B8  CONTROL: detuned WELL ABOVE the roundoff floor, the second band "
       "is a stable number over four decades -- CAN FAIL",
       b8["strong_spread"] < 1e-5,
       f"spread {b8['strong_spread']:.2e} about {b8['strong']:.9f}", "< 1e-5"))
    R(("B8  CONTROL, other side: detuned BELOW the floor it drifts back to "
       "the unperturbed value, not to the stable one -- this crossover is "
       "what proves the Gamma value is roundoff",
       b8["weak_near"] < 1e-2 and b8["weak_far"] > 1e-1,
       f"near unperturbed {b8['weak_near']:.2e}, "
       f"far from stable {b8['weak_far']:.3f}", "< 1e-2 and > 0.1"))

    print()
    print("=" * 78)
    print(f"GATE  {len(checks)} rows")
    print("=" * 78)
    for name, ok, val, crit in checks:
        print(f"  {'PASS' if ok else 'FAIL':4s}  {name:66s} {str(val):>26s} {str(crit):>18s}")

    print()
    print("  ROWS THAT EXIST ONLY TO STOP ANOTHER ROW BEING UNFALSIFIABLE:")
    print("   * B1's MATCHED-|k| control. Without it 'degenerate on the")
    print("     planes' is satisfiable by a projection that has collapsed to")
    print("     one band everywhere, which is exactly what a rank-deficient")
    print("     phase basis would produce.")
    print("   * B1's INTERIOR scan. On-plane degeneracy proves the planes are")
    print("     IN the locus, never that they are ALL of it.")
    print("   * B2's wrong-form control and its TWO-SIDED span band. A stuck")
    print("     sampler returning one k repeatedly would give a perfect fit")
    print("     to any form and a span of zero; the band's lower edge fails")
    print("     that, and the upper edge fails a divergent statistic.")
    print("   * B3's REDUCED-model row. The claim being overturned is that")
    print("     the plane is flat, so that claim is computed here rather than")
    print("     quoted -- otherwise B2 is contradicting a citation, not a")
    print("     model.")
    print("   * B4's FACE-CENTRE control. Without it 'the in-plane velocity")
    print("     is nonzero' could be a finite-difference floor rather than")
    print("     physics; at the face centre the same statistic must come out")
    print("     zero, and it does.")
    print("   * B5's a0 = -20 control, which separates the two conditions:")
    print("     the structure factor still vanishes there, and the")
    print("     degeneracy is gone anyway, so the diagonals are what the")
    print("     reference phase buys.")
    print("   * B6's no-pure-translation row. Without it the symmetry finding")
    print("     reads as ordinary bcc band folding, which would make the")
    print("     locus robust for a reason that is not true.")
    print("   * B7's LINEARITY row. 'gap = 0 at a0 = -30' is satisfied by a")
    print("     branch that is flat at zero as well as by a gapless linear")
    print("     one; only the ratio test tells them apart.")
    print("   * B7's h-SWEEP and REVERSED-SLOT-ORDER controls, and B8's")
    print("     detuned rank-two control. These are not generic hygiene. The")
    print("     record's optical gap 1.486784 was a real computation of an")
    print("     ungated quantity whose value was carried by floating-point")
    print("     accumulation order; reversing the slot iteration reproduces")
    print("     it exactly. A number that survives six decades of step size")
    print("     and an order reversal is a different kind of number, and")
    print("     these rows are what says which kind B7 has.")
    print()
    print("  ROWS DELETED RATHER THAN FIXED: a row asserting that the")
    print("  degenerate value is CONSTANT on each plane at 2/sqrt(3). It was")
    print("  written from T2's M and R readings before the scan, and it is")
    print("  false: those two points sit on the face EDGES, where the value")
    print("  is 2/sqrt(3), while the face CENTRE is sqrt(2/3). Two agreeing")
    print("  samples had been generalised to a surface. B2 replaces it.")
    print()
    print("  A ROW DELIBERATELY NOT BUILT: the real-space packet experiment.")
    print("  B4 settles the harmonic question -- there is no zero-velocity")
    print("  locus to localise on -- but 'what does a packet launched near a")
    print("  plane actually do' needs a time-domain simulation and the")
    print("  anharmonic terms, neither of which exists. Nothing here licenses")
    print("  a claim about transient icosahedral patterns in either")
    print("  direction.")
    print()
    print("  WHAT THIS FILE DOES NOT MODEL, AND WHY IT NOW MATTERS MORE: the")
    print("  lattice is still held at lambda(a0) and not allowed to breathe.")
    print("  B7 turns that from a caveat into the leading term -- the phase")
    print("  band's gap IS the breathing that was forbidden, sqrt(6) times")
    print("  |d lambda/d a|. Every wave number this epic has was measured at")
    print("  the one phase where that term vanishes. Anharmonic terms are")
    print("  still absent, and inviscid-l1d is still unfixed and still")
    print("  independent of all of this.")

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
    print("jb_bt -- the band-touching planes, re-scanned in the faithful model")
    print("=" * 78)
    print("  The reduced model put the two phase bands together wherever any")
    print("  k_i = pi/(2 lambda), and the record read that as a localisation")
    print("  mechanism. The locus survives the move to the faithful model --")
    print("  it is the simple-cubic zone boundary, and it is exact -- but the")
    print("  degenerate value DISPERSES across each plane, the bands cross")
    print("  linearly rather than touching, and no group velocity vanishes.")
    print("  There is no stall. And the whole structure belongs to a0 = -30")
    print("  alone: it is the phase where the two cells are congruent, and")
    print("  where lambda is stationary so the fixed lattice costs nothing.")
    b0 = b0_agreement()
    m = b0.pop("model")
    b1 = b1_locus(m)
    b2 = b2_inplane(m)
    b3 = b3_reduced_contrast(m)
    b4 = b4_velocity(m)
    b5 = b5_mechanism(m)
    b6 = b6_phase_specificity()
    b7 = b7_gap()
    b8 = b8_gamma_rank()

    print()
    print("-" * 78)
    print("  THE PLANE, IN NUMBERS")
    print("-" * 78)
    print(f"  face centre (X)   omega = {b2['centre']:.9f}   = sqrt(2/3)")
    print(f"  face edge  (M,R)  omega = {b2['edge']:.9f}   = 2/sqrt(3)")
    print(f"  ratio                   {b2['ratio']:.9f}   = sqrt(2)")
    print(f"  in-plane span           {b2['span']:.9f}")
    print(f"  in-plane |v|            {b4['vin_min']:.6f} .. {b4['vin_max']:.6f}"
          f"   (0 at the centre)")
    print(f"  normal   |v|            {b4['vn_min']:.6f} .. {b4['vn_max']:.6f}")
    print(f"  VERDICT                 {b4['verdict']}")
    print()
    print("  THE GAMMA GAP vs THE FORBIDDEN BREATHING")
    print(f"    {'a0':>8s} {'d lambda/d a':>16s} {'gap':>14s} "
          f"{'gap/(sqrt6*|dl|)':>20s}")
    for (a, gap, dl, r) in b7["rows"]:
        print(f"    {a:8.1f} {dl:+16.9f} {gap:14.9f} {r:20.10f}")
    print(f"    {A_REF:8.1f} {b7['dl0']:+16.9f} {b7['gap0']:14.9f} "
          f"{'(0/0: gapless)':>20s}")
    return gate(b0, b1, b2, b3, b4, b5, b6, b7, b8)


if __name__ == "__main__":
    with np.errstate(all="ignore"):
        sys.exit(main())
