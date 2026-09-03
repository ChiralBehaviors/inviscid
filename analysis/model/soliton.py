"""soliton -- the medium's NONLINEARITY, and the wave that survives it.

THE QUESTION (owner, 2026-09-02): can this medium carry a soliton?

WHY IT CANNOT BE ANSWERED BY ANYTHING MEASURED SO FAR. Every lane to date
is linearised at the working fold: the kick's scalar wave, the shove and
stir's polarized waves. In a linear medium superposition is exact --
kick_response's R10 gates it at 2e-16 -- and exact superposition is
precisely the statement that no soliton exists. A soliton is nonlinearity
balancing dispersion, so answering the question means building the
nonlinear soft-joint dynamics that T2 [23959] named as the natural next
measurement, and that is what this module is.

WHAT THE NONLINEARITY IS, AND WHY IT IS WHERE IT IS. `assembly.
lattice_constant` is STATIONARY at a = -30: dL/da = 0 there, and it is
the only angle where a coherently breathing patch's centres do not move.
That one fact runs through everything here. It is why the fold sector
decouples from translation at LINEAR order (kick_response R5), and it is
why the leading coupling between them is QUADRATIC: a fold deviation of
either sign contracts the lattice by (1/2)L''d^2, so fold amplitude
squared is a source of longitudinal strain. The medium's nonlinearity is
not a small correction bolted onto the fold law; it is the second-order
remainder of the same geometric fact that makes the fold law clean.

THE EXACT REDUCTION (R4). A PLANE-SYMMETRIC state of the bulk medium --
one that depends only on the x index -- reduces to ONE cell per plane
with three bonds: an x bond to the next plane, and two SELF bonds tying
the cell to its own transverse images at fixed lattice offset (this is
dispersion.py's "one cell per primitive cell, 7 DOF, 3 bonds" with the x
bond made non-local along the chain). Of the seven coordinates, only TWO
move: R4 measures that from any state with cells unrotated and displaced
only along x, the full seven-coordinate force has zero rotation and zero
transverse-translation components, to 1e-15. So the plane-symmetric
medium is EXACTLY a two-field lattice -- longitudinal displacement u_n
and fold gamma_n -- and this module integrates it exactly, with the true
joint separations (`weld_residual`'s own quantity), not a linearisation.

TWO EXACT STRUCTURAL FACTS, BOTH SURPRISING.

  * THE TWO FIELDS ARE THE SAME HARMONIC CHAIN (R5). Onsite 4, coupling
    -2, mass 8, for BOTH: omega = sin(q/2) exactly, sound speed 1/2, and
    the two branches are exactly degenerate. That degeneracy is the one
    shear_response R2 reports as "the longitudinal translation is
    degenerate with the fold at 1/2", here in closed form.
  * A PLANE LONGITUDINAL WAVE IS EXACTLY HARMONIC, AT EVERY AMPLITUDE
    (R7). Zero fold is an exact invariant subspace (a pure strain exerts
    no fold force, 1e-15) and the strain energy is exactly quadratic --
    V(10u) / V(u) = 100 to twelve digits. The shove wave of
    shear_response can never steepen, never break, and can carry no
    soliton. ALL of the medium's nonlinearity is the fold's.

WHAT THE FOLD'S NONLINEARITY IS. Its on-site potential is exactly
QUARTIC (R8): a uniform fold at clamped transverse spacing costs
4 delta^4 per cell with no quadratic and no cubic term -- again dL/da = 0,
so the breathe is free to first order and pays only at fourth. Between
planes the cubic terms do not vanish, and the fold-to-strain source
d3V/du_{n+1} d gamma_n^2 = STEP exactly.

THE ANSWER (S1-S4): YES -- an ENVELOPE soliton. A fold wavepacket at a
carrier wavenumber, with its amplitude and width in the right relation,
propagates at the group velocity with its envelope UNCHANGED over
hundreds of lattice steps, while the same packet in the harmonic model
disperses away. Two of them collide and come out with their amplitudes
and widths intact and a position shift -- the signature that separates a
soliton from a merely solitary wave, and the exact thing linear
superposition cannot do.

WHY AN ENVELOPE SOLITON AND NOT A KdV PULSE. At LONG wavelength the fold
and strain branches are not just degenerate but phase-matched: a fold
pulse's g^2 source travels at exactly the speed of the strain wave it
drives, so the transfer is secular and a broad fold pulse pumps itself
into a strain wave instead of holding together (measured: a 2-degree
Gaussian of width 16 loses half its height by t = 200 and goes bipolar).
At a finite carrier the second-harmonic resonance is detuned -- 2 omega(q)
- omega(2q) grows like q^3 -- and what is left is the hard on-site
quartic against the band's negative curvature: focusing, and a bright
envelope soliton. The medium selects the soliton it can hold.

UNITS as the rest of the lane: k = 1, plate mass 1 per triangle, distance
in primitive lattice steps, time as in every spectrum module. Fold
deviations are carried in RADIANS from a = -30 and reported in degrees.

SCOPE. Exact kinematics and exact soft-joint energy, no linearisation
anywhere; the approximation is the PLANE-SYMMETRIC RESTRICTION, which is
a uniaxial-strain plane wave in the bulk (transverse spacing clamped --
the correct condition for a plane wave in an infinite medium, and the
reason the breathe costs quartically here rather than nothing). Fold
excursions stay well inside the physical |a| <= 60 window. The collision
is measured to be NEARLY elastic and the radiation is quantified; this
lattice is not integrable and the row does not claim it is.
"""
from __future__ import annotations

import base64
import json
import pathlib
import sys

import numpy as np

from analysis.model import assembly as RC
from analysis.model import dispersion as OC
from analysis.model.assembly import VMASS, body, hat

A_REF = -30.0
K_JOINT = 1.0
J7, M7, BONDS = OC.periodic_cell()
BND = {e: prs for (e, prs) in BONDS}
_two, _ = RC.honeycomb_single([(0, 0, 0), (2, 0, 0)], gc=A_REF)
#: centre-to-centre distance between axis neighbours, the model's own
STEP = float(_two.ctr0[1][0])
XP = list(BND[(1, 0, 0)])
NV = 12
EX = np.array([1.0, 0.0, 0.0])
EY = np.array([0.0, 1.0, 0.0])
EZ = np.array([0.0, 0.0, 1.0])

#: the soliton this lane exhibits: carrier, envelope width, amplitude
Q0 = np.pi / 2
SOL_W = 8.0
SOL_AMP_DEG = 6.0
#: the family's invariant, amplitude x width (S2)
CONST_AW = SOL_AMP_DEG * SOL_W
DT = 0.05


# --------------------------------------------------------------------------
# the cell's closed form, vectorised over planes
# --------------------------------------------------------------------------

def basis():
    """x_v(a) = P cos(a-60) + Q sin(a-60) + S + T cos(a) -- `body`'s own closed
    form, whose coefficients are recovered by solving at four angles rather
    than re-deriving FACES here. R1 gates the result against body() itself,
    value and both derivatives."""
    angs = np.array([-70.0, -20.0, 15.0, 55.0])
    A = np.zeros((4, 4))
    for i, ad in enumerate(angs):
        r = np.radians(ad)
        A[i] = [np.cos(r - np.pi / 3), np.sin(r - np.pi / 3), 1.0, np.cos(r)]
    Y = np.array([body(ad) for ad in angs])
    return np.linalg.solve(A, Y.reshape(4, -1)).reshape(4, NV, 3)


COEF = basis()


def _trig(g):
    a = np.radians(A_REF) + np.asarray(g, float)
    return (np.cos(a - np.pi / 3), np.sin(a - np.pi / 3), np.cos(a), np.sin(a))


def verts(g, nder=1):
    """Vertex positions (and d/dg, d2/dg2) for fold deviations g in radians
    from A_REF. Shapes (n, 12, 3)."""
    c1, s1, c2, s2 = _trig(g)
    P, Q, S, T = COEF

    def e(v, w):
        return v[None, :, :] * np.asarray(w)[:, None, None]

    out = [e(P, c1) + e(Q, s1) + S[None] + e(T, c2)]
    if nder >= 1:
        out.append(-e(P, s1) + e(Q, c1) - e(T, s2))
    if nder >= 2:
        out.append(-e(P, c1) - e(Q, s1) - e(T, c2))
    return out


def _self_coeffs():
    """The y and z bonds tie a cell to its OWN transverse image, so their
    residual is a difference of two of its own vertices minus one lattice
    step -- a function of that plane's fold alone."""
    P, Q, S, T = COEF
    out = []
    for prs, e in ((BND[(0, 1, 0)], EY), (BND[(0, 0, 1)], EZ)):
        for (a, b) in prs:
            out.append((P[a] - P[b], Q[a] - Q[b], S[a] - S[b], T[a] - T[b],
                        STEP * e))
    return out


SELF = _self_coeffs()


def _mass_coeffs():
    """M_gamma(g) = sum_v m_v |dx_v/dg|^2 is a trig polynomial in g; six
    scalars, computed once. R2 gates it against the direct vertex sum."""
    P, Q, S, T = COEF

    def dot(X, Y):
        return float(np.sum(VMASS * np.einsum('vi,vi->v', X, Y)))

    return (dot(P, P), dot(Q, Q), dot(T, T),
            -2 * dot(P, Q), 2 * dot(P, T), -2 * dot(Q, T))


MC = _mass_coeffs()
MU = float(VMASS.sum())          # 8: the plane's translational mass


def mass_fold(g, deriv=False):
    """M_gamma(g), and dM/dg on request, from the closed form."""
    c1, s1, c2, s2 = _trig(g)
    A, B, C, D, E, F = MC
    m = (A * s1 * s1 + B * c1 * c1 + C * s2 * s2
         + D * s1 * c1 + E * s1 * s2 + F * c1 * s2)
    if not deriv:
        return m
    dm = (2 * A * s1 * c1 - 2 * B * c1 * s1 + 2 * C * s2 * c2
          + D * (c1 * c1 - s1 * s1) + E * (c1 * s2 + s1 * c2)
          + F * (c1 * c2 - s1 * s2))
    return m, dm


def mass_fold_direct(g):
    """The same sum taken straight over the twelve vertices."""
    _, V1, _ = verts(g, 2)
    return np.einsum('v,nvi,nvi->n', VMASS, V1, V1)


# --------------------------------------------------------------------------
# the two-field chain: exact energy, exact gradient
# --------------------------------------------------------------------------

def energy_grad(u, g, k=K_JOINT):
    """(V, dV/du, dV/dg) for a periodic ring of planes, from the TRUE joint
    separations. `np.roll`, never `np.add.at`: the partner index is a shift,
    and add.at on a permutation costs more than the whole force evaluation."""
    c1, s1, c2, s2 = _trig(g)
    V = 0.0
    dVg = np.zeros(len(u))
    for (dP, dQ, dS, dT, T) in SELF:
        r = (np.outer(c1, dP) + np.outer(s1, dQ) + dS[None, :]
             + np.outer(c2, dT) - T[None, :])
        dr = -np.outer(s1, dP) + np.outer(c1, dQ) - np.outer(s2, dT)
        V += 0.5 * k * float(np.sum(r * r))
        dVg += k * np.sum(r * dr, axis=1)
    P, Q, S, T4 = COEF
    dVu = np.zeros(len(u))
    # the ring closes on itself: plane n-1's partner is plane 0 one lattice
    # length away, and (n-1)STEP - nSTEP = -STEP is the interior formula again
    dx0 = u - np.roll(u, -1) - STEP
    for (va, vb) in XP:
        pa = (np.outer(c1, P[va]) + np.outer(s1, Q[va]) + S[va][None, :]
              + np.outer(c2, T4[va]))
        da = -np.outer(s1, P[va]) + np.outer(c1, Q[va]) - np.outer(s2, T4[va])
        pb = (np.outer(c1, P[vb]) + np.outer(s1, Q[vb]) + S[vb][None, :]
              + np.outer(c2, T4[vb]))
        db = -np.outer(s1, P[vb]) + np.outer(c1, Q[vb]) - np.outer(s2, T4[vb])
        r = pa - np.roll(pb, -1, axis=0)
        r[:, 0] += dx0
        V += 0.5 * k * float(np.sum(r * r))
        f = k * r[:, 0]
        dVu += f - np.roll(f, 1)
        dVg += k * np.sum(r * da, axis=1)
        dVg -= np.roll(k * np.sum(r * np.roll(db, -1, axis=0), axis=1), 1)
    return V, dVu, dVg


def accel(u, g, ud, gd, k=K_JOINT):
    """The exact equations of motion. The fold's inertia is configuration
    dependent -- M_gamma(g) -- so its equation carries the -(1/2)M' gdot^2
    term, which is `swing`'s own 1-DOF equation written per plane."""
    _V, dVu, dVg = energy_grad(u, g, k)
    m, dm = mass_fold(g, deriv=True)
    return -dVu / MU, (-0.5 * dm * gd * gd - dVg) / m


def energy(u, g, ud, gd, k=K_JOINT):
    V, _, _ = energy_grad(u, g, k)
    return V + 0.5 * float(np.sum(MU * ud * ud + mass_fold(g) * gd * gd))


def rk4(u, g, ud, gd, dt=DT, steps=1000, k=K_JOINT, sample=None, linear=False):
    """Fixed-step RK4 on the exact system (or the harmonic control)."""
    if linear:
        return rk4_linear(u, g, ud, gd, dt, steps, k, sample)
    n = len(u)
    y = np.concatenate([u, g, ud, gd])
    out = []

    def f(y):
        au, ag = accel(y[:n], y[n:2 * n], y[2 * n:3 * n], y[3 * n:], k)
        return np.concatenate([y[2 * n:3 * n], y[3 * n:], au, ag])

    for s in range(steps + 1):
        if sample and s % sample == 0:
            out.append((s * dt, y[:n].copy(), y[n:2 * n].copy(),
                        energy(y[:n], y[n:2 * n], y[2 * n:3 * n], y[3 * n:], k)))
        if s == steps:
            break
        k1 = f(y)
        k2 = f(y + 0.5 * dt * k1)
        k3 = f(y + 0.5 * dt * k2)
        k4 = f(y + dt * k3)
        y = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return out, (y[:n], y[n:2 * n], y[2 * n:3 * n], y[3 * n:])


def linear_blocks(k=K_JOINT):
    """(Hnn, Hn,n+1, Hn,n-1): the 2x2 Hessian blocks at rest, EXACTLY.

    At the reference every joint separation is zero, so the second-derivative
    -of-the-residual term drops and the Hessian is just k (dr/dq)^T (dr/dq)
    with dr/dq read from the same closed form the force uses. Computed rather
    than differenced because R5 compares it to dispersion.py's published bands
    and R6 uses it as the harmonic control -- a differenced Hessian carries
    ~1e-10 of its own error, which is larger than the nonlinearity R6 is
    trying to see at small amplitude.
    """
    c1, s1, c2, s2 = (float(x) for x in _trig(0.0))
    Hnn = np.zeros((2, 2))
    Hnp = np.zeros((2, 2))
    # self bonds: fold only, on-site
    for (dP, dQ, dS, dT, _T) in SELF:
        d = -dP * s1 + dQ * c1 - dT * s2
        Hnn[1, 1] += k * float(d @ d)
    # x bonds: (u_n, g_n) against (u_{n+1}, g_{n+1})
    P, Q, S, T4 = COEF
    for (va, vb) in XP:
        da = -P[va] * s1 + Q[va] * c1 - T4[va] * s2
        db = -P[vb] * s1 + Q[vb] * c1 - T4[vb] * s2
        Ja = np.stack([EX, da], axis=1)          # (3, 2) for plane n
        Jb = np.stack([-EX, -db], axis=1)        # (3, 2) for plane n+1
        Hnn += k * (Ja.T @ Ja) + k * (Jb.T @ Jb)   # bond ahead and bond behind
        Hnp += k * (Ja.T @ Jb)
    return Hnn, Hnp, Hnp.T


def rk4_linear(u, g, ud, gd, dt, steps, k=K_JOINT, sample=None):
    """The harmonic model on the same initial condition: the control every
    soliton row is measured against."""
    Hnn, Hnp, Hnm = linear_blocks(k)
    M = np.array([MU, float(mass_fold(np.zeros(1))[0])])
    n = len(u)

    def f(y):
        s = np.stack([y[:n], y[n:2 * n]], axis=1)
        # np.dot, not @ -- the Accelerate BLAS status-flag issue the model's
        # own _solve documents
        F = -(np.dot(s, Hnn.T) + np.dot(np.roll(s, -1, axis=0), Hnp.T)
              + np.dot(np.roll(s, 1, axis=0), Hnm.T))
        a = F / M
        return np.concatenate([y[2 * n:3 * n], y[3 * n:], a[:, 0], a[:, 1]])

    y = np.concatenate([u, g, ud, gd])
    out = []
    for s in range(steps + 1):
        if sample and s % sample == 0:
            out.append((s * dt, y[:n].copy(), y[n:2 * n].copy(), 0.0))
        if s == steps:
            break
        k1 = f(y)
        k2 = f(y + 0.5 * dt * k1)
        k3 = f(y + 0.5 * dt * k2)
        k4 = f(y + dt * k3)
        y = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return out, (y[:n], y[n:2 * n], y[2 * n:3 * n], y[3 * n:])


# --------------------------------------------------------------------------
# the full seven-coordinate force, for the row that licenses the reduction
# --------------------------------------------------------------------------

def _rot(w):
    th = float(np.linalg.norm(w))
    if th < 1e-14:
        return np.eye(3)
    u = w / th
    Kx = hat(u)
    return np.eye(3) + np.sin(th) * Kx + (1.0 - np.cos(th)) * (Kx @ Kx)


def full_gradient(C, W, G, k=K_JOINT):
    """dV/d(centre, rotation vector, fold) per plane, with NOTHING frozen:
    three centre components, a finite rotation, and the fold. R4 evaluates
    this on states inside the two-field subspace and measures the components
    that the reduction throws away."""
    n = len(G)
    out = np.zeros((n, 7))
    Rs = [_rot(W[i]) for i in range(n)]
    Bs = [body(A_REF + np.degrees(G[i]), 1) for i in range(n)]
    Pw = [(Rs[i] @ Bs[i][0].T).T for i in range(n)]
    Dw = [(Rs[i] @ Bs[i][1].T).T for i in range(n)]

    def add(i, r, va, sign):
        blk = np.zeros((3, 7))
        blk[:, 0:3] = sign * np.eye(3)
        blk[:, 3:6] = -sign * hat(Pw[i][va])
        blk[:, 6] = sign * Dw[i][va]
        out[i] += k * (blk.T @ r)

    for i in range(n):
        for prs, e in ((BND[(0, 1, 0)], EY), (BND[(0, 0, 1)], EZ)):
            for (a, b) in prs:
                r = Pw[i][a] - Pw[i][b] - STEP * e
                blk = np.zeros((3, 7))
                blk[:, 3:6] = -hat(Pw[i][a]) + hat(Pw[i][b])
                blk[:, 6] = Dw[i][a] - Dw[i][b]
                out[i] += k * (blk.T @ r)
        j = (i + 1) % n
        shift = (n * STEP) * EX if j == 0 else np.zeros(3)
        for (a, b) in XP:
            r = (C[i] + Pw[i][a]) - (C[j] + Pw[j][b] + shift)
            add(i, r, a, +1.0)
            add(j, r, b, -1.0)
    return out


# --------------------------------------------------------------------------
# wavepackets and how they are read
# --------------------------------------------------------------------------

def packet(N, n0, w, amp_deg, q0=Q0, kind="sech", sign=1.0, right=True):
    """A fold wavepacket: envelope `kind` of width w at carrier q0, with the
    velocity that makes it travel at the group velocity.

    `right=False` reverses it. Note that NEGATING q0 does NOT: the carrier is
    real, so a real packet holds both +q0 and -q0 and the direction lives
    entirely in the velocity field's relative phase. Flipping the sign of q0
    reproduces the same right-mover exactly -- measured, after it silently
    sent both halves of a head-on collision the same way."""
    n = np.arange(N, dtype=float)
    d = (n - n0 + N / 2.0) % N - N / 2.0
    A = sign * np.radians(amp_deg)
    if kind == "sech":
        psi = A / np.cosh(d / w)
        dpsi = -psi * np.tanh(d / w) / w
    else:
        psi = A * np.exp(-(d / w) ** 2)
        dpsi = -2.0 * d / w ** 2 * psi
    w0, vg = np.sin(q0 / 2.0), 0.5 * np.cos(q0 / 2.0)
    s = 1.0 if right else -1.0
    return (psi * np.cos(q0 * n),
            s * (-vg * dpsi * np.cos(q0 * n) + w0 * psi * np.sin(q0 * n)))


def envelope(g):
    """|analytic signal| -- the envelope, with the carrier divided out, so a
    breathing carrier cannot be mistaken for a decaying pulse."""
    n = len(g)
    F = np.fft.fft(g)
    H = np.zeros(n)
    H[0] = 1.0
    H[1:(n + 1) // 2] = 2.0
    if n % 2 == 0:
        H[n // 2] = 1.0
    return np.abs(np.fft.ifft(F * H))


def peak_width(e, frac=0.5):
    """(peak, full width at `frac` of peak, index) of the tallest structure,
    read periodically so a packet at the seam is not split."""
    n = len(e)
    i = int(np.argmax(e))
    idx = (np.arange(n) - i + n // 2) % n - n // 2
    m = e >= frac * e[i]
    return float(e[i]), float(np.abs(idx[m]).max() * 2 + 1), i


def centroid(e, i=None, window=None):
    """Envelope centroid about its peak, periodic.

    `window` restricts the average to |n - peak| < window. Pass it whenever
    more than one packet is present: taken over the whole ring, the centroid
    of a two-packet field is their MIDPOINT, which silently reported both
    halves of a collision sitting on top of each other long before they met.
    """
    n = len(e)
    if i is None:
        i = int(np.argmax(e))
    d = (np.arange(n) - i + n // 2) % n - n // 2
    p = e ** 2
    if window is not None:
        p = np.where(np.abs(d) < window, p, 0.0)
    return float((i + (p * d).sum() / p.sum()) % n)


# --------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------

def gate():
    checks, out = [], {}
    A = checks.append
    rng = np.random.default_rng(0)

    # ---- R1: the vectorised kinematics ARE the model's own -----------------
    worst = 0.0
    for ad in (-58.0, -30.0, -12.3, 7.7):
        v = verts(np.array([np.radians(ad - A_REF)]), 2)
        ref = body(ad, 2)
        for j in range(3):
            worst = max(worst, float(np.abs(v[j][0] - ref[j]).max()))
    out["r1"] = worst
    A(("R1 the vectorised closed form IS body(): positions and both "
       "derivatives, at four fold angles", worst < 1e-14,
       f"worst {worst:.2e}"))

    # ---- R2: the fold's inertia --------------------------------------------
    gg = np.linspace(-0.4, 0.4, 25)
    m, dm = mass_fold(gg, deriv=True)
    md = float(np.abs(m - mass_fold_direct(gg)).max())
    h = 1e-6
    dnum = float(np.abs(dm - (mass_fold(gg + h) - mass_fold(gg - h)) / (2 * h)).max())
    m0 = float(mass_fold(np.zeros(1))[0])
    out["r2"] = (md, dnum, m0)
    A(("R2 the fold's configuration-dependent inertia: closed form == the "
       "direct vertex sum, dM/dg == central differences, and M(0) is the "
       "kick lane's own measured m_gamma = 8",
       md < 1e-12 and dnum < 1e-7 and abs(m0 - 8.0) < 1e-12,
       f"closed form vs direct {md:.2e}, dM/dg {dnum:.2e}, M(0) = {m0:.12f}, "
       f"M_u = {MU:.6f}"))

    # ---- R3: the reference configuration is exact --------------------------
    V0, du0, dg0 = energy_grad(np.zeros(9), np.zeros(9))
    u1 = 0.03 * rng.normal(size=9)
    g1 = 0.02 * rng.normal(size=9)
    _V, du, dg = energy_grad(u1, g1)
    hh = 1e-7
    nu, ng = np.zeros(9), np.zeros(9)
    for i in range(9):
        e = np.zeros(9)
        e[i] = hh
        nu[i] = (energy_grad(u1 + e, g1)[0] - energy_grad(u1 - e, g1)[0]) / (2 * hh)
        ng[i] = (energy_grad(u1, g1 + e)[0] - energy_grad(u1, g1 - e)[0]) / (2 * hh)
    gerr = max(float(np.abs(du - nu).max()), float(np.abs(dg - ng).max()))
    out["r3"] = (V0, gerr)
    A(("R3 the reference is exact and the gradient is the energy's: V and "
       "both gradients vanish at a = -30, and away from it the analytic "
       "gradient matches central differences",
       V0 < 1e-20 and abs(du0).max() < 1e-14 and gerr < 1e-8,
       f"V(ref) {V0:.2e}, |grad(ref)| {max(abs(du0).max(), abs(dg0).max()):.2e}; "
       f"gradient vs differences {gerr:.2e}"))

    # ---- R4: THE REDUCTION IS EXACT ---------------------------------------
    n4 = 11
    u4 = 0.02 * rng.normal(size=n4)
    g4 = np.radians(4.0) * rng.normal(size=n4)
    C4 = np.zeros((n4, 3))
    C4[:, 0] = np.arange(n4) * STEP + u4
    F = full_gradient(C4, np.zeros((n4, 3)), g4)
    tor = float(np.abs(F[:, 3:6]).max())
    tra = float(np.abs(F[:, 1:3]).max())
    scale = float(np.abs(F[:, [0, 6]]).max())
    out["r4"] = (tor, tra, scale)
    A(("R4 THE PLANE-SYMMETRIC MEDIUM IS EXACTLY A TWO-FIELD CHAIN: from a "
       "state with cells unrotated and displaced only along x, the full "
       "SEVEN-coordinate force has no torque and no transverse component, so "
       "the coordinates this module drops can never be excited. Two-sided: "
       "the kept components are of order one",
       tor < 1e-12 and tra < 1e-12 and scale > 1e-3,
       f"torque {tor:.2e}, transverse force {tra:.2e}, against kept "
       f"components of size {scale:.2e}"))

    # ---- R5: the harmonic limit IS dispersion.py's own bands ---------------
    Hnn, Hnp, Hnm = linear_blocks()
    Mb = np.array([MU, m0])
    ident = (abs(Hnn[0, 0] - 4.0) < 1e-6 and abs(Hnn[1, 1] - 4.0) < 1e-6
             and abs(Hnp[0, 0] + 2.0) < 1e-6 and abs(Hnp[1, 1] + 2.0) < 1e-6
             and abs(Hnn[0, 1]) < 1e-8 and abs(Hnp[0, 1]) < 1e-8
             and abs(Hnp[1, 0]) < 1e-8)
    worst, wsin, sp = 0.0, 0.0, []
    for q in (0.02, 0.3, 0.9, 1.7, 2.6, 3.0):
        Hq = Hnn + Hnp * np.exp(1j * q) + Hnm * np.exp(-1j * q)
        Hq = (Hq + Hq.conj().T) / 2.0
        w = np.sqrt(np.clip(np.linalg.eigvalsh(np.diag(1.0 / Mb) @ Hq).real,
                            0.0, None))
        ref = OC.bands(np.array([q, 0.0, 0.0]), (J7, M7, BONDS))
        for x in w:
            worst = max(worst, float(np.min(np.abs(ref - x))))
        wsin = max(wsin, float(np.max(np.abs(w - np.sin(q / 2.0)))))
        if q == 0.02:
            sp = sorted(w / q)
    out["r5"] = (worst, wsin, sp)
    A(("R5 THE TWO FIELDS ARE THE SAME CHAIN, and it is dispersion.py's: "
       "onsite 4, coupling -2, no harmonic cross term, so both branches are "
       "omega = sin(q/2) exactly -- degenerate, sound speed 1/2 -- and both "
       "are members of the medium's own bands along [100]",
       ident and worst < 1e-9 and wsin < 1e-9,
       f"membership in dispersion.bands([q,0,0]) {worst:.2e}; both branches "
       f"vs sin(q/2) {wsin:.2e}; speeds at q = 0.02 "
       f"{sp[0]:.6f}, {sp[1]:.6f}"))

    # ---- R6: the integrator ------------------------------------------------
    n6 = 200
    errs = []
    for amp in (1e-4, 1e-2):
        g0, gd0 = packet(n6, 60.0, 8.0, amp)
        _o, fn = rk4(np.zeros(n6), g0, np.zeros(n6), gd0, 0.05, 800)
        _o, fl = rk4(np.zeros(n6), g0, np.zeros(n6), gd0, 0.05, 800, linear=True)
        # per AMPLITUDE SQUARED: the leading nonlinearity is quadratic, so
        # this is the quantity that must be the same at both amplitudes
        errs.append(float(np.abs(fn[1] - fl[1]).max()) / np.radians(amp) ** 2)
    ratio = errs[1] / errs[0] if errs[0] > 0 else 0.0
    g0, gd0 = packet(n6, 60.0, 8.0, 6.0)
    rec, _ = rk4(np.zeros(n6), g0, np.zeros(n6), gd0, 0.05, 4000, sample=1000)
    drift = max(abs(r[3] / rec[0][3] - 1.0) for r in rec)
    fine = [rk4(np.zeros(n6), g0, np.zeros(n6), gd0, d, int(round(40.0 / d)))[1][1]
            for d in (0.08, 0.04, 0.02)]
    e1 = float(np.abs(fine[0] - fine[2]).max())
    e2 = float(np.abs(fine[1] - fine[2]).max())
    order = e1 / e2 if e2 > 0 else 0.0
    out["r6"] = (errs, ratio, drift, order)
    A(("R6 the nonlinear integrator converges to the harmonic one QUADRATICALLY "
       "in amplitude (the leading nonlinearity is quadratic), conserves energy, "
       "and RK4 converges at dt^4",
       0.9 < ratio < 1.1 and drift < 1e-5 and 10.0 < order < 22.0,
       f"|g_nonlinear - g_harmonic| / amplitude^2 = {errs[0]:.4f} at 1e-4 deg "
       f"and {errs[1]:.4f} at 1e-2 deg -- the same to {abs(ratio - 1) * 100:.1f}% "
       f"over a 100x amplitude range, so the deviation is exactly quadratic; "
       f"energy drift {drift:.1e}; dt refinement ratio {order:.1f} (16 = "
       f"fourth order)"))

    # ---- R7: a plane longitudinal wave is EXACTLY harmonic -----------------
    uu = rng.normal(size=13)
    vals, gmax = [], 0.0
    for lam in (0.5, 1.0, 3.0, 10.0):
        V, _d, dgg = energy_grad(lam * uu, np.zeros(13))
        vals.append(V / lam ** 2)
        gmax = max(gmax, float(np.abs(dgg).max()) / lam)
    spread = max(abs(v / vals[0] - 1.0) for v in vals)
    out["r7"] = (spread, gmax)
    A(("R7 A PLANE LONGITUDINAL WAVE IS EXACTLY HARMONIC, AT EVERY AMPLITUDE: "
       "zero fold is an exact invariant subspace (a pure strain exerts no fold "
       "force) and the strain energy is exactly quadratic out to ten times "
       "amplitude -- shear_response's shove wave can never steepen, so every "
       "nonlinearity in this medium is the fold's",
       spread < 1e-12 and gmax < 1e-13,
       f"V/lambda^2 constant to {spread:.1e} over lambda = 0.5..10; fold force "
       f"under a pure strain {gmax:.1e}"))

    # ---- R8: the fold's on-site potential is exactly quartic ---------------
    coef, quad = [], []
    for dd in (0.25, 0.5, 1.0, 2.0):
        d = np.radians(dd)
        V, _, _ = energy_grad(np.zeros(9), np.full(9, d))
        coef.append((V / 9) / d ** 4)
        quad.append((V / 9) / d ** 2)
    lift = quad[-1] / quad[0]
    out["r8"] = (coef, lift)
    A(("R8 THE FOLD'S ON-SITE POTENTIAL IS EXACTLY QUARTIC -- 4 delta^4 per "
       "cell, no quadratic and no cubic term. That is `lattice_constant` being "
       "STATIONARY at a = -30 (dL/da = 0): the breathe is free to first order "
       "and pays only at fourth, which is the same fact that makes the fold "
       "sector decouple linearly (kick_response R5). Two-sided: a quadratic "
       "term would show as a CONSTANT V/delta^2, and instead it grows as "
       "delta^2 -- 64x across a 8x range of delta",
       all(abs(c - 4.0) < 5e-3 for c in coef) and 55.0 < lift < 70.0,
       f"V/delta^4 = {', '.join('%.4f' % c for c in coef)} over delta = "
       f"0.25..2 deg; V/delta^2 GROWS {lift:.0f}x across the same range"))


    # ---- S1: THE SOLITON --------------------------------------------------
    n1, st1 = 600, 10000
    g0, gd0 = packet(n1, 150.0, SOL_W, SOL_AMP_DEG)
    res = {}
    for lin in (False, True):
        rec, _ = rk4(np.zeros(n1), g0, np.zeros(n1), gd0, DT, st1,
                     sample=st1, linear=lin)
        e0, e1 = envelope(rec[0][2]), envelope(rec[-1][2])
        p0, w0_, i0 = peak_width(e0)
        p1, w1_, i1 = peak_width(e1)
        res[lin] = (p1 / p0, w1_ / w0_,
                    (centroid(e1, i1) - centroid(e0, i0)) / (st1 * DT))
    (kn, wn, vn), (kl, wl, vl) = res[False], res[True]
    out["s1"] = (kn, wn, vn, kl, wl, vl)
    A((f"S1 THE SOLITON: a {SOL_AMP_DEG:.0f}-degree fold wavepacket of envelope "
       f"width {SOL_W:.0f} travels {vn * st1 * DT:.0f} lattice steps with its "
       "envelope INTACT, while the SAME initial condition in the harmonic "
       "model disperses. Two-sided in the only way that matters: the control "
       "must fail, and it does",
       kn > 0.95 and wn < 1.25 and kl < 0.8 and wl > 1.5
       and abs(vn - 0.5 * np.cos(Q0 / 2)) < 0.01,
       f"nonlinear keeps {100 * kn:.0f}% of its height and {wn:.2f}x its width; "
       f"harmonic keeps {100 * kl:.0f}% and spreads to {wl:.2f}x; both travel at "
       f"the group velocity ({vn:.4f} and {vl:.4f} vs "
       f"{0.5 * np.cos(Q0 / 2):.4f})"))

    # ---- S2: the amplitude-width law --------------------------------------
    n2, st2 = 500, 6000
    fam, prod = [], []
    for amp, w in ((3.0, 16.0), (6.0, 8.0), (12.0, 4.0)):
        g0, gd0 = packet(n2, 120.0, w, amp)
        rec, _ = rk4(np.zeros(n2), g0, np.zeros(n2), gd0, DT, st2,
                     sample=st2 // 4)
        hs = [peak_width(envelope(r[2]))[0] for r in rec]
        # the extremes over the run, not the endpoint: a lattice soliton is
        # not exactly the continuum sech, so it BREATHES a little about it,
        # and an endpoint snapshot can land on either side of that breath
        fam.append((min(hs) / hs[0], max(hs) / hs[0]))
        prod.append(amp * w)
    out["s2"] = (fam, prod)
    A(("S2 THE SOLITON FAMILY IS AMPLITUDE x WIDTH = CONSTANT, which is what a "
       "focusing envelope equation predicts: three members an octave apart in "
       "amplitude, each with width scaled inversely, each holding its height",
       all(lo > 0.85 and hi < 1.25 for (lo, hi) in fam)
       and max(prod) - min(prod) < 1e-9,
       f"(A, w) = (3, 16), (6, 8), (12, 4), all with A*w = {prod[0]:.0f}; over "
       f"{st2 * DT:.0f} time units each stays within "
       f"{', '.join('%.0f-%.0f%%' % (100 * lo, 100 * hi) for (lo, hi) in fam)} "
       f"of its launch height -- breathing about the lattice soliton, not "
       f"decaying"))


    # ---- S3: THE COLLISION -------------------------------------------------
    n3, st3, win = 500, 13000, 60
    xa, xb = 130.0, 370.0
    samp = st3 // 90

    def _pair(right):
        return packet(n3, xa if right else xb, SOL_W, SOL_AMP_DEG, right=right)

    def _series(lin):
        """Per frame: how far each packet is from where it would have been had
        the other never existed. The packet BREATHES about the exact lattice
        soliton, so a single end-of-run reading lands wherever that breath
        happens to be -- measured, it ranges over about 0.2 of a cell. The
        shift is therefore the MEAN over the post-collision frames, and the
        spread is reported rather than hidden."""
        ga, va = _pair(True)
        gb, vb = _pair(False)
        z = np.zeros(n3)
        rc, _ = rk4(z, ga + gb, z, va + vb, DT, st3, sample=samp, linear=lin)
        rr, _ = rk4(z, ga, z, va, DT, st3, sample=samp, linear=lin)
        rl, _ = rk4(z, gb, z, vb, DT, st3, sample=samp, linear=lin)
        rows = []
        idx = np.arange(n3)
        for k in range(len(rc)):
            er, el = envelope(rr[k][2]), envelope(rl[k][2])
            ec = envelope(rc[k][2])
            cr = centroid(er, int(np.argmax(er)), win)
            cl = centroid(el, int(np.argmax(el)), win)
            if abs(cr - cl) < 110:
                rows.append((rc[k][0], None, None, None))
                continue
            got = []
            for c0 in (cr, cl):
                d = (idx - int(round(c0)) + n3 // 2) % n3 - n3 // 2
                i = int(np.argmax(np.where(np.abs(d) < 50, ec, -1.0)))
                got.append((centroid(ec, i, win) - c0, float(ec[i])))
            hr = float(er[int(round(cr)) % n3])
            rows.append((rc[k][0], got[0][0], got[1][0],
                         got[0][1] / hr if hr > 0 else 0.0))
        return rows

    rows = _series(False)
    tmax3 = rows[-1][0]
    post = [r for r in rows if r[1] is not None and r[0] > 0.72 * tmax3]
    pre = [r for r in rows if r[1] is not None and r[0] < 0.25 * tmax3]
    shift_r = float(np.mean([r[1] for r in post]))
    shift_l = float(np.mean([r[2] for r in post]))
    spread = float(np.std([r[1] for r in post]))
    before = float(np.mean([abs(r[1]) for r in pre]))
    keep = float(np.mean([r[3] for r in post]))
    lrows = _series(True)
    lpost = [r for r in lrows if r[1] is not None and r[0] > 0.72 * tmax3]
    lshift = float(np.mean([r[1] for r in lpost]))
    lkeep = float(np.mean([r[3] for r in lpost]))
    # the harmonic model's zero shift is an identity, not a coincidence:
    # superposition holds exactly there
    ga, va = _pair(True)
    gb, vb = _pair(False)
    z = np.zeros(n3)
    _o, fboth = rk4(z, ga + gb, z, va + vb, DT, 4000, linear=True)
    _o, fA = rk4(z, ga, z, va, DT, 4000, linear=True)
    _o, fB = rk4(z, gb, z, vb, DT, 4000, linear=True)
    sup = float(np.abs(fboth[1] - (fA[1] + fB[1])).max())
    out["s3"] = (shift_r, shift_l, spread, before, keep, lshift, sup)
    A(("S3 THE COLLISION -- what separates a soliton from a merely solitary "
       "wave. Two packets meet head-on, pass through, and come out at full "
       "height but displaced from where free propagation would have left "
       "them, each advanced along its own direction. The harmonic control "
       "cannot do this: there superposition is exact, so the two never "
       "interact and the shift is identically zero. WHAT THE SHIFT IS "
       "(envelope S1, 2026-09-02): these packets are launched COLD, with the "
       "strain at rest, and each radiates strain at c = 1/2 while its dent "
       "forms. Most of this row's shift is that launch transient: launched "
       "warm, with the dents "
       "in place, the collision's own shift is a tenth of a cell, and the "
       "envelope equation predicts it",
       keep > 0.93 and shift_r > 0.3 and shift_r * shift_l < 0
       and abs(abs(shift_l) - shift_r) < 0.25 * shift_r
       and before < 0.02 and abs(lshift) < 1e-3 and sup < 1e-12,
       f"height kept {100 * keep:.0f}% (harmonic {100 * lkeep:.0f}%); shift "
       f"{shift_r:+.3f} and {shift_l:+.3f} cells, mean over {len(post)} "
       f"post-collision frames, spread {spread:.3f} (the packet breathes); "
       f"before the meeting {before:.4f}; harmonic shift {lshift:+.5f} with "
       f"superposition residual {sup:.1e}"))


    # ---- S4/S5: how far the scaling goes, and what stops it ---------------
    def _hold(w, amp, q0=Q0, widths=12.0):
        """Run one family member for a fixed number of ITS OWN envelope widths,
        so every scale is judged on the same dimensionless clock."""
        vg = 0.5 * np.cos(q0 / 2.0)
        tmax = widths * w / vg
        nn = int(min(1400, max(200, 2.4 * vg * tmax + 14 * w)))
        st = int(round(tmax / DT))
        g0, gd0 = packet(nn, nn * 0.12, w, amp, q0=q0)
        z = np.zeros(nn)
        rec, _ = rk4(z, g0, z, gd0, DT, st, sample=max(1, st // 6))
        h = [float(np.degrees(peak_width(envelope(r[2]))[0])) for r in rec]
        return min(h) / h[0], max(h) / h[0]

    fam = {}
    for w in (2.0, 8.0, 16.0):
        fam[w] = _hold(w, CONST_AW / w)
    out["s4"] = fam
    A(("S4 THE SCALING WINDOW: the family is one wave at many sizes -- amplitude "
       "x width constant -- and it holds over that range until the ENVELOPE gets "
       "down to a few cells, where the lattice stops being ignorable. Two-sided: "
       "the wide members must hold AND the narrow one must visibly fail to. The "
       "family forces the narrow member to be TALL -- at w = 2 it is 24 degrees "
       "-- and that is still inside the physical +-60 window, so what fails "
       "there is the lattice, not the range",
       fam[8.0][0] > 0.95 and fam[16.0][0] > 0.95
       and fam[8.0][1] < 1.12 and fam[16.0][1] < 1.12
       and fam[2.0][0] < 0.9,
       "; ".join(f"w = {w:g} (A = {CONST_AW / w:.1f} deg): {lo:.2f}-{hi:.2f}"
                 for w, (lo, hi) in sorted(fam.items()))
       + " of launch height over 12 own-widths of travel"))

    car = {}
    for nm, q in (("pi/2", np.pi / 2), ("pi/8", np.pi / 8)):
        car[nm] = _hold(SOL_W, SOL_AMP_DEG, q0=q)
    mism = {nm: abs(2 * np.sin(q / 2) - np.sin(q))
            for nm, q in (("pi/2", np.pi / 2), ("pi/8", np.pi / 8))}
    out["s5"] = (car, mism)
    A(("S5 THE OTHER WALL IS ON THE CARRIER, NOT THE ENVELOPE. Widening the "
       "envelope detunes nothing, but lengthening the CARRIER drives the fold "
       "and strain branches into phase-match -- 2 omega(q) - omega(2q) falls "
       "like q^3 -- and the packet pumps itself into the strain wave it drives. "
       "Same envelope, same amplitude, only the carrier changed",
       car["pi/2"][0] > 0.95 and car["pi/8"][0] < 0.9,
       f"carrier pi/2 (4 cells, mismatch {mism['pi/2']:.2e}) keeps "
       f"{car['pi/2'][0]:.2f}; carrier pi/8 (16 cells, mismatch "
       f"{mism['pi/8']:.2e}) keeps {car['pi/8'][0]:.2f}"))

    return checks, out


# --------------------------------------------------------------------------
# frames export for the page
# --------------------------------------------------------------------------

def _quant(frames):
    """Per-run int8 quantisation with ONE scale for the whole run, so a decay
    is visible as a decay rather than being renormalised away frame by frame."""
    sc = float(np.max(np.abs(frames)))
    sc = sc if sc > 0 else 1.0
    return sc, np.clip(np.round(frames / sc * 127.0), -127, 127).astype(np.int8)


def _run(u0, g0, ud0, gd0, n, steps, sample, linear):
    rec, _ = rk4(u0, g0, ud0, gd0, DT, steps, sample=sample, linear=linear)
    return (np.array([r[0] for r in rec]),
            np.array([r[2] for r in rec]),
            np.array([r[1] for r in rec]))


def _pack(times, G, U):
    scg, qg = _quant(G)
    env = np.array([envelope(g) for g in G])
    sce, qe = _quant(env)
    return {
        "times": [round(float(t), 3) for t in times],
        "scale_deg": float(np.degrees(scg)),
        "fold_i8": base64.b64encode(qg.tobytes()).decode(),
        "env_scale_deg": float(np.degrees(sce)),
        "env_i8": base64.b64encode(qe.tobytes()).decode(),
        "peak_deg": [round(float(np.degrees(e.max())), 4) for e in env],
    }


def export(path="analysis/.pages/data/soliton.json"):
    """Both scenes, each with its harmonic control: one soliton crossing the
    lattice, and two colliding head-on."""
    data = {"step": STEP, "q0": Q0, "dt": DT, "amp_deg": SOL_AMP_DEG,
            "width": SOL_W, "vg": 0.5 * np.cos(Q0 / 2.0),
            "omega0": np.sin(Q0 / 2.0), "a_ref": A_REF}

    # ring sizes are chosen so the packet crosses nearly the whole lattice
    # WITHOUT wrapping, and so a cell is still a few pixels wide when the
    # page draws the medium at this length
    n, steps, sample = 360, 18000, 60
    g0, gd0 = packet(n, 25.0, SOL_W, SOL_AMP_DEG)
    z = np.zeros(n)
    for name, lin in (("solitonNL", False), ("solitonLin", True)):
        t, G, U = _run(z, g0, z, gd0, n, steps, sample, lin)
        data[name] = _pack(t, G, U)
    data["n_one"] = n

    n2, steps2, sample2 = 420, 12400, 41
    ga, va = packet(n2, 110.0, SOL_W, SOL_AMP_DEG, right=True)
    gb, vb = packet(n2, 310.0, SOL_W, SOL_AMP_DEG, right=False)
    z2 = np.zeros(n2)
    for name, lin in (("collideNL", False), ("collideLin", True)):
        t, G, U = _run(z2, ga + gb, z2, va + vb, n2, steps2, sample2, lin)
        data[name] = _pack(t, G, U)
        # THE SIGNATURE, frame by frame: where the right-moving packet is
        # against where it would have been had it never met the other one.
        # Undefined while the two overlap, and emitted as null there rather
        # than as a number the tracker cannot actually resolve.
        _t, GR, _u = _run(z2, ga, z2, va, n2, steps2, sample2, lin)
        _t, GL, _u = _run(z2, gb, z2, vb, n2, steps2, sample2, lin)
        off = []
        for f in range(len(t)):
            er, el = envelope(GR[f]), envelope(GL[f])
            cr = centroid(er, int(np.argmax(er)), 60)
            cl = centroid(el, int(np.argmax(el)), 60)
            sep = abs(cr - cl)
            if sep < 95:
                off.append(None)
                continue
            ec = envelope(G[f])
            lo = int(round(cr)) - 45
            w = np.zeros(n2, bool)
            w[[(lo + j) % n2 for j in range(90)]] = True
            i = int(np.argmax(np.where(w, ec, -1.0)))
            off.append(round(float(centroid(ec, i, 45) - cr), 4))
        data[name]["shift"] = off
    data["n_two"] = n2

    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data))
    return p, data


def export_ref(path="analysis/.pages/data/soliton_ref.json"):
    """A small FLOAT reference run for a page that integrates this model
    itself: the same two-field chain, in full precision, so an in-browser
    integrator can be measured against the gated one at load rather than
    trusted. Same contract as the kick lane's export-ref."""
    n, tmax, every = 120, 200.0, 20.0
    steps = int(round(tmax / DT))
    sample = int(round(every / DT))
    g0, gd0 = packet(n, 30.0, SOL_W, SOL_AMP_DEG)
    z = np.zeros(n)
    rec, _ = rk4(z, g0, z, gd0, DT, steps, sample=sample)
    data = {
        "n": n, "dt": DT, "amp_deg": SOL_AMP_DEG, "width": SOL_W,
        "q0": Q0, "x0": 30.0, "step": STEP, "a_ref": A_REF,
        "times": [round(float(r[0]), 4) for r in rec],
        "fold": [[float(f"{v:.7e}") for v in r[2]] for r in rec],
        "u": [[float(f"{v:.7e}") for v in r[1]] for r in rec],
        # the model's OWN constants, so a page that integrates this chain
        # evaluates the same energy rather than a re-derivation of it
        "coef": [[[float(x) for x in v] for v in blk] for blk in COEF],
        "self": [[[float(x) for x in c] for c in row] for row in SELF],
        "xp": [[int(a_), int(b_)] for (a_, b_) in XP],
        "mc": [float(x) for x in MC],
        "mu": float(MU),
    }
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data))
    return p, data


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "export-ref":
        p, d = export_ref(sys.argv[2] if len(sys.argv) > 2
                          else "analysis/.pages/data/soliton_ref.json")
        print(f"exported reference: {len(d['times'])} frames of {d['n']} "
              f"planes -> {p}")
        return 0
    if len(sys.argv) > 1 and sys.argv[1] == "export":
        path = sys.argv[2] if len(sys.argv) > 2 else "analysis/.pages/data/soliton.json"
        p, d = export(path)
        print(f"exported soliton frames -> {p} "
              f"({len(d['solitonNL']['times'])} + {len(d['collideNL']['times'])} frames)")
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
