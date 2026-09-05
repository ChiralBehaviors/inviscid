"""envelope -- RUNG 4: the envelope equation, and what it can and cannot see.

THE QUESTION (owner, 2026-09-02, the hierarchy page's open rung): the
soliton's carrier is four cells and its envelope is tens to hundreds. Can
the carrier be dropped and the envelope evolved on its own, and is the
result the medium's own or a fit?

THE METHOD is the project's: find the reduction, GATE it against the full
model, then run the reduced model where the full one cannot go. Here the
full model is soliton.py's exact two-field chain (itself gated against
the seven-coordinate medium, R4 there), and the reduction is the
nonlinear Schroedinger equation for the fold envelope A(X, T),

    gamma_n(t) = Re[ A(n, t) exp(i (q0 n - omega0 t)) ],
    i (A_T + v_g A_X) + P A_XX + Q |A|^2 A = 0,

with TWO coefficients. Neither is fitted.

  * P = omega''(q0) / 2 comes from the band omega = sin(q/2) that
    soliton R5 gates in closed form: P = -sin(q0/2) / 8.
  * Q is the nonlinear frequency shift, and it has TWO parts that must be
    measured separately, because a periodic ring can only see one of
    them.
      - Q_ring: a uniform plane wave on a ring at amplitude a runs at
        omega0 + Q_ring a^2 (R3). This holds the on-site quartic, the
        fold's configuration-dependent inertia, the inter-plane quartic
        terms and the strain driven at the SECOND harmonic -- everything
        except a mean strain, which a fixed-length ring cannot carry.
      - Q_dent: the fold-strain cubic is a GRADIENT coupling (R1),
        V3 = (STEP/2) sum gamma_n^2 (u_{n+1} - u_{n-1}), so a moving
        packet drives the strain chain with the gradient of its own
        |A|^2 and DRAGS A CONTRACTION along with it (R4),
            eps = -STEP <gamma^2> / (2 - M_u v_g^2),
        finite because v_g < c = 1/2 (the same mismatch that detunes the
        long-wave resonance). Sitting in its own dent the packet's
        on-site stiffness drops by 2 STEP eps, which SOFTENS the
        frequency shift: Q_dent = -STEP^2 / (16 omega0 (2 - M_u v_g^2)).
    Q = Q_ring + Q_dent. At q0 = pi/2 the ring gives 0.707 and the dent
    takes back 0.471, leaving 0.236 -- the dent is two thirds the size of
    everything else put together, and a derivation that forgot it would
    be off by a factor three in Q and 1.7 in the soliton's size.

THE GATE THAT DECIDES IT (R5). NLS has a one-parameter soliton family with
    amplitude x width = sqrt(2P/Q)
and NOTHING adjustable. With the coefficients above that is 49.6
degree-cells. The exact chain, asked directly -- at width 16, which
amplitude neither grows nor decays over a long run -- answers 48-49. The
two agree to a few percent, and the alternatives do not: the on-site
quartic alone predicts 23, the ring alone 28. The number the soliton lane
measured empirically (48, S2 there) is this prediction.

THE REDUCED MODEL (R6, R7). A split-step Fourier integrator for the
envelope on a grid of one point per two cells with a time step ten times
the chain's, checked against the chain's own envelope on the soliton
lane's scene: same height, same width, same position after 130 cells of
travel, in a hundredth of the time. The harmonic control (Q = 0) matches
the harmonic chain's dispersing packet the same way, which is P on its
own.

WHAT THE ENVELOPE MODEL CANNOT SEE (R8), stated as measurements rather
than caveats. It has no lattice, so at w = 2 it holds where the chain
loses a fifth of its height (the discreteness wall, soliton S4). And it
assumes the strain follows the packet instantly: the dent takes
2w/(c - v_g) to form, which diverges toward the resonance, so at carrier
pi/4 the chain packet runs on the bare ring coefficient for hundreds
of time units while the adiabatic Q says defocusing. That is the resonance
wall from rung 4's side: beyond it the strain is a slow field of its own,
not a coefficient. Both walls are real and both are where the
single-envelope description stops -- the same places the hierarchy page
said it would. What the chain then does at pi/4 -- focus to 1.3x and
beyond -- is NOT the half-built dent: it is an instability of the packet
to its own mean fold, which grows exponentially under it (rung 4b,
longwave.py L6, T2 [24273]); the dent story ends where R8b's row does.

THE COLLISION (S1), AND A CORRECTION TO THE SOLITON LANE. Two
counter-propagating envelopes couple only through cross-phase modulation,
Q_x |B|^2 A, with Q_x measured the same way Q is: two waves on a ring (R3
again), plus the dent the other packet carries. R3's Q_x,ring is the value
at RELATIVE CARRIER PHASE ZERO (both waves launched as cos(q n)); the
phase-averaged value is 1.47, and the difference is the umklapp four-wave
term at 4 q0 = 2 pi, which makes the pi/2 collision a phase-dependent
energy exchange (rung 4b, longwave.py L1). The envelope model's
collision displaces each soliton by 0.09 cells. The chain, with its
packets launched WARM -- their dents already in place -- gives 0.08. The
chain launched COLD, as soliton S3 launched it, gives 0.56-0.63, and that
is the number the soliton page published as the collision's signature. It
is mostly not the collision. A packet born without its dent radiates a
strain pulse at c = 1/2 while the dent forms; the pulse outruns the packet
(c > v_g) and reaches the other packet before the two meet. The strain's
initial state is the entire difference -- warm and cold launches differ in
nothing else -- and rung 4b (longwave.py L3) has the mechanism: the cold
launch detunes its own carrier through the strain transient, and the
detuned packet meets the other in the umklapp exchange regime. The
collision's own shift is a tenth of a cell, and rung 4 predicts it at
phase zero; with the umklapp term carried, the displacement is the one
number rung 4b does not yet reproduce (T2 [24240]). Every number in
soliton S3 stands as measured; what changes is what it measures.

THE CLOSED FORMS, DERIVED (R9). Lindstedt-Poincare on the exact chain's
plane wave, gamma = a cos th + a^2 g2 cos 2th + ..., u = a^2 d2 sin 2th,
th = q n - w t, w = s + a^2 w2, s = sin(q/2), with the Lagrangian's own
coefficients (each measured in R9 by central differences):
    M(g) = 8 - 4 STEP g + (16/3) g^2 + ...          (M1 = -4 STEP exactly)
    V3_fold = -(STEP/4) sum (dg)^2 (g_n + g_n+1)     (-3 STEP, STEP/2)
    V3_ug   = (STEP/2) sum g_n^2 (u_n+1 - u_n-1)
    V4      = sum 17/6 g^4 + 1/2 g_n^2 g_n+1^2 + 1/3 (g_n^3 g_n+1 + g_n g_n+1^3)
Order a^2: the DC fold source is M1 s^2/4 + STEP s^2 = 0 -- the breathe is
not driven -- and the second harmonics are g2 = STEP/16 (independent of q)
and d2 = STEP cos(q/2) cos q / (16 s^3), which VANISHES at q = pi/2: the
strain's second-harmonic source is sin 2q. Order a^3, the cos th
component that must vanish:
    16 s w2 = [mass] + [fold cubic x g2] + [strain 2q] + [quartic]
           = -(2/3)s^2 - (2/3)s^2(1 + 2c^2) + (4/3)c^2 cos^2 q / s^2
             + 12 - 8 s^2 + 4 s^4                       (c = cos(q/2))
           = 4 (1 + 4 s^2) / (3 s^2)                    (an identity)
so Q_ring = (1 + 4 s^2) / (12 s^3). The dent is R4's closed form in the
strain chain, Q_dent = -STEP^2 / (32 s^3) = -1/(6 s^3) (STEP^2 = 16/3),
and Q(q0) = (1 - 2 cos q0) / (12 s^3): it VANISHES at q0 = pi/3 and is
defocusing below it in the adiabatic limit -- but see R8b for why the
chain does not see that limit there.

UNITS as soliton.py: k = 1, distance in lattice steps, fold in radians
internally and degrees in every printed number. SCOPE: everything is at
fixed carrier q0 = pi/2 unless a row says otherwise; the plane-symmetric
restriction is inherited from the chain.
"""
from __future__ import annotations

import json
import pathlib
import re
import sys
import time

import numpy as np

from analysis.model import soliton as SL
from analysis.model.soliton import (MU, STEP, envelope, packet, peak_width,
                                    centroid, energy_grad, rk4)

Q0 = SL.Q0
SOL_W = SL.SOL_W
SOL_AMP_DEG = SL.SOL_AMP_DEG
DT = SL.DT
#: the strain chain's coupling: onsite 4, coupling -2 (soliton R5) is
#: 2 u_XX in the continuum
K_STRAIN = 2.0


# --------------------------------------------------------------------------
# the coefficients
# --------------------------------------------------------------------------

def omega(q):
    return np.sin(q / 2.0)


def vgroup(q0=Q0):
    return 0.5 * np.cos(q0 / 2.0)


def coef_p(q0=Q0):
    """P = omega''/2 from the band sin(q/2)."""
    return -np.sin(q0 / 2.0) / 8.0


def dent_depth(amp2, q0=Q0):
    """The strain a moving packet drags along: eps for <gamma^2> = amp2/2."""
    return -STEP * amp2 / (2.0 * (K_STRAIN - MU * vgroup(q0) ** 2))


def coef_q_dent(q0=Q0):
    """The frequency shift per a^2 from sitting in one's own dent."""
    return -STEP ** 2 / (16.0 * omega(q0) * (K_STRAIN - MU * vgroup(q0) ** 2))


def ring_shift(q0, N, amp, b_amp=0.0, T=200.0, dt=DT):
    """Frequency of the RIGHT-moving carrier wave of amplitude `amp` on a
    ring of N planes, minus omega0, with an optional left-mover of
    amplitude `b_amp` sharing the ring. The right mover is separated from
    the left one by combining the mode with its time derivative, and its
    frequency read from the slope of its phase."""
    w0 = omega(q0)
    nn = np.arange(N)
    g0 = (amp + b_amp) * np.cos(q0 * nn)
    gd0 = (amp - b_amp) * w0 * np.sin(q0 * nn)
    z = np.zeros(N)
    rec, _ = rk4(z, g0, z, gd0, dt, int(round(T / dt)), sample=1)
    ts = np.array([r[0] for r in rec])
    G = np.array([r[2] for r in rec])
    m = int(round(q0 * N / (2 * np.pi)))
    F = np.fft.fft(G, axis=1)[:, m] * 2.0 / N
    Fd = np.gradient(F, dt)
    w = w0
    for _ in range(3):
        R = (F + 1j * Fd / w) / 2.0
        ph = np.unwrap(np.angle(R[2:-2]))
        w = -np.polyfit(ts[2:-2], ph, 1)[0]
    return float(w - w0)


def ring_size(q0, at_least=60):
    """A ring the carrier fits on exactly: a multiple of its period."""
    per = int(round(2 * np.pi / q0))
    return per * int(np.ceil(at_least / per))


def ring_coefficient(q0=Q0, N=None, cross=False, T=200.0):
    """Q_ring (or Q_x,ring) per radian^2, extrapolated to zero amplitude
    from one- and two-degree measurements (the next correction is a^4).
    `T` must cover a few beats of the second-harmonic strain, whose period
    2 pi / (2 omega(q) - omega(2q)) lengthens as the carrier does."""
    N = ring_size(q0) if N is None else N
    a1, a2 = np.radians(1.0), np.radians(2.0)
    if not cross:
        s1 = ring_shift(q0, N, a1, T=T) / a1 ** 2
        s2 = ring_shift(q0, N, a2, T=T) / a2 ** 2
    else:
        base = ring_shift(q0, N, a1, T=T)
        s1 = (ring_shift(q0, N, a1, a1, T=T) - base) / a1 ** 2
        s2 = (ring_shift(q0, N, a1, a2, T=T) - base) / a2 ** 2
    return (4.0 * s1 - s2) / 3.0, s1, s2


def q_ring_closed(q0=Q0):
    """The ring shift, (1 + 4 s^2) / (12 s^3), s = sin(q0/2): derived in the
    docstring and R9, measured in R3."""
    s = np.sin(q0 / 2.0)
    return (1.0 + 4.0 * s * s) / (12.0 * s ** 3)


def q_total_closed(q0=Q0):
    """Ring plus dent: (4 s^2 - 1) / (12 s^3) = (1 - 2 cos q0) / (12 s^3).
    Vanishes at q0 = pi/3 and changes sign below it."""
    s = np.sin(q0 / 2.0)
    return (4.0 * s * s - 1.0) / (12.0 * s ** 3)


def q_ring_pieces(q0=Q0):
    """The four cos(theta) contributions to 16 s w2 in the Lindstedt
    expansion (docstring), from the Lagrangian coefficients as measured, so
    R9 can gate the derivation piece by piece. Returns (mass, fold cubic
    through g2, strain second harmonic, quartic), summing to 16 s Q_ring."""
    s2 = np.sin(q0 / 2.0) ** 2
    c2 = 1.0 - s2
    cq = 1.0 - 2.0 * s2
    M1, M2 = -4.0 * STEP, 32.0 / 3.0
    g2 = STEP / 16.0
    mass = s2 * (-1.5 * g2 * M1 - 0.25 * M2)
    fold3 = -2.0 * STEP * g2 * s2 * (1.0 + 2.0 * c2)
    strain2 = (4.0 / 3.0) * c2 * cq * cq / s2
    quartic = 12.0 - 8.0 * s2 + 4.0 * s2 * s2
    return mass, fold3, strain2, quartic


def nls_coefficients(q0=Q0, N=None):
    """(P, Q, v_g, omega0) in the NLS sign convention -- a uniform wave
    A = a exp(i Q a^2 T) runs at omega0 - Q a^2, so Q = -(ring + dent)."""
    qr, _s1, _s2 = ring_coefficient(q0, N)
    return coef_p(q0), -(qr + coef_q_dent(q0)), vgroup(q0), omega(q0)


def soliton_aw(P, Q):
    """The family invariant, amplitude x width, in radian-cells."""
    r = 2.0 * P / Q
    return float(np.sqrt(r)) if r > 0 else float("nan")


# --------------------------------------------------------------------------
# the reduced model: split-step Fourier NLS on a coarse grid
# --------------------------------------------------------------------------

def nls_run(A0, L, P, Q, vg, dt, steps, sample=None, B0=None, Qx=0.0):
    """Integrate i(A_T + vg A_X) + P A_XX + Q|A|^2 A (+ Qx |B|^2 A) = 0 on a
    periodic domain of length L cells with len(A0) points, Strang split:
    half a nonlinear phase, the full linear step exactly in Fourier space,
    half a nonlinear phase. `B0` adds the counter-propagating field, which
    sees -vg and the same P, Q, Qx.

    Returns (samples, A, B): samples is a list of (t, |A|, |B| or None)."""
    n = len(A0)
    dX = L / n
    k = 2.0 * np.pi * np.fft.fftfreq(n, d=dX)
    la = np.exp(-1j * (vg * k + P * k * k) * dt)
    lb = np.exp(-1j * (-vg * k + P * k * k) * dt)
    A = np.asarray(A0, complex).copy()
    B = None if B0 is None else np.asarray(B0, complex).copy()
    out = []

    def phase(h):
        nonlocal A, B
        if B is None:
            A = A * np.exp(1j * h * Q * np.abs(A) ** 2)
        else:
            a2, b2 = np.abs(A) ** 2, np.abs(B) ** 2
            A = A * np.exp(1j * h * (Q * a2 + Qx * b2))
            B = B * np.exp(1j * h * (Q * b2 + Qx * a2))

    for s in range(steps + 1):
        if sample and s % sample == 0:
            out.append((s * dt, np.abs(A), None if B is None else np.abs(B)))
        if s == steps:
            break
        phase(0.5 * dt)
        A = np.fft.ifft(np.fft.fft(A) * la)
        if B is not None:
            B = np.fft.ifft(np.fft.fft(B) * lb)
        phase(0.5 * dt)
    return out, A, B


def sech_envelope(n, L, x0, w, amp):
    """A sech envelope on the NLS grid, periodic, matching `packet`'s."""
    x = np.arange(n) * (L / n)
    d = (x - x0 + L / 2.0) % L - L / 2.0
    return amp / np.cosh(d / w)


def dent_fields(g, amp, right=True, q0=Q0):
    """The strain a packet WOULD have dragged had it been travelling forever:
    the co-moving dent of R4 under its envelope, with the velocity field of a
    profile moving at v_g. Launching with this is a WARM start; launching
    with u = 0 is a COLD one that radiates a strain pulse at c = 1/2 in both
    directions as the dent forms. On a ring the net contraction must vanish,
    so the dent is balanced by a uniform extension."""
    e = envelope(g)
    eps = dent_depth(amp * amp, q0) * (e / amp) ** 2
    eps = eps - eps.mean()
    u = np.concatenate([[0.0], np.cumsum(eps)[:-1]])
    u = u - u.mean()
    ud = -(1.0 if right else -1.0) * vgroup(q0) * (np.roll(u, -1) - np.roll(u, 1)) / 2.0
    return u, ud


# --------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------

def gate():
    checks, out = [], {}
    A = checks.append
    t_start = time.time()

    # ---- R1: the fold-strain cubic is a GRADIENT coupling ------------------
    N1, n, h = 9, 4, 1e-3

    def V(u, g):
        return energy_grad(u, g)[0]

    d3 = {}
    for m in (n - 1, n, n + 1):
        def f(du, dg):
            uu, gg = np.zeros(N1), np.zeros(N1)
            uu[m] += du
            gg[n] += dg
            return V(uu, gg)
        d3[m - n] = ((f(h, h) - 2 * f(h, 0) + f(h, -h))
                     - (f(-h, h) - 2 * f(-h, 0) + f(-h, -h))) / (2 * h ** 3)
    out["r1"] = d3
    A(("R1 THE FOLD-STRAIN CUBIC IS A GRADIENT COUPLING: d3V/du_{n+1} dg_n^2 "
       "= +STEP, d3V/du_{n-1} dg_n^2 = -STEP, and the same plane's own u "
       "does not appear -- so V3 = (STEP/2) sum g_n^2 (u_{n+1} - u_{n-1}) "
       "and a packet drives the strain with the GRADIENT of its |A|^2, "
       "which is why the strain it drags is a dent that moves with it",
       abs(d3[1] - STEP) < 1e-5 and abs(d3[-1] + STEP) < 1e-5
       and abs(d3[0]) < 1e-5,
       f"d3V/du_(n+1) = {d3[1]:.6f}, du_(n-1) = {d3[-1]:.6f}, du_n = "
       f"{d3[0]:.1e}; STEP = {STEP:.6f}"))

    # ---- R2: P from the chain's own band -----------------------------------
    Hnn, Hnp, Hnm = SL.linear_blocks()
    m0 = float(SL.mass_fold(np.zeros(1))[0])

    def band(q):
        Hq = Hnn + Hnp * np.exp(1j * q) + Hnm * np.exp(-1j * q)
        Hq = (Hq + Hq.conj().T) / 2.0
        w = np.sqrt(np.clip(np.linalg.eigvalsh(np.diag(1.0 / np.array([MU, m0])) @ Hq).real, 0, None))
        return float(np.max(w))   # both branches coincide; either
    hq = 1e-3
    wpp = (band(Q0 + hq) - 2 * band(Q0) + band(Q0 - hq)) / hq ** 2
    wp = (band(Q0 + hq) - band(Q0 - hq)) / (2 * hq)
    P = coef_p()
    out["r2"] = (wpp, wp, P)
    A(("R2 P IS THE BAND'S CURVATURE, not a fit: omega''(q0)/2 from the "
       "chain's own linear blocks matches -sin(q0/2)/8, and the group velocity "
       "matches cos(q0/2)/2",
       abs(wpp / 2.0 - P) < 1e-6 and abs(wp - vgroup()) < 1e-6,
       f"omega''/2 numerical {wpp / 2:.7f}, closed form {P:.7f}; v_g "
       f"{wp:.7f} vs {vgroup():.7f}"))

    # ---- R3: Q_ring, quadratic in amplitude, its closed form, and Q_x,ring
    qr, s1, s2 = ring_coefficient()
    qx, x1, x2 = ring_coefficient(cross=True)
    closed = {}
    for q in (np.pi / 3, 2 * np.pi / 3, np.pi / 4):
        closed[q] = (ring_coefficient(q, T=400.0)[0], q_ring_closed(q))
    cerr = max(abs(m / c - 1.0) for (m, c) in closed.values())
    out["r3"] = (qr, s1, s2, qx, x1, x2, closed)
    A(("R3 A UNIFORM WAVE ON A RING RUNS FAST BY Q_ring a^2 -- the shift per "
       "a^2 the same at one and two degrees -- and across four carriers the "
       "measurement is (1 + 4 sin^2(q/2)) / (12 sin^3(q/2)), the closed form "
       "R9 derives from the energy. A counter-propagating wave adds "
       "Q_x,ring b^2, larger than the self shift because the pair also stands "
       "a static strain at the zone boundary",
       abs(s2 / s1 - 1.0) < 0.01 and abs(x2 / x1 - 1.0) < 0.05
       and qx > qr and abs(qr / q_ring_closed() - 1.0) < 2e-3 and cerr < 6e-3,
       f"Q_ring(pi/2) = {qr:.4f} (per a^2: {s1:.4f} at 1 deg, {s2:.4f} at 2 "
       f"deg), closed form {q_ring_closed():.4f}; "
       + ", ".join(f"q = {q:.3f}: {m:.4f} vs {c:.4f}" for q, (m, c) in closed.items())
       + f"; Q_x,ring = {qx:.3f} ({x1:.3f}, {x2:.3f})"))

    # ---- R4: the dent -------------------------------------------------------
    n4, st4, smp4 = 400, 6000, 500
    a4 = np.radians(SOL_AMP_DEG)
    g0, gd0 = packet(n4, 100.0, SOL_W, SOL_AMP_DEG)
    z = np.zeros(n4)
    rec, _ = rk4(z, g0, z, gd0, DT, st4, sample=smp4)
    pred = dent_depth(a4 * a4)
    ratios, track = [], []
    for r in rec:
        if r[0] < 100.0:
            continue
        eps = np.roll(r[1], -1) - r[1]
        e = envelope(r[2])
        i = int(np.argmax(e))
        j = int(np.argmin(eps))
        ratios.append(eps[i] / pred)
        track.append(abs(((j - i + n4 // 2) % n4) - n4 // 2))
    ratios = np.array(ratios)
    out["r4"] = (pred, float(ratios.mean()), float(ratios.std()), max(track))
    A(("R4 THE PACKET DRAGS A DENT: the strain under the peak is a contraction "
       "of -STEP <g^2> / (2 - M_u v_g^2), the closed form from R1 and the "
       "strain chain, and the dent's minimum rides with the peak. It is "
       "finite because v_g < c: the same mismatch that detunes the long-wave "
       "resonance is what bounds the dent",
       abs(ratios.mean() - 1.0) < 0.12 and max(track) <= 3,
       f"predicted eps = {pred:+.5f}; measured / predicted = {ratios.mean():.3f} "
       f"+- {ratios.std():.3f} over {len(ratios)} frames after the dent forms "
       f"(the packet breathes); dent minimum within {max(track)} cells of the "
       f"peak"))

    # ---- R5: THE PREDICTION -------------------------------------------------
    qd = coef_q_dent()
    Q = -(qr + qd)
    aw_pred = np.degrees(soliton_aw(P, Q))
    aw_ring = np.degrees(soliton_aw(P, -qr))
    aw_quartic = np.degrees(soliton_aw(P, -0.75 / omega(Q0)))
    w5, n5, T5 = 16.0, 800, 1200.0
    st5 = int(round(T5 / DT))
    finals = {}
    for adeg in (2.8, 3.1, 3.4):
        g0, gd0 = packet(n5, 160.0, w5, adeg)
        z = np.zeros(n5)
        rec5, _ = rk4(z, g0, z, gd0, DT, st5, sample=st5 // 4)
        hh = [peak_width(envelope(r[2]))[0] for r in rec5]
        finals[adeg] = float(np.log(hh[-1] / hh[0]))
    ks = sorted(finals)
    lg = np.array([finals[k] for k in ks])
    # the stationary amplitude: where the log-height change crosses zero
    a_star = float(np.interp(0.0, lg, ks)) if lg[0] < 0 < lg[-1] else float("nan")
    aw_chain = a_star * w5
    out["r5"] = (aw_pred, aw_chain, finals, aw_ring, aw_quartic, qd, Q)
    A(("R5 THE ENVELOPE EQUATION PREDICTS THE SOLITON FAMILY WITH NOTHING "
       "FITTED: amplitude x width = sqrt(2P/Q), and the exact chain, asked "
       "which amplitude at width 16 neither grows nor decays over a long run, "
       "agrees to a few percent. Two-sided: the on-site quartic alone and the "
       "ring alone predict families half the size, so the dent is load-bearing",
       abs(aw_chain / aw_pred - 1.0) < 0.05 and aw_ring < 0.7 * aw_pred
       and aw_quartic < 0.6 * aw_pred,
       f"Q = -({qr:.4f} ring + {qd:+.4f} dent) = {Q:.4f}, P = {P:.5f}: "
       f"predicted A x w = {aw_pred:.1f} deg-cells; the chain's stationary "
       f"amplitude at w = 16 is {a_star:.2f} deg (log height change "
       + ", ".join(f"{k:.1f}: {finals[k]:+.3f}" for k in ks)
       + f") -> {aw_chain:.1f}. Quartic alone would say {aw_quartic:.1f}, "
       f"ring alone {aw_ring:.1f}"))

    # ---- R6: the reduced model reproduces the chain's envelope -------------
    n6, T6 = 360, 400.0
    st6 = int(round(T6 / DT))
    g0, gd0 = packet(n6, 25.0, SOL_W, SOL_AMP_DEG)
    z = np.zeros(n6)
    vg = vgroup()
    cmp6 = {}
    for lin in (False, True):
        t0 = time.time()
        rec6, _ = rk4(z, g0, z, gd0, DT, st6, sample=st6, linear=lin)
        t_chain = time.time() - t0
        e = envelope(rec6[-1][2])
        hc, wc, ic = peak_width(e)
        cc = centroid(e, ic)
        # the NLS on ONE point per cell, same dt, so the comparison is of
        # the equations and not of the grids (R7 is the coarse one)
        A0 = sech_envelope(n6, n6, 25.0, SOL_W, np.radians(SOL_AMP_DEG))
        t0 = time.time()
        _o, Af, _b = nls_run(A0, n6, P, 0.0 if lin else Q, vg, DT, st6)
        t_nls = time.time() - t0
        ea = np.abs(Af)
        ha, wa, ia = peak_width(ea)
        ca = centroid(ea, ia)
        cmp6["lin" if lin else "nl"] = (hc, wc, cc, ha, wa, ca, t_chain, t_nls)
    nl, li = cmp6["nl"], cmp6["lin"]
    dpos = lambda c: abs(((c[5] - c[2] + n6 / 2) % n6) - n6 / 2)
    out["r6"] = cmp6
    A(("R6 THE ENVELOPE MODEL REPRODUCES THE CHAIN'S OWN ENVELOPE on the "
       "soliton lane's scene after 140 cells of travel: same height, same "
       "width, same position -- and with Q = 0 it reproduces the HARMONIC "
       "chain's dispersing packet, which checks P on its own. Envelope "
       "solitons are where the envelope description is exact enough to "
       "predict the lattice's",
       abs(nl[3] / nl[0] - 1) < 0.06 and abs(nl[4] / nl[1] - 1) < 0.15
       and dpos(nl) < 2.0
       and abs(li[3] / li[0] - 1) < 0.06 and abs(li[4] / li[1] - 1) < 0.15
       and dpos(li) < 2.0,
       f"soliton: chain height {np.degrees(nl[0]):.3f} deg width {nl[1]:.0f} "
       f"at {nl[2]:.1f}; NLS {np.degrees(nl[3]):.3f} deg width {nl[4]:.0f} at "
       f"{nl[5]:.1f}. harmonic: chain {np.degrees(li[0]):.3f} deg width "
       f"{li[1]:.0f} at {li[2]:.1f}; linear Schroedinger {np.degrees(li[3]):.3f} "
       f"deg width {li[4]:.0f} at {li[5]:.1f}"))

    # ---- R7: the coarse grid, and what it costs ----------------------------
    A0 = sech_envelope(n6, n6, 25.0, SOL_W, np.radians(SOL_AMP_DEG))
    _o, Afine, _b = nls_run(A0, n6, P, Q, vg, DT, st6)
    nc, dtc = n6 // 2, 10 * DT
    A0c = sech_envelope(nc, n6, 25.0, SOL_W, np.radians(SOL_AMP_DEG))
    t0 = time.time()
    _o, Acoarse, _b = nls_run(A0c, n6, P, Q, vg, dtc, int(round(T6 / dtc)))
    t_coarse = time.time() - t0
    fine_on_coarse = np.abs(Afine)[::2]
    dev = float(np.abs(np.abs(Acoarse) - fine_on_coarse).max() / np.abs(Afine).max())
    speed = nl[6] / t_coarse
    out["r7"] = (dev, t_coarse, nl[6], speed, nc, dtc)
    A(("R7 THE COARSE GRID: one point per two cells and ten times the chain's "
       "time step give the same envelope as the fine NLS to a fraction of a "
       "percent, and the whole scene costs a hundredth of the chain -- the "
       "order of magnitude the hierarchy page promised, exceeded",
       dev < 0.01 and speed > 20.0,
       f"coarse ({nc} points, dt {dtc:.2f}) vs fine NLS: {dev:.1e} of peak; "
       f"chain {nl[6]:.2f} s, coarse envelope {t_coarse * 1000:.1f} ms -- "
       f"{speed:.0f}x"))

    # ---- R8: the two walls, seen from rung 4 -----------------------------
    def hold_chain(w, amp, q0):
        vgq = vgroup(q0)
        tmax = 12.0 * w / vgq
        nn = int(min(1400, max(200, 2.4 * vgq * tmax + 14 * w)))
        st = int(round(tmax / DT))
        g0, gd0 = packet(nn, nn * 0.12, w, amp, q0=q0)
        z = np.zeros(nn)
        rec, _ = rk4(z, g0, z, gd0, DT, st, sample=max(1, st // 6))
        hh = [peak_width(envelope(r[2]))[0] for r in rec]
        return min(hh) / hh[0], nn, tmax

    def hold_nls(w, amp, q0, nn, tmax, per_cell=1):
        Pq, Qq, vgq, _w = nls_coefficients(q0)
        A0 = sech_envelope(nn * per_cell, nn, nn * 0.12, w, np.radians(amp))
        st = int(round(tmax / dtc))
        o, _a, _b = nls_run(A0, nn, Pq, Qq, vgq, dtc, st, sample=max(1, st // 6))
        hh = [float(s[1].max()) for s in o]
        return min(hh) / hh[0]

    amp2 = aw_pred / 2.0
    kc, nn2, tm2 = hold_chain(2.0, amp2, Q0)
    kn = hold_nls(2.0, amp2, Q0, nn2, tm2, per_cell=4)
    out["r8a"] = (kc, kn)
    A(("R8a THE DISCRETENESS WALL IS INVISIBLE FROM RUNG 4, measured: the "
       "envelope model has no lattice, so at w = 2 -- its own family member, "
       "on a grid fine enough to resolve it -- it holds its height where the "
       "chain loses a fifth (soliton S4). Where the chain fails and the "
       "envelope does not is where the scale-free description stops",
       kn > 0.95 and kc < 0.9,
       f"w = 2, A = {amp2:.1f} deg: chain keeps {kc:.2f}, envelope model "
       f"{kn:.2f}"))

    # the long-carrier side: the dent grows like 1/sin^3(q/2), takes
    # 2w/(c - v_g) to form, and the adiabatic Q crosses zero at pi/3
    q4 = np.pi / 4
    qt3 = q_total_closed(np.pi / 3)
    qt4 = q_total_closed(q4)
    meas3 = ring_coefficient(np.pi / 3)[0] + coef_q_dent(np.pi / 3)
    n8, st8 = 500, 12000
    a8 = np.radians(SOL_AMP_DEG)
    g0, gd0 = packet(n8, 100.0, SOL_W, SOL_AMP_DEG, q0=q4)
    z = np.zeros(n8)
    rec8, _ = rk4(z, g0, z, gd0, DT, st8, sample=2000)
    pred4 = dent_depth(a8 * a8, q4)
    form, hts8 = [], []
    for r in rec8:
        e = envelope(r[2])
        eps = np.roll(r[1], -1) - r[1]
        wgt = e * e
        form.append(float((eps * wgt).sum() / wgt.sum()
                          / (pred4 * (wgt * wgt).sum() / wgt.sum() / (a8 * a8))))
        hts8.append(float(e.max() / a8))
    t_form2 = 2 * SOL_W / (0.5 - vgroup(Q0))
    t_form4 = 2 * SOL_W / (0.5 - vgroup(q4))
    out["r8b"] = (qt3, meas3, qt4, form, hts8, t_form2, t_form4)
    A(("R8b THE RESONANCE WALL, SEEN FROM RUNG 4: the dent grows like "
       "1/sin^3(q/2) and the adiabatic Q = (1 - 2 cos q) / (12 sin^3(q/2)) "
       "VANISHES at carrier pi/3 (ring 4/3 against dent -4/3, measured) and "
       "changes sign below it -- yet at pi/4 the chain packet does not "
       "defocus, because the dent has not formed: it takes 2w/(c - v_g), "
       "which diverges as v_g -> c, and at t = 600 it is still half-built. "
       "The single envelope assumes its long-wave fields follow it instantly, "
       "and that is exactly what fails at the resonance; beyond it the strain "
       "is a field of its own (and the focusing that follows is the mean-fold "
       "instability, longwave L6, not a dent effect)",
       abs(meas3) < 0.01 and qt4 < 0 and form[-1] < 0.75
       and max(hts8) > 1.1,
       f"Q(pi/3): closed form {qt3:.4f}, measured {meas3:+.4f}; Q(pi/4) = "
       f"{qt4:+.3f} (defocusing if adiabatic); dent formation time 2w/(c - v_g) "
       f"= {t_form2:.0f} at pi/2, {t_form4:.0f} at pi/4; at pi/4 the dent is "
       f"{form[-1]:.2f} formed at t = 600 (pi/2 reaches {out['r4'][1]:.2f} by "
       f"200) and the packet peaks at {max(hts8):.2f} of launch -- focusing"))

    # ---- R9: the closed form is DERIVED from the chain's own energy -------
    def cd(f, h):
        return (f(h) - f(-h)) / (2 * h)

    def cd2(f, h):
        return (f(h) - 2 * f(0.0) + f(-h)) / (h * h)

    def rich(fn, h):
        return (4 * fn(h / 2) - fn(h)) / 3.0

    m1 = float(SL.mass_fold(np.zeros(1), deriv=True)[1][0])
    m2 = rich(lambda h: cd2(lambda x: float(SL.mass_fold(np.array([x]))[0]), h), 1e-3)
    N9, k9 = 11, 5

    def Vg(**kw):
        g = np.zeros(N9)
        for idx, val in kw.items():
            g[k9 + int(idx[1:]) * (1 if idx[0] == "p" else -1) if idx != "n" else k9] += val
        return energy_grad(np.zeros(N9), g)[0]

    def d4_on(h):
        return (Vg(n=2 * h) - 4 * Vg(n=h) + 6 * Vg() - 4 * Vg(n=-h) + Vg(n=-2 * h)) / h ** 4

    def d4_22(h):
        def f(a, b):
            g = np.zeros(N9)
            g[k9] += a
            g[k9 + 1] += b
            return energy_grad(np.zeros(N9), g)[0]
        return (f(h, h) - 2 * f(h, 0) + f(h, -h) - 2 * (f(0, h) - 2 * f(0, 0) + f(0, -h))
                + f(-h, h) - 2 * f(-h, 0) + f(-h, -h)) / h ** 4

    def d4_31(h):
        def f(a, b):
            g = np.zeros(N9)
            g[k9] += a
            g[k9 + 1] += b
            return energy_grad(np.zeros(N9), g)[0]
        return ((f(2 * h, h) - f(2 * h, -h)) - 2 * (f(h, h) - f(h, -h))
                + 2 * (f(-h, h) - f(-h, -h)) - (f(-2 * h, h) - f(-2 * h, -h))) / (4 * h ** 4)

    def d3_nnn(h):
        return (Vg(n=2 * h) - 2 * Vg(n=h) + 2 * Vg(n=-h) - Vg(n=-2 * h)) / (2 * h ** 3)

    def d3_nnp(h):
        def f(a, b):
            g = np.zeros(N9)
            g[k9] += a
            g[k9 + 1] += b
            return energy_grad(np.zeros(N9), g)[0]
        return ((f(h, h) - 2 * f(0, h) + f(-h, h)) - (f(h, -h) - 2 * f(0, -h) + f(-h, -h))) / (2 * h ** 3)

    c40 = rich(d4_on, 2e-2) / 24.0
    c22 = rich(d4_22, 2e-2) / 4.0
    c31 = rich(d4_31, 2e-2) / 6.0
    k3n = rich(d3_nnn, 1e-2)
    k3p = rich(d3_nnp, 1e-2)
    coeffs = {"M1 + 4 STEP": m1 + 4 * STEP, "M2 - 32/3": m2 - 32.0 / 3.0,
              "c40 - 17/6": c40 - 17.0 / 6.0, "c22 - 1/2": c22 - 0.5, "c31 - 1/3": c31 - 1.0 / 3.0,
              "V_ggg + 3 STEP": k3n + 3 * STEP, "V_ggg' - STEP/2": k3p - STEP / 2}
    cworst = max(abs(v) for v in coeffs.values())
    derived, dworst = {}, 0.0
    for q, meas in [(Q0, qr)] + [(q, m) for q, (m, _c) in closed.items()]:
        pieces = q_ring_pieces(q)
        w2 = sum(pieces) / (16.0 * np.sin(q / 2.0))
        derived[q] = (w2, meas, pieces)
        dworst = max(dworst, abs(w2 / meas - 1.0))
    # the identity: the four pieces sum to 4(1 + 4 s^2)/(3 s^2)
    idw = max(abs(sum(q_ring_pieces(q)) - 4 * (1 + 4 * np.sin(q / 2) ** 2) / (3 * np.sin(q / 2) ** 2))
              for q in np.linspace(0.3, 3.0, 28))
    # and the strain's second harmonic vanishes at pi/2: its source is sin 2q
    nn9 = 64
    nn = np.arange(nn9)
    a9 = np.radians(3.0)
    rec9, _ = rk4(np.zeros(nn9), a9 * np.cos(Q0 * nn), np.zeros(nn9), a9 * omega(Q0) * np.sin(Q0 * nn),
                  DT, 2000, sample=200)
    u2q = max(float(np.abs(np.fft.fft(r[1])[nn9 // 2]) * 2 / nn9) for r in rec9) / a9 ** 2
    out["r9"] = (coeffs, derived, idw, u2q)
    pp = derived[Q0][2]
    A(("R9 THE CLOSED FORM IS DERIVED, NOT MATCHED: Lindstedt on the chain's "
       "plane wave with the Lagrangian's own coefficients -- fold inertia "
       "8 - 4 STEP g + (16/3) g^2, fold cubic -(STEP/4)(dg)^2 (g + g'), "
       "quartics 17/6, 1/2, 1/3, all measured here by central differences -- "
       "gives four cos(theta) contributions at third order whose sum is "
       "4(1 + 4 s^2)/(3 s^2), an identity, so Q_ring = (1 + 4 s^2)/(12 s^3); "
       "the DC fold source cancels (M1 = -4 STEP, the breathe is not driven) "
       "and the strain's second harmonic is sourced by sin 2q, so at pi/2 "
       "there is none",
       cworst < 2e-3 and dworst < 3e-3 and idw < 1e-12 and u2q < 1e-6,
       "coefficients off by at most " + f"{cworst:.1e}; derived vs measured "
       + ", ".join(f"q = {q:.3f}: {w2:.4f} vs {m:.4f}" for q, (w2, m, _p) in derived.items())
       + f"; pieces at pi/2 (mass, fold cubic, strain 2q, quartic) = "
       f"{pp[0]:+.4f}, {pp[1]:+.4f}, {pp[2]:+.4f}, {pp[3]:+.4f}; identity residual "
       f"{idw:.1e}; strain at 2q on the pi/2 ring {u2q:.1e} a^2"))

    # ---- S1: the collision -- envelope model, warm chain, cold chain ------
    Qx = -(qx + qd)
    n9, T9 = 420, 500.0
    xa, xb = 110.0, 310.0
    st9 = int(round(T9 / dtc))
    smp9 = st9 // 50
    A0 = sech_envelope(n9 // 2, n9, xa, SOL_W, np.radians(SOL_AMP_DEG))
    B0 = sech_envelope(n9 // 2, n9, xb, SOL_W, np.radians(SOL_AMP_DEG))
    both, _a, _b = nls_run(A0, n9, P, Q, vg, dtc, st9, sample=smp9, B0=B0, Qx=Qx)
    free, _a, _b = nls_run(A0, n9, P, Q, vg, dtc, st9, sample=smp9)
    esh, ehs = [], []
    for (t, ea, eb), (_t, ef, _n) in zip(both, free):
        ia, ib = int(np.argmax(ea)), int(np.argmax(eb))
        sep = abs(((ia - ib + n9 // 4) % (n9 // 2)) - n9 // 4) * 2.0
        if t < T9 * 0.55 or sep < 100:
            continue
        ca = centroid(ea, ia, 25) * 2.0
        cf = centroid(ef, int(np.argmax(ef)), 25) * 2.0
        esh.append(((ca - cf + n9 / 2) % n9) - n9 / 2)
        ehs.append(float(ea.max() / ef.max()))
    esh = np.array(esh)

    # the chain, launched WARM (dents in place) and COLD (soliton S3's way)
    ga, va = packet(n9, xa, SOL_W, SOL_AMP_DEG, right=True)
    gb, vb = packet(n9, xb, SOL_W, SOL_AMP_DEG, right=False)
    a9 = np.radians(SOL_AMP_DEG)
    st9c = int(round(T9 / DT))
    smp9c = 40
    z9 = np.zeros(n9)

    def chain_shift(warm):
        if warm:
            ua, uda = dent_fields(ga, a9, True)
            ub, udb = dent_fields(gb, a9, False)
        else:
            ua = uda = ub = udb = z9
        runs = [rk4(u, g, ud, gd, DT, st9c, sample=smp9c)[0] for (u, g, ud, gd) in
                ((ua + ub, ga + gb, uda + udb, va + vb), (ua, ga, uda, va), (ub, gb, udb, vb))]
        sh, hs = [], []
        for rb, rf, rl in zip(*runs):
            if rb[0] < T9 * 0.55:
                continue
            er, el = envelope(rf[2]), envelope(rl[2])
            cr = centroid(er, int(np.argmax(er)), 60)
            cl = centroid(el, int(np.argmax(el)), 60)
            if abs(cr - cl) < 95:
                continue
            ec = envelope(rb[2])
            lo = int(round(cr)) - 45
            wmask = np.zeros(n9, bool)
            wmask[[(lo + j) % n9 for j in range(90)]] = True
            i = int(np.argmax(np.where(wmask, ec, -1.0)))
            sh.append(centroid(ec, i, 45) - cr)
            hs.append(ec[i] / er.max())
        return np.array(sh), float(np.mean(hs))

    warm, warm_h = chain_shift(True)
    cold, cold_h = chain_shift(False)
    out["s1"] = (Qx, float(esh.mean()), float(esh.std()), float(np.mean(ehs)),
                 float(warm.mean()), float(warm.std()), float(cold.mean()), float(cold.std()))
    A(("S1 THE COLLISION, AND A CORRECTION. The envelope model's two "
       "counter-propagating solitons, coupled only by cross-phase modulation "
       "with Q_x from the ring plus the dent the other one carries, pass "
       "through and come out displaced by less than a tenth of a cell. The "
       "chain agrees -- WHEN ITS PACKETS ARE LAUNCHED WITH THEIR DENTS. "
       "Launched cold, as soliton S3 launched them, each packet radiates "
       "strain at c = 1/2 while its dent forms, and the 0.63 cells S3 "
       "published is mostly that: the launch transient, not the collision "
       "(the strain's initial state is the only difference). The collision's own "
       "shift is the small one, and rung 4 predicts it",
       abs(np.mean(ehs) - 1.0) < 0.05 and abs(warm_h - 1.0) < 0.05
       and abs(warm.mean() - esh.mean()) < 0.04 and cold.mean() > 0.4
       and esh.mean() < 0.15,
       f"Q_x = {Qx:.3f} (ring {qx:.3f}, dent {qd:+.3f}); envelope model "
       f"{esh.mean():+.3f} +- {esh.std():.3f} cells; chain warm "
       f"{warm.mean():+.3f} +- {warm.std():.3f}; chain cold {cold.mean():+.3f} "
       f"+- {cold.std():.3f} (S3 published +0.630 +- 0.067); heights "
       f"{np.mean(ehs):.3f} / {warm_h:.3f} of free"))

    out["elapsed"] = time.time() - t_start
    return checks, out


# --------------------------------------------------------------------------
# page export
# --------------------------------------------------------------------------

def export(path="analysis/.pages/data/envelope.json", page="docs/envelope.html"):
    """The page's data: the coefficients, the chain's envelope carpet on the
    soliton scene (the envelope model is integrated IN the page and checked
    against this), the family scan, the dent, and the collision."""
    P, Q, vg, w0 = nls_coefficients()
    qr, _s1, _s2 = ring_coefficient()
    qx, _x1, _x2 = ring_coefficient(cross=True)
    qd = coef_q_dent()
    data = {"step": STEP, "q0": Q0, "dt": DT, "amp_deg": SOL_AMP_DEG,
            "width": SOL_W, "vg": vg, "omega0": w0, "P": P, "Q": Q,
            "Qx": -(qx + qd), "q_ring": qr, "q_dent": qd, "qx_ring": qx,
            "aw_pred_deg": float(np.degrees(soliton_aw(P, Q))),
            "aw_ring_deg": float(np.degrees(soliton_aw(P, -qr))),
            "aw_quartic_deg": float(np.degrees(soliton_aw(P, -0.75 / w0)))}

    # the soliton scene: chain envelope, sampled for a carpet
    n, T, smp = 360, 900.0, 3.0
    st = int(round(T / DT))
    g0, gd0 = packet(n, 25.0, SOL_W, SOL_AMP_DEG)
    z = np.zeros(n)
    for name, lin in (("chainNL", False), ("chainLin", True)):
        rec, _ = rk4(z, g0, z, gd0, DT, st, sample=int(round(smp / DT)), linear=lin)
        E = np.array([envelope(r[2]) for r in rec])
        eps = np.array([np.roll(r[1], -1) - r[1] for r in rec])
        data[name] = {"times": [round(float(r[0]), 3) for r in rec],
                      "env_deg": np.round(np.degrees(E), 3).tolist()}
        if not lin:
            data[name]["strain"] = np.round(eps, 5).tolist()
    data["n_one"] = n
    data["x0"] = 25.0
    data["dent_pred"] = dent_depth(np.radians(SOL_AMP_DEG) ** 2)

    # the family scan at w = 16
    w5, n5, T5 = 16.0, 800, 1200.0
    st5 = int(round(T5 / DT))
    fam = {}
    for adeg in (2.8, 3.1, 3.4):
        g0, gd0 = packet(n5, 160.0, w5, adeg)
        z5 = np.zeros(n5)
        rec5, _ = rk4(z5, g0, z5, gd0, DT, st5, sample=st5 // 24)
        fam[str(adeg)] = {"times": [round(float(r[0]), 2) for r in rec5],
                          "height": [round(float(np.degrees(peak_width(envelope(r[2]))[0])), 4)
                                     for r in rec5]}
    data["family"] = fam
    data["family_w"] = w5

    # the collision, chain side, launched warm and cold: envelope carpets
    # and the position-minus-free trace of the right mover
    n2, T2 = 420, 500.0
    st2 = int(round(T2 / DT))
    smp2 = int(round(smp / DT))
    xa, xb = 110.0, 310.0
    ga, va = packet(n2, xa, SOL_W, SOL_AMP_DEG, right=True)
    gb, vb = packet(n2, xb, SOL_W, SOL_AMP_DEG, right=False)
    a2 = np.radians(SOL_AMP_DEG)
    z2 = np.zeros(n2)
    for name, warm in (("collideWarm", True), ("collideCold", False)):
        if warm:
            ua, uda = dent_fields(ga, a2, True)
            ub, udb = dent_fields(gb, a2, False)
        else:
            ua = uda = ub = udb = z2
        runs = [rk4(u, g, ud, gd, DT, st2, sample=smp2)[0] for (u, g, ud, gd) in
                ((ua + ub, ga + gb, uda + udb, va + vb), (ua, ga, uda, va), (ub, gb, udb, vb))]
        E = np.array([envelope(r[2]) for r in runs[0]])
        off = []
        for rb, rf, rl in zip(*runs):
            er, el = envelope(rf[2]), envelope(rl[2])
            cr = centroid(er, int(np.argmax(er)), 60)
            cl = centroid(el, int(np.argmax(el)), 60)
            if abs(cr - cl) < 95:
                off.append(None)
                continue
            ec = envelope(rb[2])
            lo = int(round(cr)) - 45
            wmask = np.zeros(n2, bool)
            wmask[[(lo + j) % n2 for j in range(90)]] = True
            i = int(np.argmax(np.where(wmask, ec, -1.0)))
            off.append(round(float(centroid(ec, i, 45) - cr), 4))
        data[name] = {"times": [round(float(r[0]), 3) for r in runs[0]],
                      "env_deg": np.round(np.degrees(E), 3).tolist(),
                      "shift": off}
    data["n_two"] = n2
    data["xa"], data["xb"] = xa, xb
    data["chain_shift_published"] = [0.630, 0.067]

    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data))
    inline(path, page)
    return p, data


def inline(path="analysis/.pages/data/envelope.json", page="docs/envelope.html"):
    """Embed the exported JSON in the page's data block (publish-to-docs:
    the page is one self-contained file)."""
    pg = pathlib.Path(page)
    if not pg.exists():
        return False
    txt = pathlib.Path(path).read_text()
    html = pg.read_text()
    new = re.sub(r'(<script type="application/json" id="envdata">).*?(</script>)',
                 lambda m: m.group(1) + txt + m.group(2), html, count=1, flags=re.S)
    pg.write_text(new)
    return True


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "export":
        p, d = export(*sys.argv[2:4])
        print(f"exported envelope data -> {p} ({len(d['chainNL']['times'])} frames)")
        return 0
    if len(sys.argv) > 1 and sys.argv[1] == "inline":
        print("inlined" if inline(*sys.argv[2:4]) else "no page to inline into")
        return 0
    checks, out = gate()
    fails = 0
    for name, ok, detail in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}  [{detail}]")
        fails += 0 if ok else 1
    print(f"{len(checks) - fails}/{len(checks)} rows pass  ({out['elapsed']:.0f} s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
