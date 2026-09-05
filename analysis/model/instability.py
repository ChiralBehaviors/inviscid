"""instability -- RUNG 4b: the long-carrier instability, as five fields.

THE QUESTION (bead inviscid-grl, after longwave.py L7): the gapped mean fold
carries the packet instability at a 6 and nothing at a 1.5 or pi/3, and
leaves a fold pulse on a plane-wave ring oscillating where the chain grows
it fourfold. Which fields does the chain use that the averaged
delta^2 |A|^2 theory integrated out, and does the system with them carry
every long-carrier row with the chain's own coefficients?

TWO PROCESSES, TWO MORE FIELDS. Both come from terms longwave.py already
has, resolved at the wavenumbers it averaged over.

  1. THE PLANE WAVE'S FOUR-WAVE DECAY into the second-harmonic band (T2
     [24443]): carrier(q) + carrier(q) -> (2q - k) + fold(k), resonant where
     the mismatch 2w(q) - w(2q) meets the gap sqrt(K_d) a. The fourth field
     is the second-harmonic envelope C on e^{2i(q n - w t)}: its driven part
     is the Lindstedt harmonic g2 A^2 = (STEP/16) A^2, its free part grows
     at 2q - k. The coupling is the resonant quartic A A C^* delta from V4
     with the second harmonic plus the M2 inertia, L_int = -Gc delta
     Re(A^2 C^*), Gc(q) = 9 + 2 cos q + cos 2q + (4/3) sin^2(q/2) (T2 [24444]).
     Measured on a strain-frozen ring, the seeded linear rate is the
     chain's to 10-50 % with nothing fitted (0.0059/0.0075/0.0069 vs
     0.0039/0.0066/0.0072 per t at lambda 200/100/50, pi/4 a 3), and the
     carrier scan peaks where mismatch/gap = 1.

  2. THE PACKET'S STRAIN-MEDIATED GROWTH (T2 [24465] [24472]). Under a warm
     packet the strain that grows is a packet-shaped envelope AT THE
     CARRIER'S OWN WAVENUMBER q, co-moving with the packet and growing with
     the mean fold at one rate (chain: band x500 while the mean fold's goes
     x700; 2q and mean strain bands flat). The fifth field is the strain
     envelope U on e^{i(q n - w t)}, u_n = u_slow + Re[U e^{i theta}]. It is
     resonant because the chain's two linear branches are EXACTLY
     degenerate (soliton R5: onsite 4, coupling -2, mass 8 for both, so
     omega = sin(q/2) and v_g = cos(q/2)/2 for both), and the fold-strain
     cubic V3_ug = (STEP/2) sum g_n^2 (u_{n+1} - u_{n-1}) couples
     delta x A -> U at (q, w) on resonance. The strain branch has no
     nonlinearity of its own (R7), so U is otherwise a free wave with A's
     transport. Newtonian route, g = (A e^{i th} + cc)/2 convention, ONE
     coefficient G = STEP sin q / (MU w) (0.533 at pi/4):

        U_T + v_g U_X - i P U_XX = -G delta A - (G/2) A^* C
        A_T  gets  + G delta U - (G/2) C U^*
        C_T  gets  + (G/4) A U
        MU delta_TT  gets  - STEP sin q  Im(A U^*)

     Energy |A|^2 + |U|^2 + 4 |C|^2 is conserved by these terms (masses 8,
     frequencies w, w, 2w). On a plane wave U follows delta in QUADRATURE,
     so Im(A U^*) = 0 and there is no k -> 0 force on the mean fold: the
     ring shows the drive (eps(q) = 0.30-0.46 delta0, linear in delta,
     independent of a, chain dc rows) and no growth. Under the packet the
     dent's carrier shift and the lowered gap rotate A against U and the
     loop closes -- linear in delta, threshold-free, and gone when the strain
     is frozen, all as the chain has it ([24456]). The same cubic is why live
     strain DAMPS the plane wave's four-wave decay: the seeded fold shuttles
     into U at q +- k instead of into its 2q - k partner.

  Plus the ON-SITE QUARTIC 4 delta^4 per cell (soliton R8) as the force
  -16 delta^3 on the mean fold: a contracted region with no carrier is a
  fold well (-2 STEP eps delta), and the medium BUCKLES there to
  delta^2 = -STEP eps / 8 (chain 2.9 / 6.1 / 10.8 deg for a 3 / a 6 / 3x the
  a 6 dent, [24456]); without the quartic every growing row diverges.

THE SYSTEM (single packet; B of longwave.py not carried here):

    i(A_T + v_g A_X) + P A_XX - Q_ring |A|^2 A - (STEP/8s) eps A - (K_d/s) delta^2 A
        + (STEP/4s) src A - (Gc/4s) delta A^* C - (Gx/4s) |C|^2 A
        - i G delta U + i (G/2) C U^* = 0
    i(C_T + v_2 C_X) + P_2 C_XX + Dl C + S A^2 - (Gc/8 sin q) delta A^2
        + (Gx/4 sin q) |A|^2 C + (STEP/8 sin q) eps C - i (G/4) A U = 0
    i(U_T + v_g U_X) + P U_XX + i G delta A + i (G/2) A^* C = 0
    MU eps_TT = 2 eps_XX + (STEP/2)(|A|^2)_XX + STEP (delta^2)_XX
    MU delta_TT = 2 delta_XX - MU K_d |A|^2 delta - 2 STEP eps delta + lin. source
        - Gc Re(A^2 C^*) - STEP sin q Im(A U^*) - 16 delta^3

with v_2 = w'(2q) = cos(q)/2, P_2 = -sin(q)/8, Dl = 2 w(q) - w(2q),
S = -Dl STEP/16 (so the driven C is g2 A^2), Gx = 4.5 + (cos q cos 2q +
cos q + cos 2q)/2 - (10/3) sin^2(q/2) the |A|^2 |C|^2 cross term (moves
nothing measurable, kept for completeness). `strain=False` freezes eps AND
U (G = 0), as the chain's frozen u freezes every strain band.

WHAT IT CARRIES (I1-I5), against the chain: the ring drive; the kicked
packet's purely growing fold with U co-located, and none with G off; the
packet rates at a 1.5 / 3 / 6 and pi/3 with the pi/2 null; live strain's
damping of the plane wave's seeded decay, the frozen rate the chain's;
the buckled fold in a static dent. RESIDUALS, as measured: the kicked rate
is 1.5x the chain's at a 3, the unkicked a 3 onset is ~100 t late (chain
t 50), pi/3 is 0.7x, U's amplitude is half the chain's at the same fold,
and the a 3 ring carrier scan peaks one step low (0.25 vs 0.29 pi). Open
with the bead: the pi/2 collision displacement (T2 [24240]).

UNITS as longwave.py: chain rows use soliton.py's exact chain; model rows
the split-step below (dX 2, dt 0.25). Chain constants quoted from T2 are
marked as such.
"""
from __future__ import annotations

import sys
import time

import numpy as np

from analysis.model.envelope import (DT, MU, STEP, coef_p, dent_fields, omega,
                                     q_ring_closed, sech_envelope, vgroup)
from analysis.model.longwave import DTM, DXM, dent_of, k_gap
from analysis.model.soliton import accel, packet


def g_c(q0):
    """Gc(q): the resonant A A C^* delta coupling (docstring)."""
    return 9.0 + 2.0 * np.cos(q0) + np.cos(2 * q0) + (4.0 / 3.0) * np.sin(q0 / 2.0) ** 2


def g_x(q0):
    """Gx(q): the |A|^2 |C|^2 cross term."""
    return (4.5 + (np.cos(q0) * np.cos(2 * q0) + np.cos(q0) + np.cos(2 * q0)) / 2.0
            - (10.0 / 3.0) * np.sin(q0 / 2.0) ** 2)


def g_u(q0):
    """G(q) = STEP sin q / (MU w): the one coefficient of the strain envelope."""
    return STEP * np.sin(q0) / (MU * omega(q0))


# --------------------------------------------------------------------------
# the five-field model
# --------------------------------------------------------------------------

def run(A0, L, q0, dt, steps, sample, eps0=None, epsd0=None, d0=None, dd0=None, C0=None, U0=None,
        strain=True, gu=1.0, sat=1.0):
    """Split-step integration of the system in the docstring. Linear parts
    of A, C, U (exact in k) and of eps, delta (harmonic chains, exact in k)
    go by Fourier multipliers; the nonlinear phases by RK4 on (A, C, U);
    the sources by half kicks. `strain=False` freezes eps at its initial
    value and U at zero. `gu` scales G, `sat` the on-site quartic. Returns
    [(t, A, eps, C, delta, U), ...] every `sample` steps."""
    n = len(A0)
    dX = L / n
    k = 2 * np.pi * np.fft.fftfreq(n, d=dX)
    P, vg, s = coef_p(q0), vgroup(q0), omega(q0)
    qr, kd, Gc, Gx = q_ring_closed(q0), k_gap(q0), g_c(q0), g_x(q0)
    G = (gu if strain else 0.0) * g_u(q0)
    sq = np.sin(q0)
    w2, vg2, P2 = np.sin(q0), 0.5 * np.cos(q0), -np.sin(q0) / 8.0
    Dl = 2 * s - w2
    S = -Dl * STEP / 16.0
    la = np.exp(-1j * (vg * k + P * k * k) * dt)
    lc = np.exp(-1j * (vg2 * k + P2 * k * k - Dl) * dt)
    wk = np.abs(k) / 2.0
    cw, sw = np.cos(wk * dt), np.sin(wk * dt)
    prop = np.where(wk > 0, sw / np.where(wk > 0, wk, 1.0), dt)
    coup, coupC = STEP / (8.0 * s), STEP / (8.0 * w2)
    z = np.zeros(n)
    A = np.asarray(A0, complex).copy()
    C = np.zeros(n, complex) if C0 is None else np.asarray(C0, complex).copy()
    U = np.zeros(n, complex) if U0 is None else np.asarray(U0, complex).copy()
    e = z.copy() if eps0 is None else np.asarray(eps0, float).copy()
    ed = z.copy() if epsd0 is None else np.asarray(epsd0, float).copy()
    d = z.copy() if d0 is None else np.asarray(d0, float).copy()
    dd = z.copy() if dd0 is None else np.asarray(dd0, float).copy()
    E, Ed, D, Dd = np.fft.fft(e), np.fft.fft(ed), np.fft.fft(d), np.fft.fft(dd)
    hist, out = [], []

    def lin_source(Ih, Itt):
        return STEP * Itt - (STEP / 4.0) * (-k * k * Ih)

    def d_accel(Ih, Itt, eps, dl, A_, C_, U_):
        # MU delta_TT - 2 delta_XX = this (real space)
        return (np.fft.ifft(lin_source(Ih, Itt)).real - MU * kd * np.fft.ifft(Ih).real * dl
                - 2.0 * STEP * eps * dl - Gc * (A_ * A_ * np.conj(C_)).real
                - STEP * sq * (A_ * np.conj(U_)).imag - 16.0 * sat * dl ** 3)

    def force(A_, C_, U_, eps, dl, src):
        fa = 1j * ((-qr * np.abs(A_) ** 2 - coup * eps - (kd / s) * dl * dl + (STEP / (4.0 * s)) * src
                    - (Gx / (4.0 * s)) * np.abs(C_) ** 2) * A_ - (Gc / (4.0 * s)) * dl * np.conj(A_) * C_)
        fa = fa + G * dl * U_ - (G / 2.0) * C_ * np.conj(U_)
        fc = -1j * (-S * A_ * A_ + (Gc / (8.0 * w2)) * dl * A_ * A_
                    + (Gx / (4.0 * w2)) * np.abs(A_) ** 2 * C_ + coupC * eps * C_) + (G / 4.0) * A_ * U_
        fu = -G * dl * A_ - (G / 2.0) * np.conj(A_) * C_
        return fa, fc, fu

    def nl(h, eps, dl, src):
        nonlocal A, C, U
        k1 = force(A, C, U, eps, dl, src)
        k2 = force(A + 0.5 * h * k1[0], C + 0.5 * h * k1[1], U + 0.5 * h * k1[2], eps, dl, src)
        k3 = force(A + 0.5 * h * k2[0], C + 0.5 * h * k2[1], U + 0.5 * h * k2[2], eps, dl, src)
        k4 = force(A + h * k3[0], C + h * k3[1], U + h * k3[2], eps, dl, src)
        A = A + h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
        C = C + h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
        U = U + h * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6

    def kick(h, Ih, Itt):
        nonlocal Ed, Dd
        eps, dl = np.fft.ifft(E).real, np.fft.ifft(D).real
        if strain:
            Ed = Ed + h * (-(STEP / 2.0) * k * k * Ih - STEP * k * k * np.fft.fft(dl * dl)) / MU
        Dd = Dd + h * np.fft.fft(d_accel(Ih, Itt, eps, dl, A, C, U)) / MU

    for st in range(steps + 1):
        Ih = np.fft.fft(np.abs(A) ** 2)
        hist = (hist + [Ih])[-3:]
        Itt = (hist[-1] - 2 * hist[-2] + hist[-3]) / (dt * dt) if len(hist) == 3 else np.zeros(n, complex)
        eps, dl = np.fft.ifft(E).real, np.fft.ifft(D).real
        if st % sample == 0:
            out.append((st * dt, A.copy(), eps.copy(), C.copy(), dl.copy(), U.copy()))
        if st == steps:
            break
        src = d_accel(Ih, Itt, eps, dl, A, C, U) / MU
        nl(0.5 * dt, eps, dl, src)
        A = np.fft.ifft(np.fft.fft(A) * la)
        C = np.fft.ifft(np.fft.fft(C) * lc)
        U = np.fft.ifft(np.fft.fft(U) * la)
        kick(0.5 * dt, Ih, Itt)
        if strain:
            E, Ed = E * cw + Ed * prop, -E * wk * sw + Ed * cw
        D, Dd = D * cw + Dd * prop, -D * wk * sw + Dd * cw
        kick(0.5 * dt, np.fft.fft(np.abs(A) ** 2), Itt)
        nl(0.5 * dt, np.fft.ifft(E).real, np.fft.ifft(D).real, src)
    return out


# --------------------------------------------------------------------------
# measurements
# --------------------------------------------------------------------------

def _rk4_chain(y, f, steps, every):
    """RK4 on the exact chain, yielding (step, y) every `every` steps."""
    for s_ in range(steps + 1):
        if s_ % every == 0:
            yield s_, y
        if s_ == steps:
            break
        k1 = f(y)
        k2 = f(y + 0.5 * DT * k1)
        k3 = f(y + 0.5 * DT * k2)
        k4 = f(y + DT * k3)
        y = y + (DT / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def ring_drive(q0, adeg, d0deg, N=800, T=300.0):
    """A plane wave (q0, a) on a ring with a UNIFORM fold offset d0, strain
    live from u = 0: the peak amplitude of the strain at the carrier's
    wavenumber, eps(q) = |eps_q|, per radian of d0 -- chain and model."""
    a, d0 = np.radians(adeg), np.radians(d0deg)
    nn = np.arange(N)
    z = np.zeros(N)
    mq = int(round(q0 * N / (2 * np.pi)))

    def f(y):
        u, g, ud, gd = y[:N], y[N:2 * N], y[2 * N:3 * N], y[3 * N:]
        au, ag = accel(u, g, ud, gd)
        return np.concatenate([ud, gd, au, ag])

    y0 = np.concatenate([z, a * np.cos(q0 * nn) + d0, z, a * omega(q0) * np.sin(q0 * nn)])
    ch = 0.0
    for _s, y in _rk4_chain(y0, f, int(round(T / DT)), 100):
        u = y[:N]
        ch = max(ch, 2 * abs(np.fft.fft(np.roll(u, -1) - u)[mq]) / N)
    n = N // 2
    out = run(np.full(n, a, complex), N, q0, DTM, int(round(T / DTM)), 20,
              d0=np.full(n, d0), C0=np.full(n, (STEP / 16.0) * a * a, complex))
    md = max(2 * np.sin(q0 / 2.0) * np.abs(U).mean() for *_r, U in out)
    return ch / d0, md / d0


def packet_track(q0, w, adeg, T, kick_deg=0.0, N=800, nsamp=12, gu=1.0):
    """A warm w-wide packet (C warm, U cold), optionally with a fold kick of
    kick_deg sech^2(X/10) under its peak: (t, h/a, fold at the peak in
    degrees, max |delta| in degrees, its X, the peak's X, eps(q) = 2s|U| max,
    its X) every T/nsamp."""
    a = np.radians(adeg)
    n = N // 2
    X = np.arange(n) * DXM
    A0 = sech_envelope(n, N, 100.0, w, a)
    e, ed = dent_of(np.abs(A0), a, True, q0)
    d0 = np.radians(kick_deg) / np.cosh((X - 100.0) / 10.0) ** 2 if kick_deg else None
    stm = int(round(T / DTM))
    rows = []
    for t, A, _e, _C, dl, U in run(A0, N, q0, DTM, stm, stm // nsamp, e, ed, d0=d0,
                                    C0=(STEP / 16.0) * A0 * A0, gu=gu):
        ea = np.abs(A)
        if not np.all(np.isfinite(ea)):
            break
        ip, j, ju = int(np.argmax(ea)), int(np.argmax(np.abs(dl))), int(np.argmax(np.abs(U)))
        rows.append((t, float(ea.max() / a), float(np.degrees(dl[ip])), float(np.degrees(abs(dl[j]))),
                     j * DXM, ip * DXM, float(2 * np.sin(q0 / 2.0) * np.abs(U).max()), ju * DXM))
    return rows


def _rate(rows, t0, t1):
    r = {row[0]: row for row in rows}
    return float(np.log(r[t1][3] / r[t0][3]) / (t1 - t0))


def seeded_decay(q0, adeg, mk, live, N=800, T=400.0):
    """A plane wave (q0, a) on a ring with a seeded fold mode k = 2 pi mk/N
    and its partner at 2q - k (0.02 deg each), strain live or frozen: the
    fold mode's amplitude every 40 t and its fitted rate over t > 100."""
    a, sd = np.radians(adeg), np.radians(0.02)
    n = N // 2
    X = np.arange(n) * DXM
    k = 2 * np.pi * mk / N
    d0 = sd * np.cos(k * X)
    C0 = (STEP / 16.0) * a * a + sd * np.exp(-1j * k * X)
    stm = int(round(T / DTM))
    ts, am = [], []
    for t, _A, _e, _C, dl, _U in run(np.full(n, a, complex), N, q0, DTM, stm, 160, d0=d0, C0=C0, strain=live):
        ts.append(t)
        am.append(np.degrees(np.abs(np.fft.fft(dl))[mk] * 2 / n))
    ts, la = np.array(ts), np.log(np.array(am))
    i0 = len(ts) // 4
    return float(np.polyfit(ts[i0:], la[i0:], 1)[0]), float(la[-1] - la[0]), am


def buckling(adeg, scale=1.0, q0=np.pi / 4, N=800, T=300.0):
    """The dent a packet of amplitude a would drag, frozen and SCALED, with
    NO carrier, and a 0.3-degree fold kick at its centre: the fold at the
    centre over the last third of the run (degrees, mean and spread) against
    the buckled value sqrt(-STEP eps / 8)."""
    a = np.radians(adeg)
    n = N // 2
    X = np.arange(n) * DXM
    e = scale * dent_of(sech_envelope(n, N, 100.0, 16.0, a), a, True, q0)[0]
    d0 = np.radians(0.3) / np.cosh((X - 100.0) / 10.0) ** 2
    stm = int(round(T / DTM))
    vals = [np.degrees(dl[50]) for t, _A, _e, _C, dl, _U in run(np.zeros(n, complex), N, q0, DTM, stm, 20,
                                                                 e, np.zeros(n), d0=d0, strain=False)
            if t > 2 * T / 3]
    pred = np.degrees(np.sqrt(max(-STEP * e[50] / 8.0, 0.0)))
    return float(np.mean(vals)), float(np.std(vals)), pred


# --------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------

def gate():
    checks, out = [], {}
    A = checks.append
    t_start = time.time()
    q = np.pi / 4

    # ---- I1: the drive ------------------------------------------------------
    ch, md = ring_drive(q, 3.0, 1.0)
    out["i1"] = (ch, md)
    A(("I1 A FOLD OFFSET DRIVES THE STRAIN AT THE CARRIER'S OWN WAVENUMBER: a "
       "uniform delta on a plane-wave ring, strain live, raises a strain wave "
       "at q of amplitude proportional to delta (chain 0.30-0.46 per radian "
       "across delta 0.5-2 deg and a 1.5-6, T2 [24472]) and grows nothing; the "
       "model's one coefficient G = STEP sin q/(MU w) gives it within a third",
       0.25 < md < 0.5 and 0.25 < ch < 0.6 and abs(md / ch - 1.0) < 0.4,
       f"pi/4 a 3 delta 1 deg: eps(q)/delta chain {ch:.3f}, model {md:.3f}"))

    # ---- I2: the kicked packet ----------------------------------------------
    kk = packet_track(q, 16.0, 3.0, 300.0, kick_deg=0.3)
    k0 = packet_track(q, 16.0, 3.0, 300.0, kick_deg=0.3, gu=0.0)
    pk = {r[0]: r for r in kk}
    late = [r[2] for r in kk if 100 <= r[0] <= 225]      # the exponential stage; it saturates past 250
    mono = all(b > a_ for a_, b in zip(late, late[1:])) and late[0] > 0
    rk = float(np.log(pk[225.0][2] / pk[125.0][2]) / 100.0)
    coloc = abs(pk[225.0][7] - pk[225.0][5]) < 20
    g0late = [abs(r[2]) for r in k0 if r[0] >= 100]
    out["i2"] = (kk, k0, rk)
    A(("I2 UNDER A WARM PACKET A FOLD KICK IS PURELY GROWING, NOT GAPPED, and "
       "the strain envelope at q grows with it, co-located: chain 0.03 -> 0.19 "
       "-> 0.34 -> 0.59 -> 0.87 deg at the peak over t 100..225 with no zero "
       "crossing (0.017/t, T2 [24456]); with G off the kick rings in its gap "
       "and decays, as longwave.py's did",
       mono and 0.01 < rk < 0.035 and coloc and pk[275.0][6] > 5 * pk[25.0][6] and max(g0late) < 0.2,
       f"pi/4 a 3 w 16 + 0.3 deg: fold at peak " + " ".join(f"{r[2]:.3f}" for r in kk if 100 <= r[0] <= 225)
       + f" (t 100..225, rate {rk:.3f}/t; saturates past 250); eps(q) {pk[25.0][6]:.1e} -> {pk[275.0][6]:.1e} at X "
       f"{pk[275.0][7]:.0f}, peak at {pk[275.0][5]:.0f}; G off: |fold| max {max(g0late):.3f} after t 100"))

    # ---- I3: the rates across amplitude and carrier -------------------------
    rows = {"a 1.5": (packet_track(q, 16.0, 1.5, 1000.0, nsamp=10), 400.0, 900.0, 0.0063),
            "a 3": (packet_track(q, 16.0, 3.0, 400.0, nsamp=8), 150.0, 350.0, 0.015),
            "a 6": (packet_track(q, 16.0, 6.0, 200.0, nsamp=6), 33.333333333333336, 100.0, 0.04),
            "pi/3 a 3": (packet_track(np.pi / 3, 16.0, 3.0, 1000.0, nsamp=10), 500.0, 900.0, 0.008)}
    rates = {}
    for name, (tr, t0, t1, chain) in rows.items():
        r = {row[0]: row for row in tr}
        tt0 = min(r, key=lambda t: abs(t - t0))
        tt1 = min(r, key=lambda t: abs(t - t1))
        rates[name] = (float(np.log(r[tt1][3] / r[tt0][3]) / (tt1 - tt0)), chain)
    ctl = packet_track(np.pi / 2, 16.0, 3.0, 1200.0)
    cvals = [r[3] for r in ctl if r[0] >= 100]
    out["i3"] = (rates, ctl)
    A(("I3 THE INSTABILITY IS THRESHOLD-FREE AND ORDERED BY THE CARRIER, as the "
       "chain's: a warm w 16 packet's mean fold grows at pi/4 for a 1.5, 3 and "
       "6 and at pi/3, each within a factor two of the chain's rate (0.0063 / "
       "0.015 / 0.04 / 0.008 per t, T2 [24273]), and not at all at pi/2 to "
       "t 1200. longwave.py needed a ~ 2 and had nothing at pi/3",
       all(0.5 < m / c < 2.0 for m, c in rates.values())
       and max(cvals) / min(cvals) < 3.0 and max(cvals) < 0.05,
       "; ".join(f"{k}: {m:.4f}/t (chain {c})" for k, (m, c) in rates.items())
       + f"; pi/2 a 3 mean fold {min(cvals):.4f}..{max(cvals):.4f} deg"))

    # ---- I4: live strain damps the plane wave's decay ------------------------
    dec = {("100", True): seeded_decay(q, 3.0, 8, True), ("100", False): seeded_decay(q, 3.0, 8, False),
           ("33", True): seeded_decay(q, 3.0, 24, True)}
    out["i4"] = dec
    A(("I4 ON A PLANE WAVE THE FOUR-WAVE DECAY RUNS WITH THE STRAIN FROZEN AND "
       "IS DAMPED WITH IT LIVE: a seeded fold mode with its 2q - k partner "
       "grows at the chain's frozen rate (0.0066/t at lambda 100, T2 [24444]) "
       "and, live, shuttles into the strain at q +- k instead -- less than "
       "half the growth at lambda 100, none at lambda 33 (chain: ~0 and "
       "negative, T2 [24456]). longwave.py + C alone had live = frozen",
       0.004 < dec[("100", False)][0] < 0.012 and dec[("100", True)][1] < 0.5 * dec[("100", False)][1]
       and dec[("33", True)][1] < 0.8,
       f"lambda 100: frozen {dec[('100', False)][0]:.4f}/t (ln {dec[('100', False)][1]:+.2f}), live "
       f"{dec[('100', True)][0]:.4f}/t (ln {dec[('100', True)][1]:+.2f}); lambda 33 live "
       f"{dec[('33', True)][0]:.4f}/t (ln {dec[('33', True)][1]:+.2f})"))

    # ---- I5: buckling ---------------------------------------------------------
    bk = {"a 3": buckling(3.0), "a 6": buckling(6.0), "a 6 x3": buckling(6.0, 3.0)}
    chain_bk = {"a 3": 2.9, "a 6": 6.1, "a 6 x3": 10.8}          # T2 [24456]
    out["i5"] = bk
    A(("I5 THE MEDIUM BUCKLES UNDER CONTRACTION, AND THE QUARTIC HOLDS IT: a "
       "static dent with no carrier grows a fold kick at the -2 STEP eps delta "
       "rate and saturates at delta^2 = -STEP eps/8 against the on-site "
       "4 delta^4 -- chain 2.9 / 6.1 / 10.8 deg for the a 3, a 6 and tripled a 6 "
       "dents (T2 [24456]); without the -16 delta^3 term every growing row here "
       "diverged",
       all(abs(v[0] - chain_bk[k]) < 0.3 and abs(v[0] - v[2]) < 0.4 for k, v in bk.items()),
       "; ".join(f"{k}: {v[0]:.2f} +- {v[1]:.2f} deg (buckled {v[2]:.2f}, chain {chain_bk[k]})"
                 for k, v in bk.items())))

    out["elapsed"] = time.time() - t_start
    return checks, out


def main():
    checks, out = gate()
    fails = 0
    for name, ok, detail in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}  [{detail}]")
        fails += 0 if ok else 1
    print(f"{len(checks) - fails}/{len(checks)} rows pass  ({out['elapsed']:.0f} s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
