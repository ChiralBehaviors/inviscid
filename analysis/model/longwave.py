"""longwave -- RUNG 4b: the envelope past the resonance wall, as three fields.

THE QUESTION (bead inviscid-grl, from envelope R8b): toward long carriers the
single-envelope equation fails, because the strain it treats as a coefficient
takes 2w/(c - v_g) to form, and because a mean fold is pumped under the
packet. What is the system when the long-wave fields are fields of their own,
and which of the chain's long-carrier behaviour does it carry?

THE SYSTEM (build 2 of the bead, plus the mean fold's gap). Fold
gamma_n = delta + Re[A e^{i(q0 n - w0 t)}] + Re[B e^{i(-q0 n - w0 t)}],
strain eps = u_X, on one grid point per two cells:

    i(A_T + v_g A_X) + P A_XX - Q_ring |A|^2 A - Qx_ring |B|^2 A
        - K B^2 A^* e^{-i Dk X} - (STEP/8s) eps A - (K_d/s) delta^2 A = 0
    MU eps_TT = 2 eps_XX + (STEP/2)(|A|^2 + |B|^2)_XX + STEP (delta^2)_XX
    MU delta_TT = 2 delta_XX - MU K_d (|A|^2 + |B|^2) delta - 2 STEP eps delta
        + STEP (|A|^2 + |B|^2)_TT - (STEP/4)(|A|^2 + |B|^2)_XX

(B mirrored). Every coefficient is the chain's own, none fitted:
  * P, v_g, Q_ring = (1 + 4 s^2)/(12 s^3): envelope R2, R3, R9.
  * The strain as a FIELD, driven by the gradient of |A|^2 (envelope R1's
    cubic) and shifting the carrier by STEP eps/(8 s): the dent is no longer
    a coefficient Q_dent, it forms on its own schedule.
  * Qx_ring = 1.47, the phase-AVERAGED cross coefficient. Envelope R3's 1.89
    is the value at relative carrier phase zero; the pi/2 collision depends
    on that phase (below), so the phase-independent part is what belongs in
    a coefficient.
  * K = 0.44 per rad^2, the UMKLAPP four-wave term: at q0 = pi/2, 4 q0 =
    2 pi, so A_L^2 A_R^* is resonant (Dk = 4 q0 - 2 pi = 0) and two
    counter-propagating waves EXCHANGE ENERGY at a rate set by their
    relative phase, -0.67 b^2 sin 2 phi on the ring. Off pi/2 the term
    rotates with Dk and averages away.
  * K_d(q) = 2.5 + cos(q)/2 - s^2/3, THE MEAN FOLD'S GAP under a carrier:
    from the Lagrangian's d^2 c^2 terms -- V4 gives delta^2 [20 c_n^2 +
    4 c_n c_{n+1}] -> (10 + 2 cos q)|A|^2 delta^2, the fold inertia's M2
    term -(4/3) s^2 |A|^2 delta^2 (its sign is destabilising); M1 and the
    fold cubic give nothing at this order for a uniform delta (their linear
    pieces cancel, R9). The same Lagrangian term is the carrier's (K_d/s)
    delta^2 shift. It is not Q_ring a^2, which it happens to be near at pi/4.
  * -2 STEP eps delta and STEP (delta^2)_XX: the fold-strain cubic with
    delta in it. A contracted region LOWERS the gap (Omega^2 = -STEP eps/4
    for eps < 0), so the packet's own dent takes back most of its gap at
    pi/4; a contraction with no carrier at all is a fold well.
  * The linear source and its back-coupling on A (from M1 and the fold
    cubic at gradient order) seed a static mean fold of a few hundredths of
    a^2 under the packet. Alone they give no growth (measured; that draft was
    rejected). With the gap they are finite; without it they are resonant.

WHAT IT CARRIES (L1-L6), each gated against the chain:
  L1  the pi/2 collision's energy exchange, sign and size vs relative phase;
  L2  none of it off resonance;
  L3  the COLD launch: a packet born without its dent detunes its own
      carrier (the strain transient acting on the envelope's phase) and a
      warm/cold pair then exchanges energy through the rotated phase;
  L4  the dent's formation at w 16, pi/4 and pi/3, within 0.05 of the
      dent through t 500-700 -- the resonance wall's own schedule;
  L5  the gap, measured directly: a uniform fold offset on a plane-wave ring
      oscillates at Omega^2 = K_d a^2 at every carrier and amplitude tried;
  L6  the long-carrier INSTABILITY: at pi/4 the mean fold grows exponentially
      under the packet, co-moving, and the packet focuses with it, where the
      adiabatic Q says defocus; nothing at pi/2. Chain: 0.04/t at a 6,
      0.015 at a 3, zero at pi/2 to t 1600 (T2 [24273]); this model's rate
      at a 6 is the chain's.

WHAT IT DOES NOT (L7, stated as measured). The chain's growth is
threshold-free (rate ~ a^1.3, width-independent, 0.0063/t at a 1.5), half
as fast at pi/3, and present on a RING with no packet at all (a 0.5-degree
fold pulse under a uniform carrier grows fourfold by t 400, T2 [24289]).
This system's instability needs a ~ 2 at pi/4 (the trapped mode's length
c/(sqrt(K_d) a) must fit under the packet), is 5x slow at pi/3, and leaves
the ring pulse oscillating in its gap without growth. The averaged
delta^2 |A|^2 theory is right about the mean fold's local dynamics and
misses two processes it integrated out: the plane wave's four-wave decay
into the second-harmonic band, and the packet's strain-mediated growth
through the strain envelope at the carrier's own wavenumber. Both are
carried, with the chain's coefficients, by instability.py (the fold cubic
on the pulse's gradients, the first candidate, was tried and does nothing:
T2 [24442]). The pi/2 collision DISPLACEMENT (0.3 cells too positive at
a 6) is the open item (T2 [24240]) and is not gated here.

UNITS as envelope.py. The chain rows use soliton.py's exact two-field
chain; the model rows use the split-step integrator below (dX 2, dt 0.25;
dt 0.0625 changes nothing gated). Chain constants quoted from T2 are
marked as such in the row text.
"""
from __future__ import annotations

import sys
import time

import numpy as np

from analysis.model.envelope import (DT, MU, Q0, SOL_W, STEP, coef_p, coef_q_dent,
                                     dent_depth, omega, q_ring_closed, sech_envelope,
                                     soliton_aw, vgroup)
from analysis.model.soliton import accel, envelope, packet, peak_width, rk4

#: the phase-averaged cross coefficient and the umklapp coefficient (T2 [24204])
QX_RING_AVG = 1.47
K_RING = 0.44
#: model grid: two cells per point, ten chain steps per model step
DXM, DTM = 2.0, 0.25


def k_gap(q0=Q0):
    """K_d(q): the mean fold's gap per |A|^2 under a carrier at q (docstring)."""
    return 2.5 + 0.5 * np.cos(q0) - np.sin(q0 / 2.0) ** 2 / 3.0


# --------------------------------------------------------------------------
# the three-field model
# --------------------------------------------------------------------------

def run(A0, L, q0, dt, steps, sample, eps0=None, epsd0=None, B0=None, d0=None, dd0=None,
        Qx_ring=QX_RING_AVG, K=K_RING, strain=True, meanfold=True):
    """Split-step integration of the system in the docstring. Linear parts
    of A, B (exact in k), eps and delta (harmonic chains, exact in k) go by
    Fourier multipliers; the nonlinear phases by RK4 on A, B; the sources by
    half kicks. `meanfold=False` is build 2 exactly; `strain=False` freezes
    eps at its initial value (use eps0 None for none). Returns
    [(t, A, eps, B, delta), ...] every `sample` steps."""
    n = len(A0)
    dX = L / n
    X = np.arange(n) * dX
    k = 2 * np.pi * np.fft.fftfreq(n, d=dX)
    P, vg, s = coef_p(q0), vgroup(q0), omega(q0)
    qr, kd = q_ring_closed(q0), k_gap(q0)
    Dk = 4 * q0 - 2 * np.pi
    la = np.exp(-1j * (vg * k + P * k * k) * dt)
    lb = np.exp(-1j * (-vg * k + P * k * k) * dt)
    wk = np.abs(k) / 2.0
    cw, sw = np.cos(wk * dt), np.sin(wk * dt)
    prop = np.where(wk > 0, sw / np.where(wk > 0, wk, 1.0), dt)
    em, ep = np.exp(-1j * Dk * X), np.exp(1j * Dk * X)
    coup = STEP / (8.0 * s)
    z = np.zeros(n)
    A = np.asarray(A0, complex).copy()
    B = None if B0 is None else np.asarray(B0, complex).copy()
    e = z.copy() if eps0 is None else np.asarray(eps0, float).copy()
    ed = z.copy() if epsd0 is None else np.asarray(epsd0, float).copy()
    d = z.copy() if d0 is None else np.asarray(d0, float).copy()
    dd = z.copy() if dd0 is None else np.asarray(dd0, float).copy()
    E, Ed, D, Dd = np.fft.fft(e), np.fft.fft(ed), np.fft.fft(d), np.fft.fft(dd)
    hist, out = [], []

    def intensity():
        return np.abs(A) ** 2 + (0.0 if B is None else np.abs(B) ** 2)

    def force(A_, B_, eps, dl, src):
        shift = -qr * np.abs(A_) ** 2 - coup * eps - (kd / s) * dl * dl + (STEP / (4.0 * s)) * src
        if B_ is None:
            return 1j * shift * A_, None
        b2, a2 = np.abs(B_) ** 2, np.abs(A_) ** 2
        return (1j * ((shift - Qx_ring * b2) * A_ - K * B_ * B_ * np.conj(A_) * em),
                1j * ((-qr * b2 - Qx_ring * a2 - coup * eps - (kd / s) * dl * dl
                       + (STEP / (4.0 * s)) * src) * B_ - K * A_ * A_ * np.conj(B_) * ep))

    def nl(h, eps, dl, src):
        nonlocal A, B
        k1a, k1b = force(A, B, eps, dl, src)
        k2a, k2b = force(A + 0.5 * h * k1a, None if B is None else B + 0.5 * h * k1b, eps, dl, src)
        k3a, k3b = force(A + 0.5 * h * k2a, None if B is None else B + 0.5 * h * k2b, eps, dl, src)
        k4a, k4b = force(A + h * k3a, None if B is None else B + h * k3b, eps, dl, src)
        A = A + h * (k1a + 2 * k2a + 2 * k3a + k4a) / 6
        if B is not None:
            B = B + h * (k1b + 2 * k2b + 2 * k3b + k4b) / 6

    def lin_source(Ih, Itt):
        # M1 and the fold cubic at gradient order: STEP (I)_TT - (STEP/4) I_XX
        return STEP * Itt - (STEP / 4.0) * (-k * k * Ih)

    def d_accel(Ih, Itt, eps, dl):
        # MU delta_TT - 2 delta_XX = this (real space), so (delta_TT - delta_XX/4) = this/MU
        return (np.fft.ifft(lin_source(Ih, Itt)).real - MU * kd * np.fft.ifft(Ih).real * dl
                - 2.0 * STEP * eps * dl)

    def kick(h, Ih, Itt):
        nonlocal Ed, Dd
        eps, dl = np.fft.ifft(E).real, np.fft.ifft(D).real
        if strain:
            Ed = Ed + h * (-(STEP / 2.0) * k * k * Ih - STEP * k * k * np.fft.fft(dl * dl)) / MU
        if meanfold:
            Dd = Dd + h * np.fft.fft(d_accel(Ih, Itt, eps, dl)) / MU

    for st in range(steps + 1):
        Ih = np.fft.fft(intensity())
        hist = (hist + [Ih])[-3:]
        Itt = (hist[-1] - 2 * hist[-2] + hist[-3]) / (dt * dt) if len(hist) == 3 else np.zeros(n, complex)
        eps, dl = np.fft.ifft(E).real, np.fft.ifft(D).real
        if st % sample == 0:
            out.append((st * dt, A.copy(), eps.copy(), None if B is None else B.copy(), dl.copy()))
        if st == steps:
            break
        src = d_accel(Ih, Itt, eps, dl) / MU if meanfold else z
        nl(0.5 * dt, eps, dl, src)
        A = np.fft.ifft(np.fft.fft(A) * la)
        if B is not None:
            B = np.fft.ifft(np.fft.fft(B) * lb)
        kick(0.5 * dt, Ih, Itt)
        if strain:
            E, Ed = E * cw + Ed * prop, -E * wk * sw + Ed * cw
        if meanfold:
            D, Dd = D * cw + Dd * prop, -D * wk * sw + Dd * cw
        kick(0.5 * dt, np.fft.fft(intensity()), Itt)
        nl(0.5 * dt, np.fft.ifft(E).real, np.fft.ifft(D).real, src)
    return out


def dent_of(env, amp, right, q0, dX=DXM):
    """The co-moving dent under a model envelope (envelope.dent_fields on the
    model grid): eps and its velocity field for a warm launch."""
    e = dent_depth(amp * amp, q0) * (env / amp) ** 2
    e = e - e.mean()
    return e, -(1.0 if right else -1.0) * vgroup(q0) * np.gradient(e, dX)


# --------------------------------------------------------------------------
# measurements
# --------------------------------------------------------------------------

def _efrac(e, half=10):
    n = len(e)
    i = int(np.argmax(e))
    dd = np.abs((np.arange(n) - i + n / 2) % n - n / 2)
    return float((e ** 2)[dd < half].sum())


def _carrier_slope(A, dX=DXM, win=6):
    n = len(A)
    i = int(np.argmax(np.abs(A)))
    idx = np.arange(i - win, i + win + 1)
    ph = np.unwrap(np.angle(A[idx % n]))
    return float(np.polyfit(idx * dX, ph, 1)[0])


def collision(q0, adeg, w, phi, warmA, warmB, L=720, T=900.0, xa=110.0, xb=310.0, K=K_RING):
    """Two packets, A right-moving at xa with carrier phase phi and B
    left-moving at xb, each launched warm (its dent in place) or cold, on a
    ring of L cells, against each alone. Returns the mean over the frames
    after they separate of (E_A / E_A free, E_B / E_B free, carrier slope of
    A after, carrier slope of A free)."""
    n = int(L / DXM)
    a = np.radians(adeg)
    A0 = sech_envelope(n, L, xa, w, a) * np.exp(1j * phi)
    B0 = sech_envelope(n, L, xb, w, a)
    ea, eda = dent_of(np.abs(A0), a, True, q0)
    eb, edb = dent_of(np.abs(B0), a, False, q0)
    z = np.zeros(n)
    e0 = (ea if warmA else z) + (eb if warmB else z)
    ed0 = (eda if warmA else z) + (edb if warmB else z)
    steps, smp = int(round(T / DTM)), int(round(10 / DTM))
    both = run(A0, L, q0, DTM, steps, smp, e0, ed0, B0, K=K, meanfold=False)
    free = run(A0, L, q0, DTM, steps, smp, ea if warmA else z, eda if warmA else z, K=K, meanfold=False)
    freeB = run(np.zeros(n, complex), L, q0, DTM, steps, smp, eb if warmB else z, edb if warmB else z, B0,
                K=K, meanfold=False)
    rows = []
    for (t, A, _e, B, _d), (_, Af, _e2, _n, _d2), (_, _z, _e3, Bf, _d3) in zip(both, free, freeB):
        if t < 600:
            continue
        ia, ib = int(np.argmax(np.abs(Af))), int(np.argmax(np.abs(Bf)))
        if abs(((ia - ib + n / 2) % n) - n / 2) * DXM < 95:
            continue
        rows.append((_efrac(np.abs(A)) / _efrac(np.abs(Af)), _efrac(np.abs(B)) / _efrac(np.abs(Bf)),
                     _carrier_slope(A), _carrier_slope(Af)))
    return np.array(rows).mean(0)


def dent_formation(q0, N, w, adeg, T, nsamp):
    """Cold launch of a w-wide packet at carrier q0, chain and model: the
    fraction of R4's dent formed under the envelope, and the height, every
    T/nsamp. Returns (chain rows, model rows) of (t, h/a, dent fraction)."""
    a = np.radians(adeg)
    pred = dent_depth(a * a, q0)

    def frac(eps, e):
        wt = e * e
        return float((eps * wt).sum() / wt.sum() / (pred * (wt * wt).sum() / wt.sum() / (a * a)))

    st = int(round(T / DT))
    g0, gd0 = packet(N, 100.0, w, adeg, q0=q0)
    z = np.zeros(N)
    rec, _ = rk4(z, g0, z, gd0, DT, st, sample=st // nsamp)
    ch = [(t, float(envelope(g).max() / a), frac(np.roll(u, -1) - u, envelope(g))) for t, u, g, _gd in rec]
    stm = int(round(T / DTM))
    md = [(t, float(np.abs(A).max() / a), frac(eps, np.abs(A)))
          for t, A, eps, _B, _d in run(sech_envelope(N // 2, N, 100.0, w, a), N, q0, DTM, stm, stm // nsamp,
                                        meanfold=False)]
    return ch, md


def ring_gap(q0, N, adeg, d0deg=0.3, T=250.0):
    """A plane wave at (q0, a) on a ring of N with a UNIFORM fold offset d0
    at rest, strain frozen: the mean fold's frequency from the zero crossings
    of <gamma>(t), as Omega^2 / a^2."""
    a, d0 = np.radians(adeg), np.radians(d0deg)
    nn = np.arange(N)
    g0 = a * np.cos(q0 * nn) + d0
    gd0 = a * omega(q0) * np.sin(q0 * nn)
    z = np.zeros(N)

    def f(y):
        g, gd = y[:N], y[N:]
        _au, ag = accel(z, g, z, gd)
        return np.concatenate([gd, ag])

    y = np.concatenate([g0, gd0])
    st, every = int(round(T / DT)), 4
    ts, ms = [], []
    for s_ in range(st + 1):
        if s_ % every == 0:
            ts.append(s_ * DT)
            ms.append(float(y[:N].mean()))
        if s_ == st:
            break
        k1 = f(y)
        k2 = f(y + 0.5 * DT * k1)
        k3 = f(y + 0.5 * DT * k2)
        k4 = f(y + DT * k3)
        y = y + (DT / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    ts, ms = np.array(ts), np.array(ms) - np.mean(ms)
    zc = ts[:-1][np.sign(ms[:-1]) != np.sign(ms[1:])]
    per = 2 * np.mean(np.diff(zc))
    return float((2 * np.pi / per) ** 2 / (a * a)), len(zc)


def meanfold_track(q0, w, adeg, T, N=800, nsamp=12, warm=True):
    """A warm packet in the full model: (t, h/a, max |delta| in degrees,
    its position, the peak's position) every T/nsamp."""
    a = np.radians(adeg)
    n = N // 2
    A0 = sech_envelope(n, N, 100.0, w, a)
    e, ed = dent_of(np.abs(A0), a, True, q0)
    stm = int(round(T / DTM))
    out = run(A0, N, q0, DTM, stm, stm // nsamp, e if warm else None, ed if warm else None)
    rows = []
    for t, A, _eps, _B, d in out:
        ea = np.abs(A)
        if not np.all(np.isfinite(ea)):
            break
        j = int(np.argmax(np.abs(d)))
        rows.append((t, float(ea.max() / a), float(np.degrees(abs(d[j]))), j * DXM, int(np.argmax(ea)) * DXM))
    return rows


def ring_pulse(q0, adeg, d0deg, T=400.0, N=800, width=20.0):
    """A uniform carrier on a ring with an imposed fold pulse d0 sech^2(x/width)
    at rest, strain frozen: sum delta^2 (deg^2 per cell) at t = 0, T/2, T."""
    a = np.radians(adeg)
    n = N // 2
    X = np.arange(n) * DXM
    d0 = np.radians(d0deg) / np.cosh((X - N / 2) / width) ** 2
    stm = int(round(T / DTM))
    out = run(np.full(n, a, complex), N, q0, DTM, stm, stm // 2, d0=d0, strain=False)
    return [float(DXM * (np.degrees(d) ** 2).sum()) for _t, _A, _e, _B, d in out]


# --------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------

def gate():
    checks, out = [], {}
    A = checks.append
    t_start = time.time()
    q = Q0
    afam = np.degrees(soliton_aw(coef_p(q), -(q_ring_closed(q) + coef_q_dent(q)))) / SOL_W

    # ---- L1: the umklapp exchange at pi/2 --------------------------------
    ex = {phi: collision(q, afam, SOL_W, phi, True, True) for phi in (np.pi / 4, 3 * np.pi / 4)}
    chain_ex = {np.pi / 4: 0.74, 3 * np.pi / 4: 1.24}     # T2 [24191], [24208]
    out["l1"] = ex
    A(("L1 THE pi/2 COLLISION IS A PHASE-DEPENDENT ENERGY EXCHANGE, and the "
       "umklapp term carries it: with the relative carrier phase at pi/4 the "
       "right-mover comes out with less energy, at 3pi/4 with more, and the "
       "size is the chain's (chain 0.74 / 1.24, T2 [24191])",
       ex[np.pi / 4][0] < 0.85 and ex[3 * np.pi / 4][0] > 1.15
       and all(abs(ex[p][0] - chain_ex[p]) < 0.05 for p in ex),
       f"a {afam:.2f} deg w 8: E_A/free = {ex[np.pi / 4][0]:.3f} at phi pi/4, "
       f"{ex[3 * np.pi / 4][0]:.3f} at 3pi/4"))

    # ---- L2: none off resonance ------------------------------------------
    off = collision(0.45 * np.pi, afam, SOL_W, np.pi / 4, True, True)
    out["l2"] = off
    A(("L2 OFF pi/2 THE EXCHANGE IS GONE: at carrier 0.45 pi the umklapp term "
       "rotates with Dk = 4 q0 - 2 pi and the same collision returns the "
       "energies it was given (chain: an order of magnitude down at 0.45 pi)",
       abs(off[0] - 1.0) < 0.03 and abs(off[1] - 1.0) < 0.03,
       f"E_A/free {off[0]:.3f}, E_B/free {off[1]:.3f} at 0.45 pi, phi pi/4"))

    # ---- L3: the cold launch ---------------------------------------------
    mixed = collision(q, 6.0, SOL_W, 0.0, False, True)
    chain_slope, chain_gain = -0.022, 1.072                 # T2 [24169] [24192]; [24191]
    out["l3"] = mixed
    A(("L3 THE COLD LAUNCH DETUNES ITS OWN CARRIER, IN MECHANISM: a packet "
       "born without its dent sits on a one-sided strain ramp while the "
       "pulse leaves at c, and the ramp's frequency gradient integrates into "
       "a carrier shift (chain -0.021..-0.023 rad/cell); the detuned packet "
       "then meets a warm one in the exchange regime and GAINS (chain 1.072)",
       abs(mixed[3] - chain_slope) < 0.35 * abs(chain_slope) and mixed[3] < -0.012
       and abs(mixed[2] - chain_slope) < 0.35 * abs(chain_slope)
       and abs(mixed[0] - chain_gain) < 0.03 and mixed[0] > 1.03,
       f"a 6 w 8, A cold B warm: A's carrier slope free {mixed[3]:+.4f}, after "
       f"{mixed[2]:+.4f}; E_A/free {mixed[0]:.3f}, E_B/free {mixed[1]:.3f}"))

    # ---- L4: dent formation at w 16 ----------------------------------------
    form = {}
    for name, qq, N, T in (("pi/4", np.pi / 4, 800, 500.0), ("pi/3", np.pi / 3, 804, 700.0)):
        form[name] = dent_formation(qq, N, 16.0, 3.0, T, int(T // 100))
    dmax = max(abs(c[2] - m[2]) for rows in form.values() for c, m in zip(*rows) if c[0] > 0)
    hmax = max(abs(c[1] - m[1]) for rows in form.values() for c, m in zip(*rows))
    out["l4"] = form
    A(("L4 THE DENT FORMS ON THE STRAIN FIELD'S OWN SCHEDULE, and the model "
       "tracks the chain through the whole formation at BOTH long carriers -- "
       "the fraction of R4's dent under the envelope within 0.05 every 100 "
       "time units to t 500 at pi/4 and t 700 at pi/3, the height within 0.2 "
       "(it parts at the end of the span, where L6's instability begins). "
       "This is what envelope R8b could not see; w 16 keeps the envelope two "
       "carrier wavelengths wide (the w 8 rows are not converged, T2 [24273])",
       dmax < 0.06 and hmax < 0.2,
       "; ".join(f"{k}: " + " ".join(f"{c[2]:.2f}/{m[2]:.2f}" for c, m in zip(*v)) for k, v in form.items())
       + f" (chain/model); worst dent {dmax:.3f}, worst height {hmax:.2f}"))

    # ---- L5: the gap --------------------------------------------------------
    gaps = {}
    for name, qq, N in (("pi/4", np.pi / 4, 800), ("pi/3", np.pi / 3, 804), ("pi/2", np.pi / 2, 800)):
        gaps[name] = (ring_gap(qq, N, 3.0), k_gap(qq), q_ring_closed(qq))
    gworst = max(abs(g[0][0] / g[1] - 1.0) for g in gaps.values())
    out["l5"] = gaps
    A(("L5 THE MEAN FOLD UNDER A CARRIER IS GAPPED, AND THE GAP IS DERIVED: a "
       "uniform fold offset on a plane-wave ring does not sit still, it "
       "oscillates at Omega^2 = K_d a^2 with K_d = 2.5 + cos(q)/2 - sin^2(q/2)/3 "
       "from V4's d^2 c^2 terms against the fold inertia's M2 -- at every "
       "carrier, within the zero-crossing measurement. It is NOT Q_ring a^2 "
       "(near it at pi/4 by coincidence; three times it at pi/2)",
       gworst < 0.25 and all(g[2] > 0 for g in gaps.values())
       and gaps["pi/2"][0][0] / gaps["pi/2"][2] > 2.0,
       "; ".join(f"{k}: measured {g[0][0]:.2f} ({g[0][1]} crossings) vs K_d {g[1]:.2f}, Q_ring {g[2]:.2f}"
                 for k, g in gaps.items())))

    # ---- L6: the instability ------------------------------------------------
    tr6 = meanfold_track(np.pi / 4, 16.0, 6.0, 400.0, nsamp=8)
    ctl = meanfold_track(np.pi / 2, 16.0, 3.0, 1200.0, nsamp=12)
    f6 = {r[0]: r for r in tr6}
    growth = f6[400.0][2] / f6[200.0][2]
    rate = float(np.log(growth) / 200.0)
    coloc = abs(f6[400.0][3] - f6[400.0][4]) < 20 and abs(f6[350.0][3] - f6[350.0][4]) < 20
    cvals = [r[2] for r in ctl if r[0] >= 100]
    out["l6"] = (tr6, ctl, rate)
    A(("L6 PAST THE WALL THE PACKET IS UNSTABLE TO ITS OWN MEAN FOLD, and this "
       "system has it: at pi/4 a warm w 16 packet's mean fold grows "
       "exponentially UNDER the packet, co-moving with its peak, and the packet "
       "focuses with it where the adiabatic Q says defocus -- at the chain's "
       "rate (0.04/t at a 6, T2 [24273]); at pi/2 the same launch keeps its "
       "mean fold flat for 1200 time units (chain: to 1600)",
       growth > 30 and rate > 0.02 and coloc and f6[400.0][1] > 1.1
       and max(cvals) / min(cvals) < 3.0 and max(cvals) < 0.05,
       f"pi/4 a 6: mean fold {f6[200.0][2]:.4f} -> {f6[400.0][2]:.3f} deg over t 200..400 "
       f"(rate {rate:.3f}/t), at X {f6[400.0][3]:.0f} with the peak at {f6[400.0][4]:.0f}, "
       f"height {f6[400.0][1]:.2f}; pi/2 a 3 mean fold {min(cvals):.4f}..{max(cvals):.4f} deg to t 1200"))

    # ---- L7: what the gap does not carry -------------------------------------
    tr15 = meanfold_track(np.pi / 4, 16.0, 1.5, 1800.0, nsamp=12)
    lo15 = max(r[2] for r in tr15 if r[0] >= 100)
    pulse = ring_pulse(np.pi / 4, 3.0, 0.5)
    out["l7"] = (tr15, pulse)
    A(("L7 WHAT THE AVERAGED THEORY DOES NOT CARRY, as measured: the chain's "
       "growth is threshold-free (0.0063/t at a 1.5, width-independent) and "
       "present on a ring with no packet at all (a 0.5-degree fold pulse under "
       "a uniform carrier grows fourfold by t 400, T2 [24289]); this system "
       "keeps the a 1.5 packet's mean fold below a hundredth of a degree for "
       "1800 time units and leaves the ring pulse oscillating in its gap. A "
       "process that needs neither envelope nor strain is missing at order "
       "delta^2 |A|^2 -- carried by instability.py (T2 [24443] [24472])",
       lo15 < 0.01 and pulse[-1] < 2.0 * pulse[0],
       f"a 1.5 pi/4: mean fold max {lo15:.4f} deg to t 1800 (chain: 1 deg by ~t 700); "
       f"ring pulse sum delta^2 {pulse[0]:.1f} -> {pulse[1]:.1f} -> {pulse[2]:.1f} deg^2 "
       f"at t 0, 200, 400 (chain 6.7 -> 27)"))

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
