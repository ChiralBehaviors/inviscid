"""inertial_array -- the medium, running. A time-domain honeycomb, and what it transports.

*** VOID -- OWNER DECISION 2026-08-27. ***
The strut-compliance fork is REJECTED. This file models the jitterbug honeycomb
as an ELASTIC NETWORK -- bars carrying springs of unit stiffness -- and the
model of record is a RIGID-STRUT LINKAGE. Rigid struts are what make a jitterbug
a jitterbug; putting a spring on every strut does not add compliance to the
linkage, it replaces the linkage. Every number here is a property of the spring
network, not of the medium.

The fork was never declared and never authorised. T2 inviscid 23533 section 2,
written 2026-08-26, admits it: "Every gapless result in this epic is downstream
of the first choice, which had never been declared."

Kept, not deleted, because this repo's README is explicit that these scripts are
the evidence base for its retractions and "a retraction whose evidence has
evaporated is just folklore."

Read T2 inviscid [23562] "OWNER DECISION 2026-08-27 -- the strut-compliance fork
is REJECTED" before using anything below. The live problem is finding a
potential energy V for the RIGID linkage: this repo's README states "there is no
potential energy anywhere in the model, so its six DOF are six zero-frequency
modes and no dispersion relation can exist yet." Resume at the jb_h..jb_q V
survey and bead inviscid-qvf.21 (V=0, contact).


Everything before this file is STATIC: packing, exchange, dispersion relations,
one order of anharmonic coupling. All of it linearised about a fixed
configuration. Nothing took an initial condition and evolved it, and a medium
you cannot evolve is a dispersion relation with ambitions.

THE MODEL, and every choice in it is already declared elsewhere rather than
invented here:
    STRUCTURE   the rectified cubic honeycomb as a bar-and-joint framework --
                nodes at the shared vertices, bars along the rigid triangles'
                edges (jb_hc). Finite patch, free boundary, no periodicity.
    COMPLIANCE  struts compliant, hinges free. That is the fork jb_fl's C rows
                declare, and it is why the phase field is gapless rather than
                massive. Bars are springs of unit stiffness; nodes carry unit
                mass. Neither number is physical -- omega = sqrt(K/M)*sigma, so
                every speed here is a pure number, exactly as everywhere else.
    INTEGRATOR  velocity Verlet. Symplectic, so the energy audit is a real
                check rather than a hopeful one.

WHY THE EXPERIMENT IS GATED AND NOT JUST THE PHYSICS. Four scratch prototypes
preceded this file and all four produced a transport speed that could not be
defended -- a threshold front biased by its own threshold, a cross-correlation
picking up boundary reflections, a packet centroid tracking a superposition that
had split across branches, and an eigenmode launch in which 170 of 584 nodes
silently received no initial condition at all. Not one of those was a physics
error. Each was an experiment that had not been asked whether it was valid, and
each would have died instantly against a row that could fail:
    E2  the initial condition reaches EVERY node          (kills the fourth)
    E3  the launched packet stays on ONE branch           (kills the third)
    E4  the measurement window closes before the first
        reflection returns                               (kills the second)
    E5  the front is tracked by a statistic that does not
        depend on an arbitrary threshold                 (kills the first)
The house style exists for exactly this and it was not applied, because the
experiment was being treated as scaffolding rather than as the deliverable. It
is the deliverable.

WHAT IT ESTABLISHES. One thing, and it is spectral rather than dynamical.

THE PHASE DEGREE OF FREEDOM IS NOT A NORMAL MODE. Its spectral weight is
dominant on one band but not concentrated there: 0.588 and 0.644 at the two
wavevectors probed, consistent with jb_fl's independent measurement that the
maximum weight on any single band is about 0.5. So there is no pure phase wave
to launch; a phase disturbance is a superposition by construction. This is
computed entirely in the PERIODIC UNIT CELL, never touches the patch, and
stands.

WHAT IT NO LONGER CLAIMS, AND WHY -- retracted 2026-08-26. This file previously
concluded that a phase disturbance DISPERSES in the time domain, at 0.31
against the branch's own 0.54 at k = 0.20 G and not advancing at all at
k = 0.35 G. Both numbers were real and neither meant what was claimed. Two
independent defects, either of which is sufficient on its own:

    THE INITIAL CONDITION WAS NOT THE EIGENMODE. The loop discarded the slot's
    Bravais translation nn. A corner of this honeycomb is shared between 2 and
    8 cells, so each node was written a mean of 3.66 times with a DIFFERENT
    unit cell's amplitude each time -- disagreeing by up to 0.252 against a
    mode amplitude of 0.408 -- and the last write won by the order of
    med.cells. E2 passed throughout, because E2 counts how many nodes the loop
    REACHED and a node written eight times is reached eight times. Coverage was
    gated; single-valuedness was assumed and enforced nowhere. Now gated as
    E2c, which is free: the same redundancy that corrupted the field proves it
    correct once the writes agree, and they now agree to 8e-16. Fixing this
    also dropped the energy drift by a factor of ten, from 2e-05 to 3e-06,
    which is what a genuine eigenmode should do.

    THE LAUNCH IS A STANDING WAVE. A velocity-only kick of Re(ev e^{ikR})
    excites +k and -k equally, and the centroid of two counter-propagating
    packets sits still no matter what the medium transports. Measured
    directionality 0.03 and 0.12 on a scale where 1 is one-way, using a
    statistic calibrated against synthetic packets rather than argued from a
    sign convention. E6 now states this, and CAN FAIL: make the launch one-way
    and it goes red.

AND THE PATCH CANNOT MEASURE A GROUP VELOCITY REGARDLESS. At k = 0.20 G one
wavelength is 23.1 against a 39.3 span and a 4.7 envelope. A packet narrower
than its own wavelength has no meaningful centre velocity, so no repair to the
initial condition would have rescued the measurement. Transport needs a patch
several wavelengths long, or a method that does not use a packet at all.

The linearised medium and the running medium are therefore NOT in conflict, and
never were shown to be. jb_fl's F10 found the phase spectral peak is a
continuous branch with group velocity 0.534, and that stands. What a phase
DISTURBANCE does in time is a separate question, and this file does not answer
it -- which is a smaller claim than the one it used to make, and the only one
its measurements support.
"""

from __future__ import annotations

import itertools as it
import sys

import numpy as np

from analysis.retired.strut_springs import honeycomb_waves as HC
from analysis.model import plates as Z
A_REF = -30.0
PHASE_OFFSET = 60.0

#: Patch extent, in integer honeycomb lattice steps. Long in x so a packet has
#: room to travel; thin across so the node count stays tractable.
PATCH = (34, 2, 2)

DT = 0.02                 # Verlet step
KICK = 0.02               # launch amplitude, small enough to stay linear
ENERGY_TOL = 1e-4         # relative drift, measured 7e-06 .. 2e-05
PURITY_BAND = (0.30, 0.80)   # the phase DOF is DOMINANT on one band but not
                             # PURE. Measured 0.588 and 0.644, matching jb_fl's
                             # independent max-weight-per-band of about 0.5.
COVERAGE_TOL = 0          # nodes allowed to miss the initial condition: none
SPEED_TOL = 0.12          # |measured - predicted| / predicted
PROBE_K = (0.20, 0.35)    # wavevector fractions of G to launch at


def build(nx, ny, nz, a=A_REF):
    """A finite honeycomb patch as a bar-and-joint framework.

    Free boundary and no periodic reduction: this is a real object, not a unit
    cell. Cells sit at all-even (VE) and all-odd (hole) integer sites, the hole
    running PHASE_OFFSET ahead as the exchange requires."""
    L = HC.lattice(a)
    nodes, reps, cells, bars = {}, [], [], set()

    def nid(p):
        k = tuple(np.round(p, 6))
        if k not in nodes:
            nodes[k] = len(reps)
            reps.append(np.array(k))
        return nodes[k]

    sites = [c for c in it.product(range(nx + 1), range(ny + 1), range(nz + 1))
             if all(x % 2 == 0 for x in c) or all(x % 2 != 0 for x in c)]
    for s in sites:
        even = sum(abs(x) for x in s) % 2 == 0
        ph = a if even else a + PHASE_OFFSET
        origin = L * np.array(s, dtype=float)
        X = Z.corners(ph) + origin
        h = 1e-5
        D = ((Z.corners(ph + h) - Z.corners(ph - h)) / (2 * h))
        ids = {}
        for f in range(8):
            fid = [nid(X[f][c]) for c in range(3)]
            for c in range(3):
                ids[fid[c]] = D[f][c]
            for u, v in it.combinations(fid, 2):
                bars.add((min(u, v), max(u, v)))
        cells.append(dict(site=s, even=even, ids=ids, origin=origin))
    return np.array(reps), sorted(bars), cells, L


class Medium:
    """The patch, its springs, and the phase direction of every cell."""

    def __init__(self, patch=PATCH, a=A_REF):
        self.P, B, self.cells, self.L = build(*patch, a=a)
        self.n = len(self.P)
        self.bi = np.array([b[0] for b in B])
        self.bj = np.array([b[1] for b in B])
        self.L0 = np.linalg.norm(self.P[self.bi] - self.P[self.bj], axis=1)
        self.nbars = len(B)
        pd = []
        for c in self.cells:
            v = np.zeros((self.n, 3))
            for k, d in c["ids"].items():
                v[k] = d
            pd.append(v / np.linalg.norm(v))
        self.PD = np.array(pd)
        self.cx = np.array([self.L * c["site"][0] for c in self.cells])
        self._index = {tuple(np.round(p, 6)): i for i, p in enumerate(self.P)}

    def index_of(self, p):
        return self._index.get(tuple(np.round(p, 6)))

    def forces(self, x):
        d = x[self.bi] - x[self.bj]
        ln = np.linalg.norm(d, axis=1)
        e = ln - self.L0
        f = (e / ln)[:, None] * d
        F = np.zeros_like(x)
        np.add.at(F, self.bi, -f)
        np.add.at(F, self.bj, f)
        return F, 0.5 * float(np.sum(e * e))

    def phase_field(self, vel):
        return np.einsum("cnd,nd->c", self.PD, vel)

    def fastest_speed(self):
        """Largest signal speed in the medium, from the stiffest bar mode. Sets
        how long the measurement window can be before a reflection returns."""
        return float(np.sqrt(2.0 * self.nbars / self.n))


def launch(med, kfrac, direction=(1.0, 0.0, 0.0)):
    """A narrow-band packet on the PHASE branch's own Bloch eigenmode.

    Not the phase DIRECTION -- that is not an eigenvector, so a packet built
    from it splits across branches within a few steps and its centroid stops
    meaning anything (prototype three). The eigenvector is taken from the
    periodic unit cell and mapped onto every patch node by reducing the node's
    position modulo the cubic lattice; E2 gates that the mapping reaches all of
    them, because in prototype four it silently reached 29%."""
    fw = HC.h4_framework(A_REF)
    Pu, bars_u, Au, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
    G = np.pi / Au
    u = np.array(direction, dtype=float)
    u = u / np.linalg.norm(u)
    kv = kfrac * G * u

    M = HC.bloch(Pu, bars_u, Au, kv)
    w, V = np.linalg.eigh(M.conj().T @ M)
    w = np.sqrt(np.maximum(w, 0.0))
    B = HC.phase_basis(Pu, slots, Au, kv)
    up = B[:, 0] + B[:, 1]
    up = up / np.linalg.norm(up)
    m = int(np.argmax(np.abs(V.conj().T @ up) ** 2))
    ev = V[:, m]

    # Map patch node -> unit-cell node BY CONSTRUCTION rather than by reducing
    # coordinates. Both the patch and the unit cell place a corner as
    # (cell parity, face, corner), so that triple IS the identity of the node
    # and no arithmetic is needed. Two earlier versions reduced modulo the
    # lattice and matched by rounded key or nearest neighbour; they dropped 414
    # and then 48 and then 80 nodes respectively, because a node on a cell
    # boundary reduces to either side on the last bit. Geometry was the wrong
    # tool for a bookkeeping problem.
    # AND THE SLOT'S OWN BRAVAIS TRANSLATION nn, which the first version of
    # this loop unpacked as _nn and threw away. That is not a detail. A corner
    # of this honeycomb is SHARED between 2 and 8 cells, so the loop reaches
    # each physical node several times -- mean 3.66, max 8 -- and every visit
    # wrote the amplitude of a DIFFERENT unit cell to the same node. The writes
    # disagreed by up to 0.252 against a mode amplitude of 0.408, i.e. 62%, and
    # the last one won by iteration order. The launched field was therefore
    # never the eigenmode; it was a scramble selected by the order of
    # med.cells. With nn included, every visit computes the same unit cell
    # R_cell + nn and the writes agree to 4e-16.
    #
    # The redundancy that was corrupting the initial condition IS the check
    # that it is now right, so it is measured here and gated as E2c. Nothing
    # gated it before: E2 counts how many nodes the loop REACHED, and a node
    # written eight times with eight different values is reached eight times.
    # Coverage was enforced; single-valuedness was assumed and enforced
    # nowhere.
    covered = 0
    worst = 0.0
    x0 = med.cx.max() * 0.18
    wd = med.cx.max() * 0.12
    vel = np.zeros((med.n, 3))
    Lc = med.L
    seen = {}
    for cell in med.cells:
        ci = 0 if cell["even"] else 1
        base = cell["origin"] - (0.0 if ci == 0 else Lc) * np.ones(3)
        ph_ref = A_REF if ci == 0 else A_REF + PHASE_OFFSET
        X = Z.corners(ph_ref) + cell["origin"]
        for (sci, _ph, _off, f, c, j, nn) in slots:
            if sci != ci:
                continue
            p = X[f][c]
            node = med.index_of(p)
            if node is None:
                continue
            R = base + Au * np.asarray(nn, dtype=float)
            z = ev[3 * j:3 * j + 3] * np.exp(1j * float(kv @ R))
            if node in seen:
                worst = max(worst, float(np.abs(z - seen[node]).max()))
            seen[node] = z
            env = np.exp(-((p[0] - x0) / wd) ** 2)
            vel[node] = KICK * env * np.real(z)
    covered = len(seen)
    # A node legitimately gets zero velocity where the envelope has decayed, so
    # counting nonzero rows understates the coverage. The fix for that was
    # `np.any(np.isfinite(vel))` -- which is True for every row of an array
    # initialised to zeros, so it returned med.n unconditionally and E2 could
    # not fail. Count the nodes the loop actually ASSIGNED, which is what
    # `seen` holds.
    reached = len(seen)
    purity = float(np.abs(np.vdot(ev, up)) ** 2)
    return vel, reached, float(w[m]), G, kv, purity, worst


def group_velocity(kfrac, direction=(1.0, 0.0, 0.0), h=1e-4):
    """d(omega)/dk of the phase-carrying branch, from the dispersion."""
    fw = HC.h4_framework(A_REF)
    Pu, bars_u, Au, slots = fw["P"], fw["bars"], fw["A"], fw["slots"]
    G = np.pi / Au
    u = np.array(direction, dtype=float)
    u = u / np.linalg.norm(u)

    def om(t):
        kv = t * G * u
        M = HC.bloch(Pu, bars_u, Au, kv)
        w, V = np.linalg.eigh(M.conj().T @ M)
        w = np.sqrt(np.maximum(w, 0.0))
        B = HC.phase_basis(Pu, slots, Au, kv)
        up = B[:, 0] + B[:, 1]
        up = up / np.linalg.norm(up)
        return w[int(np.argmax(np.abs(V.conj().T @ up) ** 2))]
    return (om(kfrac + h) - om(kfrac - h)) / (2 * h * G)


def run(med, vel0, steps, sample=25):
    """Velocity Verlet. Returns (times, phase-field history, energy history)."""
    x = med.P.copy()
    vel = vel0.copy()
    F, _E = med.forces(x)
    T, U, E = [], [], []
    for n in range(steps):
        vel += 0.5 * DT * F
        x += DT * vel
        F, pot = med.forces(x)
        vel += 0.5 * DT * F
        if n % sample == 0:
            T.append(n * DT)
            U.append(med.phase_field(vel))
            E.append(pot + 0.5 * float(np.sum(vel * vel)))
    return np.array(T), np.array(U), np.array(E)


def directionality(T, U, cx, nx=128):
    """Signed one-way-ness of the phase field in [-1, 1]. +1 right, 0 standing.

    Needed because a velocity-only kick of Re(ev e^{ikR}) excites +k and -k
    EQUALLY -- it is a standing wave, and a standing wave's centroid does not
    advance no matter what the medium transports. Nothing in this file
    measured that, so E6 read a launch artefact as a property of the medium.

    NOT the spatial spectrum of a single frame: for a REAL field |FFT|^2 is
    symmetric in k by construction, so a +k/-k comparison there is identically
    zero. It returned exactly 0.000 on four different runs, which is how that
    was caught. The (k, omega) spectrum obeys P(k,w) = P(-k,-w) instead, so
    right-going content sits at (k>0, w>0) and left-going at (k<0, w>0), and
    THOSE two quadrants do differ for a real field.

    SIGN CALIBRATED, not argued. np.fft uses e^{-i...}, so a physical e^{+ikx}
    lands at a negative FFT index. Rather than reason about that, the statistic
    was run on synthetic packets with known answers: right-going gave -0.950,
    left-going +0.945, standing -0.062. Hence the flip below."""
    order = np.argsort(cx)
    xs = np.linspace(cx.min(), cx.max(), nx)
    G = np.array([np.interp(xs, cx[order], fr[order]) for fr in U])
    G = G - G.mean()
    P = np.abs(np.fft.fft2(G)) ** 2
    nt, nxx = P.shape
    ht, hx = nt // 2, nxx // 2
    rgt = P[1:ht, 1:hx].sum()
    lft = P[1:ht, hx + 1:].sum()
    return float((lft - rgt) / (rgt + lft + 1e-30))


def experiment(med, kfrac):
    vel0, covered, omega, G, kv, purity, mapdist = launch(med, kfrac)
    span = med.cx.max()
    t_reflect = 2.0 * span / med.fastest_speed()
    steps = int(0.85 * t_reflect / DT)
    T, U, E = run(med, vel0, steps)

    # E3's purity comes from `launch`: |<eigenvector, phase direction>|^2 in the
    # PERIODIC unit cell, which is a genuine fraction bounded by 1. The first
    # version summed projections of the patch launch onto the per-cell phase
    # directions, which are not mutually orthogonal -- it returned 1.41 and
    # 1.55, i.e. "141% of the energy is on one branch", and passed its own
    # threshold while measuring nothing. A row whose statistic can exceed its
    # own maximum is not a row.

    # E5: threshold-free front statistic -- the energy-weighted centroid of the
    # phase field, which needs no cutoff and no arbitrary arrival criterion
    P2 = U ** 2
    tot = P2.sum(axis=1)
    good = tot > 0
    cen = (P2[good] * med.cx).sum(axis=1) / tot[good]
    tg = T[good]
    m = (tg > T[-1] * 0.25) & (tg < T[-1] * 0.75) & (cen < span * 0.75)
    speed = float(np.polyfit(tg[m], cen[m], 1)[0]) if m.sum() > 4 else float("nan")
    vg = group_velocity(kfrac)
    return dict(k=kfrac, covered=covered, n=med.n, omega=omega, purity=purity,
                mapdist=mapdist, oneway=directionality(T, U, med.cx),
                G=G, span=span, wd=span * 0.12,
                drift=abs(E[-1] - E[0]) / E[0], steps=steps,
                t_window=T[-1], t_reflect=t_reflect,
                speed=speed, vg=vg,
                err=abs(speed - vg) / abs(vg) if np.isfinite(speed) else float("inf"),
                monotone=bool(np.all(np.diff(cen[m]) > -span * 0.02)))


def gate(med, runs):
    checks = []
    R = checks.append

    R(("E1  the patch is a real honeycomb: every bar is a strut, sqrt(2)",
       med.nbars > 0
       and abs(med.L0.max() - np.sqrt(2)) < 1e-5
       and abs(med.L0.min() - np.sqrt(2)) < 1e-5,
       f"{med.nbars} bars, {med.L0.min():.6f}..{med.L0.max():.6f}",
       f"{np.sqrt(2):.6f}"))
    R(("E2  the initial condition reaches EVERY node -- prototype four "
       "silently reached 29% of them and its answer meant nothing",
       all(r["covered"] == r["n"] for r in runs),
       f"{min(r['covered'] for r in runs)}/{med.n} nodes", "all of them"))
    R(("E2c EVERY SHARED NODE GETS THE SAME VALUE from all 2-8 cells that "
       "reach it. A corner is shared, so the loop writes each node a mean of "
       "3.66 times; before nn was included those writes disagreed by 0.252 "
       "against a mode amplitude of 0.408 and the last won by iteration "
       "order. The redundancy that was corrupting the IC is the check that "
       "it is right",
       all(r["mapdist"] < 1e-9 for r in runs),
       f"worst disagreement between writes "
       f"{max(r['mapdist'] for r in runs):.2e}", "< 1e-9"))
    R(("E3  THE PHASE DEGREE OF FREEDOM IS NOT A NORMAL MODE: dominant on one "
       "band, not concentrated there, so there is no pure phase wave to launch "
       "-- TWO-SIDED, and both edges would be findings",
       all(PURITY_BAND[0] < r["purity"] < PURITY_BAND[1] for r in runs),
       f"purity {[round(r['purity'], 4) for r in runs]}", f"in {PURITY_BAND}"))
    R(("E4  the measurement window CLOSES before the first boundary reflection "
       "returns", all(r["t_window"] < r["t_reflect"] for r in runs),
       f"window {max(r['t_window'] for r in runs):.1f} vs reflection at "
       f"{min(r['t_reflect'] for r in runs):.1f}", "window < reflection"))
    R(("E5  ENERGY is conserved -- the integrator is symplectic, so this is a "
       "real check and not a hopeful one",
       all(r["drift"] < ENERGY_TOL for r in runs),
       f"worst relative drift {max(r['drift'] for r in runs):.2e}",
       f"< {ENERGY_TOL:.0e}"))
    R(("E6  THIS LAUNCH IS A STANDING WAVE, so its centroid speed is NOT a "
       "transport measurement. A velocity-only kick of Re(ev e^{ikR}) "
       "excites +k and -k equally; the centroid of two counter-propagating "
       "packets sits still whatever the medium does. CAN FAIL -- if the "
       "launch were one-way this row goes red and E7 becomes readable",
       all(abs(r["oneway"]) < 0.25 for r in runs),
       f"directionality {[round(r['oneway'], 3) for r in runs]} "
       f"(+1 right, 0 standing)", "|dir| < 0.25, i.e. standing"))
    R(("E7  AND THE PATCH CANNOT MEASURE A GROUP VELOCITY ANYWAY: the packet "
       "is narrower than its own wavelength, so it has no meaningful centre "
       "velocity. At k = 0.20 G one wavelength is 23.1 against a 39.3 span "
       "and a 4.7 envelope. PRINTED NOT GATED -- the centroid speeds below "
       "are recorded, and they are not transport speeds",
       True,
       f"lambda {[round(2 * np.pi / (r['k'] * r['G']), 1) for r in runs]} vs "
       f"envelope {runs[0]['wd']:.1f} and span {runs[0]['span']:.1f}; "
       f"centroid {[round(r['speed'], 4) for r in runs]} vs branch "
       f"{runs[0]['vg']:.4f}", "recorded, not gated"))
    R(("E8  PRINTED NOT GATED: every speed here is a pure number. omega = "
       "sqrt(K/M)*sigma and this file sets K = M = 1, exactly as jb_fl's C "
       "rows declare", True,
       f"measured {runs[0]['speed']:.6f} x sqrt(K/M)", "printed"))

    print()
    print("=" * 78)
    print(f"GATE  {len(checks)} rows")
    print("=" * 78)
    for name, ok, val, crit in checks:
        print(f"  {'PASS' if ok else 'FAIL':4s}  {name:66s} {str(val):>26s} {str(crit):>18s}")

    print()
    print("  ROWS THAT EXIST ONLY TO STOP ANOTHER ROW BEING UNFALSIFIABLE:")
    print("   * ALL OF E2 THROUGH E6. They are not hygiene, they are the four")
    print("     scratch prototypes that preceded this file, each of which")
    print("     produced a transport speed that could not be defended. E2 is")
    print("     the one where 170 of 584 nodes got no initial condition. E3 is")
    print("     the one where the packet split across branches. E4 is the one")
    print("     where the cross-correlation was reading reflections. E6 is the")
    print("     one where the centroid doubled back and the fit was fitting")
    print("     noise. Every one of them would have died in a line against")
    print("     its own row, and none of them was a physics error -- they were")
    print("     experiments nobody had asked whether they were valid.")
    print("   * E7 without E2-E6 is worthless: a speed that agrees with the")
    print("     prediction while the initial condition covers a third of the")
    print("     patch is a coincidence, and would have been reported as a")
    print("     result.")
    print()
    print("  A ROW DELIBERATELY NOT BUILT: anything about the <110> anharmonic")
    print("  sink. jb_fl measured a cubic vertex forcing a zero-restoring-force")
    print("  channel along the closest-packing directions, and the obvious next")
    print("  move is to launch along <110> here and watch for it. This file")
    print("  does not, because the amplitude is deliberately small enough to")
    print("  stay linear and a linear run cannot show a second-order effect.")
    print("  Seeing it needs finite amplitude and a longer window, and that is")
    print("  its own experiment with its own validity rows.")
    print()
    print("  WHAT THIS FILE DOES NOT MODEL: contact. Cells can pass through one")
    print("  another here. The medium's valid range is the closed interval")
    print("  [0, -60] and separation goes positive outside it (jb_hc), so a run")
    print("  driven far enough would leave the range without being stopped.")
    print("  Small-amplitude motion about a = -30 does not approach either end,")
    print("  and no row here claims otherwise.")

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
    print("jb_aa -- the medium, running")
    print("=" * 78)
    print("  A finite honeycomb patch as a bar-and-joint framework, struts")
    print("  compliant and hinges free per jb_fl's declared fork, integrated")
    print("  with velocity Verlet. A narrow-band packet is launched on the")
    print("  phase branch's own Bloch eigenmode and tracked by a")
    print("  threshold-free centroid. The question is whether the medium run")
    print("  FORWARD IN TIME transports at the speed the medium LINEARISED")
    print("  says it should.")
    med = Medium()
    print()
    print(f"  patch: {len(med.cells)} cells, {med.n} nodes, {med.nbars} bars, "
          f"{3 * med.n} DOF, length {med.cx.max():.2f}")
    runs = [experiment(med, k) for k in PROBE_K]
    print()
    print("-" * 78)
    print("  TRANSPORT")
    print("-" * 78)
    print(f"    {'k/G':>6s} {'omega':>10s} {'measured':>10s} {'predicted':>10s} "
          f"{'err':>8s} {'purity':>8s} {'drift':>9s}")
    for r in runs:
        print(f"    {r['k']:6.2f} {r['omega']:10.6f} {r['speed']:10.6f} "
              f"{r['vg']:10.6f} {r['err']:8.4f} {r['purity']:8.4f} "
              f"{r['drift']:9.2e}")
    return gate(med, runs)


if __name__ == "__main__":
    with np.errstate(all="ignore"):
        sys.exit(main())
