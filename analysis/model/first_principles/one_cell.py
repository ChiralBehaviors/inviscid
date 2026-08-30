"""one_cell -- a single jitterbug, from first principles: it passes through
the VE with its joints intact, the VE is where it chooses a sense, and the
octahedron is a dead end that reopens one way only.

THE QUESTION (owner, 2026-08-30). "Let's take this from first principles
again. Let's just use a single VE, a single jitterbug. Can you now drive a
triangle from -60 to +60?" And then: at the VE "you have a choice ... as to
which chirality you chose to transition to the octa, correct?" And then:
"let's use the current -60 as 0. then we can move -60 and +60 from that
point ... Correct?"

WHAT IS MEASURED HERE. The cell is `plates.corners(a)`: eight rigid
plates, each spinning about its own axis by sigma*(a - 60) and riding it at
height Z cos a. The twelve joints are the twelve coincident corner pairs,
read once at a generic angle and then only CHECKED, never re-read.

  R1  Driven -60 -> +60 the same twelve corner pairs coincide on every
      frame (gap ~1e-16); struts stay sqrt2; twelve distinct vertices except
      at +-60 where they sit as six (Fuller's "12 vertices congruent as 6");
      no two struts that share no joint come within 1.03 of each other.
  R2  The driven plate's spin about its own axis is a, its height Z cos a.
  R3  At -60 and at +60 the twelve joints merge onto the six points in
      DIFFERENT pairings. Which two joints share a point is the physical
      record of which way the cell closed.
  R4  As a FREE LINKAGE -- eight rigid plates, 48 dof, twelve point joints,
      no symmetric-family assumption -- the cell has six internal freedoms
      at the octahedron, the same as anywhere. The opening it came in by is
      in the tangent space and walks; the mirror-sense opening (spin the
      other way and rise) violates the joints at first order, and what
      survives projection onto them is not an opening. The octahedron is an
      END.
  R5  At the VE both senses are real motions of the free linkage: either
      way closes it. The VE is the point with the choice.

So a single cell's range is VE-centred, octahedron <- VE -> octahedron,
120 degrees with the choice in the middle; the answer to "use -60 as 0 and
move +-60 from there" is no.

T2: [23791]. Page: "One Jitterbug" (pages/export_one_cell.py).
"""
from __future__ import annotations

import sys

import numpy as np

from analysis.model import plates as Z
from analysis.model.first_principles import geometry as G

FACES = Z.faces()
AXES = np.array([f[2] for f in FACES])
SIG = np.array([f[3] for f in FACES])
DRIVEN = 0


def spin_and_height(C, f):
    """Spin of plate f about its own axis relative to the VE, and its height."""
    u = AXES[f]
    C0 = Z.corners(0.0)
    v0 = C0[f][0] - C0[f].mean(0)
    v1 = C[f][0] - C[f].mean(0)
    v0p, v1p = v0 - u * (v0 @ u), v1 - u * (v1 @ u)
    ang = float(np.degrees(np.arctan2(np.cross(v0p, v1p) @ u, v0p @ v1p)))
    return ang, float(C[f].mean(0) @ u)


# ---- the free linkage: eight rigid plates, twelve point joints ---------------
def _jac(C, J):
    M = np.zeros((36, 48))
    for r, ((f1, c1), (f2, c2)) in enumerate(J):
        for sgn, f, c in ((1, f1, c1), (-1, f2, c2)):
            p = C[f][c] - C[f].mean(0)
            M[3 * r:3 * r + 3, 6 * f:6 * f + 3] += sgn * np.eye(3)
            px = np.array([[0, -p[2], p[1]], [p[2], 0, -p[0]], [-p[1], p[0], 0]])
            M[3 * r:3 * r + 3, 6 * f + 3:6 * f + 6] += sgn * (-px)
    return M


def _gaps(C, J):
    return np.concatenate([C[f1][c1] - C[f2][c2] for ((f1, c1), (f2, c2)) in J])


def _move(C, v, eps):
    out = np.empty_like(C)
    for f in range(8):
        t, w = v[6 * f:6 * f + 3] * eps, v[6 * f + 3:6 * f + 6] * eps
        th = np.linalg.norm(w)
        c0 = C[f].mean(0)
        if th > 0:
            k = w / th
            K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
            R = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K
        else:
            R = np.eye(3)
        out[f] = (R @ (C[f] - c0).T).T + c0 + t
    return out


def _project(C, J, iters=40):
    for _ in range(iters):
        g = _gaps(C, J)
        if np.abs(g).max() < 1e-13:
            break
        C = _move(C, -np.linalg.lstsq(_jac(C, J), g, rcond=None)[0], 1.0)
    return C


def linkage_nullity(C, J):
    s = np.linalg.svd(_jac(C, J), compute_uv=False)
    return 48 - int((s > s[0] * 1e-9).sum())


def height_rate(a):
    """d(height)/da of every plate along its own axis, per radian, at a."""
    h = lambda x: float(Z.corners(x)[0].mean(0) @ AXES[0])
    return (h(a + 1e-4) - h(a - 1e-4)) / (2e-4 * np.pi / 180.0)


def symmetric_field(sense, rate):
    """Spin every plate about its own axis with the given per-plate sense
    (1 rad per unit) and move it along the axis at `rate` per unit."""
    v = np.zeros(48)
    for f in range(8):
        v[6 * f:6 * f + 3] = AXES[f] * rate
        v[6 * f + 3:6 * f + 6] = AXES[f] * sense[f]
    return v


def walk_height(C, J, v, steps=4, deg=5.0):
    """Walk the linkage along v, projecting back onto the joints; report the
    driven plate's height change and the worst joint gap."""
    h0 = float(C[DRIVEN].mean(0) @ AXES[DRIVEN])
    worst = 0.0
    for _ in range(steps):
        C = _project(_move(C, v, np.radians(deg)), J)
        worst = max(worst, float(np.abs(_gaps(C, J)).max()))
    return float(C[DRIVEN].mean(0) @ AXES[DRIVEN]) - h0, worst


def gate():
    checks, out = [], {}
    A = checks.append
    J = G.joints()

    # ---- R1: driven -60 -> +60, joints permanent -------------------------------
    gaps, nverts, struts, cross, clear = [], {}, [], [], {}
    for a in np.arange(-60.0, 60.0 + 1e-9, 2.0):
        C = Z.corners(a)
        gaps.append(G.joint_gap(C, J))
        nverts[a] = len(G.classes(C))
        struts += [float(np.linalg.norm(C[f][i] - C[f][j])) for f in range(8) for i, j in ((0, 1), (1, 2), (2, 0))]
        S = G.cell_struts(C)
        cross.append(G.crossings([s[1] for s in S], [s[2] for s in S], [s[0] for s in S]))
        if a in (-30.0, 0.0, 30.0):
            clear[a] = G.self_clearance(C)
    out["R1"] = (max(gaps), sorted(set(nverts.values())), min(struts), max(struts), max(cross), clear)
    A(("R1  DRIVEN -60 -> +60 THE TWELVE JOINTS ARE THE SAME TWELVE CORNER PAIRS on "
       "every frame; struts stay sqrt2; twelve vertices except at +-60 where they sit "
       "as six; no strut ever crosses another (pairs of struts CONVERGE onto each other "
       "toward the octahedra -- 24 edges congruent as 12 -- but never cross); mid-fold "
       "no two struts that share no joint come within 1.0.",
       max(gaps) < 1e-12 and nverts[-60.0] == 6 and nverts[60.0] == 6
       and all(nverts[a] == 12 for a in nverts if abs(a) < 60)
       and abs(min(struts) - np.sqrt(2)) < 1e-9 and abs(max(struts) - np.sqrt(2)) < 1e-9
       and max(cross) == 0 and min(clear.values()) > 1.0,
       f"max joint gap {max(gaps):.1e}, vertices {{-60: {nverts[-60.0]}, 0: {nverts[0.0]}, +60: {nverts[60.0]}}}, "
       f"strut {min(struts):.6f}..{max(struts):.6f}, crossings {max(cross)}, clearance at -30/0/30 "
       f"{clear[-30.0]:.3f}/{clear[0.0]:.3f}/{clear[30.0]:.3f}",
       "gap < 1e-12; 6 / 12 / 6; sqrt2; 0; > 1.0"))

    # ---- R2: the driven plate ------------------------------------------------
    err_spin, err_h = 0.0, 0.0
    for a in (-60.0, -30.0, 0.0, 30.0, 60.0):
        sp, h = spin_and_height(Z.corners(a), DRIVEN)
        err_spin = max(err_spin, abs(sp - a))
        err_h = max(err_h, abs(h - Z.Z * np.cos(np.radians(a))) if hasattr(Z, "Z") else 0.0)
    h60, h0 = spin_and_height(Z.corners(60.0), DRIVEN)[1], spin_and_height(Z.corners(0.0), DRIVEN)[1]
    out["R2"] = (err_spin, h60, h0)
    A(("R2  THE DRIVEN PLATE SPINS ABOUT ITS OWN AXIS BY EXACTLY a, and rides it from "
       "height 0.577 at the octahedron to 1.155 at the VE and back.",
       err_spin < 1e-9 and abs(h60 - 1 / np.sqrt(3)) < 1e-9 and abs(h0 - 2 / np.sqrt(3)) < 1e-9,
       f"spin error {err_spin:.1e}; height {h60:.4f} at 60, {h0:.4f} at 0",
       "0; 0.5774, 1.1547"))

    # ---- R3: the two octahedra pair the joints differently --------------------
    jid = {c: i for i, j in enumerate(J) for c in j}
    merged = lambda a: sorted(tuple(sorted({jid[c] for c in g})) for g in G.classes(Z.corners(a)))
    m_minus, m_plus = merged(-60.0), merged(60.0)
    out["R3"] = (m_minus, m_plus)
    A(("R3  AT -60 AND AT +60 THE TWELVE JOINTS MERGE ONTO THE SIX POINTS IN DIFFERENT "
       "PAIRINGS. Which joints share a point records which way the cell closed.",
       len(m_minus) == 6 and len(m_plus) == 6 and m_minus != m_plus
       and all(len(p) == 2 for p in m_minus + m_plus),
       f"-60: {m_minus}; +60: {m_plus}", "six pairs each, and not the same six"))

    # ---- R4: the octahedron is a dead end, as a free linkage -------------------
    C60 = Z.corners(60.0)
    n60, n30 = linkage_nullity(C60, J), linkage_nullity(Z.corners(30.0), J)
    Jm = _jac(C60, J)
    rise = -height_rate(60.0)                 # a decreasing from 60: the plates rise
    v_open = symmetric_field(-SIG, rise)      # the family's own opening, spin -sigma
    v_mirror = symmetric_field(+SIG, rise)    # spin the other way, rise the same
    first_open = np.linalg.norm(Jm @ v_open) / np.linalg.norm(v_open)
    first_mirror = np.linalg.norm(Jm @ v_mirror) / np.linalg.norm(v_mirror)
    dh_open, g_open = walk_height(C60, J, v_open)
    dh_mirror, g_mirror = walk_height(C60, J, v_mirror)
    out["R4"] = (n60, n30, first_open, first_mirror, dh_open, dh_mirror, g_open, g_mirror)
    A(("R4  AS A FREE LINKAGE THE OCTAHEDRON HAS SIX INTERNAL FREEDOMS, LIKE ANYWHERE, "
       "AND STILL OPENS ONE WAY: the opening it came in by is in the tangent space and "
       "walks (height rises > 0.3 in four 5-degree steps); the mirror-sense opening "
       "violates the joints at first order, and what survives projection is not an "
       "opening (height change < 0.02).",
       n60 == 12 and n30 == 12 and first_open < 1e-8 and first_mirror > 1.0
       and dh_open > 0.3 and abs(dh_mirror) < 0.02 and g_open < 1e-12 and g_mirror < 1e-12,
       f"nullity {n60} at 60, {n30} at 30; |Jv|/|v| open {first_open:.1e}, mirror {first_mirror:.2f}; "
       f"height change open {dh_open:+.3f}, mirror {dh_mirror:+.4f}; joint gaps {g_open:.0e}, {g_mirror:.0e}",
       "12, 12; < 1e-8, > 1; > 0.3, < 0.02; both < 1e-12"))

    # ---- R5: at the VE both senses are real -----------------------------------
    C0 = Z.corners(0.0)
    J0 = _jac(C0, J)
    res = {}
    # at the VE the height is stationary (d/da of Z cos a is 0 there): the
    # tangent is pure spin, either sense, and the sinking is second order
    rate0 = height_rate(0.0)
    for name, sense in (("one way", -SIG), ("the other", +SIG)):
        v = symmetric_field(sense, rate0)
        dh, g = walk_height(C0, J, v)
        res[name] = (np.linalg.norm(J0 @ v) / np.linalg.norm(v), dh, g)
    out["R5"] = res
    A(("R5  AT THE VE BOTH SENSES ARE REAL MOTIONS OF THE FREE LINKAGE: the tangent there "
       "is pure spin (the height is stationary at the VE), either sense is in the tangent "
       "space, and walking either way the joints hold and the cell closes (height falls "
       "at second order). The VE is the point with the choice.",
       all(r[0] < 1e-8 and r[1] < -0.04 and r[2] < 1e-12 for r in res.values()),
       "; ".join(f"{k}: |Jv|/|v| {r[0]:.1e}, height {r[1]:+.3f}, gap {r[2]:.0e}" for k, r in res.items()),
       "both: < 1e-8, < -0.04 (four 5-degree steps), < 1e-12"))
    return checks, out


def main():
    np.set_printoptions(precision=6, suppress=True)
    with np.errstate(all="ignore"):
        print("=" * 78)
        print("one_cell -- a single jitterbug: through the VE with its joints, the VE the "
              "choice, the octahedron a dead end")
        print("=" * 78)
        checks, out = gate()
        bad = 0
        for name, ok, got, want in checks:
            bad += 0 if ok else 1
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
            print(f"        got {got}")
            print(f"        want {want}")
        print()
        print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
        print("   * ONE CELL'S RANGE IS VE-CENTRED: octahedron <- VE -> octahedron, 120")
        print("     degrees, and the choice of sense is made at the VE, not at an end.")
        print("   * AN OCTAHEDRON WITH ITS JOINTS AS THEY ARE OPENS ONE WAY, even with")
        print("     six internal freedoms. 'Use -60 as 0 and move +-60 from there' is no.")
        print("   * NOTHING HERE IS THE MEDIUM. It is what one cell can do; what a tied")
        print("     array can do is ring.py and vertex_point.py.")
        print()
        print("  ALL CHECKS PASSED." if not bad
              else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
        return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
