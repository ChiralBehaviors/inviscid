"""lift -- G0: the sign-continuous SU(2) lift of a sampled SO(3) path, and
its null controls. No plate result (G1, G2, G7) is trusted until this
instrument passes its own gates.

THE QUESTION (analysis/notes/su2_boundary_conditions.md section 8.4, bead
inviscid-qjo). The G1 measurement -- does one 360-degree overdrive cycle
carry a plate's quaternion to -q -- is only as good as the lift code that
produces the sign. Section 8.4 demands null controls on the lift itself: a
trivial loop must give +q, a bare 2 pi rotation -q, a bare 4 pi rotation +q,
and the sampling must be fine enough that the nearest-sign choice is
meaningful. This module is that instrument and those controls.

WHAT THE INSTRUMENT IS. `quat_from_matrix` converts one rotation matrix to
a unit quaternion by Shepperd's max-branch method, so the conversion stays
well-conditioned at and near trace -1 (180-degree rotations -- the
inside-out region of the overdrive passes through these). `lift` walks a
sampled path R_0 .. R_n, converting each frame and flipping the sign of any
quaternion whose dot with its predecessor is negative; it REFUSES the path
(ValueError) if any successive pair of frames differs by more than
`max_step_deg` of rotation, because the nearest-sign choice degrades toward
a coin flip as the step approaches 180 degrees. `holonomy` requires the
path to close in SO(3), lifts it, and returns the sign of q_end . q_start:
+1 on a contractible loop, -1 on the nontrivial class.

  R1  The conversion round-trips: quat -> matrix -> quat -> matrix at 1e-12
      over seeded random rotations, and at exact and near 180-degree
      rotations about coordinate and skew axes (the trace ~ -1 branches).
  R2  Trivial loops lift to +q: the constant path, a one-axis excursion out
      and back, and a two-axis non-planar excursion out and back.
  R3  A bare 2 pi rotation lifts to -q about every axis tried (x, y, z,
      skew, seeded random).
  R4  A bare 4 pi rotation lifts to +q about the same axes.
  R5  The sign is a class invariant as measured: reversing the 2 pi loop
      and conjugating it frame-by-frame by a random rotation both leave
      the holonomy at -1.
  R6  The instrument refuses what it cannot measure: a 2 pi path sampled
      at 60-degree steps raises, an unclosed path raises in `holonomy`,
      and the worst step `lift` reports on a 0.5-degree path is 0.5 deg.

Vocabulary note (section 8.5): the -1 is classical Z2 holonomy of a path,
a fact about SO(3), not a spin sector and not anything the medium's joints
remember (section 8.2).

T2: [23865]. Ref: su2_boundary_conditions.md section 8 (read before 1-7).
"""
from __future__ import annotations

import sys

import numpy as np


def quat_from_matrix(R):
    """Unit quaternion (w, x, y, z) of a rotation matrix, Shepperd max-branch."""
    R = np.asarray(R, float)
    m00, m11, m22 = R[0, 0], R[1, 1], R[2, 2]
    t = m00 + m11 + m22
    if t >= max(m00, m11, m22):
        r = np.sqrt(1.0 + t)
        q = np.array([0.5 * r, (R[2, 1] - R[1, 2]) / (2 * r),
                      (R[0, 2] - R[2, 0]) / (2 * r), (R[1, 0] - R[0, 1]) / (2 * r)])
    elif m00 >= max(m11, m22):
        r = np.sqrt(1.0 + m00 - m11 - m22)
        q = np.array([(R[2, 1] - R[1, 2]) / (2 * r), 0.5 * r,
                      (R[1, 0] + R[0, 1]) / (2 * r), (R[0, 2] + R[2, 0]) / (2 * r)])
    elif m11 >= m22:
        r = np.sqrt(1.0 - m00 + m11 - m22)
        q = np.array([(R[0, 2] - R[2, 0]) / (2 * r), (R[1, 0] + R[0, 1]) / (2 * r),
                      0.5 * r, (R[2, 1] + R[1, 2]) / (2 * r)])
    else:
        r = np.sqrt(1.0 - m00 - m11 + m22)
        q = np.array([(R[1, 0] - R[0, 1]) / (2 * r), (R[0, 2] + R[2, 0]) / (2 * r),
                      (R[2, 1] + R[1, 2]) / (2 * r), 0.5 * r])
    return q / np.linalg.norm(q)


def matrix_from_quat(q):
    """Rotation matrix of a unit quaternion (w, x, y, z)."""
    w, x, y, z = np.asarray(q, float) / np.linalg.norm(q)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


#: named control axes; `rot` and `path` also take any 3-vector
AXES = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1), "skew": (1, 1, 1)}


def rot(axis, deg):
    """Rotation matrix: `deg` degrees about `axis` (a 3-vector or an AXES name)."""
    k = np.asarray(AXES.get(axis, axis) if isinstance(axis, str) else axis, float)
    k = k / np.linalg.norm(k)
    th = np.radians(deg)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def lift(Rs, max_step_deg=30.0):
    """Sign-continuous SU(2) lift of a sampled SO(3) path.

    Returns (q, worst): the (n, 4) lifted quaternions and the largest
    rotation step between successive frames, in degrees. Raises ValueError
    when any step exceeds `max_step_deg` -- the nearest-sign choice is only
    trustworthy well below the 180-degree ambiguity, so a coarse path is
    refused rather than silently lifted.
    """
    Rs = np.asarray(Rs, float)
    q = np.empty((len(Rs), 4))
    q[0] = quat_from_matrix(Rs[0])
    worst = 0.0
    for i in range(1, len(Rs)):
        qi = quat_from_matrix(Rs[i])
        d = float(qi @ q[i - 1])
        step = 2.0 * np.degrees(np.arccos(min(1.0, abs(d))))
        worst = max(worst, step)
        if step > max_step_deg:
            raise ValueError(
                f"lift step {step:.1f} deg between frames {i - 1} and {i} "
                f"exceeds {max_step_deg:g} deg -- sample the path finer")
        q[i] = -qi if d < 0 else qi
    return q, worst


def holonomy(Rs, max_step_deg=30.0, closure_tol=1e-9):
    """+1 or -1: the lifted sign at the end of a CLOSED sampled SO(3) path."""
    Rs = np.asarray(Rs, float)
    gap = float(np.abs(Rs[-1] - Rs[0]).max())
    if gap > closure_tol:
        raise ValueError(f"path is not closed in SO(3): |R_end - R_start| = {gap:.1e}")
    q, _ = lift(Rs, max_step_deg)
    s = float(q[-1] @ q[0])
    if abs(abs(s) - 1.0) > 1e-6:
        raise ValueError(f"lift end is not +-start: dot = {s:.6f}")
    return 1 if s > 0 else -1


def path(axis, degs):
    """The sampled path of rotations about one fixed axis through `degs`."""
    return np.array([rot(axis, d) for d in degs])


def gate():
    checks = []
    A = checks.append
    rng = np.random.default_rng(0)

    def rand_axis():
        v = rng.normal(size=3)
        return v / np.linalg.norm(v)

    # ---- R1: the conversion round-trips, including the trace ~ -1 branches ----
    worst_rand = 0.0
    for _ in range(200):
        q0 = rng.normal(size=4)
        q0 /= np.linalg.norm(q0)
        R = matrix_from_quat(q0)
        worst_rand = max(worst_rand, float(np.abs(matrix_from_quat(quat_from_matrix(R)) - R).max()))
    worst_pi = 0.0
    for ax in list(AXES.values()) + [rand_axis() for _ in range(4)]:
        for deg in (180.0, 180.0 - 1e-4, -180.0 + 1e-4, 179.999999):
            R = rot(ax, deg)
            worst_pi = max(worst_pi, float(np.abs(matrix_from_quat(quat_from_matrix(R)) - R).max()))
    A(("R1  THE CONVERSION ROUND-TRIPS: quat -> matrix -> quat -> matrix over 200 seeded "
       "random rotations, and at exact and near 180-degree rotations about coordinate, "
       "skew and random axes (the trace ~ -1 branches of Shepperd's method).",
       worst_rand < 1e-12 and worst_pi < 1e-9,
       f"max |R'' - R| random {worst_rand:.1e}, near-pi {worst_pi:.1e}",
       "random < 1e-12, near-pi < 1e-9"))

    # ---- R2: trivial loops lift to +q -----------------------------------------
    const = np.array([np.eye(3)] * 721)
    up = np.linspace(0, 90, 181)
    one_axis = path("x", np.concatenate([up, up[::-1][1:]]))
    leg1 = path("x", up)
    leg2 = np.array([rot("x", 90) @ rot("y", d) for d in up])
    two_axis = np.concatenate([leg1, leg2[1:], leg2[::-1][1:], leg1[::-1][1:]])
    triv = {"constant": holonomy(const), "one-axis out-and-back": holonomy(one_axis),
            "two-axis out-and-back": holonomy(two_axis)}
    A(("R2  TRIVIAL LOOPS LIFT TO +q: the constant path, a one-axis excursion out and "
       "back, and a two-axis non-planar excursion out and back all return with sign +1.",
       all(v == 1 for v in triv.values()),
       ", ".join(f"{k} {v:+d}" for k, v in triv.items()),
       "all +1"))

    # ---- R3: a bare 2 pi lifts to -q ------------------------------------------
    degs2 = np.linspace(0, 360, 721)
    two_pi = {k: holonomy(path(ax, degs2)) for k, ax in AXES.items()}
    two_pi["random"] = holonomy(path(rand_axis(), degs2))
    A(("R3  A BARE 2 PI ROTATION LIFTS TO -q about every axis tried: the sampled loop "
       "closes in SO(3) and its lift returns with sign -1 (the nontrivial class).",
       all(v == -1 for v in two_pi.values()),
       ", ".join(f"{k} {v:+d}" for k, v in two_pi.items()),
       "all -1"))

    # ---- R4: a bare 4 pi lifts back to +q -------------------------------------
    degs4 = np.linspace(0, 720, 1441)
    four_pi = {k: holonomy(path(ax, degs4)) for k, ax in AXES.items()}
    four_pi["random"] = holonomy(path(rand_axis(), degs4))
    A(("R4  A BARE 4 PI ROTATION LIFTS TO +q about the same axes: twice around the "
       "nontrivial loop is contractible and the lift says so.",
       all(v == 1 for v in four_pi.values()),
       ", ".join(f"{k} {v:+d}" for k, v in four_pi.items()),
       "all +1"))

    # ---- R5: the sign is a class invariant as measured ------------------------
    loop = path("x", degs2)
    g = rot(rand_axis(), float(rng.uniform(10, 170)))
    conj = np.einsum("ij,njk,kl->nil", g, loop, g.T)
    inv = {"reversed": holonomy(loop[::-1]), "conjugated": holonomy(conj)}
    A(("R5  THE SIGN IS A CLASS INVARIANT AS MEASURED: reversing the 2 pi loop and "
       "conjugating it frame-by-frame by a random rotation both leave the holonomy "
       "at -1 -- the lift is reading the loop's class, not its parametrisation.",
       all(v == -1 for v in inv.values()),
       ", ".join(f"{k} {v:+d}" for k, v in inv.items()),
       "both -1"))

    # ---- R6: the instrument refuses what it cannot measure --------------------
    coarse_raised = unclosed_raised = False
    try:
        lift(path("x", np.linspace(0, 360, 7)))          # 60-degree steps
    except ValueError:
        coarse_raised = True
    try:
        holonomy(path("x", np.linspace(0, 350, 701)))    # never returns
    except ValueError:
        unclosed_raised = True
    _, worst_fine = lift(path("x", degs2))               # 0.5-degree steps
    A(("R6  THE INSTRUMENT REFUSES WHAT IT CANNOT MEASURE: a 2 pi path sampled at "
       "60-degree steps raises, an unclosed path raises in holonomy, and the worst "
       "step reported on the 0.5-degree path is 0.5 deg -- the guard measures what "
       "it claims to guard.",
       coarse_raised and unclosed_raised and abs(worst_fine - 0.5) < 1e-9,
       f"coarse raised {coarse_raised}, unclosed raised {unclosed_raised}, "
       f"worst fine step {worst_fine:.6f} deg",
       "True, True, 0.5 +- 1e-9"))
    return checks


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("lift -- G0: the SU(2) lift of a sampled SO(3) path, and its null controls")
    print("=" * 78)
    checks = gate()
    bad = 0
    for name, ok, got, want in checks:
        bad += 0 if ok else 1
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")
    print()
    print("  WHAT THIS LICENSES AND WHAT IT DOES NOT.")
    print("   * THE INSTRUMENT IS CALIBRATED: G1/G2/G7 may lift plate paths with this")
    print("     module and quote the sign, at sampling the guard accepts.")
    print("   * A -q DOWNSTREAM IS A FACT ABOUT A PATH IN SO(3) -- classical Z2")
    print("     holonomy, not a spin sector, and not anything the medium's point")
    print("     joints remember (note section 8.2).")
    print("   * NOTHING HERE TOUCHES THE MODEL: no plate data was lifted. G1's")
    print("     expected -q remains a prediction (note section 8.3), not a result.")
    print()
    print("  ALL CHECKS PASSED." if not bad
          else f"  !! {bad} CHECK(S) FAILED -- this is a bug report, not a measurement.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
