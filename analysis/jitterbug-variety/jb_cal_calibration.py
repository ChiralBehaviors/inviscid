"""jb_cal -- calibrate the model before trusting anything the medium says.

*** VOID IN PART -- OWNER DECISION 2026-08-27. ***
The strut-compliance fork is REJECTED (T2 inviscid [23562]); the model of record
is a RIGID-STRUT LINKAGE, not an elastic network.

WHAT STANDS HERE: everything that is rigidity or kinematics. The unit cell, the
bar list with its lattice offsets, the Bloch phase convention, the reciprocal
lattice, and every ZERO-MODE COUNT. Nullity of a rigidity matrix does not depend
on stiffness, so the mechanism counts and the six <110> floppy lines are
statements about the rigid linkage and are exactly the right model.

WHAT IS VOID: anything read as a FREQUENCY -- elastic bands, sound speeds, phase
speeds, gaps, dispersion. Those exist only once the bars are springs.


Three transport measurements in a row failed their own validity checks, and the
common feature of all three was that NONE of them had a known answer. A number
came out, nothing could say whether it was right, and the argument moved on to
whether the physics was surprising. That is backwards. This file measures two
things whose answers are known in advance, so that a failure is unambiguous.

C1, HERE: does the model reproduce itself? A periodic SUPERCELL of N1 x N2 x N3
unit cells, diagonalised directly, must have exactly the spectrum you get by
evaluating the Bloch matrix at the N1*N2*N3 commensurate wavevectors and taking
the union. This is not an approximation that holds well or badly -- it is an
identity, true to machine precision or else something is wrong. It exercises
the unit cell, the bar list with its lattice offsets, the Bloch phase
convention, and the reciprocal lattice, which is to say almost everything the
rest of the programme stands on.

C2, SEPARATELY: does the measurement APPARATUS work? That is the monatomic
chain, and it does not live here.

WHY AN IDENTITY AND NOT A COMPARISON. A free-boundary patch was the obvious
thing to diagonalise, since that is the object jb_aa actually integrates. But a
free patch has surface modes and no exact quantisation, so the comparison
against bulk bands is a judgement call about how close is close -- and a
judgement call is what we are trying to eliminate. The supercell gives up the
free boundary and buys an exact answer. Both are worth having; only one of them
can fail cleanly.

WHAT WOULD MAKE THIS VACUOUS. An identity that holds no matter what is not
evidence, and this one has an obvious failure mode: if the supercell were built
by calling the same Bloch code, it would agree with itself for free. It is not.
The supercell is assembled directly from the bar list in real space with no
wavevector anywhere in its construction, and C3 corrupts a single bar's lattice
offset to show the comparison goes red when the two constructions disagree.
"""

from __future__ import annotations

import itertools as it
import sys

import numpy as np

import jb_hc_honeycomb as HC

A_REF = -30.0

#: Supercell extents to test. Two different shapes, because an identity that
#: holds only for a cube would be a coincidence about the cube.
SUPERCELLS = ((3, 2, 2), (4, 3, 2))

IDENTITY_TOL = 1e-9      # max |w2_supercell - w2_folded| over the whole spectrum
MUTATION_MIN = 1e-3      # a corrupted offset must break the identity by AT LEAST


def supercell(N, a=A_REF):
    """Assemble the periodic supercell in REAL SPACE. No wavevector appears.

    Returns (K, ncells, n_unit, nbars_unit). Nodes are indexed (cell, slot)
    flattened cell-major; a bar (i, j, R) in the unit cell becomes, for every
    cell c, a bar from (c, i) to (c + R mod N, j) carrying the same geometric
    direction the unit cell gives it."""
    P, bars, A, _slots = HC.unit_cell(a)
    n = len(P)
    cells = list(it.product(*(range(m) for m in N)))
    cidx = {c: k for k, c in enumerate(cells)}
    nc = len(cells)
    rows = []
    for c in cells:
        for (i, j, Rv) in bars:
            Rw = np.array(Rv, dtype=float) * A
            d = P[i] - (P[j] + Rw)
            d = d / np.linalg.norm(d)
            cj = tuple((c[t] + Rv[t]) % N[t] for t in range(3))
            row = np.zeros(3 * n * nc)
            row[3 * (cidx[c] * n + i):3 * (cidx[c] * n + i) + 3] += d
            row[3 * (cidx[cj] * n + j):3 * (cidx[cj] * n + j) + 3] -= d
            rows.append(row)
    S = np.array(rows)
    return S.T @ S, nc, n, len(bars)


def folded(N, a=A_REF):
    """The same spectrum the other way: Bloch at every commensurate wavevector.

    Allowed k for a supercell of N cells are k_d = 2*pi*m_d / (N_d * A), which
    is the reciprocal lattice divided by N."""
    P, bars, A, _slots = HC.unit_cell(a)
    out = []
    for m in it.product(*(range(t) for t in N)):
        kv = 2 * np.pi * np.array([m[t] / (N[t] * A) for t in range(3)])
        M = HC.bloch(P, bars, A, kv)
        out.append(np.linalg.eigvalsh(M.conj().T @ M))
    return np.sort(np.concatenate(out))


def mutated_folded(N, a=A_REF):
    """Control: the SAME fold with one bar's lattice offset corrupted.

    The identity must break. If it does not, the comparison is insensitive to
    the lattice offsets and proves nothing about them."""
    P, bars, A, _slots = HC.unit_cell(a)
    bad = list(bars)
    i, j, Rv = bad[0]
    bad[0] = (i, j, tuple(np.array(Rv) + np.array([1, 0, 0])))
    out = []
    for m in it.product(*(range(t) for t in N)):
        kv = 2 * np.pi * np.array([m[t] / (N[t] * A) for t in range(3)])
        M = HC.bloch(P, bad, A, kv)
        out.append(np.linalg.eigvalsh(M.conj().T @ M))
    return np.sort(np.concatenate(out))


def zero_modes(K, tol=1e-8):
    w = np.linalg.eigvalsh(K)
    return int(np.sum(w < tol)), w


def gamma_mechanism(a=A_REF):
    """WHAT the Gamma mechanism is, not just how many there are.

    Returns (n_zero, n_mech, weight of the mechanism inside the phase span,
    per-phase-direction weight inside the mechanism span). The last one is the
    two-sided part: if BOTH phase directions were soft the medium would have
    two free channels, and if neither were the soft mode would be something
    else entirely. Either would be a finding."""
    P, bars, A, slots = HC.unit_cell(a)
    M = HC.bloch(P, bars, A, np.zeros(3))
    K = (M.conj().T @ M).real
    w, V = np.linalg.eigh(K)
    null = V[:, w < 1e-8]
    T = np.zeros((3, K.shape[0]))
    for d in range(3):
        T[d, d::3] = 1.0
    T /= np.linalg.norm(T, axis=1, keepdims=True)
    Q = null - T.T @ (T @ null)
    u, sv, _ = np.linalg.svd(Q, full_matrices=False)
    mech = u[:, sv > 1e-8]
    B = np.real(HC.phase_basis(P, slots, A, np.zeros(3)))
    B = B / np.linalg.norm(B, axis=0, keepdims=True)
    Pspan = B @ np.linalg.pinv(B)
    Mspan = mech @ mech.T
    return (null.shape[1], mech.shape[1],
            [float(mech[:, i] @ Pspan @ mech[:, i]) for i in range(mech.shape[1])],
            [float(B[:, j] @ Mspan @ B[:, j]) for j in range(B.shape[1])])


def zero_locus(N, a=A_REF):
    """Where in the zone the zero modes sit. Returns [(k/(pi/A), count)]."""
    P, bars, A, _slots = HC.unit_cell(a)
    out = []
    for m in it.product(*(range(t) for t in N)):
        kv = 2 * np.pi * np.array([m[t] / (N[t] * A) for t in range(3)])
        w = np.linalg.eigvalsh(HC.bloch(P, bars, A, kv).conj().T
                               @ HC.bloch(P, bars, A, kv))
        nz = int((w < 1e-8).sum())
        if nz:
            out.append((tuple(round(float(x) / (np.pi / A), 3) for x in kv), nz))
    return out


def measure(N):
    K, nc, n, nb = supercell(N)
    w2_direct = np.sort(np.linalg.eigvalsh(K))
    w2_fold = folded(N)
    w2_mut = mutated_folded(N)
    nz, _ = zero_modes(K)
    return dict(
        N=N, ncells=nc, n_unit=n, nbars_unit=nb,
        dof=3 * n * nc, nbars=nb * nc,
        n_direct=len(w2_direct), n_fold=len(w2_fold),
        gap=float(np.abs(w2_direct - w2_fold).max()),
        mut=float(np.abs(w2_direct - w2_mut).max()),
        zero=nz,
        locus=zero_locus(N),
        wmax=float(np.sqrt(max(w2_direct.max(), 0.0))),
    )


def gate(runs):
    checks = []
    R = checks.append

    R(("C1  the supercell has exactly N1*N2*N3 copies of the unit cell -- "
       "node and bar counts, so a mis-assembled cell cannot hide",
       all(r["dof"] == 3 * r["n_unit"] * r["ncells"]
           and r["nbars"] == r["nbars_unit"] * r["ncells"] for r in runs),
       f"dof {[r['dof'] for r in runs]}, bars {[r['nbars'] for r in runs]}",
       "n_unit x ncells"))

    R(("C1b both constructions produce the SAME NUMBER of modes, which is the "
       "precondition for comparing them elementwise at all",
       all(r["n_direct"] == r["n_fold"] == r["dof"] for r in runs),
       f"direct {[r['n_direct'] for r in runs]} vs "
       f"folded {[r['n_fold'] for r in runs]}", "equal, and equal to dof"))

    R(("C2  THE IDENTITY: a real-space supercell diagonalised directly has "
       "exactly the spectrum of the Bloch matrix folded over the commensurate "
       "wavevectors. No wavevector appears anywhere in the supercell's "
       "construction, so this is two independent routes to one answer",
       all(r["gap"] < IDENTITY_TOL for r in runs),
       f"worst |w2_direct - w2_folded| = {max(r['gap'] for r in runs):.3e}",
       f"< {IDENTITY_TOL:.0e}"))

    R(("C3  MUTATION PROBE -- corrupt ONE bar's lattice offset and the same "
       "comparison must go RED. Without this, C2 could be an identity that "
       "holds however the offsets are wired, and would say nothing about them",
       all(r["mut"] > MUTATION_MIN for r in runs),
       f"smallest break {min(r['mut'] for r in runs):.3e}",
       f"> {MUTATION_MIN:.0e}"))

    R(("C4  ZERO MODES ARE COUNTED, NOT ASSUMED. A periodic framework has 3 "
       "rigid translations and no rotations; anything above that is a "
       "mechanism and is the interesting number. Reported rather than "
       "subtracted -- `nullity - 6` once returned -3 in this programme, which "
       "is nonsense that gets explained away instead of investigated",
       all(r["zero"] >= 3 for r in runs),
       f"zero modes {[r['zero'] for r in runs]} of dof "
       f"{[r['dof'] for r in runs]}; mechanisms beyond translation "
       f"{[r['zero'] - 3 for r in runs]}", ">= 3 translations"))

    nz, nm, mech_in_phase, phase_in_mech = gamma_mechanism()
    R(("C6  WHAT THE MECHANISM IS. The one Gamma mechanism beyond rigid "
       "translation lies ENTIRELY in the phase span, and exactly ONE of the "
       "two phase directions is soft while the other is not. TWO-SIDED: two "
       "soft directions would mean two free channels, none would mean the "
       "soft mode is something other than the jitterbug coordinate. Either "
       "would be a finding, and neither is what the medium does",
       nm == 1 and mech_in_phase[0] > 0.99
       and sum(1 for x in phase_in_mech if x > 0.99) == 1
       and sum(1 for x in phase_in_mech if abs(x) < 0.01) == 1,
       f"{nz} zero modes at Gamma = 3 translations + {nm} mechanism; "
       f"mechanism sits {mech_in_phase[0]:.4f} inside the phase span; "
       f"phase directions sit {[round(x, 4) for x in phase_in_mech]} inside "
       f"the mechanism", "1 mechanism, all phase, exactly one direction soft"))

    R(("C7  AND THE ONLY OTHER ZERO IS ON THE ZONE BOUNDARY, at a wavevector "
       "with two components at pi/A. That is jb_bt's band-touching locus, "
       "arrived at here by a construction with no wavevector anywhere in it",
       all(all(nzz == 4 for kk, nzz in r["locus"] if all(abs(x) < 1e-6 for x in kk))
           and all(sum(1 for x in kk if abs(abs(x) - 1.0) < 1e-6) >= 2
                   for kk, nzz in r["locus"]
                   if any(abs(x) > 1e-6 for x in kk)) for r in runs),
       "; ".join(f"N={r['N']}: " + ", ".join(f"k/(pi/A)={kk} -> {nzz}"
                                             for kk, nzz in r["locus"])
                 for r in runs),
       "4 at Gamma, the rest with >=2 components at the boundary"))

    R(("C8  PRINTED NOT GATED: the spectrum is bounded and real, so the "
       "supercell is a stable framework rather than a numerically sick one",
       True, f"max omega {[round(r['wmax'], 4) for r in runs]}", "printed"))

    return checks


def main():
    np.set_printoptions(precision=6, suppress=True)
    print("=" * 78)
    print("jb_cal -- C1: does the model reproduce itself?")
    print("=" * 78)
    print("  A periodic supercell assembled in REAL SPACE, diagonalised")
    print("  directly, against the Bloch matrix folded over the commensurate")
    print("  wavevectors. An identity, not an approximation: it holds to")
    print("  machine precision or something is wrong. Nothing here is a")
    print("  property of the medium -- the answer is known in advance, which")
    print("  is the entire point after three measurements whose answers were")
    print("  not.")
    print()

    runs = [measure(N) for N in SUPERCELLS]

    print(f"  {'supercell':>12} {'cells':>6} {'dof':>6} {'bars':>6} "
          f"{'identity':>11} {'mutated':>11} {'zero':>5}")
    for r in runs:
        print(f"  {str(r['N']):>12} {r['ncells']:6d} {r['dof']:6d} "
              f"{r['nbars']:6d} {r['gap']:11.3e} {r['mut']:11.3e} "
              f"{r['zero']:5d}")
    print()

    checks = gate(runs)
    bad = 0
    for name, ok, got, want in checks:
        tag = "PASS" if ok else "FAIL"
        bad += 0 if ok else 1
        print(f"  {tag}  {name}")
        print(f"        got {got}")
        print(f"        want {want}")
    print()
    print("  WHAT THIS DOES AND DOES NOT LICENSE.")
    print("   * It licenses the linearised machinery: unit cell, bar offsets,")
    print("     Bloch phase convention, reciprocal lattice. Two independent")
    print("     constructions agree exactly, and C3 shows the agreement is")
    print("     sensitive to the offsets rather than automatic.")
    print("   * It licenses NOTHING about transport. No time integration")
    print("     happens here and no packet is launched. The instrument that")
    print("     measures transport is calibrated separately, against a chain")
    print("     whose dispersion is known in closed form.")
    print("   * The free-boundary patch jb_aa integrates is a DIFFERENT")
    print("     object: surface modes, no exact quantisation. C2 says the")
    print("     bulk model is self-consistent, not that the patch is bulk.")
    print()
    print("  ALL CHECKS PASSED." if not bad else f"  {bad} CHECK(S) FAILED.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
