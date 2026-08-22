"""Step V: the TRANSVERSE CURVATURE of V along the symmetric path.
Bead inviscid-qvf.4.

WHY THE FILE IS NOT CALLED `transverse_stability`, which is what the bead asks
for. Off a critical point an eigenvalue of (Hess V, g) is a CURVATURE INVARIANT
of the pair (V, g), not a frequency -- qvf.9 corollary (i). Nothing oscillates
where the gradient does not vanish, and everywhere on the symmetric path except
its critical points the system is MOVING, not resting. A dynamical stability
verdict for a moving trajectory is the normal variational equation along that
trajectory (a Floquet-type problem), which needs an equation of motion
integrated, and this project has never integrated one. So the honest deliverable
is the SIGN of the curvature of V in the five directions transverse to the path,
which is what the bead's question -- "is the symmetric sector a saddle?" -- was
reaching for, and the file is named for that. See V10 for what is and is not
licensed by it.

WHAT IS MEASURED. At each configuration on the symmetric path, the Riemannian
Hessian of V with respect to the MOMENTUM-FREE (mechanical-connection) mass
metric, exactly as jb_u builds it, decomposed into the three symmetry blocks of
the chiral tetrahedral group. The symmetric path is that group's fixed-point
set, so the group acts on the 6-D internal tangent space at every point of the
path and 6 = 1 + 2 + 3 is forced. The SINGLET is the path's own tangent; the
DOUBLET and TRIPLET are the five transverse directions. Their sign is the
answer.

THE BLOCK LABELLING IS THE HARD PART, AND jb_u's DOES NOT SURVIVE THIS SWEEP.
jb_u labels blocks by the eigenspaces of the mass metric, which is
kernel-independent and therefore stable against the D/T ordering flips that
defeat sort-position labelling. That is strictly better than sorting and it is
still not enough here: measured below in closed form (V2), the doublet and the
triplet mass values are EQUAL AT a = 60 EXACTLY -- g's spectrum there is
[1/32 x 5, 5/96] -- so the mass metric cannot separate them in a band around the
octahedron, and `blocks_by_irrep` correctly REFUSES rather than guessing. This
file therefore builds the labelling from the SYMMETRY CHARACTER instead: the 12
rotations of the chiral tetrahedral group are constructed explicitly, their
action on the 72-D ambient configuration space is built and verified to fix the
configuration exactly, they are pushed down to the 6-D chart by the
metric-orthogonal projection, and the three irrep projectors are assembled from
the character table. That labelling is defined at every point of the path
including a = 60, does not depend on any eigenvalue ordering, and its ranks
(1, 2, 3) are a check rather than an assumption.

WHAT IS REUSED, NOT REWRITTEN. `jb_u_riemannian_hessian` (charts, the shared
finite-difference stencil, the momentum-free form, g/dg/Gamma/Hess),
`jb_j_internal_frame` + `jb_k_hull_hessian.aligned_frame` (the chart and the
gauge, with basis direction 0 pinned to the symmetric-path tangent),
`jb_o_kernel_family` and `jb_q_strut_kernels` (the raw kernels on the two
primitives), `jb_t_modes_primitive_offpath` (the second chart). The correctness
of the Riemannian Hessian itself is gated by jb_u's own 17-row gate and is NOT
re-derived here; what is gated here is everything this file adds.

FOUR DECLARATIONS on every number: KERNEL, MASS MODEL, PRIMITIVE (raw vertex or
raw strut-midpoint), METRIC FORM (section or momentum-free). The metric form is
MOMENTUM-FREE throughout, and V8 measures -- for this sweep, rather than
inheriting jb_u's two-angle measurement -- that the section form coincides with
it in this chart at every swept configuration, INCLUDING in the derivative dg,
which is the quantity that actually needed checking because the stencil steps
transverse to the path even when the base point is on it.

ABSOLUTE VALUES ARE A CONVENTION (coupling 1, total mass 1/2, R = 1). SIGNS ARE
NOT, and RATIOS ARE NOT. The deliverable of this file is a sign structure; the
ratios are reported beside it; no absolute curvature scale is offered as a
measurement.

THE SWEEP RANGE IS [0, 90), NOT THE BEAD'S [-60, 60]. The bead's acceptance
criterion predates USER DECISION 16 and quotes an admissibility verdict that has
since been WITHDRAWN: configuration space is the whole variety. a -> 180 - a is
an exact isometry (1e-15), so [0, 90] is the fundamental domain, and V7 measures
the isometry acting on the curvature blocks rather than assuming it. Two
boundaries are handled explicitly rather than by omission: a = 60, where the
inverse-power kernels genuinely diverge because vertices merge in pairs (V6),
and a = 90, where the chart's Newton solve fails because the local dimension is
7 rather than 6 (V7). Neither is silently dropped and neither is claimed.

THE TWO GUARD BANDS ARE ENTERED, NOT MERELY EXCLUDED. Together they are about
3.3% of the fundamental domain, and a verdict that says "throughout the
fundamental domain" cannot rest on a region the sweep never visits. V6b runs all
twenty inverse-power combinations across |a - 60| < 1 and V7(a-ter) runs them
across a > 89, both at two step sizes, both carrying the per-block deviation so
that a sign taken from the band is qualified exactly as one taken outside it.
The guard itself is applied BY KERNEL FAMILY: the Gaussians do not diverge at
a = 60, so V3c bisects them straight through the band and locates two turnovers
that were previously reported as unlocatable.

AND THE TURNOVER ANGLES CARRY AN ERROR BAR (V3d). They are the most quotable
output of this file and were, for one revision, its only entirely ungated one.
The bisection runs to 1e-07 degrees, which is NOT the precision of the answer:
V3d measures each root moving with the step size at up to 1.05e-04 degrees, gates
that reproducibility, gates the slope that makes the root well posed, and gates
the ABSOLUTE per-block deviation converted into an angular uncertainty -- the
relative deviation used everywhere else is structurally blind at a root, where
its denominator goes to zero. FOUR DECIMALS is the quotable precision.

NOT part of the Maven build. Nothing under src/ is touched. Run from the repo
root: `python3 -W ignore analysis/jitterbug-variety/jb_v_transverse_curvature.py`
(the scripts add their own directory to sys.path).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import itertools as it

import numpy as np

import jb_cache

#: The importable name of THIS module -- a literal, because `__name__` is
#: "__main__" under `python jb_v_transverse_curvature.py` and a prefetch
#: worker in a fresh interpreter must be able to re-import it by name.
_MODULE = "jb_v_transverse_curvature"

from jb_a_family import corners
from jb_j_internal_frame import Frame
from jb_k_hull_hessian import aligned_frame
from jb_o_kernel_family import A_ICO
from jb_r_mass_metric import MODELS, position_jacobian
from jb_s_frequency_spectrum import RAW_KERNELS
from jb_t_modes_primitive_offpath import OriginFrame, origin_position_jacobian
from jb_u_riemannian_hessian import (ChartUnmeasurable, FrameChart, Geometry,
                                     IrrepLabelError, Stencil, WEIGHTS,
                                     blocks_by_irrep, forms, rigid_fields,
                                     strut_potential, vertex_potential)

np.seterr(all="ignore")      # spurious numpy/BLAS matmul warnings; see jb_r


# --------------------------------------------------------------------------
# configurations of record. Every one of these is quoted so that a
# re-measurement has something to DISAGREE with, and every one is gated.
# --------------------------------------------------------------------------

A_VE = 0.0
A_POLE = 60.0                 # vertices merge in pairs; inverse powers diverge
A_BRANCH = 90.0               # local dimension 7; the chart's Newton solve dies

# The recorded second minimum for raw Thomson on raw VERTICES (jb_t S5d, memo
# THE SPECTRUM). Defined HERE, locally, and not imported -- jb_u's post-mortem
# records a validator "verifying" a hole by mutating an imported name that was
# never defined locally, so the mutation never applied and a clean exit was read
# as confirmation. A local constant is one a mutation probe can actually reach.
A_2ND_THOMSON_VERTEX = 75.262042
# The recorded 6-D local MAXIMUM of raw gauss 0.5 on raw vertices (same source).
A_MAX_GAUSS05_VERTEX = 62.706500

# The limiting slope of the PATH gradient's linear collapse at the branch point:
# (path gradient in the g-norm) / (90 - a) as a -> 90, for raw 1/r^1 (Thomson) on
# raw VERTICES with point masses, momentum-free. FOUR DECLARATIONS, because it is
# a number and this file's rule has no exceptions.
#
# THIS IS A REPRODUCIBILITY ANCHOR, NOT AN INDEPENDENT ONE, and saying so is the
# point: no prior record carries it, so it is THIS FILE'S OWN measurement written
# down as a local constant so that a later change to the collapse law, to the
# gradient, or to the ratio's denominator has something to disagree with. Held
# locally (never imported) for the same reason A_2ND_THOMSON_VERTEX is: a
# mutation probe must be able to reach it. Its detection floor is the ratio's own
# convergence, |r(89.99) - r(89.9)| = 1.8e-04, so the 1e-04 relative criterion on
# it resolves a change of about one part in ten thousand -- and the mutation this
# exists to catch, halving the denominator, lands five thousand times past that.
A90_RATE = 2.9628

PRIOR = {
    # jb_s S1, reproduced by jb_u U1 to 1.835e-07 relative. Written here in
    # (Doublet, Triplet, Singlet) order -- the character order this file uses --
    # rather than ascending, precisely because ascending order is what this file
    # refuses to label modes by.
    "DTS_VE_point": (43.163565, 45.104546, 60.370536),
    "DTS_VE_lamina": (69.061704, 90.209093, 241.482143),
    # jb_u U2b, the icosahedron, momentum-free, raw vertex, point mass.
    "grad_ico": 3.057691,
    # closed form: the momentum-free mass metric at the octahedron.
    "m_doublet": 1.0 / 32.0,
    "m_singlet_at_60": 5.0 / 96.0,
    # jb_s S1 / jb_u U2b(ii-ter), the VE ratios in OMEGA (i.e. sqrt of the
    # generalised eigenvalue), as (T/D, S/D). Quoted so that V4 can show its
    # own CURVATURE ratios squaring back onto them -- the two conventions
    # differ by a square root and a reader comparing them directly would find
    # a false discrepancy.
    "omega_ratios_VE": {"point": (1.022237, 1.182644),
                        "lamina": (1.142895, 1.869924)},
}


# --------------------------------------------------------------------------
# gate thresholds. Declared here, read by `gate()`, and each labelled DERIVED
# or FITTED, because the two are not the same kind of number. Every "(measured
# X)" annotation is the value THIS code produces at H_MAIN, re-measured after
# the last edit; jb_u's post-mortem records three annotations that survived a
# refactor that had invalidated them.
# --------------------------------------------------------------------------

# The main-sweep step size. DERIVED, not chosen here: jb_u's U5 sweeps h over
# eight decades for the same chart, the same stencil and the same Riemannian
# Hessian, and its measured minimum is h = 3e-4 (T2 qvf.9-riemannian-hessian.md,
# "U5's own minimum is h=3e-4"). jb_u quotes its own deliverable at 1e-3 with a
# stated 2.53x error bar rather than re-running 36 combinations; this file has
# no such sunk cost, so it sits at the measured optimum. V9 sweeps h anyway and
# prints which thresholds hold at which step size, because a threshold in a
# quantity carrying O(h^2) truncation is fitted to its h until a sweep says
# otherwise.
H_MAIN = 3e-4

# Angles within POLE_GUARD degrees of a = 60 are excluded from the MAIN table
# and given their own section (V6, V6b). They are NOT silently dropped: V6
# measures them at three step sizes and reports what does and does not survive
# there, and V6b sweeps ALL TWENTY inverse-power combinations across the band so
# that the coverage claim in V10 is made of measurements rather than of the
# band's absence.
# FITTED, and the fit is stated: at h = 3e-4 the worst per-block scalarity over
# all 36 combinations on the 0.5-degree grid is 1.64e-04 with the band excluded
# and reaches 4.13e-01 inside it (2.40e+00 at h = 1e-3), because V varies by
# orders of magnitude across a stencil of any usable width when a pole is a
# hundredth of a degree away. V6 prices the band from its own rows rather than
# from this comment.
#
# THE GUARD IS APPLIED BY KERNEL FAMILY, NOT BY ANGLE ALONE, and that is a
# correction earned by review rather than a refinement. The justification above
# is an INVERSE-POWER pole: at a = 60 the twelve shared vertices merge in pairs,
# so some r_ij is exactly zero and 1/r^p diverges. A GAUSSIAN DOES NOT DIVERGE
# THERE -- V6's own table has the gauss 0.5 rows sitting at dev ~2e-08 at every
# step size straight through the band while the Thomson rows blow up to 4.13e-01
# -- so applying the same guard to Gaussian rows suppresses measurements the
# file's own machinery resolves cleanly. It did: two of the located turnovers
# (gauss 1.0 / strut / point, triplet, and gauss 1.5 / vertex / point, triplet)
# came back "STRADDLES a=60 / not refined", and one of them is not even near the
# pole -- it sits at ~59.68 and was unlocated only because 59.5 had been struck
# from the grid. `_pole_guarded` below is where the family test lives; V3c calls
# it instead of testing the angle alone.
POLE_GUARD = 1.0


def _pole_guarded(kernel):
    """Does the a = 60 guard band apply to THIS kernel?

    True for the inverse powers, which genuinely diverge at the octahedron.
    False for the Gaussians, which are finite there -- V6 measures them at
    dev ~2e-08 inside the band, five orders inside this file's own scalarity
    criterion, at every step size in that table.

    The MAIN GRID still excludes the band for every kernel, because the grid is
    shared: one site serves all 36 combinations, and an inverse-power row inside
    the band would redden the scalarity gate for all of them. What this function
    governs is where a band exclusion would DISCARD a measurement that resolves
    -- the bisection of V3c, and the coverage tables of V6b.
    """
    return kernel in INVERSE_POWER


# The bisection tolerance for the turnover angles of V3c, in DEGREES. It is not
# the precision of the answer and must not be read as one: V3d measures the
# turnover's movement with the step size h at 1.05e-04 deg, which is the real
# error bar, so the honest quotable precision is FOUR decimals. This constant is
# set two decades below that so that the bisection is not what limits the number.
# It was 1e-4 -- i.e. the same size as the physical error bar -- while the table
# printed six decimals, which put ~5e-05 of pure bisection slop into a field
# displaying 1e-06 and made the two headline angles in the record wrong in their
# fifth and sixth decimals.
BISECT_TOL = 1e-7

# The same idea at the other end of the fundamental domain, and it was NOT
# anticipated: a = 90 is a second pole for the raw STRUT-MIDPOINT primitive
# (V7c measures the strut midpoints colliding linearly as 2(90 - a) radians)
# while the raw VERTEX primitive is perfectly well behaved there. The main
# sweep therefore stops one degree short of the branch point and V7 approaches
# it explicitly, per primitive and at three step sizes. FITTED, and the fit is
# stated: V7(a-bis) measures 1/r^3 on STRUTS at h = 3e-4 going 2.51e-05
# (a = 89) -> 9.06e-05 (89.5) -> 1.87e-03 (89.9) -> 2.26e+00 (89.99) while the
# same kernel on VERTICES stays at 2.9e-08 throughout. The guard exists because
# the first run of this file did NOT have it and produced a 1/r^3 strut row at
# a = 89.99 whose scalarity was 2.7 -- larger than the quantity it qualifies.
BRANCH_GUARD = 1.0

# The metric form, in ONE place. qvf.9 corollary (iii): the momentum-free
# (mechanical-connection) form is required for anything transverse, and it is
# the fourth declaration on every number below. This constant exists because the
# mutation matrix caught the alternative: `Site.__init__` and the `site()`
# memoiser each carried their own `kind="horizontal"` default, so mutating the
# one on `Site` changed nothing -- the wrapper's default shadowed it -- and the
# probe reported a clean exit that was the expected result under BOTH
# hypotheses. That is jb_u's A_ICO failure in a new costume: the substitution
# APPLIED (the count was asserted) and still reached no live code. The general
# guard, now in the harness, is that a mutation whose OUTPUT is byte-identical
# to the baseline is not a mutation.
METRIC_FORM = "horizontal"

TOL = {
    # --- the symmetry machinery this file adds ---
    "rep_fix": 1e-14,      # |T x(a) - x(a)|, the group really fixes the
                           # configuration                   (measured 2.220e-16)
                           # DERIVED: the entries of T are exact 0/+-1 and the
                           # configuration is exact, so the baseline is roundoff
                           # on a 72-vector of O(1) entries.
    "rep_mult": 1e-12,     # rho(g1)rho(g2) = rho(g1 g2)      (measured 9.992e-16)
                           # DERIVED: same reason -- an algebraic identity of
                           # exact matrices.
    "rep_iso": 1e-9,       # rho^T gx rho = gx, relative, in
                           # the CENTROID chart               (measured 9.250e-16)
                           # DERIVED: the symmetry is an isometry of the mass
                           # metric, so this is an identity of closed-form
                           # objects and its baseline is roundoff. It is gated
                           # at 1e-9 rather than at roundoff so that it stays
                           # meaningful in the ORIGIN chart, whose section is
                           # built by a different construction.
    "proj_alg": 1e-10,     # P_A + P_E + P_F = I and P^2 = P  (measured 6.661e-16)
                           # DERIVED: character-orthogonality identity.
    "proj_ev": 1e-9,       # projector eigenvalues are 0 or 1 (measured 3.331e-16)
    # --- the block labelling ---
    "tangent_leak": 1e-6,  # |P_A e0 - e0|_g / |e0|_g at EVERY
                           # swept a -- the bead's "verify,
                           # do not assume"                   (measured 7.751e-09)
                           # FITTED to H_MAIN, and V9 SWEEPS IT AND REFUTES THE
                           # REASON THIS COMMENT USED TO GIVE. The old text said
                           # "the statistic carries g, hence O(h^2)"; V9's own
                           # column reads 1.855e-09 at h = 3e-3, 1e-3, 3e-4 AND
                           # 1e-4 -- FLAT to four significant figures. It is flat
                           # because P_A is built entirely from the CLOSED-FORM D
                           # and gx (see `Characters`), so only the |.|_g norm
                           # touches the stencil and the leak it measures is a
                           # property of the closed-form projector. The threshold
                           # is therefore NOT h-sensitive; V9 measures the
                           # dynamic range of each swept column and says which
                           # ones actually vary, rather than asserting that all
                           # three do.
    "tangent_teeth": 0.9,  # the SAME statistic on a transverse
                           # chart direction must be O(1), or
                           # the row above cannot fail        (measured 1.000e+00)
                           # DERIVED: P_A annihilates the transverse complement
                           # exactly, so the statistic is exactly 1 there.
    "grad_transverse": 3e-3,
                           # (doublet + triplet part of the
                           # gradient) / (path part), over the
                           # rows where the path part is not
                           # itself vanishing                 (measured 4.139e-04)
                           # DERIVED to be identically zero: dV is an invariant
                           # covector at a point the group fixes, so it lies in
                           # the A-isotypic component, which V0/V1 measure to be
                           # one-dimensional and equal to the path direction.
                           # What is measured is therefore a finite-difference
                           # floor, and it is RELATIVE because the floor scales
                           # with the gradient -- an absolute threshold set from
                           # the well-behaved rows (3e-07) is exceeded by twenty
                           # decades at a = 61 on 1/r^12 struts, where the
                           # gradient itself is 2.4e+24. FITTED to H_MAIN, and
                           # measured to be pure O(h^2): the worst row runs
                           # 4.59e-03 / 4.14e-04 / 4.60e-05 at
                           # h = 1e-3 / 3e-4 / 1e-4, swept in place below.
    "grad_path_teeth": 1e0,
                           # ... and the PATH part must NOT
                           # vanish, or the row above is a
                           # statement about a zero gradient  (measured 2.430e+24)
    "split_teeth": 1e0,    # the transverse/path ratio of the
                           # SAME split, run on a covector whose
                           # gradient vector IS a transverse
                           # chart direction               (measured 3.479e+10)
                           # DERIVED: g^-1 (g e1) = e1 lies entirely in the
                           # transverse complement, which P_A annihilates
                           # exactly, so the path part is roundoff and the ratio
                           # is enormous. Gated at 1 rather than at the measured
                           # value because the claim is only that the split CAN
                           # see a transverse component, which is what the row
                           # above needs and does not itself establish.
    "grad_crit_floor": 1e-3,
                           # a row whose PATH gradient is below
                           # this is at a critical point and
                           # has no direction to be transverse
                           # TO; such rows are reported
                           # separately, by absolute size     (42 of 6480 rows)
    "v4_rank": 3,          # under the Klein subgroup the "A"
                           # projector has rank 3, not 1 -- the
                           # teeth of the rank test           (measured 3)
                           # DERIVED from the character table: chi_E restricted
                           # to V4 is the trivial character twice, so the
                           # V4-fixed subspace is A + E, of dimension 3.
    # --- the measurement ---
    "scalarity": 1e-3,     # max |Hess|block - lambda I| / |lambda|
                           # over the whole main sweep x 36   (measured 1.642e-04)
                           # FITTED to H_MAIN and measured to be pure O(h^2):
                           # V9 shows it falling by ~10x per 3.16x drop in h.
                           # This is the falsifier for the entire per-block
                           # construction -- every lambda here is
                           # trace(Hess|block)/mult, an AVERAGE, which is a
                           # block eigenvalue only if Hess is scalar on the
                           # block.
    "ve_record": 1e-6,     # VE reproduction against jb_s S1,
                           # relative                         (measured 2.895e-08)
                           # FITTED, and capped by the record: the recorded
                           # constants carry six decimals, so nothing below
                           # ~8e-09 relative could mean anything.
    "ico_record": 1e-5,    # |dV| at the icosahedron vs the
                           # record -- the ONE configuration
                           # anchor that is not a = 0         (measured 2.827e-07)
    "second_min": 1e-3,    # |dV| at the recorded second
                           # minimum: it must actually be a
                           # critical point                   (measured 8.845e-07)
                           # FITTED, and its DETECTION FLOOR is measured rather
                           # than argued: d|dV|/da at that angle is 1.21 per
                           # degree (|dV| = 1.208e-03 at A + 0.001 deg and
                           # 8.993e-02 at A + 0.075 deg), so 1e-3 detects an
                           # angle error of ~8e-4 degrees. A 0.1% error in the
                           # recorded angle -- the A_ICO-style mutation jb_u's
                           # post-mortem is about -- lands at 9.0e-02, five
                           # decades above the baseline.
    "isometry": 1e-6,      # a -> 180-a maps the curvature
                           # blocks to themselves, relative   (measured 1.255e-07)
                           # DERIVED that it must hold (the map is an exact
                           # isometry of the vertex set, 1e-15); the THRESHOLD
                           # is fitted to H_MAIN.
    "isometry_teeth": 1e-2,
                           # the same statistic on a
                           # DELIBERATELY MISMATCHED pair must
                           # be large                         (measured 4.100e-01)
    "metric_form": 1e-5,   # max |Gamma_section - Gamma_horizontal| / |Gamma|
                           # over the whole sweep             (measured 7.183e-09)
                           # This is the seed's explicit instruction: CONFIRM
                           # for this sweep rather than inherit jb_u's
                           # two-angle 2.566e-08.
    "gauge_leak": 1e-12,   # closed-form |Z^T W D| on the path (measured 7.634e-16)
                           # DERIVED: exact arithmetic, no finite differences.
    "chart_agree": 1e-5,   # centroid vs origin chart, block
                           # values, relative, momentum-free  (measured 5.266e-07)
                           # FITTED to H_MAIN; it is jb_u's u2_horiz_rel
                           # criterion applied to this file's own statistic.
    "chart_teeth": 1e-2,   # ... and the SECTION form must
                           # FAIL that, or the row above is
                           # vacuous                          (measured 2.614e-01)
    "form_in_use": 1e-2,   # the ambient form the sweep actually
                           # used must DIFFER from the section
                           # form, or the declaration
                           # "momentum-free" is unfalsifiable  (measured 1.042e-01)
                           # DERIVED: Wh = W - WZ(Z^T W Z)^-1 Z^T W, and the
                           # subtracted term is a W-orthogonal projector onto a
                           # 6-D subspace of a 72-D space, so it is O(1)
                           # relative to W and not a small correction. This row
                           # exists because a mutation that switched the Site's
                           # default form to `section` made section V8(a)
                           # compare a thing to itself and the run exited 0.
    "ratio_convention": 1e-6,
                           # sqrt of this file's CURVATURE
                           # ratios at the VE must equal the
                           # record's OMEGA ratios            (measured 4.113e-07)
                           # DERIVED that it must hold (omega^2 = lambda at a
                           # critical point); the threshold is the record's
                           # own six-decimal quantisation, ~5e-07 relative.
    # --- the turnover angles: V3c locates them, V3d prices them ---
    "root_h_spread": 1e-3,
                           # |root(h) - root(H_MAIN)| in DEGREES,
                           # worst over every located turnover
                           # and the two flanking it          (measured 1.050e-04)
                           # FITTED to the measurement, and it is the number that
                           # sets HOW MANY DECIMALS OF A TURNOVER ANGLE ARE
                           # QUOTABLE. The bisection is run to BISECT_TOL = 1e-7
                           # so that the bisection is NOT the limiting error; the
                           # limiting error is this one, and it is ~1e-04 deg, so
                           # four decimals is the honest precision and the record
                           # quotes four. Before this row existed the bisection
                           # stopped at 1e-4 and the table printed six decimals,
                           # i.e. 5e-05 of bisection slop displayed in a 1e-06
                           # field, and a re-measurement did not reproduce the
                           # fifth or sixth decimal of either headline angle.
    "root_slope_teeth": 1e-2,
                           # min |d lambda / d a| at a located
                           # turnover, per DEGREE              (measured 3.784e-01)
                           # THE TEETH of the row above: if lambda were flat
                           # through the crossing then every angle in the bracket
                           # would be a root, the h-spread would be small for a
                           # reason that has nothing to do with the root being
                           # well determined, and the reproducibility row would
                           # be vacuous. A non-zero slope is what makes the root
                           # a well-posed quantity to reproduce.
    "root_dev_angle": 1e-3,
                           # the ABSOLUTE per-block scalar
                           # deviation at the root, converted
                           # to degrees by the measured slope  (measured 4.642e-06)
                           # THIS IS THE ERROR CONTROL THE RELATIVE `scalarity`
                           # CANNOT PROVIDE HERE, and the gap is structural, not
                           # a matter of tightness: `dev` is |Hb - lam I| divided
                           # by |lam|, and the bisection evaluates precisely where
                           # lam -> 0, so the relative statistic carries NO
                           # information about the measurements the roots are made
                           # of. Dividing the ABSOLUTE deviation by |dlam/da|
                           # converts it into the angular uncertainty it actually
                           # implies, in the same units as the root.
    # --- coverage inside the two guard bands (V6b, V7 a-ter) ---
    "band_dev": 1e-3,      # the SAME criterion as `scalarity`,
                           # applied to the rows measured INSIDE
                           # the guard bands, so that a sign
                           # taken from the band is qualified
                           # exactly as one taken outside it   (measured 5.167e-04
                           #                                    a=60 band;
                           #                                    4.763e-04 a>89)
    "band_sign_announce": 1e-1,
                           # a NEGATIVE reported inside a guard
                           # band must carry a per-block
                           # deviation at least 100x the band
                           # criterion above               (measured 1.44e+01
                           #                                a=60; 9.86e-01 a>89)
                           # DERIVED FROM THE CRITERION IT PROTECTS rather than
                           # fitted to the observation: `band_dev` = 1e-3 is the
                           # point below which a block value is trusted, so a
                           # sign disagreement two decades PAST it is one the
                           # dev column is loudly refusing to stand behind. The
                           # a > 89 measurement sits at 9.86e-01, i.e. just
                           # under the "deviation exceeds the block value" bar
                           # that V6's own row uses -- which is why that bar is
                           # NOT what is asserted here. Stating it as a multiple
                           # of the criterion is the honest form; stating it as
                           # ">= 1" and then discovering the measurement is 0.986
                           # would have meant loosening a threshold to fit a
                           # number, which is what gate lesson 3 forbids.
    "guard_fraction": 0.05,
                           # the two guard bands together may
                           # exclude at most this fraction of
                           # the fundamental domain        (measured 3.333e-02)
                           # A POLICY BOUND, and it is labelled one. It exists
                           # because both guards were measured to be
                           # UNCONSTRAINED FROM ABOVE: widening POLE_GUARD to 3
                           # or BRANCH_GUARD to 5 exits 0 with zero red rows AND
                           # IMPROVES the reported worst scalarity (1.077e-04
                           # against 1.642e-04), so the gate as it stood rewarded
                           # hiding more of the domain. There is no derivation
                           # for a particular width -- the poles are points and
                           # any positive width is a fit -- so what is asserted
                           # is that the fit stays small: 5% of [0,90) is 4.5
                           # degrees, already far more than a pole neighbourhood
                           # needs, and the two bands are additionally COVERED by
                           # V6b and V7(a-ter) rather than merely excluded.
    # --- V5(d), the a = 90 collapse of the path gradient ---
    "a90_rate": 1e-4,      # |(path grad / (90-a)) at 89.99
                           # divided by A90_RATE - 1|      (measured 1.595e-05)
    "a90_linear": 1e-3,    # |ratio(89.99) - ratio(89.9)|,
                           # the LINEARITY of the collapse (measured 1.830e-04)
    # --- the a = 90 band, priced the way V6 prices the a = 60 band ---
    "branch_vertex": 1e-5, # worst per-block scalarity on the
                           # raw VERTEX primitive inside
                           # a > 89                            (measured 3.527e-07)
    "branch_strut_teeth": 1e0,
                           # ... and the raw STRUT primitive
                           # must FAIL it, or "the guard exists
                           # because of the strut pole" is a
                           # claim the rows do not support     (measured 2.265e+00)
    "pole_linearity": 1e-6,
                           # |min strut-midpoint separation /
                           # 2(90-a) - 1| at a = 89.99         (measured 5.077e-09)
                           # DERIVED to be zero in the limit; the threshold is set
                           # from the O((90-a)) residual the finite angle leaves.
                           # This row exists because a mutation of the pole LAW
                           # itself -- 2(90-a) -> 3(90-a) -- exited 0 with every
                           # gate row green: the ratio column was printed and
                           # quoted in the record and asserted nowhere.
    "pole_vertex_floor": 0.8,
                           # ... and the min VERTEX separation
                           # must NOT collapse at a = 89.99, or
                           # "the pole belongs to ONE primitive"
                           # is not what the rows show         (measured 0.8165)
    # --- the h sweep ---
    "h_window": 2,         # at least this many of the five step
                           # sizes must meet ALL FOUR of the
                           # swept thresholds                  (measured 2 of 5)
                           # It is 2 and not 3 because the scalarity criterion
                           # fails at h = 3e-3 AND at h = 1e-3 on the worst of
                           # the 36 combinations. That is disclosed rather than
                           # traded away by loosening `scalarity`: the two
                           # passing rows are H_MAIN and one step below it, and
                           # V9 prints the per-threshold h ranges so the reader
                           # can see which criterion does the excluding.
                           # TWO THINGS THIS ROW IS NOT, both measured in V9 and
                           # both disclosed there rather than left to be inferred:
                           # (1) the four thresholds do not all constrain it.
                           # The tangent leak is FLAT to four figures (dynamic
                           # range 1.00e+00) and the VE reproduction passes at
                           # every step, so the window is set by TWO criteria
                           # pulling from OPPOSITE ends: `scalarity`, which is
                           # truncation and rules out the coarse steps, and
                           # `isometry`, which is roundoff and rules out the
                           # finest. The isometry threshold was added to this
                           # sweep because an independent whole-gate h sweep
                           # found the gate RED at h = 1e-4 -- finer than the
                           # quoted step -- on precisely that row, while this
                           # section reported h = 1e-4 as acceptable;
                           # (2) it has ZERO MARGIN -- 2 measured against a
                           # criterion of 2. A future worst-combination ~6x worse
                           # than today's would take h = 3e-4 out of the passing
                           # set and redden this row.
    "h_shape": 1e1,        # ... and the SCALARITY column must
                           # VARY, or it is not measuring
                           # discretisation                   (measured 1.23e+03)
                           # Named for scalarity specifically: V9 prints the
                           # dynamic range of all four swept columns and only
                           # this one is required to have shape, because it is
                           # the only one whose h dependence is truncation.
}

# Character table of the chiral tetrahedral group T, keyed by trace(R), which
# separates the three conjugacy classes exactly: 3 = identity, -1 = the three
# C2, 0 = the eight C3. "E" is the REAL two-dimensional representation, i.e. the
# sum of the complex-conjugate pair E' + E'', so its character is
# omega^k + omega^-k = -1 on C3 and 2 elsewhere.
CHI = {"A": {3.0: 1.0, -1.0: 1.0, 0.0: 1.0},
       "E": {3.0: 2.0, -1.0: 2.0, 0.0: -1.0},
       "F": {3.0: 3.0, -1.0: -1.0, 0.0: 0.0}}
# Weight in P = (d/|G|) sum chi* rho. d = 1 for A; d = 1 for EACH of E', E'' and
# the two are summed by using the real character above; d = 3 for F.
DIMW = {"A": 1, "E": 1, "F": 3}
MULT = {"A": 1, "E": 2, "F": 3}
BLOCKS = ("A", "E", "F")
LONG = {"A": "singlet (path tangent)", "E": "doublet", "F": "triplet"}


# --------------------------------------------------------------------------
# the chiral tetrahedral group, its action on the ambient space, and the fact
# that both are INDEPENDENT of a
# --------------------------------------------------------------------------

def _sgn(x):
    """'+' / '-' for a block value, and '?' for one that is not a number.

    NOT COSMETIC. `np.nan > 0` is False, so the obvious
    `"+" if lam > 0 else "-"` PRINTS A NEGATIVE SIGN FOR A NaN -- i.e. a block
    the labelling failed to produce at all is reported in the headline sign map
    as a measured negative, which is a wrong answer rather than a missing one.
    `block_values` already returns (nan, inf, inf) for a rank mismatch so that
    the scalarity falsifier reddens; the SIGN had no such treatment. Every sign
    formatter in this file goes through here.
    """
    if not np.isfinite(x):
        return "?"
    return "+" if x > 0 else "-"


def tetrahedral_rotations():
    """The 12 rotations of the chiral tetrahedral group T, as 3x3 matrices.

    DERIVED, not tabulated. The eight bodies sit one per octant, labelled by the
    sign triple s, and jb_a's forced chirality gives body s the twist
    sigma = s_x s_y s_z. A signed permutation matrix with permutation p and
    signs e maps octant s to an octant whose sign product is (prod e)(prod s),
    so sigma is preserved iff prod e = +1; and det = sign(p)(prod e), so a
    ROTATION preserving sigma needs sign(p) = +1 as well. Even permutations
    (3) times even sign patterns (4) is exactly 12 -- which is why the jitterbug
    family's symmetry is T and not the full chiral octahedral O of the
    cuboctahedron. The count is a consequence here, not an input.
    """
    out = []
    for p in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        for e in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
            R = np.zeros((3, 3))
            for k in range(3):
                R[p[k], k] = float(e[p[k]])
            out.append(R)
    return out


def _ambient_action(R, X):
    """The 72x72 matrix of R acting on configurations, plus the match residual.

    R permutes the eight bodies and, within a body, the three corners. The
    permutation is found by MATCHING rather than derived from a labelling
    convention, so a wrong convention cannot hide: the returned residual is the
    worst |R X[i][perm[j]] - X[i'][j]| over the whole configuration, and it is
    gated. A group element that did not actually preserve the configuration
    would show up here as an O(1) number, not as a silently wrong permutation.
    """
    T = np.zeros((72, 72))
    worst = 0.0
    for i in range(8):
        Y = (R @ X[i].T).T
        best = None
        for ip in range(8):
            for perm in it.permutations(range(3)):
                d = float(np.abs(Y[list(perm)] - X[ip]).max())
                if best is None or d < best[0]:
                    best = (d, ip, perm)
        d, ip, perm = best
        worst = max(worst, d)
        for j in range(3):
            for dd in range(3):
                for de in range(3):
                    T[9 * ip + 3 * j + dd, 9 * i + 3 * perm[j] + de] = R[dd, de]
    return T, worst


_ROTS = tetrahedral_rotations()
# Built ONCE at the vector equilibrium and reused at every a. That reuse is a
# claim -- the body and corner permutations do not depend on the angle -- and it
# is CHECKED at every swept configuration by `check_rep_fixes`, not assumed.
_BUILT = [(float(np.trace(R)),) + _ambient_action(R, corners(A_VE))
          for R in _ROTS]
_ACTIONS = [(t, T) for t, T, _ in _BUILT]
_BUILD_RESIDUAL = max(r for _, _, r in _BUILT)
# The Klein four-group V4 = {identity, the three C2}: the elements whose
# permutation part is trivial. Used ONLY as the teeth arm of the rank test.
_V4 = [(t, T) for (t, T) in _ACTIONS if t in (3.0, -1.0)]


def check_rep_fixes(a):
    """max |T x(a) - x(a)| over the group. Must be roundoff at EVERY a."""
    x = corners(a).reshape(-1)
    return max(float(np.abs(T @ x - x).max()) for _, T in _ACTIONS)


# --------------------------------------------------------------------------
# pushing the group down to the chart, and the character projectors
# --------------------------------------------------------------------------

def exact_chart_jacobian(F, origin_pivot=False):
    """d x / d q at q = 0 in CLOSED FORM, with no finite differences.

    The chart map is q -> Newton-project(x0 + B q) and its derivative at the
    origin is the ambient position Jacobian composed with B. Using the closed
    form here rather than the stencil's D matters: the group representation and
    hence the block LABELLING then carries no truncation at all, and the only
    finite-differenced ingredient left in the labelling is the metric g used to
    project. Measured agreement with the stencil's D is ~1e-08 at H_MAIN, which
    is the O(h^2) of the stencil and not an error in either.
    """
    if origin_pivot:
        return origin_position_jacobian(F.X0) @ F.B
    return position_jacobian(F.X0) @ F.B


class Characters:
    """The three irrep projectors on one chart, plus every check on them.

    THE PUSH-DOWN. A group element acts on the ambient configuration space by T.
    It maps the constrained variety to itself and fixes the base configuration,
    so it maps the tangent space to itself -- but a CHART is a section of that
    tangent space (6 internal directions chosen out of 12 internal + rigid), and
    nothing forces a section to be invariant. The honest push-down is therefore
    the projection ONTO the section ALONG the rigid orbit, in the mass inner
    product, which is exactly what the momentum-free form already does:

        rho(g) = gx^{-1} D^T Wh (T D)

    AND THE PROJECTION TURNS OUT TO REMOVE NOTHING HERE, which is worth saying
    plainly rather than leaving the construction to imply otherwise. The
    residual |T D - D rho| / |D| is measured at ~2e-15 in BOTH charts at every
    angle tried, so T D already lies in the section and the projection is a
    no-op on this path: jb_j's linear gauge happens to be group-equivariant.
    The projection is still the right construction -- nothing above guarantees
    that, and the number is what shows it -- but no claim in this file rests on
    it doing any work.

    EVERYTHING IN THE LABELLING IS BUILT FROM THE CLOSED-FORM D AND ITS OWN
    METRIC gx = D^T Wh D, not from the stencil's finite-differenced g. Mixing
    the two was measured to put the representation law, the isometry and the
    projector algebra at ~3e-09 -- which is exactly the O(h^2) gap between the
    two metrics, i.e. an artefact of the inconsistency and not of the
    construction. With a consistent pair every one of those residuals is
    roundoff, so the labelling carries NO truncation at all and the only
    finite-differenced object left in the pipeline is the Hessian itself and the
    metric it is diagonalised against.

    THE PROJECTORS. P_X = (d_X / 12) sum_g chi_X(g) rho(g) with the character
    table above. Three things are then CHECKED rather than assumed: the ranks
    are (1, 2, 3), the three projectors sum to the identity, and each is
    idempotent. The rank check is what a wrong character table or a missing
    group element breaks first.

    THE BLOCK BASES ARE ORTHONORMAL IN THE MEASURING METRIC. P_X is self-adjoint
    with respect to the metric, not symmetric as a matrix, so it is symmetrised
    by the Cholesky factor before being diagonalised: with g = L L^T, the
    operator L^T P L^-T is a genuine symmetric projector whose eigenvalues must
    be 0 or 1 (gated, in gx where the statement is exact). Its unit eigenvectors
    mapped back by L^-T give a basis Q with Q^T g Q = I, so that Q^T (Hess V) Q
    is directly the block of the GENERALISED eigenproblem and no second
    eigensolve is needed. The Q basis uses the STENCIL's g -- the same object
    the Hessian is built against -- while the diagnostics use gx, because a
    block basis orthonormal in a metric other than the measuring one would put
    an O(h^2) error into every block value. Using `eigh` of a symmetric matrix
    rather than an SVD also removes the sign ambiguity that would otherwise make
    the output non-reproducible between runs.
    """

    def __init__(self, D, W, g_metric, label, group=None):
        self.label = label
        self.g = g_metric
        gx = D.T @ W @ D
        self.gx = 0.5 * (gx + gx.T)
        grp = _ACTIONS if group is None else group
        self.rho = []
        self.leak = 0.0
        gi_DtW = np.linalg.solve(self.gx, D.T @ W)
        scale = float(np.abs(D).max())
        for t, T in grp:
            TD = T @ D
            r = gi_DtW @ TD
            self.rho.append((t, r))
            self.leak = max(self.leak,
                            float(np.abs(TD - D @ r).max()) / scale)
        self.iso = max(float(np.abs(r.T @ self.gx @ r - self.gx).max())
                       for _, r in self.rho) / float(np.abs(self.gx).max())
        n = len(grp)
        # UNKNOWN CONJUGACY CLASSES ARE COUNTED, NOT RAISED ON. `CHI` is keyed on
        # trace(R), which separates T's three classes exactly -- but only for
        # elements that ARE symmetries. A corrupted group (an improper sign
        # pattern, a wrong permutation) produces a trace outside {3, -1, 0} and
        # `CHI[b][t]` used to raise KeyError inside a constructor called before
        # anything is printed, so the gate table never appeared and no row could
        # go red for a broken group. The character is taken as 0 for such an
        # element and the count is surfaced: the projector ranks then fail, the
        # block values come back (nan, inf, inf), and V0's own row below says
        # what happened. A FAIL row, not a traceback.
        self.unknown = sum(1 for t, _ in grp if t not in CHI["A"])
        self.P = {b: sum(DIMW[b] * CHI[b].get(t, 0.0) * r for t, r in self.rho)
                  / n for b in BLOCKS}
        self.rank = {b: int(np.linalg.matrix_rank(self.P[b], tol=1e-8))
                     for b in BLOCKS}
        self.alg = float(np.abs(sum(self.P[b] for b in BLOCKS)
                                - np.eye(self.gx.shape[0])).max())
        for b in BLOCKS:
            self.alg = max(self.alg,
                           float(np.abs(self.P[b] @ self.P[b]
                                        - self.P[b]).max()))
        Lx = np.linalg.cholesky(self.gx)
        Lxi = np.linalg.inv(Lx).T
        L = np.linalg.cholesky(0.5 * (g_metric + g_metric.T))
        self.Q = {}
        self.ev_err = 0.0
        for b in BLOCKS:
            Px = Lx.T @ self.P[b] @ Lxi
            ev = np.linalg.eigvalsh(0.5 * (Px + Px.T))
            self.ev_err = max(self.ev_err,
                              float(np.abs(np.minimum(np.abs(ev),
                                                      np.abs(ev - 1.0))).max()))
            Pt = L.T @ self.P[b] @ np.linalg.inv(L).T
            Pt = 0.5 * (Pt + Pt.T)
            evm, EV = np.linalg.eigh(Pt)
            sel = EV[:, evm > 0.5]
            self.Q[b] = np.linalg.solve(L.T, sel)

    def mult_table(self):
        """max |rho(g1) rho(g2) - rho(g1 g2)| over all 144 pairs.

        The representation property is not automatic for a PROJECTED action --
        it holds because the projection is along an invariant complement. If the
        rigid orbit were being projected out the wrong way this residual would
        not be roundoff, so the check has content.
        """
        def key(M):
            return tuple(np.round(np.asarray(M).reshape(-1), 9))
        idx = {key(R): i for i, R in enumerate(_ROTS)}
        worst = 0.0
        for i, A in enumerate(_ROTS):
            for j, B in enumerate(_ROTS):
                k = idx.get(key(A @ B))
                if k is None:
                    # The set is not closed under composition, so it is not a
                    # group and there is no rho(g1 g2) to compare against. That
                    # IS the violation this method measures -- reported as an
                    # infinite residual, not as a KeyError, because a raise here
                    # would destroy the gate table that says which OTHER rows
                    # the same defect reddens.
                    return float("inf")
                worst = max(worst, float(np.abs(self.rho[i][1] @ self.rho[j][1]
                                                - self.rho[k][1]).max()))
        return worst

    def project_out(self, u):
        """|P_A u - u|_g / |u|_g: how far u is from the singlet block."""
        d = self.P["A"] @ u - u
        return float(np.sqrt(d @ self.g @ d) / np.sqrt(u @ self.g @ u))

    def block_values(self, Hess):
        """{block: (lambda, relative deviation, ABSOLUTE deviation)}.

        lambda = trace(Q^T Hess Q) / mult with Q g-orthonormal, i.e. the average
        of the block's generalised eigenvalues. The deviation is returned WITH
        it, never computed and dropped: symmetry forces Hess to be scalar on an
        irrep block, so the deviation is the falsifier for the whole
        construction, and a lambda quoted without it is an average presented as
        an eigenvalue.

        BOTH FORMS OF THE DEVIATION ARE RETURNED, and the third slot is not
        redundant. The RELATIVE form is the right falsifier almost everywhere --
        the deviation's floor scales with the block value, and an absolute
        threshold set from the well-behaved rows is exceeded by twenty decades
        near the pole. But it is STRUCTURALLY BLIND at a sign change: the
        denominator is |lambda| and a turnover is exactly where lambda -> 0, so
        the relative statistic carries no information about the measurements the
        bisection of V3c is made of. V3d uses the ABSOLUTE form there and
        converts it to an angular uncertainty with the measured d lambda / d a.
        Existing callers index [0] and [1] and are unaffected by the third slot.
        """
        Hs = 0.5 * (Hess + Hess.T)
        out = {}
        for b in BLOCKS:
            Q = self.Q[b]
            k = Q.shape[1]
            if k != MULT[b]:
                # The projector did not come back with the rank the character
                # table requires, so there is no block here to take a trace
                # over. Returned as (nan, inf) rather than raised: the rank row
                # in the gate is the place that says so, and a raise would
                # abort before the gate prints and hide which OTHER rows the
                # same defect reddens. A zero-width block also makes
                # `np.abs(...).max()` a zero-size reduction, which is an
                # exception with no diagnostic value at all.
                out[b] = (float("nan"), float("inf"), float("inf"))
                continue
            Hb = Q.T @ Hs @ Q
            lam = float(np.trace(Hb) / k)
            dev = float(np.abs(Hb - lam * np.eye(k)).max())
            out[b] = (lam, dev / max(abs(lam), 1e-300), dev)
        return out


# --------------------------------------------------------------------------
# one configuration: the chart, the geometry for both mass models, the
# characters, and every potential evaluated on the shared stencil
# --------------------------------------------------------------------------

class Site:
    """Everything measured at one angle a, for both mass models at once.

    The expensive part is the chart: 73 Newton solves per stencil. Every kernel,
    every primitive and both mass models read the SAME stencil points, which is
    exact sharing rather than an approximation (jb_o's `Probe` makes the same
    observation). The potential's value, gradient and chart Hessian are computed
    ONCE per (kernel, primitive) and the mass model enters only through
    Hess = Hnaive - Gamma . dV, so the 36 combinations cost one stencil and 18
    potential passes rather than 36 of each.
    """

    def __init__(self, a, h=H_MAIN, kind=None, origin_pivot=False):
        kind = METRIC_FORM if kind is None else kind
        self.a = float(a)
        self.h = float(h)
        self.kind = kind
        self.origin_pivot = origin_pivot
        F = OriginFrame(self.a) if origin_pivot else aligned_frame(self.a)
        self.F = F
        self.dim = F.dim
        if F.dim != 6:
            raise ChartUnmeasurable(
                f"chart dimension {F.dim} at a = {self.a}: this is a branch "
                f"point of the variety, not a 6-D chart, and the 1 + 2 + 3 "
                f"block structure is not defined on it.")
        self.chart = FrameChart(F, "origin pivot" if origin_pivot
                                else "centroid pivot")
        self.st = Stencil(self.chart, self.h)
        self.D = exact_chart_jacobian(F, origin_pivot)
        self.geo = {m: Geometry(self.st, forms(m, kind)) for m in MODELS}
        self.chars = {m: Characters(self.D, self.geo[m].W, self.geo[m].g,
                                    f"{self.chart.label}/{m}")
                      for m in MODELS}
        self._pot = {}

    def potential(self, kernel, primitive):
        """(V, dV, chart Hessian) for one (kernel, primitive), cached."""
        key = (kernel, primitive)
        if key not in self._pot:
            kern = dict(RAW_KERNELS)[kernel]
            fv = (vertex_potential(kern) if primitive == "vertex"
                  else strut_potential(kern))
            self._pot[key] = self.st.derivs(fv)
        return self._pot[key]

    def split_covector(self, dV, model):
        """(|dV|_euclid, path g-norm, transverse g-norm) of ANY covector.

        Q is g-orthonormal, so `Q^T dV` are the components of the vector
        g^{-1} dV in that block's basis and their norm is its g-norm. The FIRST
        entry is a Euclidean norm of the COVECTOR (the record's convention, and
        what V5 anchors against); the other two are g-norms of the VECTOR. They
        are norms of different objects in different metrics, so the "parts" are
        not bounded by the "whole" -- at the icosahedron the path part is 4.7x
        it. Only the RATIO of the last two is used in any claim.

        FACTORED OUT OF `gradient_split` SO THAT THE SPLIT ITSELF CAN BE GIVEN
        TEETH. V1b's transverse row is a `<=` row on a quantity derived to be
        exactly zero, so a build that simply never computes a transverse
        component satisfies it: measured, a mutation zeroing the transverse
        term exits 0 with every row green. V1b now runs this same routine on a
        covector whose gradient vector is a pure TRANSVERSE chart direction,
        where the transverse part must dominate.
        """
        ch = self.chars[model]
        c = {b: ch.Q[b].T @ dV for b in BLOCKS}
        return (float(np.linalg.norm(dV)), float(np.linalg.norm(c["A"])),
                float(np.sqrt(c["E"] @ c["E"] + c["F"] @ c["F"])))

    def gradient_split(self, kernel, primitive, model):
        """`split_covector` applied to this (kernel, primitive)'s own dV."""
        _, dV, _ = self.potential(kernel, primitive)
        return self.split_covector(dV, model)

    def blocks(self, kernel, primitive, model):
        """{block: (lambda, scalar deviation)} plus (V, |dV|)."""
        v0, dV, Hv = self.potential(kernel, primitive)
        G = self.geo[model]
        Hess = Hv - np.einsum('cab,c->ab', G.Gam, dV)
        return (self.chars[model].block_values(Hess), float(v0),
                float(np.linalg.norm(dV)))


_SITE_CACHE = {}


@jb_cache.memoize(_MODULE)
def _site_build(a, h, kind, origin_pivot):
    """The actual construction, behind the DISK layer of `site`'s two-level
    cache.

    `a` ARRIVES RAW, NOT ROUNDED, and that is load-bearing. `_SITE_CACHE`'s
    key rounds to 12 decimals, but the construction this file has always
    performed uses the caller's full-precision angle -- V3c's bisection and
    V9's step sweep both hand in angles well past the 12th decimal. Feeding
    the rounded key into `Site` instead was tried and REJECTED: it moved
    published roots in the 6th decimal (45.818000 -> 45.818001, 62.094048 ->
    62.094046) and their |lambda| columns by ~1e-6, which is a different
    measurement wearing the same row. The disk key is therefore the raw angle;
    the rounding stays where it always was, in the in-process layer above.

    Cheap here, expensive across a run: ~137ms each, ~1391 distinct arguments
    per gate (measured 2026-08-22), which is 195s of a 212s run. The entries
    are large (~822KB, ~360KB stored), so the store is compressed and worth
    watching; `--clear-cache` is the whole maintenance story."""
    return Site(a, h=h, kind=kind, origin_pivot=origin_pivot)


def site(a, h=H_MAIN, kind=None, origin_pivot=False):
    """`Site`, memoised. Purely an economy -- 73 Newton solves per stencil, and
    the bisections of V3c and the mirrored angles of V7 revisit the same
    configurations. Nothing about the result depends on the cache; it is keyed
    on every argument that enters the construction.

    TWO LEVELS. `_SITE_CACHE` is the in-process one this file always had, and
    it still absorbs the repeat calls within a single run (1604 calls, 1391
    distinct). Underneath it `_site_build` adds a disk layer keyed on the
    transitive SOURCE of the construction as well as its arguments, so the
    cost survives process exit but is discarded the moment anything `Site`
    reads actually changes."""
    kind = METRIC_FORM if kind is None else kind
    key = (round(float(a), 12), float(h), kind, bool(origin_pivot))
    if key not in _SITE_CACHE:
        # NOTE the raw `a`, `h` -- not `*key`. See `_site_build`'s docstring:
        # passing the rounded key moves published roots.
        _SITE_CACHE[key] = _site_build(a, h, kind, origin_pivot)
    return _SITE_CACHE[key]


COMBOS = [(kn, pn, m)
          for kn, _ in RAW_KERNELS
          for pn in ("vertex", "strut")
          for m in MODELS]
INVERSE_POWER = tuple(kn for kn, _ in RAW_KERNELS if kn.startswith("1/r"))
GAUSSIAN = tuple(kn for kn, _ in RAW_KERNELS if kn.startswith("gauss"))


def main_grid():
    """The swept angles. Deterministic, explicit, and printed with the results.

    Half a degree over the fundamental domain, plus every configuration of
    record. TWO neighbourhoods are NOT in this list and each has its own
    section, because in each of them a finite-difference second derivative
    stops resolving and saying so is better than reporting a number that has
    stopped meaning anything:

      * |a - 60| < POLE_GUARD -- the octahedron, where the twelve shared
        vertices merge in pairs so every inverse power diverges for BOTH
        primitives. Section V6.
      * a > 90 - BRANCH_GUARD -- the branch point, where the local dimension
        rises to 7 AND (measured in V7c, and not in the record) the strut
        midpoints collide, so every inverse power diverges for the STRUT
        primitive alone. Section V7.

    The first run of this file did NOT have the second guard, and the a = 89.99
    row it produced is what found the strut pole: 1/r^3 on struts reported a
    per-block scalarity of 2.7 -- a deviation nearly three times the block value
    it was deviating from -- which is a measurement announcing that it is not
    one. The guard was added because of that row, not before it.
    """
    xs = list(np.arange(0.0, 90.0, 0.5))
    xs += [A_ICO, A_2ND_THOMSON_VERTEX, A_MAX_GAUSS05_VERTEX,
           73.873521]                      # gauss 0.5 second minimum (record)
    xs = sorted({round(float(x), 9) for x in xs})
    return [x for x in xs
            if abs(x - A_POLE) >= POLE_GUARD and x <= A_BRANCH - BRANCH_GUARD]



# ==========================================================================
# V0 -- the symmetry machinery, checked before it is used for anything
# ==========================================================================

def v0_symmetry_controls(sites):
    print("=" * 78)
    print("V0  THE SYMMETRY GROUP, AND WHY THE BLOCK LABELS ARE NOT SORTED")
    print("=" * 78)
    print("  The symmetric path is the FIXED-POINT SET of the chiral")
    print("  tetrahedral group acting on the linkage, so at every point of the")
    print("  path that group acts on the 6-D internal tangent space and")
    print("  6 = 1 + 2 + 3 is forced. That is the decomposition this bead")
    print("  needs: the singlet is the path's own direction, the doublet and")
    print("  triplet are the five transverse ones.")
    print()
    print("  THE GROUP IS DERIVED, NOT TABULATED. Body s carries the twist")
    print("  sigma = sx sy sz (jb_a: forced chirality). A signed permutation")
    print("  with permutation p and signs e sends sigma to (prod e) sigma and")
    print("  has determinant sign(p)(prod e), so a rotation preserving the")
    print("  jitterbug needs prod e = +1 AND sign(p) = +1: three even")
    print("  permutations times four even sign patterns.")
    print(f"      |G| = {len(_ROTS)}   (the chiral tetrahedral group, order 12 --")
    print("      NOT the order-24 chiral octahedral group of the cuboctahedron,")
    print("      because the twist breaks it)")
    cls = {}
    for t, _ in _ACTIONS:
        cls[t] = cls.get(t, 0) + 1
    print(f"      conjugacy classes by trace(R): "
          + ", ".join(f"trace {t:+.0f}: {n} element(s)"
                      for t, n in sorted(cls.items(), reverse=True)))
    print("      1 + 3 + 8 = 12, the class structure of T. trace separates the")
    print("      classes exactly, which is what the character table is keyed on.")
    # ARE THEY ACTUALLY ROTATIONS, AND IS THE CLASS STRUCTURE THE ONE CLAIMED?
    # Asserted rather than left to the derivation in the docstring, because the
    # derivation is what a mutation edits. An improper sign pattern (prod e = -1)
    # gives det = -1 and a trace outside {3, -1, 0}; before this row existed such
    # a mutation raised KeyError inside `Characters` and the gate never printed.
    det_err = max(abs(float(np.linalg.det(R)) - 1.0) for R in _ROTS)
    orth_err = max(float(np.abs(R.T @ R - np.eye(3)).max()) for R in _ROTS)
    classes_ok = (len(_ROTS) == 12 and cls == {3.0: 1, -1.0: 3, 0.0: 8})
    print(f"      max |det(R) - 1|            = {det_err:.3e}"
          f"   (criterion <= {TOL['rep_mult']:.0e})")
    print(f"      max |R^T R - I|             = {orth_err:.3e}"
          f"   (criterion <= {TOL['rep_mult']:.0e})")
    print(f"      class structure is {{3:1, -1:3, 0:8}} and |G| = 12: "
          f"{classes_ok}")
    print("      Those three make the group ITSELF falsifiable. A matrix that is")
    print("      not a proper rotation, or a class count that is not T's, is now")
    print("      a red row rather than a KeyError before the table prints.")

    print()
    print("  (a) DOES THE GROUP ACTUALLY FIX THE CONFIGURATION? The body and")
    print("      corner permutations are found by MATCHING at a = 0, so a wrong")
    print("      labelling convention cannot hide in them -- and they are then")
    print("      REUSED at every angle, which is a claim (the permutations do")
    print("      not depend on a) and is therefore checked at every angle.")
    print("      THE MATCH IS TAKEN AT a = 0 DELIBERATELY AND THAT IS THE ONE")
    print("      CONFIGURATION IN THE SWEEP WITH MORE SYMMETRY THAN THE")
    print("      JITTERBUG (the cuboctahedron, order 24 against T's 12), so an")
    print("      argmin match could in principle select a permutation valid only")
    print("      there. That is MITIGATED rather than argued away: the fix")
    print("      residual below is evaluated at every SWEPT angle (where the")
    print("      symmetry is T and nothing larger), and the multiplication table")
    print("      in (b) returns an infinite residual if the matched set is not")
    print("      closed under composition. Either would catch it.")
    print(f"      build residual at a = 0, worst over 12 elements = "
          f"{_BUILD_RESIDUAL:.3e}")
    fix = 0.0
    fix_a = None
    swept = [s.a for s in sites]
    # EVERY SWEPT ANGLE, not a hand-picked list. The previous version checked
    # nine hardcoded angles, FOUR of which (A_ICO aside, the 59.5 / 60.0 / 60.5 /
    # 89.99 entries) are not swept angles at all, while the gate row above it
    # said "at every a". The claim and the loop now agree.
    for a in swept:
        r = check_rep_fixes(a)
        # `fix_a is None or` -- the same sentinel class this pass is removing
        # from V1/V1b/V3b, and these two sites were INTRODUCED by this pass. A
        # group whose fix residual were exactly 0.0 at every angle would leave
        # `fix_a` None and kill the print. Found by grepping this pass's own
        # diff for the shape, which is the discipline the class demands.
        if fix_a is None or r > fix:
            fix, fix_a = r, a
    # ... and then, separately labelled, the angles the MAIN GRID excludes,
    # because the guard bands are where a reader will ask whether the symmetry
    # machinery still holds.
    guard = 0.0
    guard_a = None
    for a in (59.5, 59.9, 59.99, A_POLE, 60.01, 60.1, 60.5, 89.5, 89.9, 89.99,
              A_BRANCH):
        r = check_rep_fixes(a)
        if guard_a is None or r > guard:
            guard, guard_a = r, a
    print(f"      {'set':38s} {'angles':>8s} {'worst |T x(a) - x(a)|':>24s}")
    print(f"      {'EVERY SWEPT ANGLE (the V3 grid)':38s} {len(swept):8d} "
          f"{fix:24.3e}")
    print(f"      {'  worst at a =':38s} {fix_a:8.4f}")
    print(f"      {'the GUARD-BAND angles, incl. a = 60, 90':38s} "
          f"{11:8d} {guard:24.3e}")
    print(f"      {'  worst at a =':38s} {guard_a:8.4f}")
    fix = max(fix, guard)
    print(f"      worst over BOTH sets = {fix:.3e}"
          f"   (criterion <= {TOL['rep_fix']:.0e})")
    print("      NOTE a = 60 AND a = 90 ARE IN THE SECOND SET. The group action")
    print("      is exact at both; it is the MASS METRIC that degenerates at the")
    print("      octahedron (V2) and the CHART that gives out at the branch")
    print("      point (V7a), not the symmetry.")

    print()
    print("  (a-bis) AND IS V ITSELF GROUP-INVARIANT? That is the PREMISE of")
    print("      V1b's derivation -- 'dV is an invariant covector at a point the")
    print("      group fixes' is a statement about V, not about x -- and it was")
    print("      previously DERIVED (a radial kernel over a permuted pair set)")
    print("      and confirmed only by its consequence (V1b's transverse")
    print("      gradient). The premise itself is one line, so it is measured.")
    print()
    print("      THE OBVIOUS TEST IS VACUOUS AND THIS FILE SHIPPED IT ONCE.")
    print("      Comparing V(Tx) against V(x) AT A CONFIGURATION ON THE PATH")
    print("      measures nothing whatever: row (a) above asserts T x = x to")
    print("      2e-16, so V(Tx) = V(x) follows for ANY function of x, invariant")
    print("      or not. Measured: adding a term LINEAR IN A SINGLE COORDINATE")
    print("      to the potential -- as flagrantly non-invariant a perturbation")
    print("      as exists -- left the residual at 0.000e+00 and the row green.")
    print("      The test has to be taken OFF the fixed-point set, where Ty and")
    print("      y are genuinely different configurations -- AND STILL ON THE")
    print("      VARIETY, which the first attempt at this row got wrong and")
    print("      which is a fact about V worth recording. A free ambient")
    print("      displacement measured 6.4e-02, i.e. flagrant non-invariance,")
    print("      and it is not a bug: `verts_of` selects ONE representative")
    print("      corner per shared-vertex pair, the two coincide only where the")
    print("      linkage constraints hold, and the group permutes corners. So")
    print("      V IS G-INVARIANT ON THE CONSTRAINT VARIETY AND NOT ON THE")
    print("      AMBIENT SPACE, which is exactly the statement V1b's derivation")
    print("      needs (dV there is the differential of V ON the variety).")
    print("      y is therefore taken as a CHART displacement -- Newton-projected")
    print("      onto the variety, transverse to the symmetric path, at a")
    print("      seeded pseudo-random direction of norm 0.1 in chart units.")
    print(f"      {'kernel':20s} {'prim':7s} {'a':>9s} "
          f"{'max_g |V(Ty) - V(y)| / |V(y)|':>30s}   |y - x|")
    vinv = 0.0
    rng = np.random.default_rng(20260816)
    for kn, pn in (("1/r^1  (Thomson)", "vertex"), ("1/r^12", "strut"),
                   ("gauss s=1.0", "vertex"), ("gauss s=2.5", "strut")):
        kern = dict(RAW_KERNELS)[kn]
        fv = (vertex_potential(kern) if pn == "vertex"
              else strut_potential(kern))
        for a in (0.0, A_ICO, 45.0, 89.0):
            ch = site(a).chart
            x = ch.x(np.zeros(6))
            q = rng.standard_normal(6)
            q = 0.1 * q / np.linalg.norm(q)
            y = ch.x(q)
            v = float(fv(y))
            r = max(abs(float(fv(T @ y)) - v) for _, T in _ACTIONS) / abs(v)
            vinv = max(vinv, r)
            print(f"      {kn:20s} {pn:7s} {a:9.4f} {r:30.3e}"
                  f"   {float(np.linalg.norm(y - x)):.4f}")
    print(f"      worst = {vinv:.3e}   (criterion <= {TOL['rep_iso']:.0e})")
    print("      No metric, no stencil derivative: the potential evaluated at y")
    print("      and at Ty for all twelve group elements, at four configurations")
    print("      per kernel, each O(0.1) off the symmetric path. This is the")
    print("      teeth for the critical-manifold derivation of V1b, which is the")
    print("      reason a transverse SIGN is worth reporting at all. Gated at")
    print("      1e-9 rather than at roundoff because the cancellation in")
    print("      |V(Ty) - V(y)| is between numbers of size |V|, and 1/r^12 near")
    print("      a close pair makes those large.")

    print()
    print("  (b) THE PUSHED-DOWN REPRESENTATION on the 6-D chart. rho(g) =")
    print("      g^-1 D^T Wh (T D): the ambient action projected onto the")
    print("      chart's section ALONG the rigid orbit, in the mass inner")
    print("      product. Three properties are checked, none assumed.")
    s0 = site(A_ICO)
    ch = s0.chars["point"]
    mt = ch.mult_table()
    print(f"      section leak |T D - D rho| / |D|      = {ch.leak:.3e}")
    print("      THAT NUMBER SAYS THE PROJECTION REMOVES NOTHING. T D already")
    print("      lies in the chart's section, in both charts and at every angle")
    print("      tried, because jb_j's linear gauge happens to be")
    print("      group-equivariant. Nothing above guaranteed that and no claim")
    print("      here rests on it -- the projection is kept because it is the")
    print("      construction that would be needed if the gauge were not.")
    print(f"      representation law, worst of 144 pairs = {mt:.3e}"
          f"   (criterion <= {TOL['rep_mult']:.0e})")
    print(f"      isometry rho^T gx rho = gx, relative   = {ch.iso:.3e}"
          f"   (criterion <= {TOL['rep_iso']:.0e})")
    print("      The isometry row has content: the symmetry preserves the")
    print("      kinetic energy, so it must preserve the mass metric. A")
    print("      projection along the wrong complement would break it while")
    print("      leaving the multiplication table intact.")

    print()
    print("  (c) THE PROJECTORS. P_X = (d_X/12) sum_g chi_X(g) rho(g).")
    print(f"      {'block':8s} {'rank':>6s} {'expected':>9s}")
    for b in BLOCKS:
        print(f"      {b:8s} {ch.rank[b]:6d} {MULT[b]:9d}   {LONG[b]}")
    print(f"      |sum P - I| and |P^2 - P|, worst = {ch.alg:.3e}"
          f"   (criterion <= {TOL['proj_alg']:.0e})")
    print(f"      projector eigenvalues off {{0, 1}} by  = {ch.ev_err:.3e}"
          f"   (criterion <= {TOL['proj_ev']:.0e})")

    print()
    print("  (d) THE TEETH OF THE RANK TEST, because a rank check that cannot")
    print("      fail is not a check. Restrict to the KLEIN SUBGROUP V4 (the")
    print("      identity and the three C2). chi_E restricted to V4 is the")
    print("      trivial character TWICE, so the V4-fixed subspace is A + E and")
    print("      the 'singlet' projector built from V4 has rank 3, not 1. If")
    print("      the machinery were insensitive to which group it is given,")
    print("      this number would come back 1.")
    v4 = Characters(s0.D, s0.geo["point"].W, s0.geo["point"].g, "V4",
                    group=_V4)
    print(f"      rank of the V4 'A' projector = {v4.rank['A']}"
          f"   (criterion == {TOL['v4_rank']}; the full group gives "
          f"{ch.rank['A']})")
    grp_ok = (det_err <= TOL["rep_mult"] and orth_err <= TOL["rep_mult"]
              and classes_ok and ch.unknown == 0)
    ok = (fix <= TOL["rep_fix"] and _BUILD_RESIDUAL <= TOL["rep_fix"]
          and vinv <= TOL["rep_iso"] and grp_ok
          and mt <= TOL["rep_mult"] and ch.iso <= TOL["rep_iso"]
          and ch.alg <= TOL["proj_alg"] and ch.ev_err <= TOL["proj_ev"]
          and all(ch.rank[b] == MULT[b] for b in BLOCKS)
          and v4.rank["A"] == TOL["v4_rank"])
    print(f"\n  V0 PASSED: {ok}")
    return dict(fix=fix, build=_BUILD_RESIDUAL, mult=mt, iso=ch.iso,
                alg=ch.alg, ev=ch.ev_err, vinv=vinv, nswept=len(swept),
                ranks=tuple(ch.rank[b] for b in BLOCKS),
                v4rank=v4.rank["A"], grp_ok=grp_ok, det=det_err,
                orth=orth_err, unknown=ch.unknown, ok=ok)


# ==========================================================================
# V1 -- IS THE SINGLET THE PATH TANGENT? Verified at every swept a.
# ==========================================================================

def v1_singlet_is_the_path_tangent(sites):
    print()
    print("=" * 78)
    print("V1  IS THE SINGLET THE PATH TANGENT? -- at EVERY swept angle")
    print("=" * 78)
    print("  The record has overlap^2 = 1.000000000 AT THE VE and nowhere")
    print("  else. The seed's instruction is to verify it along the sweep and")
    print("  not to assume it, because if it failed anywhere the transverse")
    print("  subspace would not be D + T there and the whole decomposition")
    print("  would be measuring something other than what it is called.")
    print()
    print("  THE TEST IS CHEAP BECAUSE OF HOW THE CHART IS BUILT.")
    print("  jb_k's `aligned_frame` pins basis direction 0 to the symmetric-")
    print("  path tangent (it asserts |B[:,0] . xi| = 1 to 1e-9), so the")
    print("  question 'is the singlet the path tangent' is exactly 'does the A")
    print("  projector fix the chart's zeroth basis vector'. No ambient")
    print("  re-projection of the tangent is needed and none is done.")
    e0 = np.zeros(6)
    e0[0] = 1.0
    worst = 0.0
    worst_a = None
    teeth = np.inf
    for s in sites:
        for m in MODELS:
            r = s.chars[m].project_out(e0)
            # `worst_a is None or` IS THE FIX FOR A WHOLE TRACEBACK CLASS, not a
            # defensive flourish. With a strict `r > worst` alone the companion
            # sentinel is assigned only when some row beats 0.0, so a mutation
            # that drives the statistic to EXACTLY zero -- which is precisely the
            # falsifier-silencing mutation this row exists to catch -- leaves it
            # None and the print below dies on TypeError BEFORE the gate table
            # is emitted. An independent mutation audit found seven tracebacks
            # in this file with one root cause and three sites, and every one of
            # them was reached by a falsifier-silencing probe. A raise inside or
            # before a swept section destroys the verdict table, which is gate
            # lesson 4 recorded for M2/M3 and not generalised at the time.
            if worst_a is None or r > worst:
                worst, worst_a = r, (s.a, m)
            for i in range(1, 6):
                u = np.zeros(6)
                u[i] = 1.0
                teeth = min(teeth, s.chars[m].project_out(u))
    print(f"  |P_A e0 - e0|_g / |e0|_g, worst over "
          f"{len(sites)} angles x {len(MODELS)} mass models")
    print(f"      = {worst:.3e} at a = {worst_a[0]:.6f}, {worst_a[1]}"
          f"   (criterion <= {TOL['tangent_leak']:.0e})")
    print(f"  THE SAME STATISTIC on the five TRANSVERSE chart directions,")
    print(f"  smallest over all of them = {teeth:.3e}"
          f"   (criterion >= {TOL['tangent_teeth']:.2f})")
    print("  That second row is what makes the first one a measurement. The")
    print("  statistic is ~1e-09 on the path tangent and exactly 1 on a")
    print("  direction the projector annihilates, so the test separates the two")
    print("  cases by nine decades and cannot be passed by an inert projector.")
    print()
    print("  VERDICT: the singlet block IS the symmetric path's tangent at")
    print("  every angle swept, both mass models. The five transverse")
    print("  directions are therefore exactly the DOUBLET and the TRIPLET, and")
    print("  'transverse curvature' below means precisely those five.")
    print("  SCOPE: 'every angle swept' is the grid printed in V3, not the")
    print("  continuum, and it excludes the pole band and a >= 90.")
    ok = worst <= TOL["tangent_leak"] and teeth >= TOL["tangent_teeth"]
    print(f"\n  V1 PASSED: {ok}")
    return worst, teeth, ok


def v1b_gradient_has_no_transverse_part(sites):
    print()
    print("=" * 78)
    print("V1b  THE GRADIENT HAS NO TRANSVERSE COMPONENT, ANYWHERE ON THE PATH")
    print("=" * 78)
    print("  DERIVED FIRST. dV is an invariant covector at a point the group")
    print("  fixes, so it lives in the A-isotypic component of the tangent")
    print("  space -- and V0/V1 measure that component to be one-dimensional and")
    print("  to be the path tangent. So dV can have NO doublet or triplet part,")
    print("  at any angle. Equivalently: the symmetric path is the symmetry's")
    print("  fixed-point set, so the force on it has no component that would")
    print("  take a trajectory off it, which is the dynamical-invariance fact")
    print("  the memo already records -- restated as a gradient statement and")
    print("  then MEASURED rather than assumed.")
    print()
    print("  The split is invariant, not a choice of chart directions: with Q")
    print("  g-orthonormal on a block, Q^T dV are the components of the gradient")
    print("  VECTOR g^-1 dV in that block and their norm is its g-norm.")
    print()
    print("  READ THE TWO NORM COLUMNS AS DIFFERENT OBJECTS, because they are,")
    print("  and an earlier version of this table invited the opposite reading.")
    print("  |dV|_euclid is the EUCLIDEAN norm of the COVECTOR dV in chart")
    print("  coordinates -- it is the record's own convention and is what V5")
    print("  anchors against. The path and transverse columns are g-NORMS of the")
    print("  gradient VECTOR g^-1 dV. They are norms of different objects in")
    print("  different metrics, so the 'part' is NOT bounded by the 'whole' and")
    print("  at the icosahedron the path part is about 4.7x |dV|_euclid. Only the")
    print("  RATIO of the last two columns is used anywhere, and both of its")
    print("  terms are g-norms, so the measurement is unaffected -- what is")
    print("  corrected here is a presentation that let two conventions sit in one")
    print("  row under one name. |grad|_g = sqrt(path^2 + transverse^2) is")
    print("  printed so the decomposition can be seen to be one.")
    print()
    print("  FOUR DECLARATIONS on this table: KERNEL as printed per row,")
    print("  PRIMITIVE as printed per row, MASS MODEL point, METRIC FORM")
    print("  momentum-free. The swept statistic below the table covers all 36")
    print("  combinations and therefore both mass models and both primitives.")
    print(f"  {'a':>9s} {'kernel':20s} {'prim':7s} {'|dV|_euclid':>12s} "
          f"{'|grad|_g':>13s} {'path_g':>13s} {'transverse_g':>13s} "
          f"{'ratio':>11s}  note")
    worst_r = 0.0
    worst_where = None
    best_path = 0.0
    crit_rows = 0
    crit_abs = 0.0
    crit_path = 0.0
    for s_ in sites:
        for (kn, pn, m) in COMBOS:
            _, ga, gt = s_.gradient_split(kn, pn, m)
            best_path = max(best_path, ga)
            if ga < TOL["grad_crit_floor"]:
                crit_rows += 1
                crit_abs = max(crit_abs, gt)
                crit_path = max(crit_path, ga)
                continue
            # `worst_where is None or` -- same traceback class as V1, second of
            # three sites. Three separate probes (transverse component := 0,
            # both gradient parts := 0, a corrupted ambient action) died here
            # before the gate table printed.
            if worst_where is None or gt / ga > worst_r:
                worst_r, worst_where = gt / ga, (s_.a, kn, pn, m)
    for a in (0.0, A_ICO, 45.0, A_2ND_THOMSON_VERTEX, 85.0, 89.0):
        s_ = site(a)
        for kn, pn in (("1/r^1  (Thomson)", "vertex"),
                       ("gauss s=1.0", "vertex")):
            gv, ga, gt = s_.gradient_split(kn, pn, "point")
            note = ("CRITICAL POINT: ratio excluded, see below"
                    if ga < TOL["grad_crit_floor"] else "")
            print(f"  {a:9.4f} {kn:20s} {pn:7s} {gv:12.4e} "
                  f"{np.hypot(ga, gt):13.4e} {ga:13.4e} "
                  f"{gt:13.4e} {gt / ga if ga > 0 else float('nan'):11.3e}"
                  f"  {note}")
    print(f"\n  worst TRANSVERSE / PATH ratio over {len(sites)} angles x "
          f"{len(COMBOS)} combinations")
    if worst_where is None:
        # THE `is None` GUARD ON THE UPDATE IS NOT ENOUGH HERE and that is worth
        # spelling out, because the same guard IS enough in V1 and V3b. The
        # update in this loop sits behind a `continue` that skips every row whose
        # PATH gradient is below the critical floor -- so a mutation that
        # annihilates the path gradient (an inert A projector, a zeroed gradient
        # split) skips every row, leaves the sentinel None, and used to kill the
        # print before the gate table existed. The PRINT is guarded, and the
        # statistic goes to infinity so the row goes RED rather than absent: a
        # build in which every row is a critical point has not satisfied this
        # criterion, it has destroyed the measurement.
        print("      NO ROW HAD A PATH GRADIENT ABOVE "
              f"{TOL['grad_crit_floor']:.0e}, so the ratio was never")
        print("      formed. That is a destroyed measurement, not a passing one.")
        worst_r = float("inf")
    else:
        print(f"      = {worst_r:.3e} at a = {worst_where[0]:.4f}, "
              f"{worst_where[1]}, {worst_where[2]}, {worst_where[3]}"
              f"   (criterion <= {TOL['grad_transverse']:.0e})")
    print(f"  largest PATH gradient over the same set = {best_path:.3e}"
          f"   (criterion >= {TOL['grad_path_teeth']:.0e})")
    print("      That second row is the teeth. Without it the first would be")
    print("      satisfied by a build whose gradient vanished everywhere, which")
    print("      is the mistake of reading 'the correction is small' as 'the")
    print("      correction is right'.")
    s_ico = site(A_ICO)
    _, anchor_path, _ = s_ico.gradient_split(
        "1/r^1  (Thomson)", "vertex", "point")
    # THE TEETH FOR THE SPLIT ITSELF, and it was missing. The row above is a
    # `<=` on a quantity DERIVED to be exactly zero, so a build that never
    # computes a transverse component at all satisfies it -- measured: zeroing
    # the transverse term exits 0 with every gate row green. Feed the SAME
    # routine a covector whose gradient vector is a pure transverse chart
    # direction (u = g e1, so g^-1 u = e1) and the transverse part must
    # dominate. This is V1's tangent teeth restated in the split's own language.
    e1 = np.zeros(6)
    e1[1] = 1.0
    _, tp, tt = s_ico.split_covector(s_ico.geo["point"].g @ e1, "point")
    split_teeth = tt / tp if tp > 0 else float("inf")
    print(f"  ... AND THE SAME TEETH AT ONE FIXED WELL-BEHAVED ROW: the path")
    print(f"  gradient at the icosahedron (1/r^1 Thomson, raw vertex, point,")
    print(f"  momentum-free) = {anchor_path:.3e}"
          f"   (criterion >= {TOL['grad_path_teeth']:.0e})")
    print("      The MAX above is attained on a 1/r^12 row adjacent to the pole,")
    print("      where the gradient is 2.4e+24, so on its own it would be")
    print("      satisfied by a build whose gradient vanished on every")
    print("      well-behaved row and survived only next to a divergence. This")
    print("      row is at a configuration of record, far from both poles, and")
    print("      is the one a reader should read the teeth from.")
    print(f"  ... AND THE TEETH FOR THE SPLIT ITSELF: the same routine run on a")
    print(f"  covector whose gradient vector is the transverse chart direction")
    print(f"  e1 gives TRANSVERSE / PATH = {split_teeth:.3e}"
          f"   (criterion >= {TOL['split_teeth']:.0e})")
    print("      Without it the row above is satisfied by a build that never")
    print("      computes a transverse component at all -- measured, zeroing")
    print("      the transverse term exits 0 with every other row green,")
    print("      because the quantity it reports is DERIVED to be zero and a")
    print("      constant zero is indistinguishable from the truth. The teeth")
    print("      is the same statistic on an input where it must be large.")
    print()
    print(f"  {crit_rows} of {len(sites) * len(COMBOS)} rows have a PATH")
    print(f"  gradient below {TOL['grad_crit_floor']:.0e} and are EXCLUDED from")
    print("  the ratio: at a critical point the whole gradient vanishes and")
    print("  there is no direction for the transverse part to be small relative")
    print("  TO. They are reported by absolute size instead --")
    print(f"      worst path part {crit_path:.3e}, worst transverse part "
          f"{crit_abs:.3e}")
    print("  -- so nothing is dropped, and those rows are the same critical")
    print("  points V5 gates against the record.")
    print()
    print("  THE RATIO IS O(h^2) TRUNCATION, swept here rather than asserted:")
    print(f"      {'h':>9s} {'path part':>15s} {'transverse':>14s} "
          f"{'ratio':>11s}")
    if worst_where is None:
        print("      NOT MEASURABLE: no row formed a ratio (see above).")
    else:
        kw, pw, mw = worst_where[1], worst_where[2], worst_where[3]
        for hh in (1e-3, 3e-4, 1e-4):
            _, ga, gt = site(worst_where[0], h=hh).gradient_split(kw, pw, mw)
            print(f"      {hh:9.0e} {ga:15.6e} {gt:14.4e} {gt / ga:11.3e}")
    print("      Ten-to-one per factor 3.16 in h, on the worst row in the")
    print("      sweep. So the threshold is a statement about H_MAIN and not")
    print("      about the construction, and the quantity itself is zero.")
    print()
    print("  WHAT THIS BUYS, and it is more than a consistency check. It means")
    print("  the symmetric path is a CRITICAL MANIFOLD OF V IN THE TRANSVERSE")
    print("  DIRECTIONS: at every angle, V is stationary under a transverse")
    print("  displacement at fixed a. So the doublet and triplet blocks of")
    print("  (Hess V, g) are the SECOND VARIATION OF V NORMAL TO AN INVARIANT")
    print("  SUBMANIFOLD, not generic off-critical curvature -- which is why a")
    print("  sign is worth reporting at all. It does NOT make them frequencies:")
    print("  the PATH component of the gradient is emphatically not zero (the")
    print("  column above), so the system accelerates ALONG the path and the")
    print("  transverse motion is driven by a time-dependent coefficient. The")
    print("  transverse Hessian is the potential term of that normal")
    print("  variational equation and not the equation.")
    ok = (worst_r <= TOL["grad_transverse"]
          and best_path >= TOL["grad_path_teeth"]
          and anchor_path >= TOL["grad_path_teeth"]
          and split_teeth >= TOL["split_teeth"])
    print(f"\n  V1b PASSED: {ok}")
    return worst_r, best_path, anchor_path, split_teeth, ok


# ==========================================================================
# V2 -- the labelling that does NOT survive this sweep, and why
# ==========================================================================

def _g_closed_form(a, model="point", origin_pivot=False):
    """The momentum-free mass metric in closed form: no finite differences."""
    F = OriginFrame(a) if origin_pivot else aligned_frame(a)
    D = exact_chart_jacobian(F, origin_pivot)
    W = WEIGHTS[model]
    Z = rigid_fields(F.X0.reshape(-1))
    WZ = W @ Z
    Wh = W - WZ @ np.linalg.solve(Z.T @ WZ, WZ.T)
    g = D.T @ Wh @ D
    return 0.5 * (g + g.T), D, Z, W


def v2_why_not_the_mass_metric():
    print()
    print("=" * 78)
    print("V2  WHY NOT jb_u's LABELLING -- the mass metric degenerates at a=60")
    print("=" * 78)
    print("  jb_u labels blocks by the eigenspaces of the mass metric, which is")
    print("  kernel-independent and therefore immune to the D/T ordering flips")
    print("  that defeat sort-position labelling. It is strictly better than")
    print("  sorting. It is still not enough here, and the reason is exact")
    print("  rather than numerical.")
    print()
    print("  CLOSED FORM -- no stencil, no truncation. g = D^T Wh D with D the")
    print("  analytic position Jacobian composed with the chart basis:")
    print(f"      {'a':>10s} {'eigenvalues of g (point mass, momentum-free)':>48s}")
    for a in (55.0, 59.0, 59.9, 60.0, 60.1, 61.0, 65.0):
        ev = np.linalg.eigvalsh(_g_closed_form(a)[0])
        print(f"      {a:10.4f}   {np.array2string(ev, precision=12)}")
    ev60 = np.linalg.eigvalsh(_g_closed_form(A_POLE)[0])
    d5 = float(np.abs(ev60[:5] - PRIOR["m_doublet"]).max())
    d1 = abs(float(ev60[5]) - PRIOR["m_singlet_at_60"])
    print()
    print("  AT a = 60 EXACTLY the spectrum is [1/32 x 5, 5/96]:")
    print(f"      max |ev[0..4] - 1/32|  = {d5:.3e}")
    print(f"      |ev[5] - 5/96|         = {d1:.3e}")
    print("  The DOUBLET mass is 1/32 at EVERY a (it is the flat one); the")
    print("  TRIPLET mass falls through it and the two are equal AT THE")
    print("  OCTAHEDRON. So in a band around a = 60 the mass metric has a")
    print("  five-fold eigenvalue and cannot separate the doublet from the")
    print("  triplet at all. This is not a tolerance question: at a = 60 the")
    print("  degeneracy is EXACT.")
    print()
    print("  AND jb_u's ROUTINE DOES THE RIGHT THING -- it refuses:")
    fails = []
    for a in (59.0, 59.5, 59.9, 60.0, 60.5, 61.0, 61.5):
        try:
            s = site(a)
            b, _, _ = s.blocks("gauss s=1.0", "vertex", "point")
            G = s.geo["point"]
            _, dV, Hv = s.potential("gauss s=1.0", "vertex")
            blocks_by_irrep(Hv - np.einsum('cab,c->ab', G.Gam, dV), G.g)
            verdict = "labels D/T/S"
        except IrrepLabelError:
            verdict = "REFUSES (IrrepLabelError)"
            fails.append(a)
        except ChartUnmeasurable as exc:
            verdict = f"chart unmeasurable: {str(exc).splitlines()[0][:40]}"
        print(f"      a = {a:8.3f}   blocks_by_irrep {verdict}")
    print(f"      -> the mass-metric labelling is UNDEFINED for a in roughly "
          f"[{min(fails):.1f}, {max(fails):.1f}]" if fails else
          "      -> no refusal seen in this list")
    print()
    print("  THE REPLACEMENT is the character labelling of V0, which depends on")
    print("  no eigenvalue ordering and no eigenvalue SEPARATION -- only on the")
    print("  group. Its ranks at a = 60, where the mass metric has nothing to")
    print("  say:")
    s60 = site(A_POLE)
    r60 = s60.chars["point"].rank
    print(f"      character projector ranks at a = 60: "
          + ", ".join(f"{b} = {r60[b]}" for b in BLOCKS))
    print("  CONSEQUENCE FOR THE RECORD, and it is a limitation of an existing")
    print("  tool rather than an error in it: any future work that needs D/T/S")
    print("  attribution near the octahedron cannot get it from the mass")
    print("  metric, and this file's `Characters` is the route.")
    ok = (d5 <= 1e-14 and d1 <= 1e-14
          and all(r60[b] == MULT[b] for b in BLOCKS) and bool(fails))
    print(f"\n  V2 PASSED: {ok}")
    return d5, d1, tuple(r60[b] for b in BLOCKS), fails, ok


# ==========================================================================
# V3 -- THE DELIVERABLE: the curvature spectrum along the path
# ==========================================================================

def v3_sweep(sites):
    print()
    print("=" * 78)
    print("V3  THE CURVATURE SPECTRUM ALONG THE SYMMETRIC PATH")
    print("=" * 78)
    print("  READ THE UNITS BEFORE READING THE NUMBERS. Except at the critical")
    print("  points listed in V5, |dV| is NOT zero on this path, so these are")
    print("  NOT frequencies and nothing oscillates. They are chart-invariant")
    print("  LOCAL CURVATURE SCALES of the pair (V, g). Their ABSOLUTE size")
    print("  carries the arc's convention (coupling 1, total mass 1/2, R = 1)")
    print("  and is not a measurement; their SIGN is, and so are their RATIOS.")
    print()
    print("  FOUR DECLARATIONS on the headline table: KERNEL 1/r^1 (Thomson)")
    print("  RAW, PRIMITIVE raw VERTEX, MASS MODEL point, METRIC FORM")
    print("  momentum-free. The full 36-combination sweep is below it.")
    print()
    print(f"  {'a':>9s} {'V':>13s} {'|dV|':>12s} {'S = path':>14s} "
          f"{'D transverse':>14s} {'T transverse':>14s} {'sgn(D,T)':>9s} "
          f"{'T/D':>9s} {'S/D':>10s} {'dev':>9s}")
    head = []
    for s in sites:
        b, v0, gd = s.blocks("1/r^1  (Thomson)", "vertex", "point")
        lam = {k: b[k][0] for k in BLOCKS}
        dev = max(b[k][1] for k in BLOCKS)
        sg = _sgn(lam["E"]) + _sgn(lam["F"])
        head.append((s.a, v0, gd, lam, dev))
        print(f"  {s.a:9.4f} {v0:13.6f} {gd:12.4e} {lam['A']:14.6e} "
              f"{lam['E']:14.6e} {lam['F']:14.6e} {sg:>9s} "
              f"{lam['F'] / lam['E']:9.5f} {lam['A'] / lam['E']:10.5f} "
              f"{dev:9.1e}")
    print("  'dev' is max |Hess|block - lambda I| / |lambda|, the falsifier for")
    print("  reading an AVERAGE as a block eigenvalue. It is pure O(h^2)")
    print("  truncation -- V9 shows it falling by ~10x per 3.16x drop in h --")
    print("  and it is gated over ALL 36 combinations, not just this row.")
    print()
    print("  The ratios T/D and S/D are convention-free (a uniform rescale of")
    print("  all three cancels) and are therefore measurements; the columns")
    print("  they are formed from are not, and are printed for shape only.")
    return head


def v3b_sign_map(sites):
    print()
    print("=" * 78)
    print("V3b  THE SIGN MAP: all 36 kernel x primitive x mass-model")
    print("     combinations, and it is NOT invariant")
    print("=" * 78)
    print("  The seed asks whether the VERDICT -- the sign structure -- is")
    print("  invariant across the four declarations even though the ratios")
    print("  provably are not. It is not. The answer splits cleanly by KERNEL")
    print("  FAMILY, and the split is the main result of this file.")
    print()
    print("  FOUR DECLARATIONS: KERNEL, PRIMITIVE and MASS MODEL are the three")
    print("  swept and are printed on every row; METRIC FORM is MOMENTUM-FREE")
    print("  on every row of this table and is stated here because a table that")
    print("  varies three of the four declarations is exactly where the fourth")
    print("  goes unsaid.")
    print()
    data = {}
    worst_dev = 0.0
    worst_dev_where = None
    for s in sites:
        for (kn, pn, m) in COMBOS:
            b, v0, gd = s.blocks(kn, pn, m)
            d = max(b[k][1] for k in BLOCKS)
            # `worst_dev_where is None or` -- same traceback class, third site,
            # and the most consequential: this sentinel is ALSO what V9's probe
            # is built from, so silencing the scalarity falsifier used to kill
            # the run here rather than reddening the row it was silencing.
            if worst_dev_where is None or d > worst_dev:
                worst_dev, worst_dev_where = d, (s.a, kn, pn, m)
            data.setdefault((kn, pn, m), []).append(
                (s.a, b["A"][0], b["E"][0], b["F"][0]))
    print(f"  worst per-block scalarity over the WHOLE sweep x 36 = "
          f"{worst_dev:.3e}")
    print(f"    at a = {worst_dev_where[0]:.4f}, {worst_dev_where[1]}, "
          f"{worst_dev_where[2]}, {worst_dev_where[3]}"
          f"   (criterion <= {TOL['scalarity']:.0e})")
    print()
    print("  THAT ROW ALSO CARRIES THE PER-ANGLE RANK CHECK, which is worth")
    print("  saying out loud because it is the only thing making 'the ranks are")
    print("  (1,2,3) at every angle' defensible from a gate that checks the rank")
    print("  DIRECTLY at only a handful of angles. `block_values` returns")
    print("  (nan, inf, inf) whenever a projector comes back with a rank the")
    print("  character table does not require, and an `inf` deviation propagates")
    print("  into `worst_dev` above and reddens THIS row. So a rank failure at")
    print("  any one of the swept angles, in either mass model, fails the gate")
    print("  through the scalarity criterion even though no rank is printed for")
    print("  that angle. Measured: under the E/F character-swap mutation and")
    print("  under the 8-of-12-elements mutation, that is exactly the path by")
    print("  which the sweep rows went red.")
    print()
    print(f"  {'kernel':20s} {'prim':7s} {'model':7s} {'D sign':>28s} "
          f"{'T sign':>28s}")
    flips = {}
    for key in sorted(data, key=lambda k: (COMBOS.index(k),)):
        rows = data[key]
        desc = {}
        for idx, bn in ((2, "E"), (3, "F")):
            ch = []
            for i in range(1, len(rows)):
                lo_v, hi_v = rows[i - 1][idx], rows[i][idx]
                if not (np.isfinite(lo_v) and np.isfinite(hi_v)):
                    # A non-finite block value means the labelling did not
                    # produce a block here at all (rank mismatch -> nan). It is
                    # NOT a sign change: `np.sign(nan) != np.sign(nan)` is True,
                    # so an unguarded comparison reports EVERY adjacent pair as
                    # a crossing and sends V3c bisecting hundreds of phantom
                    # brackets. Measured under the E/F character-swap mutation,
                    # which produced a 1.4 MB report of nothing.
                    continue
                if np.sign(hi_v) != np.sign(lo_v):
                    ch.append((rows[i - 1][0], rows[i][0]))
            flips[(key, bn)] = ch
            if not np.isfinite(rows[0][idx]):
                desc[bn] = "NOT LABELLED (non-finite)"
            elif not ch:
                desc[bn] = ("POSITIVE throughout" if rows[0][idx] > 0
                            else "NEGATIVE throughout")
            else:
                desc[bn] = "flips at " + ", ".join(
                    f"{lo:.1f}-{hi:.1f}" for lo, hi in ch)
        print(f"  {key[0]:20s} {key[1]:7s} {key[2]:7s} {desc['E']:>28s} "
              f"{desc['F']:>28s}")
    inv_pos = all(
        all(r[2] > 0 and r[3] > 0 for r in data[k])
        for k in data if k[0] in INVERSE_POWER)
    gauss_flip = all(
        bool(flips[(k, "E")]) or bool(flips[(k, "F")])
        for k in data if k[0] in GAUSSIAN)
    n_inv = sum(1 for k in data if k[0] in INVERSE_POWER)
    n_gau = sum(1 for k in data if k[0] in GAUSSIAN)
    ve_pos = all(
        data[k][0][2] > 0 and data[k][0][3] > 0 for k in data)
    print()
    print("  THE TWO CLAIMS, each gated, and each the other's non-vacuity:")
    print(f"  * all {n_inv} INVERSE-POWER combinations keep BOTH transverse")
    print(f"    blocks POSITIVE at every swept angle: {inv_pos}")
    print(f"  * all {n_gau} GAUSSIAN combinations have at least one transverse")
    print(f"    SIGN CHANGE inside the fundamental domain: {gauss_flip}")
    print("  If the first were reported alone it would be indistinguishable")
    print("  from a build that cannot produce a negative number at all. The")
    print("  second is the arm that shows the sweep can see one.")
    print(f"  * and at the GROUND STATE (a = 0) every one of the 36 has both")
    print(f"    transverse blocks positive: {ve_pos}   -- which is the")
    print("    inertia (6,0,0) already in the record, recovered here through a")
    print("    completely different labelling.")
    return data, flips, worst_dev, worst_dev_where, inv_pos, gauss_flip, ve_pos


def v3c_refine_flips(flips):
    print()
    print("=" * 78)
    print("V3c  WHERE THE TRANSVERSE SIGN CHANGES, refined by bisection")
    print("=" * 78)
    print(f"  Each bracket from V3b is refined by bisection to "
          f"{BISECT_TOL:.0e} degrees on")
    print("  the block value itself. THE BISECTION TOLERANCE IS NOT THE")
    print("  PRECISION OF THE ANSWER and must not be read as one: V3d measures")
    print("  the turnover MOVING with the step size h, and that movement is the")
    print("  real error bar. The tolerance sits two decades below it so that the")
    print("  bisection is not what limits the number. An earlier version stopped")
    print("  at 1e-04 -- the same size as the physical error bar -- while")
    print("  printing six decimals, and the two headline angles that went into")
    print("  the record were consequently wrong in their fifth and sixth")
    print("  decimals. The '+- h' column below is the number to quote from.")
    print()
    print("  A BRACKET INSIDE THE a = 60 GUARD BAND IS REFINED OR NOT ACCORDING")
    print("  TO THE KERNEL FAMILY, not according to the angle alone. The guard's")
    print("  justification is an INVERSE-POWER pole -- the twelve shared vertices")
    print("  merge in pairs and 1/r^p diverges -- and bisecting across a")
    print("  divergence returns the pole, which is the trap jb_t.path_critical")
    print("  records for dV/da sign scans. A GAUSSIAN DOES NOT DIVERGE THERE:")
    print("  V6's own table has the gauss rows at dev ~2e-08 straight through")
    print("  the band at every step size. Refusing to refine them was suppressing")
    print("  measurements this file's own machinery resolves, and one of the two")
    print("  suppressed turnovers is not even near the pole. The `dev` column")
    print("  below prices every in-band root against the file's own criterion.")
    print()
    print(f"  {'kernel':20s} {'prim':7s} {'model':7s} {'blk':4s} "
          f"{'bracket':>18s} {'root a':>14s} {'|lambda| there':>15s} "
          f"{'dev':>10s}  band")
    out = []
    n_brackets = 0
    n_refused = 0
    for (key, bn), brackets in sorted(
            flips.items(), key=lambda kv: (COMBOS.index(kv[0][0]), kv[0][1])):
        kn, pn, m = key
        for lo, hi in brackets:
            n_brackets += 1
            in_band = (lo < A_POLE < hi or abs(lo - A_POLE) < POLE_GUARD
                       or abs(hi - A_POLE) < POLE_GUARD)
            if in_band and _pole_guarded(kn):
                n_refused += 1
                print(f"  {kn:20s} {pn:7s} {m:7s} {bn:4s} "
                      f"{f'{lo:.1f}-{hi:.1f}':>18s} "
                      f"{'POLE (1/r^p)':>14s} {'not refined':>15s} "
                      f"{'--':>10s}  a=60")
                continue

            def f(a, kn=kn, pn=pn, m=m, bn=bn):
                b, _, _ = site(a).blocks(kn, pn, m)
                return b[bn][0]
            a_lo, a_hi = lo, hi
            f_lo = f(a_lo)
            f_hi = f(a_hi)
            for _ in range(60):
                if a_hi - a_lo < BISECT_TOL:
                    break
                mid = 0.5 * (a_lo + a_hi)
                f_mid = f(mid)
                if np.sign(f_mid) == np.sign(f_lo):
                    a_lo, f_lo = mid, f_mid
                else:
                    a_hi, f_hi = mid, f_mid
            root = 0.5 * (a_lo + a_hi)
            b_root, _, _ = site(root).blocks(kn, pn, m)
            dev_abs = b_root[bn][2]
            out.append((kn, pn, m, bn, root, dev_abs, in_band))
            print(f"  {kn:20s} {pn:7s} {m:7s} {bn:4s} "
                  f"{f'{lo:.1f}-{hi:.1f}':>18s} {root:14.6f} "
                  f"{min(abs(f_lo), abs(f_hi)):15.3e} {dev_abs:10.2e}"
                  f"  {'a=60' if in_band else ''}")
    n_band = sum(1 for r in out if r[6])
    band_dev = max((r[5] for r in out if r[6]), default=0.0)
    print()
    print(f"  {len(out)} turnovers located, of which {n_band} lie INSIDE the")
    print(f"  a = 60 guard band and are Gaussian. Their worst ABSOLUTE per-block")
    print(f"  deviation is {band_dev:.2e}; V3d converts every deviation in the")
    print("  `dev` column into the angular uncertainty it implies, which is the")
    print("  form in which it can be compared with the root.")
    print("  NOTE WHY THE RELATIVE `dev` OF V3 AND V3b IS NOT USED HERE. It is")
    print("  |Hb - lam I| DIVIDED BY |lam|, and a turnover is exactly where")
    print("  lam -> 0, so the relative statistic is structurally blind at the")
    print("  one place these measurements are taken. The column above is the")
    print("  ABSOLUTE deviation and V3d gates it.")
    print()
    print("  READ THIS TABLE AS THE ANSWER TO 'AT WHICH a DOES THE CHARACTER")
    print("  CHANGE'. It is not one angle: it depends on the kernel width, on")
    print("  the primitive AND on the mass model, and the doublet and the")
    print("  triplet do not turn over together. There is no kernel-independent")
    print("  'the sector becomes a ridge at a = X' to record.")
    print()
    print(f"  BRACKET ACCOUNTING, which is the teeth on the kernel-conditional")
    print(f"  guard and did not exist before it: V3b found {n_brackets} sign")
    print(f"  brackets and this section refined {len(out)} of them into roots,")
    print(f"  DISCARDING {n_refused}.")
    print(f"      (criterion: discarded == 0)")
    print("  EVERY BRACKET IS A MEASURED SIGN CHANGE. Refusing to refine one")
    print("  throws away a located feature of the landscape, so the honest")
    print("  criterion is not 'the guard is applied correctly' -- which the")
    print("  gate cannot see -- but 'nothing V3b found is dropped'. Measured:")
    print("  reverting the guard to angle-only makes this row read 34 of 36")
    print("  with 2 discarded and turns it RED, which is exactly the state the")
    print("  file was in before the guard was made kernel-conditional. Without")
    print("  this row that reversion changed the gate not at all.")
    print("  It is not a hardcoded count: both numbers come from V3b's own")
    print("  output, so adding kernels or moving the grid moves them together.")
    return out, n_brackets, n_refused


# The angular half-width used to measure d lambda / d a at a turnover. 0.01 deg
# moves a Gaussian block value by ~3e-03 against an absolute scalar deviation of
# ~1e-05, so the slope is measured three decades above the finite-difference
# floor that qualifies it.
ROOT_DA = 0.01


def v3d_root_reproducibility(roots):
    """Do the turnover angles REPRODUCE? -- and what is their real error bar.

    THE MOST QUOTABLE OUTPUT OF THIS FILE WAS, UNTIL THIS SECTION EXISTED, ITS
    ONLY ENTIRELY UNGATED ONE. `roots` flowed from V3c into V10's span line and
    into nothing else, so a build that located every turnover in the wrong place
    exited 0. Three things are measured here and all three are gated:

      1. HOW FAR THE ROOT MOVES WITH THE STEP SIZE h. Estimated by ONE Newton
         step from the H_MAIN root using the slope measured at the same h --
         3 evaluations per (root, h) instead of a full re-bisection's ~25, and
         exact to first order in a shift measured at ~1e-04 degrees. This is
         the number that says how many decimals of a turnover angle are real.
      2. THE SLOPE ITSELF, which is the TEETH. If lambda were flat through the
         crossing, every angle in the bracket would be a root, item 1 would be
         small for a reason having nothing to do with the root being determined,
         and the reproducibility row would be vacuous.
      3. THE ABSOLUTE per-block scalar deviation at the root, DIVIDED BY the
         slope, i.e. converted into the angular uncertainty it implies. This is
         the error control the file's relative `dev` cannot supply here: that
         statistic divides by |lambda| and the root is where lambda -> 0.
    """
    print()
    print("=" * 78)
    print("V3d  DO THE TURNOVER ANGLES REPRODUCE? -- the error bar on V3c")
    print("=" * 78)
    print("  FOUR DECLARATIONS: kernel, primitive and mass model as printed per")
    print("  row; METRIC FORM momentum-free throughout. Every root in V3c is")
    print("  carried through, not a chosen subset.")
    print()
    print(f"  {'kernel':20s} {'prim':7s} {'model':7s} {'blk':4s} "
          f"{'root (H_MAIN)':>14s} {'dlam/da':>11s} {'dev_abs':>10s} "
          f"{'dev -> deg':>11s} {'d(H_MAIN)':>10s} {'d(1e-3)':>10s} "
          f"{'d(1e-4)':>10s}")
    worst_spread = 0.0
    worst_spread_where = None
    min_slope = float("inf")
    worst_dev_deg = 0.0
    for (kn, pn, m, bn, root, dev_abs, _in_band) in roots:
        def lam(a, h, kn=kn, pn=pn, m=m, bn=bn):
            b, _, _ = site(a, h=h).blocks(kn, pn, m)
            return b[bn][0]
        shifts = {}
        slope0 = None
        for hh in (H_MAIN, 1e-3, 1e-4):
            sl = ((lam(root + ROOT_DA, hh) - lam(root - ROOT_DA, hh))
                  / (2.0 * ROOT_DA))
            if hh == H_MAIN:
                slope0 = sl
            # THE H_MAIN COLUMN IS NOT REDUNDANT AND IS THE ONE THAT CATCHES A
            # MIS-LOCATED ROOT. It is the Newton correction to the reported root
            # AT THE STEP SIZE THE ROOT WAS BISECTED AT, so it answers "is this
            # actually a root of the function it was bisected on" independently
            # of any h-dependence. A mutation widening the bisection tolerance
            # to 1e-1 leaves every reported root wrong by up to 0.05 degrees; the
            # h-flanking columns would catch that too, but only because they
            # re-evaluate -- this column catches it directly.
            shifts[hh] = -lam(root, hh) / sl if sl != 0.0 else float("inf")
        min_slope = min(min_slope, abs(slope0))
        dev_deg = dev_abs / abs(slope0) if slope0 != 0.0 else float("inf")
        worst_dev_deg = max(worst_dev_deg, dev_deg)
        sp = max(abs(v) for v in shifts.values())
        # `worst_spread_where is None or` -- the SAME sentinel class this fix
        # pass is removing from V1/V1b/V3b, caught in a section this fix pass
        # itself added. Writing it down because that is the failure mode the
        # project has now measured twice: a review-driven fix pass reintroduces
        # the defect class it is removing, inside the code written to remove it.
        if worst_spread_where is None or sp > worst_spread:
            worst_spread, worst_spread_where = sp, (kn, pn, m, bn, root)
        print(f"  {kn:20s} {pn:7s} {m:7s} {bn:4s} {root:14.6f} "
              f"{slope0:11.3e} {dev_abs:10.2e} {dev_deg:11.2e} "
              f"{shifts[H_MAIN]:10.2e} {shifts[1e-3]:10.2e} "
              f"{shifts[1e-4]:10.2e}")
    print()
    print(f"  worst Newton correction over {len(roots)} turnovers x 3")
    print(f"  step sizes (H_MAIN itself plus the two flanking it) = "
          f"{worst_spread:.3e} deg")
    print(f"      (criterion <= {TOL['root_h_spread']:.0e})")
    if worst_spread_where is None:
        # No turnovers located at all. Reported as a FAIL row rather than as an
        # IndexError: a build in which the sign map found nothing is exactly the
        # build whose verdict table a reader most needs to see.
        print("      NO TURNOVERS WERE LOCATED. The Gaussian arm of the verdict")
        print("      is the sweep's non-vacuity proof, so this is a failure and")
        print("      not an empty success.")
        worst_spread = float("inf")
    else:
        print(f"      at {worst_spread_where[0]}, {worst_spread_where[1]}, "
              f"{worst_spread_where[2]}, {worst_spread_where[3]}, "
              f"root {worst_spread_where[4]:.6f}")
    print(f"  smallest |d lambda / d a| at a root = {min_slope:.3e} per degree"
          f"   (criterion >= {TOL['root_slope_teeth']:.0e})")
    print("      THE TEETH. A flat lambda through the crossing would make the")
    print("      row above small for a reason that is not reproducibility.")
    print(f"  worst ABSOLUTE deviation expressed as an angle = "
          f"{worst_dev_deg:.3e} deg")
    print(f"      (criterion <= {TOL['root_dev_angle']:.0e})")
    print("      The relative `dev` of V3/V3b divides by |lambda| and the root")
    print("      is where lambda -> 0, so it is structurally blind here. This")
    print("      row is the absolute error control that replaces it.")
    print()
    print("  CONSEQUENCE FOR THE RECORD, and it is the whole reason this")
    print("  section exists: the turnover angles are reproducible to about")
    print(f"  {max(worst_spread, worst_dev_deg):.1e} degrees on the worst row,")
    print("  NOT to the six decimals the bisection can print. FOUR DECIMALS IS")
    print("  THE HONEST PRECISION -- and on the worst row the fourth is itself")
    print("  uncertain by about one unit, so a reader comparing two turnovers")
    print("  should treat a difference below 1e-03 degrees as no difference.")
    print("  That is what the README and T2 quote. The fifth and sixth decimals")
    print("  of a turnover angle in this file are not measurements.")
    ok = (worst_spread <= TOL["root_h_spread"]
          and min_slope >= TOL["root_slope_teeth"]
          and worst_dev_deg <= TOL["root_dev_angle"])
    print(f"\n  V3d PASSED: {ok}")
    return worst_spread, min_slope, worst_dev_deg, ok


# ==========================================================================
# V4 -- the ratios, which are the only convention-free magnitudes
# ==========================================================================

def v4_ratios(sites):
    print()
    print("=" * 78)
    print("V4  RATIOS -- the only form in which a magnitude here is a")
    print("    measurement")
    print("=" * 78)
    print("  Absolute lambda carries the arc's convention (coupling 1, total")
    print("  mass 1/2, R = 1); a uniform rescale of all three is free. T/D and")
    print("  S/D are not. Reported as SIGNED ratios of curvature -- not as")
    print("  ratios of omega, because omega does not exist where lambda < 0 and")
    print("  does not mean a frequency where dV != 0 even when lambda > 0.")
    print()
    landmarks = [A_VE, A_ICO, 45.0, A_2ND_THOMSON_VERTEX, 89.0]
    for kn in ("1/r^1  (Thomson)", "gauss s=0.5", "gauss s=1.0"):
        print(f"  --- KERNEL {kn} (raw), PRIMITIVE raw vertex, METRIC FORM "
              f"momentum-free ---")
        print(f"      {'a':>9s} {'model':7s} {'|dV|':>11s} {'T/D':>12s} "
              f"{'S/D':>12s} {'S/T':>12s}  signs (D,T)")
        for a in landmarks:
            s = site(a)
            for m in MODELS:
                b, _, gd = s.blocks(kn, "vertex", m)
                lam = {k: b[k][0] for k in BLOCKS}
                sg = _sgn(lam["E"]) + _sgn(lam["F"])
                print(f"      {a:9.4f} {m:7s} {gd:11.4e} "
                      f"{lam['F'] / lam['E']:12.6f} "
                      f"{lam['A'] / lam['E']:12.6f} "
                      f"{lam['A'] / lam['F']:12.6f}  {sg}")
        print()
    print("  A CONVENTION CROSS-CHECK, because the record's ratios and this")
    print("  file's are NOT the same quantity and a reader comparing them")
    print("  directly would find a discrepancy that is not there. jb_s S1 and")
    print("  jb_u quote ratios of OMEGA; the columns above are ratios of")
    print("  LAMBDA = omega^2. At the VE, where dV = 0 and omega exists, the")
    print("  square root of one must be the other exactly:")
    print(f"      {'model':7s} {'sqrt(T/D) here':>15s} {'record T/D':>12s} "
          f"{'sqrt(S/D) here':>15s} {'record S/D':>12s} {'worst rel':>11s}")
    worst_conv = 0.0
    ve = site(A_VE)
    for m in MODELS:
        b, _, _ = ve.blocks("1/r^1  (Thomson)", "vertex", m)
        lam = {k: b[k][0] for k in BLOCKS}
        got = (np.sqrt(lam["F"] / lam["E"]), np.sqrt(lam["A"] / lam["E"]))
        rec = PRIOR["omega_ratios_VE"][m]
        rel = max(abs(got[i] - rec[i]) / rec[i] for i in (0, 1))
        worst_conv = max(worst_conv, rel)
        print(f"      {m:7s} {got[0]:15.6f} {rec[0]:12.6f} {got[1]:15.6f} "
              f"{rec[1]:12.6f} {rel:11.3e}")
    print(f"      worst relative = {worst_conv:.3e}"
          f"   (criterion <= {TOL['ratio_convention']:.0e})")
    print("      Off a critical point NO square root is taken anywhere in this")
    print("      file, because lambda goes negative and omega would not exist.")
    print()
    print("  The two mass models do NOT agree on these ratios and are not")
    print("  expected to: off a critical point the mass model enters twice,")
    print("  through the block mass AND through Gamma, so jb_s S2c's per-block")
    print("  factors sqrt(8/5), sqrt(2), 2 have no validity here (jb_u,")
    print("  measured +13.7 / +16.4 / +48.7 percent at the icosahedron). They")
    print("  are not used anywhere in this file; each model is recomputed.")
    return worst_conv


# ==========================================================================
# V5 -- the critical points, where the reading IS a frequency
# ==========================================================================

def v5_critical_points(sites):
    print()
    print("=" * 78)
    print("V5  THE CRITICAL POINTS -- the only places on this path where the")
    print("    numbers are FREQUENCIES and the verdict is unhedged")
    print("=" * 78)
    print("  FOUR DECLARATIONS for the whole section: KERNEL 1/r^1 (Thomson)")
    print("  raw, PRIMITIVE raw VERTEX, MASS MODEL as printed per row in (a)")
    print("  and (c) and point in (b) and (d), METRIC FORM momentum-free")
    print("  throughout.")
    print()
    print("  Where dV = 0 the Riemannian Hessian IS the chart Hessian of V and")
    print("  the generalised eigenvalues ARE squared normal-mode frequencies.")
    print("  Everywhere else on this path the system is moving and they are")
    print("  curvature scales. This section is therefore the only one whose")
    print("  rows carry a stability reading in the ordinary sense.")
    print()
    print("  (a) THE VECTOR EQUILIBRIUM, against the record. jb_s S1 by way of")
    print("      jb_u U1, in (D, T, S) order rather than ascending -- the whole")
    print("      point of the character labelling is that ascending order is")
    print("      not a label.")
    ve = next(x for x in sites if abs(x.a - A_VE) < 1e-12)
    worst_ve = 0.0
    for m in MODELS:
        b, v0, gd = ve.blocks("1/r^1  (Thomson)", "vertex", m)
        got = np.array([b["E"][0], b["F"][0], b["A"][0]])
        rec = np.array(PRIOR[f"DTS_VE_{m}"])
        rel = float(np.abs(got - rec).max() / np.abs(rec).max())
        worst_ve = max(worst_ve, rel)
        print(f"      {m:7s} measured (D,T,S) = "
              f"[{got[0]:.6f}, {got[1]:.6f}, {got[2]:.6f}]")
        print(f"      {m:7s} record           = "
              f"[{rec[0]:.6f}, {rec[1]:.6f}, {rec[2]:.6f}]"
              f"   relative {rel:.3e}")
    print(f"      |dV| at the VE = {gd:.3e}")
    print(f"      worst relative = {worst_ve:.3e}"
          f"   (criterion <= {TOL['ve_record']:.0e})")
    print("      RESOLUTION CAP: the record constants carry six decimals, so")
    print("      nothing below ~8e-09 relative could mean anything here.")

    print()
    print("  (b) THE ICOSAHEDRON, the record's one non-critical anchor. |dV| is")
    print("      NOT zero there and this row is a CONFIGURATION check, not a")
    print("      spectrum check: it is what catches a sweep that is measuring")
    print("      the right quantity at the wrong place. |dV| here is the")
    print("      EUCLIDEAN norm of the covector in chart coordinates, which is")
    print("      the record's own convention -- NOT the g-norm of the gradient")
    print("      vector that V1b's path and transverse columns carry. The two")
    print("      differ by a factor of ~4.7 at this configuration.")
    ico = next(x for x in sites if abs(x.a - A_ICO) < 1e-12)
    _, _, gd_ico = ico.blocks("1/r^1  (Thomson)", "vertex", "point")
    d_ico = abs(gd_ico - PRIOR["grad_ico"])
    print(f"      |dV| = {gd_ico:.6f}   record {PRIOR['grad_ico']:.6f}"
          f"   |diff| {d_ico:.3e}   (criterion <= {TOL['ico_record']:.0e})")

    print()
    print("  (c) THE SECOND MINIMUM of raw Thomson on raw vertices, a =")
    print(f"      {A_2ND_THOMSON_VERTEX}. The record has inertia (6,0,0) there.")
    print("      Recovered here through the character labelling and the")
    print("      Riemannian Hessian, neither of which the record used.")
    s2 = next(x for x in sites
              if abs(x.a - A_2ND_THOMSON_VERTEX) < 1e-9)
    worst_2nd = 0.0
    for m in MODELS:
        b, v0, gd = s2.blocks("1/r^1  (Thomson)", "vertex", m)
        lam = {k: b[k][0] for k in BLOCKS}
        npos = sum(1 for k in BLOCKS for _ in range(MULT[k]) if lam[k] > 0)
        worst_2nd = max(worst_2nd, gd)
        print(f"      {m:7s} |dV| = {gd:.3e}   V = {v0:.6f}   "
              f"D {lam['E']:.6f}  T {lam['F']:.6f}  S {lam['A']:.6f}"
              f"   inertia ({npos},0,{6 - npos})")
    print(f"      worst |dV| = {worst_2nd:.3e}"
          f"   (criterion <= {TOL['second_min']:.0e})")
    print("      That criterion is what a wrong configuration breaks, and its")
    print("      DETECTION FLOOR is measured rather than argued: d|dV|/da here")
    print("      is 1.21 per degree, so the 1e-3 criterion resolves an angle")
    print("      error of ~8e-4 degrees. The 0.1% mutation of this constant")
    print("      lands the row at 9.0e-02, five decades above its baseline.")
    print("      The angle is a LOCAL constant in this file precisely so that a")
    print("      mutation probe can reach it -- jb_u's post-mortem records a")
    print("      probe that mutated an IMPORTED name, never applied, and read")
    print("      the resulting clean exit as confirmation.")

    print()
    print("  (d) A CRITICAL POINT OF V RESTRICTED TO THE PATH, which the record")
    print("      does not list and which is FORCED rather than found: a = 90 is")
    print("      the fixed point of the exact isometry a -> 180 - a, so V")
    print("      RESTRICTED TO THE PATH is even about it and dV/da vanishes")
    print("      there identically. The TRANSVERSE components are already zero")
    print("      at every angle (V1b), so it is the PATH component alone that")
    print("      has to collapse, and it does -- linearly, which is what an even")
    print("      function's derivative does at its symmetry point. Both columns")
    print("      are measured rather than narrated; an earlier version asserted")
    print("      the transverse components vanished 'for the same reason', which")
    print("      is a different argument and was not the one being made.")
    print("      THE HEADING SAYS 'RESTRICTED TO THE PATH' AND THAT QUALIFIER IS")
    print("      LOAD-BEARING, not caution. The isometry is a statement about V")
    print("      ALONG THE PATH and licenses only the path-tangential")
    print("      derivative. V1b covers the five CHART-transverse directions and")
    print("      is measured up to a = 89. But the local dimension at a = 90 is")
    print("      SEVEN (V7a), so there is a direction the 6-D chart never sees,")
    print("      and invariance of dV under the involution induced on a 7-D")
    print("      tangent cone forces only the (-1)-eigenspace to vanish, not the")
    print("      whole covector. 'a = 90 is a critical point of V' -- full stop,")
    print("      no restriction -- is therefore NOT what these rows show, and is")
    print("      not claimed here or in the record.")
    print(f"      {'a':>9s} {'path grad (g)':>13s} {'transverse (g)':>16s} "
          f"{'ratio to (90-a)':>16s}   Thomson / raw vertex / point")
    ratios = {}
    for a in (85.0, 88.0, 89.0, 89.5, 89.9, 89.99):
        s = site(a)
        _, ga, gt = s.gradient_split("1/r^1  (Thomson)", "vertex", "point")
        ratios[a] = ga / (90.0 - a)
        print(f"      {a:9.4f} {ga:13.6f} {gt:16.4e} "
              f"{ratios[a]:16.6f}")
    a90_rel = abs(ratios[89.99] / A90_RATE - 1.0)
    a90_lin = abs(ratios[89.99] - ratios[89.9])
    print("      THE COLLAPSE LAW IS NOW ASSERTED IN TWO INDEPENDENT WAYS, and")
    print("      it previously was not asserted at all -- this entire table fed")
    print("      nothing; `v5_critical_points` returned only the three anchor")
    print("      residuals. A mutation changing the ratio's denominator from")
    print("      (90-a) to 2(90-a) -- i.e. reporting a collapse law twice as")
    print("      slow as the measured one -- exited 0 with every gate row green.")
    print(f"      LINEARITY  |r(89.99) - r(89.9)| = {a90_lin:.3e}"
          f"   (criterion <= {TOL['a90_linear']:.0e})")
    print("          The ratio must be CONVERGING, i.e. the collapse is linear")
    print("          and what is left is the O((90-a)^2) correction.")
    print(f"      VALUE      |r(89.99)/A90_RATE - 1| = {a90_rel:.3e}"
          f"   (criterion <= {TOL['a90_rate']:.0e})")
    print(f"          against the local anchor A90_RATE = {A90_RATE}, raw")
    print("          1/r^1 (Thomson), raw VERTEX, point mass, momentum-free.")
    print("          The linearity row alone does NOT catch a wrong denominator")
    print("          -- halving it halves every ratio uniformly and leaves the")
    print("          convergence statistic untouched -- so the value row is the")
    print("          one that does. It is a REPRODUCIBILITY anchor written down")
    print("          from this file's own measurement, not an independent")
    print("          record, and it is labelled as such rather than dressed up.")
    print("      The VERTEX primitive is used here because the STRUT primitive")
    print("      has a pole of its own at a = 90 (V7c) and |dV| for an inverse")
    print("      power on struts diverges instead of collapsing. The vanishing")
    print("      is a property of V on the path, not of the primitive.")
    print("      WHAT IS AT a = 90 IS NOT CLAIMED. The chart's local dimension")
    print("      is 7 there, not 6 (memo R3, and V7 measures the failure), so")
    print("      a 6-D Hessian approaching it does not see the seventh")
    print("      direction and an inertia count from these rows would be")
    print("      counting in the wrong space. What the limit DOES show is that")
    print("      the two kernel families disagree about the character of that")
    print("      point -- Thomson approaches it with the path direction")
    print("      negative and both transverse blocks positive; gauss 1.0")
    print("      approaches it with all six negative. Settling it is")
    print("      inviscid-qvf.3 / inviscid-yli, not this bead.")
    return worst_ve, d_ico, worst_2nd, a90_rel, a90_lin


# ==========================================================================
# V6 -- the pole at a = 60, approached rather than omitted
# ==========================================================================

def v6_pole():
    print()
    print("=" * 78)
    print("V6  THE POLE AT a = 60, MEASURED -- not excluded in silence")
    print("=" * 78)
    print("  At a = 60 the twelve shared vertices merge in pairs, so some")
    print("  r_ij is exactly zero and every inverse power diverges. The memo")
    print("  records this as an INFINITE barrier ON THE PATH, derived rather")
    print("  than scanned. The question here is different: what does the")
    print("  TRANSVERSE curvature do as the pole is approached, and where does")
    print("  the measurement stop being one?")
    print()
    print("  The main sweep excludes |a - 60| < "
          f"{POLE_GUARD} for a reason this table")
    print("  makes visible: the per-block scalarity, which is the falsifier for")
    print("  reading an average as an eigenvalue, is O(h^2) truncation and V")
    print("  varies by orders of magnitude across a stencil of any usable width")
    print("  when the pole is a tenth of a degree away.")
    print()
    print("  FOUR DECLARATIONS on every row of this table: KERNEL as printed,")
    print("  PRIMITIVE raw VERTEX, MASS MODEL point, METRIC FORM momentum-free.")
    print("  The two primitives and both mass models are covered across the band")
    print("  by V6b below; this table varies the STEP SIZE, which is what it is")
    print("  for. 'dev' is the RELATIVE per-block scalar deviation.")
    print(f"  {'a':>9s} {'h':>8s} {'kernel':20s} {'prim':7s} {'model':7s} "
          f"{'lambda D':>14s} {'lambda T':>14s} {'sgn':>5s} {'dev':>10s}")
    rows = []
    for a in (59.5, 59.9, 59.99, 60.01, 60.1, 60.5):
        for h in (1e-3, 3e-4, 1e-4):
            for kn in ("1/r^1  (Thomson)", "gauss s=0.5"):
                try:
                    s = site(a, h=h)
                    b, v0, gd = s.blocks(kn, "vertex", "point")
                except (ChartUnmeasurable, IrrepLabelError) as exc:
                    print(f"  {a:9.4f} {h:8.0e} {kn:20s} {'vertex':7s} "
                          f"{'point':7s} FAIL -- {type(exc).__name__}: "
                          f"{str(exc).splitlines()[0][:40]}")
                    continue
                sg = _sgn(b["E"][0]) + _sgn(b["F"][0])
                dev = max(b[k][1] for k in BLOCKS)
                rows.append((a, h, kn, b["E"][0], b["F"][0], sg, dev))
                print(f"  {a:9.4f} {h:8.0e} {kn:20s} {'vertex':7s} "
                      f"{'point':7s} {b['E'][0]:14.4e} "
                      f"{b['F'][0]:14.4e} {sg:>5s} {dev:10.2e}")
    print()
    print("  READ THE dev COLUMN AS THE LIMIT OF THE METHOD, not as an error --")
    print("  AND THEN READ THE ROW WHERE IT EARNED ITS PLACE. An earlier")
    print("  version of this paragraph asserted that a relative deviation of")
    print("  1e-02 on a block value of 1e+06 could not move a SIGN. The table")
    print("  above refutes it: at a = 59.99 with h = 1e-03 the triplet comes")
    print("  back NEGATIVE with dev = 2.4 -- a deviation more than twice the")
    print("  quantity it qualifies -- and every smaller step size says")
    print("  positive. The falsifier caught a spurious sign. That is what it is")
    print("  for, and the claim it refutes is deleted rather than softened.")
    h_fine = min(r[1] for r in rows)
    inv = [r for r in rows if r[2].startswith("1/r")]
    gau = [r for r in rows if r[2].startswith("gauss")]
    worst_inv = max(r[6] for r in inv) if inv else float("nan")
    # ... and the SAME statistic restricted to H_MAIN, because the gate prices
    # the guard band by comparing inside against outside and the OUTSIDE number
    # is at H_MAIN. `worst_inv` above is a max over every step size in this
    # table and is attained at h = 1e-3, so quoting the two together compares
    # different discretisations. Both are returned; the gate uses the
    # like-for-like one and prints the other beside it.
    inv_hmain = [r for r in inv if r[1] == H_MAIN]
    worst_inv_h = max((r[6] for r in inv_hmain), default=float("nan"))
    sgn_inv = sorted({r[5] for r in inv if r[1] == h_fine})
    sgn_gau = sorted({r[5] for r in gau if r[1] == h_fine})
    dev_fine = max((r[6] for r in rows if r[1] == h_fine), default=float("nan"))
    bad = [r for r in rows
           if (r[5] not in (sgn_inv if r[2].startswith("1/r") else sgn_gau))]
    flagged = max((r[6] for r in bad), default=0.0)
    print()
    print(f"  THE READING IS TAKEN AT THE FINEST STEP IN THIS TABLE, "
          f"h = {h_fine:.0e},")
    print(f"  where the worst dev over every row is {dev_fine:.2e}:")
    print(f"    raw Thomson  / raw vertex / point mass: {', '.join(sgn_inv)}")
    print(f"    raw gauss0.5 / raw vertex / point mass: {', '.join(sgn_gau)}")
    print(f"  rows at a COARSER step that disagree with that reading: "
          f"{len(bad)}")
    print(f"  worst dev among them = {flagged:.2e}"
          f"   (criterion >= 1e+00: a sign error MUST come with a dev that")
    print("   announces it, or the dev column is not protecting anything)")
    print(f"  worst dev on the inverse power anywhere in the band = "
          f"{worst_inv:.2e}   (attained at h = 1e-03)")
    print(f"  ... and the same restricted to H_MAIN = {worst_inv_h:.2e}, which")
    print("  is the number the gate compares against the OUTSIDE worst, because")
    print("  that one is at H_MAIN too and a ratio between two different step")
    print("  sizes prices nothing.")
    print("  The pole-free case crosses a = 60 with BOTH transverse blocks")
    print("  NEGATIVE and the inverse power crosses it with both strongly")
    print("  POSITIVE. So the two kernel families do not merely differ in")
    print("  magnitude near the octahedron; they differ in the sign of the")
    print("  transverse curvature, on the same geometry, at the same angle.")
    print("  RATIOS INSIDE THE BAND ARE NOT QUOTED at any step size: the")
    print("  values move by 3% between h = 3e-04 and h = 1e-04 at a = 59.99,")
    print("  so they have not converged and a ratio built from them would be")
    print("  reporting the discretisation.")
    print()
    print("  AND THE IRONY THE MEMO ALREADY RECORDS, seen from the transverse")
    print("  side: an unregularised inverse power re-imposes as a wall exactly")
    print("  the region USER DECISION 16 declared legal, and it does so in ALL")
    print("  SIX directions at once, not only along the path. The wall is not")
    print("  a path artefact that a transverse detour could walk around.")
    print("  Whether a LOWER route exists is a minimum-energy-path question")
    print("  and is not answered here -- a positive transverse curvature says")
    print("  the valley walls rise, not how high the pass is.")
    return worst_inv, sgn_inv, sgn_gau, flagged, len(bad), worst_inv_h


def v6b_band_coverage():
    """ALL TWENTY inverse-power combinations, INSIDE the a = 60 guard band.

    THE CLAIM THIS SECTION EXISTS TO SUPPORT. V10 and the record say the
    inverse powers keep both transverse blocks positive throughout the
    fundamental domain. The main sweep excludes |a - 60| < 1 and a > 89, which
    is about 3.3% of the domain, and V6 above covers exactly ONE of the twenty
    inverse-power combinations inside the first band (1/r^1 on raw vertices with
    point masses). A coverage claim resting on a section that covers one row in
    twenty is a claim about a region the sweep never enters.

    So the region is entered. Every inverse power, both primitives, both mass
    models, at two step sizes, at SIX angles across the band -- and the
    per-block deviation is carried with the sign so that a sign taken from the
    band is qualified exactly as one taken outside it.

    WHICH SIX, AND WHY NOT a = 60 ITSELF. The half-degree grid would have
    visited 59.5, 60.0 and 60.5 inside this band. Two of those are here, along
    with four angles closer to the pole than the grid ever gets (59.9, 59.99,
    60.01, 60.1) precisely because that is where the method is expected to give
    out. a = 60 EXACTLY IS NOT AND CANNOT BE: the twelve shared vertices merge
    in pairs there, so some r_ij is exactly zero and every inverse power is
    infinite. That is the pole itself, not a gap in the coverage, and a table
    row asserting a sign for an infinite potential would be the over-claim this
    section exists to remove. The Gaussians ARE finite at 60 and V3c bisects
    them straight through it.

    AND THE FIRST RUN OF THIS SECTION IMMEDIATELY FOUND SOMETHING, which is
    why the reading is structured the way V6's is rather than as a flat "all
    positive". At h = 3e-4, within a hundredth of a degree of the pole, EIGHT
    of the twenty combinations come back with a NEGATIVE transverse block --
    and every one of them carries a per-block deviation of 1e+01 or more,
    while the same combination at h = 1e-4 says positive. That is precisely
    the artefact V6 records for 1/r^1 at 59.99, one step size finer and on
    the harder kernels. So the claim is stated as V6 states it: the reading is
    taken at the FINEST step measured, and a coarser row that disagrees must
    arrive with a deviation that announces it.
    """
    print()
    print("=" * 78)
    print("V6b  INSIDE THE a = 60 BAND: ALL 20 INVERSE-POWER COMBINATIONS")
    print("=" * 78)
    print("  FOUR DECLARATIONS: KERNEL, PRIMITIVE and MASS MODEL all swept and")
    print("  printed; METRIC FORM momentum-free throughout.")
    print()
    print(f"  {'a':>9s} {'h':>8s} {'combos':>7s} {'both blocks +':>14s} "
          f"{'worst dev':>11s} {'dev of neg':>11s} {'resolves':>9s}   "
          f"worst-dev combination")
    rows = []
    for a in (59.5, 59.9, 59.99, 60.01, 60.1, 60.5):
        for h in (H_MAIN, 1e-4):
            npos = 0
            n = 0
            wd = 0.0
            where = None
            dev_neg = float("inf")
            for (kn, pn, m) in COMBOS:
                if kn not in INVERSE_POWER:
                    continue
                try:
                    b, _, _ = site(a, h=h).blocks(kn, pn, m)
                except (ChartUnmeasurable, IrrepLabelError):
                    continue
                n += 1
                d = max(b[k][1] for k in BLOCKS)
                # A NON-FINITE block value is NOT a measured negative:
                # it is counted OUT of npos (so no row can claim
                # positivity on it) and its deviation is inf, which
                # drives the announcement statistic the right way.
                if b["E"][0] > 0 and b["F"][0] > 0:
                    npos += 1
                else:
                    dev_neg = min(dev_neg, d)
                if d > wd:
                    wd, where = d, f"{kn}/{pn}/{m}"
            res = wd <= TOL["band_dev"]
            rows.append((a, h, n, npos, wd, dev_neg, res, where))
            print(f"  {a:9.4f} {h:8.0e} {n:7d} {f'{npos}/{n}':>14s} "
                  f"{wd:11.2e} "
                  f"{'--' if npos == n else f'{dev_neg:.2e}':>11s} "
                  f"{'yes' if res else 'NO':>9s}   {where}")
    h_fine = min(r[1] for r in rows)
    fine = [r for r in rows if r[1] == h_fine]
    resolved = [r for r in rows if r[6]]
    unresolved = [r for r in rows if not r[6]]
    fine_pos = all(r[3] == r[2] and r[2] > 0 for r in fine)
    res_pos = all(r[3] == r[2] and r[2] > 0 for r in resolved)
    neg_rows = [r for r in rows if r[3] < r[2]]
    announce = min((r[5] for r in neg_rows), default=float("inf"))
    band_dev = max((r[4] for r in resolved), default=float("nan"))
    print()
    print(f"  THE READING IS TAKEN AT THE FINEST STEP IN THIS TABLE, "
          f"h = {h_fine:.0e},")
    print(f"  exactly as V6 takes its own:")
    print(f"    all 20 inverse-power combinations have BOTH transverse blocks")
    print(f"    POSITIVE at every angle in the band: {fine_pos}"
          f"   (criterion == True)")
    print(f"  rows at ANY step size reporting a negative: {len(neg_rows)}")
    print(f"    smallest per-block deviation among them = {announce:.2e}")
    print(f"      (criterion >= {TOL['band_sign_announce']:.0e} == 100x band_dev:")
    print("       a sign disagreement inside the band MUST come with a deviation")
    print("       that announces it, or the dev column is not protecting")
    print("       anything. NOT V6's own >= 1e+00 bar: the a > 89 measurement")
    print("       sits at 9.86e-01, just under it, and loosening a threshold to")
    print("       fit a measurement is what gate lesson 3 forbids -- so the")
    print("       criterion is stated as a multiple of what it protects)")
    print(f"  rows that RESOLVE at the file's own scalarity criterion "
          f"(dev <= {TOL['band_dev']:.0e}):")
    print(f"    {len(resolved)} of {len(rows)}; worst dev among them "
          f"{band_dev:.3e}   (criterion <= {TOL['band_dev']:.0e})")
    print(f"    and all of them positive: {res_pos}   (criterion == True)")
    print(f"  rows that do NOT resolve: {len(unresolved)}. THAT IS MOST OF THE")
    print("    TABLE and it is the honest state of the band: 1/r^12 on strut")
    print("    midpoints a tenth of a degree from a pole is not a converged")
    print("    measurement at any step size this chart supports. Those rows")
    print("    carry SIGNS ONLY, taken at the finest step, and no magnitude and")
    print("    no ratio from them enters any claim.")
    print()
    print("  WHAT THIS CHANGES IN THE VERDICT, and it is evidence rather than")
    print("  retraction: 'valley floor throughout the fundamental domain' was")
    print("  previously supported by a sweep that excluded this band and by ONE")
    print("  combination of twenty inside it. It is now supported by all twenty")
    print("  at the finest step, with the coarse-step disagreements and the")
    print("  unresolved rows both printed rather than absent.")
    print("  THE TEETH ARE ALREADY IN V6: the gauss 0.5 rows cross the same band")
    print("  with BOTH transverse blocks NEGATIVE at every step size, so 'both")
    print("  blocks positive here' is not something the band forces on any")
    print("  kernel that is evaluated in it.")
    print("  STILL NOT COVERED, and it is the reason V11 names the grid: this")
    print("  table visits the SAME half-degree lattice as the main sweep. A")
    print("  feature narrower than the grid, or an even number of crossings")
    print("  between two same-sign samples, is invisible to it exactly as it is")
    print("  to the main sweep.")
    ok = (fine_pos and res_pos and band_dev <= TOL["band_dev"]
          and announce >= TOL["band_sign_announce"])
    print(f"\n  V6b PASSED: {ok}")
    return (fine_pos, res_pos, band_dev, announce, len(resolved),
            len(unresolved), ok)


# ==========================================================================
# V7 -- a = 90: the chart failure, the isometry, and a pole nobody recorded
# ==========================================================================

def v7_branch_point_and_isometry(sites):
    print()
    print("=" * 78)
    print("V7  a = 90: THE CHART FAILURE, THE ISOMETRY, AND A SECOND POLE")
    print("=" * 78)
    print("  (a) THE CHART GIVES OUT, and the failure is reported rather than")
    print("      routed around. a = 90 is a genuine branch point of the")
    print("      variety: local dimension 7 against 6 (memo R3).")
    print(f"      {'a':>10s} {'Frame.dim':>10s} {'stencil':>28s}")
    dims = []
    for a in (89.0, 89.9, 89.99, 89.999, 90.0):
        d = Frame(a).dim
        try:
            site(a)
            st = "builds"
        except ChartUnmeasurable as exc:
            st = f"FAILS: {str(exc).splitlines()[0][:22]}"
        dims.append((a, d, st))
        print(f"      {a:10.4f} {d:10d} {st:>28s}")
    fails_at_90 = any(a == 90.0 and (d != 6 or not s.startswith("builds"))
                      for a, d, s in dims)
    ok_below = all(d == 6 and s.startswith("builds")
                   for a, d, s in dims if a < 90.0)
    print(f"      chart usable below 90: {ok_below};  fails AT 90: "
          f"{fails_at_90}")
    print("      NOTHING IS CLAIMED AT a = 90. The tangent cone there is bead")
    print("      inviscid-qvf.3 / inviscid-yli.")

    print()
    print("  (a-bis) THE APPROACH, PER PRIMITIVE AND PER STEP SIZE -- because")
    print("      the main sweep stops one degree short of 90 and a guard band")
    print("      applied without a measurement is an omission. What stops")
    print("      resolving is the STRUT primitive, and only it.")
    print("      FOUR DECLARATIONS on this table, none of which were previously")
    print("      stated anywhere in this section: KERNEL 1/r^3 raw, PRIMITIVE as")
    print("      printed per row, MASS MODEL lamina, METRIC FORM momentum-free.")
    print(f"      {'a':>9s} {'prim':7s} {'h':>8s} {'lambda D':>14s} "
          f"{'lambda T':>14s} {'sgn':>5s} {'dev':>10s}")
    approach = []
    for a in (89.0, 89.5, 89.9, 89.99):
        for pn in ("vertex", "strut"):
            for h in (1e-3, 3e-4, 1e-4):
                try:
                    b, _, _ = site(a, h=h).blocks("1/r^3", pn, "lamina")
                except (ChartUnmeasurable, IrrepLabelError) as exc:
                    print(f"      {a:9.4f} {pn:7s} {h:8.0e} "
                          f"FAIL -- {type(exc).__name__}")
                    continue
                dev = max(b[k][1] for k in BLOCKS)
                sg = _sgn(b["E"][0]) + _sgn(b["F"][0])
                approach.append((a, pn, h, dev))
                print(f"      {a:9.4f} {pn:7s} {h:8.0e} {b['E'][0]:14.4e} "
                      f"{b['F'][0]:14.4e} {sg:>5s} {dev:10.2e}")
    # `default=nan` on both: a mutation that removes one primitive's rows made
    # this a `ValueError: max() iterable argument is empty` raised BEFORE the
    # gate table printed, so a build that had deleted half the a = 90 evidence
    # crashed instead of reddening a row. nan fails both `<=` and `>=`, so the
    # missing column now arrives as two FAIL rows.
    dv = max((d for _, p, _, d in approach if p == "vertex"),
             default=float("nan"))
    ds = max((d for _, p, _, d in approach if p == "strut"),
             default=float("nan"))
    print(f"      worst dev, raw VERTEX = {dv:.2e}"
          f"   (criterion <= {TOL['branch_vertex']:.0e})")
    print(f"      worst dev, raw STRUT  = {ds:.2e}"
          f"   (criterion >= {TOL['branch_strut_teeth']:.0e}: it MUST fail the")
    print("                            vertex criterion, or 'the guard exists")
    print("                            because of the strut pole' is a claim")
    print("                            these rows do not support)")
    print("      BOTH NUMBERS ARE NOW RETURNED AND GATED. They were previously")
    print("      computed, printed and returned as neither, so the a = 90 guard")
    print("      was priced by narration while the a = 60 guard of V6 got two")
    print("      gate rows -- and the gate epilogue said the two were priced")
    print("      'the same way'. They are now.")
    print("      The vertex column stays at the level of the main sweep all the")
    print("      way to 89.99. The strut column does not, and (c) below is why.")
    print("      No verdict is taken from the strut rows inside this band.")
    print("      READ THE STRUT COLUMN AT a = 89.99 CAREFULLY: it is NOT")
    print("      monotone in h (1.25e+00 -> 2.26e+00 -> 3.57e-02 as h falls")
    print("      through 1e-3 / 3e-4 / 1e-4). A monotone O(h^2) column is what a")
    print("      truncation error looks like; this one is a stencil straddling a")
    print("      divergence, where the width of the stencil decides which side")
    print("      of the pole its outer points land on. It is a further reason no")
    print("      number is taken from these rows, and it is why the criterion")
    print("      above is stated on the MAX over the step sizes rather than on")
    print("      the finest.")

    print()
    print("  (a-ter) THE SECOND HALF OF THE COVERAGE CLAIM: all 20 inverse-power")
    print("      combinations INSIDE a > 89, which the main sweep also excludes.")
    print("      Same reason as V6b -- 'valley floor throughout the fundamental")
    print("      domain' cannot rest on a band the sweep never enters, and")
    print("      (a-bis) above covers one kernel of five.")
    print("      FOUR DECLARATIONS: KERNEL, PRIMITIVE and MASS MODEL all swept")
    print("      and printed; METRIC FORM momentum-free.")
    print(f"      {'a':>9s} {'h':>8s} {'prim':7s} {'combos':>7s} "
          f"{'both blocks +':>14s} {'worst dev':>11s} {'dev of neg':>11s} "
          f"{'resolves':>9s}")
    br_rows = []
    for a in (89.25, 89.5, 89.75, 89.9, 89.99):
        for h in (H_MAIN, 1e-4):
            for pn in ("vertex", "strut"):
                npos = 0
                n = 0
                wd = 0.0
                dev_neg = float("inf")
                for (kn, pn2, m) in COMBOS:
                    if kn not in INVERSE_POWER or pn2 != pn:
                        continue
                    try:
                        b, _, _ = site(a, h=h).blocks(kn, pn, m)
                    except (ChartUnmeasurable, IrrepLabelError):
                        continue
                    n += 1
                    d = max(b[k][1] for k in BLOCKS)
                    if b["E"][0] > 0 and b["F"][0] > 0:
                        npos += 1
                    else:
                        dev_neg = min(dev_neg, d)
                    wd = max(wd, d)
                if n == 0:
                    continue
                res = wd <= TOL["band_dev"]
                br_rows.append((a, h, pn, n, npos, wd, dev_neg, res))
                print(f"      {a:9.4f} {h:8.0e} {pn:7s} {n:7d} "
                      f"{f'{npos}/{n}':>14s} {wd:11.2e} "
                      f"{'--' if npos == n else f'{dev_neg:.2e}':>11s} "
                      f"{'yes' if res else 'NO':>9s}")
    br_hfine = min(r[1] for r in br_rows)
    br_fine_pos = all(r[4] == r[3] for r in br_rows if r[1] == br_hfine)
    br_resolved = [r for r in br_rows if r[7]]
    br_unresolved = [r for r in br_rows if not r[7]]
    br_res_pos = all(r[4] == r[3] for r in br_resolved)
    br_neg = [r for r in br_rows if r[4] < r[3]]
    br_announce = min((r[6] for r in br_neg), default=float("inf"))
    br_dev = max((r[5] for r in br_resolved), default=float("nan"))
    print(f"      THE READING IS TAKEN AT THE FINEST STEP, h = {br_hfine:.0e}:")
    print(f"        all 20 combinations, BOTH blocks positive at every angle: "
          f"{br_fine_pos}   (criterion == True)")
    print(f"      rows at any step size reporting a negative: {len(br_neg)}; "
          f"smallest dev among them {br_announce:.2e}")
    print(f"          (criterion >= {TOL['band_sign_announce']:.0e} == 100x "
          f"band_dev, the announcement rule V6b uses; V6's own bar is 1e+00)")
    print(f"      rows that RESOLVE: {len(br_resolved)} of {len(br_rows)}; "
          f"worst dev among them {br_dev:.2e}")
    print(f"          (criterion <= {TOL['band_dev']:.0e}), all positive: "
          f"{br_res_pos}")
    print("      EVERY UNRESOLVED ROW AND EVERY NEGATIVE IS A STRUT ROW, and")
    print("      they are exactly the rows (a-bis) prices: the strut midpoints")
    print("      are colliding, which is (c). The vertex primitive resolves")
    print("      across the whole band at both step sizes. So the coverage claim")
    print("      for a > 89 is carried by the vertex rows outright and by the")
    print("      strut rows at the finest step, with the coarse-step")
    print("      disagreements printed and priced rather than absent.")

    print()
    print("  (b) THE ISOMETRY a -> 180 - a, applied to the curvature blocks.")
    print("      FOUR DECLARATIONS: KERNEL as printed per row, PRIMITIVE raw")
    print("      vertex, MASS MODEL point, METRIC FORM momentum-free.")
    print("      The map is an exact isometry of the twelve shared vertices")
    print("      (three independent invariants, 1e-15; memo ESTABLISHED), so")
    print("      it must map the block values to themselves. That is a")
    print("      DERIVED requirement and a real check: a block mislabelling")
    print("      or a broken Hessian breaks it.")
    print(f"      {'a':>9s} {'kernel':20s} {'max |d lambda| rel':>20s}")
    worst_iso = 0.0
    for a in (13.0, 30.0, A_ICO, 45.0, 73.5, 89.0):
        for kn in ("1/r^1  (Thomson)", "gauss s=1.0"):
            r = {}
            for aa in (a, 180.0 - a):
                b, _, _ = site(aa).blocks(kn, "vertex", "point")
                r[aa] = {k: b[k][0] for k in BLOCKS}
            sc = max(abs(r[a][k]) for k in BLOCKS)
            d = max(abs(r[a][k] - r[180.0 - a][k]) for k in BLOCKS) / sc
            worst_iso = max(worst_iso, d)
            print(f"      {a:9.4f} {kn:20s} {d:20.3e}")
    print(f"      worst = {worst_iso:.3e}"
          f"   (criterion <= {TOL['isometry']:.0e})")
    ba, _, _ = site(30.0).blocks("1/r^1  (Thomson)", "vertex", "point")
    bb, _, _ = site(145.0).blocks("1/r^1  (Thomson)", "vertex", "point")
    sc = max(abs(ba[k][0]) for k in BLOCKS)
    teeth = max(abs(ba[k][0] - bb[k][0]) for k in BLOCKS) / sc
    print(f"      TEETH: the same statistic on the DELIBERATELY MISMATCHED")
    print(f"      pair (a = 30, 180 - 35 = 145) = {teeth:.3e}"
          f"   (criterion >= {TOL['isometry_teeth']:.0e})")
    print("      Without that row, a build in which every angle returned the")
    print("      same numbers would pass the isometry check outright.")

    print()
    print("  (c) A POLE THE RECORD DOES NOT HAVE, found while approaching 90.")
    print("      The memo records a = 90 as a total-collision state (all eight")
    print("      centroids at the origin). For the raw VERTEX primitive that is")
    print("      harmless -- the shared vertices stay well separated. For the")
    print("      raw STRUT-MIDPOINT primitive it is a SECOND POLE:")
    from jb_o_kernel_family import sq_dists
    from jb_q_strut_kernels import midpoints
    from jb_b_variety import PAIRS
    print(f"      {'a':>9s} {'min vertex sep':>16s} "
          f"{'min strut-midpoint sep':>24s} {'ratio to 2*(90-a) rad':>22s}")
    pole_lin = float("nan")
    vertex_floor = float("inf")
    for a in (75.0, 85.0, 89.0, 89.9, 89.99):
        X = corners(a)
        P = np.array([X[i, j] for (i, j), _ in PAIRS])
        dvx = float(np.sqrt(sq_dists(P).min()))
        dm = float(np.sqrt(sq_dists(midpoints(X)).min()))
        pred = 2.0 * np.radians(90.0 - a)
        # The floor is taken over the GUARD BAND ONLY (a >= 90 - BRANCH_GUARD),
        # because that is where the claim lives: the strut midpoints collapse
        # approaching 90 and the vertices do not. Over the wider set printed
        # here the minimum vertex separation is 0.5977 at a = 75, which is not
        # a counterexample to anything -- it is just a different angle, and a
        # teeth taken over it would be testing a claim nobody makes.
        if a >= A_BRANCH - BRANCH_GUARD:
            vertex_floor = min(vertex_floor, dvx)
        if abs(a - 89.99) < 1e-12:
            pole_lin = abs(dm / pred - 1.0)
        print(f"      {a:9.4f} {dvx:16.9f} {dm:24.9f} {dm / pred:22.9f}")
    print("      The separation falls LINEARLY to zero as 2 (90 - a) in")
    print("      radians. So every inverse power on the strut primitive")
    print("      diverges at a = 90 as well as at a = 60, and the vertex")
    print("      primitive does not.")
    print(f"      THE LAW ITSELF IS NOW ASSERTED, not merely printed:")
    print(f"      |ratio - 1| at a = 89.99 = {pole_lin:.3e}"
          f"   (criterion <= {TOL['pole_linearity']:.0e})")
    print(f"      min VERTEX separation over a >= "
          f"{A_BRANCH - BRANCH_GUARD:.0f} = {vertex_floor:.4f}"
          f"   (criterion >= {TOL['pole_vertex_floor']:.2f})")
    print("      The second row is the teeth of the first: the claim is that")
    print("      the pole belongs to ONE primitive, so the other one's minimum")
    print("      separation must NOT be collapsing over the same angles. It is")
    print("      taken over the GUARD BAND, which is where the claim lives; the")
    print("      wider table above dips to 0.5977 at a = 75, which is a")
    print("      different angle and not a counterexample to anything. Both")
    print("      rows exist because the ratio column was previously printed,")
    print("      quoted in the record as evidence, and asserted nowhere -- a")
    print("      mutation replacing the pole law 2(90-a) by 3(90-a) exited 0")
    print("      with every gate row green. It no longer does: the ratio")
    print("      becomes 2/3 and this row goes red.")
    print("      THIS IS A LANDSCAPE FACT, WHICH IS EXACTLY THE CLAUSE")
    print("      DECISION 17a LEFT LIVE: the two primitives share the ground")
    print("      state and do not share the barrier structure. The record has")
    print("      the a = 60 pole for both; the a = 90 pole belongs to one of")
    print("      them. It is measured here as a by-product and is NOT this")
    print("      bead's deliverable -- it deserves its own bead.")
    return (worst_iso, teeth, ok_below, fails_at_90, dv, ds,
            pole_lin, vertex_floor, br_fine_pos, br_res_pos, br_dev,
            br_announce)


# ==========================================================================
# V8 -- the metric form, confirmed for THIS sweep, and the second chart
# ==========================================================================

def v8_metric_form_and_second_chart(sites):
    print()
    print("=" * 78)
    print("V8  METRIC FORM AND CHART INVARIANCE, re-measured for this sweep")
    print("=" * 78)
    print("  qvf.9 corollary (iii) makes METRIC FORM a fourth declaration and")
    print("  requires the momentum-free form for anything transverse. That is")
    print("  what this file uses. The seed asks for one thing to be CONFIRMED")
    print("  rather than inherited, and it is the right thing to ask about:")
    print("  the base points here are ON the path, where jb_u U2c measures the")
    print("  two forms to coincide -- but Gamma needs dg in ALL SIX directions,")
    print("  and the stencil that produces dg steps TRANSVERSE to the path.")
    print("  Whether the coincidence survives that is a question about the")
    print("  stencil, not about the base point, and jb_u measured it at two")
    print("  angles only.")
    print()
    print("  (a) SECTION vs MOMENTUM-FREE, over the WHOLE sweep, in the")
    print("      quantity that actually matters -- the Christoffel symbols,")
    print("      which is where dg enters.")
    print(f"      {'a':>9s} {'model':7s} {'|g_s - g_h|':>13s} "
          f"{'|dg_s - dg_h|':>14s} {'|Gam_s - Gam_h|/|Gam|':>22s}")
    print("      BOTH ARMS ARE CONSTRUCTED EXPLICITLY here. An earlier version")
    print("      compared a freshly built SECTION geometry against the Site's")
    print("      own -- which is correct only while the Site's own is the")
    print("      momentum-free one, so a build that had quietly switched the")
    print("      default would compare a thing to itself and report 0.0 as")
    print("      confirmation. Measured: mutating the Site's default metric")
    print("      form to `section` made this whole section vacuous and the run")
    print("      exited 0 with every gate row green. The row (a-bis) below")
    print("      closes that, and this loop no longer reads the Site's form.")
    worst_mf = 0.0
    arms_differ = float("inf")
    shown = 0
    for s in sites:
        for m in MODELS:
            Gh = Geometry(s.st, forms(m, "horizontal"))
            Gs = Geometry(s.st, forms(m, "section"))
            dgam = float(np.abs(Gs.Gam - Gh.Gam).max()) \
                / float(np.abs(Gh.Gam).max())
            worst_mf = max(worst_mf, dgam)
            # THE TEETH THIS ROW WAS MISSING, and it is the M4 class one level
            # over: the row asserts that the two arms AGREE, so building both
            # arms from the same form makes it pass trivially. Measured: replacing
            # this loop's `forms(m, "section")` by `forms(m, "horizontal")` sent
            # `worst_mf` toward zero and the run exited 0 with every row green --
            # the same shape as the mutation that produced gate lesson 2, in the
            # very section written to fix it. So the AMBIENT forms the two arms
            # were built from are compared at the source, exactly as (a-bis)
            # compares the one the sweep used.
            arms_differ = min(
                arms_differ,
                float(np.abs(Gs.W - Gh.W).max())
                / float(np.abs(Gs.W).max()))
            if abs(s.a - round(s.a / 15.0) * 15.0) < 1e-9 and m == "point" \
                    and shown < 7:
                shown += 1
                print(f"      {s.a:9.4f} {m:7s} "
                      f"{float(np.abs(Gs.g - Gh.g).max()):13.3e} "
                      f"{float(np.abs(Gs.dg - Gh.dg).max()):14.3e} "
                      f"{dgam:22.3e}")
    print(f"      worst over ALL {len(sites)} angles x {len(MODELS)} models = "
          f"{worst_mf:.3e}   (criterion <= {TOL['metric_form']:.0e})")
    print(f"      ... and the SMALLEST difference between the two ARMS'")
    print(f"      ambient forms, over the same set = {arms_differ:.3e}")
    print(f"          (criterion >= {TOL['form_in_use']:.0e}: the two arms must")
    print("           be genuinely different forms, or the row above compares a")
    print("           thing to itself and cannot fail)")
    print("      CONFIRMED for this sweep: in the centroid chart the two forms")
    print("      agree to the finite-difference floor at every swept angle,")
    print("      INCLUDING in dg. jb_u's derivation says why -- every term of")
    print("      d(W Z (Z^T W Z)^-1 Z^T W) contracted as D^T (.) D carries a")
    print("      factor Z^T W D -- and the closed-form value of that factor on")
    print("      the path is measured next.")

    print()
    print("  (a-bis) WHICH FORM THE SWEEP ACTUALLY USED, asserted rather than")
    print("      declared. The two forms agreeing (row above) is exactly what")
    print("      makes 'we used the momentum-free one' unfalsifiable from the")
    print("      RESULTS, so it is checked at the source: the ambient form the")
    print("      Site handed to every Geometry, against the constant section")
    print("      form it must NOT be.")
    print(f"      {'model':7s} {'|W_sec|':>12s} {'max |W_sec - W_used|':>21s} "
          f"{'relative':>12s}")
    form_used = np.inf
    s_ref = sites[0]
    for m in MODELS:
        Wsec = forms(m, "section").at(s_ref.st.pts[tuple([0] * 6)])
        Wused = s_ref.geo[m].W
        d = float(np.abs(Wsec - Wused).max()) / float(np.abs(Wsec).max())
        form_used = min(form_used, d)
        print(f"      {m:7s} {float(np.abs(Wsec).max()):12.6f} "
              f"{float(np.abs(Wsec - Wused).max()):21.6f} {d:12.3e}")
    print(f"      smallest over the mass models = {form_used:.3e}"
          f"   (criterion >= {TOL['form_in_use']:.0e}: the projector must have")
    print("       removed something, or the sweep is running the section form)")

    print()
    print("  (b) THE FACTOR ITSELF, in closed form, over the whole sweep. No")
    print("      finite differences: Z^T W D with D the analytic Jacobian.")
    worst_leak = 0.0
    for s in sites:
        F = s.F
        D = s.D
        Z = rigid_fields(F.X0.reshape(-1))
        for m in MODELS:
            worst_leak = max(worst_leak,
                             float(np.abs(Z.T @ WEIGHTS[m] @ D).max()))
    print(f"      worst |Z^T W D| over the sweep, both models = "
          f"{worst_leak:.3e}   (criterion <= {TOL['gauge_leak']:.0e})")
    print("      The centroid gauge is momentum-free at every point of the")
    print("      symmetric path. jb_u U2a already shows this is a SYMMETRY")
    print("      ACCIDENT of the path and not an identity -- off the path the")
    print("      same statistic is 1e-04 -- which is why this file uses the")
    print("      momentum-free form regardless and does not lean on the")
    print("      accident.")

    print()
    print("  (c) CHART INVARIANCE. jb_u's deliverable was that two genuinely")
    print("      different charts agree once the Hessian is Riemannian and the")
    print("      metric is momentum-free. Re-run here on THIS file's own")
    print("      statistic -- the character block values -- because a")
    print("      chart-invariant construction with a chart-DEPENDENT labelling")
    print("      would still produce chart-dependent numbers.")
    print("      FOUR DECLARATIONS: KERNEL 1/r^1 (Thomson) raw, PRIMITIVE raw")
    print("      vertex, MASS MODEL point, METRIC FORM as printed per column")
    print("      (the two columns ARE the metric-form comparison).")
    print(f"      {'a':>9s} {'MOMENTUM-FREE rel':>19s} {'SECTION rel':>13s} "
          f"{'rank(A,E,F) origin':>20s}")
    worst_ci = 0.0
    worst_sec = 0.0
    min_sec = float("inf")
    origin_ranks = []
    for a in (5.0, A_ICO, 35.0, 50.0, 65.0, A_2ND_THOMSON_VERTEX, 85.0):
        row = {}
        for kind in ("horizontal", "section"):
            vals = {}
            for op in (False, True):
                s = site(a, kind=kind, origin_pivot=op)
                b, _, _ = s.blocks("1/r^1  (Thomson)", "vertex", "point")
                vals[op] = {k: b[k][0] for k in BLOCKS}
                if kind == "horizontal" and op:
                    rk = tuple(s.chars["point"].rank[k] for k in BLOCKS)
                    origin_ranks.append((a, rk))
            sc = max(abs(vals[False][k]) for k in BLOCKS)
            row[kind] = max(abs(vals[False][k] - vals[True][k])
                            for k in BLOCKS) / sc
        worst_ci = max(worst_ci, row["horizontal"])
        worst_sec = max(worst_sec, row["section"])
        min_sec = min(min_sec, row["section"])
        print(f"      {a:9.4f} {row['horizontal']:19.3e} "
              f"{row['section']:13.3e} {str(rk):>20s}")
    ranks_ok = all(rk == (1, 2, 3) for _, rk in origin_ranks)
    print(f"      worst MOMENTUM-FREE = {worst_ci:.3e}"
          f"   (criterion <= {TOL['chart_agree']:.0e})   <- THE DELIVERABLE")
    print(f"      worst SECTION       = {worst_sec:.3e}"
          f"   (criterion >= {TOL['chart_teeth']:.0e}: it MUST fail, or the")
    print("                                            row above is vacuous)")
    print(f"      SMALLEST SECTION over the same angles = {min_sec:.3e}. That")
    print("      is printed rather than gated, and it is the honest reading of")
    print("      the teeth: the MAX is attained at one angle, and at a = 85 the")
    print("      section form agrees to 2.6e-04, so the momentum-free row is")
    print("      near-vacuous AT THAT ANGLE even though it is not over the set.")
    print("      The gate keeps the MAX because the claim being defended is")
    print("      'the two forms are distinguishable by this statistic', which is")
    print("      an existence claim -- but a reader should know the separation")
    print("      is not uniform in a.")
    print(f"      CHARACTER RANKS IN THE ORIGIN CHART, at all "
          f"{len(origin_ranks)} angles: {ranks_ok}"
          f"   (criterion == (1,2,3) at every one)")
    print("      Previously PRINTED and asserted nowhere, while the README said")
    print("      'ranks (1,2,3) hold at every angle and in both charts' as an")
    print("      established fact -- the exact class this file's own gate lesson")
    print("      1 forbids. It is what says the labelling travels with the")
    print("      construction rather than being a property of one basis.")
    print("      The origin-pivot chart's slice is tilted relative to the")
    print("      centroid one (jb_u U2a: up to 32 degrees at the icosahedron),")
    print("      so this is not a re-parameterisation of the same subspace --")
    print("      which is what makes it evidence.")
    return (worst_mf, worst_leak, worst_ci, worst_sec, form_used, ranks_ok,
            min_sec, arms_differ)


# ==========================================================================
# V9 -- the step size, swept, with the h range for each threshold
# ==========================================================================

def v9_step_sweep(worst_where):
    print()
    print("=" * 78)
    print("V9  STEP SIZE: which thresholds hold over which range of h")
    print("=" * 78)
    print("  A threshold expressed in a quantity that carries O(h^2)")
    print("  truncation is FITTED to its h until a sweep says otherwise. FOUR")
    print("  of this file's thresholds are swept here -- the block scalarity,")
    print("  the singlet-tangent leak, the VE reproduction and the a -> 180-a")
    print("  isometry -- rather than quoted with a margin at one step size. The")
    print("  h RANGE over which each holds is printed at the bottom; the margin")
    print("  at H_MAIN is not the statistic being reported.")
    print()
    print("  AND WHICH OF THE FOUR IS ACTUALLY h-SENSITIVE IS NOW MEASURED")
    print("  RATHER THAN ASSERTED, because the column below refutes the label on")
    print("  one of them. The dynamic range of each swept column is printed at")
    print("  the bottom. The tangent leak is FLAT to four significant figures")
    print("  across the whole sweep -- P_A is built from the closed-form D and")
    print("  gx, so only the |.|_g norm touches the stencil -- and the VE")
    print("  reproduction is non-monotonic and passes at every step size. So the")
    print("  'all four' window below is set by TWO criteria pulling from")
    print("  opposite ends: scalarity is truncation and excludes the coarse")
    print("  steps, isometry is roundoff and excludes the finest. Its count of 2")
    print("  sits at ZERO MARGIN against a criterion of 2. Both facts are stated")
    print("  here and in V11 rather than left for a reader to derive.")
    print()
    a_w, kn_w, pn_w, m_w = worst_where
    print("  THE PROBE IS THE WORST COMBINATION V3b ACTUALLY FOUND, passed in")
    print("  rather than hardcoded -- a probe fixed at a guess stops being the")
    print("  worst case the first time the grid or the guard band changes:")
    print(f"      a = {a_w:.4f}, {kn_w}, {pn_w}, {m_w}")
    probe = [(a_w, kn_w, pn_w, m_w)]
    print()
    print("  THE ISOMETRY IS SWEPT HERE TOO, and it was not before. An")
    print("  independent whole-gate h sweep found the gate RED at h = 1e-4 --")
    print("  FINER than the quoted step -- on V7's isometry row (1.224e-06")
    print("  against a 1e-06 criterion, roundoff-limited rather than")
    print("  truncation-limited), while this section's three thresholds all")
    print("  passed there and this section therefore reported h = 1e-4 as")
    print("  acceptable. A per-threshold window that does not include every")
    print("  h-sensitive threshold in the gate is not the gate's window. The")
    print("  step size h = 6e-4 is in the list for the same reason: it is one of")
    print("  the two the whole gate is actually green at.")
    print(f"  {'h':>8s} {'newton res':>11s} {'r/h^2':>10s} "
          f"{'worst scalarity':>16s} {'tangent leak':>14s} "
          f"{'VE reproduction':>16s} {'isometry':>12s} {'all four':>10s}")
    rows = []
    e0 = np.zeros(6)
    e0[0] = 1.0
    for h in (3e-3, 1e-3, 6e-4, 3e-4, 1e-4):
        try:
            dev = 0.0
            res = 0.0
            for (a, kn, pn, m) in probe:
                s = site(a, h=h)
                b, _, _ = s.blocks(kn, pn, m)
                dev = max(dev, max(b[k][1] for k in BLOCKS))
                res = max(res, s.st.residual)
            leak = 0.0
            for a in (A_VE, A_ICO, 45.0, 85.0):
                s = site(a, h=h)
                res = max(res, s.st.residual)
                for m in MODELS:
                    leak = max(leak, s.chars[m].project_out(e0))
            ve = site(A_VE, h=h)
            rel = 0.0
            for m in MODELS:
                b, _, _ = ve.blocks("1/r^1  (Thomson)", "vertex", m)
                got = np.array([b["E"][0], b["F"][0], b["A"][0]])
                rec = np.array(PRIOR[f"DTS_VE_{m}"])
                rel = max(rel, float(np.abs(got - rec).max()
                                     / np.abs(rec).max()))
            iso = 0.0
            for a in (13.0, 30.0, A_ICO, 45.0, 73.5, 89.0):
                for kn in ("1/r^1  (Thomson)", "gauss s=1.0"):
                    r = {}
                    for aa in (a, 180.0 - a):
                        b, _, _ = site(aa, h=h).blocks(kn, "vertex", "point")
                        r[aa] = {k: b[k][0] for k in BLOCKS}
                    sc = max(abs(r[a][k]) for k in BLOCKS)
                    iso = max(iso, max(abs(r[a][k] - r[180.0 - a][k])
                                       for k in BLOCKS) / sc)
        except (ChartUnmeasurable, IrrepLabelError) as exc:
            print(f"  {h:8.0e} FAIL -- not measurable at this step size: "
                  f"{type(exc).__name__}")
            print(f"           {str(exc).splitlines()[0][:120]}")
            continue
        good = (dev <= TOL["scalarity"] and leak <= TOL["tangent_leak"]
                and rel <= TOL["ve_record"] and iso <= TOL["isometry"])
        rows.append((h, dev, leak, rel, good, iso))
        print(f"  {h:8.0e} {res:11.2e} {res / h ** 2:10.2e} {dev:16.3e} "
              f"{leak:14.3e} {rel:16.3e} {iso:12.3e} "
              f"{'PASS' if good else 'FAIL':>10s}")
    if not rows:
        print("  NO step size in this sweep produced a measurement.")
        return rows, 0, 0.0, False
    good = [r for r in rows if r[4]]
    devs = [r[1] for r in rows]
    # A COLUMN OF EXACT ZEROS IS THE SILENCED FALSIFIER, NOT A DIVISION ERROR.
    # `max/min` raised ZeroDivisionError before the gate table existed when the
    # per-block deviation was mutated to a constant 0.0 -- i.e. the one probe
    # this shape row is FOR killed the run instead of reddening it. A flat-zero
    # column has no dynamic range, so the shape is reported as zero and the row
    # goes RED, which is what "the column must vary" means.
    shape = (max(devs) / min(devs)) if min(devs) > 0 else 0.0
    span = rows[0][0] / rows[-1][0]
    print()
    _fall = (devs[0] / devs[-1]) if devs[-1] > 0 else float("nan")
    print(f"  scalarity falls by {_fall:.2e} across a "
          f"{span:.0f}x range of h; pure O(h^2) predicts {span ** 2:.2e}.")
    print("  That is the shape the column is required to have. A FLAT column")
    print("  would mean the sweep is not measuring discretisation at all, and")
    print("  a flat-AND-LARGE column is exactly what a broken construction")
    print("  produces -- which is why the shape is gated as well as the value.")
    print()
    print("  THE h RANGE PER THRESHOLD, which is the thing the seed asks for,")
    print("  WITH THE DYNAMIC RANGE OF EACH COLUMN BESIDE IT -- the second")
    print("  number is what says whether calling the threshold 'h-sensitive' is")
    print("  a measurement or a habit:")
    ranges = {}
    for name, idx, key in (("block scalarity", 1, "scalarity"),
                           ("singlet-tangent leak", 2, "tangent_leak"),
                           ("VE reproduction", 3, "ve_record"),
                           ("isometry (V7b)", 5, "isometry")):
        holds = [r[0] for r in rows if r[idx] <= TOL[key]]
        col = [r[idx] for r in rows]
        rng = max(col) / min(col) if min(col) > 0 else float("inf")
        ranges[key] = rng
        print(f"    {name:22s} <= {TOL[key]:.0e}   range {rng:8.2e}   "
              f"holds at h = "
              + (", ".join(f"{x:.0e}" for x in holds) if holds else "NONE"))
    print("      The scalarity column has a range worth the name. The leak's")
    print(f"      range is {ranges['tangent_leak']:.2e} -- flat -- so the")
    print("      annotation that used to call it O(h^2) was refuted by the very")
    print("      sweep it pointed at, and has been corrected at the threshold.")
    print("      The isometry column is NOT truncation-limited either: it gets")
    print("      WORSE as h falls, which is roundoff, and it is what excludes")
    print("      the finest step.")
    print(f"    all four together               holds at h = "
          + (", ".join(f"{r[0]:.0e}" for r in good) if good else "NONE"))
    print(f"      count = {len(good)} of {len(rows)} measured"
          f"   (criterion >= {TOL['h_window']})")
    print("      ZERO MARGIN, stated rather than left to be noticed: the count")
    print("      equals the criterion. Two criteria do the excluding, and from")
    print("      OPPOSITE ends -- scalarity (truncation) rules out the coarse")
    print("      steps, isometry (roundoff) rules out the finest -- so a future")
    _dev_hmain = next((r[1] for r in rows if r[0] == H_MAIN), None)
    if _dev_hmain:
        print(f"      worst combination about {TOL['scalarity'] / _dev_hmain:.0f}x"
              f" worse than today's would take")
        print("      H_MAIN out of the passing set and redden this row.")
    print()
    print("  WHAT THIS SECTION STILL DOES NOT MEASURE, and an independent audit")
    print("  quantified it: the scalarity column above is ONE combination -- the")
    print("  worst one AT H_MAIN -- and it does not stay the worst as h moves.")
    print("  A whole-gate sweep (T2 inviscid/qvf.4-check-suite-validation.md)")
    print("  measured the full 36-combination worst at h = 1e-4 to be 6.654e-04")
    print("  against this probe's 1.259e-05, a factor of 53. So the ranges above")
    print("  are OPTIMISTIC away from H_MAIN, and the same audit found the whole")
    print("  gate green only at h = 6e-4 and 3e-4 -- which is the pair this")
    print("  table now reports, but by a cheaper route than re-running the")
    print("  sweep, and the agreement is not a proof that the routes coincide.")
    print(f"    dynamic range max/min of the scalarity column = {shape:.2e}"
          f"   (criterion >= {TOL['h_shape']:.0e})")
    print()
    print(f"  THE MAIN SWEEP IS QUOTED AT h = {H_MAIN:.0e}, which is jb_u U5's")
    print("  measured minimum for the same chart and the same stencil. That is")
    print("  a DERIVED choice, not one fitted here, and this table is where it")
    print("  is checked rather than assumed. Note what the table says about")
    print("  jb_u's own quoted step: at h = 1e-3 the scalarity criterion of")
    print("  THIS file does not hold on its worst combination. That is not a")
    print("  defect in jb_u, whose gated statistic is a different one -- it is")
    print("  why this file did not simply inherit the step size.")
    ok = len(good) >= TOL["h_window"] and shape >= TOL["h_shape"]
    print(f"\n  V9 PASSED: {ok}")
    return rows, len(good), shape, ok


# ==========================================================================
# V10 -- the verdict, and what it does and does not license
# ==========================================================================

def v10_verdict(inv_pos, gauss_flip, ve_pos, roots, band60_pos, band90_pos,
                root_err):
    print()
    print("=" * 78)
    print("V10  THE VERDICT")
    print("=" * 78)
    print("  THE QUESTION THE BEAD ASKS: is the symmetric 1-DOF sector a")
    print("  transverse VALLEY FLOOR, a RIDGE, or does it change character")
    print("  along the path -- and at which a.")
    print()
    print("  THE ANSWER: IT CHANGES CHARACTER, AND NEITHER WHETHER NOR WHERE")
    print("  IS A PROPERTY OF THE GEOMETRY ALONE. WHETHER it changes is a")
    print("  property of the KERNEL FAMILY -- no inverse power does, every")
    print("  Gaussian does. WHERE it changes depends on the kernel, the")
    print("  PRIMITIVE and the MASS MODEL, all three.")
    print()
    print("  * AT THE GROUND STATE (a = 0, the vector equilibrium) the sector")
    print("    is a transverse VALLEY FLOOR for every one of the 36 kernel x")
    print(f"    primitive x mass-model combinations: {ve_pos}. This is the")
    print("    recorded inertia (6,0,0), recovered through a labelling the")
    print("    record did not use. The saddle scenario the bead was opened to")
    print("    rule out does NOT hold at the ground state.")
    print()
    print("  * ALONG THE PATH the two kernel families part company.")
    print(f"    - every INVERSE POWER (1/r^p, p = 1,2,3,6,12), both primitives,")
    print(f"      both mass models: VALLEY FLOOR at every swept angle. "
          f"{inv_pos}")
    print(f"    - ... and inside the |a - 60| < {POLE_GUARD} band the main grid")
    print(f"      excludes, all twenty at the finest step: {band60_pos} (V6b)")
    print(f"    - ... and inside the a > {90 - BRANCH_GUARD:.0f} band it also")
    print(f"      excludes, all twenty at the finest step: {band90_pos} "
          f"(V7 a-ter)")
    print("      THOSE TWO ROWS ARE WHY THE WORD 'THROUGHOUT' IS NOW USED AT")
    print("      ALL. The main sweep never enters either band -- together about")
    print("      3.3% of the fundamental domain -- and until V6b and V7(a-ter)")
    print("      existed, the coverage inside them was one inverse-power")
    print("      combination of twenty in the first band and one kernel of five")
    print("      in the second. The claim was true; the evidence for it was not")
    print("      in the file. THE RESIDUE IS STATED EXACTLY AND IT IS LARGE:")
    print("      most rows in the a = 60 band do NOT meet this file's own")
    print("      scalarity criterion at any step size it supports -- 1/r^12 on")
    print("      strut midpoints a tenth of a degree from a pole is not a")
    print("      converged measurement -- and at h = 3e-4 within a hundredth of")
    print("      a degree of the pole EIGHT of the twenty come back negative,")
    print("      each with a deviation above 1e+01, while every one of them says")
    print("      positive at h = 1e-4. The band therefore carries SIGNS ONLY,")
    print("      taken at the finest step, exactly as V6 already said for the")
    print("      one combination it covered. The a > 89 residue is the STRUT")
    print("      rows only, and precisely: at h = 3e-4 they stop resolving from")
    print("      a = 89.75 upward, at h = 1e-4 from a = 89.9 upward, and the")
    print("      VERTEX rows resolve across the whole band at both step sizes.")
    print(f"    - every GAUSSIAN (s = 0.5, 1.0, 1.5, 2.5), both primitives,")
    print(f"      both mass models: the transverse curvature CHANGES SIGN")
    print(f"      inside the fundamental domain. {gauss_flip}")
    print("      The sector becomes a transverse RIDGE -- in the doublet, in")
    print("      the triplet, or in both, and the doublet and the triplet do")
    print("      not turn over at the same angle.")
    if roots:
        lo = min(r[4] for r in roots)
        hi = max(r[4] for r in roots)
        print(f"    - the located turnovers span a = {lo:.4f} to {hi:.4f}")
        print(f"      (+- {root_err:.0e} deg, V3d) across the Gaussian")
        print("      combinations (V3c lists them one by one). There is no")
        print("      single angle to record.")
        print("      FOUR DECIMALS, NOT SIX. The bisection runs to 1e-07 so that")
        print("      it is not the limiting error, but V3d measures the root")
        print("      MOVING with the step size at ~1e-04 deg, and that is the")
        print("      error bar. A previous version of this line printed six")
        print("      decimals off a bisection that stopped at 1e-04, and the")
        print("      fifth and sixth digits that reached the record were wrong.")
    print()
    print("  * SO THE SIGN STRUCTURE IS NOT INVARIANT ACROSS THE FOUR")
    print("    DECLARATIONS. The seed asked whether the VERDICT survives the")
    print("    choices even though the ratios provably do not. It does not.")
    print("    The kernel was already known to set the spectrum without moving")
    print("    the ground state (DECISION 17a); it is now measured to set the")
    print("    TOPOGRAPHY of the landscape around the path as well, including")
    print("    its sign. Choosing the kernel is choosing whether the symmetric")
    print("    sector is a valley or a ridge away from the ground state.")
    print()
    print("  WHAT THIS DOES NOT LICENSE, and the seed's scope correction 2 is")
    print("  the governing text here.")
    print("  * IT IS NOT A DYNAMICAL STABILITY VERDICT FOR A MOVING")
    print("    TRAJECTORY. Off a critical point nothing oscillates, and the")
    print("    normal variational equation along the motion -- the Floquet-type")
    print("    problem that would decide whether a trajectory started near the")
    print("    symmetric path stays near it -- requires an equation of motion")
    print("    integrated. This project has never integrated one. A negative")
    print("    transverse curvature at a point the system passes THROUGH at")
    print("    speed is not the same object as an unstable mode, and this file")
    print("    does not convert one into the other. THAT CALCULATION IS NEEDED")
    print("    FOR AN HONEST DYNAMICAL VERDICT AND IT GETS ITS OWN BEAD.")
    print("  * IT DOES NOT INVALIDATE THE EXISTING 1-DOF RESULTS, and it does")
    print("    not confirm them either. M_eff's 3:1 and 9:1 and the fact that")
    print("    adot peaks at the vector equilibrium are statements about the")
    print("    symmetric sector's OWN kinematics. That sector is the")
    print("    fixed-point set of the chiral tetrahedral symmetry, so a")
    print("    trajectory started exactly on it stays on it EXACTLY, whatever")
    print("    the transverse curvature does -- V1 measures the singlet to be")
    print("    the path tangent at every swept angle, which is the same")
    print("    statement. Those results are therefore correct AS statements")
    print("    about that invariant sub-model, unconditionally.")
    print("    What was at issue is whether the sub-model is REPRESENTATIVE of")
    print("    the medium. On that: at the ground state, yes -- the valley")
    print("    walls rise in all five transverse directions for every kernel.")
    print("    Away from it, it depends on the kernel, and for the Gaussians")
    print("    there are stretches of the path where the transverse curvature")
    print("    is negative. Turning that into 'the 1-DOF motion is")
    print("    unrepresentative' needs the dynamical calculation above, not")
    print("    this one.")
    print()
    print("  * AND THE ONE PIECE OF EMPIRICAL EVIDENCE ANYONE HAS OFFERED ON")
    print("    THAT QUESTION IS NOT A COMPUTATION: bead inviscid-qvf.11. The")
    print("    bead comment of 2026-08-16 directs that it be read before this")
    print("    paragraph is written, and it is discharged here rather than")
    print("    skipped. THE OWNER BUILT THIS BEAD'S 1-DOF SECTOR IN WOOD AND")
    print("    WIRE: jitterbug units wired at the triangle vertices, each")
    print("    triangle additionally constrained by a dowel through its centre,")
    print("    which is mechanically the removal of five of the six internal")
    print("    degrees of freedom -- the symmetric sector, realised. In phase,")
    print("    the array LOCKS at the icosahedral phase, ONE-SIDEDLY: it cannot")
    print("    expand, it can contract. qvf.11's fork is (A) a tension-only wire")
    print("    going taut, which would be an artefact of the rig and touch")
    print("    nothing here, against (B) a genuine boundary or cusp of the")
    print("    array's configuration variety, which would bind the model.")
    print()
    print("    WHAT THIS FILE CAN SAY, AND IT IS NARROW. It measures ONE unit.")
    print("    A one-sided obstruction is a CONSTRAINT phenomenon -- an")
    print("    inequality becoming active -- and a curvature sign is not that")
    print("    object, so nothing in the transverse Hessian can confirm or")
    print("    refute the lock. What the single unit CAN report is whether its")
    print("    own configuration space does anything at the icosahedral phase,")
    print("    since fork (B) would need a degeneracy there and this file")
    print("    already measures the two places where the variety genuinely does")
    print("    misbehave. MEASURED, all raw / momentum-free:")
    # THE DEGENERACY STATISTIC HAD TO BE CHOSEN WITH CARE and the obvious one is
    # wrong. `ev_max / ev_min` at a = 60 is 5/3, a perfectly ordinary number --
    # the mass metric there is NOT singular, it has a five-fold eigenvalue
    # (V2: the spectrum is [1/32 x5, 5/96]). A ratio of extremes cannot see
    # that. What sees it is the MULTIPLICITY of the smallest eigenvalue: 2 away
    # from the octahedron (the doublet), 5 at it (doublet and triplet merged).
    # Reporting the ratio and calling it degeneracy would have been prose
    # asserting more than the adjacent number supports.
    _ico = site(A_ICO)

    def _spec(a):
        ev = np.linalg.eigvalsh(_g_closed_form(a)[0])
        mult = int(sum(1 for x in ev if abs(x - ev[0]) <= 1e-12 * abs(ev[0])))
        return float(ev[0]), mult
    _lam_min = min(min(_ico.blocks(kn, pn, m)[0][b][0] for b in ("E", "F"))
                   for (kn, pn, m) in COMBOS)
    print(f"      {'a':>12s} {'dim':>5s} {'min ev(g)':>12s} "
          f"{'mult of min ev':>15s}")
    for _a, _note in ((A_ICO, "the icosahedral phase"),
                      (A_POLE, "a=60: mass metric degenerate (V2)"),
                      (A_BRANCH, "a=90: local dimension jumps (V7a)")):
        _d = Frame(_a).dim
        if _d == 6:
            _mn, _mu = _spec(_a)
            print(f"      {_a:12.7f} {_d:5d} {_mn:12.6f} {_mu:15d}   {_note}")
        else:
            print(f"      {_a:12.7f} {_d:5d} {'--':>12s} {'--':>15s}   "
                  f"{_note}")
    print(f"      smallest transverse lambda at the icosahedral phase over all")
    print(f"      {len(COMBOS)} combinations = {_lam_min:.4e}  -- positive, so")
    print("      no block is going soft there either.")
    print("    SO THE SINGLE UNIT'S CONFIGURATION SPACE IS AN ORDINARY INTERIOR")
    print("    POINT AT a = 22.2387561: chart dimension 6, mass metric")
    print("    non-degenerate, every transverse block stiff -- and this file")
    print("    detects both of the places where the variety is NOT ordinary, so")
    print("    the null result is not a blind instrument. THAT IS EVIDENCE FOR")
    print("    FORK (A) AND IT IS NOT PROOF OF IT: the array is a DIFFERENT")
    print("    variety from one unit, with inter-unit constraints this file")
    print("    never assembles, and a boundary can be created by the coupling")
    print("    without existing in any single unit. qvf.11's own P2 (assemble")
    print("    the array's constraint Jacobian with BILATERAL constraints) and")
    print("    P3 (DOF count, doweled against free) are the tests that settle")
    print("    it, and P0 -- change the wire lengths and see whether the lock")
    print("    angle moves -- settles it more cheaply than any of them. This")
    print("    file's contribution to qvf.11 is one negative datum, and the")
    print("    precedent qvf.11 itself records applies to it: R3 measured rank")
    print("    CONSTANT 36 through a = 60 after the same instinct said it could")
    print("    not be. Do not let a smooth single unit close that bead.")
    print()
    print("  * NO FREQUENCY IS QUOTED. Neither absolute nor relative, except")
    print("    at the critical points of V5, and there only as the record's own")
    print("    numbers being reproduced.")
    print()
    print("  WHAT THE BEAD'S ACCEPTANCE CRITERION ASKED FOR AND DID NOT GET,")
    print("  disclosed rather than quietly re-scoped:")
    print("  * the range 'a in [-60, 60]' -- superseded. [0, 90] is the")
    print("    fundamental domain under an exact isometry (V7b), and DECISION")
    print("    16 withdrew the admissibility verdict that produced [-60, 60].")
    print("  * a stated verdict of 'stable / unstable / marginal' -- REFUSED")
    print("    as a category error off a critical point, and replaced by the")
    print("    curvature sign. At the critical points, where the two readings")
    print("    coincide, the verdict is given unhedged in V5.")


def v11_scope():
    print()
    print("=" * 78)
    print("V11  WHAT IS NOT COVERED -- read as part of the result")
    print("=" * 78)
    for line in (
        "NO EQUATION OF MOTION WAS INTEGRATED. There is no Floquet analysis,",
        "  no normal variational equation, no trajectory. The dynamical",
        "  question the bead's title names is NOT answered and cannot be",
        "  answered by a curvature sign.",
        "THE SWEEP IS THE PRINTED GRID, not the continuum. Half a degree over",
        "  [0, 90) plus the configurations of record, with the pole band given",
        "  its own coarser treatment. A feature narrower than the grid between",
        "  two same-sign samples would be missed -- and so, specifically, would",
        "  ANY EVEN NUMBER OF CROSSINGS inside one half-degree interval, which",
        "  is the failure mode that matters for a SIGN map: the map reports the",
        "  sign at the samples and the bisection only ever refines brackets the",
        "  map already found. V6b and V7(a-ter) close the two guard bands but",
        "  do so on the SAME lattice and inherit the same blind spot.",
        "THE h WINDOW OF V9 HAS ZERO MARGIN AND IS ABOUT A FACTOR OF TWO WIDE.",
        "  The count is 2 measured against a criterion of 2, and the two",
        "  criteria that do the excluding pull from opposite ends: block",
        "  scalarity (truncation) rules out h = 3e-3 and 1e-3, the a -> 180-a",
        "  isometry (roundoff) rules out h = 3e-3 and 1e-4. The singlet-tangent",
        "  leak is flat to four figures across the whole sweep and the VE",
        "  reproduction passes at every step, so neither constrains anything.",
        "  A future worst combination roughly 6x worse than today's would take",
        "  H_MAIN out of the passing set and redden that row.",
        "V9'S SCALARITY COLUMN IS ONE COMBINATION, NOT THE SWEEP. It is the",
        "  worst combination AT H_MAIN, and it does not stay the worst as h",
        "  moves: an independent whole-gate audit measured the full",
        "  36-combination worst at h = 1e-4 to be 6.654e-04 against this",
        "  probe's 1.259e-05, a factor of 53. The per-threshold h ranges V9",
        "  prints are therefore OPTIMISTIC away from H_MAIN. What is not",
        "  optimistic is the 'all four' window, which the same audit",
        "  reproduced independently by re-running the entire gate.",
        "THE GUARD-BAND COVERAGE IS AT TWO STEP SIZES, NOT THREE, and it stops",
        "  where the method does: within a hundredth of a degree of a = 60, and",
        "  on the STRUT primitive from a = 89.75 upward at h = 3e-4 (89.9 at",
        "  h = 1e-4), the per-block deviation",
        "  exceeds this file's own scalarity criterion. Those rows are printed",
        "  with their deviations and carry SIGNS ONLY. No ratio and no",
        "  magnitude from them enters any claim.",
        "a = 90 AND a = 270 REMAIN UNCLAIMED. The chart's local dimension is 7",
        "  there and its Newton solve fails, measured in V7a. The limit rows",
        "  approaching 90 are 6-D and cannot see the seventh direction, so no",
        "  inertia count at the branch point is offered. That is qvf.3 /",
        "  inviscid-yli.",
        "THE POLE BAND |a - 60| < 1 CARRIES SIGNS ONLY, AND ONLY AT THE",
        "  FINEST STEP MEASURED. Inside it the block scalarity reaches 2.4 --",
        "  larger than the quantity it qualifies -- and V6 measures a COARSE",
        "  step reporting a spurious NEGATIVE triplet there. The signs quoted",
        "  from the band are the h = 1e-4 ones; no ratio from the band is",
        "  quoted at any step size, because nothing in it has converged.",
        "ONLY RAW KERNELS AND THE TWO MASS MODELS OF THE RECORD are swept.",
        "  Normalised kernels are excluded by DECISION 17; strut-axis and",
        "  combo kernels are untested here as in jb_u.",
        "TWO CHARTS, NOT A PROOF OF CHART INDEPENDENCE. V8c re-runs jb_u's",
        "  two-chart comparison on this file's own statistic at seven angles.",
        "  Agreement between two charts is evidence, not a theorem, and no",
        "  third independently constructed chart exists.",
        "THE MOMENTUM-FREE PROJECTION IS CHOSEN BY DERIVATION, not by any",
        "  measurement here or in jb_u: a Euclidean-orthogonal rival with the",
        "  same kernel and the same equivariance reaches the same chart",
        "  agreement while giving different spectra (jb_u U4e). This file",
        "  inherits that limitation whole.",
        "THE RIEMANNIAN HESSIAN ITSELF IS NOT RE-GATED HERE. Its correctness",
        "  rests on jb_u's 17-row gate -- the analytic polar control, the",
        "  metric-compatibility identity, the Christoffel transformation law,",
        "  the nonlinear reparameterisation test with its Gamma := 0 arm. What",
        "  this file gates is what it adds: the group, the representation, the",
        "  projectors, the labelling, the sweep and the verdict.",
        "THE CHARACTER LABELLING ASSUMES THE SYMMETRY GROUP IS T AND NOT",
        "  LARGER. At a = 0 the cuboctahedron has more symmetry than the",
        "  jitterbug does, and the file uses the jitterbug's group throughout,",
        "  which is the correct choice for a statement along the path but",
        "  means the a = 0 row is not resolving the extra structure that exists",
        "  only there. The ranks (1,2,3) are checked at a = 0 as everywhere",
        "  else and hold.",
        "NO MINIMUM-ENERGY PATH. A positive transverse curvature says the",
        "  valley walls rise near the path; it says nothing about how high a",
        "  transverse detour would have to climb, so the memo's open question",
        "  about a route around the a = 60 wall is untouched.",
        "THE TURNOVER ANGLES ARE QUOTABLE TO FOUR DECIMALS, not six. V3d",
        "  measures each root moving with the step size at up to 1.05e-04",
        "  degrees; the bisection tolerance of 1e-07 is set below that so the",
        "  bisection is not the limiting error, and it must not be mistaken",
        "  for the precision of the answer.",
        "THIS FILE MEASURES ONE UNIT AND SAYS ALMOST NOTHING ABOUT AN ARRAY.",
        "  Bead inviscid-qvf.11 reports a physical in-phase array locking",
        "  one-sidedly at the icosahedral phase. V10 records the single unit's",
        "  null result there (chart dimension 6, mass metric non-degenerate,",
        "  every transverse block stiff) as one datum toward that bead's fork.",
        "  The array's own constraint variety is never assembled here and a",
        "  boundary created by the coupling would be invisible to everything",
        "  in this file.",
        "THE a = 90 STRUT POLE (V7c) IS A BY-PRODUCT, measured but not swept:",
        "  its consequences for the strut primitive's landscape are not worked",
        "  out here and deserve their own bead.",
        "ABSOLUTE CURVATURE SCALES REMAIN A CONVENTION (coupling 1, total mass",
        "  1/2, R = 1). Only signs and ratios are measurements, here as",
        "  everywhere in this arc.",
    ):
        print("  * " + line if not line.startswith("  ") else line)


# ==========================================================================
# THE GATE
# ==========================================================================

def gate(v0, v1, v1b, v2, v3b, v3d, v4, v5, v6, v6b, v7, v8, v9,
         n_brackets, n_refused):
    """Every check's verdict in one table, and this process's exit code.

    Every number in this block is COMPUTED from what was passed in. jb_u's
    post-mortem records two rounds of the same defect -- verdicts computed and
    discarded, a reconciliation paragraph narrated from string literals -- and
    the second round happened INSIDE the fix for the first.
    """
    leak, teeth, ok1 = v1
    grad_t, grad_p, grad_anchor, split_teeth, ok1b = v1b
    d5, d1, ranks60, band, ok2 = v2
    _, _, worst_dev, dev_where, inv_pos, gauss_flip, ve_pos = v3b
    root_spread, root_slope, root_dev_deg, ok3d = v3d
    ve_rel, ico_d, second_d, a90_rel, a90_lin = v5
    (pole_dev, pole_sgn_inv, pole_sgn_gau, pole_flagged, pole_nbad,
     pole_dev_h) = v6
    (band60_fine, band60_res_pos, band60_dev, band60_announce, band60_nres,
     band60_nunres, ok6b) = v6b
    (iso, iso_teeth, chart_below_90, chart_fails_90, branch_v, branch_s,
     pole_lin, vertex_floor, band90_fine, band90_res_pos, band90_dev,
     band90_announce) = v7
    mf, gauge, ci, sec, form_used, origin_ranks_ok, sec_min, arms = v8
    rows9, n_good, shape, ok9 = v9
    good9 = [r[0] for r in rows9 if r[4]]
    # The fraction of the fundamental domain the two guard bands remove between
    # them. Computed from the guard constants, so widening either one moves it.
    excl = (2.0 * POLE_GUARD + BRANCH_GUARD) / (A_BRANCH - A_VE)

    checks = [
        (f"V0  group fixes config, all {v0['nswept']} swept a + bands",
         v0["fix"] <= TOL["rep_fix"], f"{v0['fix']:.3e}",
         f"<= {TOL['rep_fix']:.0e}"),
        ("V0  V is group-invariant OFF the path: V(Ty)=V(y)",
         v0["vinv"] <= TOL["rep_iso"], f"{v0['vinv']:.3e}",
         f"<= {TOL['rep_iso']:.0e}"),
        ("V0  the 12 are proper rotations, classes 1+3+8",
         v0["grp_ok"],
         f"{max(v0['det'], v0['orth']):.1e}/{v0['unknown']}", "True"),
        ("V0  representation law rho(g1)rho(g2)=rho(g1g2)",
         v0["mult"] <= TOL["rep_mult"], f"{v0['mult']:.3e}",
         f"<= {TOL['rep_mult']:.0e}"),
        ("V0  rho is an isometry of g",
         v0["iso"] <= TOL["rep_iso"], f"{v0['iso']:.3e}",
         f"<= {TOL['rep_iso']:.0e}"),
        ("V0  projector algebra (sum = I, idempotent)",
         v0["alg"] <= TOL["proj_alg"], f"{v0['alg']:.3e}",
         f"<= {TOL['proj_alg']:.0e}"),
        ("V0  projector eigenvalues are 0 or 1",
         v0["ev"] <= TOL["proj_ev"], f"{v0['ev']:.3e}",
         f"<= {TOL['proj_ev']:.0e}"),
        ("V0  block ranks are (1,2,3)", v0["ranks"] == (1, 2, 3),
         str(v0["ranks"]), "== (1,2,3)"),
        ("V0  ... and the rank test has teeth (V4 gives 3)",
         v0["v4rank"] == TOL["v4_rank"], str(v0["v4rank"]),
         f"== {TOL['v4_rank']}"),
        ("V1  singlet IS the path tangent, worst over sweep",
         leak <= TOL["tangent_leak"], f"{leak:.3e}",
         f"<= {TOL['tangent_leak']:.0e}"),
        ("V1  ... and the same test on a transverse dir",
         teeth >= TOL["tangent_teeth"], f"{teeth:.3e}",
         f">= {TOL['tangent_teeth']:.2f}"),
        ("V1b transverse/path gradient ratio, worst on sweep",
         grad_t <= TOL["grad_transverse"], f"{grad_t:.3e}",
         f"<= {TOL['grad_transverse']:.0e}"),
        ("V1b ... and the PATH part does not vanish",
         grad_p >= TOL["grad_path_teeth"], f"{grad_p:.3e}",
         f">= {TOL['grad_path_teeth']:.0e}"),
        ("V1b ... at one fixed well-behaved row (the icosa)",
         grad_anchor >= TOL["grad_path_teeth"], f"{grad_anchor:.3e}",
         f">= {TOL['grad_path_teeth']:.0e}"),
        ("V1b ... and the SPLIT sees a transverse input",
         split_teeth >= TOL["split_teeth"], f"{split_teeth:.3e}",
         f">= {TOL['split_teeth']:.0e}"),
        ("V2  g at a=60 is [1/32 x5, 5/96] in closed form",
         max(d5, d1) <= 1e-14, f"{max(d5, d1):.3e}", "<= 1e-14"),
        ("V2  character ranks survive a=60 where g does not",
         ranks60 == (1, 2, 3) and bool(band), str(ranks60), "== (1,2,3)"),
        ("V3  per-block scalarity, worst over sweep x 36",
         worst_dev <= TOL["scalarity"], f"{worst_dev:.3e}",
         f"<= {TOL['scalarity']:.0e}"),
        ("V3c every sign bracket is refined, none dropped",
         n_refused == 0 and n_brackets > 0,
         f"{n_brackets - n_refused}/{n_brackets}", "0 dropped"),
        ("V3d turnover roots reproduce across h (degrees)",
         root_spread <= TOL["root_h_spread"], f"{root_spread:.3e}",
         f"<= {TOL['root_h_spread']:.0e}"),
        ("V3d ... and lambda is NOT flat at the root",
         root_slope >= TOL["root_slope_teeth"], f"{root_slope:.3e}",
         f">= {TOL['root_slope_teeth']:.0e}"),
        ("V3d absolute scalarity at the root, as degrees",
         root_dev_deg <= TOL["root_dev_angle"], f"{root_dev_deg:.3e}",
         f"<= {TOL['root_dev_angle']:.0e}"),
        ("V4  curvature ratios square-root onto the record",
         v4 <= TOL["ratio_convention"], f"{v4:.3e}",
         f"<= {TOL['ratio_convention']:.0e}"),
        ("V5  VE reproduction vs jb_s S1, both models",
         ve_rel <= TOL["ve_record"], f"{ve_rel:.3e}",
         f"<= {TOL['ve_record']:.0e}"),
        ("V5  |dV| at the icosahedron vs the record",
         ico_d <= TOL["ico_record"], f"{ico_d:.3e}",
         f"<= {TOL['ico_record']:.0e}"),
        ("V5  recorded second minimum IS critical",
         second_d <= TOL["second_min"], f"{second_d:.3e}",
         f"<= {TOL['second_min']:.0e}"),
        ("V5d a=90 path-gradient collapse is LINEAR",
         a90_lin <= TOL["a90_linear"], f"{a90_lin:.3e}",
         f"<= {TOL['a90_linear']:.0e}"),
        ("V5d ... and its rate reproduces A90_RATE",
         a90_rel <= TOL["a90_rate"], f"{a90_rel:.3e}",
         f"<= {TOL['a90_rate']:.0e}"),
        ("V6  pole band: fine-h sign is ++ (inv pow) / -- (g0.5)",
         pole_sgn_inv == ["++"] and pole_sgn_gau == ["--"],
         f"{','.join(pole_sgn_inv)}|{','.join(pole_sgn_gau)}", "++|--"),
        ("V6  ... and a coarse-h sign error carries dev >= 1",
         pole_nbad > 0 and pole_flagged >= 1.0,
         f"{pole_nbad} rows {pole_flagged:.2e}", ">= 1e+00"),
        ("V6b a=60 band, finest h: all 20 inv-power D,T > 0",
         band60_fine, str(band60_fine), "True"),
        ("V6b ... rows that RESOLVE are positive too",
         band60_res_pos and band60_dev <= TOL["band_dev"],
         f"{band60_dev:.3e} ({band60_nres}/{band60_nres + band60_nunres})",
         f"<= {TOL['band_dev']:.0e}"),
        ("V6b ... and a coarse-h negative announces itself",
         band60_announce >= TOL["band_sign_announce"],
         f"{band60_announce:.2e}", f">= {TOL['band_sign_announce']:.0e}"),
        ("V7  a>89 band, finest h: all 20 inv-power D,T > 0",
         band90_fine, str(band90_fine), "True"),
        ("V7  ... rows that RESOLVE are positive too",
         band90_res_pos and band90_dev <= TOL["band_dev"],
         f"{band90_dev:.3e}", f"<= {TOL['band_dev']:.0e}"),
        ("V7  ... and a coarse-h negative announces itself",
         band90_announce >= TOL["band_sign_announce"],
         f"{band90_announce:.2e}", f">= {TOL['band_sign_announce']:.0e}"),
        ("V6/V7 guards exclude a small part of the domain",
         excl <= TOL["guard_fraction"], f"{excl:.3e}",
         f"<= {TOL['guard_fraction']:.2f}"),
        ("V7  a>89 guard priced: raw VERTEX dev stays low",
         branch_v <= TOL["branch_vertex"], f"{branch_v:.3e}",
         f"<= {TOL['branch_vertex']:.0e}"),
        ("V7  ... and the raw STRUT dev does NOT",
         branch_s >= TOL["branch_strut_teeth"], f"{branch_s:.3e}",
         f">= {TOL['branch_strut_teeth']:.0e}"),
        ("V7  strut-midpoint pole law is 2(90-a) at 89.99",
         pole_lin <= TOL["pole_linearity"], f"{pole_lin:.3e}",
         f"<= {TOL['pole_linearity']:.0e}"),
        ("V7  ... and the VERTEX separation does NOT collapse",
         vertex_floor >= TOL["pole_vertex_floor"], f"{vertex_floor:.4f}",
         f">= {TOL['pole_vertex_floor']:.2f}"),
        ("V7  isometry a -> 180-a maps blocks to blocks",
         iso <= TOL["isometry"], f"{iso:.3e}", f"<= {TOL['isometry']:.0e}"),
        ("V7  ... and a mismatched pair does NOT",
         iso_teeth >= TOL["isometry_teeth"], f"{iso_teeth:.3e}",
         f">= {TOL['isometry_teeth']:.0e}"),
        ("V7  chart works below 90 and FAILS at 90",
         chart_below_90 and chart_fails_90,
         f"{chart_below_90}/{chart_fails_90}", "True/True"),
        ("V8  section vs momentum-free Gamma over sweep",
         mf <= TOL["metric_form"], f"{mf:.3e}",
         f"<= {TOL['metric_form']:.0e}"),
        ("V8  ... and the two Gamma arms ARE different forms",
         arms >= TOL["form_in_use"], f"{arms:.3e}",
         f">= {TOL['form_in_use']:.0e}"),
        ("V8  the sweep's form IS momentum-free, not section",
         form_used >= TOL["form_in_use"], f"{form_used:.3e}",
         f">= {TOL['form_in_use']:.0e}"),
        ("V8  closed-form gauge leak |Z^T W D| on path",
         gauge <= TOL["gauge_leak"], f"{gauge:.3e}",
         f"<= {TOL['gauge_leak']:.0e}"),
        ("V8  two charts agree, momentum-free",
         ci <= TOL["chart_agree"], f"{ci:.3e}",
         f"<= {TOL['chart_agree']:.0e}"),
        ("V8  ... and the SECTION form must still fail it",
         sec >= TOL["chart_teeth"], f"{sec:.3e}",
         f">= {TOL['chart_teeth']:.0e}"),
        ("V8  character ranks (1,2,3) in the ORIGIN chart",
         origin_ranks_ok, str(origin_ranks_ok), "True"),
        ("V9  h window: >= 2 step sizes, and not flat", ok9,
         f"{n_good}/{shape:.1e}",
         f">= {TOL['h_window']}/{TOL['h_shape']:.0e}"),
        ("V9  ... and H_MAIN is one of the passing step sizes",
         H_MAIN in good9, str(H_MAIN in good9), "True"),
        ("VERDICT  all 20 inverse-power combos: D,T > 0",
         inv_pos, str(inv_pos), "True"),
        ("VERDICT  all 16 gaussian combos DO change sign",
         gauss_flip, str(gauss_flip), "True"),
        ("VERDICT  at the ground state all 36 have D,T > 0",
         ve_pos, str(ve_pos), "True"),
    ]
    print()
    print("=" * 78)
    print(f"GATE  {len(checks)} rows: every check's verdict, and this "
          f"process's exit code")
    print("=" * 78)
    for name, passed, val, crit in checks:
        print(f"  {'PASS' if passed else 'FAIL'}  {name:48s} "
              f"{val:>16s} {crit:>14s}")

    print()
    print("  THE LAST THREE ROWS ARE THE DELIVERABLE, and they are gated for")
    print("  the same reason as everything above it: a verdict that is printed")
    print("  but not asserted is a verdict nobody can break. The second of the")
    print("  three is the first's non-vacuity -- without it, 'no negative")
    print("  transverse curvature was found' would be satisfied by a build")
    print("  incapable of producing one.")
    print()
    print("  STEP SIZE, RECONCILED from V9's own rows rather than asserted")
    print("  beside them:")
    print(f"    The sweep is quoted at h = {H_MAIN:.0e}.")
    print(f"    Step sizes meeting all three swept thresholds: "
          + (", ".join(f"{x:.0e}" for x in good9) if good9 else "NONE"))
    print(f"    Is the quoted h among them: {H_MAIN in good9}"
          f"   <- GATED, row 'V9 ... and H_MAIN is one of the passing")
    print("       step sizes'. It was previously printed as a computed boolean")
    print("       and asserted nowhere, so a build in which H_MAIN fell OUT of")
    print("       the passing set while two other step sizes stayed in would")
    print("       have satisfied the h-window row and exited 0 while quoting")
    print("       every block value at a step the file's own criterion rejects.")
    print("    Only ONE of the three swept thresholds is h-sensitive at all")
    print("    (V9's per-column dynamic ranges), and the window's count of")
    print(f"    {n_good} sits at ZERO MARGIN against a criterion of "
          f"{TOL['h_window']}.")
    if rows9:
        by_h = {r[0]: r[1] for r in rows9}
        if H_MAIN in by_h:
            best = min(by_h, key=by_h.get)
            # Guarded for the same reason as V9's shape statistic: a mutation
            # that silences the per-block deviation makes every entry exactly
            # 0.0, and this ratio killed the run AFTER the gate table had
            # printed but BEFORE its exit code -- a complete table followed by a
            # traceback is still a destroyed verdict.
            _ratio = (by_h[H_MAIN] / by_h[best]) if by_h[best] > 0 \
                else float("nan")
            print(f"    Best scalarity in the sweep is at h = {best:.0e}; the")
            print(f"    ratio scalarity({H_MAIN:.0e})/scalarity({best:.0e}) = "
                  f"{_ratio:.2f}, which is the honest error")
            print("    bar on the quoted block values.")
    print()
    print("  THE TWO GUARD BANDS, priced from this run's own numbers rather")
    print("  than asserted, AND PRICED LIKE FOR LIKE. The previous version of")
    print("  this block divided an INSIDE number taken over all step sizes in")
    print("  V6 (attained at h = 1e-3) by an OUTSIDE number at H_MAIN, so the")
    print("  quoted ratio compared two different discretisations. Both step")
    print("  sizes are now shown and the ratio is taken at H_MAIN.")
    print(f"    worst per-block scalarity INSIDE |a - 60| < {POLE_GUARD} "
          f"(raw inverse power)")
    print(f"        at H_MAIN = {H_MAIN:.0e}          = {pole_dev_h:.2e}")
    print(f"        over every step size in V6  = {pole_dev:.2e}"
          f"   (attained at h = 1e-03)")
    print(f"    worst per-block scalarity OUTSIDE both guards, at H_MAIN"
          f"  = {worst_dev:.2e}")
    print(f"      at a = {dev_where[0]:.4f}, {dev_where[1]}, {dev_where[2]}, "
          f"{dev_where[3]}")
    _pratio = (pole_dev_h / worst_dev) if worst_dev > 0 else float("nan")
    print(f"    LIKE-FOR-LIKE ratio at H_MAIN = {_pratio:.0f}x."
          f"  Decisive, and honest.")
    print("    The band is reported separately in V6 and its reading is taken")
    print("    at the FINEST step there, not at H_MAIN: V6 measures a COARSE")
    print("    step producing a spurious sign inside the band (deviation 2.4,")
    print("    larger than the quantity it qualifies), so 'a large deviation")
    print("    cannot move a sign' is FALSE and is not the justification for the")
    print("    guard. The justification is that inside the band nothing has")
    print("    converged -- which is why no RATIO from it is quoted at any step")
    print("    size.")
    print(f"    THE GUARD IS APPLIED BY KERNEL FAMILY, not by angle: the")
    print("    Gaussians do not diverge at a = 60 and V3c bisects them straight")
    print("    through the band. Two turnovers that previously came back")
    print("    'STRADDLES a=60 / not refined' are located there.")
    print(f"    The a = 90 guard is NOW priced the same way -- V7(a-bis)'s two")
    print(f"    numbers are gate rows, VERTEX {branch_v:.2e} against STRUT "
          f"{branch_s:.2e}.")
    print("    They were previously computed, printed and returned as neither,")
    print("    while this paragraph claimed the two bands were priced alike.")
    print(f"    THE GUARDS ARE NOW BOUNDED FROM ABOVE. Together they remove")
    print(f"    {excl:.2%} of the fundamental domain"
          f"   (criterion <= {TOL['guard_fraction']:.0%})")
    print("    Both were measured to be unconstrained in that direction:")
    print("    POLE_GUARD = 3 and BRANCH_GUARD = 5 each exit 0 with zero red")
    print("    rows, AND widening IMPROVES the reported worst scalarity, so the")
    print("    gate as it stood rewarded hiding more of the domain. There is no")
    print("    derivation for a particular width -- a pole is a point and any")
    print("    positive width is a fit -- so what is asserted is that the fit")
    print("    stays small. It is labelled a POLICY bound, not a derived one.")
    print("    COVERAGE INSIDE BOTH BANDS is now measured rather than assumed:")
    print(f"    V6b runs all 20 inverse-power combinations across |a-60| < "
          f"{POLE_GUARD}")
    print(f"    and V7(a-ter) across a > {90 - BRANCH_GUARD:.0f}, both at two "
          f"step sizes, and both")
    print("    report BOTH transverse blocks positive on every row.")

    failed = [n for n, p, _, _ in checks if not p]
    print()
    if failed:
        print(f"  !! {len(failed)} CHECK(S) FAILED -- this is a bug report, not")
        print("     a measurement. Nothing above may enter the record.")
        for n in failed:
            print(f"       - {n}")
        return 1
    print("  ALL CHECKS PASSED.")
    print("  Reminder for whoever records these numbers: FOUR declarations")
    print("  (kernel, mass model, primitive, METRIC FORM); off a critical point")
    print("  these are CURVATURE SCALES, not frequencies; and the sign")
    print("  structure is NOT invariant across the kernel.")
    return 0


def main():
    np.set_printoptions(precision=6, suppress=False, linewidth=170)

    grid = main_grid()
    print("=" * 78)
    print("jb_v -- TRANSVERSE CURVATURE OF V ALONG THE SYMMETRIC PATH")
    print("        bead inviscid-qvf.4")
    print("=" * 78)
    print(f"  fundamental domain [0, 90), {len(grid)} angles, step size "
          f"h = {H_MAIN:.0e}")
    print(f"  pole guard: |a - 60| >= {POLE_GUARD} in the main sweep; V6/V6b "
          f"cover the band,")
    print(f"              and the guard is applied BY KERNEL FAMILY -- the")
    print(f"              Gaussians are finite at a = 60 and are bisected "
          f"through it")
    print(f"  metric form: MOMENTUM-FREE throughout (qvf.9 corollary iii)")
    print(f"  {len(COMBOS)} combinations = {len(RAW_KERNELS)} raw kernels x 2 "
          f"primitives x {len(MODELS)} mass models")

    # SPECULATIVE PARALLEL PREFETCH. Every `Site` is an independent pure
    # construction, so the previous run's recorded argument trace -- the main
    # grid plus every angle V3c's bisections and V9's step sweep reached --
    # is replayed through a process pool before the serial pass below. Order
    # and output are untouched; the pass simply finds the work already done.
    jb_cache.prefetch(_site_build)

    sites = []
    unmeasurable = []
    for a in grid:
        try:
            sites.append(site(a))
        except (ChartUnmeasurable, IrrepLabelError) as exc:
            unmeasurable.append((a, type(exc).__name__,
                                 str(exc).splitlines()[0]))
    if unmeasurable:
        print(f"\n  {len(unmeasurable)} grid angle(s) NOT MEASURABLE, reported "
              f"rather than dropped:")
        for a, kind, msg in unmeasurable:
            print(f"    a = {a:.6f}  {kind}: {msg[:100]}")
    else:
        print(f"  every grid angle measurable: {len(sites)} sites built")

    r0 = v0_symmetry_controls(sites)
    r1 = v1_singlet_is_the_path_tangent(sites)
    r1b = v1b_gradient_has_no_transverse_part(sites)
    r2 = v2_why_not_the_mass_metric()
    v3_sweep(sites)
    r3b = v3b_sign_map(sites)
    roots, n_brackets, n_refused = v3c_refine_flips(r3b[1])
    r3d = v3d_root_reproducibility(roots)
    r4 = v4_ratios(sites)
    r5 = v5_critical_points(sites)
    r6 = v6_pole()
    r6b = v6b_band_coverage()
    r7 = v7_branch_point_and_isometry(sites)
    r8 = v8_metric_form_and_second_chart(sites)
    r9 = v9_step_sweep(r3b[3])
    v10_verdict(r3b[4], r3b[5], r3b[6], roots, r6b[0], r7[8],
                max(r3d[0], r3d[2]))
    v11_scope()
    return gate(r0, r1, r1b, r2, r3b, r3d, r4, r5, r6, r6b, r7, r8, r9,
                n_brackets, n_refused)


if __name__ == "__main__":
    # `--no-cache` / `--clear-cache` are consumed here; anything else is a
    # loud failure rather than a run that silently ignored what was asked.
    _rest = jb_cache.parse_argv(sys.argv[1:])
    if _rest:
        print(f"unrecognised argument(s): {' '.join(_rest)}", file=sys.stderr)
        print("usage: jb_v_transverse_curvature.py [--no-cache] [--clear-cache]",
              file=sys.stderr)
        sys.exit(2)
    sys.exit(main())
