/**
 * Copyright (c) 2016 Chiral Behaviors, LLC, all rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.chiralbehaviors.inviscid.automaton.measure;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.automaton.Necronomata;
import com.chiralbehaviors.inviscid.automaton.lga.CollisionSweep;
import com.chiralbehaviors.inviscid.automaton.lga.ContactAtlasGenerator;
import com.chiralbehaviors.inviscid.automaton.lga.ContactPredicate;
import com.chiralbehaviors.inviscid.automaton.lga.ContactScan;
import com.chiralbehaviors.inviscid.automaton.lga.FccNeighborhood;
import com.chiralbehaviors.inviscid.automaton.lga.HybridAutomaton;
import com.chiralbehaviors.inviscid.automaton.lga.MemberGeometry;
import com.chiralbehaviors.inviscid.automaton.lga.QuantaExchangeRule;

/**
 * B.5 (bead inviscid-0nx.10): the isotropy discriminator. Measures, per
 * probe direction {@code d} in {@link StructureFactor.Direction} ({@code
 * X100}/{@code X110}/{@code X111}), TWO INDEPENDENT estimators of
 * transport magnitude along {@code d}, on the SAME underlying
 * configuration (one automaton run per seed), and reports the anisotropy
 * ratio {@code A = max_d/min_d} for each -- both as a naive per-seed
 * statistic AND, per the stacked-review correction below, as a properly
 * pooled/null-calibrated statistic.
 *
 * <h2>The cardinal constraint (do not weaken this class's honesty
 * contract)</h2>
 * This class MEASURES; it does not decide an isotropy posture. Every
 * degenerate/no-signal condition (see "Non-fabrication contract" below)
 * MUST surface as {@link OptionalDouble#empty()}, never a silently
 * fabricated {@code 1.0} or any other number. Disagreement between the
 * two estimators is reported as-is (both {@link EstimatorResult}s live
 * side by side in {@link SeedResult}) -- this class never averages them
 * together.
 *
 * <h2>STACKED-REVIEW CORRECTION (T2 {@code
 * inviscid/critique-anisotropy-probe-inviscid-0nx10.md}, T3 {@code
 * critique-pattern-max-min-ratio-order-statistic-bias}) -- read before
 * trusting {@link #bootstrapCi}</h2>
 * The FIRST version of this class summarized a Phase A campaign as "A =
 * 1.293, 95% CI [1.179,1.415], excludes 1.0" using {@link #bootstrapCi}
 * over the list of PER-SEED {@code max_d/min_d} ratios (a "mean of
 * ratios" statistic). That statistic is bounded BELOW by exactly 1.0 (max
 * &gt;= min always, by construction of a ratio of the same three
 * quantities) and is systematically upward-biased by seed-to-seed noise
 * -- an order-statistic artifact, not evidence of physical anisotropy.
 * The seed-pooled "ratio of means" on the SAME campaign data (average
 * each direction's magnitude across all 8 seeds FIRST, then take
 * max/min) was only 1.063, a 6x smaller effect, fully consistent with a
 * back-of-envelope order-statistic null check (E[max]/E[min] ~ 1.4 for 3
 * iid samples at the campaign's ~20% per-direction CV, from noise
 * alone). The "winning" (max) direction also flipped essentially at
 * random across seeds -- the signature of noise, not a stable
 * crystallographic axis.
 *
 * <p>{@link #bootstrapCi} and the raw per-seed-ratio list therefore
 * remain in this class as a DIAGNOSTIC (useful for eyeballing per-seed
 * spread), but are NOT the significance statistic -- {@link
 * #pooledEstimate} is. It computes (a) the seed-pooled ratio-of-means
 * with a proper resample-then-aggregate bootstrap CI (resample SEED
 * INDICES with replacement, recompute per-direction MEANS from the
 * resample, THEN take max/min -- never resample the pre-collapsed ratio
 * list), and (b) a permutation/null-calibration test: shuffle which
 * magnitude is labeled X100/X110/X111 WITHIN each seed (preserving each
 * seed's own noise realization, destroying only the direction-label
 * information), recompute the pooled ratio-of-means under the shuffle,
 * repeat {@link #PERMUTATION_COUNT} times, and report the empirical
 * p-value (fraction of permuted statistics &gt;= the observed one) plus
 * the null distribution's 95th percentile for context. A small p-value
 * is the actual evidence a stable, direction-linked effect exists; "CI
 * excludes 1.0" on the naive per-seed-ratio statistic is not.
 *
 * <h2>TRANSPORT estimator -- exact definition</h2>
 * "Quanta distribution" here means: {@code Necronomata.frequency} is the
 * conserved per-member quanta count (design memo, "Structural insight");
 * a spatially LOCALIZED excess of quanta (a "packet") seeded at one
 * lattice cell is the transport probe. The estimator's input is a
 * per-cell scalar field of quanta ({@link
 * StructureFactor#coarseGrainedField(com.chiralbehaviors.inviscid.automaton.QuantaField)},
 * reused directly), and it weights each cell by {@code Math.abs(...)} of
 * that cell's value.
 *
 * <h3>The discrete maximum principle is RETIRED, not re-attributed
 * (beads inviscid-0nx.28 / inviscid-o24; two earlier versions of this
 * javadoc were wrong about it in two different ways)</h3>
 * {@link QuantaExchangeRule} moves ONE quantum from a strictly-higher
 * member to a strictly-lower one. Exactly ONE invariant survives that,
 * and it holds for EVERY initial condition, signed or not:
 * <ul>
 * <li><b>(I1) exact integer conservation.</b> Every firing is
 * {@code (-1,+1)} and every {@code CollisionRule.Delta} is zero-sum by
 * construction, so the lattice-wide member total is constant however many
 * contacts resolve in a tick. Tick totals are read from {@link
 * ConservationAudit#ledger()}, never hand-derived.</li>
 * </ul>
 *
 * <p><b>(I2) interval invariance -- "the maximum principle" -- is FALSE
 * for this automaton, under every initial condition.</b> The proof that
 * used to be given here ("a firing requires {@code q_a >= q_b + 1};
 * afterwards {@code q_a - 1 >= q_b >= m0} ... and no other member is
 * touched") is a proof for SEQUENTIAL single-firing dynamics. {@link
 * CollisionSweep} does not implement sequential firing. It implements
 * SNAPSHOT resolution (bead inviscid-72s, that class's own javadoc):
 * every contact in a tick resolves against the FROZEN pre-tick
 * {@code frequency} array, and the resulting deltas ACCUMULATE additively
 * ({@code deltaF[i] += delta.deltaX()}), so multiple contacts touching the
 * same member in one tick ALL contribute. "No other member is touched" is
 * exactly the hypothesis snapshot resolution denies.
 *
 * <p>Counterexample, in one tick: member {@code a} at {@code q_a = 1} with
 * two same-tick partners {@code b}, {@code c} at {@code 0}. Both contacts
 * see the frozen {@code q_a = 1 > 0} and both fire, so {@code a}
 * accumulates {@code -1} twice and lands at {@code q_a = -1}, below
 * {@code m0 = 0}. The mirror construction pushes a member above
 * {@code M0}. This is not an edge case being tolerated: {@link
 * CollisionSweep} records that {@code ~18.8%} of ticks exhibit a
 * same-member multi-contact at typical seed densities, and that quanta are
 * signed {@code long}s with no floor, so "overdrawing a member below zero
 * is LEGAL BY DESIGN -- there was never anything to protect against". The
 * principle therefore never held here; the first version of this javadoc
 * attributed a false claim to the {@code {0, packetQuanta}} initial
 * condition and the second promoted the same false claim to universality.
 * Both are withdrawn. Every downstream claim that rested on {@code (I2)}
 * is withdrawn with it:
 * <ul>
 * <li>"the coarse-grained per-cell field is therefore always
 * non-negative" -- withdrawn as a THEOREM on every path, not merely under
 * a signed background (see the measured status of the
 * {@code {0, packetQuanta}} path immediately below);</li>
 * <li>"the field literally IS the quanta-deviation from the zero
 * background -- no separate background subtraction is needed" -- FALSE
 * under a signed background, and it is exactly the claim §D-A's
 * matched-pair subtraction replaces. See "MATCHED-PAIR transport"
 * below;</li>
 * <li>"invariant interval {@code [-q, max(q, packetQuanta)]} under a
 * signed background" -- withdrawn; there is no invariant interval, which
 * bears on any absorbing-state / {@code rho_c} characterization built on
 * T2 {@code design-seeding-radius.md} §D-A and must not be assumed
 * there.</li>
 * </ul>
 *
 * <p><b>The {@code {0, packetQuanta}} path: MEASURED, not proved.</b>
 * Because {@code (I2)} is false, nothing guarantees the Phase A field is
 * non-negative, and so nothing guarantees {@code Math.abs} was the
 * identity on it. That question is settled by census rather than by
 * theorem. Running the exact Phase A campaign parameters (8 seeds
 * {@code 42..49}, {@code 128} ticks, {@code 8^3}, the
 * {@code {0, packetQuanta=30}} initial condition) and counting at every
 * one of the 1024 snapshots gives <b>zero</b> negative MEMBERS and
 * <b>zero</b> negative CELLS; the minimum member value ever observed is
 * {@code 0} and the minimum cell value is {@code 0.0}. {@code Math.abs}
 * WAS the identity on every Phase A snapshot, the committed {@code
 * anisotropy-report-phaseA.tsv} weights no negative cell positively, and
 * every Phase A number is bit-for-bit what it was. The mechanism is the
 * sparsity of the regime, not an invariant: driving a member negative
 * needs more same-tick firing contacts on one member than it holds quanta,
 * and at Phase A's contact density the frontier members that sit at
 * {@code q = 1} are effectively never multiply contacted in one tick.
 * <b>This is an empirical property of these parameters and does not
 * extend to longer or denser runs.</b>
 *
 * <p><b>What is pinned in the suite, and how little of the campaign that
 * is.</b> The full census costs minutes, so {@code
 * AnisotropyProbeTest#phaseARegimeIsNonNegativeAndFarFromSaturationByCensus}
 * runs it at reduced scope -- seeds {@code {42,43}} x ticks {@code 0..23}
 * only -- which is a TRIPWIRE ON THE REGIME, not a proof about the
 * campaign, and its coverage is small and BIASED LOW. Measured over the
 * full campaign: the window is 48 of 1024 snapshots ({@code 4.69%}), 1482
 * of 39869 occupied ({@code q>0}) member observations ({@code 3.72%}), and
 * -- the number that matters, since only a member at {@code q == 1} can be
 * driven negative by a single extra same-tick contact -- just 42 of 5842
 * at-risk observations, {@code 0.72%}. Worse, that 0.72% is not a uniform
 * sample: the at-risk population by 16-tick bucket across all 8 seeds runs
 * {@code 159 / 471 / 875 / 1012 / 1013 / 875 / 738 / 699}, so the tripwire
 * sits in the FIRST bucket, at the very bottom of a risk curve that
 * plateaus ~6x higher at ticks 48-79. <b>The window samples the regime
 * where negativity is LEAST likely.</b> A parameter change that first
 * drives a member negative after tick 24, or on seeds 44-49, passes it
 * silently. The full census is reproducible from this tree on demand --
 * {@link #main} with the {@code --census} argument -- so the headline
 * number can be re-derived rather than taken on trust.
 *
 * <p>Per tick {@code t} and direction {@code d}, the estimator computes
 * the mass-weighted SECOND MOMENT of displacement from a FIXED reference
 * point (the packet's own seed cell, never a per-tick recomputed
 * centroid):
 * <pre>
 *   M_d(t) = sum_cell |field(cell,t)| * proj_d(cell - origin)^2
 *            / sum_cell |field(cell,t)|
 * </pre>
 * where {@code proj_d} matches {@link StructureFactor}'s own
 * cross-direction-comparable convention ({@code X100 -> dx}, {@code
 * X110 -> (dx+dy)/sqrt(2)}, {@code X111 -> (dx+dy+dz)/sqrt(3)}) so both
 * estimators treat "distance along [110]" identically -- required for
 * their disagreement to be a meaningful finding rather than an artifact
 * of inconsistent axis conventions. Using a FIXED origin (a
 * mean-squared-displacement definition), not a re-centered variance, is
 * deliberate: a packet that drifts ballistically along one axis without
 * spreading would report LOW variance-about-its-own-centroid despite
 * being a striking anisotropy signal; MSD-from-origin captures both
 * drift and spread and is harder to fool into reporting isotropy (the
 * bead's "falsify, not illustrate" instruction).
 *
 * <p>{@code D_hat(d) = |OLS slope of M_d(t) against t|} (group
 * transport-rate estimate; the OLS machinery mirrors {@link
 * StructureFactor#extractRidge}'s pattern applied to (tick, moment)
 * pairs instead of (k, omega) pairs).
 *
 * <h2>MATCHED-PAIR transport under a signed background (bead
 * inviscid-0nx.28; T2 {@code design-seeding-radius.md} §D-A)</h2>
 * §D-A anchors the decontaminated transport campaign ABOVE the absorbing
 * transition on a uniform signed background, and extracts the signal by
 * MATCHED PAIRS: two runs at the same seed and therefore the same RNG
 * background, one with the packet and one without, subtracted. Write
 * {@code f_P} and {@code f_C} for the two coarse-grained fields and
 * {@code D = f_P - f_C} for the difference. <b>{@code D}, not
 * {@code f_P}, is the estimator's input</b> ({@link
 * #transportEstimateMatchedPair}); every statement in this section is
 * about {@code D}.
 *
 * <h3>What survives the subtraction, and what does not</h3>
 * <ul>
 * <li><b>(I1) survives, and becomes the load-bearing invariant.</b> Each
 * run conserves its own member total exactly, so
 * {@code S := sum_cell D(cell,t) = T_P - T_C} equals the INJECTED EXCESS
 * ({@code 30 * packetQuanta}) and is constant for all {@code t}, exactly,
 * as an integer. {@link #transportEstimateMatchedPair} asserts that
 * tick-invariance: besides being the physics, it is the sharpest
 * available check that the caller really did run both halves of the pair
 * on the SAME background.</li>
 * <li><b>Nothing bounds {@code D}'s values.</b> There is no {@code (I2)}
 * to lose -- it was never true (above) -- and subtraction would destroy it
 * even if it had been: the exchange rule is nonlinear (both whether a
 * firing occurs and its direction depend on {@code sign(q_a - q_b)},
 * which the packet perturbs), so the two runs take genuinely different
 * firing sequences and evolution does not commute with subtraction.
 * {@code D} is genuinely signed: one quantum the packet run moved and the
 * control run did not appears as {@code +1} at the destination cell and
 * {@code -1} at the source.</li>
 * <li><b>No upper bound on {@code ||D||_1} survives.</b> At tick 0
 * {@code D} is {@code +S} at the origin cell and exactly zero elsewhere,
 * so {@code ||D||_1 = |S|}; thereafter it GROWS as the two runs
 * decorrelate, bounded only trivially by
 * {@code ||f_P||_1 + ||f_C||_1}. Localization of {@code D} is an
 * EMPIRICAL, early-time property, not a theorem -- which is precisely why
 * the wrap-safety precondition below has to actually bite.</li>
 * </ul>
 *
 * <h3>What replaces the maximum principle as the licence for
 * {@code Math.abs} weighting</h3>
 * Not another sign guarantee: conservation plus the triangle inequality.
 * Since {@code |S| = |sum_cell D| <= sum_cell |D| = ||D||_1} and
 * {@code S} is tick-invariant,
 * <pre>
 *   ||D(.,t)||_1 &gt;= |S| = 30 * packetQuanta &gt; 0   for every tick t
 * </pre>
 * so {@code mu_t(cell) := |D(cell,t)| / ||D(.,t)||_1} is a well-defined
 * PROBABILITY MEASURE on cells at every tick, and {@code M_d(t)} is
 * exactly its second moment about the fixed origin. The moment definition
 * never required its weights to come from a non-negative PHYSICAL field;
 * it required them non-negative and summable with a strictly positive
 * total, which {@code |D|} is by construction plus the bound above.
 * {@code Math.abs}, which the first version of this javadoc described as
 * merely "defensive", is therefore DEFINITIONAL on the matched-pair path.
 *
 * <p>Total variation is the right weight, not merely a convenient one:
 * <ul>
 * <li><b>the signed alternative is degenerate.</b>
 * {@code sum_cell D * proj_d^2} is not a spread -- it cancels precisely
 * in the case of interest, because the elementary transport event here is
 * a single hop contributing {@code +1} and {@code -1} at two cells, which
 * the signed form records as a DIFFERENCE of squared displacements while
 * total variation records it as two units of response at two places: the
 * correct accounting of "where the perturbation has reached";</li>
 * <li><b>it reduces EXACTLY to the pre-existing estimator.</b> On the
 * Phase A zero-background substrate the control run holds no quanta at
 * all, so no member is ever strictly greater than another, nothing ever
 * transfers, {@code f_C == 0} identically, and {@code D == f_P >= 0}.
 * The matched-pair estimator on Phase A data is the SAME numbers, not
 * merely a compatible estimator -- pinned by {@code
 * AnisotropyProbeTest#matchedPairOnZeroBackgroundReproducesSingleFieldEstimate};</li>
 * <li><b>it is the norm the wrap-safety precondition is already stated
 * in,</b> so weight and precondition speak about one quantity rather than
 * two.</li>
 * </ul>
 *
 * <h3>Degeneracy of the ratio-of-means under signed cancellation</h3>
 * Two channels, which must not be conflated:
 * <ul>
 * <li><b>Normalization degeneracy</b> ({@code secondMoment}'s
 * {@code totalMass > 0} branch) CANNOT fire for a matched pair with a
 * nonzero injected excess: by the bound above, signed cancellation can
 * shrink {@code ||D||_1} at most down to {@code |S|}, never to zero. The
 * only way to reach it is {@code S == 0} -- subtracting a run from
 * itself -- where {@code D == 0} identically, every moment is
 * {@code 0.0}, every slope is exactly {@code 0.0}, and {@link
 * #ratio(Map)} correctly returns {@link OptionalDouble#empty()}. {@link
 * #transportEstimateMatchedPair} rejects a zero injected excess up front
 * rather than returning that vacuous empty.</li>
 * <li><b>Ratio degeneracy</b> ({@link #RATIO_DEGENERATE_EPSILON}) is
 * applied to the per-direction OLS SLOPE, and to pooled means of slopes,
 * never to {@code totalMass} -- it is dimensionally unrelated to L1
 * cancellation and needs no change. What DOES change is its DIAGNOSIS. In
 * the {@code {0, packetQuanta}} regime an exactly-zero slope meant "no
 * collision ever fired" (the K=0 baseline). On the matched-pair path it
 * can ALSO mean "the response spread, but with no net second-moment
 * trend". Reporting empty stays correct and non-fabricating either way;
 * telling the two apart is what {@link
 * MatchedPairTransport#responseL1First()} / {@link
 * MatchedPairTransport#responseL1Last()} are for --
 * {@code responseL1Last == |S|} means the response never left the seed
 * cell, whereas {@code responseL1Last >> |S|} with a flat slope means it
 * spread without a trend.</li>
 * </ul>
 *
 * <p>None of this touches significance: {@link
 * PooledResult#permutationPValue()} remains the authoritative statistic,
 * and no CI-excludes-1.0 reasoning is reintroduced anywhere.
 *
 * <h2>The periodic-wrap precondition: two limbs, one of them provably
 * vacuous on campaign geometry (FIX 6, then re-derived by bead
 * inviscid-0nx.28)</h2>
 * <b>Which field.</b> The precondition applies to the field the moment is
 * actually computed from: the raw quanta field on the
 * {@code {0, packetQuanta}} path, and the DIFFERENCE field {@code D} on
 * the matched-pair path -- never to a raw signed-background field.
 * Asserting localization on a raw signed field is not merely weak, it is
 * meaningless: such a field has full support by construction, so the
 * criterion is either satisfied vacuously (any campaign origin -- see
 * W1's unreachability below) or violated by background quanta that have
 * never moved (an off-center origin). The precondition is a statement
 * about the RESPONSE's reach, and the response is {@code D}.
 *
 * <h3>Limb W1 -- exactness ({@link #assertWrapSafe}, criterion
 * unchanged)</h3>
 * {@code proj_d} above is an UNWRAPPED Cartesian displacement -- valid
 * only while the naive {@code |coord-origin|} still equals the TRUE
 * periodic minimum-image distance {@code min(|coord-origin|,
 * extentAxis-|coord-origin|)}. Those two agree EXACTLY when
 * {@code |coord-origin| <= extentAxis/2} (both paths are the same length
 * at the tie point {@code extentAxis/2}) and diverge -- the naive value
 * silently overestimates the true distance -- once
 * {@code |coord-origin| > extentAxis/2}. {@link #assertWrapSafe}
 * enforces exactly this ORIGIN-RELATIVE criterion, per axis, for every
 * snapshot: any mass whose per-axis displacement from the SUPPLIED
 * {@code originCell} exceeds {@code extentAxis/2} trips {@link
 * IllegalStateException}. This is sound for ANY {@code originCell}
 * {@link #transportEstimate} is called with -- the first version of this
 * guard checked literal array-boundary coordinates ({@code 0} or
 * {@code extent-1}) regardless of {@code originCell}, which is only a
 * correct proxy for the centered case and silently under-protects an
 * off-center origin.
 *
 * <p><b>W1 is PROVABLY UNREACHABLE for every legal campaign geometry, on
 * every field -- signed or not, localized or not.</b> {@link
 * FccNeighborhood} requires every extent axis to be EVEN (and at least
 * 4), and {@link #nearestEvenParityCenter} places the origin at
 * {@code L/2}, or on {@code z} only at {@code L/2 - 1}; in BOTH cases
 * {@code max_i |i - origin| = max(origin, L-1-origin) = L/2} exactly,
 * while the criterion is a STRICT {@code >}. (Enumerated over all legal
 * extents with each axis even in {@code 4..16}, non-cubic included, in
 * {@code
 * AnisotropyProbeTest#wrapSafetyLimbW1IsStructurallyUnreachableOnEveryCampaignGeometry}.)
 * So W1 is a real, reachable backstop for the general any-origin public
 * API -- see {@code
 * AnisotropyProbeTest#transportEstimateFailsLoudlyWhenMassExceedsHalfPeriodFromOrigin}
 * -- and a STRUCTURAL no-op for the campaign.
 *
 * <p><b>Provenance, stated accurately: this is a RE-EVALUATION, not a
 * discovery.</b> The unreachability argument was already derived verbatim
 * in this javadoc before bead inviscid-0nx.28 touched it, where it was
 * framed as a VIRTUE ("not vacuous by accident, but exact by
 * construction"). The mathematics was right and documented, and it had
 * passed review in that framing. What bead inviscid-0nx.28 changed is the
 * verdict, not the derivation: a guard that provably cannot fire on the
 * geometry it is deployed on certifies nothing there, and the question
 * review should have asked -- "then what protects the campaign?" -- had no
 * answer. What is genuinely new is the even-extent enumeration test and
 * limb W2. (The bead's own stated premise, that a uniform signed
 * background makes W1 throw at tick 0 on {@code 8^3}, is separately FALSE
 * and was never satisfiable; the plan audit that recorded it as "confirmed
 * in source" was wrong. Verified by exhaustive execution -- see the
 * enumeration test.)
 *
 * <p>W1's structural no-op was harmless while the estimator's input was a
 * localized non-negative packet. It is a liability once the input is a
 * difference field whose localization is empirical rather than provable,
 * because W1 alone then certifies nothing whatsoever about a campaign run.
 * Hence W2.
 *
 * <h3>Limb W2 -- saturation ({@link #assertResponseLocalized};
 * matched-pair path only)</h3>
 * The failure mode a centered-origin run actually faces is one W1 is
 * structurally blind to. The torus bounds the achievable displacement, so
 * {@code M_d(t)} cannot grow without limit: it relaxes onto the second
 * moment of the UNIFORM distribution about the same origin, after which
 * the OLS slope UNDERSTATES transport, worse the longer the run --
 * silently, with no exception, no {@code NaN}, and no empty {@link
 * OptionalDouble}.
 *
 * <p>W2 measures that directly. Define, per direction, the <b>moment
 * saturation ratio</b>
 * <pre>
 *   sat_d(t) = M_d(t) / M_d^unif,
 *   M_d^unif = (1/N) * sum_cell proj_d(cell - origin)^2
 * </pre>
 * ({@code M_d^unif} is computed by running the SAME {@code secondMoment}
 * code over an all-ones field, so the two quantities cannot drift apart).
 * W2 asserts {@code max_d sat_d(t) <=} {@link
 * #RESPONSE_MOMENT_SATURATION_TOLERANCE} at every snapshot. A response
 * still concentrated near the origin has {@code sat ~ 0}; a response that
 * has filled the box has {@code sat ~ 1}.
 *
 * <p><b>What this replaces, and why.</b> The first version of W2 measured
 * the fraction of {@code ||D||_1} sitting on the half-period SHELL (some
 * axis at {@code |coord-origin| >= extentAxis/2}), against a 1% tolerance
 * justified by the shell's 33.0% cell fraction under complete
 * delocalization ({@code (L^3-(L-1)^3)/L^3} = {@code 169/512} at
 * {@code L=8} -- that arithmetic was correct). The measured SET was not.
 * With {@code L=8} and the origin at {@code 4}, {@code |i-4| >= 4} holds
 * only at {@code i=0}, so the "shell" is a THREE-PLANE SKIN, and a
 * difference field spread uniformly over the interior {@code (L-1)^3} with
 * zero on those three planes measures <b>exactly 0.0</b> -- W2 stayed
 * SILENT on a maximally delocalized field. That is the same species of
 * defect as W1's: a precondition that cannot fire on the case it exists
 * for. Under the saturation ratio the same field measures
 * {@code 4.0/5.5 = 0.727} and fires; it is pinned as a regression by
 * {@code
 * AnisotropyProbeTest#saturationLimbW2FiresOnTheInteriorUniformFieldTheShellMeasureMissed}.
 * The shell measure was also only ever a LAGGING proxy -- it reports
 * nothing until mass reaches the very edge -- whereas the saturation ratio
 * tracks the biased quantity itself.
 *
 * <p><b>The tolerance is calibrated against SLOPE BIAS, which is what
 * actually matters -- but the calibration is FAMILY-CONDITIONAL, and that
 * qualifier is load-bearing.</b> The predecessor's 1% was calibrated
 * against distinguishability from full delocalization, not against the
 * estimator's error. Calibrating against the error instead: place a
 * wrapped (periodic-image-summed) Gaussian of variance {@code s} on
 * {@code L} sites about a centered origin -- the exact 1-D problem the
 * {@code X100} projection solves -- and compare {@code dM/ds} against the
 * free-diffusion value of 1. The fractional understatement of the
 * instantaneous rate at {@code sat = 0.25} is <b>0.53% at {@code L=6},
 * 0.81% at {@code L=8}, 1.36% at {@code L=40}</b>: nearly
 * {@code L}-independent because both sides scale with the box. {@link
 * #RESPONSE_MOMENT_SATURATION_TOLERANCE} is therefore set to
 * {@code 0.25}, and the calibration is re-derived in-test by {@code
 * AnisotropyProbeTest#momentSaturationToleranceBoundsTheSlopeUnderstatement}
 * rather than asserted in prose. Note the contrast with the old number: 1%
 * of MASS pinned at {@code (L/2)^2 = 16} contributes {@code 0.16} cells^2
 * to {@code M_d}, which against an early-time {@code M_d ~ 1} is a ~16%
 * distortion -- the old tolerance admitted an order of magnitude MORE bias
 * than its "1%" suggested.
 *
 * <h3>"~1.4% slope understatement" is a statement ABOUT THAT FAMILY, not
 * about an arbitrary response, and must not be quoted as one</h3>
 * An earlier version of this javadoc read the {@code 0.25} tolerance as
 * meaning, flatly, "at most ~1.4% instantaneous slope understatement at
 * the worst snapshot". That is an overclaim in two independent ways.
 * <ul>
 * <li><b>{@code sat} is a LEVEL; the estimator's error is driven by a
 * RATE.</b> The wrapped Gaussian is a ONE-PARAMETER family -- level and
 * rate are both functions of the single variance {@code s}, hence tied to
 * each other, so bounding one bounds the other. Nothing ties them for a
 * general two-component response: a packet core plus a decorrelation halo
 * of L1 mass fraction {@code f}. Writing
 * {@code M = (1-f)s + f*M^unif} gives
 * {@code dM/dt = (1-f)ds/dt + (df/dt)(M^unif - s)}. {@code sat} bounds
 * {@code f}; it says nothing whatever about {@code df/dt}. At {@code L=8},
 * {@code M^unif - s} is of order {@code 4.5} while {@code ds/dt} is set by
 * a {@code ~1.1e-4}/tick contact rate ({@link QuantaExchangeRule}'s own
 * javadoc), so the halo term is not a small correction -- and it carries
 * the OPPOSITE sign from the failure mode described above. Solving
 * {@code (1-f)*1 + f*5.5 = 0.25*5.5} puts {@code f = 0.083} exactly at
 * the tolerance boundary: a third of the measured second moment could be
 * halo rather than transport while W2 reads "compliant". The halo is
 * structural, not hypothetical -- {@link Necronomata}'s
 * {@code deltaA = QUANTUM_RATE * frequency} turns a quanta difference into
 * an ANGLE difference and {@code ContactScan}'s predicate is angle-driven,
 * so divergence propagates through the contact graph independently of
 * quanta transport. That is what limb W3 below is for, and why
 * {@code f_halo} is ASSERTED rather than merely reported.</li>
 * <li><b>The calibration validates the {@code X100} limb only.</b> Its
 * family is 1-D and its uniform reference is literally
 * {@code uniformSecondMoment(..., X100)}. It says nothing about a response
 * whose spread is not aligned with the probe-direction set -- see the
 * blind spot below.</li>
 * </ul>
 *
 * <h3>The MODEL-FREE reading of {@code sat}, which holds for any
 * response</h3>
 * {@code M_d} is the second moment of the probability measure
 * {@code mu_t = |D|/||D||_1}, so Markov's inequality applies verbatim: for
 * any threshold {@code P},
 * <pre>
 *   mu_t{ cell : proj_d(cell-origin)^2 &gt;= P }  &lt;=  M_d(t) / P
 * </pre>
 * Taking {@code P = M_d^unif} states it in the tolerance's own units:
 * <b>at most {@code sat_d} of the response's L1 MASS can sit at or beyond
 * uniform-typical squared displacement along {@code d}</b> -- so at most a
 * quarter of it at the {@code 0.25} tolerance. Taking {@code P} at W1's
 * ceiling gives a sharper one: W1 admits only {@code |d_a| <= L/2} per
 * axis, so {@code proj_X111^2 <= (3*(L/2))^2/3 = 3*(L/2)^2}, which is
 * {@code 48} at {@code L=8} against {@code M_X111^unif = 6.0}; at
 * {@code sat = 0.25} that is {@code M_X111 <= 1.5}, hence <b>at most
 * {@code 3.125%} of the mass at the extreme wrap-safe displacement</b>.
 * Both are distribution-free -- no family, no fit -- and are what
 * {@code sat} bounds unconditionally. Pinned by {@code
 * AnisotropyProbeTest#momentSaturationBoundsDelocalizedMassDistributionFree}.
 *
 * <h3>W2's blind spot: a response confined to the LEAST-PROBED axis
 * (round-2 review; MEASURED, UNBOUNDED, and OPEN)</h3>
 * The guard is weakest precisely in the anisotropic regime this instrument
 * exists to detect. {@link StructureFactor.Direction}'s set is
 * {@code {dx, (dx+dy)/sqrt(2), (dx+dy+dz)/sqrt(3)}} -- no probe isolates
 * {@code y} or {@code z}. A {@code z}-only displacement therefore reaches
 * {@code X111} alone, divided by 3, while {@code M^unif} is the 3-D
 * uniform moment with all three axes contributing. Measured on
 * {@code 8^3} about {@code (4,4,4)}, where {@code M^unif} is
 * {@code (5.5, 5.75, 6.0)} for {@code (X100, X110, X111)}:
 * <ul>
 * <li>a response uniform along the origin's {@code z} column reaches
 * {@code M_X111 = 5.5/3 = 1.833}, i.e. {@code sat = 0.3056} -- which is
 * the ENTIRE attainable range for a {@code z}-confined response. The
 * {@code 0.25} tolerance therefore sits at {@code 81.8%} of that range,
 * not at 25% of it;</li>
 * <li>compared AT MATCHED SATURATION, which is the only fair comparison:
 * a separable wrapped Gaussian with per-axis variance rates
 * {@code (1,1,20)} understates the true <b>INSTANTANEOUS</b> {@code X111}
 * rate by <b>19.6164% at {@code sat = 0.20} and 34.7264% at
 * {@code sat = 0.24}</b> -- W2 SILENT at both, since both are under the
 * {@code 0.25} tolerance. The isotropic {@code (1,1,1)} field at exactly
 * those saturations understates by {@code 0.1314%} and {@code 0.5411%}.
 * The calibration above is honoured for the isotropic response and missed
 * by a factor of {@code 150} and {@code 64} respectively for the
 * anisotropic one. <b>"Instantaneous" is load-bearing and must be carried
 * whenever these two numbers are quoted:</b> the corresponding
 * WINDOW-AVERAGED understatement -- an OLS slope fitted across a window
 * ending at that saturation, rather than {@code dM/ds} at it -- is only
 * ~{@code 7.6%}, a 2-3x smaller figure for the same field. W2 is a
 * PER-SNAPSHOT trip, so the instantaneous convention is the correct one
 * here; omitting the word is what produced a cross-harness disagreement
 * of exactly that 2-3x factor in round 3 (see below);</li>
 * <li>statically: a field uniform over the interior {@code z} line
 * ({@code k=1..7}) measures {@code sat = 0.2222} and PASSES; adding a
 * compact isotropic blob at the origin worth 2%/5%/10% of its mass LOWERS
 * that to {@code 0.2179}/{@code 0.2116}/{@code 0.2020}, still silent. The
 * {@code y}-pencil equivalent measures {@code 0.3478} and FIRES -- the
 * blindness is axis-specific, and {@code z} is the worst axis by exactly
 * the {@code 3/2} the two normalizations differ by.</li>
 * </ul>
 * Measured, not argued, by {@code
 * AnisotropyProbeTest#saturationLimbW2IsBlindToAResponseConfinedToTheLeastProbedAxis}.
 *
 * <p><b>The gap is UNBOUNDED, not the 20-35% effect those rows suggest --
 * and this is the reason a small measured {@code sat} certifies nothing.</b>
 * {@code momentSaturation} is a MASS-WEIGHTED MEAN, so an escaping mass
 * fraction {@code f} can be hidden arbitrarily well by parking the
 * remaining {@code (1-f)} on the origin: the {@code X111}-rate
 * understatement stays at {@code 99.81%} while {@code sat} scales LINEARLY
 * down with {@code f}. At {@code f = 0.001} such a field measures
 * {@code sat = 0.00089} -- <b>5.7x BELOW Phase A's own worst-case census
 * value of {@code 0.00512}</b>. The "49x headroom" recorded below therefore
 * bounds how delocalized the observed regime's mass is; it does NOT certify
 * that the regime's transport estimate carries low bias, and must never be
 * quoted as if it did.
 *
 * <p><b>Two ORTHOGONAL defect classes, and the named repair closes only
 * one.</b> The direction-set class (no probe isolates {@code y} or
 * {@code z}) is repaired by a direction set, or a per-axis reference, that
 * probes them -- a change to {@link StructureFactor}'s shared
 * cross-direction-comparable convention, which BOTH estimators consume and
 * on which E.2's composition depends. The ESCAPED-FRACTION class is not: the
 * witness field above reads {@code sat = 0.0029} on an {@code X001} probe
 * and stays silent, so adding {@code y}/{@code z} directions leaves it
 * undetected. It is not a {@link StructureFactor.Direction} change at all --
 * it needs a measure other than a mass-weighted mean (a tail/quantile
 * statistic, or an escaped-mass bound). <b>This is NOT fixed here.</b> Both
 * classes are tracked as bead {@code inviscid-fii}, resized to cover both,
 * to land BEFORE E.2 composes on W2. What makes it non-blocking TODAY is not
 * the census margin but that the matched-pair path has no campaign consumer
 * yet, and the only campaign that has run this code is Phase A, whose
 * difference field IS the packet field ({@code f_halo == 0} at all 1024
 * snapshots) rather than an adversarially-shaped response.
 *
 * <p><b>Round-3 cross-harness dispute, and how it resolved</b> -- recorded
 * because the next reader may otherwise re-derive the withdrawn numbers.
 * A second harness initially reported {@code -1.9%} (opposite sign) for the
 * isotropic case and {@code 9.7%}/{@code 11.3%} for the anisotropic one.
 * <b>The figures above stand; the competing ones were withdrawn.</b> Root
 * cause was a truth-model convention error on the other side: the isotropic
 * sign flip came from comparing an OLS slope over a window against the
 * ANALYTIC CONTINUUM rate of a sub-cell Gaussian -- lattice discretization
 * suppresses measured variance at small {@code s} and does so
 * non-uniformly across the window, corrupting the slope by 1.2-2.1%, which
 * was the entire {@code -1.9%}. The central-difference construction used
 * here does not have that transient, and a convention-free {@code 96^3}
 * big-box comparand agrees with it ({@code +0.006%} to {@code +0.043%}); three
 * independent constructions now put the isotropic bias at ~{@code 0}. The
 * 2-3x gap on the anisotropic figures was window-averaged vs instantaneous,
 * per the bullet above.
 *
 * <p><b>What W2 still does not detect, beyond that.</b> It bounds the bias
 * from the response filling the box. It says nothing about drift of the
 * response's CENTROID away from the fixed origin, which raises
 * {@code M_d} for a reason other than spread -- outside the periodic-wrap
 * question this precondition is about. Decorrelation of the two halves,
 * previously listed here as undetected with {@code ||D||_1} growth
 * "reported, not asserted", is now limb W3. The measured maximum is
 * RETURNED ({@link MatchedPairTransport#maxMomentSaturation()}), not
 * merely used as a trip-wire -- a guard whose only output is "did not
 * throw" cannot be told apart from a vacuous one.
 *
 * <h3>Limb W3 -- halo fraction ({@link #assertResponseCoreDominates};
 * matched-pair path only)</h3>
 * W2 asks where the response's mass SITS. W3 asks what that mass IS.
 * Define the <b>decorrelation halo fraction</b>
 * <pre>
 *   f_halo(t) = 1 - |S| / ||D(.,t)||_1
 * </pre>
 * Both terms are already computed on this path, and the quantity is exact
 * and model-free -- a definition, not a fit. By the
 * {@code |S| &lt;= ||D||_1} bound above it lies in {@code [0,1)}. It is
 * {@code 0} exactly when every nonzero cell of {@code D} shares a sign,
 * i.e. when the response is pure transported excess; it rises strictly
 * above zero as the two halves' firing sequences diverge and {@code D}
 * acquires cancelling {@code +/-} pairs, each of which adds to
 * {@code ||D||_1} while leaving {@code S} untouched. W3 asserts
 * {@code f_halo(t) <=} {@link #RESPONSE_HALO_FRACTION_TOLERANCE} at every
 * snapshot; the maximum is returned as {@link
 * MatchedPairTransport#maxHaloFraction()}.
 *
 * <p><b>What a LEVEL bound on {@code f_halo} does and does not buy --
 * stated precisely, because this limb exists to repair exactly that
 * confusion in W2's own justification.</b> It does NOT bound
 * {@code df/dt}, and so does not by itself bound the halo's contribution
 * to the OLS slope. What makes it bite anyway is the ANCHOR: precondition
 * (3) of {@link #transportEstimateMatchedPair} asserts
 * {@code ||D(.,0)||_1 == |S|}, i.e. {@code f_halo(0) = 0} -- to within
 * {@link #INJECTED_EXCESS_TOLERANCE} ({@code 1e-6}), which is what that
 * precondition actually enforces, NOT bit-exactly. Bit-exactness holds only
 * through the runner's support check, which independently establishes that
 * {@code D(.,0)} is single-signed.
 * Decomposing the measure into core and halo parts with conditional
 * moments {@code m_core}, {@code m_halo},
 * <pre>
 *   M_d(t) - M_d(0) = [m_core(t) - m_core(0)]
 *                     + f_halo(t) * [m_halo(t) - m_core(t)]
 * </pre>
 * the second term is the WHOLE of the halo's contribution to the measured
 * change, and its magnitude is at most
 * {@code f_halo(t) * max_cell proj_d^2} -- itself bounded by W1. So the
 * tolerance bounds the halo's total moment EXCURSION over the window,
 * which a level bound alone could not do and which the
 * {@code f_halo(0) = 0} anchor supplies for free. Pinned by {@code
 * AnisotropyProbeTest#haloFractionIsAnchoredAtZeroAndBoundsTheHaloMomentExcursion}.
 *
 * <p><b>What that excursion bound is WORTH, numerically -- because the
 * algebra is correct but its magnitude is the point, and quoting the
 * algebra without the magnitude over-advertises this limb.</b> At
 * {@code L=8} under W1's ceiling ({@code |d_a| <= L/2}): for {@code X111},
 * {@code f_halo * max proj^2 = 0.05 * 48 = 2.4}, against the
 * {@code M_X111 <= 1.5} that W2 already admits at its own tolerance -- the
 * excursion bound EXCEEDS the entire admissible range and so constrains
 * nothing W2 does not. For {@code X100} it is {@code 0.05 * 16 = 0.8}
 * against {@code 1.375}, i.e. it still permits the halo to be ~{@code 58%}
 * of the measured moment change. That is structural, not an accident of
 * these numbers: it is the SAME quantity as the share bound {@link
 * #RESPONSE_HALO_FRACTION_TOLERANCE}'s own javadoc already concedes
 * vacuous, taken over a smaller denominator, hence a fortiori weaker.
 * <b>W3's defence is therefore the MEASURED LICENSED WINDOW recorded on
 * that constant</b> -- which genuinely refuses {@code 8^3} x 128 ticks on a
 * signed background -- and not this excursion argument, which should be
 * read as an existence statement rather than as a useful bound.
 *
 * <p>W3 reads EXACTLY {@code 0} on the {@code {0, packetQuanta}} path, at
 * every one of the 1024 Phase A campaign snapshots (measured): the control
 * run holds no quanta, so {@code D == f_P >= 0} and
 * {@code ||D||_1 == |S|} identically. Like W2 it is applied only to the
 * matched-pair path, for the same bit-for-bit-preservation reason.
 *
 * <p>W2 is deliberately NOT applied to the {@code {0, packetQuanta}}
 * path. That path is pinned bit-for-bit ({@code SeamGoldenCompatTest},
 * {@code SubstrateFactorySeamTest}); adding a new failure mode to it
 * would be a behaviour change, not a re-derivation. The campaign path
 * ({@link #runCampaign} -&gt; {@link #runOneSeed}) therefore still carries
 * only the proven-vacuous W1, which is a real and acknowledged gap. It is
 * bounded by MEASUREMENT rather than left open -- but by two measurements
 * of very different strength, which must not be conflated:
 * <ul>
 * <li><b>the campaign number, measured out-of-suite:</b> the worst
 * saturation over all 1024 Phase A snapshots (8 seeds x 128 ticks) is
 * {@code 0.00512}, a 49x margin under the {@code 0.25} tolerance.
 * Reproducible from this tree via {@link #main} {@code --census};</li>
 * <li><b>the number actually pinned in the suite:</b> a 2-seed, 24-tick
 * reduction of that census. Saturation grows monotonically with ticks, so
 * what the suite asserts is a LOWER BOUND on the 128-tick worst case, not
 * the worst case -- a tripwire on the regime, with the coverage and the
 * low bias quantified under "the {@code {0, packetQuanta}} path" above.
 * An earlier version of this paragraph said the campaign's own worst-case
 * saturation was "pinned", which the suite does not do and which
 * contradicted that test's own (correct) javadoc.</li>
 * </ul>
 * Both read the saturation off the packet field directly, which is exact
 * here: on the zero-background substrate the control holds no quanta, so
 * the difference field IS the packet field.
 *
 * <p>This class still chooses SHORT RUNS (a small, budget-bounded tick
 * count) over a periodic-aware circular-moment technique as the overall
 * strategy: (a) a circular (von Mises-style) moment does not compose
 * cleanly across the three simultaneously-probed directions with their
 * different {@code sqrt(Nd)} normalizations; (b) {@link
 * QuantaExchangeRule} is a single-quantum, sparse-contact
 * (~1.1e-4/tick, {@code QuantaExchangeRule}'s own javadoc) process, so
 * spread accumulates slowly and a short-tick regime is physically
 * appropriate for early-time transport characterization regardless of
 * the exact-correctness guard's reach.
 *
 * <h2>SPECTRAL estimator</h2>
 * Per direction {@code d}: {@code points = structureFactor.spectrum
 * (fieldByTick, d)} (public, real-field overload); {@code ridge =
 * structureFactor.extractRidge(points)} over the RAW, UNFILTERED points
 * list -- per bead inviscid-0nx.9's stacked-review final state (comments
 * on this bead), {@code Ridge.slope()} is now genuinely
 * cross-direction-comparable ({@code sqrt(Nd)}-scaled) and real-field
 * mirror pairs REINFORCE rather than cancel in {@code extractRidge}'s
 * unweighted OLS -- no manual pre-filtering or normalization is applied
 * here, matching the reviewed-safe consumption pattern. {@code
 * magnitude_d = Math.abs(ridge.slope())} (a speed, not a signed
 * propagation direction, so the max/min ratio answers "how much faster",
 * not "which way"). The {@code [111]} probe's narrower/fewer-point range
 * (real FCC physics, not a bug -- see {@link StructureFactor}'s class
 * javadoc) is surfaced directly: {@link DirectionMagnitude#sampleSize()}
 * carries {@code points.size()} so a report reader sees the asymmetric
 * precision, not just prose about it.
 *
 * <p><b>On the campaign's fully-degenerate spectral result (stacked
 * review, FIX 3): this is the EXPECTED signature of purely diffusive
 * dynamics, not an instrument malfunction.</b> Diffusion has no
 * propagating branch (the archetypal diffusive dispersion relation is
 * {@code omega ~ i*D*k^2}, overdamped -- no real, nonzero temporal
 * frequency dominates any {@code k} the way a propagating wave's
 * {@code omega = c*k} would), so a temporal-FFT ridge-slope estimator
 * legitimately finding no dominant nonzero-frequency peak at any
 * {@code k}, in any direction, is that estimator correctly reporting
 * "no propagating collective mode here" -- consistent with, not
 * contradicting, the TRANSPORT estimator's real-space diffusive signal.
 * The two estimators measure DIFFERENT PHYSICS (propagating-mode speed
 * vs. real-space spread rate); a diffusive system with no propagating
 * branch produces exactly this pattern by construction, not
 * "disagreement" in the pejorative sense.
 *
 * <h2>Non-fabrication contract</h2>
 * {@link #ratio(Map)} (and {@link #pooledEstimate}'s internal
 * ratio-of-means) is the choke point every ratio computation in this
 * class reduces through: if the smaller of the direction magnitudes (or
 * pooled means) is at or below {@link #RATIO_DEGENERATE_EPSILON}, the
 * ratio is {@link OptionalDouble#empty()} -- covers the K=0
 * (collision-free) baseline exactly (a field that never changes produces
 * an OLS slope of EXACTLY {@code 0.0} in every direction, by
 * construction), without a separate "is this K=0" special case.
 *
 * @author halhildebrand
 */
public final class AnisotropyProbe {

    /**
     * One direction's measured magnitude: {@code magnitude} is
     * {@code |D_hat(d)|} (transport) or {@code |ridge.slope()|}
     * (spectral); {@code sampleSize} is the point/tick count the
     * magnitude was fit from (informational -- surfaces the [111]
     * fewer-points caveat for the spectral estimator).
     */
    public record DirectionMagnitude(StructureFactor.Direction direction,
                                      double magnitude, int sampleSize) {
    }

    /**
     * One estimator's full per-direction result plus the NAIVE per-seed
     * anisotropy ratio -- {@link OptionalDouble#empty()} means
     * DEGENERATE (see class javadoc, "Non-fabrication contract"), never
     * a fabricated number. See class javadoc, "STACKED-REVIEW
     * CORRECTION": this per-seed ratio is a diagnostic, not the
     * significance statistic -- use {@link #pooledEstimate} for that.
     */
    public record EstimatorResult(Map<StructureFactor.Direction, DirectionMagnitude> perDirection,
                                   OptionalDouble ratio) {
        public EstimatorResult {
            perDirection = Map.copyOf(perDirection);
        }
    }

    /**
     * Both estimators' results for ONE seed, computed from the SAME
     * {@code fieldByTick} snapshot sequence (bead's "same configuration"
     * requirement) -- disagreement between {@link #transport()} and
     * {@link #spectral()} is a finding, read directly off this record,
     * never averaged away. {@code totalCollisions}/{@code
     * effectiveCollisions} (FIX 2, stacked review) are this seed's {@link
     * CollisionStatistics} counts over the full run -- surfaced so a
     * reader can judge whether the OLS window sampled enough real
     * transfer events to be in a genuinely diffusive regime, rather than
     * a handful-of-discrete-hops small-N regime.
     */
    public record SeedResult(long seed, EstimatorResult transport,
                              EstimatorResult spectral, long totalCollisions,
                              long effectiveCollisions) {
    }

    /**
     * The matched-pair TRANSPORT result (bead inviscid-0nx.28; T2 {@code
     * design-seeding-radius.md} §D-A): the estimator's own output on the
     * DIFFERENCE field, plus the three quantities the re-derivation
     * showed a reader needs in order to interpret it.
     *
     * @param transport             the estimator result, computed on
     *                              {@code D = f_P - f_C}
     * @param injectedExcess        {@code S = sum_cell D}, the exactly
     *                              tick-invariant injected packet total
     *                              (invariant (I1)); {@code ||D||_1} is
     *                              bounded below by {@code |S|} at every
     *                              tick
     * @param responseL1First       {@code ||D||_1} at the first snapshot.
     *                              ASSERTED, not merely reported, to equal
     *                              {@code |S|} -- see {@link
     *                              #transportEstimateMatchedPair}, which
     *                              refuses the pair otherwise. This is the
     *                              SHARP same-background check; the
     *                              runner's {@code S == 30*packetQuanta}
     *                              equality is the blunt one
     * @param responseL1Last        {@code ||D||_1} at the last snapshot.
     *                              Read WITH the ratio: a degenerate
     *                              (empty) ratio at
     *                              {@code responseL1Last == |S|} means the
     *                              response never left the seed cell,
     *                              whereas the same empty ratio at
     *                              {@code responseL1Last >> |S|} means it
     *                              spread with no net second-moment trend
     *                              -- see the class javadoc, "Degeneracy
     *                              of the ratio-of-means", and {@code
     *                              AnisotropyProbeTest#responseL1LastSeparatesTheTwoDegenerateRatioDiagnoses}
     *                              which exercises exactly that
     *                              disambiguation
     * @param maxMomentSaturation   the largest over-ticks, over-directions
     *                              value of limb W2's measured moment
     *                              saturation ratio
     *                              {@code M_d(t)/M_d^unif}. Returned, not
     *                              merely asserted, so the guard has an
     *                              observable output rather than only "did
     *                              not throw"
     * @param maxHaloFraction       the largest over-ticks value of limb
     *                              W3's measured DECORRELATION HALO
     *                              FRACTION {@code 1 - |S|/||D(.,t)||_1}
     *                              -- the share of the L1 mass the second
     *                              moment is computed over that is signed
     *                              structure created by the two halves
     *                              decorrelating rather than transported
     *                              injected excess. Exactly {@code 0} on
     *                              the {@code {0, packetQuanta}} path.
     *                              Reported for the same reason as {@code
     *                              maxMomentSaturation}, and because it is
     *                              the quantity that tells a flat-slope
     *                              result "the packet did not move" apart
     *                              from "the pair stopped being a
     *                              measurement of the packet"
     */
    public record MatchedPairTransport(EstimatorResult transport,
                                        double injectedExcess,
                                        double responseL1First,
                                        double responseL1Last,
                                        double maxMomentSaturation,
                                        double maxHaloFraction) {
    }

    /**
     * A percentile bootstrap confidence interval over the per-seed ratios
     * that were present (non-degenerate). Diagnostic only -- see class
     * javadoc, "STACKED-REVIEW CORRECTION": this statistic is bounded
     * BELOW by exactly 1.0 (every individual sample is itself
     * {@code max_d/min_d >= 1.0} by construction), is upward-biased by
     * seed-to-seed noise (an order-statistic artifact, T3 {@code
     * critique-pattern-max-min-ratio-order-statistic-bias}), and "CI
     * excludes 1.0" on THIS statistic is not evidence of a real
     * direction-linked effect. Use {@link #pooledEstimate}'s {@code
     * pooledRatio} + permutation p-value for that. {@code
     * nSeedsDegenerate} counts seeds whose ratio was {@link
     * OptionalDouble#empty()} -- excluded from the resample, but
     * reported, not silently dropped.
     */
    public record BootstrapCi(double mean, double lower, double upper,
                               int nSeedsUsed, int nSeedsDegenerate) {
    }

    /**
     * One direction's seed-pooled statistic: the mean magnitude across
     * all seeds, with a resample-then-aggregate bootstrap CI (seed
     * indices resampled with replacement, per-direction mean recomputed
     * from the resample -- see {@link #pooledEstimate}).
     */
    public record PooledDirectionStats(StructureFactor.Direction direction,
                                        double mean, double ciLower,
                                        double ciUpper, int nSeeds) {
    }

    /**
     * THE significance statistic for a Phase A campaign (see class
     * javadoc, "STACKED-REVIEW CORRECTION"): the seed-pooled
     * ratio-of-means ({@code pooledRatio}, with its own
     * resample-then-aggregate bootstrap CI {@code pooledRatioCiLower}/
     * {@code pooledRatioCiUpper}), plus a permutation/null-calibration
     * test ({@code permutationPValue} = (countGe+1)/(permutationCount+1),
     * the standard +1 continuity-corrected fraction of direction-label-
     * shuffled resamples whose pooled ratio meets or exceeds the
     * observed one -- the observed statistic is exchangeable with the
     * null draws under H0, so it counts as one of its own reference set,
     * and a finite permutation count can never license a fabricated
     * exact-zero p-value; {@code permutationNull95} = the null
     * distribution's 95th percentile, for context; {@code
     * permutationCount} = how many permutation draws were actually
     * usable, i.e. non-degenerate).
     */
    public record PooledResult(Map<StructureFactor.Direction, PooledDirectionStats> perDirection,
                                OptionalDouble pooledRatio,
                                double pooledRatioCiLower,
                                double pooledRatioCiUpper,
                                double permutationPValue,
                                double permutationNull95,
                                int permutationCount) {
        public PooledResult {
            perDirection = Map.copyOf(perDirection);
        }
    }

    /** The full Phase A campaign result: every seed's raw result plus both estimators' pooled/naive statistics. */
    public record Report(List<SeedResult> perSeed, BootstrapCi transportCi,
                          BootstrapCi spectralCi, PooledResult pooledTransport,
                          PooledResult pooledSpectral, Point3i extent,
                          int ticks, Point3i originCell, long[] seeds,
                          int packetQuanta) {
        public Report {
            perSeed = List.copyOf(perSeed);
            seeds = seeds.clone();
        }

        @Override
        public long[] seeds() {
            return seeds.clone();
        }
    }

    /**
     * Below this, a direction's magnitude (or pooled mean) is treated as
     * zero for ratio purposes -- see class javadoc, "Non-fabrication
     * contract". Chosen far below any genuine transport/group-velocity
     * scale in this domain, and the K=0 baseline produces an EXACT
     * {@code 0.0}, not merely a small one, so the precise epsilon value
     * is not load-bearing for that case.
     */
    static final double RATIO_DEGENERATE_EPSILON = 1e-9;

    /**
     * Wrap-safety limb W2's tolerance (bead inviscid-0nx.28, fix round):
     * the largest MOMENT SATURATION RATIO {@code max_d M_d(t)/M_d^unif}
     * any snapshot of the difference field may reach before {@link
     * #assertResponseLocalized} refuses it. See the class javadoc, "Limb
     * W2 -- saturation", for the full derivation and for why the
     * predecessor's shell-mass-fraction measure was retired (it read
     * exactly {@code 0.0} on a maximally delocalized interior-uniform
     * field).
     *
     * <p>Calibrated against SLOPE BIAS, not against distinguishability
     * from full delocalization: for a wrapped Gaussian on {@code L} sites
     * about a centered origin, {@code sat = 0.25} corresponds to an
     * instantaneous-rate understatement of 0.53% at {@code L=6}, 0.81% at
     * {@code L=8} and 1.36% at {@code L=40} -- nearly {@code L}-independent
     * because both sides scale with the box. Re-derived in-test by {@code
     * AnisotropyProbeTest#momentSaturationToleranceBoundsTheSlopeUnderstatement}.
     */
    static final double RESPONSE_MOMENT_SATURATION_TOLERANCE = 0.25;

    /**
     * Wrap-safety limb W3's tolerance (bead inviscid-0nx.28, round-2 fix):
     * the largest DECORRELATION HALO FRACTION
     * {@code f_halo(t) = 1 - |S|/||D(.,t)||_1} any snapshot of the
     * difference field may reach before {@link
     * #assertResponseCoreDominates} refuses it. See the class javadoc,
     * "Limb W3 -- halo fraction", for what it does and does not buy.
     *
     * <p><b>This is a DECLARED REGIME BOUND with measured headroom, not a
     * derived constant</b> -- said plainly, because W2's neighbouring
     * tolerance IS derived and the two must not be read the same way.
     * There is no derivation available: a small {@code f_halo} does not
     * imply a small halo contribution to {@code M_d} (the halo's own
     * conditional moment is unbounded below the W1 ceiling, so the
     * rigorous share bound {@code f_halo * max proj^2 / M_d} goes vacuous
     * whenever the core is tightly localized). {@code f_halo} is asserted
     * because it is the only MODEL-FREE handle on the composition of the
     * mass the moment is computed over, not because a bound on it bounds
     * the bias.
     *
     * <p>The value is set so that it bites well before the halo can carry
     * the share of the second moment that {@link
     * #RESPONSE_MOMENT_SATURATION_TOLERANCE} alone would admit: a halo
     * fraction of {@code 0.083} makes a third of {@code M_d} halo rather
     * than transport at {@code sat = 0.25} (the round-2 critique that
     * motivated this limb), so the tolerance must sit below that, and
     * {@code 0.05} does.
     *
     * <p><b>The licensed window, measured this round rather than assumed.</b>
     * On the real signed-quanta substrate ({@code packetQuanta = 30},
     * seeds 42/43/44), worst-over-ticks and worst-over-seeds
     * {@code f_halo}:
     * <pre>
     *   4^3,  32 ticks, amp  60 : 0.0153   (the in-suite fixture)  PASSES
     *   4^3,  64 ticks, amp  60 : 0.0239                           PASSES
     *   6^3,  64 ticks, amp  60 : 0.0239                           PASSES
     *   8^3,  64 ticks, amp  60 : 0.0175                           PASSES
     *   4^3, 128 ticks, amp  60 : 0.0426                           PASSES
     *   8^3, 128 ticks, amp  60 : 0.0683                           REFUSED
     *   8^3, 128 ticks, amp 200 : 0.0586                           REFUSED
     * </pre>
     * <b>The last two rows are the point, not an oversight.</b> This limb
     * genuinely narrows the licensed window: a 128-tick matched pair on a
     * signed {@code 8^3} background is NOT licensed at this tolerance, and
     * that was previously invisible -- those runs' W2 saturations are
     * {@code 0.014}-{@code 0.019}, two orders under W2's tolerance, so W2
     * passes them without comment while two thirds of a tenth of the L1
     * mass the moment is normalized by is decorrelation structure rather
     * than packet. A campaign that needs 128 ticks on a signed background
     * must shorten the window, or raise this tolerance <b>on evidence that
     * the halo does not bias its estimate</b> -- not because the number was
     * inconvenient.
     *
     * <p>On the Phase A {@code {0, packetQuanta}} substrate {@code f_halo}
     * is EXACTLY {@code 0} at every one of the 1024 campaign snapshots
     * (measured): the control run holds no quanta, so {@code D == f_P >= 0}
     * and {@code ||D||_1 == |S|} identically.
     */
    static final double RESPONSE_HALO_FRACTION_TOLERANCE = 0.05;

    /**
     * Absolute tolerance on the tick-to-tick drift of the difference
     * field's signed total (the injected excess {@code S}, invariant
     * (I1)). Quanta are integer-valued and every per-tick sum is taken in
     * the same order over {@code double}s, so the comparison is exact in
     * practice for any campaign that fits {@code 2^53}; the tolerance
     * exists so the check reports a physics violation rather than a
     * last-bit artifact.
     */
    static final double INJECTED_EXCESS_TOLERANCE = 1e-6;

    static final int  BOOTSTRAP_RESAMPLES = 5000;
    static final long BOOTSTRAP_RNG_SEED  = 1_000_003L;

    /** >= 1000 required by the stacked-review fix; 2000 for extra headroom at negligible cost. */
    static final int  PERMUTATION_COUNT    = 2000;
    static final long PERMUTATION_RNG_SEED = 7_000_001L;

    /**
     * Mean effective-collision count per seed below which the campaign
     * header flags itself small-N/early-time (FIX 2, stacked review) --
     * a documented, not arbitrary-and-hidden, threshold: below this, an
     * OLS fit over the recorded window is fitting mostly zero-change
     * ticks punctuated by a handful of discrete single-quantum jumps,
     * not a settled diffusive trend.
     */
    static final double SMALL_N_EFFECTIVE_COLLISIONS_THRESHOLD = 50.0;

    public static final Point3i DEFAULT_EXTENT       = new Point3i(8, 8, 8);
    public static final int     DEFAULT_TICKS        = 128;
    public static final int     DEFAULT_PACKET_QUANTA = 30;
    /** Literal seed list, per the bead's instruction -- never derived/generated. */
    public static final long[]  DEFAULT_SEEDS        = { 42L, 43L, 44L, 45L,
                                                           46L, 47L, 48L, 49L };

    private static final String GOLDEN_RELATIVE_PATH = "src/test/resources/lga/anisotropy-report-phaseA.tsv";

    private AnisotropyProbe() {
    }

    // ------------------------------------------------------------------
    // Data-agnostic estimator core -- operates on any fieldByTick
    // sequence in StructureFactor.coarseGrainedField's layout, real or
    // synthetic.
    // ------------------------------------------------------------------

    /**
     * The TRANSPORT estimator -- see class javadoc for the exact
     * definition and the origin-relative wrap-safety precondition this
     * method asserts (FIX 6).
     *
     * @param fieldByTick snapshots in {@link
     *                    StructureFactor#coarseGrainedField(com.chiralbehaviors.inviscid.automaton.QuantaField)}'s
     *                    layout; length &gt;= 2 required (an OLS fit
     *                    needs at least two distinct {@code t} values)
     * @param extent      the periodic-wrap extent every snapshot is
     *                    shaped for
     * @param originCell  the FIXED reference cell (the packet's seed
     *                    location) every direction's displacement is
     *                    measured from -- may be any cell, not
     *                    necessarily centered; see class javadoc for the
     *                    exact origin-relative correctness criterion this
     *                    implies
     * @throws IllegalStateException if any snapshot's mass has, on any
     *                                axis, moved past the exact
     *                                half-period distance from {@code
     *                                originCell} -- see class javadoc
     */
    public static EstimatorResult transportEstimate(double[][] fieldByTick,
                                                      Point3i extent,
                                                      Point3i originCell) {
        if (fieldByTick == null || fieldByTick.length < 2) {
            throw new IllegalArgumentException("fieldByTick must have at least 2 snapshots, had "
                                                + (fieldByTick == null ? 0
                                                                        : fieldByTick.length));
        }
        int expected = extent.x * extent.y * extent.z;
        int t = fieldByTick.length;
        for (int i = 0; i < t; i++) {
            if (fieldByTick[i] == null || fieldByTick[i].length != expected) {
                throw new IllegalArgumentException("fieldByTick[" + i
                                                    + "] must have length "
                                                    + expected);
            }
            assertWrapSafe(fieldByTick[i], extent, originCell, i);
        }
        return momentSlopes(fieldByTick, extent, originCell);
    }

    /**
     * The estimator's arithmetic core, shared verbatim by {@link
     * #transportEstimate} and {@link #transportEstimateMatchedPair} so the
     * two paths cannot drift: per direction, the mass-weighted second
     * moment per snapshot, then {@code |OLS slope|} against tick index.
     * Assumes its caller has already validated shapes and asserted the
     * wrap-safety precondition appropriate to ITS field (see the class
     * javadoc, "The periodic-wrap precondition").
     */
    private static EstimatorResult momentSlopes(double[][] fieldByTick,
                                                 Point3i extent,
                                                 Point3i originCell) {
        int t = fieldByTick.length;
        double[] xs = new double[t];
        for (int i = 0; i < t; i++) {
            xs[i] = i;
        }

        Map<StructureFactor.Direction, DirectionMagnitude> perDirection = new EnumMap<>(StructureFactor.Direction.class);
        for (StructureFactor.Direction d : StructureFactor.Direction.values()) {
            double[] moments = new double[t];
            for (int i = 0; i < t; i++) {
                moments[i] = secondMoment(fieldByTick[i], extent, originCell,
                                           d);
            }
            double slope = Math.abs(olsSlope(xs, moments));
            perDirection.put(d, new DirectionMagnitude(d, slope, t));
        }
        return new EstimatorResult(perDirection, ratio(perDirection));
    }

    /**
     * The MATCHED-PAIR transport estimator (bead inviscid-0nx.28; T2
     * {@code design-seeding-radius.md} §D-A) -- the signed-background
     * entry point. Subtracts the control run's field from the packet
     * run's field snapshot by snapshot and runs the SAME estimator core on
     * the difference {@code D = f_P - f_C}. See the class javadoc,
     * "MATCHED-PAIR transport under a signed background", for the full
     * re-derivation this method enforces.
     *
     * <p>Preconditions asserted here, in order:
     * <ol>
     * <li>shapes agree (both sequences the same length &gt;= 2, every
     * snapshot of length {@code extent.x*extent.y*extent.z});</li>
     * <li>the injected excess {@code S = sum_cell D} is NONZERO -- a zero
     * excess is a run subtracted from itself, for which {@code D == 0}
     * identically and every magnitude is a vacuous {@code 0.0};</li>
     * <li><b>{@code ||D(.,0)||_1 == |S|} -- the SHARED-BACKGROUND
     * check.</b> Before any tick runs the two halves can differ ONLY by
     * the seeded packet, so {@code D(.,0)} is supported on the seed cell
     * alone and its {@code L1} norm equals {@code |S|}. Any background
     * that failed to cancel shows up as extra {@code L1} mass. This runs
     * BEFORE the per-snapshot limbs deliberately: on a genuinely
     * mismatched pair {@code D} carries the whole background difference,
     * which is delocalized, so W2 would otherwise fire first and
     * misdiagnose a broken pair as wrap-saturation -- exactly on the input
     * this check exists to catch;</li>
     * <li>{@code S} is TICK-INVARIANT (invariant (I1)); a drift means the
     * two runs were not driven on the same conserved background, which is
     * a broken matched pair rather than a physics result;</li>
     * <li>per snapshot, wrap-safety limbs W1 and W2 on {@code D} (never on
     * {@code f_P} or {@code f_C}) -- {@link #assertResponseLocalized} --
     * then limb W3, the halo fraction, {@link
     * #assertResponseCoreDominates}.</li>
     * </ol>
     *
     * <p><b>Exactly what check (3) catches, stated without overclaim.</b>
     * By the triangle inequality {@code |S| <= ||D(.,0)||_1} always, with
     * equality iff every nonzero term shares a sign. So the check fires
     * for every background difference that is not single-signed relative
     * to the packet -- which is every difference an RNG-drawn background
     * produces. It does NOT fire for an adversarial difference that is
     * single-signed off-origin and cancels against the origin cell's own
     * difference (e.g. {@code -5} at the origin and {@code +5} at one
     * other cell leaves both {@code ||D||_1} and {@code S} untouched).
     * Neither this check nor the runner's {@code S == 30*packetQuanta}
     * equality closes that class, and their conjunction does not either --
     * both are NORM checks, and that class is constructed to leave both
     * norms fixed.
     *
     * <p><b>That class IS closed, on the campaign path, by a third check
     * that looks at the SUPPORT instead of at a norm</b> (round-2 fix):
     * {@link #assertTickZeroSupportIsTheSeedCell}, called by {@link
     * #runOneSeedMatchedPair} before this method. It requires every
     * non-origin cell of {@code D(.,0)} to be bit-exactly {@code 0.0},
     * which is a total-precision check with no tolerance to tune -- see
     * that method for why the equality is exact rather than approximate.
     * It lives on the runner rather than here because seed-cell support is
     * a property of how the runner builds its pair, not of this method's
     * contract: this method's public API accepts any two field sequences,
     * including a synthetic packet that already has spatial extent at tick
     * 0. So the class is shut for every campaign run, and remains open only
     * for a caller that assembles the two sequences itself.
     *
     * <p>Compare the runner's {@code S == 30*packetQuanta} equality alone,
     * which a mismatched background passes whenever the two background
     * TOTALS happen to coincide -- a per-seed probability of order
     * {@code 1e-3}, i.e. ~1% over an 8-seed campaign. That is why check (3)
     * here, not that one, is the sharp NORM check.
     *
     * @param packetFieldByTick  the run WITH the packet
     * @param controlFieldByTick the run on the SAME RNG background WITHOUT
     *                           the packet
     * @throws IllegalArgumentException if shapes disagree or the injected
     *                                   excess is zero
     * @throws IllegalStateException    if the injected excess drifts, or a
     *                                   snapshot fails W1, W2 or W3
     */
    public static MatchedPairTransport transportEstimateMatchedPair(double[][] packetFieldByTick,
                                                                      double[][] controlFieldByTick,
                                                                      Point3i extent,
                                                                      Point3i originCell) {
        if (packetFieldByTick == null || packetFieldByTick.length < 2) {
            throw new IllegalArgumentException("packetFieldByTick must have at least 2 snapshots, had "
                                                + (packetFieldByTick == null
                                                    ? 0
                                                    : packetFieldByTick.length));
        }
        if (controlFieldByTick == null
            || controlFieldByTick.length != packetFieldByTick.length) {
            throw new IllegalArgumentException("controlFieldByTick must have the same snapshot count as packetFieldByTick ("
                                                + packetFieldByTick.length
                                                + "), had "
                                                + (controlFieldByTick == null
                                                    ? 0
                                                    : controlFieldByTick.length));
        }
        int expected = extent.x * extent.y * extent.z;
        int t = packetFieldByTick.length;
        double[][] difference = new double[t][];
        for (int i = 0; i < t; i++) {
            difference[i] = differenceField(packetFieldByTick[i],
                                             controlFieldByTick[i], expected,
                                             i);
        }

        double injectedExcess = signedTotal(difference[0]);
        if (Math.abs(injectedExcess) <= INJECTED_EXCESS_TOLERANCE) {
            throw new IllegalArgumentException("matched pair has a ZERO injected excess (sum of the difference field at tick 0 is "
                                                + injectedExcess
                                                + "): the two runs carry the same quanta total, so the difference field is identically zero and every"
                                                + " direction magnitude would be a vacuous 0.0 - seed the packet run with a nonzero packetQuanta");
        }

        // The SHARP same-background check, asserted BEFORE the per-tick
        // W1/W2 loop so that a mismatched pair is diagnosed as a
        // mismatched pair rather than misdiagnosed as wrap-saturation
        // (which is what a delocalized background difference looks like to
        // W2). See this method's javadoc for exactly what it does and does
        // not catch.
        double firstL1 = responseL1(difference[0]);
        if (Math.abs(firstL1 - Math.abs(injectedExcess)) > INJECTED_EXCESS_TOLERANCE) {
            throw new IllegalStateException("matched pair is not on a shared background: at tick 0 ||D||_1 = "
                                             + firstL1
                                             + " but |S| = "
                                             + Math.abs(injectedExcess)
                                             + " - before any tick runs the two runs can differ ONLY by the seeded packet, so D must be supported on the"
                                             + " seed cell alone and ||D||_1 must equal |S| exactly; the surplus is background that did NOT cancel, i.e."
                                             + " the two halves were not run on the same RNG background");
        }

        double maxSaturation = 0;
        double maxHalo = 0;
        for (int i = 0; i < t; i++) {
            double total = signedTotal(difference[i]);
            if (Math.abs(total - injectedExcess) > INJECTED_EXCESS_TOLERANCE) {
                throw new IllegalStateException("injected excess is not tick-invariant: tick 0 had "
                                                 + injectedExcess + ", tick "
                                                 + i + " has " + total
                                                 + " - each run conserves its own quanta total exactly (invariant (I1)), so the difference of the two"
                                                 + " totals cannot drift; the two runs were not driven on the same conserved background, which is a"
                                                 + " broken matched pair rather than a measurement");
            }
            maxSaturation = Math.max(maxSaturation,
                                      assertResponseLocalized(difference[i],
                                                               extent,
                                                               originCell, i));
            maxHalo = Math.max(maxHalo,
                                assertResponseCoreDominates(difference[i],
                                                             injectedExcess,
                                                             i));
        }

        return new MatchedPairTransport(momentSlopes(difference, extent,
                                                      originCell),
                                         injectedExcess, firstL1,
                                         responseL1(difference[t - 1]),
                                         maxSaturation, maxHalo);
    }

    /**
     * Elementwise {@code packet - control}, with both snapshots' shapes
     * validated against {@code expected}.
     */
    static double[] differenceField(double[] packet, double[] control,
                                     int expected, int tick) {
        if (packet == null || packet.length != expected) {
            throw new IllegalArgumentException("packetFieldByTick[" + tick
                                                + "] must have length "
                                                + expected);
        }
        if (control == null || control.length != expected) {
            throw new IllegalArgumentException("controlFieldByTick[" + tick
                                                + "] must have length "
                                                + expected);
        }
        double[] difference = new double[expected];
        for (int i = 0; i < expected; i++) {
            difference[i] = packet[i] - control[i];
        }
        return difference;
    }

    /**
     * {@code sum_cell D} -- the injected excess {@code S}, invariant (I1).
     * SIGNED, deliberately: it is the conserved quantity, and its
     * magnitude is the lower bound on {@link #responseL1}.
     */
    static double signedTotal(double[] field) {
        double sum = 0;
        for (double v : field) {
            sum += v;
        }
        return sum;
    }

    /**
     * {@code ||D||_1 = sum_cell |D(cell)|} -- the response's total
     * variation, the weight the matched-pair second moment normalizes by.
     * Bounded below by {@code |signedTotal|} at every tick (class javadoc,
     * "What replaces the maximum principle").
     */
    static double responseL1(double[] field) {
        double sum = 0;
        for (double v : field) {
            sum += Math.abs(v);
        }
        return sum;
    }

    /**
     * The second moment of the UNIFORM distribution over all cells about
     * {@code origin}, along {@code d} -- the value {@code M_d(t)} relaxes
     * onto once the response has filled the box, and hence the denominator
     * of limb W2's saturation ratio. Computed by running {@link
     * #secondMoment} itself over an all-ones field so the numerator and
     * denominator cannot drift apart under any future change to the
     * projection convention.
     *
     * <p>MEMOIZED on {@code (extent, origin, direction)}, which is its
     * complete dependency set -- it does not depend on the field. {@link
     * #momentSaturation} calls it once per direction per snapshot, so
     * without the cache every matched-pair snapshot allocated and swept
     * three extra full {@code L^3} arrays purely to recompute constants.
     * No correctness impact: the cache is keyed on exactly the arguments,
     * and {@code
     * AnisotropyProbeTest#uniformSecondMomentMemoizationIsKeyedOnItsFullDependencySet}
     * checks that varying each of the three independently still yields the
     * value a fresh all-ones sweep gives.
     */
    static double uniformSecondMoment(Point3i extent, Point3i origin,
                                       StructureFactor.Direction d) {
        return UNIFORM_MOMENT_CACHE.computeIfAbsent(new UniformMomentKey(extent.x,
                                                                          extent.y,
                                                                          extent.z,
                                                                          origin.x,
                                                                          origin.y,
                                                                          origin.z,
                                                                          d),
                                                     key -> {
                                                         double[] ones = new double[key.ex()
                                                                                     * key.ey()
                                                                                     * key.ez()];
                                                         Arrays.fill(ones, 1.0);
                                                         return secondMoment(ones,
                                                                              new Point3i(key.ex(),
                                                                                           key.ey(),
                                                                                           key.ez()),
                                                                              new Point3i(key.ox(),
                                                                                           key.oy(),
                                                                                           key.oz()),
                                                                              key.d());
                                                     });
    }

    /** The complete dependency set of {@link #uniformSecondMoment}. */
    private record UniformMomentKey(int ex, int ey, int ez, int ox, int oy,
                                     int oz, StructureFactor.Direction d) {
    }

    private static final Map<UniformMomentKey, Double> UNIFORM_MOMENT_CACHE = new ConcurrentHashMap<>();

    /**
     * Wrap-safety limb W2's measured quantity: the largest over-directions
     * MOMENT SATURATION RATIO {@code M_d(field)/M_d^unif}. {@code 0} for a
     * response still concentrated at the origin, {@code ~1} for one spread
     * uniformly over the box. See the class javadoc, "Limb W2 --
     * saturation", for the derivation, for the calibration of {@link
     * #RESPONSE_MOMENT_SATURATION_TOLERANCE} against slope bias, and for
     * why this replaced the original shell-mass-fraction measure (which
     * read exactly {@code 0.0} on a maximally delocalized interior-uniform
     * field, i.e. was blind to the very case W2 exists for).
     *
     * <p>Returns {@code 0.0} for an all-zero field -- which {@link
     * #transportEstimateMatchedPair} can never pass it, since a zero
     * injected excess is rejected first.
     */
    static double momentSaturation(double[] field, Point3i extent,
                                    Point3i origin) {
        double worst = 0;
        for (StructureFactor.Direction d : StructureFactor.Direction.values()) {
            double uniform = uniformSecondMoment(extent, origin, d);
            if (uniform <= 0) {
                continue;
            }
            worst = Math.max(worst,
                              secondMoment(field, extent, origin, d) / uniform);
        }
        return worst;
    }

    /**
     * Wrap-safety LIMB W3's measured quantity: the DECORRELATION HALO
     * FRACTION
     * <pre>
     *   f_halo(t) = 1 - |S| / ||D(.,t)||_1
     * </pre>
     * the fraction of the response's L1 mass that is NOT accounted for by
     * the conserved injected excess. It is model-free and exact -- a
     * definition, not a fit -- and both of its terms are already computed
     * on this path.
     *
     * <p>{@code f_halo = 0} iff {@code ||D||_1 = |S|}, i.e. every nonzero
     * cell of {@code D} shares a sign: the response is pure transported
     * excess. It rises strictly above zero exactly when the two halves'
     * firing sequences diverge and {@code D} acquires cancelling
     * {@code +/-} structure, since each such pair adds to {@code ||D||_1}
     * while leaving {@code S} untouched. See the class javadoc, "Limb W3",
     * for why a bound on this LEVEL is not a bound on its RATE, and for
     * what the {@code f_halo(0) = 0} anchor makes it nonetheless buy.
     *
     * @param injectedExcess {@code S}; its magnitude is the {@code ||D||_1}
     *                       lower bound, so the returned value is always in
     *                       {@code [0,1)}
     */
    static double haloFraction(double[] difference, double injectedExcess) {
        double l1 = responseL1(difference);
        return l1 > 0 ? 1.0 - Math.abs(injectedExcess) / l1 : 0.0;
    }

    /**
     * Wrap-safety LIMB W3 (halo fraction) -- see {@link #haloFraction} and
     * the class javadoc. Separate from {@link #assertResponseLocalized}
     * deliberately: W1 and W2 are statements about the response's GEOMETRY
     * (how far from the origin its mass sits) and need only the field,
     * whereas W3 is a statement about the response's COMPOSITION (how much
     * of its mass is conserved excess versus decorrelation halo) and needs
     * the conserved {@code S} as well.
     *
     * @return the measured halo fraction, so the caller can report it
     *         rather than only observing that nothing threw
     */
    static double assertResponseCoreDominates(double[] difference,
                                               double injectedExcess,
                                               int tick) {
        double halo = haloFraction(difference, injectedExcess);
        if (halo > RESPONSE_HALO_FRACTION_TOLERANCE) {
            throw new IllegalStateException("halo-fraction bound violated at tick "
                                             + tick
                                             + ": f_halo = 1 - |S|/||D||_1 has reached "
                                             + halo + " (|S| = "
                                             + Math.abs(injectedExcess)
                                             + ", ||D||_1 = "
                                             + responseL1(difference)
                                             + ", tolerance "
                                             + RESPONSE_HALO_FRACTION_TOLERANCE
                                             + ") - that fraction of the L1 mass the second moment is computed over is signed structure created by the two"
                                             + " halves' firing sequences diverging, not transported injected excess, so the moment is no longer measuring"
                                             + " the packet's spread alone; shorten the run or reduce the background amplitude");
        }
        return halo;
    }

    /**
     * The matched-pair path's GEOMETRIC wrap-safety limbs: limb W1
     * (exactness, {@link #assertWrapSafe}, unchanged) THEN limb W2
     * (saturation, {@link #momentSaturation}). Applied to the DIFFERENCE
     * field only -- see the class javadoc, "The periodic-wrap
     * precondition", for why applying either limb to a raw
     * signed-background field is meaningless, and for the proof that W1
     * alone is structurally unreachable on every legal campaign geometry
     * (hence W2's existence).
     *
     * <p><b>This is NOT the full precondition.</b> The matched-pair path's
     * full wrap-safety precondition is W1 THEN W2 THEN W3. Limb W3 (halo
     * fraction) is a statement about the response's COMPOSITION rather than
     * its geometry, needs the conserved {@code S} that this method is not
     * given, and is asserted separately by {@link
     * #assertResponseCoreDominates}. A caller invoking only this method has
     * checked two of the three limbs.
     *
     * @return the measured moment saturation ratio, so the caller can
     *         report it rather than only observing that nothing threw
     */
    static double assertResponseLocalized(double[] difference, Point3i extent,
                                           Point3i origin, int tick) {
        assertWrapSafe(difference, extent, origin, tick);
        double saturation = momentSaturation(difference, extent, origin);
        if (saturation > RESPONSE_MOMENT_SATURATION_TOLERANCE) {
            throw new IllegalStateException("wrap-saturation bound violated at tick "
                                             + tick
                                             + ": the response's second moment has reached "
                                             + saturation
                                             + " of its fully-delocalized value about origin "
                                             + origin + " (extent " + extent
                                             + ", tolerance "
                                             + RESPONSE_MOMENT_SATURATION_TOLERANCE
                                             + ") - the torus bounds the achievable displacement, so the second moment is relaxing onto its uniform-field"
                                             + " asymptote and the OLS slope now UNDERSTATES transport; reduce tick count or enlarge extent");
        }
        return saturation;
    }

    /**
     * The SPECTRAL estimator -- see class javadoc. Consumes {@code sf}'s
     * public, real-field {@link StructureFactor#spectrum} overload
     * directly, per bead inviscid-0nx.9's final-review-verified safe
     * pattern (no manual pre-filtering).
     */
    public static EstimatorResult spectralEstimate(StructureFactor sf,
                                                     double[][] fieldByTick) {
        Map<StructureFactor.Direction, DirectionMagnitude> perDirection = new EnumMap<>(StructureFactor.Direction.class);
        for (StructureFactor.Direction d : StructureFactor.Direction.values()) {
            List<StructureFactor.DispersionPoint> points = sf.spectrum(fieldByTick,
                                                                         d);
            StructureFactor.Ridge ridge = sf.extractRidge(points);
            perDirection.put(d, new DirectionMagnitude(d,
                                                         Math.abs(ridge.slope()),
                                                         points.size()));
        }
        return new EstimatorResult(perDirection, ratio(perDirection));
    }

    /**
     * The shared max/min choke point both estimators reduce through --
     * see class javadoc, "Non-fabrication contract".
     */
    static OptionalDouble ratio(Map<StructureFactor.Direction, DirectionMagnitude> perDirection) {
        double max = Double.NEGATIVE_INFINITY;
        double min = Double.POSITIVE_INFINITY;
        for (DirectionMagnitude dm : perDirection.values()) {
            max = Math.max(max, dm.magnitude());
            min = Math.min(min, dm.magnitude());
        }
        if (min <= RATIO_DEGENERATE_EPSILON) {
            return OptionalDouble.empty();
        }
        return OptionalDouble.of(max / min);
    }

    /**
     * Percentile bootstrap over the per-seed ratios that were present
     * (non-degenerate). Deterministic: the resampling RNG is seeded with
     * the literal {@link #BOOTSTRAP_RNG_SEED}, never wall-clock, per this
     * project's determinism rule. DIAGNOSTIC ONLY -- see {@link
     * BootstrapCi}'s javadoc and class javadoc "STACKED-REVIEW
     * CORRECTION": use {@link #pooledEstimate} for the significance
     * claim.
     *
     * @param presentRatios non-degenerate per-seed ratios
     * @param totalSeeds    the full seed-list size (for reporting {@code
     *                      nSeedsDegenerate = totalSeeds - presentRatios.size()})
     */
    static BootstrapCi bootstrapCi(List<Double> presentRatios, int totalSeeds) {
        int n = presentRatios.size();
        if (n == 0) {
            return new BootstrapCi(Double.NaN, Double.NaN, Double.NaN, 0,
                                    totalSeeds);
        }
        double mean = 0;
        for (double r : presentRatios) {
            mean += r;
        }
        mean /= n;

        Random random = new Random(BOOTSTRAP_RNG_SEED);
        double[] resampleMeans = new double[BOOTSTRAP_RESAMPLES];
        for (int b = 0; b < BOOTSTRAP_RESAMPLES; b++) {
            double sum = 0;
            for (int i = 0; i < n; i++) {
                sum += presentRatios.get(random.nextInt(n));
            }
            resampleMeans[b] = sum / n;
        }
        Arrays.sort(resampleMeans);
        int lowerIdx = (int) (0.025 * BOOTSTRAP_RESAMPLES);
        int upperIdx = Math.min((int) (0.975 * BOOTSTRAP_RESAMPLES),
                                 BOOTSTRAP_RESAMPLES - 1);
        return new BootstrapCi(mean, resampleMeans[lowerIdx],
                                resampleMeans[upperIdx], n,
                                totalSeeds - n);
    }

    // ------------------------------------------------------------------
    // FIX 1 (stacked review): the actual significance statistic --
    // seed-pooled ratio-of-means, resample-then-aggregate bootstrap CI,
    // and a permutation/null-calibration test.
    // ------------------------------------------------------------------

    /**
     * @param perSeedMagnitudes one entry per seed: that seed's magnitude
     *                          for every direction (e.g. {@code
     *                          EstimatorResult.perDirection()} values
     *                          mapped to their {@code magnitude()}).
     * @return the seed-pooled ratio-of-means, its resample-then-aggregate
     *         bootstrap CI, per-direction pooled stats, and a permutation
     *         null-calibration p-value -- see class javadoc,
     *         "STACKED-REVIEW CORRECTION".
     */
    public static PooledResult pooledEstimate(List<Map<StructureFactor.Direction, Double>> perSeedMagnitudes) {
        if (perSeedMagnitudes == null || perSeedMagnitudes.isEmpty()) {
            throw new IllegalArgumentException("perSeedMagnitudes must be non-empty");
        }
        int n = perSeedMagnitudes.size();
        StructureFactor.Direction[] dirs = StructureFactor.Direction.values();

        Map<StructureFactor.Direction, Double> observedMeans = meanPerDirection(perSeedMagnitudes);
        OptionalDouble observedRatio = ratioOfMeans(observedMeans);

        // Resample-then-aggregate bootstrap: resample SEED INDICES
        // jointly across directions (preserves each seed's own
        // cross-direction correlation), recompute per-direction means
        // from the resample, THEN take max/min.
        Random random = new Random(BOOTSTRAP_RNG_SEED);
        Map<StructureFactor.Direction, double[]> resampleMeans = new EnumMap<>(StructureFactor.Direction.class);
        for (StructureFactor.Direction d : dirs) {
            resampleMeans.put(d, new double[BOOTSTRAP_RESAMPLES]);
        }
        double[] resampleRatios = new double[BOOTSTRAP_RESAMPLES];
        int validRatios = 0;
        for (int b = 0; b < BOOTSTRAP_RESAMPLES; b++) {
            double[] sums = new double[dirs.length];
            for (int i = 0; i < n; i++) {
                Map<StructureFactor.Direction, Double> seedMap = perSeedMagnitudes.get(random.nextInt(n));
                for (int di = 0; di < dirs.length; di++) {
                    sums[di] += seedMap.get(dirs[di]);
                }
            }
            double max = Double.NEGATIVE_INFINITY;
            double min = Double.POSITIVE_INFINITY;
            for (int di = 0; di < dirs.length; di++) {
                double mean = sums[di] / n;
                resampleMeans.get(dirs[di])[b] = mean;
                max = Math.max(max, mean);
                min = Math.min(min, mean);
            }
            if (min > RATIO_DEGENERATE_EPSILON) {
                resampleRatios[validRatios++] = max / min;
            }
        }

        Map<StructureFactor.Direction, PooledDirectionStats> perDirectionStats = new EnumMap<>(StructureFactor.Direction.class);
        for (StructureFactor.Direction d : dirs) {
            double[] arr = resampleMeans.get(d).clone();
            Arrays.sort(arr);
            double lower = arr[(int) (0.025 * BOOTSTRAP_RESAMPLES)];
            double upper = arr[Math.min((int) (0.975 * BOOTSTRAP_RESAMPLES),
                                         BOOTSTRAP_RESAMPLES - 1)];
            perDirectionStats.put(d, new PooledDirectionStats(d,
                                                                observedMeans.get(d),
                                                                lower, upper,
                                                                n));
        }

        double pooledRatioCiLower = Double.NaN;
        double pooledRatioCiUpper = Double.NaN;
        if (validRatios > 0) {
            double[] ratios = Arrays.copyOf(resampleRatios, validRatios);
            Arrays.sort(ratios);
            pooledRatioCiLower = ratios[(int) (0.025 * ratios.length)];
            pooledRatioCiUpper = ratios[Math.min((int) (0.975 * ratios.length),
                                                  ratios.length - 1)];
        }

        double permutationPValue = Double.NaN;
        double permutationNull95 = Double.NaN;
        int permutationCount = 0;
        if (observedRatio.isPresent()) {
            double[] nulls = permutationNullDistribution(perSeedMagnitudes,
                                                           PERMUTATION_COUNT,
                                                           PERMUTATION_RNG_SEED);
            permutationCount = nulls.length;
            // guards the +1 correction too (inviscid-0sn): without this,
            // permutationCount==0 would compute (0+1)/(0+1)=1.0, a NEW
            // fabricated-p=1.0 failure mode the correction could
            // introduce, not just the original divide-by-zero.
            if (permutationCount > 0) {
                double observed = observedRatio.getAsDouble();
                long countGe = countGe(nulls, observed);
                // +1 continuity correction (inviscid-0sn): under H0 the
                // observed statistic is exchangeable with the null draws,
                // so it counts as one of its own reference set -- avoids
                // a fabricated exact-zero p-value that a finite
                // permutation count can never actually prove. Maximum
                // shift from the correction is bounded by
                // 1/(permutationCount+1) (~5e-4 at N=2000), well inside
                // the matched-noise control tests' margins to their 0.05
                // threshold (~0.36 and ~0.048 -- see
                // permutationTestIsNonSignificantForMatchedNoiseIsotropicControl
                // / permutationTestDetectsGenuineTwoFoldAnisotropyAtMatchedNoise).
                permutationPValue = (countGe + 1) / (double) (permutationCount + 1);
                double[] sorted = nulls.clone();
                Arrays.sort(sorted);
                int idx95 = Math.min((int) (0.95 * sorted.length),
                                      sorted.length - 1);
                permutationNull95 = sorted[idx95];
            }
        }

        return new PooledResult(perDirectionStats, observedRatio,
                                 pooledRatioCiLower, pooledRatioCiUpper,
                                 permutationPValue, permutationNull95,
                                 permutationCount);
    }

    /**
     * The null-calibration draw: for each of {@code permutations} trials,
     * independently shuffle EACH seed's own 3-direction magnitude triple
     * (Fisher-Yates on 3 elements), pool the shuffled values into
     * per-direction means across all seeds, and record {@code max/min}.
     * Shuffling WITHIN each seed (not across seeds) is deliberate: it
     * destroys exactly the "does this magnitude belong to X100, X110, or
     * X111" information the significance test is about, while preserving
     * each seed's own noise realization (its multiset of 3 values) --
     * the correct null for "is there a stable per-direction effect", not
     * a test of "do the seeds differ from each other" (a different
     * question). Trials whose shuffled pooled minimum is degenerate
     * (<= {@link #RATIO_DEGENERATE_EPSILON}) are dropped, not counted as
     * zero or infinity -- {@link PooledResult#permutationCount()}
     * reports how many trials were actually usable.
     */
    static double[] permutationNullDistribution(List<Map<StructureFactor.Direction, Double>> perSeedMagnitudes,
                                                  int permutations, long rngSeed) {
        Random random = new Random(rngSeed);
        int n = perSeedMagnitudes.size();
        StructureFactor.Direction[] dirs = StructureFactor.Direction.values();
        double[] nulls = new double[permutations];
        int kept = 0;
        for (int p = 0; p < permutations; p++) {
            double[] sums = new double[dirs.length];
            for (Map<StructureFactor.Direction, Double> seedMap : perSeedMagnitudes) {
                double[] vals = new double[dirs.length];
                for (int i = 0; i < dirs.length; i++) {
                    vals[i] = seedMap.get(dirs[i]);
                }
                shuffle(vals, random);
                for (int i = 0; i < dirs.length; i++) {
                    sums[i] += vals[i];
                }
            }
            double max = Double.NEGATIVE_INFINITY;
            double min = Double.POSITIVE_INFINITY;
            for (double s : sums) {
                double mean = s / n;
                max = Math.max(max, mean);
                min = Math.min(min, mean);
            }
            if (min > RATIO_DEGENERATE_EPSILON) {
                nulls[kept++] = max / min;
            }
        }
        return Arrays.copyOf(nulls, kept);
    }

    /**
     * How many of {@code nulls} meet or exceed {@code observed} -- the
     * numerator ingredient of the permutation p-value's +1 continuity
     * correction (inviscid-0sn). Ties count TOWARD the null (conservative
     * standard practice, consistent with the +1 correction's own
     * exchangeability rationale): an exact tie is treated as "at least as
     * extreme as observed", not excluded. Extracted as its own
     * package-private method so the counting rule is directly testable
     * without needing to force an exact tie through the full
     * shuffle-based {@link #permutationNullDistribution} + campaign
     * fixture (ties have ~0 probability under continuous noise).
     */
    static long countGe(double[] nulls, double observed) {
        long countGe = 0;
        for (double v : nulls) {
            if (v >= observed) {
                countGe++;
            }
        }
        return countGe;
    }

    private static void shuffle(double[] arr, Random random) {
        for (int i = arr.length - 1; i > 0; i--) {
            int j = random.nextInt(i + 1);
            double tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    private static Map<StructureFactor.Direction, Double> meanPerDirection(List<Map<StructureFactor.Direction, Double>> perSeedMagnitudes) {
        Map<StructureFactor.Direction, Double> sums = new EnumMap<>(StructureFactor.Direction.class);
        for (StructureFactor.Direction d : StructureFactor.Direction.values()) {
            sums.put(d, 0.0);
        }
        for (Map<StructureFactor.Direction, Double> seedMap : perSeedMagnitudes) {
            for (StructureFactor.Direction d : StructureFactor.Direction.values()) {
                sums.merge(d, seedMap.get(d), Double::sum);
            }
        }
        int n = perSeedMagnitudes.size();
        Map<StructureFactor.Direction, Double> means = new EnumMap<>(StructureFactor.Direction.class);
        for (Map.Entry<StructureFactor.Direction, Double> e : sums.entrySet()) {
            means.put(e.getKey(), e.getValue() / n);
        }
        return means;
    }

    private static OptionalDouble ratioOfMeans(Map<StructureFactor.Direction, Double> means) {
        double max = Double.NEGATIVE_INFINITY;
        double min = Double.POSITIVE_INFINITY;
        for (double v : means.values()) {
            max = Math.max(max, v);
            min = Math.min(min, v);
        }
        if (min <= RATIO_DEGENERATE_EPSILON) {
            return OptionalDouble.empty();
        }
        return OptionalDouble.of(max / min);
    }

    private static Map<StructureFactor.Direction, Double> magnitudesOf(EstimatorResult result) {
        Map<StructureFactor.Direction, Double> magnitudes = new EnumMap<>(StructureFactor.Direction.class);
        for (Map.Entry<StructureFactor.Direction, DirectionMagnitude> e : result.perDirection()
                                                                                 .entrySet()) {
            magnitudes.put(e.getKey(), e.getValue().magnitude());
        }
        return magnitudes;
    }

    private static double olsSlope(double[] xs, double[] ys) {
        int n = xs.length;
        double sumX = 0;
        double sumY = 0;
        for (int i = 0; i < n; i++) {
            sumX += xs[i];
            sumY += ys[i];
        }
        double xMean = sumX / n;
        double yMean = sumY / n;
        double num = 0;
        double den = 0;
        for (int i = 0; i < n; i++) {
            double dx = xs[i] - xMean;
            num += dx * (ys[i] - yMean);
            den += dx * dx;
        }
        return den > 1e-12 ? num / den : 0.0;
    }

    private static double secondMoment(double[] field, Point3i extent,
                                        Point3i origin,
                                        StructureFactor.Direction d) {
        double totalMass = 0;
        double weighted = 0;
        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    double mass = Math.abs(field[(i * extent.y + j) * extent.z
                                                  + k]);
                    if (mass == 0.0) {
                        continue;
                    }
                    double proj = projectionAlong(i - origin.x, j - origin.y,
                                                   k - origin.z, d);
                    totalMass += mass;
                    weighted += mass * proj * proj;
                }
            }
        }
        return totalMass > 0 ? weighted / totalMass : 0.0;
    }

    private static double projectionAlong(int dx, int dy, int dz,
                                           StructureFactor.Direction d) {
        switch (d) {
        case X100:
            return dx;
        case X110:
            return (dx + dy) / Math.sqrt(2);
        case X111:
            return (dx + dy + dz) / Math.sqrt(3);
        default:
            throw new IllegalArgumentException("unhandled direction: " + d);
        }
    }

    /**
     * Wrap-safety LIMB W1 (exactness) -- see the class javadoc, "The
     * periodic-wrap precondition: two limbs". ORIGIN-RELATIVE and
     * mathematically exact: a cell contributes "violating mass" iff its
     * per-axis displacement from {@code origin} STRICTLY exceeds that
     * axis's half-period ({@code extentAxis/2}) -- the exact point past
     * which the naive unwrapped distance would overestimate the true
     * periodic minimum-image distance. The criterion is unchanged by bead
     * inviscid-0nx.28's re-derivation; what that bead established is its
     * REACH.
     *
     * <p><b>Reach, stated honestly.</b> This limb is a real, reachable
     * backstop for the general any-origin public API ({@link
     * #transportEstimate} accepts any {@code originCell}), and is
     * PROVABLY UNREACHABLE for every legal campaign geometry, on every
     * field -- signed or not, localized or not -- because {@link
     * FccNeighborhood} forces even extents and {@link
     * #nearestEvenParityCenter} then places the origin so that the maximum
     * attainable {@code |coord-origin|} is exactly {@code extentAxis/2},
     * never more. It therefore certifies NOTHING about a campaign run on
     * its own, and the matched-pair path pairs it with limb W2 ({@link
     * #assertResponseLocalized}), which is the limb that can actually
     * bite there.
     *
     * <p>Deliberately NOT strengthened in place: this method is on the
     * bit-for-bit-pinned {@code {0, packetQuanta}} path ({@code
     * SeamGoldenCompatTest}, {@code SubstrateFactorySeamTest}), so W2 is
     * added alongside it rather than folded into it.
     */
    static void assertWrapSafe(double[] field, Point3i extent, Point3i origin,
                                int tick) {
        double totalMass = 0;
        double violatingMass = 0;
        int halfX = extent.x / 2;
        int halfY = extent.y / 2;
        int halfZ = extent.z / 2;
        for (int i = 0; i < extent.x; i++) {
            int dx = Math.abs(i - origin.x);
            for (int j = 0; j < extent.y; j++) {
                int dy = Math.abs(j - origin.y);
                for (int k = 0; k < extent.z; k++) {
                    double mass = Math.abs(field[(i * extent.y + j)
                                                  * extent.z + k]);
                    if (mass == 0.0) {
                        continue;
                    }
                    totalMass += mass;
                    int dz = Math.abs(k - origin.z);
                    boolean violates = dx > halfX || dy > halfY || dz > halfZ;
                    if (violates) {
                        violatingMass += mass;
                    }
                }
            }
        }
        if (totalMass > 0 && violatingMass > 0) {
            throw new IllegalStateException("wrap-safety bound violated at tick "
                                             + tick + ": " + violatingMass
                                             + "/" + totalMass
                                             + " of quanta mass has, on at least one axis, moved past the exact half-period distance from origin "
                                             + origin + " (extent " + extent
                                             + ") - the naive Cartesian second-moment is invalid past this point;"
                                             + " reduce tick count, enlarge extent, or use a centered origin so the packet stays within the half-period bound");
        }
    }

    // ------------------------------------------------------------------
    // Real-automaton harness (the Phase A campaign) -- reuses the exact
    // HybridAutomaton/CollisionSweep/QuantaExchangeRule/ConservationAudit/
    // AuditedRun wiring ContactAtlasGenerator#runDynamicReachability
    // established, not reimplemented.
    // ------------------------------------------------------------------

    /**
     * Drives one seed's automaton run: random angles (seeded), a
     * localized quanta packet at {@code originCell}, {@code ticks}
     * recorded {@link StructureFactor#coarseGrainedField} snapshots (the
     * first BEFORE any tick runs, matching {@code
     * BaselineSpectrumHarness}'s "first recorded value is before any
     * ticks" convention), then both estimators computed from that SAME
     * snapshot sequence, plus (FIX 2) this seed's total/effective
     * collision counts.
     */
    public static SeedResult runOneSeed(Point3i extent, long seed, int ticks,
                                         int packetQuanta, Point3i originCell) {
        return runOneSeed(extent, seed, ticks, packetQuanta, originCell,
                           AnisotropyProbe::phaseAHybridSubstrate);
    }

    /**
     * Substrate-injectable overload (bead inviscid-ckn / inviscid-0nx.21,
     * T2 design-ckn-lattice-seam.md §4): identical to the 5-arg {@link
     * #runOneSeed(Point3i, long, int, int, Point3i)} except the substrate
     * is built by {@code factory} instead of being hardwired to the
     * Phase A hybrid.
     */
    public static SeedResult runOneSeed(Point3i extent, long seed, int ticks,
                                         int packetQuanta, Point3i originCell,
                                         SubstrateFactory factory) {
        SubstrateFactory.Substrate substrate = factory.create(extent, seed,
                                                                packetQuanta,
                                                                originCell);

        double[][] fieldByTick = recordFieldByTick(substrate, ticks);

        StructureFactor sf = new StructureFactor(extent);
        EstimatorResult transport = transportEstimate(fieldByTick, extent,
                                                        originCell);
        EstimatorResult spectral = spectralEstimate(sf, fieldByTick);
        return new SeedResult(seed, transport, spectral,
                               substrate.statistics().totalCollisions(),
                               substrate.statistics().effectiveCollisions());
    }

    /**
     * The snapshot loop, extracted verbatim from {@link #runOneSeed} so
     * the matched-pair runner cannot drift from it: the first snapshot is
     * recorded BEFORE any tick runs (matching {@code
     * BaselineSpectrumHarness}'s convention), then one per tick.
     */
    private static double[][] recordFieldByTick(SubstrateFactory.Substrate substrate,
                                                 int ticks) {
        double[][] fieldByTick = new double[ticks][];
        fieldByTick[0] = StructureFactor.coarseGrainedField(substrate.field());
        for (int tick = 1; tick < ticks; tick++) {
            substrate.run().tick(tick - 1);
            fieldByTick[tick] = StructureFactor.coarseGrainedField(substrate.field());
        }
        return fieldByTick;
    }

    /**
     * Drives ONE matched pair (bead inviscid-0nx.28; T2 {@code
     * design-seeding-radius.md} §D-A): the same {@code factory}, the same
     * {@code seed}, run twice -- once with {@code packetQuanta} and once
     * with {@code 0} -- and the two snapshot sequences fed to {@link
     * #transportEstimateMatchedPair}.
     *
     * <p><b>Why {@code packetQuanta = 0} is the control.</b> {@link
     * #phaseAHybridSubstrate} draws randomness ONLY in {@code
     * seedRandomAngles}; seeding the packet consumes no RNG at all. The
     * control run therefore replays the packet run's background draw for
     * draw, which is exactly the "same RNG background, one run with the
     * packet and one without" of §D-A. This method requires the same of
     * any injected {@link SubstrateFactory}: <b>a factory's background
     * draws must not depend on {@code packetQuanta}</b>. That requirement
     * is not left to trust -- the assertion below is its check.
     *
     * <p><b>The same-background assertion, and its ORDER.</b> The injected
     * excess must come out as exactly {@code 30 * packetQuanta} (30
     * members per cell, the {@code Necronomata} layout {@code
     * StructureFactor.coarseGrainedField} sums over). If a factory draws a
     * different background for the two halves, both runs still conserve
     * their own totals, so the tick-invariance check inside {@link
     * #transportEstimateMatchedPair} would pass while the pair is silently
     * mismatched. This equality is evaluated <b>before</b> {@link
     * #transportEstimateMatchedPair} is called at all: run the other way
     * round (as the first version of this method did), a mismatched pair
     * -- whose difference field carries the whole delocalized background
     * difference -- trips wrap-safety limb W2 first and is reported as
     * "reduce tick count or enlarge extent", the wrong root cause on
     * precisely the input the check exists for. Pinned by {@code
     * AnisotropyProbeTest#matchedPairRunnerDiagnosesAMismatchedBackgroundBeforeWrapSaturation}.
     *
     * <p>This equality is the BLUNTEST of the three same-background checks
     * the runner path now applies, and each catches something the others
     * cannot:
     * <ul>
     * <li><b>this one</b> -- two different backgrounds whose totals
     * coincide pass it (order {@code 1e-3} per seed). What it alone adds
     * is the packet MAGNITUDE: it catches a factory that shares a
     * background correctly but injects the wrong excess;</li>
     * <li><b>{@code ||D(.,0)||_1 == |S|}</b> inside {@link
     * #transportEstimateMatchedPair} -- the sharp NORM check, which fires
     * on every background difference that is not single-signed relative to
     * the packet;</li>
     * <li><b>{@link #assertTickZeroSupportIsTheSeedCell}</b>, called below
     * -- the SUPPORT check, which closes the single-signed zero-sum class
     * ({@code -5} at the origin, {@code +5} elsewhere) that is constructed
     * to leave both norms untouched and so defeats both of the above.</li>
     * </ul>
     *
     * @throws IllegalStateException if the two runs did not share a
     *                                background (injected excess is not
     *                                {@code 30 * packetQuanta}, or the
     *                                tick-0 difference field carries
     *                                uncancelled background in either its
     *                                norm or its support), or if a
     *                                snapshot fails wrap-safety limb W1,
     *                                W2 or W3
     */
    public static MatchedPairTransport runOneSeedMatchedPair(Point3i extent,
                                                               long seed,
                                                               int ticks,
                                                               int packetQuanta,
                                                               Point3i originCell,
                                                               SubstrateFactory factory) {
        double[][] packetFieldByTick = recordFieldByTick(factory.create(extent,
                                                                         seed,
                                                                         packetQuanta,
                                                                         originCell),
                                                          ticks);
        double[][] controlFieldByTick = recordFieldByTick(factory.create(extent,
                                                                          seed,
                                                                          0,
                                                                          originCell),
                                                           ticks);
        // FIRST, before the estimator runs: a mismatched pair's difference
        // field is the whole (delocalized) background difference, which
        // would trip wrap-safety limb W2 and be misreported as saturation.
        double expectedExcess = 30.0 * packetQuanta;
        double[] tickZero = differenceField(packetFieldByTick[0],
                                              controlFieldByTick[0],
                                              extent.x * extent.y * extent.z,
                                              0);
        double actualExcess = signedTotal(tickZero);
        if (Math.abs(actualExcess - expectedExcess) > INJECTED_EXCESS_TOLERANCE) {
            throw new IllegalStateException("matched pair is not on a shared background: injected excess was "
                                             + actualExcess
                                             + ", expected 30*packetQuanta="
                                             + expectedExcess
                                             + " - the substrate factory's background draws must not depend on packetQuanta (see this method's javadoc"
                                             + " and SubstrateFactory's RNG draw-order contract)");
        }
        assertTickZeroSupportIsTheSeedCell(tickZero, extent, originCell);
        return transportEstimateMatchedPair(packetFieldByTick,
                                             controlFieldByTick, extent,
                                             originCell);
    }

    /**
     * The SHARPEST same-background check available on the runner path, and
     * the one that closes the adversarial class the two norm-based checks
     * cannot (bead inviscid-0nx.28, round-2 fix): <b>the tick-0 difference
     * field's nonzero support must be exactly the seed cell.</b>
     *
     * <p>{@code ||D(.,0)||_1 == |S|} fires only on a background difference
     * that is not single-signed relative to the packet; the runner's
     * {@code S == 30*packetQuanta} fires only when the two background
     * TOTALS differ. A difference of {@code -5} at the origin and
     * {@code +5} at one other cell leaves BOTH untouched. This check
     * refuses it, because it looks at the support rather than at either
     * norm.
     *
     * <p><b>Why it is a TOTAL-PRECISION check with no tolerance to tune.</b>
     * Before any tick runs, the two halves come from the same factory at
     * the same seed and differ only by {@link #seedPacket}, which touches
     * the origin cell's 30 members and nothing else. {@link
     * StructureFactor#coarseGrainedField} then sums, per cell, in a fixed
     * order, the SAME 30 {@code float}s on both sides of every non-origin
     * cell -- so every non-origin entry of {@code D(.,0)} is a
     * {@code double} minus itself, i.e. BIT-EXACTLY {@code 0.0}. The
     * comparison is {@code != 0.0}, not an epsilon.
     *
     * <p>It lives on the RUNNER rather than inside {@link
     * #transportEstimateMatchedPair} because it is a property of how the
     * runner builds its pair, not of the estimator's contract: the
     * estimator's public API accepts any two field sequences, including a
     * synthetic packet that already has spatial extent at tick 0 (the
     * {@code 40^3} Gaussian fixture in {@code
     * AnisotropyProbeTest#matchedPairRecoversPacketSpreadThroughASignedDriftingBackground}
     * has 64000 strictly-nonzero tick-0 cells, measured), for which
     * seed-cell support is false by construction and correctly so.
     *
     * @throws IllegalStateException if any non-origin cell of the tick-0
     *                                difference is nonzero
     */
    static void assertTickZeroSupportIsTheSeedCell(double[] tickZeroDifference,
                                                    Point3i extent,
                                                    Point3i originCell) {
        int seed = (originCell.x * extent.y + originCell.y) * extent.z
                    + originCell.z;
        for (int i = 0; i < tickZeroDifference.length; i++) {
            if (i != seed && tickZeroDifference[i] != 0.0) {
                int k = i % extent.z;
                int j = (i / extent.z) % extent.y;
                int x = i / (extent.y * extent.z);
                throw new IllegalStateException("matched pair is not on a shared background: at tick 0 the difference field is nonzero at cell ("
                                                 + x + ", " + j + ", " + k
                                                 + ") = " + tickZeroDifference[i]
                                                 + ", away from the seed cell "
                                                 + originCell
                                                 + " - before any tick runs the two halves differ ONLY by the seeded packet, so every non-origin cell must be"
                                                 + " bit-exactly 0.0; a nonzero one is background that did not cancel. This catches the single-signed,"
                                                 + " zero-sum background difference (e.g. -5 at the origin and +5 elsewhere) that leaves both S and ||D||_1"
                                                 + " untouched and so passes every norm-based check");
            }
        }
    }

    /**
     * The Phase A default {@link SubstrateFactory}: holds the prior
     * lines 875-889 VERBATIM, same order (bead inviscid-ckn /
     * inviscid-0nx.21, T2 design-ckn-lattice-seam.md §4) -- so Phase A
     * reproducibility is preserved by construction, not by
     * re-measurement. <b>RNG draw order is part of the contract:</b>
     * {@link #seedRandomAngles} then {@link #seedPacket}, exactly as
     * before -- a factory that reorders those draws changes every
     * seeded trajectory even at the same seed (pinned by {@code
     * SeamGoldenCompatTest#runOneSeedThroughTheSeamMatchesPinnedPhaseANumerics}
     * and {@code
     * SubstrateFactorySeamTest#fiveArgRunOneSeedDelegatesToPhaseAHybridSubstrateFactoryUnchanged}).
     */
    public static SubstrateFactory.Substrate phaseAHybridSubstrate(Point3i extent,
                                                                     long seed,
                                                                     int packetQuanta,
                                                                     Point3i originCell) {
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, seed);
        seedPacket(automaton, originCell, packetQuanta);

        FccNeighborhood neighborhood = new FccNeighborhood(automaton);
        ContactPredicate predicate = new ContactPredicate(new MemberGeometry(ContactAtlasGenerator.GEOMETRY_RESOLUTION,
                                                                               ContactAtlasGenerator.RADIUS));
        ContactScan scan = new ContactScan(automaton, neighborhood, predicate);
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);
        ConservationAudit audit = new ConservationAudit(automaton);
        AuditedRun run = new AuditedRun(hybrid, audit);

        return new SubstrateFactory.Substrate(automaton, run, statistics);
    }

    /**
     * The Phase A campaign: one {@link #runOneSeed} per seed, aggregated
     * into both estimators' naive per-seed {@link BootstrapCi} (diagnostic)
     * AND the pooled/null-calibrated {@link PooledResult} (the
     * significance statistic -- see class javadoc, "STACKED-REVIEW
     * CORRECTION"). Wall-time-budgeted for manual/main() invocation --
     * NOT run inside surefire (see {@link #main(String[])}).
     */
    public static Report runCampaign(Point3i extent, long[] seeds, int ticks,
                                      int packetQuanta) {
        return runCampaign(extent, seeds, ticks, packetQuanta,
                            AnisotropyProbe::phaseAHybridSubstrate);
    }

    /**
     * Substrate-injectable overload (bead inviscid-ckn / inviscid-0nx.21,
     * T2 design-ckn-lattice-seam.md §4) -- same overload pattern as
     * {@link #runOneSeed(Point3i, long, int, int, Point3i, SubstrateFactory)}.
     * {@code nearestEvenParityCenter(extent)} origin selection stays
     * outside the factory (campaign geometry, not substrate
     * construction, per the design memo).
     */
    public static Report runCampaign(Point3i extent, long[] seeds, int ticks,
                                      int packetQuanta,
                                      SubstrateFactory factory) {
        Point3i origin = nearestEvenParityCenter(extent);
        List<SeedResult> perSeed = new ArrayList<>(seeds.length);
        List<Double> transportRatios = new ArrayList<>();
        List<Double> spectralRatios = new ArrayList<>();
        List<Map<StructureFactor.Direction, Double>> transportMagnitudes = new ArrayList<>();
        List<Map<StructureFactor.Direction, Double>> spectralMagnitudes = new ArrayList<>();
        for (long seed : seeds) {
            SeedResult result = runOneSeed(extent, seed, ticks, packetQuanta,
                                            origin, factory);
            perSeed.add(result);
            result.transport().ratio().ifPresent(transportRatios::add);
            result.spectral().ratio().ifPresent(spectralRatios::add);
            transportMagnitudes.add(magnitudesOf(result.transport()));
            spectralMagnitudes.add(magnitudesOf(result.spectral()));
        }
        BootstrapCi transportCi = bootstrapCi(transportRatios, seeds.length);
        BootstrapCi spectralCi = bootstrapCi(spectralRatios, seeds.length);
        PooledResult pooledTransport = pooledEstimate(transportMagnitudes);
        PooledResult pooledSpectral = pooledEstimate(spectralMagnitudes);
        return new Report(perSeed, transportCi, spectralCi, pooledTransport,
                           pooledSpectral, extent, ticks, origin, seeds,
                           packetQuanta);
    }

    /**
     * @return the even-parity cell nearest {@code extent}'s geometric
     *         center -- {@code floor(extent/2)} per axis, with {@code z}
     *         decremented by one if that lands on an odd-parity index
     *         (guaranteed room since {@code extent} axes are each &gt;=4
     *         per {@link FccNeighborhood}'s precondition).
     */
    static Point3i nearestEvenParityCenter(Point3i extent) {
        int cx = extent.x / 2;
        int cy = extent.y / 2;
        int cz = extent.z / 2;
        if (((cx + cy + cz) & 1) != 0) {
            cz -= 1;
        }
        return new Point3i(cx, cy, cz);
    }

    /**
     * Seeds the transport packet at {@code originCell}.
     *
     * <p><b>ACCUMULATES ({@code +=}); it does not assign.</b> This is a
     * deliberate, load-bearing decision taken by bead inviscid-0nx.28's
     * fix round, not an incidental style choice, and E.2 builds directly
     * on it.
     *
     * <p>On the Phase A {@code {0, packetQuanta}} substrate the two are
     * indistinguishable: {@link #seedRandomAngles} writes ANGLES only, a
     * fresh {@link Necronomata}'s {@code frequency} array is all zeros, and
     * {@code 0 + packetQuanta == packetQuanta}. The change is therefore
     * bit-for-bit inert on the pinned path (which {@code
     * SeamGoldenCompatTest} and {@code SubstrateFactorySeamTest} enforce
     * without modification). They diverge the moment a factory supplies a
     * NONZERO background -- §D-A's uniform signed background, the next
     * thing to be built:
     * <ul>
     * <li><b>{@code +=} (chosen)</b> makes the packet an EXCESS on top of
     * the background. The control half then differs from the packet half
     * by exactly {@code 30 * packetQuanta} at the origin cell and by
     * nothing anywhere else, which is precisely what the matched-pair
     * derivation assumes and what both same-background checks
     * ({@code ||D(.,0)||_1 == |S|} and {@code S == 30*packetQuanta}) are
     * written against.</li>
     * <li><b>{@code =} (rejected)</b> would OVERWRITE the background at
     * the origin cell, making the packet half "background with a hole at
     * the origin, refilled with packetQuanta". The injected excess would
     * then be {@code 30*packetQuanta - background(origin)}, so the runner's
     * equality would fire spuriously on a perfectly well-formed pair, and
     * the physical object being measured would not be the localized excess
     * §D-A specifies.</li>
     * </ul>
     * Pinned on a pre-loaded (nonzero-background) automaton by {@code
     * AnisotropyProbeTest#seedPacketAccumulatesOntoTheBackgroundRatherThanOverwritingIt}
     * -- which is why this method is package-private rather than private.
     *
     * <p>NOTE for E.2: {@code PhaseCMeasurement.seedPacket} carries the
     * same {@code =} today. It is equally inert there right now (that
     * factory's {@code seedRandomPhases} writes phases only, so its quanta
     * also start at zero) and it is out of this bead's two-file scope, but
     * it is the same latent trap and must be converted in lock-step with
     * any signed-background substrate.
     */
    static void seedPacket(Necronomata automaton, Point3i originCell,
                            int packetQuanta) {
        int base = automaton.indexOfCell(originCell);
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int m = 0; m < 30; m++) {
                frequency[base + m] += packetQuanta;
            }
        });
    }

    /**
     * Mirrors {@code ContactAtlasGenerator}'s private {@code
     * seedRandomAngles} (that method is package-private to {@code lga}
     * and cannot be imported into {@code measure} without widening its
     * visibility -- mirrored, not reused, per the relay's instruction).
     */
    private static void seedRandomAngles(Necronomata automaton,
                                          Point3i extent, long seed) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[length];
        for (int i = 0; i < length; i++) {
            angles[i] = random.nextFloat() * (float) (2 * Math.PI);
        }
        automaton.process((angleArray, frequency, deltaA,
                            deltaF) -> System.arraycopy(angles, 0, angleArray,
                                                         0, length));
    }

    // ------------------------------------------------------------------
    // Provenance / golden-artifact convention -- mirrors
    // BaselineSpectrumHarness/ContactAtlasGenerator.
    // ------------------------------------------------------------------

    /**
     * Regenerates the Phase A report with the default parameters and
     * overwrites the committed artifact. Run manually (IDE/classpath
     * invocation -- no exec plugin is configured in this project); the
     * regenerated file must then be reviewed and committed by hand. NOT
     * invoked by surefire -- see {@code AnisotropyProbeTest}'s
     * committed-artifact structural validation for what surefire actually
     * checks.
     */
    public static void main(String[] args) throws IOException {
        if (args.length > 0 && "--census".equals(args[0])) {
            printPhaseANonNegativityCensus();
            return;
        }
        long start = System.nanoTime();
        Report report = runCampaign(DEFAULT_EXTENT, DEFAULT_SEEDS,
                                     DEFAULT_TICKS, DEFAULT_PACKET_QUANTA);
        double wallSeconds = (System.nanoTime() - start) / 1e9;
        String tsv = toTsv(report, resolveGitCommit());
        Path path = Paths.get(GOLDEN_RELATIVE_PATH);
        Files.createDirectories(path.getParent());
        Files.write(path, tsv.getBytes(StandardCharsets.UTF_8));
        System.out.println("Wrote " + path.toAbsolutePath() + " in "
                            + wallSeconds + "s");
        System.out.println("transport (naive per-seed, DIAGNOSTIC ONLY): "
                            + report.transportCi());
        System.out.println("spectral  (naive per-seed, DIAGNOSTIC ONLY): "
                            + report.spectralCi());
        System.out.println("transport POOLED (significance statistic): "
                            + report.pooledTransport());
        System.out.println("spectral  POOLED (significance statistic): "
                            + report.pooledSpectral());
    }

    /**
     * The FULL Phase A non-negativity census -- the out-of-suite
     * measurement the class javadoc's {@code {0, packetQuanta}} paragraph
     * and limb W2's campaign-headroom claim both rest on, made
     * REPRODUCIBLE FROM THIS TREE (round-2 fix: previously it existed only
     * as a number quoted in prose, with no artifact, no hook and no
     * {@code @Ignore}d test, so a reader had no way to re-derive it).
     * Costs minutes -- which is why the suite pins a reduced tripwire
     * instead -- so it is a {@link #main} argument rather than a test,
     * following the same precedent as the committed TSV artifact.
     *
     * <p>Counts, over all {@code seeds x ticks} snapshots: negative MEMBER
     * and negative CELL observations, the minima of both, the worst W2
     * moment saturation, the worst W3 halo fraction (exactly {@code 0}
     * here -- on this substrate the control run holds no quanta, so
     * {@code D} is the packet field and is non-negative), and the coverage
     * the in-suite tripwire achieves, including the at-risk
     * ({@code q == 1}) population broken out by 16-tick bucket so that the
     * DIRECTION of the tripwire's bias is visible rather than asserted.
     *
     * <p>Prints only -- it asserts nothing and writes no artifact. The
     * assertions live in {@code
     * AnisotropyProbeTest#phaseARegimeIsNonNegativeAndFarFromSaturationByCensus}
     * at reduced scope; this is the instrument that produced the numbers
     * those javadocs quote.
     */
    static void printPhaseANonNegativityCensus() {
        Point3i extent = DEFAULT_EXTENT;
        Point3i origin = nearestEvenParityCenter(extent);
        int ticks = DEFAULT_TICKS;
        long[] seeds = DEFAULT_SEEDS;
        long[] tripwireSeeds = { 42L, 43L };
        int tripwireTicks = 24;
        int bucket = 16;

        long snapshots = 0, memberObs = 0, cellObs = 0;
        long negativeMembers = 0, negativeCells = 0;
        long minMember = Long.MAX_VALUE;
        double minCell = Double.MAX_VALUE;
        double worstSaturation = 0, worstHalo = 0;
        long occupied = 0, atRisk = 0;
        long snapshotsInWindow = 0, occupiedInWindow = 0, atRiskInWindow = 0;
        long[] atRiskByBucket = new long[(ticks + bucket - 1) / bucket];

        long start = System.nanoTime();
        for (long seed : seeds) {
            SubstrateFactory.Substrate substrate = phaseAHybridSubstrate(extent,
                                                                          seed,
                                                                          DEFAULT_PACKET_QUANTA,
                                                                          origin);
            boolean tripwireSeed = Arrays.stream(tripwireSeeds)
                                          .anyMatch(s -> s == seed);
            for (int tick = 0; tick < ticks; tick++) {
                if (tick > 0) {
                    substrate.run().tick(tick - 1);
                }
                boolean inWindow = tripwireSeed && tick < tripwireTicks;
                snapshots++;
                if (inWindow) {
                    snapshotsInWindow++;
                }
                int slots = substrate.field().slotCount();
                memberObs += slots;
                for (int s = 0; s < slots; s++) {
                    long q = substrate.field().quantaAt(s);
                    if (q < 0) {
                        negativeMembers++;
                    }
                    minMember = Math.min(minMember, q);
                    if (q > 0) {
                        occupied++;
                        if (inWindow) {
                            occupiedInWindow++;
                        }
                    }
                    if (q == 1) {
                        atRisk++;
                        atRiskByBucket[tick / bucket]++;
                        if (inWindow) {
                            atRiskInWindow++;
                        }
                    }
                }
                double[] cells = StructureFactor.coarseGrainedField(substrate.field());
                cellObs += cells.length;
                for (double v : cells) {
                    if (v < 0) {
                        negativeCells++;
                    }
                    minCell = Math.min(minCell, v);
                }
                worstSaturation = Math.max(worstSaturation,
                                            momentSaturation(cells, extent,
                                                              origin));
                worstHalo = Math.max(worstHalo,
                                      haloFraction(cells, signedTotal(cells)));
            }
        }

        System.out.println("Phase A non-negativity census: extent=" + extent
                            + " ticks=" + ticks + " seeds="
                            + Arrays.toString(seeds) + " packetQuanta="
                            + DEFAULT_PACKET_QUANTA);
        System.out.println("  snapshots                 = " + snapshots);
        System.out.println("  member observations       = " + memberObs);
        System.out.println("  cell observations         = " + cellObs);
        System.out.println("  negative member obs       = " + negativeMembers);
        System.out.println("  negative cell obs         = " + negativeCells);
        System.out.println("  min member / min cell     = " + minMember + " / "
                            + minCell);
        System.out.println("  worst W2 saturation       = " + worstSaturation
                            + " (tolerance "
                            + RESPONSE_MOMENT_SATURATION_TOLERANCE + ")");
        System.out.println("  worst W3 halo fraction    = " + worstHalo
                            + " (tolerance "
                            + RESPONSE_HALO_FRACTION_TOLERANCE + ")");
        System.out.println("  occupied (q>0) member obs = " + occupied);
        System.out.println("  at-risk  (q==1) member obs= " + atRisk);
        System.out.println("  IN-SUITE TRIPWIRE COVERAGE (seeds "
                            + Arrays.toString(tripwireSeeds) + " x ticks 0.."
                            + (tripwireTicks - 1) + "):");
        System.out.println("    snapshots  " + snapshotsInWindow + "/"
                            + snapshots + " = "
                            + pct(snapshotsInWindow, snapshots));
        System.out.println("    occupied   " + occupiedInWindow + "/" + occupied
                            + " = " + pct(occupiedInWindow, occupied));
        System.out.println("    AT-RISK    " + atRiskInWindow + "/" + atRisk
                            + " = " + pct(atRiskInWindow, atRisk));
        System.out.println("    at-risk by " + bucket + "-tick bucket = "
                            + Arrays.toString(atRiskByBucket)
                            + "  <- the tripwire covers the FIRST bucket only,"
                            + " i.e. the regime where negativity is LEAST likely");
        System.out.println("  elapsed s = "
                            + (System.nanoTime() - start) / 1e9);
    }

    private static String pct(long part, long whole) {
        return whole > 0
                ? String.format(Locale.ROOT, "%.4f%%", 100.0 * part / whole)
                : "n/a";
    }

    static String toTsv(Report report, String gitCommit) {
        StringBuilder sb = new StringBuilder();
        sb.append("# AnisotropyProbe Phase A measurement report\n");
        sb.append("# bead=inviscid-0nx.10\n");
        sb.append("# generator=").append(AnisotropyProbe.class.getName())
          .append('\n');
        sb.append("# gitCommit=").append(gitCommit).append('\n');
        sb.append("# extent=").append(report.extent().x).append(',')
          .append(report.extent().y).append(',').append(report.extent().z)
          .append('\n');
        sb.append("# ticks=").append(report.ticks()).append('\n');
        sb.append("# originCell=").append(report.originCell().x).append(',')
          .append(report.originCell().y).append(',')
          .append(report.originCell().z).append('\n');
        sb.append("# packetQuanta=").append(report.packetQuanta())
          .append('\n');
        StringBuilder seeds = new StringBuilder();
        for (long seed : report.seeds()) {
            if (seeds.length() > 0) {
                seeds.append(',');
            }
            seeds.append(seed);
        }
        sb.append("# seeds=").append(seeds).append('\n');
        sb.append("# transportEstimatorDefinition=abs(OLS slope of mass-weighted mean-squared-displacement-from-origin(packet seed cell) projected along d, vs tick t)\n");
        sb.append("# spectralEstimatorDefinition=abs(StructureFactor.extractRidge(StructureFactor.spectrum(fieldByTick,d)).slope()), raw unfiltered points\n");
        sb.append("# spectralZeroSlopeFraming=an all-zero spectral ridge slope is the EXPECTED signature of purely diffusive dynamics (no propagating branch, omega~i*D*k^2) - NOT an instrument malfunction, and not \"disagreement\" in a pejorative sense; TRANSPORT and SPECTRAL measure different physics (real-space spread rate vs. propagating-mode speed)\n");
        sb.append("# ratioDegenerateEpsilon=").append(RATIO_DEGENERATE_EPSILON)
          .append('\n');
        sb.append("# naivePerSeedRatioCaveat=SUMMARY rows below (mean of per-seed max/min ratios) are a DIAGNOSTIC, bounded below by 1.0 by construction, upward-biased by seed noise (order-statistic artifact, T3 critique-pattern-max-min-ratio-order-statistic-bias) - the significance statistic is POOLED_SUMMARY (seed-pooled ratio-of-means + permutation null calibration)\n");
        sb.append("# bootstrapResamples=").append(BOOTSTRAP_RESAMPLES)
          .append('\n');
        sb.append("# bootstrapRngSeed=").append(BOOTSTRAP_RNG_SEED)
          .append('\n');
        sb.append("# permutationCount=").append(PERMUTATION_COUNT)
          .append('\n');
        sb.append("# permutationRngSeed=").append(PERMUTATION_RNG_SEED)
          .append('\n');
        sb.append("# permutationDefinition=within each seed, shuffle which magnitude is labeled X100/X110/X111, pool into per-direction means across seeds, recompute ratio-of-means; empirical p-value = fraction of permuted ratios >= observed pooled ratio\n");
        sb.append("# note111=X111 probes half the k-range of X100/X110 (real FCC physics, see StructureFactor) - fewer points, higher variance\n");

        double meanEffective = 0;
        double meanTotal = 0;
        for (SeedResult sr : report.perSeed()) {
            meanEffective += sr.effectiveCollisions();
            meanTotal += sr.totalCollisions();
        }
        int nSeeds = report.perSeed().size();
        meanEffective /= nSeeds;
        meanTotal /= nSeeds;
        boolean smallN = meanEffective < SMALL_N_EFFECTIVE_COLLISIONS_THRESHOLD;
        sb.append("# smallNEarlyTimeFlag=").append(smallN).append(" (mean effective collisions/seed=")
          .append(formatPrecise(meanEffective)).append(", mean total collisions/seed=")
          .append(formatPrecise(meanTotal)).append(", threshold=")
          .append(SMALL_N_EFFECTIVE_COLLISIONS_THRESHOLD)
          .append(smallN
                  ? " - FEW real transfer events observed per seed; the OLS-fit diffusive-window assumption is NOT independently verified at this campaign scale, treat as early-time/small-N"
                  : " - collision counts comfortably above the small-N threshold")
          .append('\n');
        sb.append("# precision=%.9e\n");
        sb.append("# columns(DIRECTION rows)=recordType\tseed\testimator\tdirection\tmagnitude\tsampleSize\n");
        sb.append("# columns(COLLISIONS rows)=recordType\tseed\ttotalCollisions\teffectiveCollisions\n");
        sb.append("# columns(SUMMARY rows, DIAGNOSTIC per-seed-ratio, see naivePerSeedRatioCaveat)=recordType\testimator\tratio\tciLower\tciUpper\tnSeedsUsed\tnSeedsDegenerate\n");
        sb.append("# columns(POOLED_DIRECTION rows)=recordType\testimator\tdirection\tmean\tciLower\tciUpper\tnSeeds\n");
        sb.append("# columns(POOLED_SUMMARY rows, THE SIGNIFICANCE STATISTIC)=recordType\testimator\tpooledRatio\tpooledRatioCiLower\tpooledRatioCiUpper\tpermutationPValue\tpermutationNull95\tpermutationCount\n");
        for (SeedResult seedResult : report.perSeed()) {
            appendDirectionRows(sb, seedResult.seed(), "TRANSPORT",
                                 seedResult.transport());
            appendDirectionRows(sb, seedResult.seed(), "SPECTRAL",
                                 seedResult.spectral());
            sb.append("COLLISIONS\t").append(seedResult.seed()).append('\t')
              .append(seedResult.totalCollisions()).append('\t')
              .append(seedResult.effectiveCollisions()).append('\n');
        }
        appendSummaryRow(sb, "TRANSPORT", report.transportCi());
        appendSummaryRow(sb, "SPECTRAL", report.spectralCi());
        appendPooledDirectionRows(sb, "TRANSPORT", report.pooledTransport());
        appendPooledDirectionRows(sb, "SPECTRAL", report.pooledSpectral());
        appendPooledSummaryRow(sb, "TRANSPORT", report.pooledTransport());
        appendPooledSummaryRow(sb, "SPECTRAL", report.pooledSpectral());
        return sb.toString();
    }

    private static void appendDirectionRows(StringBuilder sb, long seed,
                                             String estimator,
                                             EstimatorResult result) {
        for (StructureFactor.Direction d : StructureFactor.Direction.values()) {
            DirectionMagnitude dm = result.perDirection().get(d);
            sb.append("DIRECTION\t").append(seed).append('\t')
              .append(estimator).append('\t').append(d).append('\t')
              .append(formatPrecise(dm.magnitude())).append('\t')
              .append(dm.sampleSize()).append('\n');
        }
    }

    private static void appendSummaryRow(StringBuilder sb, String estimator,
                                          BootstrapCi ci) {
        sb.append("SUMMARY\t").append(estimator).append('\t')
          .append(formatPrecise(ci.mean())).append('\t')
          .append(formatPrecise(ci.lower())).append('\t')
          .append(formatPrecise(ci.upper())).append('\t')
          .append(ci.nSeedsUsed()).append('\t').append(ci.nSeedsDegenerate())
          .append('\n');
    }

    private static void appendPooledDirectionRows(StringBuilder sb,
                                                    String estimator,
                                                    PooledResult pooled) {
        for (StructureFactor.Direction d : StructureFactor.Direction.values()) {
            PooledDirectionStats stats = pooled.perDirection().get(d);
            sb.append("POOLED_DIRECTION\t").append(estimator).append('\t')
              .append(d).append('\t').append(formatPrecise(stats.mean()))
              .append('\t').append(formatPrecise(stats.ciLower())).append('\t')
              .append(formatPrecise(stats.ciUpper())).append('\t')
              .append(stats.nSeeds()).append('\n');
        }
    }

    private static void appendPooledSummaryRow(StringBuilder sb,
                                                 String estimator,
                                                 PooledResult pooled) {
        sb.append("POOLED_SUMMARY\t").append(estimator).append('\t')
          .append(formatPrecise(pooled.pooledRatio().orElse(Double.NaN)))
          .append('\t').append(formatPrecise(pooled.pooledRatioCiLower()))
          .append('\t').append(formatPrecise(pooled.pooledRatioCiUpper()))
          .append('\t').append(formatPrecise(pooled.permutationPValue()))
          .append('\t').append(formatPrecise(pooled.permutationNull95()))
          .append('\t').append(pooled.permutationCount()).append('\n');
    }

    private static String formatPrecise(double v) {
        return String.format(Locale.ROOT, "%.9e", v);
    }

    /**
     * Mirrors {@code ContactAtlasGenerator#resolveGitCommit} exactly
     * (that method is private to {@code lga} and cannot be imported into
     * {@code measure} -- mirrored per the relay's instruction). Runs
     * {@code git rev-parse HEAD}, appending {@code "-dirty"} if {@code
     * git status --porcelain} reports uncommitted changes; falls back to
     * {@code "UNKNOWN"} (never throws) if {@code git} is unavailable.
     */
    static String resolveGitCommit() {
        String sha = runGit("rev-parse", "HEAD");
        if (sha == null || sha.isBlank()) {
            return "UNKNOWN";
        }
        return isDirty() ? sha + "-dirty" : sha;
    }

    private static boolean isDirty() {
        String status = runGit("status", "--porcelain");
        return status != null && !status.isBlank();
    }

    private static String runGit(String... args) {
        try {
            List<String> command = new ArrayList<>();
            command.add("git");
            command.addAll(List.of(args));
            Process process = new ProcessBuilder(command).redirectErrorStream(true)
                                                           .start();
            String output = new String(process.getInputStream().readAllBytes(),
                                        StandardCharsets.UTF_8).trim();
            int exit = process.waitFor();
            return exit == 0 ? output : null;
        } catch (IOException e) {
            return null;
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return null;
        }
    }
}
