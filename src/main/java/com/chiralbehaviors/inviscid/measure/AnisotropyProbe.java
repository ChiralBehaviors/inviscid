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

package com.chiralbehaviors.inviscid.measure;

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

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.lga.CollisionSweep;
import com.chiralbehaviors.inviscid.lga.ContactAtlasGenerator;
import com.chiralbehaviors.inviscid.lga.ContactPredicate;
import com.chiralbehaviors.inviscid.lga.ContactScan;
import com.chiralbehaviors.inviscid.lga.FccNeighborhood;
import com.chiralbehaviors.inviscid.lga.HybridAutomaton;
import com.chiralbehaviors.inviscid.lga.MemberGeometry;
import com.chiralbehaviors.inviscid.lga.QuantaExchangeRule;

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
 * lattice cell, against an all-zero background, is the natural transport
 * probe -- {@link QuantaExchangeRule} moves quanta one unit at a time
 * from a higher-quanta member to a lower one, so starting from
 * {@code {0, packetQuanta}} the field obeys a discrete maximum principle
 * (every value stays in {@code [0, packetQuanta]} for all time -- no
 * value can ever go negative or exceed the packet's own initial value,
 * since a transfer only ever moves a single unit from strictly-higher to
 * strictly-lower). The coarse-grained per-cell field ({@link
 * StructureFactor#coarseGrainedField(Necronomata)}, reused directly) is
 * therefore itself always non-negative for this specific seeding, and
 * literally IS the "quanta-deviation" from the zero background -- no
 * separate background subtraction is needed. {@link #transportEstimate}
 * still takes {@code Math.abs(...)} of every mass value defensively (the
 * design brief's own phrase is "the {@code |quanta-deviation|}
 * distribution"), so the estimator itself does not silently assume every
 * caller obeys the maximum-principle precondition (e.g. a future
 * collision rule, or a synthetic test field, might not).
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
 * <h2>The periodic-wrap subtlety, the choice made, and the exact
 * origin-relative correctness criterion (FIX 6, code-review round)</h2>
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
 * IllegalStateException}. This is now sound for ANY {@code originCell}
 * {@link #transportEstimate} is called with -- the first version of this
 * guard checked literal array-boundary coordinates ({@code 0} or
 * {@code extent-1}) regardless of {@code originCell}, which is only a
 * correct proxy for the centered case and silently under-protects an
 * off-center origin (a caller-supplied off-center origin could have mass
 * reach a genuinely wrap-invalidated cell well before touching the
 * literal array boundary, undetected by the old check).
 *
 * <p><b>Consequence, proved and documented rather than left as a silent
 * surprise:</b> for the EXACT-CENTER origin every current caller uses
 * ({@link #nearestEvenParityCenter}), the maximum possible
 * {@code |coord-origin|} over the WHOLE domain is exactly
 * {@code extentAxis/2} (attained only at the single antipodal index) --
 * i.e. a centered-origin run can mathematically never exceed this
 * criterion, for any extent, any tick count, any amount of spread. The
 * guard is therefore a real, reachable correctness backstop for the
 * general (any-origin) public API, and a PROVABLY-safe no-op for the
 * centered-origin campaign usage that motivated it -- not vacuous by
 * accident, but exact by construction. See {@code
 * AnisotropyProbeTest#transportEstimateFailsLoudlyWhenMassExceedsHalfPeriodFromOrigin}
 * for the off-center case this guard actually protects.
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
     *                    StructureFactor#coarseGrainedField(Necronomata)}'s
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
     * See class javadoc, "The periodic-wrap subtlety, the choice made,
     * and the exact origin-relative correctness criterion (FIX 6)".
     * ORIGIN-RELATIVE, mathematically exact: a cell contributes
     * "violating mass" iff its per-axis displacement from {@code origin}
     * strictly exceeds that axis's half-period ({@code extentAxis/2}) --
     * the exact point past which the naive unwrapped distance would
     * overestimate the true periodic minimum-image distance.
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

        double[][] fieldByTick = new double[ticks][];
        fieldByTick[0] = StructureFactor.coarseGrainedField(automaton);
        for (int tick = 1; tick < ticks; tick++) {
            run.tick(tick - 1);
            fieldByTick[tick] = StructureFactor.coarseGrainedField(automaton);
        }

        StructureFactor sf = new StructureFactor(extent);
        EstimatorResult transport = transportEstimate(fieldByTick, extent,
                                                        originCell);
        EstimatorResult spectral = spectralEstimate(sf, fieldByTick);
        return new SeedResult(seed, transport, spectral,
                               statistics.totalCollisions(),
                               statistics.effectiveCollisions());
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
        Point3i origin = nearestEvenParityCenter(extent);
        List<SeedResult> perSeed = new ArrayList<>(seeds.length);
        List<Double> transportRatios = new ArrayList<>();
        List<Double> spectralRatios = new ArrayList<>();
        List<Map<StructureFactor.Direction, Double>> transportMagnitudes = new ArrayList<>();
        List<Map<StructureFactor.Direction, Double>> spectralMagnitudes = new ArrayList<>();
        for (long seed : seeds) {
            SeedResult result = runOneSeed(extent, seed, ticks, packetQuanta,
                                            origin);
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

    private static void seedPacket(Necronomata automaton, Point3i originCell,
                                    int packetQuanta) {
        int base = automaton.indexOfCell(originCell);
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int m = 0; m < 30; m++) {
                frequency[base + m] = packetQuanta;
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
