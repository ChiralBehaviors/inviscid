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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.automaton.Necronomata;
import com.chiralbehaviors.inviscid.automaton.lga.CollisionSweep;
import com.chiralbehaviors.inviscid.automaton.lga.ContactAtlasGenerator;
import com.chiralbehaviors.inviscid.automaton.lga.ContactPredicate;
import com.chiralbehaviors.inviscid.automaton.lga.ContactScan;
import com.chiralbehaviors.inviscid.automaton.lga.FccNeighborhood;
import com.chiralbehaviors.inviscid.automaton.lga.HybridAutomaton;
import com.chiralbehaviors.inviscid.automaton.lga.MemberGeometry;
import com.chiralbehaviors.inviscid.automaton.lga.QuantaExchangeRule;
import com.chiralbehaviors.inviscid.automaton.measure.AnisotropyProbe.BootstrapCi;
import com.chiralbehaviors.inviscid.automaton.measure.AnisotropyProbe.DirectionMagnitude;
import com.chiralbehaviors.inviscid.automaton.measure.AnisotropyProbe.EstimatorResult;
import com.chiralbehaviors.inviscid.automaton.measure.StructureFactor.Direction;

/**
 * B.5 (bead inviscid-0nx.10): the isotropy discriminator's own tests. Tests
 * 1/2/4 exercise the data-agnostic estimator core against SYNTHETIC fields
 * (never {@code Necronomata}) per the bead's own text; test 3 exercises the
 * real automaton at the K=0 (collision-free) baseline, matching {@code
 * BaselineSpectrumHarness}'s established convention.
 *
 * @author halhildebrand
 */
public class AnisotropyProbeTest {

    private static final String RESOURCE_PATH = "lga/anisotropy-report-phaseA.tsv";

    // ------------------------------------------------------------------
    // Test 1: positive control.
    // ------------------------------------------------------------------

    /**
     * An analytically isotropic diffusion field (synthetic, NOT the
     * automaton): a separable product of per-axis Gaussian-shaped weights
     * whose variance grows at the SAME rate on x/y/z. Because the weight
     * is separable, x/y/z are exactly independent under the discrete
     * mass distribution (no residual covariance term to worry about), so
     * {@code Var([110]) = (Var(x)+Var(y))/2 = Var(x)} and {@code
     * Var([111]) = (Var(x)+Var(y)+Var(z))/3 = Var(x)} EXACTLY (up to
     * negligible tail truncation) whenever the per-axis rates are equal -
     * the anisotropy ratio A must be 1.0.
     */
    @Test
    public void syntheticIsotropicDiffusionGivesRatioOne() {
        double rate = 0.3;
        double[][] fieldByTick = gaussianPacketField(EXTENT, ORIGIN, TICKS,
                                                       SIGMA0_SQ, rate, rate,
                                                       rate, 42L, 0.0);
        EstimatorResult result = AnisotropyProbe.transportEstimate(fieldByTick,
                                                                     EXTENT,
                                                                     ORIGIN);
        assertTrue("expected a non-degenerate ratio for a genuinely spreading isotropic packet",
                   result.ratio().isPresent());
        double a = result.ratio().getAsDouble();
        assertEquals("isotropic synthetic field must report A very close to 1.0",
                     1.0, a, 0.02);

        List<Double> perSeedRatios = new ArrayList<>();
        for (long seed : SEEDS) {
            double[][] jittered = gaussianPacketField(EXTENT, ORIGIN, TICKS,
                                                        SIGMA0_SQ, rate, rate,
                                                        rate, seed, 0.05);
            perSeedRatios.add(AnisotropyProbe.transportEstimate(jittered,
                                                                  EXTENT,
                                                                  ORIGIN)
                                              .ratio().getAsDouble());
        }
        BootstrapCi ci = AnisotropyProbe.bootstrapCi(perSeedRatios,
                                                       SEEDS.length);
        // NOTE: A = max_d/min_d is mathematically bounded below by 1.0 (by
        // construction, max>=min always) - a genuinely isotropic signal
        // therefore cannot straddle 1.0 the way a two-sided statistic
        // would; the correct "within CI of 1.0" check for this one-sided-
        // bounded ratio is that the CI stays tight against that floor
        // (both bounds close to 1.0), not that it spans values below 1.0
        // (mathematically impossible here).
        assertTrue("expected the 95% CI to sit tight against the 1.0 floor for an isotropic field, was ["
                   + ci.lower() + "," + ci.upper() + "]",
                   ci.lower() >= 1.0 && ci.upper() <= 1.02);
    }

    // ------------------------------------------------------------------
    // Test 2: negative control (discriminating power).
    // ------------------------------------------------------------------

    /**
     * The SAME synthetic generator, but with the x-axis rate set to 4x
     * the y/z rate. Derivation (see class javadoc / bead handback):
     * {@code D_[100]=4R}, {@code D_[110]=(4R+R)/2=2.5R}, {@code
     * D_[111]=(4R+R+R)/3=2R} -&gt; {@code max/min = D_[100]/D_[111] =
     * 4R/2R = 2.0} EXACTLY in the continuum limit. This is the
     * discriminating-power proof: the CI must EXCLUDE 1.0.
     */
    @Test
    public void syntheticAnisotropicFieldIsDetected() {
        double rateY = 0.2;
        double rateX = 4 * rateY;
        double[][] fieldByTick = gaussianPacketField(EXTENT, ORIGIN, TICKS,
                                                       SIGMA0_SQ, rateX,
                                                       rateY, rateY, 42L, 0.0);
        EstimatorResult result = AnisotropyProbe.transportEstimate(fieldByTick,
                                                                     EXTENT,
                                                                     ORIGIN);
        assertTrue(result.ratio().isPresent());
        double a = result.ratio().getAsDouble();
        assertEquals("2x-faster-axis synthetic field must report A close to 2.0",
                     2.0, a, 0.05);
        assertEquals("max direction must be X100 (the fast axis)", Direction.X100,
                     maxDirection(result));
        assertEquals("min direction must be X111 (mixes in the two slow axes most)",
                     Direction.X111, minDirection(result));

        List<Double> perSeedRatios = new ArrayList<>();
        for (long seed : SEEDS) {
            double[][] jittered = gaussianPacketField(EXTENT, ORIGIN, TICKS,
                                                        SIGMA0_SQ, rateX,
                                                        rateY, rateY, seed,
                                                        0.05);
            perSeedRatios.add(AnisotropyProbe.transportEstimate(jittered,
                                                                  EXTENT,
                                                                  ORIGIN)
                                              .ratio().getAsDouble());
        }
        BootstrapCi ci = AnisotropyProbe.bootstrapCi(perSeedRatios,
                                                       SEEDS.length);
        assertTrue("expected the 95% CI to EXCLUDE 1.0 for a genuinely anisotropic field, was ["
                   + ci.lower() + "," + ci.upper() + "]", ci.lower() > 1.0);
    }

    // ------------------------------------------------------------------
    // Test 3: degenerate baseline, never fabricated.
    // ------------------------------------------------------------------

    /**
     * The K=0 (collision-free) baseline, real {@code Necronomata}: no
     * {@link com.chiralbehaviors.inviscid.automaton.lga.HybridAutomaton},
     * no collision resolution at all -- only {@link Necronomata#step()}
     * is called each tick, exactly matching {@code
     * BaselineSpectrumHarness}'s established K=0 convention ({@code
     * Necronomata.process(Point3i)} stays a no-op, nothing ever writes
     * {@code deltaF}). The seeded quanta packet therefore never moves:
     * every recorded field snapshot is bit-identical, so both estimators
     * MUST report DEGENERATE ({@link OptionalDouble#empty()}), never a
     * fabricated {@code 1.0} -- the direct bead inviscid-3is /
     * inviscid-0nx.7 non-vacuity criterion, restated at the anisotropy
     * layer.
     */
    @Test
    public void k0DegenerateCaseIsReportedNotFabricated() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int base = automaton.indexOfCell(origin);
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int m = 0; m < 30; m++) {
                frequency[base + m] = 5f;
            }
        });

        int t = 8;
        double[][] fieldByTick = new double[t][];
        fieldByTick[0] = StructureFactor.coarseGrainedField(automaton);
        for (int i = 1; i < t; i++) {
            automaton.step();
            fieldByTick[i] = StructureFactor.coarseGrainedField(automaton);
        }
        for (int i = 1; i < t; i++) {
            assertTrue("K=0 baseline must never change the quanta field (no collision handling ran)",
                       java.util.Arrays.equals(fieldByTick[0], fieldByTick[i]));
        }

        EstimatorResult transport = AnisotropyProbe.transportEstimate(fieldByTick,
                                                                        extent,
                                                                        origin);
        assertFalse("TRANSPORT estimator must report DEGENERATE (empty), never a fabricated ratio, on the K=0 baseline",
                    transport.ratio().isPresent());
        for (Direction d : Direction.values()) {
            assertEquals("K=0 baseline transport magnitude must be exactly 0.0 along "
                         + d, 0.0, transport.perDirection().get(d).magnitude(),
                         0.0);
        }

        StructureFactor sf = new StructureFactor(extent);
        EstimatorResult spectral = AnisotropyProbe.spectralEstimate(sf,
                                                                      fieldByTick);
        assertFalse("SPECTRAL estimator must report DEGENERATE (empty), never a fabricated ratio, on the K=0 baseline",
                    spectral.ratio().isPresent());
    }

    // ------------------------------------------------------------------
    // Test 4: bootstrap CI narrower than the claimed effect size.
    // ------------------------------------------------------------------

    /**
     * A mildly-anisotropic-but-seed-stable synthetic generator (fixed
     * rates, only a small per-seed jitter on the nuisance parameter
     * sigma0) over the 8 literal seeds 42..49: the bootstrap CI's width
     * must be narrower than the claimed effect size ({@code |2.0-1.0| =
     * 1.0}), i.e. the CI machinery correctly discriminates "stable
     * signal" from "noise swamping the estimate".
     */
    @Test
    public void estimatorsAreSeedStable() {
        double rateY = 0.2;
        double rateX = 4 * rateY;
        List<Double> perSeedRatios = new ArrayList<>();
        for (long seed : SEEDS) {
            double[][] fieldByTick = gaussianPacketField(EXTENT, ORIGIN,
                                                           TICKS, SIGMA0_SQ,
                                                           rateX, rateY,
                                                           rateY, seed, 0.05);
            EstimatorResult result = AnisotropyProbe.transportEstimate(fieldByTick,
                                                                         EXTENT,
                                                                         ORIGIN);
            assertTrue(result.ratio().isPresent());
            perSeedRatios.add(result.ratio().getAsDouble());
        }
        assertEquals("expected exactly 8 seeds", 8, SEEDS.length);

        BootstrapCi ci = AnisotropyProbe.bootstrapCi(perSeedRatios,
                                                       SEEDS.length);
        double effectSize = 1.0; // |trueRatio(2.0) - null(1.0)|
        double width = ci.upper() - ci.lower();
        assertTrue("expected CI width (" + width
                   + ") narrower than the claimed effect size (" + effectSize
                   + ") - seed-to-seed noise swamps the signal otherwise",
                   width < effectSize);
        assertEquals(0, ci.nSeedsDegenerate());
        assertEquals(8, ci.nSeedsUsed());
    }

    // ------------------------------------------------------------------
    // Bonus: the wrap-safety assertion itself (FIX 6, code-review round:
    // the guard is now ORIGIN-RELATIVE and mathematically exact -
    // |coord-origin| > extentAxis/2 on any axis is the precise point past
    // which the naive Cartesian distance overestimates the true periodic
    // minimum-image distance). A perfectly CENTERED origin can never
    // exceed this (proved in the class javadoc: max |coord-origin| over
    // the whole domain is exactly extentAxis/2 for a centered origin) -
    // so this test uses a deliberately OFF-CENTER origin, exactly the
    // case the fix targets (the public transportEstimate API accepts any
    // originCell, not just centered ones).
    // ------------------------------------------------------------------

    @Test
    public void transportEstimateFailsLoudlyWhenMassExceedsHalfPeriodFromOrigin() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = new Point3i(1, 4, 4); // deliberately off-center on x
        int cellCount = extent.x * extent.y * extent.z;
        double[][] fieldByTick = new double[2][cellCount];
        fieldByTick[0][(origin.x * extent.y + origin.y) * extent.z
                        + origin.z] = 10.0;
        // Tick 1: place mass at x=6 - |6-1|=5 > extent.x/2=4, a genuine
        // half-period violation relative to this off-center origin, well
        // short of the literal array boundary (x=7) the OLD,
        // origin-independent check would have required.
        fieldByTick[1][(6 * extent.y + 4) * extent.z + 4] = 10.0;

        try {
            AnisotropyProbe.transportEstimate(fieldByTick, extent, origin);
            fail("expected IllegalStateException: half-period-exceeding mass invalidates the Cartesian moment");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().toLowerCase().contains("wrap"));
        }
    }

    /**
     * Companion regression, proving the class javadoc's claim: a
     * perfectly CENTERED origin can place mass at the literal array
     * boundary (the old check's trigger point) WITHOUT tripping the new,
     * mathematically-exact guard, because that boundary is exactly the
     * half-period tie point (naive distance == true periodic distance
     * there, not a violation).
     */
    @Test
    public void transportEstimateAcceptsBoundaryMassForACenteredOrigin() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int cellCount = extent.x * extent.y * extent.z;
        double[][] fieldByTick = new double[2][cellCount];
        fieldByTick[0][(origin.x * extent.y + origin.y) * extent.z
                        + origin.z] = 10.0;
        fieldByTick[1][0] = 10.0; // literal (0,0,0) - the old check's trigger, now a tie, not a violation

        EstimatorResult result = AnisotropyProbe.transportEstimate(fieldByTick,
                                                                     extent,
                                                                     origin);
        assertTrue("centered-origin boundary mass at the exact half-period tie point must NOT trip the guard",
                   result != null);
    }

    // ------------------------------------------------------------------
    // FIX 5 (stacked review): test the null-calibration machinery itself
    // with a MATCHED-NOISE isotropic control (sparse-discrete-jump-like
    // CV comparable to the real campaign, not the earlier tests' smooth
    // low-jitter Gaussian) and confirm the 2x-anisotropic case is still
    // flagged. These two together make the null calibration falsifiable.
    // ------------------------------------------------------------------

    /**
     * Per-seed-per-direction magnitudes drawn directly (bypassing field
     * generation entirely) from a noisy, ISOTROPIC generator at ~22% CV
     * -- comparable to the real Phase A campaign's measured per-direction
     * CV (21/25/18%), unlike the earlier smooth-Gaussian-field positive
     * control's negligible jitter. The permutation test MUST NOT flag
     * this as significant (p &gt; 0.05) -- this is the machinery's own
     * false-positive-rate check at realistic noise.
     */
    @Test
    public void permutationTestIsNonSignificantForMatchedNoiseIsotropicControl() {
        double mean = 1.5e-4;
        double cv = 0.22;
        List<Map<Direction, Double>> perSeedMagnitudes = new ArrayList<>();
        for (long seed : SEEDS) {
            Random seedRandom = new Random(seed * 31 + 7);
            Map<Direction, Double> magnitudes = new EnumMap<>(Direction.class);
            for (Direction d : Direction.values()) {
                magnitudes.put(d, noisyMagnitude(seedRandom, mean, cv));
            }
            perSeedMagnitudes.add(magnitudes);
        }

        AnisotropyProbe.PooledResult pooled = AnisotropyProbe.pooledEstimate(perSeedMagnitudes);
        assertTrue("expected a non-degenerate pooled ratio for genuinely noisy magnitudes",
                   pooled.pooledRatio().isPresent());
        assertTrue("expected the permutation p-value (" + pooled.permutationPValue()
                   + ") to be NON-significant (>0.05) for a matched-noise isotropic control",
                   pooled.permutationPValue() > 0.05);
    }

    /**
     * The companion negative control: the SAME matched-noise generator,
     * but with the X100 mean set to 2x the X110/X111 mean -- the
     * permutation test MUST flag this as significant (p &lt; 0.05),
     * proving the machinery has discriminating power at this noise level,
     * not just a tendency to always report "not significant".
     */
    @Test
    public void permutationTestDetectsGenuineTwoFoldAnisotropyAtMatchedNoise() {
        double meanMinor = 1.0e-4;
        double meanMajor = 2 * meanMinor;
        double cv = 0.22;
        List<Map<Direction, Double>> perSeedMagnitudes = new ArrayList<>();
        for (long seed : SEEDS) {
            Random seedRandom = new Random(seed * 31 + 7);
            Map<Direction, Double> magnitudes = new EnumMap<>(Direction.class);
            magnitudes.put(Direction.X100, noisyMagnitude(seedRandom, meanMajor, cv));
            magnitudes.put(Direction.X110, noisyMagnitude(seedRandom, meanMinor, cv));
            magnitudes.put(Direction.X111, noisyMagnitude(seedRandom, meanMinor, cv));
            perSeedMagnitudes.add(magnitudes);
        }

        AnisotropyProbe.PooledResult pooled = AnisotropyProbe.pooledEstimate(perSeedMagnitudes);
        assertTrue(pooled.pooledRatio().isPresent());
        assertEquals("pooled ratio-of-means should recover close to the true 2.0x signal",
                     2.0, pooled.pooledRatio().getAsDouble(), 0.3);
        assertTrue("expected the permutation p-value (" + pooled.permutationPValue()
                   + ") to be significant (<0.05) for a genuine 2x anisotropic signal",
                   pooled.permutationPValue() < 0.05);
    }

    /** A non-negative, mean/cv-parameterized noise draw (reflected Gaussian) - see the two tests above. */
    private static double noisyMagnitude(Random random, double mean, double cv) {
        return Math.max(1e-12, mean * (1 + cv * random.nextGaussian()));
    }

    // ------------------------------------------------------------------
    // inviscid-0sn: the +1 continuity correction on the permutation
    // p-value. Under H0 the observed statistic is exchangeable with the
    // null draws, so it must be counted as one of its own reference set:
    // p = (countGe+1)/(permutations+1), never the uncorrected
    // countGe/permutations (which can report a fabricated exact-zero
    // p-value that a finite permutation count can never actually prove).
    // This test does NOT independently reimplement null-distribution
    // generation: it calls {@link AnisotropyProbe#permutationNullDistribution}
    // directly, the SAME shuffle-based helper {@link AnisotropyProbe#pooledEstimate}
    // uses internally (same fixture, same PERMUTATION_COUNT, same
    // PERMUTATION_RNG_SEED -&gt; bit-identical null array). What IS
    // independently recomputed is only the layer on top of that array --
    // the countGe tally and the +1-correction formula -- and THAT
    // independent computation is what gets pinned against the SUT's
    // reported p-value, not against the old (uncorrected) formula. See
    // {@code countGeIncludesExactTieWithObserved} below for the
    // dedicated, forced-tie test of the countGe boundary condition this
    // test's continuous-noise fixture cannot exercise (ties have ~0
    // probability here).
    // ------------------------------------------------------------------

    @Test
    public void permutationPValueUsesPlusOneContinuityCorrection() {
        double mean = 1.5e-4;
        double cv = 0.22;
        List<Map<Direction, Double>> perSeedMagnitudes = new ArrayList<>();
        for (long seed : SEEDS) {
            Random seedRandom = new Random(seed * 31 + 7);
            Map<Direction, Double> magnitudes = new EnumMap<>(Direction.class);
            for (Direction d : Direction.values()) {
                magnitudes.put(d, noisyMagnitude(seedRandom, mean, cv));
            }
            perSeedMagnitudes.add(magnitudes);
        }

        AnisotropyProbe.PooledResult pooled = AnisotropyProbe.pooledEstimate(perSeedMagnitudes);
        assertTrue("expected a non-degenerate pooled ratio for this fixture",
                   pooled.pooledRatio().isPresent());
        double observed = pooled.pooledRatio().getAsDouble();

        double[] nulls = AnisotropyProbe.permutationNullDistribution(perSeedMagnitudes,
                                                                       AnisotropyProbe.PERMUTATION_COUNT,
                                                                       AnisotropyProbe.PERMUTATION_RNG_SEED);
        assertEquals("null draw count must match the SUT's reported permutationCount",
                     nulls.length, pooled.permutationCount());
        long countGe = AnisotropyProbe.countGe(nulls, observed);
        double expectedP = (countGe + 1) / (double) (nulls.length + 1);
        assertEquals("permutation p-value must use the +1 continuity correction: (countGe+1)/(permutationCount+1)",
                     expectedP, pooled.permutationPValue(), 1e-12);
    }

    // ------------------------------------------------------------------
    // inviscid-0sn / stacked-review follow-up: the countGe tie-handling
    // boundary. Real fixtures use continuous synthetic noise, so an
    // exact tie between a null draw and the observed statistic has ~0
    // probability of ever occurring -- mutation testing confirmed
    // flipping the countGe loop's &gt;= to &gt; left every existing test
    // green. This forces the tie directly against the extracted
    // {@link AnisotropyProbe#countGe} method (bypassing the shuffle
    // machinery entirely) to close that gap.
    // ------------------------------------------------------------------

    @Test
    public void countGeIncludesExactTieWithObserved() {
        double[] nulls = { 1.0, 2.0, 3.0 };
        double observed = 2.0;
        // ties count TOWARD the null (conservative standard practice,
        // consistent with the +1 correction's exchangeability
        // rationale): 2.0 (the exact tie) and 3.0 both count, so 2 - not
        // 1, which is what an (incorrect) strict "&gt;" would report.
        assertEquals("an exact tie between a null draw and the observed statistic must count toward countGe",
                     2L, AnisotropyProbe.countGe(nulls, observed));
    }

    // ==================================================================
    // Bead inviscid-0nx.28 (E.1) / inviscid-o24: the signed-background
    // re-derivation. See AnisotropyProbe's class javadoc, "MATCHED-PAIR
    // transport under a signed background" and "The periodic-wrap
    // precondition: two limbs".
    // ==================================================================

    /**
     * <b>The corrected red test.</b> Bead inviscid-0nx.28 was written on
     * the premise that a uniform SIGNED background makes the existing
     * {@code assertWrapSafe} throw at tick 0 on the default 8^3 campaign
     * geometry. It does not, and cannot: {@code FccNeighborhood} forces
     * every extent axis EVEN, and {@code nearestEvenParityCenter} then
     * places the origin at {@code L/2} (or, on z, {@code L/2 - 1}), for
     * both of which {@code max_i |i-origin| = L/2} exactly -- while the
     * criterion is a STRICT {@code >}. So limb W1 is not "wrong under
     * signed background", it is STRUCTURALLY UNREACHABLE on every legal
     * campaign geometry, for every field. That is the real defect (a
     * guard that certifies nothing where it is needed), and it is what
     * motivates limb W2.
     *
     * <p><b>This is a re-EVALUATION, not a discovery.</b> The
     * unreachability argument was already derived verbatim in {@code
     * AnisotropyProbe}'s javadoc before this bead, framed there as a
     * virtue ("exact by construction") and reviewed in that framing. What
     * changed is the verdict; what is new is this enumeration.
     *
     * <p>Exhaustive rather than illustrative: EVERY legal extent with each
     * axis even in {@code 4..16} (343 combinations, NON-CUBIC included),
     * with all three axes asserted -- an earlier version covered only
     * cubic {@code L} and only the x/z axes, which would not have caught
     * an axis-asymmetric regression in {@code nearestEvenParityCenter}.
     *
     * <p>Non-vacuity, two ways: the fixture is asserted to genuinely have
     * full support with BOTH signs present (so the guard cannot be passing
     * by skipping zero cells), and the SAME fixture is asserted to make
     * the guard THROW for an off-center origin (so the guard is not
     * globally dead code -- it is unreachable on campaign geometry
     * specifically).
     */
    @Test
    public void wrapSafetyLimbW1IsStructurallyUnreachableOnEveryCampaignGeometry() {
        int[] axes = { 4, 6, 8, 10, 12, 14, 16 };
        int combinations = 0;
        for (int lx : axes) {
            for (int ly : axes) {
                for (int lz : axes) {
                    combinations++;
                    Point3i extent = new Point3i(lx, ly, lz);
                    Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
                    double[] uniformSigned = uniformSignedField(extent);

                    int positives = 0;
                    int negatives = 0;
                    for (double v : uniformSigned) {
                        assertTrue("fixture must have FULL support - a zero cell would let the guard skip it",
                                   v != 0.0);
                        if (v > 0) {
                            positives++;
                        } else {
                            negatives++;
                        }
                    }
                    assertTrue("fixture must carry BOTH signs at " + extent,
                               positives > 0 && negatives > 0);

                    // The geometric reason, asserted directly on ALL THREE
                    // axes rather than implied.
                    assertEquals("max |i-origin| must be exactly half the period on x at "
                                 + extent, lx / 2,
                                 Math.max(origin.x, lx - 1 - origin.x));
                    assertEquals("...on y at " + extent, ly / 2,
                                 Math.max(origin.y, ly - 1 - origin.y));
                    assertEquals("...and on z, where the parity decrement lands, at "
                                 + extent, lz / 2,
                                 Math.max(origin.z, lz - 1 - origin.z));

                    // Therefore: no throw, on a fully-supported signed field.
                    AnisotropyProbe.assertWrapSafe(uniformSigned, extent,
                                                    origin, 0);
                }
            }
        }
        assertEquals("the enumeration must actually cover every even-axis combination in 4..16",
                     343, combinations);

        // ...and the guard is NOT globally dead: the same fixture trips it
        // as soon as the origin is off-center.
        Point3i extent = new Point3i(8, 8, 8);
        try {
            AnisotropyProbe.assertWrapSafe(uniformSignedField(extent), extent,
                                            new Point3i(1, 1, 1), 0);
            fail("expected IllegalStateException: with origin (1,1,1) the maximum |i-origin| is 6 > 8/2, so limb W1 IS reachable off-center");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("wrap-safety bound violated"));
        }
    }

    /**
     * Limb W2 fires exactly where limb W1 is blind, proving the two are
     * not redundant. Same 8^3 centered-origin campaign geometry as above:
     * a difference field with a quarter of its L1 mass parked on the
     * antipodal shell passes W1 (nothing exceeds the half-period, so the
     * naive Cartesian distance is still exact) yet has driven the second
     * moment well up its relaxation curve, which is what W2 refuses.
     *
     * <p>The saturation ratio is asserted to an EXACT hand-computed value,
     * not merely to "&gt; tolerance": {@code M_[100] = (30*0 + 10*4^2)/40 =
     * 4.0} and {@code M^unif_[100] = (1/8)*sum_{i=0..7}(i-4)^2 = 5.5}, so
     * {@code sat = 4/5.5 = 8/11}. A hardcoded return value cannot pass
     * this.
     */
    @Test
    public void saturationLimbW2FiresWhereExactnessLimbW1IsBlind() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        double[] difference = new double[extent.x * extent.y * extent.z];
        difference[index(extent, origin.x, origin.y, origin.z)] = 30.0;
        difference[index(extent, 0, 4, 4)] = 10.0; // |0-4| == 4 == half: shell, not a W1 violation

        // W1 alone: silent.
        AnisotropyProbe.assertWrapSafe(difference, extent, origin, 7);

        assertEquals("the uniform-field [100] second moment about a centered origin on L=8 is (1/8)*sum(i-4)^2 = 5.5",
                     5.5,
                     AnisotropyProbe.uniformSecondMoment(extent, origin,
                                                          Direction.X100),
                     1e-12);
        assertEquals("moment saturation must be M_[100]/M^unif_[100] = ((30*0 + 10*16)/40) / 5.5 = 8/11",
                     8.0 / 11.0,
                     AnisotropyProbe.momentSaturation(difference, extent,
                                                       origin),
                     1e-12);

        try {
            AnisotropyProbe.assertResponseLocalized(difference, extent, origin,
                                                     7);
            fail("expected IllegalStateException: the response's second moment has reached 8/11 of its fully-delocalized value, far above the "
                 + AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE
                 + " tolerance");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("wrap-saturation"));
        }
    }

    /**
     * <b>The regression that retired limb W2's first measure.</b> W2
     * originally measured the fraction of {@code ||D||_1} sitting on the
     * half-period SHELL ({@code |coord-origin| >= extentAxis/2} on some
     * axis). At {@code L=8} with the origin at {@code 4} that set is
     * {@code {i=0} u {j=0} u {k=0}} -- a three-plane skin. A difference
     * field spread UNIFORMLY over the interior {@code (L-1)^3}, with zero
     * on exactly those three planes, therefore measured {@code 0.0} and
     * passed silently while being maximally delocalized: the same species
     * of defect as limb W1's proven vacuity, a precondition that cannot
     * fire on the case it exists for.
     *
     * <p>Under the moment saturation ratio the same field measures
     * {@code M_[100]/M^unif_[100] = ((1/7)*sum_{i=1..7}(i-4)^2)/5.5 =
     * 4.0/5.5 = 8/11} and fires. Both the old measure's blindness (via the
     * shell mass, recomputed here so the claim is checked rather than
     * recounted) and the new measure's detection are asserted.
     */
    @Test
    public void saturationLimbW2FiresOnTheInteriorUniformFieldTheShellMeasureMissed() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        double[] difference = new double[extent.x * extent.y * extent.z];
        for (int i = 1; i < extent.x; i++) {
            for (int j = 1; j < extent.y; j++) {
                for (int k = 1; k < extent.z; k++) {
                    difference[index(extent, i, j, k)] = 1.0;
                }
            }
        }

        // The retired measure, recomputed locally: exactly zero.
        assertEquals("the retired shell measure read exactly 0.0 on this field - that blindness is why it was replaced",
                     0.0, legacyShellMassFraction(difference, extent, origin),
                     0.0);
        // W1 is likewise silent (nothing STRICTLY exceeds the half-period).
        AnisotropyProbe.assertWrapSafe(difference, extent, origin, 0);

        assertEquals("the replacement measure must see this field as 8/11 saturated",
                     8.0 / 11.0,
                     AnisotropyProbe.momentSaturation(difference, extent,
                                                       origin),
                     1e-12);
        try {
            AnisotropyProbe.assertResponseLocalized(difference, extent, origin,
                                                     0);
            fail("expected IllegalStateException: a field uniform over the interior (L-1)^3 is maximally delocalized and must not pass W2");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("wrap-saturation"));
        }
    }

    /**
     * <b>W2's blind spot, MEASURED and pinned as a known limitation rather
     * than left to be rediscovered</b> (bead inviscid-0nx.28, round-2
     * review). The guard is weakest precisely in the anisotropic regime the
     * instrument exists to detect.
     *
     * <p>Mechanism, checked here by arithmetic rather than asserted:
     * {@link Direction}'s set is
     * {@code {dx, (dx+dy)/sqrt(2), (dx+dy+dz)/sqrt(3)}}, so no probe
     * isolates {@code y} or {@code z}. A {@code z}-only displacement
     * reaches {@code X111} alone and is divided by 3 there, while
     * {@code M^unif_X111} is the 3-D uniform moment with all three axes
     * contributing -- so a {@code z}-confined response's saturation is
     * capped far below 1, and the tolerance sits near the TOP of its
     * attainable range instead of at a quarter of it.
     *
     * <p>Four measurements, all on {@code 8^3}:
     * <ol>
     * <li>the CEILING: a response uniform along the origin's {@code z}
     * column reaches {@code sat = 5.5/(3*6.0) = 0.3056}, so {@code 0.25}
     * is {@code 81.8%} of everything such a response can ever reach;</li>
     * <li>STATIC blindness: the interior {@code z} line reads
     * {@code 0.2222} and PASSES, and diluting it with a compact isotropic
     * blob at the origin only LOWERS that;</li>
     * <li>AXIS-SPECIFICITY, which is what makes this a probe-set defect
     * rather than a tolerance being loose: the {@code y}-pencil equivalent
     * reads {@code 0.3478} and FIRES;</li>
     * <li>DYNAMIC blindness, the one that matters for the estimator: at
     * MATCHED saturation a {@code 20:1} {@code z}-preferential wrapped
     * Gaussian understates the true {@code X111} rate by two orders of
     * magnitude more than the isotropic field the tolerance was calibrated
     * on -- with W2 silent for both.</li>
     * </ol>
     * The comparison in (4) is at matched saturation deliberately: comparing
     * at matched TIME would confound the probe-set defect with the two
     * fields simply having spread by different amounts.
     */
    @Test
    public void saturationLimbW2IsBlindToAResponseConfinedToTheLeastProbedAxis() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        double tolerance = AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE;
        int cells = extent.x * extent.y * extent.z;

        // The uniform references the ratio is taken against.
        assertEquals("M^unif[100]", 5.5,
                     AnisotropyProbe.uniformSecondMoment(extent, origin,
                                                          Direction.X100),
                     1e-12);
        assertEquals("M^unif[110]", 5.75,
                     AnisotropyProbe.uniformSecondMoment(extent, origin,
                                                          Direction.X110),
                     1e-12);
        assertEquals("M^unif[111] is the 3-D uniform moment, with all three axes contributing - which is exactly why a z-confined response cannot approach it",
                     6.0,
                     AnisotropyProbe.uniformSecondMoment(extent, origin,
                                                          Direction.X111),
                     1e-9);

        // (1) THE CEILING. A z-confined response's whole attainable range.
        double[] zColumn = new double[cells];
        for (int k = 0; k < extent.z; k++) {
            zColumn[index(extent, origin.x, origin.y, k)] = 1.0;
        }
        double ceiling = AnisotropyProbe.momentSaturation(zColumn, extent,
                                                           origin);
        assertEquals("a response uniform along the origin's z column saturates at 5.5/(3*6.0) - E[dz^2] is the uniform 5.5, but X111 divides it by 3",
                     5.5 / (3 * 6.0), ceiling, 1e-9);
        assertTrue("the tolerance must sit above 80% of a z-confined response's ENTIRE attainable range - that is the defect, stated as a number: "
                   + (tolerance / ceiling), tolerance / ceiling > 0.80);

        // (2) STATIC blindness: maximally delocalized ALONG z, yet silent.
        double[] zLine = new double[cells];
        for (int k = 1; k < extent.z; k++) {
            zLine[index(extent, origin.x, origin.y, k)] = 1.0;
        }
        double zLineSat = AnisotropyProbe.momentSaturation(zLine, extent,
                                                            origin);
        assertEquals("the interior z line reads (4/3)/6.0", 4.0 / 3.0 / 6.0,
                     zLineSat, 1e-9);
        assertTrue("...and therefore PASSES W2, which is the blindness: "
                   + zLineSat, zLineSat < tolerance);
        AnisotropyProbe.assertResponseLocalized(zLine, extent, origin, 0);
        double previous = zLineSat;
        for (double blob : new double[] { 0.02, 0.05, 0.10 }) {
            double[] diluted = zLine.clone();
            diluted[index(extent, origin.x, origin.y,
                           origin.z)] += blob * (extent.z - 1);
            double sat = AnisotropyProbe.momentSaturation(diluted, extent,
                                                           origin);
            assertTrue("adding a compact isotropic blob at the origin can only LOWER the measured saturation, so it cannot rescue the guard (blob "
                       + blob + " -> " + sat + ")", sat < previous);
            AnisotropyProbe.assertResponseLocalized(diluted, extent, origin, 0);
            previous = sat;
        }

        // (3) AXIS-SPECIFICITY. Same construction on y, which X110 probes.
        double[] yLine = new double[cells];
        for (int j = 1; j < extent.y; j++) {
            yLine[index(extent, origin.x, j, origin.z)] = 1.0;
        }
        double yLineSat = AnisotropyProbe.momentSaturation(yLine, extent,
                                                            origin);
        assertEquals("the y pencil reads 2.0/5.75 through X110 - y is divided by 2, z by 3, which is the whole difference",
                     2.0 / 5.75, yLineSat, 1e-9);
        assertTrue("the y pencil FIRES where the z pencil is silent: the blindness is axis-specific, i.e. a defect of the probe SET",
                   yLineSat > tolerance);

        // (4) DYNAMIC blindness at MATCHED saturation.
        for (double target : new double[] { 0.20, 0.24 }) {
            double isotropic = x111RateUnderstatementAtSaturation(extent, origin,
                                                                    1, 1, 1,
                                                                    target);
            double zPreferential = x111RateUnderstatementAtSaturation(extent,
                                                                        origin,
                                                                        1, 1, 20,
                                                                        target);
            assertTrue("W2 must be SILENT at the matched saturation " + target
                       + ", else the comparison is not about a blind spot",
                       target < tolerance);
            assertTrue("the calibration must be honoured for the ISOTROPIC family it was derived on: understatement at sat="
                       + target + " was " + isotropic,
                       Math.abs(isotropic) < 0.02);
            assertTrue("...and MISSED for a 20:1 z-preferential response at the SAME saturation: understatement was "
                       + zPreferential + " at sat=" + target,
                       zPreferential > 0.15);
            assertTrue("the miss must be an order of magnitude or more, else this is tolerance slack rather than a probe-set defect (iso="
                       + isotropic + ", z-pref=" + zPreferential + ")",
                       zPreferential > 10 * Math.abs(isotropic));
        }
    }

    /**
     * <b>The MODEL-FREE reading of the saturation ratio</b> (bead
     * inviscid-0nx.28, round-2 review). The wrapped-Gaussian calibration
     * below bounds SLOPE bias only inside its own one-parameter family;
     * what {@code sat} bounds for an ARBITRARY response is a MASS fraction,
     * by Markov's inequality applied to the probability measure
     * {@code mu = |D|/||D||_1} whose second moment {@code M_d} is:
     * <pre>
     *   mu{ cell : proj_d^2 &gt;= P } &lt;= M_d / P     for every P &gt; 0
     * </pre>
     * Checked here at both thresholds the javadoc quotes -- {@code P =
     * M_d^unif} (giving "at most {@code sat_d} of the mass at or beyond
     * uniform-typical displacement") and {@code P} at W1's own ceiling
     * (giving the sharper {@code 3.125%} at {@code L=8}) -- on random
     * signed fields, and with a NON-VACUITY limb: a construction that makes
     * the {@code P = M^unif} bound nearly tight, so the test could not pass
     * merely because the inequality is slack on everything.
     */
    @Test
    public void momentSaturationBoundsDelocalizedMassDistributionFree() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int cells = extent.x * extent.y * extent.z;
        Random random = new Random(20260809L);

        for (int trial = 0; trial < 100; trial++) {
            double[] field = new double[cells];
            for (int i = 0; i < cells; i++) {
                field[i] = random.nextInt(7) - 3;
            }
            if (AnisotropyProbe.responseL1(field) == 0) {
                continue;
            }
            for (Direction d : Direction.values()) {
                double uniform = AnisotropyProbe.uniformSecondMoment(extent,
                                                                      origin, d);
                double moment = uniform
                                 * saturationAlong(field, extent, origin, d);
                assertTrue("Markov at P = M^unif along " + d
                           + ": the mass fraction at or beyond uniform-typical squared displacement must not exceed sat",
                           massFractionBeyond(field, extent, origin, d, uniform)
                           <= moment / uniform + 1e-12);
            }
        }

        // The SHARPER threshold the javadoc quotes: W1 admits |d_a| <= L/2,
        // so proj_X111^2 <= (3*4)^2/3 = 48 on 8^3. At sat = 0.25 that is
        // M_X111 <= 1.5, i.e. at most 1.5/48 = 3.125% of the mass out there.
        double pMax = (3.0 * (extent.x / 2)) * (3.0 * (extent.x / 2)) / 3.0;
        assertEquals("W1's ceiling on proj_X111^2 at L=8", 48.0, pMax, 1e-12);
        double atTolerance = AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE
                              * AnisotropyProbe.uniformSecondMoment(extent,
                                                                     origin,
                                                                     Direction.X111);
        assertEquals("the model-free ceiling reading at the W2 tolerance",
                     0.03125, atTolerance / pMax, 1e-9);

        // NON-VACUITY: a field whose mass sits exactly at uniform-typical
        // X100 displacement makes the P = M^unif bound nearly tight, so the
        // loop above is not merely exploiting a slack inequality.
        double[] tight = new double[cells];
        for (int j = 0; j < extent.y; j++) {
            for (int k = 0; k < extent.z; k++) {
                // |dx| = 2 gives proj_X100^2 = 4; the nearest lattice
                // displacement at or above M^unif_X100 = 5.5 is |dx| = 3.
                tight[index(extent, origin.x - 3, j, k)] = 1.0;
            }
        }
        double tightSat = saturationAlong(tight, extent, origin,
                                           Direction.X100);
        double tightFraction = massFractionBeyond(tight, extent, origin,
                                                   Direction.X100,
                                                   AnisotropyProbe.uniformSecondMoment(extent,
                                                                                        origin,
                                                                                        Direction.X100));
        assertEquals("all of this field's mass sits at proj^2 = 9 >= M^unif = 5.5",
                     1.0, tightFraction, 1e-12);
        assertEquals("...and its saturation is 9/5.5, so the Markov bound reads 1.636 against an actual 1.0 - within a factor of 2, i.e. genuinely binding",
                     9.0 / 5.5, tightSat, 1e-9);
        assertTrue("the bound must be respected on the near-tight construction too",
                   tightFraction <= tightSat + 1e-12);
    }

    /**
     * <b>The tolerance is re-derived here, not asserted in prose.</b> The
     * predecessor tolerance (1% of L1 mass on the shell) was calibrated
     * against distinguishability from full delocalization -- the wrong
     * reference quantity, since what the estimator actually suffers is
     * OLS SLOPE BIAS. This test performs the calibration {@code
     * AnisotropyProbe}'s javadoc quotes: place a wrapped
     * (periodic-image-summed) Gaussian of variance {@code s} on {@code L}
     * sites about a centered origin -- the exact 1-D problem the {@code
     * X100} projection solves -- and compare {@code dM/ds} against the
     * free-diffusion value of {@code 1}.
     *
     * <p>Asserts (a) that at the chosen tolerance the instantaneous rate
     * is understated by less than 2% for every {@code L} in {@code 6..40},
     * and (b) NON-VACUITY: that the bias is not simply negligible
     * everywhere -- at twice the tolerance it is materially larger. Without
     * (b) the test would pass for a tolerance of {@code 0.001} or
     * {@code 0.9} alike.
     */
    @Test
    public void momentSaturationToleranceBoundsTheSlopeUnderstatement() {
        double tolerance = AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE;
        for (int l : new int[] { 6, 8, 16, 40 }) {
            double atTolerance = wrappedGaussianSlopeUnderstatement(l,
                                                                     tolerance);
            assertTrue("at sat=" + tolerance + " the instantaneous-rate understatement on L="
                       + l + " must stay below 2%, was " + atTolerance,
                       atTolerance < 0.02);

            double atDouble = wrappedGaussianSlopeUnderstatement(l,
                                                                  2 * tolerance);
            assertTrue("non-vacuity: at twice the tolerance the understatement on L="
                       + l
                       + " must be materially worse (else the tolerance is not binding anything), was "
                       + atDouble, atDouble > 5 * atTolerance);
        }
    }

    /**
     * The local {@link #secondMomentOf} replica exists only because
     * {@code AnisotropyProbe.secondMoment} is private. A replica can drift
     * from the thing it replicates, so it is validated against production
     * output rather than trusted: {@link
     * AnisotropyProbe#uniformSecondMoment} IS production's
     * {@code secondMoment} run over an all-ones field, so agreement on that
     * field in every direction pins the projection convention and the
     * normalization together.
     */
    @Test
    public void secondMomentReplicaAgreesWithProduction() {
        for (Point3i extent : new Point3i[] { new Point3i(6, 6, 6),
                                               new Point3i(8, 8, 8),
                                               new Point3i(4, 6, 8) }) {
            Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
            double[] ones = new double[extent.x * extent.y * extent.z];
            Arrays.fill(ones, 1.0);
            for (Direction d : Direction.values()) {
                assertEquals("the test-local second-moment replica must reproduce production exactly on "
                             + extent + " along " + d,
                             AnisotropyProbe.uniformSecondMoment(extent, origin,
                                                                  d),
                             secondMomentOf(ones, extent, origin, d), 1e-12);
            }
        }
    }

    /**
     * {@link AnisotropyProbe#uniformSecondMoment} is memoized (it depends
     * only on {@code (extent, origin, direction)}, never on the field, and
     * {@code momentSaturation} calls it once per direction per snapshot).
     * A cache keyed on too little is a silent wrong-answer bug, so each of
     * the three key components is varied INDEPENDENTLY here and the cached
     * value compared against a fresh, cache-free recomputation.
     *
     * <p>The non-obvious component is {@code origin}: two different origins
     * on the same extent generally give the same value on a symmetric axis
     * and different values otherwise, so a cache keyed on extent alone
     * would pass a naive same-origin test.
     */
    @Test
    public void uniformSecondMomentMemoizationIsKeyedOnItsFullDependencySet() {
        Point3i[] extents = { new Point3i(6, 6, 6), new Point3i(8, 8, 8),
                              new Point3i(4, 6, 8) };
        Point3i[] origins = { new Point3i(2, 2, 2), new Point3i(3, 3, 3),
                              new Point3i(1, 2, 3) };
        boolean sawDifference = false;
        double first = Double.NaN;
        for (Point3i extent : extents) {
            for (Point3i origin : origins) {
                if (origin.x >= extent.x || origin.y >= extent.y
                    || origin.z >= extent.z) {
                    continue;
                }
                double[] ones = new double[extent.x * extent.y * extent.z];
                Arrays.fill(ones, 1.0);
                for (Direction d : Direction.values()) {
                    double cached = AnisotropyProbe.uniformSecondMoment(extent,
                                                                         origin,
                                                                         d);
                    // ...and again, so a first-call-only bug would show.
                    assertEquals("repeated calls must agree", cached,
                                 AnisotropyProbe.uniformSecondMoment(extent,
                                                                      origin, d),
                                 0.0);
                    assertEquals("the memoized value must equal a fresh all-ones sweep for extent "
                                 + extent + " origin " + origin + " direction "
                                 + d, secondMomentOf(ones, extent, origin, d),
                                 cached, 1e-12);
                    if (Double.isNaN(first)) {
                        first = cached;
                    } else if (Math.abs(cached - first) > 1e-9) {
                        sawDifference = true;
                    }
                }
            }
        }
        assertTrue("NON-VACUITY: the probed key space must actually produce differing values, else a cache keyed on nothing at all would pass",
                   sawDifference);
    }

    /**
     * <b>Limb W3's anchor and what it buys</b> (bead inviscid-0nx.28,
     * round-2 review). The round-2 critique's finding was that W2's
     * tolerance bounds a LEVEL while the estimator's error is driven by a
     * RATE, and that bounding {@code f_halo}'s level does not bound
     * {@code df/dt} either. What rescues the level bound is the ANCHOR:
     * {@code transportEstimateMatchedPair} asserts
     * {@code ||D(.,0)||_1 == |S|}, i.e. {@code f_halo(0) = 0} exactly, so a
     * level bound IS an excursion bound.
     *
     * <p>Both halves are checked by RECOMPUTATION on a constructed
     * core-plus-halo field, not by restating the algebra:
     * <ul>
     * <li>the decomposition identity
     * {@code M(t)-M(0) = [m_core(t)-m_core(0)] + f(t)*[m_halo(t)-m_core(t)]}
     * holds exactly, so the halo's contribution to the measured CHANGE is
     * exactly that second term;</li>
     * <li>that term's magnitude is at most
     * {@code f_halo * max_cell proj_d^2}, which W1 bounds -- so the
     * tolerance does bound the halo's total moment excursion over the
     * window even though it bounds no rate.</li>
     * </ul>
     */
    @Test
    public void haloFractionIsAnchoredAtZeroAndBoundsTheHaloMomentExcursion() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int cells = extent.x * extent.y * extent.z;
        int seedIdx = index(extent, origin.x, origin.y, origin.z);
        double s = 900.0;

        // THE ANCHOR: a tick-0 response is pure core, so f_halo == 0 EXACTLY
        // (not approximately -- ||D||_1 and |S| are the same sum).
        double[] tickZero = new double[cells];
        tickZero[seedIdx] = s;
        assertEquals("f_halo(0) must be EXACTLY zero for a response supported on the seed cell alone - this is the anchor the excursion bound rests on",
                     0.0, AnisotropyProbe.haloFraction(tickZero, s), 0.0);

        // A later snapshot: core still at the origin, plus a cancelling
        // +/-h halo pair placed away from it. S is untouched (the pair is
        // zero-sum); ||D||_1 grows by 2h; f_halo = 2h/(s+2h).
        double h = 12.0;
        double[] later = new double[cells];
        later[seedIdx] = s;
        later[index(extent, origin.x + 2, origin.y, origin.z)] = h;
        later[index(extent, origin.x - 2, origin.y, origin.z)] = -h;
        assertEquals("the halo pair must leave the conserved excess untouched",
                     s, AnisotropyProbe.signedTotal(later), 1e-12);
        double f = AnisotropyProbe.haloFraction(later, s);
        assertEquals("f_halo = 2h/(|S|+2h)", 2 * h / (s + 2 * h), f, 1e-12);

        for (Direction d : Direction.values()) {
            double m0 = secondMomentOf(tickZero, extent, origin, d);
            double m1 = secondMomentOf(later, extent, origin, d);
            // The two parts' CONDITIONAL moments, measured separately.
            double[] coreOnly = new double[cells];
            coreOnly[seedIdx] = s;
            double[] haloOnly = new double[cells];
            haloOnly[index(extent, origin.x + 2, origin.y, origin.z)] = h;
            haloOnly[index(extent, origin.x - 2, origin.y, origin.z)] = -h;
            double mCore = secondMomentOf(coreOnly, extent, origin, d);
            double mHalo = secondMomentOf(haloOnly, extent, origin, d);

            assertEquals("the core/halo decomposition of M must be EXACT along "
                         + d, (1 - f) * mCore + f * mHalo, m1, 1e-9);
            assertEquals("...so the halo's contribution to the measured CHANGE is exactly f*(m_halo - m_core), along "
                         + d, (m1 - m0), f * (mHalo - mCore), 1e-9);

            // The excursion bound, from the LEVEL alone.
            double projMax = maxProjectionSquared(extent, origin, d);
            assertTrue("the halo's moment excursion must be bounded by f_halo * max proj^2 along "
                       + d + ": |" + (m1 - m0) + "| vs " + (f * projMax),
                       Math.abs(m1 - m0) <= f * projMax + 1e-9);
            assertTrue("NON-VACUITY: the excursion must be genuinely nonzero along "
                       + d + ", else the bound is trivially satisfied",
                       Math.abs(m1 - m0) > 1e-6);
        }
    }

    /**
     * <b>Limb W3 fires where W1 and W2 are both silent</b> -- the property
     * that makes it a third limb rather than a restatement of the second.
     * A response whose core never leaves the seed cell but which has
     * acquired a large cancelling halo NEARBY has a low second moment
     * (W2 silent) and no mass past the half period (W1 silent), yet most of
     * the L1 mass the moment is normalized by is decorrelation structure
     * rather than transported excess.
     *
     * <p>Also pins the measured halo fraction on the REAL signed-quanta
     * substrate, so the tolerance's headroom is a number in the suite
     * rather than a claim in a javadoc.
     */
    @Test
    public void haloFractionLimbW3FiresWhereW1AndW2AreBothSilent() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int cells = extent.x * extent.y * extent.z;
        double s = 30.0;

        // Core at the origin; a cancelling +/-40 pair ONE cell away.
        double[] difference = new double[cells];
        difference[index(extent, origin.x, origin.y, origin.z)] = s;
        difference[index(extent, origin.x + 1, origin.y, origin.z)] = 40.0;
        difference[index(extent, origin.x - 1, origin.y, origin.z)] = -40.0;

        assertEquals("the conserved excess is untouched by a zero-sum halo", s,
                     AnisotropyProbe.signedTotal(difference), 1e-12);

        // W1 and W2 are both SILENT on it -- asserted, not assumed.
        AnisotropyProbe.assertWrapSafe(difference, extent, origin, 7);
        double saturation = AnisotropyProbe.momentSaturation(difference, extent,
                                                               origin);
        assertTrue("W2 must be silent here (the halo is adjacent to the origin, so the second moment stays low): "
                   + saturation,
                   saturation < AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE);
        AnisotropyProbe.assertResponseLocalized(difference, extent, origin, 7);

        double halo = AnisotropyProbe.haloFraction(difference, s);
        assertEquals("f_halo = 80/110", 80.0 / 110.0, halo, 1e-12);
        try {
            AnisotropyProbe.assertResponseCoreDominates(difference, s, 7);
            fail("expected IllegalStateException: 73% of the L1 mass is decorrelation halo, and W1/W2 both passed it");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("halo-fraction"));
            assertTrue("the diagnosis must name the tick", expected.getMessage()
                                                                   .contains("tick 7"));
        }

        // NON-VACUITY of the tolerance itself: shrink the halo until it is
        // inside the bound and confirm the limb goes quiet, so the test is
        // not merely asserting that some large number trips some threshold.
        double[] modest = new double[cells];
        modest[index(extent, origin.x, origin.y, origin.z)] = s;
        modest[index(extent, origin.x + 1, origin.y, origin.z)] = 0.35;
        modest[index(extent, origin.x - 1, origin.y, origin.z)] = -0.35;
        double modestHalo = AnisotropyProbe.haloFraction(modest, s);
        assertEquals("f_halo = 0.7/30.7", 0.7 / 30.7, modestHalo, 1e-12);
        assertTrue("the modest halo must be inside the tolerance",
                   modestHalo < AnisotropyProbe.RESPONSE_HALO_FRACTION_TOLERANCE);
        assertEquals("...and the limb must return it rather than throw",
                     modestHalo,
                     AnisotropyProbe.assertResponseCoreDominates(modest, s, 7),
                     0.0);
    }

    /**
     * The measured halo fraction on the REAL signed-quanta substrate, at
     * the same parameters as {@link
     * #matchedPairOnARealSignedBackgroundSubstrateDecorrelatesButStaysConserved}.
     * Asserted to its exact value, not merely "below tolerance": the point
     * of {@link AnisotropyProbe#RESPONSE_HALO_FRACTION_TOLERANCE} being a
     * declared regime bound is that the regime's actual number is visible.
     *
     * <p>{@code ||D||_1} reaches {@code 910} against {@code |S| = 900}, so
     * {@code f_halo = 10/910}. That is a genuine decorrelation halo
     * produced by the real dynamics, not a synthetic construction.
     */
    @Test
    public void haloFractionOnTheRealSignedSubstrateIsTheDeclaredRegimeNumber() {
        Point3i extent = new Point3i(4, 4, 4);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);

        AnisotropyProbe.MatchedPairTransport matched = AnisotropyProbe.runOneSeedMatchedPair(extent,
                                                                                              42L,
                                                                                              32,
                                                                                              30,
                                                                                              origin,
                                                                                              signedBackgroundFactory(60,
                                                                                                                       false));

        assertEquals("the halo fraction on the real substrate is exactly 10/910 - ||D||_1 reaches 910 against |S| = 900",
                     10.0 / 910.0, matched.maxHaloFraction(), 1e-12);
        assertTrue("...which must clear the tolerance with room, else the declared regime is not actually licensed: "
                   + matched.maxHaloFraction(),
                   matched.maxHaloFraction()
                   < 0.5 * AnisotropyProbe.RESPONSE_HALO_FRACTION_TOLERANCE);
        assertTrue("NON-VACUITY: the halo must be strictly positive, i.e. the runs really did decorrelate",
                   matched.maxHaloFraction() > 0);
    }

    /**
     * <b>The adversarial background class both NORM checks are blind to,
     * now closed</b> (bead inviscid-0nx.28, round-2 review -- previously
     * documented as "left open knowingly", which was declining a free fix).
     *
     * <p>A background difference of {@code -5} at the origin cell and
     * {@code +5} at one other cell leaves the injected excess
     * {@code S = 900} untouched (zero-sum) AND leaves
     * {@code ||D(.,0)||_1 = 900} untouched (every term shares a sign), so
     * it passes the runner's {@code S == 30*packetQuanta} equality and the
     * estimator's {@code ||D(.,0)||_1 == |S|} identity alike. Both of those
     * non-vacuity facts are ASSERTED here, so the test proves the new check
     * is what caught it rather than assuming so.
     *
     * <p>{@link AnisotropyProbe#assertTickZeroSupportIsTheSeedCell} closes
     * it by looking at the SUPPORT instead: every non-origin cell of
     * {@code D(.,0)} must be bit-exactly {@code 0.0}, since before any tick
     * runs both halves sum the same 30 {@code float}s per non-origin cell.
     */
    @Test
    public void matchedPairRunnerRejectsASingleSignedZeroSumBackgroundDifference() {
        Point3i extent = new Point3i(4, 4, 4);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        Point3i decoy = new Point3i(0, 0, 0);
        int packetQuanta = 30;

        SubstrateFactory adversarial = adversarialZeroSumBackgroundFactory(origin,
                                                                            decoy,
                                                                            5);

        // NON-VACUITY, measured on the very fields the runner will build:
        // both norm-based checks pass on this pair.
        double[] packetField = StructureFactor.coarseGrainedField(adversarial.create(extent,
                                                                                       42L,
                                                                                       packetQuanta,
                                                                                       origin)
                                                                              .field());
        double[] controlField = StructureFactor.coarseGrainedField(adversarial.create(extent,
                                                                                        42L,
                                                                                        0,
                                                                                        origin)
                                                                               .field());
        double[] difference = AnisotropyProbe.differenceField(packetField,
                                                               controlField,
                                                               extent.x * extent.y
                                                                           * extent.z,
                                                               0);
        assertEquals("the runner's blunt equality passes: S is exactly 30*packetQuanta",
                     30.0 * packetQuanta,
                     AnisotropyProbe.signedTotal(difference), 1e-9);
        assertEquals("the estimator's SHARP norm check passes too: every term shares a sign, so ||D||_1 == |S|",
                     30.0 * packetQuanta,
                     AnisotropyProbe.responseL1(difference), 1e-9);
        assertEquals("...and the halo fraction is therefore exactly zero, so W3 is blind to it as well",
                     0.0, AnisotropyProbe.haloFraction(difference,
                                                        30.0 * packetQuanta),
                     0.0);

        // The SUPPORT check is what catches it.
        try {
            AnisotropyProbe.runOneSeedMatchedPair(extent, 42L, 8, packetQuanta,
                                                   origin, adversarial);
            fail("expected IllegalStateException: a single-signed zero-sum background difference passes both norm checks and must be caught by the support check");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("nonzero at cell"));
            assertTrue("the diagnosis must name the offending cell",
                       expected.getMessage().contains("(0, 0, 0)"));
        }

        // ...and a well-formed pair on the same wiring still passes, so the
        // check is not simply rejecting everything this factory builds.
        AnisotropyProbe.runOneSeedMatchedPair(extent, 42L, 8, packetQuanta,
                                               origin,
                                               adversarialZeroSumBackgroundFactory(origin,
                                                                                    decoy,
                                                                                    0));
    }

    /**
     * The L1 lower-bound theorem that replaces the maximum principle as
     * the licence for {@code Math.abs} weighting: since the injected
     * excess {@code S = sum_cell D} is exactly tick-invariant (each run
     * conserves its own total) and {@code |S| &lt;= ||D||_1} by the
     * triangle inequality, signed cancellation in the difference field can
     * shrink {@code ||D||_1} at most down to {@code |S|} -- never to zero,
     * so the second moment's normalizer is never degenerate.
     *
     * <p>The inequality on its own is close to a tautology and, on random
     * fields, is never anywhere near tight -- asserting only that would
     * prove nothing about the case the theorem is invoked for. So this
     * also asserts (a) that the bound is APPROACHED ARBITRARILY CLOSELY by
     * a deliberately near-cancelling construction (which is the regime the
     * theorem must survive), and (b) that even at maximal cancellation the
     * normalizer stays at exactly {@code |S|} and the estimator therefore
     * still produces a finite second moment rather than a
     * divide-by-almost-zero.
     */
    @Test
    public void responseL1IsBoundedBelowByTheInjectedExcess() {
        Random random = new Random(20260809L);
        for (int trial = 0; trial < 200; trial++) {
            double[] difference = new double[512];
            for (int i = 0; i < difference.length; i++) {
                // Heavy cancellation on purpose: symmetric, zero-mean draws.
                difference[i] = Math.round(random.nextGaussian() * 50);
            }
            double signedTotal = AnisotropyProbe.signedTotal(difference);
            double l1 = AnisotropyProbe.responseL1(difference);
            assertTrue("||D||_1 (" + l1
                       + ") must be >= |sum D| (" + Math.abs(signedTotal)
                       + ") - the triangle inequality is what keeps the second-moment normalizer positive",
                       l1 >= Math.abs(signedTotal));
        }

        // TIGHTNESS: +m at the origin, and n cancelling +/-c pairs
        // elsewhere. ||D||_1 = m + 2nc while S = m, so the ratio
        // ||D||_1/|S| is dialled arbitrarily -- and at n=0 the bound is
        // attained with EQUALITY, which the random fields above never come
        // close to doing.
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        double m = 30.0;
        for (int pairs = 0; pairs <= 3; pairs++) {
            double[] difference = new double[extent.x * extent.y * extent.z];
            difference[index(extent, origin.x, origin.y, origin.z)] = m;
            for (int p = 0; p < pairs; p++) {
                difference[index(extent, 3, 3, 3 + p)] = 7.0;
                difference[index(extent, 5, 5, 3 + p)] = -7.0;
            }
            assertEquals("S must be unaffected by exactly-cancelling pairs",
                         m, AnisotropyProbe.signedTotal(difference), 1e-12);
            assertEquals("||D||_1 must be exactly |S| + 2*pairs*7 - at pairs=0 the triangle-inequality bound is ATTAINED, not merely respected",
                         m + 2 * pairs * 7.0,
                         AnisotropyProbe.responseL1(difference), 1e-12);
        }
    }

    /**
     * <b>The headline test.</b> A known, analytically-exact packet
     * ({@code D_[100]=4R}, {@code D_[110]=2.5R}, {@code D_[111]=2R} -&gt;
     * {@code A = 2.0}, the same algebra as {@link
     * #syntheticAnisotropicFieldIsDetected}) is buried under a large
     * SIGNED, tick-DRIFTING background that is present, identically, in
     * both halves of a matched pair. The raw packet-plus-background field
     * is what the pre-existing single-field estimator would consume; the
     * difference field is what §D-A mandates.
     *
     * <p>Two independent failure modes of the raw field are asserted, both
     * derivable rather than merely observed:
     * <ul>
     * <li><b>dilution</b> -- the raw second moment normalizes by
     * {@code 1 + ||background||_1}, which suppresses the packet's own
     * slope by that factor (thousands here);</li>
     * <li><b>corruption</b> -- the background's own drift contributes a
     * spurious, direction-dependent trend, so the raw ratio is not the
     * packet's 2.0 either.</li>
     * </ul>
     * The matched pair recovers both the per-direction slope and the ratio
     * because the background cancels EXACTLY.
     */
    @Test
    public void matchedPairRecoversPacketSpreadThroughASignedDriftingBackground() {
        double rateY = 0.2;
        double rateX = 4 * rateY;
        double[][] packetOnly = normalizedGaussianPacketField(EXTENT, ORIGIN,
                                                                TICKS,
                                                                SIGMA0_SQ,
                                                                rateX, rateY,
                                                                rateY);
        double[][] background = signedDriftingBackground(EXTENT, ORIGIN, TICKS,
                                                           4242L);
        double[][] packetRun = new double[TICKS][];
        double[][] controlRun = new double[TICKS][];
        for (int t = 0; t < TICKS; t++) {
            controlRun[t] = background[t].clone();
            packetRun[t] = background[t].clone();
            for (int i = 0; i < packetRun[t].length; i++) {
                packetRun[t][i] += packetOnly[t][i];
            }
        }

        AnisotropyProbe.MatchedPairTransport matched = AnisotropyProbe.transportEstimateMatchedPair(packetRun,
                                                                                                      controlRun,
                                                                                                      EXTENT,
                                                                                                      ORIGIN);

        assertEquals("the injected excess must be the packet's own (normalized) total mass",
                     1.0, matched.injectedExcess(), 1e-9);
        assertEquals("at tick 0 the response IS the packet, so ||D||_1 == |S|",
                     Math.abs(matched.injectedExcess()),
                     matched.responseL1First(), 1e-9);
        // W2's measured value, asserted as a BRACKET rather than "~0": the
        // packet genuinely does spread (sigma^2 reaches 1 + 7*rateX = 6.6
        // against a uniform-field moment of 133.5 on this 40^3 box), so
        // the correct expected value is 6.6/133.5 = 0.0494, NOT zero. An
        // implementation that hardcoded 0.0 on the success path -- which
        // the whole suite would otherwise accept -- fails here.
        assertEquals("W2's measured saturation must be the packet's own M_[100]/M^unif_[100] at the last tick",
                     (SIGMA0_SQ + (TICKS - 1) * rateX) / 133.5,
                     matched.maxMomentSaturation(), 5e-3);
        assertTrue("...and must still sit well inside the tolerance, was "
                   + matched.maxMomentSaturation(),
                   matched.maxMomentSaturation() < AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE);

        assertEquals("matched pair must recover the packet's exact per-axis spread rate along [100]",
                     rateX,
                     matched.transport().perDirection().get(Direction.X100)
                            .magnitude(),
                     rateX * 0.01);
        assertTrue(matched.transport().ratio().isPresent());
        assertEquals("matched pair must recover the packet's analytic anisotropy ratio",
                     2.0, matched.transport().ratio().getAsDouble(), 0.05);

        // What the pre-existing single-field estimator sees on the same
        // data. The PRIMARY claim here is an INVARIANCE, not a distance:
        // the raw estimate is (to 1 part in 1000) the SAME whether the
        // packet is present or not, i.e. it is not measuring the packet at
        // all. Stated this way it cannot pass or fail by coincidence -- an
        // earlier version of this test asserted only "raw is far from the
        // truth" and the background's own spurious trend happened to land
        // within 4% of the true rate, which would have read as success.
        EstimatorResult raw = AnisotropyProbe.transportEstimate(packetRun,
                                                                  EXTENT,
                                                                  ORIGIN);
        EstimatorResult rawWithoutPacket = AnisotropyProbe.transportEstimate(controlRun,
                                                                               EXTENT,
                                                                               ORIGIN);
        for (Direction d : Direction.values()) {
            double withPacket = raw.perDirection().get(d).magnitude();
            double withoutPacket = rawWithoutPacket.perDirection().get(d)
                                                    .magnitude();
            assertTrue("raw single-field estimate along " + d
                       + " must be blind to the packet (background dilution): with=" + withPacket
                       + " without=" + withoutPacket,
                       Math.abs(withPacket - withoutPacket)
                       < 1e-3 * Math.abs(withoutPacket));
        }
        // And it is not merely diluted but CORRUPTED: the background's own
        // drift supplies a spurious trend. Stated as an invariance again
        // rather than as a distance-from-truth -- the raw ratio equals the
        // PACKET-FREE background's ratio to 1 part in 1000, so whatever
        // number it reports is a property of the background alone. A
        // distance assertion here ("raw ratio is far from 2.0") would be
        // coincidence-prone in exactly the way the tolerance note above
        // records: the background's spurious [100] trend already lands
        // within 4% of the true rate by accident.
        assertTrue(raw.ratio().isPresent());
        assertTrue(rawWithoutPacket.ratio().isPresent());
        assertEquals("the raw single-field RATIO is a property of the background, not of the packet",
                     rawWithoutPacket.ratio().getAsDouble(),
                     raw.ratio().getAsDouble(),
                     1e-3 * rawWithoutPacket.ratio().getAsDouble());
        assertTrue("...and the matched pair must therefore report a DIFFERENT ratio from the raw one, else subtraction changed nothing",
                   Math.abs(matched.transport().ratio().getAsDouble()
                            - raw.ratio().getAsDouble()) > 0.1);
    }

    /**
     * The {@code {0, packetQuanta}} path is preserved BY CONSTRUCTION, and
     * this pins it at the substrate level rather than asserting it in
     * prose: on the Phase A zero-background substrate the control run
     * ({@code packetQuanta = 0}) holds no quanta at all, so no member is
     * ever strictly greater than another, nothing transfers, its field is
     * identically zero, and the difference field IS the packet run's
     * field. The matched-pair estimator must therefore agree with {@link
     * AnisotropyProbe#runOneSeed} BIT FOR BIT (delta 0.0, not a
     * tolerance). Parameters mirror {@code
     * SeamGoldenCompatTest#runOneSeedThroughTheSeamMatchesPinnedPhaseANumerics}
     * so this stays inside the same in-surefire cost envelope.
     */
    @Test
    public void matchedPairOnZeroBackgroundReproducesSingleFieldEstimate() {
        Point3i extent = new Point3i(6, 6, 6);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int ticks = 16;
        int packetQuanta = 100;

        EstimatorResult single = AnisotropyProbe.runOneSeed(extent, 42L, ticks,
                                                              packetQuanta,
                                                              origin)
                                                 .transport();
        AnisotropyProbe.MatchedPairTransport matched = AnisotropyProbe.runOneSeedMatchedPair(extent,
                                                                                              42L,
                                                                                              ticks,
                                                                                              packetQuanta,
                                                                                              origin,
                                                                                              AnisotropyProbe::phaseAHybridSubstrate);

        for (Direction d : Direction.values()) {
            assertEquals("matched-pair and single-field magnitudes must be BIT-IDENTICAL on a zero background, direction "
                         + d, single.perDirection().get(d).magnitude(),
                         matched.transport().perDirection().get(d).magnitude(),
                         0.0);
        }
        assertEquals("injected excess must be exactly 30 members x packetQuanta",
                     30.0 * packetQuanta, matched.injectedExcess(), 0.0);
        assertEquals("a zero-background control means D == packet field, so ||D||_1 at tick 0 == |S|",
                     30.0 * packetQuanta, matched.responseL1First(), 0.0);
    }

    /**
     * Subtracting a run from itself is the ONE way to reach the
     * second-moment normalizer's degenerate branch (the L1 bound forbids
     * every other route). It must be refused up front rather than yielding
     * a vacuous all-zero magnitude set that {@code ratio()} would then
     * report as a legitimate-looking DEGENERATE.
     */
    @Test
    public void matchedPairRejectsAZeroInjectedExcess() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        double[][] run = new double[4][extent.x * extent.y * extent.z];
        for (int t = 0; t < run.length; t++) {
            run[t][index(extent, origin.x, origin.y, origin.z)] = 30.0;
        }

        try {
            AnisotropyProbe.transportEstimateMatchedPair(run, run, extent,
                                                          origin);
            fail("expected IllegalArgumentException for a zero injected excess (a run subtracted from itself)");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("ZERO injected excess"));
        }
    }

    /**
     * Invariant (I1) as an enforced precondition, not a comment: each run
     * conserves its own quanta total exactly, so the difference of the two
     * totals cannot drift. A drift means the pair was not driven on a
     * shared conserved background -- a broken measurement, which must fail
     * loudly rather than produce a plausible-looking number.
     */
    @Test
    public void matchedPairRejectsADriftingInjectedExcess() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int cells = extent.x * extent.y * extent.z;
        double[][] packetRun = new double[3][cells];
        double[][] controlRun = new double[3][cells];
        int seedIdx = index(extent, origin.x, origin.y, origin.z);
        for (int t = 0; t < 3; t++) {
            packetRun[t][seedIdx] = 30.0 + t; // gains quanta from nowhere
        }

        try {
            AnisotropyProbe.transportEstimateMatchedPair(packetRun, controlRun,
                                                          extent, origin);
            fail("expected IllegalStateException: the injected excess drifted from 30.0 to 32.0 across ticks");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("not tick-invariant"));
        }
    }

    /**
     * The matched-pair runner's same-background check. A factory whose
     * background draws depend on {@code packetQuanta} silently breaks the
     * pair: both halves still conserve their own totals, so the
     * tick-invariance precondition passes, and only the {@code
     * S == 30*packetQuanta} equality catches it.
     */
    @Test
    public void matchedPairRunnerRejectsAFactoryWhoseBackgroundDependsOnPacketQuanta() {
        Point3i extent = new Point3i(4, 4, 4);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        // Seeds the control half with quanta of its own instead of the
        // requested zero -- exactly the "different background per half"
        // mistake. Both halves still conserve their own totals exactly, so
        // the tick-invariance precondition passes and only the
        // S == 30*packetQuanta equality can catch this.
        SubstrateFactory rogue = (ext, seed, packetQuanta,
                                   originCell) -> AnisotropyProbe.phaseAHybridSubstrate(ext,
                                                                                          seed,
                                                                                          packetQuanta == 0
                                                                                              ? 5
                                                                                              : packetQuanta,
                                                                                          originCell);

        try {
            AnisotropyProbe.runOneSeedMatchedPair(extent, 42L, 4, 10, origin,
                                                   rogue);
            fail("expected IllegalStateException: the two halves of the pair were not on a shared background");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("not on a shared background"));
        }
    }

    /**
     * <b>The ordering regression.</b> The first version of {@code
     * runOneSeedMatchedPair} ran the estimator FIRST and checked
     * {@code S == 30*packetQuanta} afterwards. On the realistic broken
     * pair -- the two halves carrying genuinely DIFFERENT quanta
     * backgrounds -- the difference field is the whole background
     * difference, which is delocalized, so wrap-safety limb W2 fired first
     * and the caller was told to "reduce tick count or enlarge extent":
     * the wrong root cause, reported on exactly the input the
     * shared-background check exists to catch.
     *
     * <p>Non-vacuity is the point of this test, so it is asserted rather
     * than argued: the tick-0 difference field is measured to exceed the
     * W2 tolerance, i.e. W2 genuinely WOULD have fired here, and the
     * message is nevertheless the shared-background one.
     */
    @Test
    public void matchedPairRunnerDiagnosesAMismatchedBackgroundBeforeWrapSaturation() {
        Point3i extent = new Point3i(6, 6, 6);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        // A signed-background factory whose BACKGROUND DRAWS depend on
        // packetQuanta -- precisely the contract violation
        // runOneSeedMatchedPair's javadoc forbids. The tick-0 difference
        // is then the whole (delocalized) background difference.
        SubstrateFactory rogue = signedBackgroundFactory(3, true);

        // Non-vacuity: the tick-0 difference really does trip W2.
        double[] packetFirst = StructureFactor.coarseGrainedField(rogue.create(extent,
                                                                                42L,
                                                                                10,
                                                                                origin)
                                                                        .field());
        double[] controlFirst = StructureFactor.coarseGrainedField(rogue.create(extent,
                                                                                 42L,
                                                                                 0,
                                                                                 origin)
                                                                         .field());
        double[] difference = AnisotropyProbe.differenceField(packetFirst,
                                                               controlFirst,
                                                               extent.x
                                                                       * extent.y
                                                                       * extent.z,
                                                               0);
        assertTrue("fixture must genuinely trip W2, else this test cannot show the ordering matters - saturation was "
                   + AnisotropyProbe.momentSaturation(difference, extent,
                                                       origin),
                   AnisotropyProbe.momentSaturation(difference, extent, origin)
                   > AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE);

        try {
            AnisotropyProbe.runOneSeedMatchedPair(extent, 42L, 4, 10, origin,
                                                   rogue);
            fail("expected IllegalStateException for a mismatched background");
        } catch (IllegalStateException expected) {
            assertTrue("the diagnosis must be the shared-background one, not wrap-saturation: "
                       + expected.getMessage(),
                       expected.getMessage().contains("not on a shared background"));
            assertFalse("the wrap-saturation message is the WRONG root cause here: "
                        + expected.getMessage(),
                        expected.getMessage().contains("wrap-saturation"));
        }
    }

    /**
     * The SHARP same-background check, asserted at the estimator level:
     * before any tick runs the two halves can differ only by the seeded
     * packet, so {@code D(.,0)} must be supported on the seed cell alone
     * and {@code ||D(.,0)||_1} must equal {@code |S|}. Uncancelled
     * background shows up as surplus L1 mass.
     *
     * <p>This is the check the derivation memo originally described in the
     * record javadoc and then merely RETURNED instead of asserting, while
     * the enforcement rested on the runner's much blunter signed-sum
     * equality. The fixture here is built so the blunt check would NOT
     * fire: the background difference sums to exactly zero, so
     * {@code S == 30*packetQuanta} still holds and only the L1 identity
     * catches it. As with the runner test, W2 is measured to confirm it
     * would otherwise have pre-empted the diagnosis.
     */
    @Test
    public void matchedPairRejectsAnUncancelledBackgroundAtTickZero() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int cells = extent.x * extent.y * extent.z;
        double[][] packetRun = new double[3][cells];
        double[][] controlRun = new double[3][cells];
        int seedIdx = index(extent, origin.x, origin.y, origin.z);
        Random random = new Random(31337L);
        // A zero-sum background difference: the signed-sum check is blind
        // to it by construction.
        double[] delta = new double[cells];
        double sum = 0;
        for (int i = 0; i < cells - 1; i++) {
            delta[i] = random.nextInt(21) - 10;
            sum += delta[i];
        }
        delta[cells - 1] = -sum;
        for (int t = 0; t < 3; t++) {
            for (int i = 0; i < cells; i++) {
                packetRun[t][i] = delta[i];
            }
            packetRun[t][seedIdx] += 900.0;
        }

        assertEquals("fixture must keep the signed sum at exactly 30*packetQuanta, so ONLY the L1 identity can catch it",
                     900.0,
                     AnisotropyProbe.signedTotal(AnisotropyProbe.differenceField(packetRun[0],
                                                                                  controlRun[0],
                                                                                  cells,
                                                                                  0)),
                     1e-9);
        assertTrue("fixture must genuinely trip W2, else the ordering claim is untested",
                   AnisotropyProbe.momentSaturation(packetRun[0], extent,
                                                     origin)
                   > AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE);

        try {
            AnisotropyProbe.transportEstimateMatchedPair(packetRun, controlRun,
                                                          extent, origin);
            fail("expected IllegalStateException: the tick-0 difference field carries uncancelled background");
        } catch (IllegalStateException expected) {
            assertTrue("the diagnosis must be the shared-background one, not wrap-saturation: "
                       + expected.getMessage(),
                       expected.getMessage().contains("not on a shared background"));
            assertFalse(expected.getMessage(),
                        expected.getMessage().contains("wrap-saturation"));
        }
    }

    /**
     * {@code responseL1Last} is the designated disambiguator for the two
     * DIFFERENT physical situations that both surface as a degenerate
     * (empty) ratio on the matched-pair path, and it was previously
     * asserted by no test and read by no caller -- a diagnostic nobody had
     * ever shown could diagnose anything. Both branches are exercised
     * here:
     * <ul>
     * <li>a response that NEVER LEFT the seed cell: every slope is exactly
     * {@code 0.0}, the ratio is empty, and
     * {@code responseL1Last == |S|};</li>
     * <li>a response that SPREAD but with no net second-moment trend
     * (it moves out and back, ending where it started): every slope is
     * again {@code 0.0} and the ratio is again empty, but
     * {@code responseL1Last >> |S|}.</li>
     * </ul>
     * The two are indistinguishable from the {@link
     * AnisotropyProbe.EstimatorResult} alone; {@code responseL1Last} is
     * what separates them.
     */
    @Test
    public void responseL1LastSeparatesTheTwoDegenerateRatioDiagnoses() {
        Point3i extent = new Point3i(8, 8, 8);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int cells = extent.x * extent.y * extent.z;
        int seedIdx = index(extent, origin.x, origin.y, origin.z);
        int ticks = 4;

        double[][] control = new double[ticks][cells];

        // Branch 1: never leaves the seed cell.
        double[][] frozen = new double[ticks][cells];
        for (int t = 0; t < ticks; t++) {
            frozen[t][seedIdx] = 30.0;
        }
        AnisotropyProbe.MatchedPairTransport still = AnisotropyProbe.transportEstimateMatchedPair(frozen,
                                                                                                    control,
                                                                                                    extent,
                                                                                                    origin);
        assertFalse("a frozen response has an exactly-zero slope in every direction, so the ratio must be degenerate",
                    still.transport().ratio().isPresent());
        assertEquals("...and responseL1Last must equal |S| exactly: the response never left the seed cell",
                     Math.abs(still.injectedExcess()), still.responseL1Last(),
                     0.0);

        // Branch 2: genuinely spreads, but along z only. Every cell it
        // occupies has dx == dy == 0, so proj_[100] and proj_[110] are
        // identically zero there and M_[100](t) == M_[110](t) == 0 at every
        // tick -- exactly-zero slopes, so ratio() is degenerate for the
        // same reason a frozen response is. The +/-3 pair keeps S fixed
        // while ||D||_1 grows by 6, and tick 0 stays origin-only so the
        // shared-background identity still holds.
        //
        // The injected excess here is 150, not 30, BECAUSE OF LIMB W3: a
        // +/-3 halo against |S| = 30 reads f_halo = 6/36 = 0.167, four
        // times the tolerance, and W3 would (correctly) refuse the fixture
        // as one whose moment is mostly decorrelation structure. At
        // |S| = 150 the same halo reads 6/156 = 0.0385, inside the licensed
        // window -- asserted below so the interaction is visible rather
        // than a magic constant.
        double[][] spreading = new double[ticks][cells];
        for (int t = 0; t < ticks; t++) {
            spreading[t][seedIdx] = 150.0;
            if (t > 0) {
                spreading[t][index(extent, origin.x, origin.y,
                                    origin.z + 2)] = 3.0;
                spreading[t][index(extent, origin.x, origin.y,
                                    origin.z - 2)] = -3.0;
            }
        }
        AnisotropyProbe.MatchedPairTransport moved = AnisotropyProbe.transportEstimateMatchedPair(spreading,
                                                                                                    control,
                                                                                                    extent,
                                                                                                    origin);
        assertFalse("a response whose spread contributes nothing to the [100]/[110] moments leaves those slopes exactly zero, so the ratio is degenerate here too",
                    moved.transport().ratio().isPresent());
        assertEquals("...but responseL1Last is 150 + 3 + 3 = 156, above |S| - which is exactly how the two cases are told apart",
                     156.0, moved.responseL1Last(), 1e-12);
        assertTrue("the two diagnoses must be separated by responseL1Last, since the EstimatorResult alone cannot tell them apart",
                   moved.responseL1Last() > still.responseL1Last());
        assertEquals("a frozen response has NO halo at all", 0.0,
                     still.maxHaloFraction(), 0.0);
        assertEquals("the spread response's halo is 6/156, inside W3's licensed window - which is why this fixture uses |S| = 150",
                     6.0 / 156.0, moved.maxHaloFraction(), 1e-12);
        assertTrue("...and W3's own limb is what makes maxHaloFraction the second discriminator alongside responseL1Last",
                   moved.maxHaloFraction()
                   < AnisotropyProbe.RESPONSE_HALO_FRACTION_TOLERANCE);
    }

    /**
     * <b>The matched pair on a real signed-background SUBSTRATE, in the
     * NONLINEAR regime.</b> Closes the gap the headline decontamination
     * test leaves open: that test adds a synthetic background IDENTICALLY
     * to both halves, so the background cancels exactly BY CONSTRUCTION --
     * the LINEAR case, which the derivation's own §2 says does not obtain.
     * Every other matched-pair test in this class uses either a synthetic
     * additive field or the zero-background Phase A substrate, so until
     * this test the matched-pair path had never seen a nonzero quanta
     * background at all, and the decorrelation risk limb W2 exists for was
     * untested.
     *
     * <p>Here both halves run the SAME uniform signed background through
     * the real {@code QuantaExchangeRule} dynamics, and the packet
     * perturbs which contacts fire. Three claims are checked:
     * <ul>
     * <li><b>(I1) survives a signed background exactly.</b> {@code S}
     * comes out at exactly {@code 30*packetQuanta} and does not drift by a
     * single quantum across 32 ticks -- asserted with delta {@code 0.0},
     * not a tolerance.</li>
     * <li><b>The pair is well-formed on a NONZERO background.</b>
     * {@code ||D(.,0)||_1 == |S|}, i.e. the background really did cancel
     * at tick 0 -- which is what {@code seedPacket}'s ACCUMULATE semantics
     * buys and which ASSIGN would have broken.</li>
     * <li><b>The runs genuinely DECORRELATE.</b> {@code ||D||_1} rises
     * strictly above {@code |S|} by the last tick: evolution does not
     * commute with subtraction, the difference field acquires signed
     * structure away from the seed cell, and no upper bound on
     * {@code ||D||_1} survives. This is the derivation's central
     * nonlinearity claim, measured rather than argued -- and it is exactly
     * why localization of {@code D} has to be asserted (limb W2) rather
     * than assumed.</li>
     * </ul>
     */
    @Test
    public void matchedPairOnARealSignedBackgroundSubstrateDecorrelatesButStaysConserved() {
        Point3i extent = new Point3i(4, 4, 4);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int ticks = 32;
        int packetQuanta = 30;

        AnisotropyProbe.MatchedPairTransport matched = AnisotropyProbe.runOneSeedMatchedPair(extent,
                                                                                              42L,
                                                                                              ticks,
                                                                                              packetQuanta,
                                                                                              origin,
                                                                                              signedBackgroundFactory(60,
                                                                                                                       false));

        assertEquals("(I1) must hold EXACTLY on a signed background: S == 30*packetQuanta, no drift over 32 ticks",
                     30.0 * packetQuanta, matched.injectedExcess(), 0.0);
        assertEquals("the background must cancel exactly at tick 0, so D is supported on the seed cell alone",
                     30.0 * packetQuanta, matched.responseL1First(), 0.0);
        assertTrue("the two runs must genuinely decorrelate on a signed background -- ||D||_1 must exceed |S| by the last tick, was "
                   + matched.responseL1Last() + " against |S| = "
                   + matched.responseL1First(),
                   matched.responseL1Last() > matched.responseL1First());
        assertTrue("...and W2 must still have headroom in this regime, saturation was "
                   + matched.maxMomentSaturation(),
                   matched.maxMomentSaturation() < AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE);
        assertTrue("the estimator must produce a usable ratio on a signed-background substrate",
                   matched.transport().ratio().isPresent());
    }

    /**
     * The {@code +=} vs {@code =} decision in {@code
     * AnisotropyProbe.seedPacket}, pinned rather than left as a comment
     * (bead inviscid-0nx.28 fix round; E.2 builds directly on it). On the
     * Phase A zero-quanta substrate the two are indistinguishable, which is
     * why the decision needs a test on a NONZERO background to have any
     * content at all.
     *
     * <p>ACCUMULATE is the chosen semantics: the packet is an EXCESS on top
     * of whatever background a factory supplied, so the matched pair's
     * difference field is exactly {@code 30*packetQuanta} at the origin and
     * zero elsewhere. ASSIGN would have made the packet half "background
     * with a hole at the origin", whose injected excess is
     * {@code 30*packetQuanta - background(origin)} -- tripping the runner's
     * own equality on a well-formed pair.
     */
    @Test
    public void seedPacketAccumulatesOntoTheBackgroundRatherThanOverwritingIt() {
        Point3i extent = new Point3i(4, 4, 4);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int slots = 30 * extent.x * extent.y * extent.z;

        // NOTE: Necronomata stores the frequency array BY REFERENCE, so
        // each automaton gets its own copy of the background.
        float[] background = new float[slots];
        Arrays.fill(background, 7.0f);
        Necronomata automaton = new Necronomata(new float[slots], extent,
                                                  background.clone());

        AnisotropyProbe.seedPacket(automaton, origin, 30);

        int base = automaton.indexOfCell(origin);
        for (int m = 0; m < 30; m++) {
            assertEquals("the packet must ACCUMULATE onto the background (7 + 30), not overwrite it (30), member "
                         + m, 37L, automaton.quantaAt(base + m));
        }
        // ...and nothing outside the origin cell is touched.
        for (int s = 0; s < slots; s++) {
            if (s < base || s >= base + 30) {
                assertEquals("seedPacket must touch only the origin cell, slot "
                             + s, 7L, automaton.quantaAt(s));
            }
        }

        // The consequence that matters: the matched-pair difference at the
        // origin is exactly 30*packetQuanta, which is what both
        // same-background checks are written against.
        Necronomata controlAutomaton = new Necronomata(new float[slots], extent,
                                                        background.clone());
        AnisotropyProbe.seedPacket(controlAutomaton, origin, 0);
        double[] packetField = StructureFactor.coarseGrainedField(automaton);
        double[] controlField = StructureFactor.coarseGrainedField(controlAutomaton);
        double[] difference = AnisotropyProbe.differenceField(packetField,
                                                               controlField,
                                                               extent.x
                                                                       * extent.y
                                                                       * extent.z,
                                                               0);
        assertEquals("the injected excess on a NONZERO background must still be exactly 30*packetQuanta",
                     900.0, AnisotropyProbe.signedTotal(difference), 1e-9);
        assertEquals("...and the difference must be supported on the seed cell alone, so ||D||_1 == |S|",
                     900.0, AnisotropyProbe.responseL1(difference), 1e-9);
    }

    /**
     * <b>The census that replaces the retired maximum principle.</b>
     * Because {@code (I2)} interval invariance is FALSE for this automaton
     * (snapshot resolution accumulates same-tick deltas additively, so a
     * member at {@code q=1} contacted twice in one tick lands at
     * {@code -1}), nothing PROVES the {@code {0, packetQuanta}} field is
     * non-negative, and so nothing proves {@code Math.abs} was the
     * identity on it. Green tests could not settle that: no assertion in
     * the suite had ever looked.
     *
     * <p>The full Phase A campaign (8 seeds {@code 42..49}, {@code 128}
     * ticks, {@code 8^3}) censuses to ZERO negative members and ZERO
     * negative cells over all 1024 snapshots, with a minimum member value
     * of {@code 0} -- so the committed {@code
     * anisotropy-report-phaseA.tsv} weights no negative cell positively.
     * That run costs minutes, so it is not performed here; it is
     * REPRODUCIBLE FROM THIS TREE via {@link AnisotropyProbe#main} with
     * {@code --census}, which prints every number this javadoc quotes.
     *
     * <p><b>What is pinned HERE is a reduced-scope tripwire, and its
     * coverage is biased LOW.</b> Seeds {@code {42,43}} x ticks
     * {@code 0..23} is 48 of 1024 snapshots ({@code 4.69%}), 1482 of 39869
     * occupied member observations ({@code 3.72%}), and -- the number that
     * matters, since only a member at {@code q == 1} can be driven negative
     * by one extra same-tick contact -- 42 of 5842 at-risk observations,
     * {@code 0.72%}. The at-risk population by 16-tick bucket across all
     * seeds runs {@code 159 / 471 / 875 / 1012 / 1013 / 875 / 738 / 699},
     * so this window is the FIRST bucket: the bottom of a risk curve that
     * plateaus about six times higher at ticks 48-79. A regression that
     * first drives a member negative after tick 24, or on seeds 44-49,
     * passes here silently. This is a tripwire on the regime, not a proof
     * about the campaign, and the reduction is toward the SAFE end.
     *
     * <p>Also measured in the same pass, at zero extra automaton cost: the
     * W2 moment saturation. On the Phase A zero-background substrate the
     * control run holds no quanta, so the difference field IS the packet
     * field and its saturation can be read directly. This bears on the
     * acknowledged gap in the campaign path -- {@link
     * AnisotropyProbe#runCampaign} still carries only the proven-vacuous
     * limb W1. The full campaign's worst-case saturation is {@code 0.00512}
     * (49x under the tolerance), but that is the {@code --census} number;
     * what THIS test asserts is the 24-tick value, and saturation grows
     * monotonically with ticks, so the assertion below is a LOWER BOUND on
     * the campaign worst case rather than the worst case itself.
     */
    @Test
    public void phaseARegimeIsNonNegativeAndFarFromSaturationByCensus() {
        Point3i extent = AnisotropyProbe.DEFAULT_EXTENT;
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        int packetQuanta = AnisotropyProbe.DEFAULT_PACKET_QUANTA;
        int ticks = 24;

        long negativeMembers = 0;
        long negativeCells = 0;
        long minMember = 0;
        double worstSaturation = 0;
        double totalQuantaSeen = 0;

        for (long seed : new long[] { 42L, 43L }) {
            SubstrateFactory.Substrate substrate = AnisotropyProbe.phaseAHybridSubstrate(extent,
                                                                                           seed,
                                                                                           packetQuanta,
                                                                                           origin);
            for (int tick = 0; tick < ticks; tick++) {
                if (tick > 0) {
                    substrate.run().tick(tick - 1);
                }
                for (int s = 0; s < substrate.field().slotCount(); s++) {
                    long q = substrate.field().quantaAt(s);
                    if (q < 0) {
                        negativeMembers++;
                    }
                    minMember = Math.min(minMember, q);
                }
                double[] cells = StructureFactor.coarseGrainedField(substrate.field());
                for (double v : cells) {
                    if (v < 0) {
                        negativeCells++;
                    }
                    totalQuantaSeen += Math.abs(v);
                }
                worstSaturation = Math.max(worstSaturation,
                                            AnisotropyProbe.momentSaturation(cells,
                                                                              extent,
                                                                              origin));
            }
        }

        // Non-vacuity: the census must actually have had quanta to look at.
        assertEquals("the census must see the conserved packet total at every one of the 2*24 snapshots",
                     2 * ticks * 30.0 * packetQuanta, totalQuantaSeen, 1e-6);

        assertEquals("Math.abs must be the IDENTITY on the {0,packetQuanta} path: no member may go negative",
                     0L, negativeMembers);
        assertEquals("...and therefore no coarse-grained cell may go negative",
                     0L, negativeCells);
        assertEquals("the minimum member value on this path must be 0", 0L,
                     minMember);

        assertTrue("the Phase A regime must sit far below limb W2's tolerance, was "
                   + worstSaturation,
                   worstSaturation < 0.1
                                      * AnisotropyProbe.RESPONSE_MOMENT_SATURATION_TOLERANCE);
    }

    // ------------------------------------------------------------------
    // Committed Phase A artifact -- structural validation only (the
    // artifact is generated by AnisotropyProbe.main(), NOT regenerated
    // in surefire -- see CommittedContactAtlasTest's precedent).
    // ------------------------------------------------------------------

    @Test
    public void committedPhaseAArtifactIsWellFormed() throws IOException {
        URL resource = AnisotropyProbeTest.class.getClassLoader()
                                                  .getResource(RESOURCE_PATH);
        if (resource == null) {
            fail("regenerate with AnisotropyProbe.main() and review the diff: "
               + "committed report src/test/resources/lga/anisotropy-report-phaseA.tsv is missing");
        }
        Path path = Paths.get(resource.getPath());
        List<String> lines = Files.readAllLines(path);

        long directionRows = lines.stream()
                                   .filter(l -> l.startsWith("DIRECTION\t"))
                                   .count();
        long summaryRows = lines.stream()
                                 .filter(l -> l.startsWith("SUMMARY\t"))
                                 .count();
        long collisionsRows = lines.stream()
                                    .filter(l -> l.startsWith("COLLISIONS\t"))
                                    .count();
        long pooledDirectionRows = lines.stream()
                                         .filter(l -> l.startsWith("POOLED_DIRECTION\t"))
                                         .count();
        long pooledSummaryRows = lines.stream()
                                       .filter(l -> l.startsWith("POOLED_SUMMARY\t"))
                                       .count();
        assertEquals("expected 8 seeds x 2 estimators x 3 directions DIRECTION rows",
                     8 * 2 * 3, directionRows);
        assertEquals("expected 2 SUMMARY rows (TRANSPORT, SPECTRAL) - diagnostic, see naivePerSeedRatioCaveat header",
                     2, summaryRows);
        assertEquals("expected 8 COLLISIONS rows (one per seed, FIX 2)", 8,
                     collisionsRows);
        assertEquals("expected 2 estimators x 3 directions POOLED_DIRECTION rows (FIX 1)",
                     2 * 3, pooledDirectionRows);
        assertEquals("expected 2 POOLED_SUMMARY rows (TRANSPORT, SPECTRAL) - the significance statistic (FIX 1)",
                     2, pooledSummaryRows);

        boolean hasGitCommit = lines.stream()
                                     .anyMatch(l -> l.startsWith("# gitCommit=")
                                                     && l.length() > "# gitCommit=".length());
        assertTrue("gitCommit provenance header must be populated", hasGitCommit);
        boolean hasSeeds = lines.stream()
                                 .anyMatch(l -> l.startsWith("# seeds=42,43,44,45,46,47,48,49"));
        assertTrue("expected the literal 8-seed list 42..49 in the provenance header",
                   hasSeeds);
        boolean hasSmallNFlag = lines.stream()
                                      .anyMatch(l -> l.startsWith("# smallNEarlyTimeFlag="));
        assertTrue("expected the FIX 2 small-N/early-time framing header",
                   hasSmallNFlag);
        boolean hasNaiveCaveat = lines.stream()
                                       .anyMatch(l -> l.startsWith("# naivePerSeedRatioCaveat="));
        assertTrue("expected the order-statistic-bias caveat header pointing at POOLED_SUMMARY as the significance statistic",
                   hasNaiveCaveat);
        boolean hasSpectralFraming = lines.stream()
                                           .anyMatch(l -> l.startsWith("# spectralZeroSlopeFraming="));
        assertTrue("expected the FIX 3 spectral-zero-slope-is-expected-diffusive-signature framing header",
                   hasSpectralFraming);
    }

    // ------------------------------------------------------------------
    // Synthetic field generator + helpers.
    // ------------------------------------------------------------------

    private static final Point3i EXTENT    = new Point3i(40, 40, 40);
    private static final Point3i ORIGIN    = new Point3i(20, 20, 20);
    private static final int     TICKS     = 8;
    private static final double  SIGMA0_SQ = 1.0;
    private static final long[]  SEEDS     = { 42L, 43L, 44L, 45L, 46L, 47L,
                                                48L, 49L };

    /**
     * A synthetic mass distribution: a separable product of per-axis
     * Gaussian-shaped weights centered at {@code origin}, with per-axis
     * variance {@code sigma0Sq + t*rateAxis} at tick {@code t} (plus a
     * small deterministic seeded jitter on {@code sigma0Sq} only, never
     * on the rates -- a nuisance-parameter perturbation, not a signal
     * perturbation). Because the weight is separable
     * ({@code w(i,j,k)=fx(i)*fy(j)*fz(k)}), x/y/z are EXACTLY independent
     * under the discrete mass distribution this induces (no cross term),
     * so the projected variances combine by simple averaging - see the
     * two tests above for the derivation this relies on. {@code extent}
     * (40 per axis, margin from a max sigma of a few units) is chosen
     * generously so the packet's moment saturation stays far below {@link
     * AnisotropyProbe#RESPONSE_MOMENT_SATURATION_TOLERANCE} - the
     * wrap-safety precondition would otherwise (correctly) reject this
     * fixture.
     */
    private static double[][] gaussianPacketField(Point3i extent,
                                                    Point3i origin, int ticks,
                                                    double sigma0Sq,
                                                    double rateX,
                                                    double rateY,
                                                    double rateZ, long seed,
                                                    double sigma0Jitter) {
        Random random = new Random(seed);
        double jitteredSigma0Sq = sigma0Sq
                                   + (random.nextDouble() * 2 - 1)
                                     * sigma0Jitter;
        double[][] fieldByTick = new double[ticks][extent.x * extent.y
                                                    * extent.z];
        for (int t = 0; t < ticks; t++) {
            double sx2 = jitteredSigma0Sq + t * rateX;
            double sy2 = jitteredSigma0Sq + t * rateY;
            double sz2 = jitteredSigma0Sq + t * rateZ;
            double[] field = fieldByTick[t];
            for (int i = 0; i < extent.x; i++) {
                double dx = i - origin.x;
                double wx = Math.exp(-(dx * dx) / (2 * sx2));
                for (int j = 0; j < extent.y; j++) {
                    double dy = j - origin.y;
                    double wy = Math.exp(-(dy * dy) / (2 * sy2));
                    for (int k = 0; k < extent.z; k++) {
                        double dz = k - origin.z;
                        double wz = Math.exp(-(dz * dz) / (2 * sz2));
                        field[(i * extent.y + j) * extent.z
                              + k] = wx * wy * wz;
                    }
                }
            }
        }
        return fieldByTick;
    }

    private static int index(Point3i extent, int i, int j, int k) {
        return (i * extent.y + j) * extent.z + k;
    }

    /**
     * A {@link SubstrateFactory} on a uniform SIGNED quanta background --
     * the §D-A regime, which no in-tree factory supplies yet. Mirrors
     * {@code AnisotropyProbe.phaseAHybridSubstrate}'s wiring exactly, with
     * one change: the background is loaded into the {@code Necronomata}
     * BEFORE the {@link ConservationAudit} is constructed, since the audit
     * snapshots the lattice total at construction and a background applied
     * afterwards would (correctly) be reported as a conservation violation.
     *
     * <p>Draw order is {@code angles -> background -> packet}, and the
     * background draws come from their own generator so they do NOT depend
     * on {@code packetQuanta} -- the contract {@code
     * AnisotropyProbe.runOneSeedMatchedPair} requires. Passing
     * {@code rogue = true} VIOLATES that contract on purpose (the control
     * half's background is shifted by {@code +1} per member), which is what
     * the ordering regression needs.
     *
     * @param amplitude per-member background draws are uniform on
     *                  {@code [-amplitude, +amplitude]}
     */
    private static SubstrateFactory signedBackgroundFactory(int amplitude,
                                                              boolean rogue) {
        return (extent, seed, packetQuanta, originCell) -> {
            int slots = 30 * extent.x * extent.y * extent.z;
            float[] angles = new float[slots];
            Random angleRandom = new Random(seed);
            for (int i = 0; i < slots; i++) {
                angles[i] = angleRandom.nextFloat() * (float) (2 * Math.PI);
            }
            float[] frequency = new float[slots];
            Random backgroundRandom = new Random(seed ^ 0x5EEDL);
            int bias = rogue && packetQuanta == 0 ? 1 : 0;
            for (int i = 0; i < slots; i++) {
                frequency[i] = backgroundRandom.nextInt(2 * amplitude + 1)
                                - amplitude + bias;
            }
            Necronomata automaton = new Necronomata(angles, extent, frequency);
            AnisotropyProbe.seedPacket(automaton, originCell, packetQuanta);

            FccNeighborhood neighborhood = new FccNeighborhood(automaton);
            ContactPredicate predicate = new ContactPredicate(new MemberGeometry(ContactAtlasGenerator.GEOMETRY_RESOLUTION,
                                                                                   ContactAtlasGenerator.RADIUS));
            ContactScan scan = new ContactScan(automaton, neighborhood,
                                                predicate);
            CollisionStatistics statistics = new CollisionStatistics();
            CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                       new QuantaExchangeRule(),
                                                       statistics);
            HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);
            ConservationAudit audit = new ConservationAudit(automaton);
            return new SubstrateFactory.Substrate(automaton,
                                                   new AuditedRun(hybrid,
                                                                   audit),
                                                   statistics);
        };
    }

    /**
     * Limb W2's RETIRED first measure, kept here (and only here) so the
     * regression that retired it checks the blindness claim by
     * recomputation rather than by recounting the argument: the fraction of
     * a field's L1 mass at or beyond the half-period shell, i.e. with some
     * axis at {@code |coord-origin| >= extentAxis/2}. Deliberately NOT in
     * production -- it is not a guard any more.
     */
    private static double legacyShellMassFraction(double[] field,
                                                    Point3i extent,
                                                    Point3i origin) {
        double total = 0;
        double shell = 0;
        int halfX = extent.x / 2;
        int halfY = extent.y / 2;
        int halfZ = extent.z / 2;
        for (int i = 0; i < extent.x; i++) {
            int dx = Math.abs(i - origin.x);
            for (int j = 0; j < extent.y; j++) {
                int dy = Math.abs(j - origin.y);
                for (int k = 0; k < extent.z; k++) {
                    double mass = Math.abs(field[index(extent, i, j, k)]);
                    if (mass == 0.0) {
                        continue;
                    }
                    total += mass;
                    int dz = Math.abs(k - origin.z);
                    if (dx >= halfX || dy >= halfY || dz >= halfZ) {
                        shell += mass;
                    }
                }
            }
        }
        return total > 0 ? shell / total : 0.0;
    }

    /**
     * {@code max_cell proj_d(cell-origin)^2} over the cells W1 admits --
     * i.e. those with {@code |d_a| <= extentAxis/2} on every axis. The
     * ceiling the halo's moment excursion is bounded by.
     */
    private static double maxProjectionSquared(Point3i extent, Point3i origin,
                                                 Direction d) {
        double worst = 0;
        for (int i = 0; i < extent.x; i++) {
            if (Math.abs(i - origin.x) > extent.x / 2) {
                continue;
            }
            for (int j = 0; j < extent.y; j++) {
                if (Math.abs(j - origin.y) > extent.y / 2) {
                    continue;
                }
                for (int k = 0; k < extent.z; k++) {
                    if (Math.abs(k - origin.z) > extent.z / 2) {
                        continue;
                    }
                    double p = projection(i - origin.x, j - origin.y,
                                           k - origin.z, d);
                    worst = Math.max(worst, p * p);
                }
            }
        }
        return worst;
    }

    /**
     * A {@link SubstrateFactory} whose two halves differ by a SINGLE-SIGNED
     * ZERO-SUM background perturbation -- the adversarial class that both
     * norm-based same-background checks are structurally blind to (see
     * {@link #matchedPairRunnerRejectsASingleSignedZeroSumBackgroundDifference}).
     *
     * <p>The CONTROL half (the one built with {@code packetQuanta == 0})
     * gets {@code +offset} on one member of the origin cell and
     * {@code -offset} on one member of {@code decoy}. The tick-0 difference
     * is therefore {@code 30*packetQuanta - offset} at the origin and
     * {@code +offset} at {@code decoy}: same sign, same sum, same L1 norm
     * as a well-formed pair. {@code offset == 0} yields a well-formed
     * factory on identical wiring, which the test uses as its control.
     */
    private static SubstrateFactory adversarialZeroSumBackgroundFactory(Point3i originCell,
                                                                          Point3i decoy,
                                                                          int offset) {
        return (extent, seed, packetQuanta, unusedOrigin) -> {
            int slots = 30 * extent.x * extent.y * extent.z;
            float[] angles = new float[slots];
            Random angleRandom = new Random(seed);
            for (int i = 0; i < slots; i++) {
                angles[i] = angleRandom.nextFloat() * (float) (2 * Math.PI);
            }
            float[] frequency = new float[slots];
            Arrays.fill(frequency, 3.0f);
            Necronomata automaton = new Necronomata(angles, extent, frequency);
            if (packetQuanta == 0) {
                int originBase = automaton.indexOfCell(originCell);
                int decoyBase = automaton.indexOfCell(decoy);
                frequency[originBase] += offset;
                frequency[decoyBase] -= offset;
            }
            AnisotropyProbe.seedPacket(automaton, originCell, packetQuanta);

            FccNeighborhood neighborhood = new FccNeighborhood(automaton);
            ContactPredicate predicate = new ContactPredicate(new MemberGeometry(ContactAtlasGenerator.GEOMETRY_RESOLUTION,
                                                                                   ContactAtlasGenerator.RADIUS));
            ContactScan scan = new ContactScan(automaton, neighborhood,
                                                predicate);
            CollisionStatistics statistics = new CollisionStatistics();
            CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                       new QuantaExchangeRule(),
                                                       statistics);
            HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);
            ConservationAudit audit = new ConservationAudit(automaton);
            return new SubstrateFactory.Substrate(automaton,
                                                   new AuditedRun(hybrid,
                                                                   audit),
                                                   statistics);
        };
    }

    /** {@code proj_d(cell - origin)}, matching the estimator's convention. */
    private static double projection(int dx, int dy, int dz, Direction d) {
        return switch (d) {
            case X100 -> dx;
            case X110 -> (dx + dy) / Math.sqrt(2);
            case X111 -> (dx + dy + dz) / Math.sqrt(3);
        };
    }

    /**
     * {@code M_d(field)} -- the mass-weighted second moment about
     * {@code origin} along {@code d}. Recomputed here rather than reached
     * for in production ({@code AnisotropyProbe.secondMoment} is private):
     * validated against {@link AnisotropyProbe#uniformSecondMoment} on the
     * all-ones field by {@link #secondMomentReplicaAgreesWithProduction()}.
     */
    private static double secondMomentOf(double[] field, Point3i extent,
                                           Point3i origin, Direction d) {
        double mass = 0;
        double weighted = 0;
        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    double m = Math.abs(field[index(extent, i, j, k)]);
                    if (m == 0.0) {
                        continue;
                    }
                    double p = projection(i - origin.x, j - origin.y,
                                           k - origin.z, d);
                    mass += m;
                    weighted += m * p * p;
                }
            }
        }
        return mass > 0 ? weighted / mass : 0.0;
    }

    /** {@code sat_d = M_d / M_d^unif} for ONE direction (production takes the max). */
    private static double saturationAlong(double[] field, Point3i extent,
                                            Point3i origin, Direction d) {
        return secondMomentOf(field, extent, origin, d)
               / AnisotropyProbe.uniformSecondMoment(extent, origin, d);
    }

    /**
     * The L1 mass fraction sitting at {@code proj_d^2 >= threshold} -- the
     * left-hand side of the Markov bound {@link
     * #momentSaturationBoundsDelocalizedMassDistributionFree} checks.
     */
    private static double massFractionBeyond(double[] field, Point3i extent,
                                               Point3i origin, Direction d,
                                               double threshold) {
        double total = 0;
        double beyond = 0;
        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    double m = Math.abs(field[index(extent, i, j, k)]);
                    if (m == 0.0) {
                        continue;
                    }
                    double p = projection(i - origin.x, j - origin.y,
                                           k - origin.z, d);
                    total += m;
                    if (p * p >= threshold) {
                        beyond += m;
                    }
                }
            }
        }
        return total > 0 ? beyond / total : 0.0;
    }

    /**
     * A separable wrapped (periodic-image-summed) Gaussian on the lattice
     * with per-axis variances -- the 3-D generalization of {@link
     * #wrappedMoment}'s 1-D family, used to exhibit W2's anisotropic blind
     * spot dynamically.
     */
    private static double[] wrappedGaussian3d(Point3i extent, Point3i origin,
                                                double sx, double sy,
                                                double sz) {
        double[] wx = wrappedAxisWeights(extent.x, origin.x, sx);
        double[] wy = wrappedAxisWeights(extent.y, origin.y, sy);
        double[] wz = wrappedAxisWeights(extent.z, origin.z, sz);
        double[] field = new double[extent.x * extent.y * extent.z];
        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    field[index(extent, i, j, k)] = wx[i] * wy[j] * wz[k];
                }
            }
        }
        return field;
    }

    private static double[] wrappedAxisWeights(int l, int o, double s) {
        double[] w = new double[l];
        for (int i = 0; i < l; i++) {
            double sum = 0;
            for (int n = -4; n <= 4; n++) {
                double d = i - o + n * (double) l;
                sum += Math.exp(-(d * d) / (2 * s));
            }
            w[i] = sum;
        }
        return w;
    }

    /**
     * Drives the wrapped-Gaussian family with per-axis variance RATES
     * {@code (rx,ry,rz)} until its measured moment saturation equals
     * {@code targetSaturation}, then returns the fractional understatement
     * of the instantaneous {@code X111} rate there: {@code 1 - measured /
     * (rx+ry+rz)/3}, the free-diffusion value being the unwrapped
     * {@code X111} rate of a separable field with those per-axis rates.
     *
     * <p>Bisecting on SATURATION rather than on time is the point: it is
     * what lets an isotropic and an anisotropic response be compared at the
     * saturation W2 would report for each, which is the comparison the
     * tolerance's claim is about.
     */
    private static double x111RateUnderstatementAtSaturation(Point3i extent,
                                                               Point3i origin,
                                                               double rx,
                                                               double ry,
                                                               double rz,
                                                               double targetSaturation) {
        double s0 = 0.25;
        double lo = 1e-4;
        double hi = 50.0;
        for (int it = 0; it < 200; it++) {
            double mid = 0.5 * (lo + hi);
            double[] field = wrappedGaussian3d(extent, origin, s0 + rx * mid,
                                                 s0 + ry * mid, s0 + rz * mid);
            if (AnisotropyProbe.momentSaturation(field, extent, origin)
                < targetSaturation) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double t = 0.5 * (lo + hi);
        double h = 1e-5 * t;
        double[] plus = wrappedGaussian3d(extent, origin, s0 + rx * (t + h),
                                            s0 + ry * (t + h),
                                            s0 + rz * (t + h));
        double[] minus = wrappedGaussian3d(extent, origin, s0 + rx * (t - h),
                                             s0 + ry * (t - h),
                                             s0 + rz * (t - h));
        double measured = (secondMomentOf(plus, extent, origin, Direction.X111)
                           - secondMomentOf(minus, extent, origin,
                                             Direction.X111))
                          / (2 * h);
        return 1.0 - measured / ((rx + ry + rz) / 3.0);
    }

    /**
     * The wrapped-Gaussian calibration behind {@link
     * AnisotropyProbe#RESPONSE_MOMENT_SATURATION_TOLERANCE}. Places a
     * periodic-image-summed Gaussian of variance {@code s} on {@code l}
     * sites about a centered origin -- the exact 1-D problem the {@code
     * X100} projection solves -- finds the {@code s} at which the moment
     * saturation ratio {@code M(s)/M^unif} equals {@code targetSaturation},
     * and returns the fractional understatement of the instantaneous rate
     * there, {@code 1 - dM/ds} (free diffusion would give {@code dM/ds = 1}
     * exactly).
     */
    private static double wrappedGaussianSlopeUnderstatement(int l,
                                                               double targetSaturation) {
        // Integer origin at l/2, matching nearestEvenParityCenter's
        // convention on the x axis.
        double uniform = 0;
        for (int i = 0; i < l; i++) {
            uniform += (double) (i - l / 2) * (i - l / 2);
        }
        uniform /= l;

        double lo = 1e-4;
        double hi = 100.0 * uniform;
        for (int it = 0; it < 200; it++) {
            double mid = 0.5 * (lo + hi);
            if (wrappedMoment(l, mid) / uniform < targetSaturation) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double s = 0.5 * (lo + hi);
        double h = s * 1e-5;
        double slope = (wrappedMoment(l, s + h) - wrappedMoment(l, s - h))
                        / (2 * h);
        return 1.0 - slope;
    }

    /**
     * Second moment about the centered origin of a wrapped Gaussian of
     * variance {@code s} on {@code l} periodic sites. Images out to
     * {@code +/-4L} are summed, which is far past numerical relevance for
     * every {@code s} this calibration probes.
     */
    private static double wrappedMoment(int l, double s) {
        int o = l / 2;
        double num = 0;
        double den = 0;
        for (int i = 0; i < l; i++) {
            double w = 0;
            for (int n = -4; n <= 4; n++) {
                double d = i - o + n * (double) l;
                w += Math.exp(-(d * d) / (2 * s));
            }
            num += w * (double) (i - o) * (i - o);
            den += w;
        }
        return num / den;
    }

    /**
     * A fully-supported uniform SIGNED field -- every cell nonzero, both
     * signs present, magnitudes comparable. Stands in for the §D-A signed
     * background at the field level (the exact per-member draw is
     * {@code PhaseCMeasurement.sharedInitialCondition}'s business, not
     * this estimator's).
     */
    private static double[] uniformSignedField(Point3i extent) {
        double[] field = new double[extent.x * extent.y * extent.z];
        Random random = new Random(90210L);
        for (int i = 0; i < field.length; i++) {
            int q = random.nextInt(9) - 4; // [-4,+4]
            field[i] = q == 0 ? 3 : q;     // full support: never zero
        }
        return field;
    }

    /**
     * The same separable Gaussian as {@link #gaussianPacketField}, but
     * NORMALIZED to unit total at every tick, so the packet's own quanta
     * total is tick-invariant exactly as a conserved automaton packet's
     * would be -- which is what invariant (I1) requires of a matched
     * pair's difference field. Per-axis variance is
     * {@code sigma0Sq + t*rateAxis}, so {@code M_[100](t)} has slope
     * exactly {@code rateX}.
     */
    private static double[][] normalizedGaussianPacketField(Point3i extent,
                                                              Point3i origin,
                                                              int ticks,
                                                              double sigma0Sq,
                                                              double rateX,
                                                              double rateY,
                                                              double rateZ) {
        double[][] fieldByTick = gaussianPacketField(extent, origin, ticks,
                                                      sigma0Sq, rateX, rateY,
                                                      rateZ, 42L, 0.0);
        for (double[] field : fieldByTick) {
            double total = 0;
            for (double v : field) {
                total += v;
            }
            for (int i = 0; i < field.length; i++) {
                field[i] /= total;
            }
        }
        return fieldByTick;
    }

    /**
     * A large, signed, tick-DRIFTING background: a deterministic signed
     * base per cell whose magnitude is progressively tilted along the x
     * axis as ticks advance. Two properties matter for {@link
     * #matchedPairRecoversPacketSpreadThroughASignedDriftingBackground}:
     * its L1 mass dwarfs the unit-normalized packet's (so a raw
     * single-field estimate is diluted by orders of magnitude), and its
     * own spatial second moment has a direction-dependent trend (so a raw
     * estimate is corrupted in the RATIO too, not merely scaled down). It
     * is a pure function of (cell, tick), so both halves of the matched
     * pair receive bit-identical values and it cancels exactly.
     */
    private static double[][] signedDriftingBackground(Point3i extent,
                                                         Point3i origin,
                                                         int ticks, long seed) {
        Random random = new Random(seed);
        int cells = extent.x * extent.y * extent.z;
        double[] base = new double[cells];
        for (int i = 0; i < cells; i++) {
            int q = random.nextInt(9) - 4;
            base[i] = q == 0 ? 3 : q;
        }
        double[][] byTick = new double[ticks][cells];
        for (int t = 0; t < ticks; t++) {
            for (int i = 0; i < extent.x; i++) {
                double tilt = 1.0
                              + 0.5 * t * Math.abs(i - origin.x) / extent.x;
                for (int j = 0; j < extent.y; j++) {
                    for (int k = 0; k < extent.z; k++) {
                        int idx = index(extent, i, j, k);
                        byTick[t][idx] = base[idx] * tilt;
                    }
                }
            }
        }
        return byTick;
    }

    private static Direction maxDirection(EstimatorResult result) {
        Direction best = null;
        double bestMagnitude = Double.NEGATIVE_INFINITY;
        for (DirectionMagnitude dm : result.perDirection().values()) {
            if (dm.magnitude() > bestMagnitude) {
                bestMagnitude = dm.magnitude();
                best = dm.direction();
            }
        }
        return best;
    }

    private static Direction minDirection(EstimatorResult result) {
        Direction best = null;
        double bestMagnitude = Double.POSITIVE_INFINITY;
        for (DirectionMagnitude dm : result.perDirection().values()) {
            if (dm.magnitude() < bestMagnitude) {
                bestMagnitude = dm.magnitude();
                best = dm.direction();
            }
        }
        return best;
    }
}
