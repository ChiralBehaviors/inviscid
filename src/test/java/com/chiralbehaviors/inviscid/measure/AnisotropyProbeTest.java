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
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.measure.AnisotropyProbe.BootstrapCi;
import com.chiralbehaviors.inviscid.measure.AnisotropyProbe.DirectionMagnitude;
import com.chiralbehaviors.inviscid.measure.AnisotropyProbe.EstimatorResult;
import com.chiralbehaviors.inviscid.measure.StructureFactor.Direction;

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
     * {@link com.chiralbehaviors.inviscid.lga.HybridAutomaton},
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
     * generously so the boundary-adjacent mass fraction is far below
     * {@link AnisotropyProbe#WRAP_SAFETY_MASS_TOLERANCE} - {@link
     * AnisotropyProbe#transportEstimate} would otherwise (correctly)
     * reject this fixture.
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
