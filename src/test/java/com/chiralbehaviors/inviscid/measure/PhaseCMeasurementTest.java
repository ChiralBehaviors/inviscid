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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.Consumer;

import javax.vecmath.Point3i;

import org.junit.BeforeClass;
import org.junit.Test;

import com.chiralbehaviors.inviscid.QuantaField;

/**
 * C.5 (bead inviscid-0nx.22): the five named acceptance tests, validating
 * the COMMITTED {@code measurement-report-phaseC.tsv} artifact -- mirrors
 * the established precedent ({@code AnisotropyProbeTest
 * .committedPhaseAArtifactIsWellFormed}, {@code CommittedContactAtlasTest}):
 * the campaign that generates this artifact is expensive (multiple LGA/
 * hybrid runs at real tick counts) and is driven by {@link
 * PhaseCMeasurement#main(String[])}, NOT by surefire -- these tests read
 * and numerically/structurally validate what was committed.
 *
 * @author halhildebrand
 */
public class PhaseCMeasurementTest {

    private static final String RESOURCE_PATH = "lga/measurement-report-phaseC.tsv";

    private static List<String> LINES;

    @BeforeClass
    public static void loadArtifact() throws IOException {
        URL resource = PhaseCMeasurementTest.class.getClassLoader()
                                                     .getResource(RESOURCE_PATH);
        if (resource == null) {
            fail("regenerate with PhaseCMeasurement.main() and review the diff: "
               + "committed report src/test/resources/lga/measurement-report-phaseC.tsv is missing");
        }
        Path path = Paths.get(resource.getPath());
        LINES = Files.readAllLines(path);
    }

    // ------------------------------------------------------------------
    // Test 1: lgaSpectrumDiffersFromK0Baseline.
    // ------------------------------------------------------------------

    /**
     * The LGA's per-member power spectrum must be MEASURABLY BROADER
     * (lower peak-bin concentration) than the committed K=0 golden
     * baseline's regime -- a stated effect size, not eyeballed. If the
     * collisions are not doing physics, this is the finding, and the
     * assertion below fails loudly rather than silently passing on a
     * near-zero effect.
     */
    @Test
    public void lgaSpectrumDiffersFromK0Baseline() {
        double[] k0 = spectralSummaryRow("K0");
        double[] lga = spectralSummaryRow("LGA");
        double k0MeanPeakFraction = k0[0];
        double lgaMeanPeakFraction = lga[0];
        double drop = k0MeanPeakFraction - lgaMeanPeakFraction;

        assertTrue("K0 baseline mean peak fraction should be a near-pure-tone concentration (>0.9), was "
                   + k0MeanPeakFraction, k0MeanPeakFraction > 0.9);
        assertTrue("expected the LGA's collision-broadened spectrum to have a LOWER mean peak-bin "
                   + "concentration than the K0 pure-tone baseline (K0=" + k0MeanPeakFraction
                   + ", LGA=" + lgaMeanPeakFraction + ", drop=" + drop
                   + ") -- a stated, non-trivial effect size (> 1e-3), not eyeballed",
                   drop > 1e-3);

        boolean hasVerdict = LINES.stream()
                                   .anyMatch(l -> l.startsWith("# spectralProgressionVerdict="));
        assertTrue("expected the spectralProgressionVerdict provenance line stating the three-way progression finding",
                   hasVerdict);

        // User decision D3 (pre-registered): BOTH absolute (rad/tick) and
        // fractional linewidths reported, neither promoted to sole
        // headline -- verify both columns are present and non-degenerate
        // for the collision-bearing substrates.
        boolean hasLinewidthDefinition = LINES.stream()
                                                .anyMatch(l -> l.startsWith("# linewidthDefinition="));
        assertTrue("expected the D3 dual-linewidth definition provenance line",
                   hasLinewidthDefinition);
        double lgaAbsoluteLinewidth = lga[4];
        double lgaFractionalLinewidth = lga[5];
        double k0AbsoluteLinewidth = k0[4];
        assertTrue("expected the LGA's ABSOLUTE linewidth (rad/tick) to be measurably greater than the K0 baseline's ("
                   + "LGA=" + lgaAbsoluteLinewidth + ", K0=" + k0AbsoluteLinewidth + ")",
                   lgaAbsoluteLinewidth > k0AbsoluteLinewidth);
        assertTrue("expected a defined (non-NaN) mean FRACTIONAL linewidth for the LGA",
                   !Double.isNaN(lgaFractionalLinewidth));
    }

    // ------------------------------------------------------------------
    // Test 2: structureFactorShowsAtLeastOneRidge (or asserted absence).
    // ------------------------------------------------------------------

    /**
     * The report must state, unambiguously, whether the LGA's dynamic
     * structure factor shows a propagating ridge (nonzero SPECTRAL ridge
     * slope in at least one direction) or its absence -- and the verdict
     * line must be CONSISTENT with the underlying per-direction data, not
     * an unchecked assertion floating independently of it.
     */
    @Test
    public void structureFactorShowsAtLeastOneRidge() {
        String verdictLine = LINES.stream()
                                   .filter(l -> l.startsWith("# structureFactorRidgeVerdict="))
                                   .findFirst()
                                   .orElse(null);
        assertTrue("expected exactly one structureFactorRidgeVerdict provenance line",
                   verdictLine != null);
        boolean claimsPresent = verdictLine.contains("RIDGE PRESENT");
        boolean claimsAbsent = verdictLine.contains("RIDGE ABSENT");
        assertTrue("verdict must claim EXACTLY ONE of RIDGE PRESENT / RIDGE ABSENT -- an ambiguous result must fail this test",
                   claimsPresent ^ claimsAbsent);

        List<String> spectralDirectionRows = LINES.stream()
                                                     .filter(l -> l.startsWith("LGA_POOLED_DIRECTION\tSPECTRAL\t"))
                                                     .toList();
        assertEquals("expected 3 SPECTRAL pooled-direction rows (X100/X110/X111)",
                     3, spectralDirectionRows.size());
        boolean anyNonzero = false;
        for (String row : spectralDirectionRows) {
            String[] parts = row.split("\t");
            double mean = Double.parseDouble(parts[3]);
            if (Math.abs(mean) > 1e-12) {
                anyNonzero = true;
            }
        }
        assertEquals("verdict's RIDGE PRESENT/ABSENT claim must match whether any pooled SPECTRAL direction actually has a nonzero ridge slope",
                     claimsPresent, anyNonzero);
    }

    // ------------------------------------------------------------------
    // Test 3: conservationHoldsAcrossTheMeasurementRun (exact).
    // ------------------------------------------------------------------

    @Test
    public void conservationHoldsAcrossTheMeasurementRun() {
        long ticksAudited = headerLong("conservationTicksAudited");
        long violations = headerLong("conservationViolations");
        // Post-C2-fix (Critical): the exact expected total, DERIVED from
        // the campaign's own shape constants (never a bare magic number)
        // -- anisotropy (8 seeds x 127 ticks) + K0 ((256-1)*225) +
        // hybrid/LGA spectral (2 x (32-1)*225) + long-run hybrid/LGA
        // (2 x 2000). This must equal the sum of every driver's
        // ConservationAudit#ledger() size (production reads the ledger,
        // never hand-derives this count) -- cross-checks to 76341.
        long expected = expectedConservationTicksAudited();
        assertEquals("76341", 76341L, expected);
        assertEquals("expected conservationTicksAudited to be the EXACT sum of every driver's "
                     + "ConservationAudit#ledger().size() (STRICT audits every driver, every "
                     + "sub-measurement -- the header's own claim), derived here from the campaign's "
                     + "own shape constants rather than assumed",
                     expected, ticksAudited);
        assertEquals("expected EXACTLY zero conservation violations across the entire measurement campaign",
                     0L, violations);
        boolean strictMode = LINES.stream()
                                   .anyMatch(l -> l.startsWith("# conservationMode=STRICT"));
        assertTrue("expected the report to state strict-mode auditing (every driver, every sub-measurement)",
                   strictMode);
    }

    /**
     * The expected {@code conservationTicksAudited} total, derived from the
     * campaign's own shape constants -- NOT a hand-derived assumption about
     * any driver's internal loop structure (that was C2's bug). Every term
     * here is the number of {@code auditTick} calls a STRICT
     * {@link com.chiralbehaviors.inviscid.measure.ConservationAudit} makes
     * for that driver, matching {@code recordStridedPhaseSeries}'s
     * {@code (samples-1)*stride} advance count for the spectral drivers,
     * {@code AnisotropyProbe}'s {@code (ticks-1)} per-seed tick count, and
     * the long run's {@code LONG_RUN_TICKS} per substrate.
     */
    private static long expectedConservationTicksAudited() {
        long anisotropy = (long) AnisotropyProbe.DEFAULT_SEEDS.length
                           * (AnisotropyProbe.DEFAULT_TICKS - 1);
        long k0 = (long) (BaselineSpectrumHarness.FFT_LENGTH - 1)
                  * BaselineSpectrumHarness.STRIDE;
        long spectralBearing = 2L * (PhaseCMeasurement.SPECTRAL_FFT_LENGTH - 1)
                                * PhaseCMeasurement.SPECTRAL_STRIDE;
        long longRun = 2L * PhaseCMeasurement.LONG_RUN_TICKS;
        return anisotropy + k0 + spectralBearing + longRun;
    }

    // ------------------------------------------------------------------
    // Test 4: anisotropyIsReportedForBothEstimators.
    // ------------------------------------------------------------------

    @Test
    public void anisotropyIsReportedForBothEstimators() {
        List<String> summaryRows = LINES.stream()
                                          .filter(l -> l.startsWith("LGA_POOLED_SUMMARY\t"))
                                          .toList();
        assertEquals("expected exactly 2 LGA_POOLED_SUMMARY rows: TRANSPORT and SPECTRAL, neither silently dropped",
                     2, summaryRows.size());
        boolean hasTransport = summaryRows.stream()
                                            .anyMatch(l -> l.startsWith("LGA_POOLED_SUMMARY\tTRANSPORT\t"));
        boolean hasSpectral = summaryRows.stream()
                                           .anyMatch(l -> l.startsWith("LGA_POOLED_SUMMARY\tSPECTRAL\t"));
        assertTrue("expected a TRANSPORT pooled-summary row", hasTransport);
        assertTrue("expected a SPECTRAL pooled-summary row", hasSpectral);

        List<String> directionRows = LINES.stream()
                                            .filter(l -> l.startsWith("LGA_POOLED_DIRECTION\t"))
                                            .toList();
        assertEquals("expected 2 estimators x 3 directions = 6 LGA_POOLED_DIRECTION rows",
                     6, directionRows.size());

        boolean hasSmallNFlag = LINES.stream()
                                       .anyMatch(l -> l.startsWith("# lgaSmallNEarlyTimeFlag="));
        assertTrue("expected the small-N/early-time framing header, same non-fabrication discipline as Phase A",
                   hasSmallNFlag);

        // Post-critique fix (Critical): per-seed per-direction diagnostic
        // rows, matching Phase A's DIRECTION/SUMMARY format, must be
        // present -- this is the data the winning-direction-stability
        // diagnostic (and any future re-analysis) depends on.
        long perSeedDirectionRows = LINES.stream()
                                           .filter(l -> l.startsWith("LGA_DIRECTION\t"))
                                           .count();
        assertEquals("expected 8 seeds x 2 estimators x 3 directions = 48 LGA_DIRECTION rows",
                     48, perSeedDirectionRows);
        long naiveSummaryRows = LINES.stream()
                                       .filter(l -> l.startsWith("LGA_SUMMARY\t"))
                                       .count();
        assertEquals("expected 2 LGA_SUMMARY rows (naive per-seed diagnostic, TRANSPORT + SPECTRAL)",
                     2, naiveSummaryRows);

        // Winning-direction-stability diagnostic: TRANSPORT only (SPECTRAL
        // is degenerate -- see winningDirectionSpectralNote).
        long winningDirectionRows = LINES.stream()
                                           .filter(l -> l.startsWith("WINNING_DIRECTION\tTRANSPORT\t"))
                                           .count();
        assertEquals("expected 8 per-seed WINNING_DIRECTION rows for TRANSPORT",
                     8, winningDirectionRows);
        String stabilityRow = LINES.stream()
                                     .filter(l -> l.startsWith("WINNING_DIRECTION_STABILITY\tTRANSPORT\t"))
                                     .findFirst()
                                     .orElse(null);
        assertTrue("expected a WINNING_DIRECTION_STABILITY row for TRANSPORT",
                   stabilityRow != null);
        String[] stabilityParts = stabilityRow.split("\t");
        int modeCount = Integer.parseInt(stabilityParts[3]);
        int totalSeeds = Integer.parseInt(stabilityParts[4]);
        assertEquals("expected totalSeeds=8 in the stability row", 8, totalSeeds);
        assertTrue("expected modeCount in [1,8]", modeCount >= 1 && modeCount <= 8);
        boolean hasSpectralNote = LINES.stream()
                                         .anyMatch(l -> l.startsWith("# winningDirectionSpectralNote="));
        assertTrue("expected the SPECTRAL-degenerate explanatory note (not a silently-omitted diagnostic)",
                   hasSpectralNote);

        // CI-vs-permutation reconciliation and power recommendation
        // (relay batch items 2 and 4) must be present and non-empty.
        boolean hasReconciliation = LINES.stream()
                                           .anyMatch(l -> l.startsWith("# ciVsPermutationReconciliation=")
                                                          && l.contains("Phase A")
                                                          && l.contains("Phase C"));
        assertTrue("expected the CI-vs-permutation two-campaign reconciliation line naming both campaigns",
                   hasReconciliation);
        boolean hasPowerRecommendation = LINES.stream()
                                                .anyMatch(l -> l.startsWith("# powerRecommendationForGate="));
        assertTrue("expected the power-recommendation-for-.23-gate line",
                   hasPowerRecommendation);
    }

    // ------------------------------------------------------------------
    // Test 5: reportArtifactIsCompleteAndProvenanced.
    // ------------------------------------------------------------------

    @Test
    public void reportArtifactIsCompleteAndProvenanced() {
        assertHeaderPresent("# bead=inviscid-0nx.22");
        assertHeaderStartsWithNonEmptyValue("# gitCommit=");
        assertHeaderStartsWithNonEmptyValue("# atlasVersion=");
        assertHeaderStartsWithNonEmptyValue("# nLga=");
        assertHeaderStartsWithNonEmptyValue("# subBinSteps=");
        assertHeaderStartsWithNonEmptyValue("# anisotropyExtent=");
        assertHeaderStartsWithNonEmptyValue("# anisotropySeeds=");
        assertHeaderStartsWithNonEmptyValue("# anisotropyTicks=");
        assertHeaderStartsWithNonEmptyValue("# spectralExtent=");
        assertHeaderStartsWithNonEmptyValue("# spectralSeed=");
        assertHeaderStartsWithNonEmptyValue("# spectralStride=");
        assertHeaderStartsWithNonEmptyValue("# spectralFftLength=");
        assertHeaderStartsWithNonEmptyValue("# nyquistQuantaBound=");
        assertHeaderStartsWithNonEmptyValue("# collisionOpportunitiesPerRevolutionAtLongRunQuantaBound=");
        assertHeaderStartsWithNonEmptyValue("# longRunExtent=");
        assertHeaderStartsWithNonEmptyValue("# longRunSeed=");
        assertHeaderStartsWithNonEmptyValue("# longRunTicks=");

        boolean hasSeeds = LINES.stream()
                                 .anyMatch(l -> l.startsWith("# anisotropySeeds=42,43,44,45,46,47,48,49"));
        assertTrue("expected the literal 8-seed list 42..49, IDENTICAL to Phase A", hasSeeds);

        boolean hasEscalation = LINES.stream()
                                       .anyMatch(l -> l.startsWith("# ESCALATION="));
        assertTrue("expected the mandatory ESCALATION line -- no posture selected",
                   hasEscalation);
        for (String posture : List.of("# posture(i)_acceptAndCharacterize_evidence=",
                                       "# posture(ii)_fchcProjection_evidence=",
                                       "# posture(iii)_orientationalStateRestoresIsotropy_evidence=")) {
            assertHeaderStartsWithNonEmptyValue(posture);
        }
        // Explicitly must NOT contain a posture-selection statement.
        boolean selectsAPosture = LINES.stream()
                                         .anyMatch(l -> l.toLowerCase(Locale.ROOT)
                                                          .contains("selected posture")
                                                        || l.toLowerCase(Locale.ROOT)
                                                            .contains("verdict: posture"));
        assertTrue("the report must NOT select an isotropy posture", !selectsAPosture);

        long fieldComparisonRows = LINES.stream()
                                          .filter(l -> l.startsWith("COLLISION_FIELD\t"))
                                          .count();
        assertEquals("expected all 5 CollisionStatistics fields reported (totalCollisions, "
                     + "effectiveCollisions, collisionsPerDirection, transferMagnitudeHistogram, meanFreePathProxy)",
                     5, fieldComparisonRows);

        long histogramRows = LINES.stream()
                                    .filter(l -> l.startsWith("QUANTA_HISTOGRAM\t"))
                                    .count();
        assertEquals("expected before/after quanta-distribution histogram rows",
                     2, histogramRows);

        // Post-critique fix (Significant #2/#3, relay items 3 and 5): full
        // 4-window series in the verdict lines (not first-vs-last), and
        // the anisotropy campaign's own within-window homogeneity data.
        String gapVerdict = LINES.stream()
                                   .filter(l -> l.startsWith("# gapMechanismVerdict="))
                                   .findFirst()
                                   .orElse(null);
        assertTrue("expected a gapMechanismVerdict line", gapVerdict != null);
        assertTrue("expected the FULL relative-rate-gap series (not just first/last) in gapMechanismVerdict",
                   gapVerdict.contains("relativeRateGapByWindow=[")
                   && gapVerdict.chars().filter(c -> c == ',').count() >= 3);
        String noOpVerdict = LINES.stream()
                                    .filter(l -> l.startsWith("# equalityConcentrationEmpiricalTrend="))
                                    .findFirst()
                                    .orElse(null);
        assertTrue("expected an equalityConcentrationEmpiricalTrend line", noOpVerdict != null);
        assertTrue("expected the FULL no-op-fraction series (not just first/last) plus an explicit one-window-lag statement",
                   noOpVerdict.contains("lgaNoOpFractionByWindow=[")
                   && (noOpVerdict.contains("ONE WINDOW LATER")
                       || noOpVerdict.contains("different window")));

        // Post-C1-fix (Critical): the lag DIRECTION must match the two
        // series' own argmax windows, computed HERE independently of the
        // production code (not a tautological contains()-over-two-strings
        // check) -- this is exactly the check that would have caught the
        // inverted subject/order bug (guard noOpPeakWindow==gapPeakWindow+1
        // means the NO-OP trend peaks LATER, but the pre-fix prose made the
        // GAP series the subject and called it "later").
        int computedGapPeak = argmaxOf(parseSeries(gapVerdict, "relativeRateGapByWindow="));
        int computedNoOpPeak = argmaxOf(parseSeries(noOpVerdict, "lgaNoOpFractionByWindow="));
        if (computedNoOpPeak == computedGapPeak + 1) {
            assertTrue("expected the artifact to state the NO-OP trend's own peak (window "
                       + computedNoOpPeak + ") occurs ONE WINDOW LATER than Section 3's gap-trend "
                       + "peak (window " + computedGapPeak + ") -- subject must be the no-op trend, "
                       + "not the reverse. Got: " + noOpVerdict,
                       noOpVerdict.contains("(window " + computedNoOpPeak
                                             + ") occurs exactly ONE WINDOW LATER than")
                       && noOpVerdict.contains("(window " + computedGapPeak + ")"));
        } else {
            assertTrue("expected the artifact to state the NO-OP trend's own peak (window "
                       + computedNoOpPeak + ") occurs at a DIFFERENT window than Section 3's "
                       + "gap-trend peak (window " + computedGapPeak + "). Got: " + noOpVerdict,
                       noOpVerdict.contains("(window " + computedNoOpPeak + ") occurs")
                       && noOpVerdict.contains("at a different window than")
                       && noOpVerdict.contains("(window " + computedGapPeak + ")"));
        }

        long anisotropyWithinWindowRows = LINES.stream()
                                                  .filter(l -> l.startsWith("ANISOTROPY_WITHIN_WINDOW\t"))
                                                  .count();
        assertEquals("expected 4 ANISOTROPY_WITHIN_WINDOW quartile rows",
                     4, anisotropyWithinWindowRows);
        boolean hasHomogeneityVerdict = LINES.stream()
                                               .anyMatch(l -> l.startsWith("# anisotropyWithinWindowHomogeneityVerdict="));
        assertTrue("expected the anisotropy campaign's own within-window homogeneity verdict",
                   hasHomogeneityVerdict);

        boolean hasSdb = LINES.stream()
                               .anyMatch(l -> l.startsWith("# sdbClosedForm="));
        assertTrue("expected the semi-detailed-balance closed-form parity-floor line",
                   hasSdb);
        boolean hasEqualityTrend = LINES.stream()
                                          .anyMatch(l -> l.startsWith("# equalityConcentrationEmpiricalTrend="));
        assertTrue("expected the equality-concentration empirical trend line",
                   hasEqualityTrend);
    }

    // ------------------------------------------------------------------
    // Fix-round Critical 1 (Section 4B: ANISOTROPY_CAMPAIGN_EQUALITY).
    // ------------------------------------------------------------------

    @Test
    public void anisotropyCampaignEqualityIsCharacterizedOnItsOwnCampaign() {
        long perSeedRows = LINES.stream()
                                  .filter(l -> l.startsWith("ANISOTROPY_CAMPAIGN_EQUALITY_SEED\t"))
                                  .count();
        assertEquals("expected 8 per-seed ANISOTROPY_CAMPAIGN_EQUALITY_SEED rows",
                     8, perSeedRows);
        long pooledRows = LINES.stream()
                                 .filter(l -> l.startsWith("ANISOTROPY_CAMPAIGN_EQUALITY_POOLED\t"))
                                 .count();
        assertEquals("expected exactly 1 pooled ANISOTROPY_CAMPAIGN_EQUALITY_POOLED row",
                     1, pooledRows);

        String characterization = LINES.stream()
                                         .filter(l -> l.startsWith("# anisotropyCampaignEqualityCharacterization="))
                                         .findFirst()
                                         .orElse(null);
        assertTrue("expected the anisotropyCampaignEqualityCharacterization line",
                   characterization != null);
        assertTrue("expected the characterization to state the absorbing-initial-condition finding",
                   characterization.contains("ABSORBING")
                   && characterization.contains("INITIAL CONDITION"));

        // Cross-check: every seed's excitedMembersInitial must be exactly
        // AnisotropyProbe.DEFAULT_PACKET_QUANTA's cell size (30 members,
        // the packet's own cell) -- the campaign's background is
        // all-zero, so ONLY the packet cell is excited at t=0, for every
        // seed, deterministically.
        List<String> perSeedLines = LINES.stream()
                                           .filter(l -> l.startsWith("ANISOTROPY_CAMPAIGN_EQUALITY_SEED\t"))
                                           .toList();
        for (String line : perSeedLines) {
            String[] parts = line.split("\t");
            // recordType seed totalCollisions effectiveCollisions noOpFraction
            // excitedMembersInitial cellsTouchedInitial excitedMembersFinal
            // cellsTouchedFinal collisionsAtFirstTick collisionsAtLastTick
            long excitedInitial = Long.parseLong(parts[5]);
            int cellsTouchedInitial = Integer.parseInt(parts[6]);
            long excitedFinal = Long.parseLong(parts[7]);
            assertEquals("expected excitedMembersInitial==30 (exactly one packet cell) for line: "
                         + line, 30L, excitedInitial);
            assertEquals("expected cellsTouchedInitial==1 (exactly one packet cell) for line: "
                         + line, 1, cellsTouchedInitial);
            assertTrue("expected excitedMembersFinal >= excitedMembersInitial (quanta can only "
                       + "spread from the packet, never vanish) for line: " + line,
                       excitedFinal >= excitedInitial);
        }
    }

    /**
     * R2-1: the pooled ANISOTROPY_CAMPAIGN_EQUALITY_POOLED denominator
     * must be scoped to the WHOLE campaign (8 seeds' worth of even-parity
     * member slots), not one lattice -- the numerator is already summed
     * across all 8 seeds. Pins the corrected value (61440 = 256 * 30 * 8)
     * and the corrected narrative fraction (240/61440), replacing the
     * pre-fix 240/7680 reading that understated the contamination 8x.
     */
    @Test
    public void anisotropyCampaignEqualityPooledDenominatorSpansAllSeeds() {
        String pooledRow = LINES.stream()
                                  .filter(l -> l.startsWith("ANISOTROPY_CAMPAIGN_EQUALITY_POOLED\t"))
                                  .findFirst()
                                  .orElseThrow(() -> new AssertionError("missing ANISOTROPY_CAMPAIGN_EQUALITY_POOLED row"));
        String[] parts = pooledRow.split("\t");
        // recordType totalCollisions effectiveCollisions noOpFraction
        // activeMemberSlotsTotal evenParityCellsTotal
        // excitedMembersInitialTotal excitedMembersFinalTotal
        long activeMemberSlotsTotal = Long.parseLong(parts[4]);
        long evenParityCellsTotal = Long.parseLong(parts[5]);
        long excitedInitialTotal = Long.parseLong(parts[6]);
        assertEquals("expected evenParityCellsTotal to stay the SINGLE-lattice value (unaffected by "
                     + "the R2-1 fix, a different column)", 256L, evenParityCellsTotal);
        assertEquals("expected activeMemberSlotsTotal to be evenParityCellsTotal * 30 * 8 seeds "
                     + "(R2-1 fix)", evenParityCellsTotal * 30L * 8L, activeMemberSlotsTotal);
        assertEquals(61440L, activeMemberSlotsTotal);
        assertEquals("expected excitedMembersInitialTotal to remain 30*8=240 (unaffected numerator)",
                     240L, excitedInitialTotal);

        String characterization = LINES.stream()
                                         .filter(l -> l.startsWith("# anisotropyCampaignEqualityCharacterization="))
                                         .findFirst()
                                         .orElseThrow(() -> new AssertionError("missing characterization line"));
        assertTrue("expected the narrative to read 240/61440 at t=0 (corrected denominator)",
                   characterization.contains("240/61440"));
        assertTrue("expected the narrative to read the final count over 61440",
                   characterization.contains("/61440 after the campaign"));
        assertFalse("expected the pre-fix single-lattice 240/7680 reading to be gone",
                    characterization.contains("240/7680"));
    }

    /**
     * R2-2: "collisionsPerTick is FLAT per seed" was an unconditional
     * literal even though seed 47 goes 23-&gt;21. The narrative must now
     * be data-branched via {@link PhaseCMeasurement#tickFlatnessStatement}.
     */
    @Test
    public void tickFlatnessStatementIsDataBranchedNotAHardcodedLiteral() {
        assertEquals("EXACTLY FLAT per seed (first-tick == last-tick for every seed)",
                     PhaseCMeasurement.tickFlatnessStatement(new long[] { 5, 5, 5 },
                                                              new long[] { 5, 5, 5 }));
        String notFlat = PhaseCMeasurement.tickFlatnessStatement(new long[] { 23, 31 },
                                                                   new long[] { 21, 32 });
        assertTrue(notFlat, notFlat.contains("NOT EXACTLY FLAT"));
        assertTrue(notFlat, notFlat.contains("max |first-tick - last-tick| delta across seeds = 2"));

        String characterization = LINES.stream()
                                         .filter(l -> l.startsWith("# anisotropyCampaignEqualityCharacterization="))
                                         .findFirst()
                                         .orElseThrow(() -> new AssertionError("missing characterization line"));
        assertFalse("expected the unconditional FLAT literal to be gone from the artifact",
                    characterization.contains("collisionsPerTick is FLAT per seed"));
        assertTrue("expected the artifact to state the measured flatness (or lack of it)",
                   characterization.contains("collisionsPerTick is EXACTLY FLAT per seed")
                   || characterization.contains("collisionsPerTick is NOT EXACTLY FLAT per seed"));
    }

    /**
     * Reviewer R2-1 / MINOR-5: 4B's per-seed totalCollisions/
     * effectiveCollisions derive from the SAME CollisionStatistics sink
     * as 1B's LGA_COLLISIONS rows -- 1B already cross-checks against it
     * (contactDirectionCensusIsEmittedForAnisotropyCampaignAndLongRun);
     * 4B did not. Added per critic MINOR-5.
     */
    @Test
    public void anisotropyCampaignEqualityReconcilesWithLgaCollisions() {
        Map<Long, long[]> lgaCollisionsTotals = new LinkedHashMap<>();
        for (String line : LINES) {
            if (line.startsWith("LGA_COLLISIONS\t")) {
                String[] parts = line.split("\t");
                lgaCollisionsTotals.put(Long.parseLong(parts[1]),
                                         new long[] { Long.parseLong(parts[2]),
                                                       Long.parseLong(parts[3]) });
            }
        }
        assertEquals(8, lgaCollisionsTotals.size());
        List<String> perSeedRows = LINES.stream()
                                          .filter(l -> l.startsWith("ANISOTROPY_CAMPAIGN_EQUALITY_SEED\t"))
                                          .toList();
        assertEquals(8, perSeedRows.size());
        for (String row : perSeedRows) {
            String[] parts = row.split("\t");
            long seed = Long.parseLong(parts[1]);
            long totalCollisions = Long.parseLong(parts[2]);
            long effectiveCollisions = Long.parseLong(parts[3]);
            long[] expected = lgaCollisionsTotals.get(seed);
            assertTrue("expected a LGA_COLLISIONS row for seed " + seed, expected != null);
            assertEquals("expected ANISOTROPY_CAMPAIGN_EQUALITY_SEED totalCollisions to reconcile "
                         + "EXACTLY with LGA_COLLISIONS for seed " + seed, expected[0], totalCollisions);
            assertEquals("expected ANISOTROPY_CAMPAIGN_EQUALITY_SEED effectiveCollisions to reconcile "
                         + "EXACTLY with LGA_COLLISIONS for seed " + seed, expected[1], effectiveCollisions);
        }
    }

    // ------------------------------------------------------------------
    // Fix-round Critical 2 (Section 1B: FCC contact-direction census).
    // ------------------------------------------------------------------

    @Test
    public void contactDirectionCensusIsEmittedForAnisotropyCampaignAndLongRun() {
        List<String> perSeedRows = LINES.stream()
                                          .filter(l -> l.startsWith("LGA_CONTACT_DIRECTION_SEED\t"))
                                          .toList();
        assertEquals("expected 8 per-seed LGA_CONTACT_DIRECTION_SEED rows",
                     8, perSeedRows.size());

        String pooledRow = LINES.stream()
                                  .filter(l -> l.startsWith("LGA_CONTACT_DIRECTION_POOLED\t"))
                                  .findFirst()
                                  .orElse(null);
        assertTrue("expected exactly one LGA_CONTACT_DIRECTION_POOLED row", pooledRow != null);
        String[] pooledParts = pooledRow.split("\t");
        // recordType d1 d2 d3 d4 d5 d6 total chiSquare df pValue
        assertEquals(11, pooledParts.length);
        int df = Integer.parseInt(pooledParts[9]);
        assertEquals("expected df=5 (6 directions - 1)", 5, df);
        double pValue = Double.parseDouble(pooledParts[10]);
        assertTrue("expected a valid p-value in [0,1], was " + pValue,
                   pValue >= 0.0 && pValue <= 1.0);

        // Cross-check (relay item 2's explicit instruction): per-direction
        // sums must reconcile EXACTLY with the committed LGA_COLLISIONS
        // totals, for every seed.
        java.util.Map<Long, Long> lgaCollisionsTotals = new java.util.HashMap<>();
        for (String line : LINES) {
            if (line.startsWith("LGA_COLLISIONS\t")) {
                String[] parts = line.split("\t");
                lgaCollisionsTotals.put(Long.parseLong(parts[1]), Long.parseLong(parts[2]));
            }
        }
        assertEquals(8, lgaCollisionsTotals.size());
        for (String row : perSeedRows) {
            String[] parts = row.split("\t");
            // recordType seed d1 d2 d3 d4 d5 d6 total
            long seed = Long.parseLong(parts[1]);
            long directionSum = 0;
            for (int i = 2; i <= 7; i++) {
                directionSum += Long.parseLong(parts[i]);
            }
            long reportedTotal = Long.parseLong(parts[8]);
            assertEquals("expected the row's own total column to equal the sum of its 6 direction "
                         + "columns for seed " + seed, directionSum, reportedTotal);
            Long expected = lgaCollisionsTotals.get(seed);
            assertTrue("expected a LGA_COLLISIONS row for seed " + seed, expected != null);
            assertEquals("expected per-direction sum to reconcile EXACTLY with the committed "
                         + "LGA_COLLISIONS total for seed " + seed, (long) expected, directionSum);
        }

        // Long-run census + hopping tensor, both substrates.
        List<String> longRunRows = LINES.stream()
                                          .filter(l -> l.startsWith("LONG_RUN_CONTACT_DIRECTION\t"))
                                          .toList();
        assertEquals("expected 2 LONG_RUN_CONTACT_DIRECTION rows (HYBRID, LGA)",
                     2, longRunRows.size());
        List<String> tensorRows = LINES.stream()
                                         .filter(l -> l.startsWith("HOPPING_TENSOR\t"))
                                         .toList();
        assertEquals("expected 2 HOPPING_TENSOR rows (HYBRID, LGA)", 2, tensorRows.size());

        // Non-tautological cross-check: recompute the tensor from the
        // artifact's own LONG_RUN_CONTACT_DIRECTION weights and compare
        // against the artifact's own HOPPING_TENSOR row.
        for (String longRunRow : longRunRows) {
            String[] parts = longRunRow.split("\t");
            // recordType substrate d1 d2 d3 d4 d5 d6 total
            String substrate = parts[1];
            java.util.Map<Integer, Long> weights = new java.util.LinkedHashMap<>();
            for (int d = 1; d <= 6; d++) {
                weights.put(d, Long.parseLong(parts[1 + d]));
            }
            double[] expectedTensor = PhaseCMeasurement.hoppingTensor(weights);
            String tensorRow = tensorRows.stream()
                                           .filter(l -> l.startsWith("HOPPING_TENSOR\t" + substrate + "\t"))
                                           .findFirst()
                                           .orElseThrow(() -> new AssertionError("missing HOPPING_TENSOR row for "
                                                                                  + substrate));
            String[] tensorParts = tensorRow.split("\t");
            for (int i = 0; i < 6; i++) {
                double actual = Double.parseDouble(tensorParts[2 + i]);
                assertEquals("HOPPING_TENSOR component " + i + " for " + substrate
                             + " must match a live recomputation from the artifact's own weights",
                             expectedTensor[i], actual, 1e-3);
            }
        }

        boolean hasNote = LINES.stream()
                                 .anyMatch(l -> l.startsWith("# contactDirectionCensusNote="));
        assertTrue("expected the contactDirectionCensusNote framing line", hasNote);
    }

    /**
     * Round-3 item 1 (critic S-NEW-1 / reviewer R2-6): Section 5's posture
     * lines must carry a NEUTRAL bare pointer to Section 1B's census/
     * tensor evidence -- selecting/advocating nothing.
     */
    @Test
    public void section5CrossReferencesSection1BCensus() {
        String line = LINES.stream()
                             .filter(l -> l.startsWith("# section1BCrossReference="))
                             .findFirst()
                             .orElse(null);
        assertTrue("expected a section1BCrossReference line", line != null);
        assertTrue("expected the pointer to name LGA_CONTACT_DIRECTION_",
                   line.contains("LGA_CONTACT_DIRECTION_"));
        assertTrue("expected the pointer to name HOPPING_TENSOR",
                   line.contains("HOPPING_TENSOR"));
        assertTrue("expected the pointer to explicitly disclaim advocacy",
                   line.contains("no posture is selected or advocated"));

        // The line must live within Section 5 (after its header, before
        // ESCALATION), not merely exist anywhere in the file.
        int section5Index = LINES.indexOf("# === SECTION 5: ISOTROPY POSTURES -- DATA ONLY, NO POSTURE SELECTED (escalated to user/orchestrator) ===");
        int crossRefIndex = LINES.indexOf(line);
        int escalationIndex = -1;
        for (int i = 0; i < LINES.size(); i++) {
            if (LINES.get(i).startsWith("# ESCALATION=")) {
                escalationIndex = i;
                break;
            }
        }
        assertTrue("expected SECTION 5 header to be found", section5Index >= 0);
        assertTrue("expected ESCALATION line to be found", escalationIndex >= 0);
        assertTrue("expected the cross-reference to sit between the SECTION 5 header and ESCALATION",
                   crossRefIndex > section5Index && crossRefIndex < escalationIndex);
    }

    /**
     * Round-3 item 2a (critic S-NEW-2a): the chi-squared's UNIFORM-null
     * premise (all six directions symmetry-equivalent with identical
     * per-tick exposure) must be stated in-artifact.
     */
    @Test
    public void chiSquarePooledRowStatesUniformNullPremise() {
        String line = LINES.stream()
                             .filter(l -> l.startsWith("# uniformNullPremise="))
                             .findFirst()
                             .orElse(null);
        assertTrue("expected a uniformNullPremise line", line != null);
        assertTrue("expected the premise to name the <110>-type / O_h-equivalence fact",
                   line.contains("<110>-type") && line.contains("O_h-equivalent"));
        assertTrue("expected the premise to name the equal-per-tick-exposure fact",
                   line.contains("exactly once per tick"));
    }

    /**
     * Round-3 item 2b (critic S-NEW-2b): a pooled p-value that underflows
     * to a bare 0.000000000e+00 must never be printed without an
     * annotation -- red-first, the pre-fix artifact has exactly this bare
     * zero with no companion note.
     */
    @Test
    public void chiSquarePooledPValueIsNotABareUnderflowZero() {
        String pooledRow = LINES.stream()
                                  .filter(l -> l.startsWith("LGA_CONTACT_DIRECTION_POOLED\t"))
                                  .findFirst()
                                  .orElseThrow(() -> new AssertionError("missing LGA_CONTACT_DIRECTION_POOLED row"));
        String[] parts = pooledRow.split("\t");
        double pValue = Double.parseDouble(parts[10]);
        if (pValue == 0.0) {
            String note = LINES.stream()
                                 .filter(l -> l.startsWith("# pValueUnderflowNote="))
                                 .findFirst()
                                 .orElse(null);
            assertTrue("expected a pValueUnderflowNote line when the pooled pValue underflows to "
                       + "an exact 0.0 -- a bare 0.000000000e+00 is ambiguous between \"computed and "
                       + "zero\" and \"not computed\"", note != null);
            assertTrue("expected the note to state this is a double underflow, not a computed exact zero",
                       note.contains("DOUBLE UNDERFLOW") && note.contains("not a computed exact zero"));
            assertTrue("expected the note to state the true value is below a stated bound",
                       note.contains("p < 1.0e-300"));
        }
    }

    /**
     * Round-3 item 2c (critic S-NEW-2d): effective-transfer-only
     * per-direction counts (the population that actually moves quanta)
     * must be emitted alongside the raw recorded-contact census, and must
     * reconcile with Section 4B's pooled effectiveCollisions total (same
     * underlying CollisionStatistics#effectiveCollisions(), broken out by
     * direction here) -- additive only, non-tautological cross-check.
     */
    @Test
    public void effectiveTransferPerDirectionCountsAreEmittedAndReconcile() {
        List<String> perSeedRows = LINES.stream()
                                          .filter(l -> l.startsWith("LGA_CONTACT_DIRECTION_EFFECTIVE_SEED\t"))
                                          .toList();
        assertEquals("expected 8 per-seed LGA_CONTACT_DIRECTION_EFFECTIVE_SEED rows",
                     8, perSeedRows.size());

        String pooledRow = LINES.stream()
                                  .filter(l -> l.startsWith("LGA_CONTACT_DIRECTION_EFFECTIVE_POOLED\t"))
                                  .findFirst()
                                  .orElseThrow(() -> new AssertionError("missing LGA_CONTACT_DIRECTION_EFFECTIVE_POOLED row"));
        String[] pooledParts = pooledRow.split("\t");
        // recordType d1 d2 d3 d4 d5 d6 total
        assertEquals(8, pooledParts.length);
        long pooledEffectiveTotal = Long.parseLong(pooledParts[7]);

        long summedFromSeeds = 0;
        for (String row : perSeedRows) {
            String[] parts = row.split("\t");
            long rowTotal = Long.parseLong(parts[8]);
            long directionSum = 0;
            for (int i = 2; i <= 7; i++) {
                directionSum += Long.parseLong(parts[i]);
            }
            assertEquals("expected the row's own total column to equal the sum of its 6 direction columns",
                         directionSum, rowTotal);
            summedFromSeeds += rowTotal;
        }
        assertEquals("expected the pooled row's total to equal the sum of per-seed totals",
                     summedFromSeeds, pooledEffectiveTotal);

        String equalityPooledRow = LINES.stream()
                                          .filter(l -> l.startsWith("ANISOTROPY_CAMPAIGN_EQUALITY_POOLED\t"))
                                          .findFirst()
                                          .orElseThrow(() -> new AssertionError("missing ANISOTROPY_CAMPAIGN_EQUALITY_POOLED row"));
        long pooledEffectiveFromEquality = Long.parseLong(equalityPooledRow.split("\t")[2]);
        assertEquals("expected the effective-transfer-per-direction pooled total to reconcile EXACTLY "
                     + "with ANISOTROPY_CAMPAIGN_EQUALITY_POOLED's effectiveCollisions -- same underlying "
                     + "CollisionStatistics#effectiveCollisions(), broken out by direction here",
                     pooledEffectiveFromEquality, pooledEffectiveTotal);
    }

    // ------------------------------------------------------------------
    // Fix-round item 5 (QUANTA_HISTOGRAM_BIN legend + actual histogram).
    // ------------------------------------------------------------------

    @Test
    public void quantaHistogramCarriesAnActualPerValueHistogram() {
        boolean hasSummaryLegend = LINES.stream()
                                          .anyMatch(l -> l.startsWith("# columns(QUANTA_HISTOGRAM)="));
        boolean hasBinLegend = LINES.stream()
                                      .anyMatch(l -> l.startsWith("# columns(QUANTA_HISTOGRAM_BIN)="));
        assertTrue("expected a columns() legend for QUANTA_HISTOGRAM", hasSummaryLegend);
        assertTrue("expected a columns() legend for QUANTA_HISTOGRAM_BIN", hasBinLegend);

        for (String label : List.of("before", "after")) {
            String summaryRow = LINES.stream()
                                       .filter(l -> l.startsWith("QUANTA_HISTOGRAM\t" + label + "\t"))
                                       .findFirst()
                                       .orElseThrow(() -> new AssertionError("missing QUANTA_HISTOGRAM row for "
                                                                              + label));
            long n = Long.parseLong(summaryRow.split("\t")[2]);
            List<String> binRows = LINES.stream()
                                          .filter(l -> l.startsWith("QUANTA_HISTOGRAM_BIN\t" + label + "\t"))
                                          .toList();
            assertTrue("expected at least one QUANTA_HISTOGRAM_BIN row for " + label,
                       !binRows.isEmpty());
            long binSum = 0;
            for (String binRow : binRows) {
                binSum += Long.parseLong(binRow.split("\t")[3]);
            }
            assertEquals("expected QUANTA_HISTOGRAM_BIN counts to sum EXACTLY to n for " + label,
                         n, binSum);
        }
    }

    /**
     * Round-3 item 3 (critic S-NEW-3 / reviewer R2-5): the bins-boundedness
     * claim rested on a FALSE whole-system maximum-principle argument
     * (removed from {@code QuantaHistogramSummary}'s Javadoc, see the
     * corrected version citing {@code LatticeGasAutomaton}'s
     * {@code checkExactnessCeiling} guard instead). What IS true and now
     * asserted here: the bins reconcile EXACTLY with their own
     * QUANTA_HISTOGRAM summary row (n, mean, variance), and total quanta
     * is conserved across the run (recomputed bins-weighted sum, before
     * == after) -- deliberately NOT a per-tick +/-1 bound, which is false.
     */
    @Test
    public void quantaHistogramBinsReconcileWithSummaryRowAndConserveAcrossRun() {
        Map<String, Long> weightedSumByLabel = new LinkedHashMap<>();
        for (String label : List.of("before", "after")) {
            String summaryRow = LINES.stream()
                                       .filter(l -> l.startsWith("QUANTA_HISTOGRAM\t" + label + "\t"))
                                       .findFirst()
                                       .orElseThrow(() -> new AssertionError("missing QUANTA_HISTOGRAM row for "
                                                                              + label));
            String[] summaryParts = summaryRow.split("\t");
            // recordType label n mean variance min max
            long n = Long.parseLong(summaryParts[2]);
            double summaryMean = Double.parseDouble(summaryParts[3]);
            double summaryVariance = Double.parseDouble(summaryParts[4]);

            List<String> binRows = LINES.stream()
                                          .filter(l -> l.startsWith("QUANTA_HISTOGRAM_BIN\t" + label + "\t"))
                                          .toList();
            long weightedSum = 0;
            long weightedSqSum = 0;
            long recomputedN = 0;
            for (String binRow : binRows) {
                String[] binParts = binRow.split("\t");
                long value = Long.parseLong(binParts[2]);
                long count = Long.parseLong(binParts[3]);
                weightedSum += value * count;
                weightedSqSum += value * value * count;
                recomputedN += count;
            }
            assertEquals("expected bins-recomputed n to equal the summary row's n for " + label,
                         n, recomputedN);
            double recomputedMean = (double) weightedSum / recomputedN;
            assertEquals("expected bins to reconcile with the QUANTA_HISTOGRAM summary row's mean for "
                         + label, summaryMean, recomputedMean, 1e-6);
            double recomputedVariance = (double) weightedSqSum / recomputedN
                                         - recomputedMean * recomputedMean;
            assertEquals("expected bins to reconcile with the QUANTA_HISTOGRAM summary row's variance "
                         + "for " + label, summaryVariance, recomputedVariance, 1e-6);
            weightedSumByLabel.put(label, weightedSum);
        }
        assertEquals("expected total quanta to be conserved across the run (bins-recomputed sum "
                     + "before == after)",
                     weightedSumByLabel.get("before"), weightedSumByLabel.get("after"));
    }

    // ------------------------------------------------------------------
    // Fix-round Critical 1 correction: powerRecommendationForGate must
    // subordinate the seeds/ticks arithmetic to the absorbing-IC seeding
    // caveat, framed as a USER decision.
    // ------------------------------------------------------------------

    @Test
    public void powerRecommendationStatesSeedingCaveatAsAUserDecision() {
        String line = LINES.stream()
                             .filter(l -> l.startsWith("# powerRecommendationForGate="))
                             .findFirst()
                             .orElse(null);
        assertTrue("expected a powerRecommendationForGate line", line != null);
        assertTrue("expected the line to state scaling alone cannot fix the absorbing-IC contamination",
                   line.contains("ABSORBING") && line.contains("CANNOT"));
        assertTrue("expected the line to frame the seeding remediation as a USER decision",
                   line.contains("USER DESIGN DECISION"));
        assertTrue("expected the seeds/ticks power arithmetic to still be present",
                   line.contains("seeds 8->24") && line.contains("ticks 128->400-500"));
    }

    /**
     * Round-3 item 7 (reviewer R2-4 / critic MINOR-1): the seeding
     * caveat's label previously inverted priority ("SUBORDINATE to the
     * arithmetic above") when the caveat is actually what GOVERNS the
     * arithmetic's validity.
     */
    @Test
    public void seedingCaveatGovernsTheArithmeticNotSubordinateToIt() {
        String line = LINES.stream()
                             .filter(l -> l.startsWith("# powerRecommendationForGate="))
                             .findFirst()
                             .orElse(null);
        assertTrue("expected a powerRecommendationForGate line", line != null);
        assertFalse("expected the inverted-priority SUBORDINATE label to be gone",
                    line.contains("SUBORDINATE to the arithmetic"));
        assertTrue("expected the caveat to state it GOVERNS the arithmetic above",
                   line.contains("GOVERNS the arithmetic above"));
    }

    @Test
    public void postureIiiEvidenceCarriesTheContaminationCaveat() {
        String line = LINES.stream()
                             .filter(l -> l.startsWith("# posture(iii)_orientationalStateRestoresIsotropy_evidence="))
                             .findFirst()
                             .orElse(null);
        assertTrue("expected a posture(iii) evidence line", line != null);
        assertTrue("expected the posture(iii) evidence to carry its own contamination caveat",
                   line.contains("CONTAMINATION CAVEAT"));
    }

    @Test
    public void collisionFieldTotalCollisionsNoteMatchesWhatTest8ActuallyAsserts() {
        String row = LINES.stream()
                            .filter(l -> l.startsWith("COLLISION_FIELD\ttotalCollisions\t"))
                            .findFirst()
                            .orElse(null);
        assertTrue("expected a COLLISION_FIELD totalCollisions row", row != null);
        String[] parts = row.split("\t");
        assertEquals("expected verifiedComparable=true now that a real assertion backs it",
                     "true", parts[4]);
        assertTrue("expected the note to name the actual asserting test method",
                   parts[5].contains("aggregateStatisticsAgreeBeyondDivergence"));
    }

    // ------------------------------------------------------------------
    // Pure-helper unit tests (fix-round item 3 / Important #2): direct
    // synthetic-input coverage of PhaseCMeasurement's narrative-generating
    // pure helpers, with NO campaign/artifact dependency -- exactly the
    // coverage gap the reviewer named as what let C1 (inverted lag
    // direction) and C3 (hardcoded conclusions) ship undetected. Each
    // covers every branch the corresponding narrative literal takes.
    // ------------------------------------------------------------------

    @Test
    public void lagNoteStatesNoOpTrendAsSubjectWhenOneWindowLater() {
        String note = PhaseCMeasurement.lagNote(2, 1);
        assertTrue(note, note.contains("this trend's own peak (window 2) occurs exactly ONE WINDOW LATER than"));
        assertTrue(note, note.contains("Section 3's relativeRateGapByWindow peak (window 1)"));
    }

    @Test
    public void lagNoteStatesDifferentWindowWhenNotOneWindowLater() {
        String same = PhaseCMeasurement.lagNote(1, 1);
        assertTrue(same, same.contains("this trend's own peak (window 1) occurs at a different window than"));
        assertTrue(same, same.contains("Section 3's relativeRateGapByWindow peak (window 1)"));

        // The direction-inversion case the pre-fix bug would have gotten
        // wrong: the no-op peak is EARLIER than the gap peak, not later.
        String reversed = PhaseCMeasurement.lagNote(0, 2);
        assertFalse("must NOT claim ONE WINDOW LATER when the no-op peak is actually earlier: " + reversed,
                    reversed.contains("ONE WINDOW LATER"));
        assertTrue(reversed, reversed.contains("at a different window than"));
    }

    @Test
    public void trendShapeStatementBranchesOnFirstVsLast() {
        assertEquals("OVERALL RISING", PhaseCMeasurement.trendShapeStatement(0.1, 0.9));
        assertEquals("OVERALL FALLING", PhaseCMeasurement.trendShapeStatement(0.9, 0.1));
        assertEquals("OVERALL FLAT (first == last)", PhaseCMeasurement.trendShapeStatement(0.5, 0.5));
    }

    @Test
    public void sdbCorroborationClauseBranchesOnFirstVsLast() {
        assertTrue(PhaseCMeasurement.sdbCorroborationClause(0.1, 0.9)
                                    .contains("dynamical corroboration"));
        assertTrue(PhaseCMeasurement.sdbCorroborationClause(0.9, 0.1)
                                    .contains("NOT corroborating"));
        assertTrue(PhaseCMeasurement.sdbCorroborationClause(0.5, 0.5)
                                    .contains("NOT corroborating"));
    }

    @Test
    public void ciExcludesOneBranchesOnBounds() {
        assertTrue(PhaseCMeasurement.ciExcludesOne(1.05, 1.5));
        assertTrue(PhaseCMeasurement.ciExcludesOne(0.5, 0.9));
        assertFalse(PhaseCMeasurement.ciExcludesOne(0.9, 1.1));
    }

    @Test
    public void isSignificantBranchesOnPValue() {
        assertTrue(PhaseCMeasurement.isSignificant(0.01));
        assertTrue(PhaseCMeasurement.isSignificant(0.05));
        assertFalse(PhaseCMeasurement.isSignificant(0.06));
        assertFalse(PhaseCMeasurement.isSignificant(0.225));
    }

    @Test
    public void ciExclusionAndSignificanceClausesMatchBooleans() {
        assertEquals("(excludes 1.0)", PhaseCMeasurement.ciExclusionClause(true));
        assertEquals("(includes 1.0)", PhaseCMeasurement.ciExclusionClause(false));
        assertEquals("(significant)", PhaseCMeasurement.significanceClause(true));
        assertEquals("(NOT significant)", PhaseCMeasurement.significanceClause(false));
    }

    @Test
    public void ciVsPermutationConclusionCoversAllFourCombinations() {
        String both = PhaseCMeasurement.ciVsPermutationConclusion(true, false, true, false);
        assertTrue(both, both.contains("in BOTH campaigns"));

        String neither = PhaseCMeasurement.ciVsPermutationConclusion(false, true, false, true);
        assertTrue(neither, neither.contains("NEITHER campaign"));

        String aOnly = PhaseCMeasurement.ciVsPermutationConclusion(true, false, false, true);
        assertTrue(aOnly, aOnly.contains("DISAGREE"));
        assertTrue(aOnly, aOnly.contains("Phase A shows it, Phase C does not"));

        String cOnly = PhaseCMeasurement.ciVsPermutationConclusion(false, true, true, false);
        assertTrue(cOnly, cOnly.contains("DISAGREE"));
        assertTrue(cOnly, cOnly.contains("Phase C shows it, Phase A does not"));
    }

    @Test
    public void belowThresholdClauseCoversAllFourCombinations() {
        assertEquals("BOTH fall below", PhaseCMeasurement.belowThresholdClause(true, true));
        assertEquals("Phase A falls below (but Phase C does NOT fall below)",
                     PhaseCMeasurement.belowThresholdClause(true, false));
        assertEquals("Phase C falls below (but Phase A does NOT fall below)",
                     PhaseCMeasurement.belowThresholdClause(false, true));
        assertEquals("NEITHER falls below", PhaseCMeasurement.belowThresholdClause(false, false));
    }

    @Test
    public void ridgeVerdictTextCoversAllThreeBranches() {
        assertTrue(PhaseCMeasurement.ridgeVerdictText(true, false, false).startsWith("INSTRUMENT ANOMALY"));
        assertTrue(PhaseCMeasurement.ridgeVerdictText(true, true, true).startsWith("INSTRUMENT ANOMALY"));
        assertTrue(PhaseCMeasurement.ridgeVerdictText(false, true, false).startsWith("RIDGE PRESENT"));
        assertTrue(PhaseCMeasurement.ridgeVerdictText(false, true, true).startsWith("RIDGE PRESENT"));
        assertTrue(PhaseCMeasurement.ridgeVerdictText(false, false, false).startsWith("RIDGE ABSENT"));
        assertTrue(PhaseCMeasurement.ridgeVerdictText(false, false, true).startsWith("RIDGE ABSENT"));
    }

    /**
     * Fix-round S3: the RIDGE ABSENT branch must NOT assert a dynamics
     * claim ("purely diffusive dynamics") the measurement alone cannot
     * support when the campaign is low-signal -- and must NOT silently
     * omit the low-signal alternative when it applies.
     */
    @Test
    public void ridgeVerdictTextDistinguishesMeasuredAbsenceFromDynamicsClaim() {
        String lowSignal = PhaseCMeasurement.ridgeVerdictText(false, false, true);
        assertTrue(lowSignal, lowSignal.contains("does NOT distinguish"));
        assertTrue(lowSignal, lowSignal.contains("insufficient signal"));
        assertTrue(lowSignal, lowSignal.contains("no dynamics claim is made"));

        String wellPowered = PhaseCMeasurement.ridgeVerdictText(false, false, false);
        assertFalse(wellPowered, wellPowered.contains("does NOT distinguish"));
        assertTrue(wellPowered, wellPowered.contains("purely-diffusive-dynamics signature"));
    }

    // ------------------------------------------------------------------
    // Fix-round item 1/2/5 pure-helper unit tests: additive-instrumentation
    // helpers, each covering every branch/edge with synthetic inputs, no
    // campaign/artifact dependency.
    // ------------------------------------------------------------------

    @Test
    public void evenParityCellCountIsXTimesYTimesZOverTwo() {
        assertEquals(256L, PhaseCMeasurement.evenParityCellCount(new Point3i(8, 8, 8)));
        assertEquals(32L, PhaseCMeasurement.evenParityCellCount(new Point3i(4, 4, 4)));
        assertEquals(16L, PhaseCMeasurement.evenParityCellCount(new Point3i(4, 4, 2)));
    }

    @Test
    public void activeMemberCensusCountsExcitedMembersAndTouchedCellsOnly() {
        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = new Point3i(1, 1, 0);
        long[] quanta = new long[60];
        // cellA (slots 0..29): all zero -- untouched.
        // cellB (slots 30..59): exactly 2 nonzero members.
        quanta[30] = 5L;
        quanta[45] = -3L;
        QuantaField field = fakeField(List.of(cellA, cellB), quanta);

        PhaseCMeasurement.ActiveCensus census = PhaseCMeasurement.activeMemberCensus(field);
        assertEquals(2L, census.excitedMembers());
        assertEquals(1, census.cellsTouched());
    }

    @Test
    public void activeMemberCensusReportsZeroWhenEverythingIsZero() {
        Point3i cellA = new Point3i(0, 0, 0);
        QuantaField field = fakeField(List.of(cellA), new long[30]);
        PhaseCMeasurement.ActiveCensus census = PhaseCMeasurement.activeMemberCensus(field);
        assertEquals(0L, census.excitedMembers());
        assertEquals(0, census.cellsTouched());
    }

    @Test
    public void pooledDirectionCountsSumsAcrossSeedsForDirectionsOneToSix() {
        List<PhaseCMeasurement.ContactDirectionCensus> perSeed = new ArrayList<>();
        perSeed.add(new PhaseCMeasurement.ContactDirectionCensus("s1",
                                                                   directionMap(1L, 2L, 3L, 4L, 5L, 6L)));
        perSeed.add(new PhaseCMeasurement.ContactDirectionCensus("s2",
                                                                   directionMap(10L, 20L, 30L, 40L, 50L, 60L)));
        Map<Integer, Long> pooled = PhaseCMeasurement.pooledDirectionCounts(perSeed);
        assertEquals(6, pooled.size());
        assertEquals(11L, (long) pooled.get(1));
        assertEquals(22L, (long) pooled.get(2));
        assertEquals(66L, (long) pooled.get(6));
    }

    @Test
    public void chiSquareStatisticIsZeroForPerfectlyUniformCounts() {
        assertEquals(0.0, PhaseCMeasurement.chiSquareStatistic(new long[] { 100, 100, 100, 100, 100, 100 }),
                     1e-9);
    }

    @Test
    public void chiSquareStatisticIsPositiveAndLargerForMoreSkewedCounts() {
        double mild = PhaseCMeasurement.chiSquareStatistic(new long[] { 90, 100, 110, 100, 100, 100 });
        double severe = PhaseCMeasurement.chiSquareStatistic(new long[] { 10, 100, 190, 100, 100, 100 });
        assertTrue(mild > 0.0);
        assertTrue(severe > mild);
    }

    /**
     * Cross-check: the critic's own pooled anisotropy-campaign counts,
     * RAW (non-independence-adjusted) chi2=5496.07 on 5 df. Round-3 item
     * 8a (reviewer MINOR-3) correction: the raw statistic on these exact
     * counts is ~5496, NOT ~43 -- chi2~43 was the critic's own SEPARATE
     * independence-DEFLATED estimate (an effective-sample-size adjustment
     * this code deliberately does not compute, per {@code
     * chiSquareCaveat}'s framing), and does not belong to this
     * computation.
     */
    @Test
    public void chiSquareStatisticReproducesCriticMeasuredPooledCounts() {
        double stat = PhaseCMeasurement.chiSquareStatistic(new long[] { 1670, 6167, 7535, 5605, 2432, 4347 });
        assertTrue("expected the raw (non-independence-adjusted) pooled chi-squared to be large, was "
                   + stat, stat > 1000.0);
        assertEquals("expected the raw pooled chi-squared to be ~5496.07 (NOT the critic's separate "
                     + "~43 independence-deflated estimate)", 5496.07, stat, 0.5);
    }

    @Test
    public void chiSquarePValueIsOneWhenStatisticIsZero() {
        assertEquals(1.0, PhaseCMeasurement.chiSquarePValue(0.0, 5), 1e-12);
    }

    /** Textbook chi-squared critical values (df, x, p~0.05 or p~0.01). */
    @Test
    public void chiSquarePValueMatchesTextbookCriticalValues() {
        assertEquals(0.05, PhaseCMeasurement.chiSquarePValue(3.841459, 1), 1e-3);
        assertEquals(0.05, PhaseCMeasurement.chiSquarePValue(11.0705, 5), 1e-3);
        assertEquals(0.01, PhaseCMeasurement.chiSquarePValue(15.0863, 5), 1e-3);
    }

    @Test
    public void chiSquarePValueIsMonotonicallyDecreasingInStatistic() {
        double p1 = PhaseCMeasurement.chiSquarePValue(1.0, 5);
        double p2 = PhaseCMeasurement.chiSquarePValue(10.0, 5);
        double p3 = PhaseCMeasurement.chiSquarePValue(50.0, 5);
        assertTrue(p1 > p2);
        assertTrue(p2 > p3);
    }

    @Test
    public void chiSquarePValueRejectsInvalidInputs() {
        try {
            PhaseCMeasurement.chiSquarePValue(-1.0, 5);
            fail("expected IllegalArgumentException for negative chiSquare");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            PhaseCMeasurement.chiSquarePValue(1.0, 0);
            fail("expected IllegalArgumentException for non-positive df");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    /**
     * Regression cross-check against the critic's own hand-derived
     * long-run LGA hopping-tensor numbers (independently verified by
     * hand: D = sum_d w_d*outer(e_d,e_d) over the six positive FCC
     * directions).
     */
    @Test
    public void hoppingTensorReproducesCriticMeasuredLongRunLgaWeights() {
        Map<Integer, Long> weights = directionMap(2439L, 5183L, 15537L, 12220L, 4579L, 7752L);
        double[] tensor = PhaseCMeasurement.hoppingTensor(weights);
        assertEquals(18974.0, tensor[0], 1e-6); // Dxx
        assertEquals(14867.0, tensor[1], 1e-6); // Dyy
        assertEquals(13869.0, tensor[2], 1e-6); // Dzz
        assertEquals(6549.0, tensor[3], 1e-6);  // Dxz
        assertEquals(-302.0, tensor[4], 1e-6);  // Dyz
        assertEquals(2234.0, tensor[5], 1e-6);  // Dxy
    }

    @Test
    public void hoppingTensorIsIsotropicForEqualDirectionWeights() {
        Map<Integer, Long> weights = directionMap(100L, 100L, 100L, 100L, 100L, 100L);
        double[] tensor = PhaseCMeasurement.hoppingTensor(weights);
        assertEquals(tensor[0], tensor[1], 1e-9);
        assertEquals(tensor[1], tensor[2], 1e-9);
        assertEquals(0.0, tensor[3], 1e-9);
        assertEquals(0.0, tensor[4], 1e-9);
        assertEquals(0.0, tensor[5], 1e-9);
    }

    @Test
    public void histogramOfComputesExactBinsReconcilingToN() {
        long[] quanta = { -2, -2, 0, 0, 0, 1, 2 };
        PhaseCMeasurement.QuantaHistogramSummary h = PhaseCMeasurement.histogramOf("test", quanta);
        assertEquals(7L, h.n());
        assertEquals(-2L, h.min());
        assertEquals(2L, h.max());
        assertEquals(5, h.bins().size()); // -2,-1,0,1,2
        assertEquals(2L, (long) h.bins().get(-2L));
        assertEquals(0L, (long) h.bins().get(-1L));
        assertEquals(3L, (long) h.bins().get(0L));
        assertEquals(1L, (long) h.bins().get(1L));
        assertEquals(1L, (long) h.bins().get(2L));
        long binSum = 0;
        for (long v : h.bins().values()) {
            binSum += v;
        }
        assertEquals(h.n(), binSum);
    }

    private static Map<Integer, Long> directionMap(long d1, long d2, long d3, long d4,
                                                     long d5, long d6) {
        Map<Integer, Long> m = new LinkedHashMap<>();
        m.put(1, d1);
        m.put(2, d2);
        m.put(3, d3);
        m.put(4, d4);
        m.put(5, d5);
        m.put(6, d6);
        return m;
    }

    /** Minimal in-test {@link QuantaField} double -- identity-keyed cell lookup, no equals() reliance. */
    private static QuantaField fakeField(List<Point3i> cells, long[] quanta) {
        Map<Point3i, Integer> baseByCell = new java.util.IdentityHashMap<>();
        for (int i = 0; i < cells.size(); i++) {
            baseByCell.put(cells.get(i), i * 30);
        }
        return new QuantaField() {
            @Override
            public Point3i extent() {
                return new Point3i(2, 2, 2);
            }

            @Override
            public int slotCount() {
                return quanta.length;
            }

            @Override
            public long quantaAt(int slot) {
                return quanta[slot];
            }

            @Override
            public boolean isExactAt(int slot) {
                return true;
            }

            @Override
            public float phaseAt(int slot) {
                return 0f;
            }

            @Override
            public int phaseResolution() {
                return 3600;
            }

            @Override
            public void forEachCell(Consumer<? super Point3i> action) {
                for (Point3i cell : cells) {
                    action.accept(cell);
                }
            }

            @Override
            public int indexOfCell(Point3i cell) {
                return baseByCell.get(cell);
            }
        };
    }

    // ------------------------------------------------------------------
    // Helpers.
    // ------------------------------------------------------------------

    private static double[] spectralSummaryRow(String substrate) {
        String row = LINES.stream()
                            .filter(l -> l.startsWith("SPECTRAL_SUMMARY\t" + substrate
                                                       + "\t"))
                            .findFirst()
                            .orElseThrow(() -> new AssertionError("missing SPECTRAL_SUMMARY row for "
                                                                   + substrate));
        String[] parts = row.split("\t");
        // recordType substrate extent fftLength stride nMembersChecked meanPeakFraction
        // minPeakFraction maxPeakFraction meanEntropy meanAbsoluteLinewidthRadPerTick
        // meanFractionalLinewidth nMembersFractionalDefined
        // index:     0          1         2      3         4      5                6
        //            7               8              9            10                        11                   12
        return new double[] { Double.parseDouble(parts[6]),
                               Double.parseDouble(parts[7]),
                               Double.parseDouble(parts[8]),
                               Double.parseDouble(parts[9]),
                               Double.parseDouble(parts[10]),
                               Double.parseDouble(parts[11]) };
    }

    private static void assertHeaderPresent(String prefix) {
        assertTrue("expected header line starting with \"" + prefix + "\"",
                   LINES.stream().anyMatch(l -> l.equals(prefix) || l.startsWith(prefix)));
    }

    private static void assertHeaderStartsWithNonEmptyValue(String prefix) {
        String line = LINES.stream()
                             .filter(l -> l.startsWith(prefix))
                             .findFirst()
                             .orElse(null);
        assertTrue("expected header line starting with \"" + prefix + "\"", line != null);
        assertTrue("expected non-empty value after \"" + prefix + "\", was: " + line,
                   line.length() > prefix.length());
    }

    /** Extracts the {@code [a,b,c,...]} series following {@code key} in {@code line}. */
    private static List<Double> parseSeries(String line, String key) {
        int keyStart = line.indexOf(key);
        int open = line.indexOf('[', keyStart);
        int close = line.indexOf(']', open);
        String inner = line.substring(open + 1, close);
        List<Double> values = new ArrayList<>();
        for (String part : inner.split(",")) {
            values.add(Double.parseDouble(part.trim()));
        }
        return values;
    }

    private static int argmaxOf(List<Double> values) {
        int best = 0;
        for (int i = 1; i < values.size(); i++) {
            if (values.get(i) > values.get(best)) {
                best = i;
            }
        }
        return best;
    }

    private static long headerLong(String key) {
        String prefix = "# " + key + "=";
        String line = LINES.stream()
                             .filter(l -> l.startsWith(prefix))
                             .findFirst()
                             .orElseThrow(() -> new AssertionError("missing header key " + key));
        String value = line.substring(prefix.length());
        int space = value.indexOf(' ');
        if (space >= 0) {
            value = value.substring(0, space);
        }
        return Long.parseLong(value.trim());
    }
}
