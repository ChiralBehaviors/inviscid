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
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Locale;

import org.junit.BeforeClass;
import org.junit.Test;

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
        assertTrue("expected a nonzero number of ticks to have been strictly conservation-audited across the campaign",
                   ticksAudited > 0);
        assertEquals("expected EXACTLY zero conservation violations across the entire measurement campaign",
                     0L, violations);
        boolean strictMode = LINES.stream()
                                   .anyMatch(l -> l.startsWith("# conservationMode=STRICT"));
        assertTrue("expected the report to state strict-mode auditing (every driver, every sub-measurement)",
                   strictMode);
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
