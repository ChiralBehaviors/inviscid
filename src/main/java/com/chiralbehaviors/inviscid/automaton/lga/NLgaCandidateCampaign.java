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

package com.chiralbehaviors.inviscid.automaton.lga;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.automaton.lga.ContactComboCache.Combo;

/**
 * The N_lga candidate measurement campaign (bead inviscid-0nx.16's
 * "THE MEASUREMENT CAMPAIGN"): quantifies, for each candidate {@code nLga
 * in {8, 12, 16, 24}}, how well that bin grid resolves the actual
 * contact-angle structure {@link ContactPredicate} exposes - the data the
 * USER's N_lga decision (reserved to the epic, never made by this class)
 * is informed by. Read-only / analysis-only: never writes the committed
 * atlas.
 *
 * <h2>Angle-quantization alignment fix (bead inviscid-gyt follow-up,
 * corrects inviscid-0nx.16's original campaign numbers)</h2>
 * This class used to reconstruct fine-grid angles with its own {@code
 * step * 2*Math.PI/resolution} (double, LEFT-EDGE) formula and bin them
 * with a matching double-precision {@code floor(angle/binWidth)} - the
 * exact bug class discovered and fixed in {@link ContactComboCache} and
 * {@link ContactAtlasGenerator} during bead inviscid-gyt (any-overlap
 * transcription semantics): {@link MemberGeometry#stepOf} quantizes with
 * {@code Constants.TWO_PI} (float precision, compounded float-modulo then
 * division), so a double-precision left-edge reconstruction round-trips
 * incorrectly for roughly 38% of the 360 steps - it silently evaluates the
 * WRONG step's geometry for over a third of the fine sweep. This class now
 * REUSES the corrected, canonical helpers instead of maintaining a third
 * copy: {@link ContactComboCache#angleOf} (step-CENTER, {@code
 * Constants.TWO_PI} float) for fine-grid angle reconstruction, and {@link
 * ContactAtlasGenerator#binOfStep} (pure integer {@code step*nLga/
 * geometryResolution}, no floating point) for step-to-bin mapping - both
 * package-private in this same package, no visibility widening required.
 * See T2 {@code inviscid/analysis-nlga-candidates.md}'s CORRECTION-2
 * section for the re-run campaign numbers this fix produced.
 *
 * <h2>Two independently-measured inputs</h2>
 * <ol>
 * <li>Combo discovery - delegates to {@link ContactComboCache#combosFor},
 * the same EXHAUSTIVE (native {@link #GEOMETRY_RESOLUTION}=360-step, the
 * provable ceiling - {@code MemberGeometry.stepOf} floors every continuous
 * angle onto one of exactly 360 discrete LUT steps before any geometry is
 * computed, so no angle grid finer than 360 steps can ever reveal a
 * combination invisible at 360 steps) cached discovery every other
 * atlas-generation caller uses, rather than a duplicate discovery sweep of
 * this class's own. See {@link ContactComboCache}'s class Javadoc for the
 * full "why precomputed, not live" / "why a coarse pre-scan was rejected"
 * argument, and its own history of a prior 24-step under-count (bead
 * inviscid-0nx.16.1) this class's earlier version independently
 * reproduced: a non-exhaustive discovery step count is a silent
 * sampling-bias risk for exactly the "how many combos are dangerously
 * narrow" question this campaign exists to answer.</li>
 * <li>{@link #measureWidth} - for each discovered combo, a {@link
 * #FINE_STEPS}x{@link #FINE_STEPS} grid (one sweep serves both discovery
 * cross-checking and width measurement), from which {@link #circularWidth}
 * extracts the smallest circular arc enclosing every contacting {@code
 * angleA} (resp. {@code angleB}) value - the "angular width of the contact
 * region" the bead's campaign asks for.</li>
 * </ol>
 * {@link #analyze(int)} then compares each candidate's bin width ({@code
 * 2*pi/nLga}) against those measured widths (fraction-of-bin-width and
 * "narrower than one bin" risk), and against a bin-center transcription
 * error rate computed directly from the same fine grids. {@link
 * #tableCost(int)} separately drives a real {@link
 * ContactAtlasGenerator#generate} to report row count, wall time, and
 * dynamic-coverage fraction.
 *
 * @author halhildebrand
 */
public final class NLgaCandidateCampaign {

    /**
     * Bead inviscid-eho. The {@code geometryResolution % nLga == 0}
     * invariant these candidates used to be checked against is NO LONGER
     * ENFORCED ANYWHERE - it died with the phase-quantizer class bead
     * inviscid-0nx.27 retired. {@link ContactTable}'s class Javadoc
     * carries the full account; the short version is that {@code
     * nLga=16} below already violates it (measured: {@code 360 % 16 ==
     * 8}), and nothing checks a future candidate. Its call site in this
     * class is {@code binOfFineIndex} -> {@code
     * ContactAtlasGenerator.binOfStep(fineIndex, nLga, FINE_STEPS)} ->
     * {@code analyze}'s {@code transcriptionErrorRate}.
     *
     * <p>The SURVIVING invariant, {@code phaseResolution %
     * geometryResolution == 0} in {@link LatticeGasAutomaton}'s
     * constructor, excludes MORE of these candidates than the retired one
     * does. {@code ContactAtlasGenerator.generate} stamps {@code
     * subBinSteps = 150} into every header unconditionally, so {@code
     * phaseResolution = nLga * 150}, and {@code phaseResolution % 360} is
     * {@code 120} at {@code nLga=8}, {@code 0} at {@code 12}, {@code 240}
     * at {@code 16}, {@code 0} at {@code 24}: only {@code 12} and {@code
     * 24} are LGA-constructible; {@code 8} and {@code 16} throw from that
     * constructor.
     */
    public static final int[]   CANDIDATES          = { 8, 12, 16, 24 };
    public static final double  RADIUS              = ContactAtlasGenerator.RADIUS;
    public static final int     GEOMETRY_RESOLUTION = ContactAtlasGenerator.GEOMETRY_RESOLUTION;
    /** Matches {@code ContactPredicateTest}'s own documented "fine sweep" resolution. */
    public static final int     FINE_STEPS          = 360;
    public static final Point3i EXTENT              = ContactAtlasGenerator.DEFAULT_EXTENT;
    public static final long    SEED                = ContactAtlasGenerator.DEFAULT_SEED;
    public static final int     TICKS               = ContactAtlasGenerator.DEFAULT_TICKS;

    private static final double TWO_PI = 2 * Math.PI;

    private NLgaCandidateCampaign() {
    }

    /**
     * A discovered combo's fine-grid contact region: {@code widthA}/{@code
     * widthB} are the smallest circular arcs (radians) enclosing every
     * contacting {@code angleA}/{@code angleB} respectively; {@code grid}
     * is the raw {@code FINE_STEPS x FINE_STEPS} contact verdict, retained
     * for the per-candidate transcription-error computation.
     */
    public record ComboWidth(Combo combo, double widthA, double widthB,
                              int trueCellCount, boolean[][] grid) {
    }

    public record CandidateStats(int nLga, double binWidthRadians,
                                  double minWidthRadians,
                                  double medianWidthRadians,
                                  double meanWidthRadians,
                                  double maxWidthRadians,
                                  double transcriptionErrorRate,
                                  double fractionNarrowerThanBin,
                                  int atlasRows, long generationWallMs,
                                  double dynamicCoverageFraction) {
    }

    private static ContactPredicate newPredicate() {
        return new ContactPredicate(new MemberGeometry(GEOMETRY_RESOLUTION,
                                                        RADIUS));
    }

    public static ComboWidth measureWidth(ContactPredicate predicate,
                                           Combo combo) {
        boolean[][] grid = new boolean[FINE_STEPS][FINE_STEPS];
        boolean[]   aHasContact = new boolean[FINE_STEPS];
        boolean[]   bHasContact = new boolean[FINE_STEPS];
        int         trueCount = 0;

        for (int a = 0; a < FINE_STEPS; a++) {
            float angleA = ContactComboCache.angleOf(a, FINE_STEPS);
            for (int b = 0; b < FINE_STEPS; b++) {
                float angleB = ContactComboCache.angleOf(b, FINE_STEPS);
                boolean contact = predicate.contacts(combo.cubeA(),
                                                      combo.memberA(), angleA,
                                                      combo.cubeB(),
                                                      combo.memberB(), angleB,
                                                      combo.direction());
                grid[a][b] = contact;
                if (contact) {
                    trueCount++;
                    aHasContact[a] = true;
                    bHasContact[b] = true;
                }
            }
        }

        double widthA = circularWidth(aHasContact, FINE_STEPS);
        double widthB = circularWidth(bHasContact, FINE_STEPS);
        return new ComboWidth(combo, widthA, widthB, trueCount, grid);
    }

    /**
     * The smallest circular arc (radians) enclosing every {@code true}
     * index of {@code marks} (a {@code resolution}-step discretization of
     * {@code [0, 2*pi)}) - i.e. {@code 2*pi} minus the largest circular
     * gap between consecutive marked indices. Degenerate cases: no marks
     * -> {@code 0.0}; every index marked -> {@code 2*pi}; exactly one mark
     * -> {@code 0.0} (a single fine-grid cell is at or below this sweep's
     * own resolution floor - not a genuine "full circle", which the naive
     * single-point gap formula would otherwise produce).
     */
    static double circularWidth(boolean[] marks, int resolution) {
        List<Integer> trueIndices = new ArrayList<>();
        for (int i = 0; i < resolution; i++) {
            if (marks[i]) {
                trueIndices.add(i);
            }
        }
        if (trueIndices.isEmpty() || trueIndices.size() == 1) {
            return 0.0;
        }
        if (trueIndices.size() == resolution) {
            return TWO_PI;
        }
        double binWidth = TWO_PI / resolution;
        int maxGap = 0;
        for (int i = 0; i < trueIndices.size(); i++) {
            int current = trueIndices.get(i);
            int next = trueIndices.get((i + 1) % trueIndices.size());
            int gap = Math.floorMod(next - current, resolution);
            maxGap = Math.max(maxGap, gap);
        }
        int coveredSteps = resolution - maxGap;
        return coveredSteps * binWidth;
    }

    private static double[] percentileStats(List<Double> values) {
        List<Double> sorted = new ArrayList<>(values);
        sorted.sort(Double::compareTo);
        double min = sorted.get(0);
        double max = sorted.get(sorted.size() - 1);
        double median = sorted.size() % 2 == 0
                         ? (sorted.get(sorted.size() / 2 - 1)
                            + sorted.get(sorted.size() / 2)) / 2.0
                         : sorted.get(sorted.size() / 2);
        double sum = 0.0;
        for (double v : sorted) {
            sum += v;
        }
        double mean = sum / sorted.size();
        return new double[] { min, median, mean, max };
    }

    /**
     * @return the {@code nLga} bin index fine-grid index {@code
     *         fineIndex} (out of {@link #FINE_STEPS}) falls into - PURE
     *         INTEGER arithmetic via {@link
     *         ContactAtlasGenerator#binOfStep}, the same step-to-bin
     *         mapping {@link ContactAtlasGenerator#sweepOverlapAndCenter}
     *         uses for its own fine sweep. Deliberately never routed
     *         through a reconstructed floating-point angle - see this
     *         class's Javadoc "Angle-quantization alignment fix" section
     *         for why that reconstruction is unreliable.
     */
    private static int binOfFineIndex(int fineIndex, int nLga) {
        return ContactAtlasGenerator.binOfStep(fineIndex, nLga, FINE_STEPS);
    }

    /**
     * @return the per-candidate width/quantization-fidelity statistics for
     *         {@code nLga}, computed from {@code widths} (already
     *         measured, once, and shared across every candidate).
     */
    public static CandidateStats analyze(int nLga, List<ComboWidth> widths) {
        double binWidth = TWO_PI / nLga;

        List<Double> allWidths = new ArrayList<>(widths.size() * 2);
        int narrowerThanBin = 0;
        long mismatches = 0;
        long totalFineCells = 0;

        for (ComboWidth cw : widths) {
            allWidths.add(cw.widthA());
            allWidths.add(cw.widthB());
            if (Math.min(cw.widthA(), cw.widthB()) < binWidth) {
                narrowerThanBin++;
            }

            boolean[][] binVerdict = new boolean[nLga][nLga];
            ContactPredicate predicate = newPredicate();
            for (int binA = 0; binA < nLga; binA++) {
                float angleA = (float) ContactAtlasGenerator.binCenter(binA,
                                                                        nLga);
                for (int binB = 0; binB < nLga; binB++) {
                    float angleB = (float) ContactAtlasGenerator.binCenter(binB,
                                                                            nLga);
                    binVerdict[binA][binB] = predicate.contacts(cw.combo()
                                                                   .cubeA(),
                                                                 cw.combo()
                                                                   .memberA(),
                                                                 angleA,
                                                                 cw.combo()
                                                                   .cubeB(),
                                                                 cw.combo()
                                                                   .memberB(),
                                                                 angleB,
                                                                 cw.combo()
                                                                   .direction());
                }
            }
            for (int a = 0; a < FINE_STEPS; a++) {
                int binA = binOfFineIndex(a, nLga);
                for (int b = 0; b < FINE_STEPS; b++) {
                    int binB = binOfFineIndex(b, nLga);
                    totalFineCells++;
                    if (cw.grid()[a][b] != binVerdict[binA][binB]) {
                        mismatches++;
                    }
                }
            }
        }

        double[] stats = percentileStats(allWidths);
        double transcriptionErrorRate = totalFineCells == 0 ? 0.0
                                                              : (double) mismatches
                                                                / totalFineCells;
        double fractionNarrowerThanBin = widths.isEmpty() ? 0.0
                                                            : (double) narrowerThanBin
                                                              / widths.size();

        TableCost cost = tableCost(nLga);

        return new CandidateStats(nLga, binWidth, stats[0], stats[1],
                                   stats[2], stats[3], transcriptionErrorRate,
                                   fractionNarrowerThanBin, cost.rows,
                                   cost.wallMs, cost.dynamicCoverageFraction);
    }

    private record TableCost(int rows, long wallMs,
                              double dynamicCoverageFraction) {
    }

    private static TableCost tableCost(int nLga) {
        long start = System.nanoTime();
        ContactAtlas atlas = ContactAtlasGenerator.generate(nLga, EXTENT,
                                                             SEED, TICKS);
        long wallMs = (System.nanoTime() - start) / 1_000_000;

        long contactingBins = 0;
        long dynamicallyReached = 0;
        for (ContactAtlas.Row row : atlas.rows()) {
            if (row.contact()) {
                contactingBins++;
                if (row.observedCount() > 0) {
                    dynamicallyReached++;
                }
            }
        }
        double coverage = contactingBins == 0 ? 0.0
                                                : (double) dynamicallyReached
                                                  / contactingBins;
        return new TableCost(atlas.rows().size(), wallMs, coverage);
    }

    /**
     * A CRUDE, single-data-point extrapolation of {@code ticksObserved}
     * needed to reach a target dynamic-coverage fraction (bead
     * inviscid-0nx.16.1's "quantify the ticksObserved advice" ask):
     * assumes coverage saturates as a Poisson/coupon-collector process,
     * {@code coverage(n) ~= 1 - exp(-k*n)}, and backs {@code k} out of
     * the single {@code (ticksObserved=}{@link #TICKS}{@code ,
     * dynamicCoverageFraction)} point each {@link CandidateStats} already
     * carries - NOT validated against a second empirical run at a
     * different tick count, so treat the result as an order-of-magnitude
     * guide, not a precise budget.
     */
    private static String formatTicksAdvice(List<CandidateStats> results) {
        StringBuilder sb = new StringBuilder();
        sb.append("Crude ticksObserved extrapolation (coupon-collector/Poisson saturation approximation - single-point fit, NOT empirically re-validated at a second tick count; order-of-magnitude guide only):\n");
        sb.append("nLga\tmeasuredCoverageAt").append(TICKS)
          .append("\tticksFor50pctCoverage\tticksFor80pctCoverage\tticksFor95pctCoverage\n");
        for (CandidateStats s : results) {
            double coverage = Math.min(s.dynamicCoverageFraction(), 0.999999);
            double k = -Math.log(1 - coverage) / TICKS;
            sb.append(s.nLga()).append('\t')
              .append(String.format(Locale.ROOT, "%.4f", coverage))
              .append('\t').append(ticksFor(k, 0.50)).append('\t')
              .append(ticksFor(k, 0.80)).append('\t')
              .append(ticksFor(k, 0.95)).append('\n');
        }
        return sb.toString();
    }

    private static long ticksFor(double k, double targetCoverage) {
        if (k <= 0) {
            return -1;
        }
        return Math.round(-Math.log(1 - targetCoverage) / k);
    }

    private static String formatReport(List<CandidateStats> results) {
        StringBuilder sb = new StringBuilder();
        sb.append("nLga\tbinWidthDeg\tminWidthDeg\tmedianWidthDeg\tmeanWidthDeg\tmaxWidthDeg\tminWidthOverBin\tmedianWidthOverBin\ttranscriptionErrorRate\tfractionNarrowerThanBin\tatlasRows\tgenerationWallMs\tdynamicCoverageFraction\n");
        for (CandidateStats s : results) {
            double toDeg = 180.0 / Math.PI;
            sb.append(String.format(Locale.ROOT,
                                     "%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.6f\t%.4f\t%d\t%d\t%.4f%n",
                                     s.nLga(), s.binWidthRadians() * toDeg,
                                     s.minWidthRadians() * toDeg,
                                     s.medianWidthRadians() * toDeg,
                                     s.meanWidthRadians() * toDeg,
                                     s.maxWidthRadians() * toDeg,
                                     s.minWidthRadians()
                                     / s.binWidthRadians(),
                                     s.medianWidthRadians()
                                     / s.binWidthRadians(),
                                     s.transcriptionErrorRate(),
                                     s.fractionNarrowerThanBin(),
                                     s.atlasRows(), s.generationWallMs(),
                                     s.dynamicCoverageFraction()));
        }
        return sb.toString();
    }

    /**
     * Runs the full campaign: discover combos (via {@link
     * ContactComboCache#combosFor}), measure fine-grid widths once, then
     * analyze every candidate in {@link #CANDIDATES} against that shared
     * measurement. Prints a TSV report to stdout and writes the same
     * report under {@code target/} - never {@code src/test/resources/lga/}
     * (this class never picks or commits a final N_lga; see class
     * Javadoc).
     */
    public static void main(String[] args) throws IOException {
        ContactPredicate predicate = newPredicate();

        long discoverStart = System.nanoTime();
        List<Combo> combos = ContactComboCache.combosFor(predicate,
                                                           GEOMETRY_RESOLUTION,
                                                           RADIUS);
        long discoverMs = (System.nanoTime() - discoverStart) / 1_000_000;
        System.out.println("Loaded " + combos.size()
                            + " ever-contacting combos via ContactComboCache in "
                            + discoverMs
                            + "ms (cache-backed resource load when (geometryResolution="
                            + GEOMETRY_RESOLUTION + ", memberRadius=" + RADIUS
                            + ") matches the checked-in header; falls back to a live exhaustive "
                            + GEOMETRY_RESOLUTION + "x" + GEOMETRY_RESOLUTION
                            + " sweep otherwise - the native MemberGeometry LUT ceiling)");

        long widthStart = System.nanoTime();
        List<ComboWidth> widths = new ArrayList<>(combos.size());
        for (Combo combo : combos) {
            widths.add(measureWidth(predicate, combo));
        }
        long widthMs = (System.nanoTime() - widthStart) / 1_000_000;
        System.out.println("Measured fine-grid (" + FINE_STEPS + "x"
                            + FINE_STEPS + ") widths for " + widths.size()
                            + " combos in " + widthMs + "ms");

        List<CandidateStats> results = new ArrayList<>(CANDIDATES.length);
        for (int nLga : CANDIDATES) {
            long candidateStart = System.nanoTime();
            CandidateStats stats = analyze(nLga, widths);
            long candidateMs = (System.nanoTime() - candidateStart)
                                / 1_000_000;
            System.out.println("nLga=" + nLga + " analyzed in "
                                + candidateMs + "ms");
            results.add(stats);
        }

        String report = formatReport(results);
        String ticksAdvice = formatTicksAdvice(results);
        System.out.println();
        System.out.println(report);
        System.out.println(ticksAdvice);

        Path out = Path.of("target", "nlga-candidate-campaign.tsv");
        Files.createDirectories(out.getParent());
        Files.writeString(out, report + "\n" + ticksAdvice);
        System.out.println("Report written to " + out);
    }
}
