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
import java.util.List;
import java.util.Locale;
import java.util.Random;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.QuantaField;
import com.chiralbehaviors.inviscid.lga.CollisionSweep;
import com.chiralbehaviors.inviscid.lga.CollisionTable;
import com.chiralbehaviors.inviscid.lga.CollisionTable.SemiDetailedBalanceReport;
import com.chiralbehaviors.inviscid.lga.ContactAtlas;
import com.chiralbehaviors.inviscid.lga.ContactPredicate;
import com.chiralbehaviors.inviscid.lga.ContactScan;
import com.chiralbehaviors.inviscid.lga.FccNeighborhood;
import com.chiralbehaviors.inviscid.lga.HybridAutomaton;
import com.chiralbehaviors.inviscid.lga.LatticeGasAutomaton;
import com.chiralbehaviors.inviscid.lga.MemberGeometry;
import com.chiralbehaviors.inviscid.lga.QuantaExchangeRule;
import com.chiralbehaviors.inviscid.lga.TickDriver;
import com.chiralbehaviors.inviscid.measure.AnisotropyProbe.PooledResult;
import com.chiralbehaviors.inviscid.measure.AnisotropyProbe.SeedResult;
import com.chiralbehaviors.inviscid.measure.SpectrumAnalyzer.WindowFunction;

/**
 * C.5 (bead inviscid-0nx.22): re-runs the Phase B instrument suite against
 * {@link LatticeGasAutomaton} and produces the three-way K=0 / hybrid / LGA
 * comparison plus the isotropy verdict DATA (not a posture selection -- see
 * class Javadoc "Isotropy verdict" below and the committed report's own
 * POSTURE section).
 *
 * <h2>Five sub-measurements, each reusing an already-reviewed instrument
 * rather than reinventing one</h2>
 * <ol>
 * <li><b>Anisotropy</b> -- {@link AnisotropyProbe#runCampaign(Point3i, long[], int, int, SubstrateFactory)}
 * with a new LGA {@link SubstrateFactory}, at the IDENTICAL parameters
 * {@link AnisotropyProbe} used for Phase A ({@code DEFAULT_EXTENT/SEEDS/TICKS/
 * PACKET_QUANTA}) -- genuine apples-to-apples, both estimators, both the
 * naive-per-seed diagnostic and the pooled/null-calibrated significance
 * statistic, exactly as Phase A reports them (per-seed per-direction
 * {@code LGA_DIRECTION} rows for BOTH estimators, the naive-ratio {@code
 * LGA_SUMMARY} row, alongside the pooled rows) -- PLUS a winning-direction-
 * stability diagnostic Phase A's own report does not carry (T3 {@code
 * review-pattern-order-statistic-bias-max-min-ratio-2026-08-08}): which
 * direction has the largest TRANSPORT magnitude in each seed, and how
 * often that winner repeats across seeds. A STABLE winner across seeds is
 * evidence consistent with a real, direction-linked effect; a winner that
 * flips seed-to-seed is the signature of order-statistic noise -- see
 * {@code appendWinningDirectionSection}. This is also where the
 * structure-factor "ridge or absence" verdict (bead's named test 2) comes
 * from: the SPECTRAL estimator's per-direction ridge slopes ARE the ridge
 * data (see {@link AnisotropyProbe}'s own "SPECTRAL estimator" javadoc
 * section) -- a nonzero slope is a ridge, an all-zero slope is the
 * already-precedented diffusive signature, and either way this class states
 * which, never leaving the finding ambiguous.
 * <li><b>Spectral broadening</b> (three-way K=0 -&gt; hybrid -&gt; LGA) --
 * K=0 reuses {@link BaselineSpectrumHarness#run()} in-process (bit-identical
 * to the committed golden per {@code BaselineK0SpectrumTest}); hybrid and
 * LGA are driven from a SHARED initial condition (mirrors {@code
 * HybridVsLgaConsistencyTest}'s pattern) at the SAME cadence ({@link
 * SpectralCadence#aligned}'s stride, {@code phaseResolution=3600}) as the
 * K=0 baseline -- the {@code fftLength} is smaller purely for wall-time
 * budget (stated explicitly in the report header, never silently), which is
 * a sample-COUNT choice, not a resampling of already-collected data (user
 * decision D1/D2's "resampling PROHIBITED" is about not stretching/
 * compressing a recorded series, not about choosing how many aligned
 * samples to collect up front).
 * <li><b>Conservation</b> -- every driver in every sub-measurement below is
 * wrapped in a STRICT {@link ConservationAudit} (constructor arg {@code
 * true}), so any violation throws immediately rather than being checked
 * after the fact; ticks audited and the zero-violation count are recorded
 * in the report, exact.
 * <li><b>Collision-statistics comparability + gap mechanism + equilibrium
 * characterization</b> -- ONE shared long run reusing {@code
 * HybridVsLgaConsistencyTest}'s EXACT configuration (extent (4,4,4), seed
 * 42L, quanta bound 6, 2000 ticks) so this report's numbers cross-check
 * directly against that already-reviewed test's pinned 26.9/23.9
 * collision-rate values. Windowed rates surface the boot-transient
 * question (bead .21 critique carryover item 3); all 5 {@link
 * CollisionStatistics} fields are compared field-by-field with an explicit
 * comparable/not-verified verdict per field (bead .21 critique carryover
 * item 2) rather than a blanket claim; the no-op-fraction trend is the
 * dynamical corroboration of {@link CollisionTable#checkSemiDetailedBalance}'s
 * already-proven closed-form preimage profile (3:2:1 at output-diff
 * 0/&plusmn;1/else) -- the formal parity-floor fact, reported alongside its
 * empirical corroboration, per the SDB equilibrium caveat (user decision
 * 6): equilibrium statistics here do NOT assume semi-detailed balance --
 * the report states the accepted violation directly.
 * </ol>
 *
 * <h2>Isotropy verdict -- reporting rules (bead's own text, restated)</h2>
 * This class MEASURES. {@link #toTsv} writes the measured anisotropy ratio
 * and CI for each estimator, then the THREE POSTURES from the locked design
 * (accept/characterize; FCHC-style projection; member-orientational-state-
 * restores-isotropy) each WITH the evidence this campaign produced for it --
 * and explicitly does NOT select one. A caller reading the committed
 * artifact sees numbers and an escalation note, not a conclusion.
 *
 * @author halhildebrand
 */
public final class PhaseCMeasurement {

    // ------------------------------------------------------------------
    // Shared fixtures.
    // ------------------------------------------------------------------

    private static final String ATLAS_RELATIVE_PATH  = "src/test/resources/lga/contact-atlas-v2.tsv";
    private static final String REPORT_RELATIVE_PATH = "src/test/resources/lga/measurement-report-phaseC.tsv";
    private static final String PHASE_A_RELATIVE_PATH = "src/test/resources/lga/anisotropy-report-phaseA.tsv";

    // Sub-measurement 2/4: spectral broadening cadence -- SAME stride /
    // phaseResolution as the committed K=0 golden (BaselineSpectrumHarness),
    // fftLength reduced purely for wall-time (stated in the report).
    public static final Point3i SPECTRAL_EXTENT      = new Point3i(4, 4, 4);
    public static final long    SPECTRAL_SEED         = 42L;
    public static final int     SPECTRAL_MIN_QUANTA   = -5;
    public static final int     SPECTRAL_MAX_QUANTA   = 5;
    public static final int     SPECTRAL_STRIDE       = BaselineSpectrumHarness.STRIDE;
    public static final int     SPECTRAL_FFT_LENGTH   = 32;

    // Sub-measurement 5: collision-statistics / gap-mechanism / equilibrium
    // long run -- IDENTICAL config to HybridVsLgaConsistencyTest, so this
    // report's numbers directly cross-check against that already-reviewed
    // test's pinned 26.9/23.9 rate values.
    public static final Point3i LONG_RUN_EXTENT       = new Point3i(4, 4, 4);
    public static final long    LONG_RUN_SEED          = 42L;
    public static final int     LONG_RUN_QUANTA_BOUND  = 6;
    public static final int     LONG_RUN_TICKS         = 2000;
    public static final int     LONG_RUN_WINDOW        = 500;
    public static final long    SDB_WINDOW             = 10L;

    private PhaseCMeasurement() {
    }

    // ------------------------------------------------------------------
    // Record types.
    // ------------------------------------------------------------------

    /**
     * One substrate's aggregate spectral-broadening summary. Carries BOTH
     * linewidth conventions (user decision D3, pre-registered): {@code
     * meanAbsoluteLinewidthRadPerTick} (half-power/FWHM full-width, radians
     * per tick, via {@link SpectralCadence#omegaRadPerTick}) and {@code
     * meanFractionalLinewidth} (that width divided by the member's own
     * peak-bin center frequency -- dimensionless, a Q-factor-style
     * relative measure). NEITHER is promoted to sole headline -- both are
     * reported side by side, exactly as D3 requires. {@code
     * nMembersFractionalDefined} is strictly {@code <= nMembersChecked}: a
     * member whose peak bin is DC (zero center frequency) has an
     * undefined fractional linewidth (division by zero) and is excluded
     * from the fractional mean, not silently coerced to zero/NaN-then-
     * averaged-in -- the non-fabrication convention this codebase applies
     * throughout ({@link AnisotropyProbe}'s {@code RATIO_DEGENERATE_EPSILON}
     * pattern).
     */
    public record SpectralSummary(String substrate, Point3i extent, int fftLength,
                                   int stride, int nMembersChecked,
                                   double meanPeakFraction, double minPeakFraction,
                                   double maxPeakFraction, double meanEntropy,
                                   double meanAbsoluteLinewidthRadPerTick,
                                   double meanFractionalLinewidth,
                                   int nMembersFractionalDefined) {
    }

    /** One member's spectral linewidth, both conventions -- see {@link SpectralSummary}. */
    private record Linewidth(double absoluteRadPerTick, java.util.OptionalDouble fractional) {
    }

    /** One 500-tick window's collision rate/effective-ratio for both substrates. */
    public record CollisionWindow(int windowStart, int windowEnd,
                                   long hybridCollisions, long hybridEffective,
                                   long lgaCollisions, long lgaEffective) {

        public double hybridRatePerTick() {
            return (double) hybridCollisions / (windowEnd - windowStart);
        }

        public double lgaRatePerTick() {
            return (double) lgaCollisions / (windowEnd - windowStart);
        }

        public double hybridEffectiveFraction() {
            return hybridCollisions == 0 ? Double.NaN
                                          : (double) hybridEffective
                                            / hybridCollisions;
        }

        public double lgaEffectiveFraction() {
            return lgaCollisions == 0 ? Double.NaN
                                       : (double) lgaEffective / lgaCollisions;
        }
    }

    /** Field-by-field {@link CollisionStatistics} comparability verdict. */
    public record FieldComparison(String field, String hybridValue,
                                   String lgaValue, boolean verifiedComparable,
                                   String note) {
    }

    /** Quanta-value population summary (distribution-shape data point). */
    public record QuantaHistogramSummary(String label, long n, double mean,
                                          double variance, long min, long max) {
    }

    /**
     * Critic Significant #3 / relay item 5: within-window collision-rate
     * homogeneity for the anisotropy campaign's OWN 128-tick window
     * (summed across all 8 seeds), tested against the same
     * SDB-settling-decay signature the SEPARATE long run's {@link
     * CollisionWindow} series checks. {@code totalCollisionsAcrossSeeds}
     * is TOTAL (not effective-only) collisions -- {@link
     * CollisionStatistics} does not track effective-vs-no-op per tick,
     * and adding that would be a non-additive change to shared,
     * already-locked production code, out of scope for this diagnostic
     * (see {@code computeAnisotropyWithinWindowHomogeneity}'s Javadoc).
     */
    public record TickQuartileHomogeneity(int quartileIndex, int tickStart,
                                            int tickEndExclusive,
                                            long totalCollisionsAcrossSeeds) {
    }

    public record Report(AnisotropyProbe.Report lgaAnisotropy,
                          SpectralSummary k0Spectral,
                          SpectralSummary hybridSpectral,
                          SpectralSummary lgaSpectral,
                          List<CollisionWindow> windows,
                          List<TickQuartileHomogeneity> anisotropyWithinWindowHomogeneity,
                          List<FieldComparison> fieldComparisons,
                          QuantaHistogramSummary quantaBefore,
                          QuantaHistogramSummary quantaAfter,
                          SemiDetailedBalanceReport sdb,
                          long conservationTicksAudited,
                          long conservationViolations,
                          ContactAtlas.Header atlasHeader) {
    }

    // ------------------------------------------------------------------
    // LGA substrate factory (anisotropy campaign reuse) -- follows {@link
    // SubstrateFactory}'s RNG draw-ORDER contract (that interface's own
    // Javadoc): the SAME draw SEQUENCE as {@code
    // AnisotropyProbe.phaseAHybridSubstrate} -- random field state, THEN
    // the localized packet -- using substrate-APPROPRIATE types (this
    // factory draws {@code int} phase steps via {@code nextInt(3600)};
    // the hybrid factory draws {@code float} angles via {@code
    // nextFloat()*2*PI}). This is NOT a bit-identical RNG stream between
    // the two substrates (their per-draw consumption differs by type),
    // only the same NUMBER and ORDER of logical draws -- sufficient for
    // the contract's own reproducibility purpose (a substrate's OWN seed
    // reproduces its OWN trajectory), not a claim of cross-substrate
    // trajectory equality.
    // ------------------------------------------------------------------

    /**
     * @param perSeedStatsSink OPTIONAL (nullable) side-channel: if
     *                         non-null, every {@link CollisionStatistics}
     *                         this factory constructs is ALSO stashed
     *                         here keyed by seed, purely for POST-campaign
     *                         observation (e.g. {@link
     *                         CollisionStatistics#collisionsPerTick()}
     *                         within-window homogeneity checks) -- this is
     *                         additive instrumentation only: it changes
     *                         nothing about what {@link AnisotropyProbe}
     *                         computes or how, it only lets a caller keep
     *                         a reference to state {@link AnisotropyProbe}
     *                         already builds and mutates internally but
     *                         does not itself expose per-tick.
     */
    private static SubstrateFactory lgaSubstrateFactory(ContactAtlas atlas,
                                                          CollisionTable collisions,
                                                          java.util.Map<Long, CollisionStatistics> perSeedStatsSink) {
        return (extent, seed, packetQuanta, originCell) -> {
            CollisionStatistics statistics = new CollisionStatistics();
            if (perSeedStatsSink != null) {
                perSeedStatsSink.put(seed, statistics);
            }
            LatticeGasAutomaton lga = new LatticeGasAutomaton(extent, atlas,
                                                                collisions,
                                                                statistics);
            seedRandomPhases(lga, extent, seed);
            seedPacket(lga, originCell, packetQuanta);
            ConservationAudit audit = new ConservationAudit(lga, true);
            AuditedRun run = new AuditedRun(lga, audit);
            return new SubstrateFactory.Substrate(lga, run, lga.statistics());
        };
    }

    private static void seedRandomPhases(LatticeGasAutomaton lga, Point3i extent,
                                          long seed) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        int[] phases = new int[length];
        for (int i = 0; i < length; i++) {
            phases[i] = random.nextInt(3600);
        }
        lga.process((phase, quanta) -> System.arraycopy(phases, 0, phase, 0,
                                                          length));
    }

    private static void seedPacket(LatticeGasAutomaton lga, Point3i originCell,
                                    int packetQuanta) {
        int base = lga.indexOfCell(originCell);
        lga.process((phase, quanta) -> {
            for (int m = 0; m < 30; m++) {
                quanta[base + m] = packetQuanta;
            }
        });
    }

    // ------------------------------------------------------------------
    // Sub-measurement 2/4: spectral broadening, three-way.
    // ------------------------------------------------------------------

    private record SharedInitialCondition(int[] phase, long[] quanta) {
    }

    private static SharedInitialCondition sharedInitialCondition(Point3i extent,
                                                                    long seed,
                                                                    int minQuanta,
                                                                    int maxQuanta) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        int range = maxQuanta - minQuanta + 1;
        int[] phase = new int[length];
        long[] quanta = new long[length];
        for (int i = 0; i < length; i++) {
            phase[i] = random.nextInt(3600);
            quanta[i] = random.nextInt(range) + minQuanta;
        }
        return new SharedInitialCondition(phase, quanta);
    }

    private static HybridAutomaton newHybrid(ContactAtlas atlas,
                                              SharedInitialCondition ic,
                                              CollisionStatistics statistics,
                                              Point3i extent) {
        Necronomata automaton = new Necronomata(extent);
        int length = ic.phase().length;
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < length; i++) {
                angle[i] = (float) (2.0 * Math.PI * ic.phase()[i] / 3600.0);
                frequency[i] = ic.quanta()[i];
            }
        });
        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactPredicate predicate = new ContactPredicate(new MemberGeometry(atlas.header()
                                                                                    .geometryResolution(),
                                                                               atlas.header()
                                                                                    .memberRadius()));
        ContactScan scan = new ContactScan(automaton, neighborhood, predicate);
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        return new HybridAutomaton(automaton, sweep);
    }

    private static LatticeGasAutomaton newLga(ContactAtlas atlas,
                                               CollisionTable collisions,
                                               SharedInitialCondition ic,
                                               CollisionStatistics statistics,
                                               Point3i extent) {
        LatticeGasAutomaton lga = new LatticeGasAutomaton(extent, atlas,
                                                            collisions,
                                                            statistics);
        lga.process((phase, quanta) -> {
            System.arraycopy(ic.phase(), 0, phase, 0, ic.phase().length);
            System.arraycopy(ic.quanta(), 0, quanta, 0, ic.quanta().length);
        });
        return lga;
    }

    /** Every even-parity-cell member's flat slot index. */
    private static int[] activeMemberIndices(QuantaField field) {
        List<Integer> indices = new ArrayList<>();
        field.forEachCell(cell -> {
            int base = field.indexOfCell(cell);
            for (int local = 0; local < 30; local++) {
                indices.add(base + local);
            }
        });
        int[] result = new int[indices.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = indices.get(i);
        }
        return result;
    }

    /**
     * Records {@code samples} phase snapshots for every member in {@code
     * memberIndices}, advancing the substrate {@code stride} ticks between
     * samples via {@code advanceOneTick} (a single synchronized pass --
     * mirrors {@code BaselineSpectrumHarness.recordStridedAngleSeries}
     * generalized over {@link QuantaField} instead of {@code Necronomata},
     * and over an arbitrary advance callback instead of {@code
     * Necronomata::step} specifically, so K0/hybrid/LGA all share this one
     * recorder).
     */
    private static float[][] recordStridedPhaseSeries(Runnable advanceOneTick,
                                                        QuantaField field,
                                                        int[] memberIndices,
                                                        int samples, int stride) {
        float[][] series = new float[memberIndices.length][samples];
        for (int s = 0; s < samples; s++) {
            for (int m = 0; m < memberIndices.length; m++) {
                series[m][s] = field.phaseAt(memberIndices[m]);
            }
            if (s < samples - 1) {
                for (int t = 0; t < stride; t++) {
                    advanceOneTick.run();
                }
            }
        }
        return series;
    }

    /**
     * Half-power (FWHM) linewidth of {@code power}'s dominant peak, BOTH
     * conventions (user decision D3) -- see {@link SpectralSummary}'s
     * Javadoc. {@code cadence} converts the bin-domain full-width and the
     * peak bin's own center frequency to radians/tick.
     */
    private static Linewidth computeLinewidth(double[] power, int n,
                                                SpectralCadence cadence) {
        int peak = SpectrumAnalyzer.peakBin(power);
        int fwhmBins = stepsToHalfPower(power, peak, 1, n)
                        + stepsToHalfPower(power, peak, -1, n);
        double absolute = cadence.omegaRadPerTick(2 * Math.PI * fwhmBins / n);
        double centerOmega = cadence.omegaRadPerTick(2 * Math.PI
                                                       * signedIndex(peak, n)
                                                       / n);
        java.util.OptionalDouble fractional = Math.abs(centerOmega) > 1e-12
                                               ? java.util.OptionalDouble.of(absolute
                                                                              / Math.abs(centerOmega))
                                               : java.util.OptionalDouble.empty();
        return new Linewidth(absolute, fractional);
    }

    /**
     * Steps from {@code peak} in {@code direction} (+1 or -1, circular over
     * {@code n} bins) until the first bin below half the peak's power --
     * returns that step count, or {@code n} if the walk returns to {@code
     * peak} without ever dropping below half power (a degenerate,
     * essentially-flat spectrum).
     */
    private static int stepsToHalfPower(double[] power, int peak, int direction,
                                          int n) {
        double halfPower = power[peak] / 2.0;
        int idx = peak;
        for (int steps = 1; steps <= n; steps++) {
            idx = Math.floorMod(idx + direction, n);
            if (idx == peak) {
                return n;
            }
            if (power[idx] < halfPower) {
                return steps;
            }
        }
        return n;
    }

    /** Mirrors {@code StructureFactor}'s private signed-bin convention exactly. */
    private static int signedIndex(int idx, int period) {
        return idx > period / 2 ? idx - period : idx;
    }

    private record SpectralAggregate(double sumFraction, double minFraction,
                                       double maxFraction, double sumEntropy,
                                       double sumAbsoluteLinewidth,
                                       double sumFractionalLinewidth,
                                       int nFractionalDefined, int n) {
    }

    private static SpectralAggregate aggregate(float[][] series, int n,
                                                 SpectralCadence cadence) {
        double sumFraction = 0;
        double minFraction = Double.POSITIVE_INFINITY;
        double maxFraction = Double.NEGATIVE_INFINITY;
        double sumEntropy = 0;
        double sumAbsoluteLinewidth = 0;
        double sumFractionalLinewidth = 0;
        int nFractionalDefined = 0;
        int counted = 0;
        for (float[] memberSeries : series) {
            double[] power = SpectrumAnalyzer.powerSpectrum(memberSeries,
                                                              WindowFunction.RECTANGULAR);
            int peakBin = SpectrumAnalyzer.peakBin(power);
            double total = 0;
            for (double p : power) {
                total += p;
            }
            if (total <= 0.0) {
                continue;
            }
            double peakFraction = power[peakBin] / total;
            double entropy = 0.0;
            for (double p : power) {
                if (p <= 0.0) {
                    continue;
                }
                double pi = p / total;
                entropy -= pi * Math.log(pi);
            }
            Linewidth lw = computeLinewidth(power, n, cadence);
            sumFraction += peakFraction;
            minFraction = Math.min(minFraction, peakFraction);
            maxFraction = Math.max(maxFraction, peakFraction);
            sumEntropy += entropy;
            sumAbsoluteLinewidth += lw.absoluteRadPerTick();
            if (lw.fractional().isPresent()) {
                sumFractionalLinewidth += lw.fractional().getAsDouble();
                nFractionalDefined++;
            }
            counted++;
        }
        return new SpectralAggregate(sumFraction, minFraction, maxFraction,
                                      sumEntropy, sumAbsoluteLinewidth,
                                      sumFractionalLinewidth,
                                      nFractionalDefined, counted);
    }

    private static SpectralSummary summarizeK0() {
        BaselineSpectrumHarness.Result officialResult = BaselineSpectrumHarness.run();

        // Recompute the SAME deterministic K0 recipe locally (extent,
        // seed, stride, fftLength, quanta range all identical to
        // BaselineSpectrumHarness.run()'s documented defaults) purely to
        // get the raw power spectra linewidth needs -- BaselineSpectrumHarness
        // .MemberSpectrum intentionally carries only summary stats (peakFraction
        // /entropy), not the raw array, per its own "keeps the artifact
        // small and diffable" convention, so this class cannot read
        // linewidth off it directly.
        Point3i extent = BaselineSpectrumHarness.DEFAULT_EXTENT;
        Necronomata automaton = new Necronomata(extent);
        int[] memberIndices = activeMemberIndices(automaton);
        java.util.Random random = new java.util.Random(BaselineSpectrumHarness.DEFAULT_SEED);
        int range = BaselineSpectrumHarness.MAX_QUANTA
                     - BaselineSpectrumHarness.MIN_QUANTA + 1;
        float[] quanta = new float[memberIndices.length];
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int m = 0; m < memberIndices.length; m++) {
                int q = random.nextInt(range) + BaselineSpectrumHarness.MIN_QUANTA;
                frequency[memberIndices[m]] = q;
                quanta[m] = q;
            }
        });
        ConservationAudit audit = new ConservationAudit(automaton, true);
        int[] k0TickCounter = { 0 };
        float[][] series = recordStridedPhaseSeries(() -> {
            automaton.step();
            audit.auditTick(k0TickCounter[0]++);
        }, automaton, memberIndices, BaselineSpectrumHarness.FFT_LENGTH,
                                                       BaselineSpectrumHarness.STRIDE);

        // Filter to nonzero-seeded quanta only, matching
        // BaselineSpectrumHarness's own row-count-reconciliation convention
        // (a zero-quanta member is a trivial non-rotating DC case).
        int nonzeroCount = 0;
        for (float q : quanta) {
            if (q != 0f) {
                nonzeroCount++;
            }
        }
        float[][] filteredSeries = new float[nonzeroCount][];
        int fi = 0;
        for (int m = 0; m < series.length; m++) {
            if (quanta[m] != 0f) {
                filteredSeries[fi++] = series[m];
            }
        }

        SpectralCadence cadence = new SpectralCadence(automaton.phaseResolution(),
                                                        BaselineSpectrumHarness.STRIDE);
        SpectralAggregate agg = aggregate(filteredSeries,
                                            BaselineSpectrumHarness.FFT_LENGTH,
                                            cadence);

        // Cross-check: this locally-recomputed peakFraction distribution
        // must agree with the official, already-golden-verified
        // BaselineSpectrumHarness result (same seed/params -- must be
        // bit-identical modulo the empty-window edge case).
        double officialSumFraction = 0;
        int officialN = 0;
        for (BaselineSpectrumHarness.MemberSpectrum m : officialResult.members) {
            if (m.quanta == 0f) {
                continue;
            }
            officialSumFraction += m.peakFraction;
            officialN++;
        }
        if (officialN != agg.n() || Math.abs(officialSumFraction / officialN
                                              - agg.sumFraction() / agg.n()) > 1e-6) {
            throw new IllegalStateException("K0 local recomputation diverged from BaselineSpectrumHarness.run() -- "
                                             + "officialN=" + officialN + " localN=" + agg.n()
                                             + " officialMeanFraction=" + (officialSumFraction / officialN)
                                             + " localMeanFraction=" + (agg.sumFraction() / agg.n()));
        }

        return new SpectralSummary("K0", extent, BaselineSpectrumHarness.FFT_LENGTH,
                                    BaselineSpectrumHarness.STRIDE, agg.n(),
                                    agg.sumFraction() / agg.n(), agg.minFraction(),
                                    agg.maxFraction(), agg.sumEntropy() / agg.n(),
                                    agg.sumAbsoluteLinewidth() / agg.n(),
                                    agg.nFractionalDefined() == 0 ? Double.NaN
                                                                   : agg.sumFractionalLinewidth()
                                                                     / agg.nFractionalDefined(),
                                    agg.nFractionalDefined());
    }

    private static SpectralSummary summarizeCollisionBearing(String label,
                                                                AuditedRun run,
                                                                QuantaField field,
                                                                Point3i extent,
                                                                int[] tickCounter) {
        int[] memberIndices = activeMemberIndices(field);
        float[][] series = recordStridedPhaseSeries(() -> run.tick(tickCounter[0]++),
                                                       field, memberIndices,
                                                       SPECTRAL_FFT_LENGTH,
                                                       SPECTRAL_STRIDE);
        SpectralCadence cadence = new SpectralCadence(field.phaseResolution(),
                                                        SPECTRAL_STRIDE);
        SpectralAggregate agg = aggregate(series, SPECTRAL_FFT_LENGTH, cadence);
        return new SpectralSummary(label, extent, SPECTRAL_FFT_LENGTH,
                                    SPECTRAL_STRIDE, agg.n(), agg.sumFraction()
                                                               / agg.n(),
                                    agg.minFraction(), agg.maxFraction(),
                                    agg.sumEntropy() / agg.n(),
                                    agg.sumAbsoluteLinewidth() / agg.n(),
                                    agg.nFractionalDefined() == 0 ? Double.NaN
                                                                   : agg.sumFractionalLinewidth()
                                                                     / agg.nFractionalDefined(),
                                    agg.nFractionalDefined());
    }

    // ------------------------------------------------------------------
    // Sub-measurement 5: collision-statistics comparability, gap
    // mechanism, boot-transient, equilibrium characterization.
    // ------------------------------------------------------------------

    private static QuantaHistogramSummary histogramOf(String label, long[] quanta) {
        long n = quanta.length;
        double mean = 0;
        long min = Long.MAX_VALUE;
        long max = Long.MIN_VALUE;
        for (long q : quanta) {
            mean += q;
            min = Math.min(min, q);
            max = Math.max(max, q);
        }
        mean /= n;
        double variance = 0;
        for (long q : quanta) {
            double d = q - mean;
            variance += d * d;
        }
        variance /= n;
        return new QuantaHistogramSummary(label, n, mean, variance, min, max);
    }

    private static long[] snapshotQuanta(LatticeGasAutomaton lga) {
        long[][] box = new long[1][];
        lga.process((phase, quanta) -> box[0] = quanta.clone());
        return box[0];
    }

    /**
     * Critic Significant #3 / relay item 5: splits the anisotropy
     * campaign's own {@code (ticks-1)}-tick window (127 actual {@code
     * tick()} calls per seed, tick args {@code 0..126}) into 4 roughly-
     * equal quartiles and sums {@link CollisionStatistics#collisionsPerTick()}
     * (TOTAL, not effective-only -- see {@link TickQuartileHomogeneity}'s
     * Javadoc) across all 8 seeds within each quartile. Tests whether the
     * SAME within-run rate decline {@link CollisionWindow}'s series shows
     * for the separate 2000-tick long run also shows up here, at the much
     * shorter 128-tick scale the OLS transport-rate fit actually uses --
     * additive-only instrumentation: {@code perSeedStats} is populated by
     * {@link #lgaSubstrateFactory}'s side channel, so this reads data
     * {@link AnisotropyProbe} already computes and mutates internally
     * without any change to that class.
     */
    private static List<TickQuartileHomogeneity> computeAnisotropyWithinWindowHomogeneity(java.util.Map<Long, CollisionStatistics> perSeedStats,
                                                                                             int ticks) {
        int actualTickCount = ticks - 1;
        int quartiles = 4;
        List<TickQuartileHomogeneity> result = new ArrayList<>();
        for (int q = 0; q < quartiles; q++) {
            int start = (int) ((long) q * actualTickCount / quartiles);
            int end = q == quartiles - 1 ? actualTickCount
                                          : (int) ((long) (q + 1)
                                                    * actualTickCount
                                                    / quartiles);
            long sum = 0;
            for (CollisionStatistics stats : perSeedStats.values()) {
                for (java.util.Map.Entry<Integer, Long> e : stats.collisionsPerTick()
                                                                   .entrySet()) {
                    int tick = e.getKey();
                    if (tick >= start && tick < end) {
                        sum += e.getValue();
                    }
                }
            }
            result.add(new TickQuartileHomogeneity(q, start, end, sum));
        }
        return result;
    }

    // ------------------------------------------------------------------
    // Campaign assembly.
    // ------------------------------------------------------------------

    public static Report runCampaign() throws IOException {
        Path atlasPath = Paths.get(ATLAS_RELATIVE_PATH);
        ContactAtlas atlas = ContactAtlas.read(atlasPath);
        CollisionTable collisions = CollisionTable.buildFromPhaseARule(new QuantaExchangeRule());

        long conservationTicks = 0;
        long conservationViolations = 0;

        // --- Sub-measurement 1: LGA anisotropy, Phase-A-identical params ---
        // perSeedStatsSink: additive-only side channel (see
        // lgaSubstrateFactory's Javadoc) -- lets this class read
        // CollisionStatistics#collisionsPerTick() after the campaign
        // without AnisotropyProbe exposing it itself and without any
        // change to AnisotropyProbe.java.
        java.util.Map<Long, CollisionStatistics> perSeedStatsSink = new java.util.LinkedHashMap<>();
        AnisotropyProbe.Report lgaAnisotropy = AnisotropyProbe.runCampaign(AnisotropyProbe.DEFAULT_EXTENT,
                                                                             AnisotropyProbe.DEFAULT_SEEDS,
                                                                             AnisotropyProbe.DEFAULT_TICKS,
                                                                             AnisotropyProbe.DEFAULT_PACKET_QUANTA,
                                                                             lgaSubstrateFactory(atlas,
                                                                                                  collisions,
                                                                                                  perSeedStatsSink));
        // AnisotropyProbe.runOneSeed calls tick() exactly (ticks-1) times
        // per seed (tick 0 is the pre-run snapshot, no tick() call) -- see
        // that method's loop, "for (int tick = 1; tick < ticks; tick++)".
        conservationTicks += (long) AnisotropyProbe.DEFAULT_SEEDS.length
                              * (AnisotropyProbe.DEFAULT_TICKS - 1);
        List<TickQuartileHomogeneity> anisotropyHomogeneity = computeAnisotropyWithinWindowHomogeneity(perSeedStatsSink,
                                                                                                          AnisotropyProbe.DEFAULT_TICKS);

        // --- Sub-measurement 2/4: spectral broadening, three-way ---
        SpectralSummary k0 = summarizeK0();

        SharedInitialCondition spectralIc = sharedInitialCondition(SPECTRAL_EXTENT,
                                                                     SPECTRAL_SEED,
                                                                     SPECTRAL_MIN_QUANTA,
                                                                     SPECTRAL_MAX_QUANTA);
        CollisionStatistics hybridSpectralStats = new CollisionStatistics();
        HybridAutomaton hybridSpectral = newHybrid(atlas, spectralIc,
                                                     hybridSpectralStats,
                                                     SPECTRAL_EXTENT);
        AuditedRun hybridSpectralRun = new AuditedRun(hybridSpectral,
                                                        new ConservationAudit(hybridSpectral.automaton(),
                                                                               true));
        int[] hybridTickCounter = { 0 };
        SpectralSummary hybridSummary = summarizeCollisionBearing("HYBRID",
                                                                     hybridSpectralRun,
                                                                     hybridSpectral.field(),
                                                                     SPECTRAL_EXTENT,
                                                                     hybridTickCounter);
        conservationTicks += hybridTickCounter[0];

        CollisionStatistics lgaSpectralStats = new CollisionStatistics();
        LatticeGasAutomaton lgaSpectral = newLga(atlas, collisions, spectralIc,
                                                   lgaSpectralStats,
                                                   SPECTRAL_EXTENT);
        AuditedRun lgaSpectralRun = new AuditedRun(lgaSpectral,
                                                     new ConservationAudit(lgaSpectral,
                                                                            true));
        int[] lgaTickCounter = { 0 };
        SpectralSummary lgaSummary = summarizeCollisionBearing("LGA",
                                                                  lgaSpectralRun,
                                                                  lgaSpectral,
                                                                  SPECTRAL_EXTENT,
                                                                  lgaTickCounter);
        conservationTicks += lgaTickCounter[0];

        // --- Sub-measurement 5: collision-statistics / gap / equilibrium ---
        SharedInitialCondition longIc = sharedInitialCondition(LONG_RUN_EXTENT,
                                                                  LONG_RUN_SEED,
                                                                  -LONG_RUN_QUANTA_BOUND,
                                                                  LONG_RUN_QUANTA_BOUND);
        CollisionStatistics hybridLongStats = new CollisionStatistics();
        HybridAutomaton hybridLong = newHybrid(atlas, longIc, hybridLongStats,
                                                 LONG_RUN_EXTENT);
        AuditedRun hybridLongRun = new AuditedRun(hybridLong,
                                                    new ConservationAudit(hybridLong.automaton(),
                                                                           true));

        CollisionStatistics lgaLongStats = new CollisionStatistics();
        LatticeGasAutomaton lgaLong = newLga(atlas, collisions, longIc,
                                               lgaLongStats, LONG_RUN_EXTENT);
        AuditedRun lgaLongRun = new AuditedRun(lgaLong,
                                                 new ConservationAudit(lgaLong,
                                                                        true));

        long[] quantaBefore = snapshotQuanta(lgaLong);

        List<CollisionWindow> windows = new ArrayList<>();
        long prevHybridTotal = 0;
        long prevHybridEffective = 0;
        long prevLgaTotal = 0;
        long prevLgaEffective = 0;
        for (int windowStart = 0; windowStart < LONG_RUN_TICKS; windowStart += LONG_RUN_WINDOW) {
            int windowEnd = Math.min(windowStart + LONG_RUN_WINDOW,
                                      LONG_RUN_TICKS);
            for (int t = windowStart; t < windowEnd; t++) {
                hybridLongRun.tick(t);
                lgaLongRun.tick(t);
            }
            long hybridTotal = hybridLongStats.totalCollisions();
            long hybridEffective = hybridLongStats.effectiveCollisions();
            long lgaTotal = lgaLongStats.totalCollisions();
            long lgaEffective = lgaLongStats.effectiveCollisions();
            windows.add(new CollisionWindow(windowStart, windowEnd,
                                             hybridTotal - prevHybridTotal,
                                             hybridEffective - prevHybridEffective,
                                             lgaTotal - prevLgaTotal,
                                             lgaEffective - prevLgaEffective));
            prevHybridTotal = hybridTotal;
            prevHybridEffective = hybridEffective;
            prevLgaTotal = lgaTotal;
            prevLgaEffective = lgaEffective;
        }
        conservationTicks += 2L * LONG_RUN_TICKS;

        long[] quantaAfter = snapshotQuanta(lgaLong);

        List<FieldComparison> fieldComparisons = new ArrayList<>();
        fieldComparisons.add(new FieldComparison("totalCollisions",
                                                   Long.toString(hybridLongStats.totalCollisions()),
                                                   Long.toString(lgaLongStats.totalCollisions()),
                                                   true,
                                                   "verified comparable, bead .21 HybridVsLgaConsistencyTest test 8 (rate gap ~11%)"));
        fieldComparisons.add(new FieldComparison("effectiveCollisions",
                                                   Long.toString(hybridLongStats.effectiveCollisions()),
                                                   Long.toString(lgaLongStats.effectiveCollisions()),
                                                   true,
                                                   "verified comparable, bead .21 HybridVsLgaConsistencyTest test 8 (effective-ratio gap ~0.011)"));
        fieldComparisons.add(new FieldComparison("collisionsPerDirection",
                                                   hybridLongStats.collisionsPerDirection()
                                                                   .toString(),
                                                   lgaLongStats.collisionsPerDirection()
                                                                .toString(),
                                                   false,
                                                   "NOT independently verified comparable prior to this bead -- reported as raw data for the first time, not asserted equal"));
        fieldComparisons.add(new FieldComparison("transferMagnitudeHistogram",
                                                   hybridLongStats.transferMagnitudeHistogram()
                                                                   .toString(),
                                                   lgaLongStats.transferMagnitudeHistogram()
                                                                .toString(),
                                                   false,
                                                   "NOT independently verified comparable prior to this bead -- QuantaExchangeRule only ever transfers magnitude 0 or 1 on both substrates by construction, so this is expected trivially degenerate, not a genuine cross-substrate distributional check"));
        fieldComparisons.add(new FieldComparison("meanFreePathProxy",
                                                   Double.toString(hybridLongStats.meanFreePathProxy()),
                                                   Double.toString(lgaLongStats.meanFreePathProxy()),
                                                   false,
                                                   "meanFreePathProxy's own Javadoc: a coarse tick-span/collision-count proxy, not a physical mean free path -- comparable only as a coarse collision-cadence indicator, not a spatial-transport statistic"));

        SemiDetailedBalanceReport sdb = collisions.checkSemiDetailedBalance(SDB_WINDOW);

        // conservationViolations stays exactly 0: every driver above is
        // wrapped in a STRICT ConservationAudit, which throws immediately
        // on any violation -- reaching this line is itself the proof.
        return new Report(lgaAnisotropy, k0, hybridSummary, lgaSummary, windows,
                           anisotropyHomogeneity, fieldComparisons,
                           histogramOf("before", quantaBefore),
                           histogramOf("after", quantaAfter), sdb,
                           conservationTicks, conservationViolations,
                           atlas.header());
    }

    // ------------------------------------------------------------------
    // main() -- regenerates the committed artifact. Not run by surefire.
    // ------------------------------------------------------------------

    public static void main(String[] args) throws IOException {
        long start = System.nanoTime();
        Report report = runCampaign();
        double wallSeconds = (System.nanoTime() - start) / 1e9;
        String tsv = toTsv(report, resolveGitCommit());
        Path path = Paths.get(REPORT_RELATIVE_PATH);
        Files.createDirectories(path.getParent());
        Files.write(path, tsv.getBytes(StandardCharsets.UTF_8));
        System.out.println("Wrote " + path.toAbsolutePath() + " in "
                            + wallSeconds + "s");
        System.out.println("LGA anisotropy pooled TRANSPORT: "
                            + report.lgaAnisotropy().pooledTransport());
        System.out.println("LGA anisotropy pooled SPECTRAL: "
                            + report.lgaAnisotropy().pooledSpectral());
        System.out.println("K0 spectral: " + report.k0Spectral());
        System.out.println("HYBRID spectral: " + report.hybridSpectral());
        System.out.println("LGA spectral: " + report.lgaSpectral());
        System.out.println("conservationTicksAudited="
                            + report.conservationTicksAudited()
                            + " violations=" + report.conservationViolations());
    }

    // ------------------------------------------------------------------
    // TSV serialization.
    // ------------------------------------------------------------------

    static String toTsv(Report report, String gitCommit) {
        StringBuilder sb = new StringBuilder();
        sb.append("# PhaseCMeasurement report -- C.5 Phase C measurement re-run and isotropy verdict data\n");
        sb.append("# bead=inviscid-0nx.22\n");
        sb.append("# generator=").append(PhaseCMeasurement.class.getName())
          .append('\n');
        sb.append("# gitCommit=").append(gitCommit).append('\n');
        sb.append("# atlasVersion=").append(report.atlasHeader().atlasVersion())
          .append('\n');
        sb.append("# nLga=").append(report.atlasHeader().phaseResolutionNLga())
          .append('\n');
        sb.append("# subBinSteps=").append(report.atlasHeader().subBinSteps())
          .append('\n');
        sb.append("# phaseResolution=3600 (both substrates, cadence 2A -- see LatticeGasAutomaton class Javadoc)\n");
        sb.append("# anisotropyExtent=").append(AnisotropyProbe.DEFAULT_EXTENT.x)
          .append(',').append(AnisotropyProbe.DEFAULT_EXTENT.y).append(',')
          .append(AnisotropyProbe.DEFAULT_EXTENT.z)
          .append(" (IDENTICAL to committed anisotropy-report-phaseA.tsv, for genuine apples-to-apples)\n");
        StringBuilder seeds = new StringBuilder();
        for (long seed : AnisotropyProbe.DEFAULT_SEEDS) {
            if (seeds.length() > 0) {
                seeds.append(',');
            }
            seeds.append(seed);
        }
        sb.append("# anisotropySeeds=").append(seeds).append('\n');
        sb.append("# anisotropyTicks=").append(AnisotropyProbe.DEFAULT_TICKS)
          .append('\n');
        sb.append("# anisotropyPacketQuanta=")
          .append(AnisotropyProbe.DEFAULT_PACKET_QUANTA).append('\n');
        sb.append("# spectralExtent=").append(SPECTRAL_EXTENT.x).append(',')
          .append(SPECTRAL_EXTENT.y).append(',').append(SPECTRAL_EXTENT.z)
          .append('\n');
        sb.append("# spectralSeed=").append(SPECTRAL_SEED).append('\n');
        sb.append("# spectralStride=").append(SPECTRAL_STRIDE)
          .append(" (== BaselineSpectrumHarness.STRIDE, exact-alignment stride for phaseResolution=3600 -- SAME cadence as the K0 golden, see SpectralCadence)\n");
        sb.append("# spectralFftLength=").append(SPECTRAL_FFT_LENGTH)
          .append(" (K0 golden used 256; reduced here purely for wall-time budget -- a sample-COUNT choice at a MATCHED cadence, not a resampling of collected data; see class Javadoc)\n");
        // T2 analysis-73v-spectral-conversion-and-cadence.md §4.2's posture
        // 2-ii metadata requirement, extending acceptance criterion 5:
        // phaseResolution/subBinSteps/stride are already above;
        // nyquistQuantaBound and collisionOpportunitiesPerRevolution close
        // the list. Both substrates share phaseResolution=3600 under
        // cadence 2A (pinned), so these are substrate-independent.
        int nyquistQuantaBound = 3600 / (2 * SPECTRAL_STRIDE);
        sb.append("# nyquistQuantaBound=").append(nyquistQuantaBound)
          .append(" (phaseResolution/(2*spectralStride) -- seeded quanta range [")
          .append(SPECTRAL_MIN_QUANTA).append(',').append(SPECTRAL_MAX_QUANTA)
          .append("] is within this bound)\n");
        sb.append("# collisionOpportunitiesPerRevolutionAtLongRunQuantaBound=")
          .append(formatPrecise(3600.0 / LONG_RUN_QUANTA_BOUND))
          .append(" (P/q for q=longRunQuantaBound=").append(LONG_RUN_QUANTA_BOUND)
          .append(", the dimensionless per-revolution collision-opportunity count T2 73v §4.2 names -- one collision opportunity per tick on both substrates under cadence 2A, so this is directly comparable, no rescaling)\n");
        sb.append("# longRunExtent=").append(LONG_RUN_EXTENT.x).append(',')
          .append(LONG_RUN_EXTENT.y).append(',').append(LONG_RUN_EXTENT.z)
          .append('\n');
        sb.append("# longRunSeed=").append(LONG_RUN_SEED).append('\n');
        sb.append("# longRunQuantaBound=").append(LONG_RUN_QUANTA_BOUND)
          .append('\n');
        sb.append("# longRunTicks=").append(LONG_RUN_TICKS)
          .append(" (IDENTICAL config to HybridVsLgaConsistencyTest -- numbers cross-check against that test's pinned 26.9/23.9 rate values)\n");
        sb.append("# longRunWindowSize=").append(LONG_RUN_WINDOW).append('\n');
        sb.append("# conservationMode=STRICT (every driver, every sub-measurement)\n");
        sb.append("# conservationTicksAudited=").append(report.conservationTicksAudited())
          .append('\n');
        sb.append("# conservationViolations=").append(report.conservationViolations())
          .append(" (exact -- strict ConservationAudit throws immediately on any violation, so a nonzero value here is structurally impossible unless this campaign itself crashed)\n");
        sb.append("# precision=%.9e\n");
        sb.append('\n');

        appendAnisotropySection(sb, report.lgaAnisotropy());
        appendSpectralSection(sb, report);
        appendCollisionSection(sb, report);
        appendEquilibriumSection(sb, report);
        appendPostureSection(sb, report);

        return sb.toString();
    }

    private static void appendAnisotropySection(StringBuilder sb,
                                                  AnisotropyProbe.Report lga) {
        sb.append("# === SECTION 1: ANISOTROPY (LGA, both estimators; POOLED_SUMMARY is the significance statistic) ===\n");
        sb.append("# columns(LGA_DIRECTION)=recordType\tseed\testimator\tdirection\tmagnitude\tsampleSize\n");
        sb.append("# columns(LGA_SUMMARY)=recordType\testimator\tratio\tciLower\tciUpper\tnSeedsUsed\tnSeedsDegenerate (DIAGNOSTIC, naive per-seed ratio -- see naivePerSeedRatioCaveat)\n");
        sb.append("# columns(LGA_POOLED_DIRECTION)=recordType\testimator\tdirection\tmean\tciLower\tciUpper\tnSeeds\n");
        sb.append("# columns(LGA_POOLED_SUMMARY)=recordType\testimator\tpooledRatio\tpooledRatioCiLower\tpooledRatioCiUpper\tpermutationPValue\tpermutationNull95\tpermutationCount\n");
        sb.append("# columns(LGA_COLLISIONS)=recordType\tseed\ttotalCollisions\teffectiveCollisions\n");
        sb.append("# columns(WINNING_DIRECTION)=recordType\testimator\tseed\twinningDirection (argmax magnitude across X100/X110/X111; degenerate/tied estimators are noted, not silently tie-broken)\n");
        sb.append("# columns(WINNING_DIRECTION_STABILITY)=recordType\testimator\tmodeDirection\tmodeCount\ttotalSeeds\tstabilityFraction (T3 review-pattern-order-statistic-bias-max-min-ratio-2026-08-08: a STABLE winner across seeds is evidence consistent with a real direction-linked effect; a flipping winner is order-statistic noise)\n");
        sb.append("# naivePerSeedRatioCaveat=LGA_SUMMARY rows (mean of per-seed max/min ratios) are a DIAGNOSTIC, bounded below by 1.0 by construction, upward-biased by seed noise (order-statistic artifact, T3 critique-pattern-max-min-ratio-order-statistic-bias) -- the significance statistic is LGA_POOLED_SUMMARY's permutationPValue, see ciVsPermutationReconciliation below\n");

        for (SeedResult sr : lga.perSeed()) {
            appendDirectionRows(sb, sr.seed(), "TRANSPORT", sr.transport());
            appendDirectionRows(sb, sr.seed(), "SPECTRAL", sr.spectral());
        }
        appendNaiveSummaryRow(sb, "TRANSPORT", lga.transportCi());
        appendNaiveSummaryRow(sb, "SPECTRAL", lga.spectralCi());
        appendPooledDirectionRows(sb, "TRANSPORT", lga.pooledTransport());
        appendPooledDirectionRows(sb, "SPECTRAL", lga.pooledSpectral());
        appendPooledSummaryRow(sb, "TRANSPORT", lga.pooledTransport());
        appendPooledSummaryRow(sb, "SPECTRAL", lga.pooledSpectral());
        double meanEffective = 0;
        double meanTotal = 0;
        for (SeedResult sr : lga.perSeed()) {
            sb.append("LGA_COLLISIONS\t").append(sr.seed()).append('\t')
              .append(sr.totalCollisions()).append('\t')
              .append(sr.effectiveCollisions()).append('\n');
            meanEffective += sr.effectiveCollisions();
            meanTotal += sr.totalCollisions();
        }
        int nSeeds = lga.perSeed().size();
        meanEffective /= nSeeds;
        meanTotal /= nSeeds;

        appendWinningDirectionSection(sb, lga);

        boolean smallN = meanEffective < AnisotropyProbe.SMALL_N_EFFECTIVE_COLLISIONS_THRESHOLD;
        sb.append("# lgaSmallNEarlyTimeFlag=").append(smallN)
          .append(" (mean effective collisions/seed=").append(formatPrecise(meanEffective))
          .append(", mean total collisions/seed=").append(formatPrecise(meanTotal))
          .append(", threshold=").append(AnisotropyProbe.SMALL_N_EFFECTIVE_COLLISIONS_THRESHOLD)
          .append(smallN
                  ? " -- FEW real transfer events per seed, same small-N caveat as the committed Phase A report; the LGA anisotropy CI/p-value below is NOT independently power-verified at this campaign scale"
                  : " -- collision counts comfortably above the small-N threshold")
          .append('\n');
        appendCiVsPermutationReconciliation(sb, lga);
        appendPowerRecommendation(sb, meanEffective);
        sb.append("# structureFactorRidgeVerdict=");
        boolean anyNonzeroRidge = false;
        for (var e : lga.pooledSpectral().perDirection().entrySet()) {
            if (Math.abs(e.getValue().mean()) > 1e-12) {
                anyNonzeroRidge = true;
                break;
            }
        }
        sb.append(anyNonzeroRidge
                  ? "RIDGE PRESENT -- at least one pooled SPECTRAL direction has a nonzero ridge slope"
                  : "RIDGE ABSENT -- every pooled SPECTRAL direction's ridge slope is exactly zero, the same signature AnisotropyProbe's own Javadoc documents as EXPECTED for purely diffusive dynamics (no propagating branch, omega~i*D*k^2), not an instrument malfunction")
          .append('\n');
        sb.append('\n');
    }

    private static void appendDirectionRows(StringBuilder sb, long seed,
                                              String estimator,
                                              AnisotropyProbe.EstimatorResult result) {
        for (var d : StructureFactor.Direction.values()) {
            var dm = result.perDirection().get(d);
            sb.append("LGA_DIRECTION\t").append(seed).append('\t')
              .append(estimator).append('\t').append(d).append('\t')
              .append(formatPrecise(dm.magnitude())).append('\t')
              .append(dm.sampleSize()).append('\n');
        }
    }

    private static void appendNaiveSummaryRow(StringBuilder sb, String estimator,
                                                 AnisotropyProbe.BootstrapCi ci) {
        sb.append("LGA_SUMMARY\t").append(estimator).append('\t')
          .append(formatPrecise(ci.mean())).append('\t')
          .append(formatPrecise(ci.lower())).append('\t')
          .append(formatPrecise(ci.upper())).append('\t')
          .append(ci.nSeedsUsed()).append('\t').append(ci.nSeedsDegenerate())
          .append('\n');
    }

    /**
     * Critic Critical fix (item 1, relay batch): the winning-direction-
     * stability diagnostic. Computed for TRANSPORT only -- SPECTRAL's
     * per-seed magnitudes are EXACTLY 0.0 in every direction for every
     * seed (the RIDGE ABSENT finding), so an argmax over them would
     * report a spurious 100%-stable "winner" driven purely by iteration-
     * order tie-breaking, not a real signal; this is stated explicitly
     * rather than silently computed and left to mislead a reader.
     */
    private static void appendWinningDirectionSection(StringBuilder sb,
                                                         AnisotropyProbe.Report lga) {
        java.util.Map<StructureFactor.Direction, Integer> modeCounts = new java.util.EnumMap<>(StructureFactor.Direction.class);
        for (StructureFactor.Direction d : StructureFactor.Direction.values()) {
            modeCounts.put(d, 0);
        }
        for (SeedResult sr : lga.perSeed()) {
            StructureFactor.Direction winner = winningDirection(sr.transport());
            sb.append("WINNING_DIRECTION\tTRANSPORT\t").append(sr.seed())
              .append('\t').append(winner).append('\n');
            modeCounts.merge(winner, 1, Integer::sum);
        }
        StructureFactor.Direction modeDirection = null;
        int modeCount = -1;
        for (var e : modeCounts.entrySet()) {
            if (e.getValue() > modeCount) {
                modeCount = e.getValue();
                modeDirection = e.getKey();
            }
        }
        int totalSeeds = lga.perSeed().size();
        sb.append("WINNING_DIRECTION_STABILITY\tTRANSPORT\t").append(modeDirection)
          .append('\t').append(modeCount).append('\t').append(totalSeeds)
          .append('\t').append(formatPrecise((double) modeCount / totalSeeds))
          .append('\n');
        sb.append("# winningDirectionSpectralNote=omitted for SPECTRAL -- every per-seed SPECTRAL magnitude is EXACTLY 0.0 in all 3 directions (RIDGE ABSENT), so an argmax winner would be a spurious iteration-order artifact (100% \"stable\"), not a real signal\n");
    }

    private static StructureFactor.Direction winningDirection(AnisotropyProbe.EstimatorResult result) {
        StructureFactor.Direction best = null;
        double bestMagnitude = Double.NEGATIVE_INFINITY;
        for (var dm : result.perDirection().values()) {
            if (dm.magnitude() > bestMagnitude) {
                bestMagnitude = dm.magnitude();
                best = dm.direction();
            }
        }
        return best;
    }

    /**
     * Reviewer Important #2 / critic Significant #2 (item 2, relay
     * batch): an explicit in-artifact statement of which statistic is
     * authoritative, plus the cross-campaign pattern -- read live from
     * the committed Phase A artifact rather than hardcoded, so this
     * cannot silently drift out of sync with that file.
     */
    private static void appendCiVsPermutationReconciliation(StringBuilder sb,
                                                                AnisotropyProbe.Report lga) {
        double[] phaseA = readPhaseATransportPooledSummary();
        var pcTransport = lga.pooledTransport();
        sb.append("# ciVsPermutationReconciliation=the pooled TRANSPORT ratio's resample-then-aggregate bootstrap CI (bounded below by 1.0 by construction -- a max/min-of-3 order statistic, upward-biased by seed noise, AnisotropyProbe's own class Javadoc \"STACKED-REVIEW CORRECTION\") is NOT the significance statistic; permutationPValue IS, per that same locked convention. TWO-CAMPAIGN PATTERN: Phase A pooled TRANSPORT ratio=")
          .append(formatPrecise(phaseA[0])).append(" CI=[").append(formatPrecise(phaseA[1]))
          .append(',').append(formatPrecise(phaseA[2])).append("] (excludes 1.0) permutationP=")
          .append(formatPrecise(phaseA[3]))
          .append(" (NOT significant); Phase C pooled TRANSPORT ratio=")
          .append(formatPrecise(pcTransport.pooledRatio().orElse(Double.NaN)))
          .append(" CI=[").append(formatPrecise(pcTransport.pooledRatioCiLower()))
          .append(',').append(formatPrecise(pcTransport.pooledRatioCiUpper()))
          .append("] (excludes 1.0) permutationP=").append(formatPrecise(pcTransport.permutationPValue()))
          .append(" (NOT significant) -- in BOTH campaigns the CI-excludes-1.0 signal does NOT correspond to permutation significance, consistent with this CI being anti-conservative at this campaign scale, not evidence of a real per-seed direction effect. A reader must not treat \"CI excludes 1.0\" as significance evidence here.\n");
    }

    /**
     * Critic relay item 4: a concrete, non-decisional recommendation for
     * .23 -- framed as an input to the gate, not a posture choice.
     */
    private static void appendPowerRecommendation(StringBuilder sb,
                                                     double meanEffectivePerSeed) {
        sb.append("# powerRecommendationForGate=both Phase A (mean 29.25 effective collisions/seed) and Phase C (mean ")
          .append(formatPrecise(meanEffectivePerSeed))
          .append(", worse) fall below the informational small-N threshold (")
          .append(AnisotropyProbe.SMALL_N_EFFECTIVE_COLLISIONS_THRESHOLD)
          .append("). RECOMMENDATION (an input to the .23 gate, NOT a decision made here): a properly-powered follow-up anisotropy campaign, seeds 8->24 and ticks 128->400-500 -- this campaign's full wall time (~177s) leaves ample headroom for that scale-up. This is a DISTINCT decision point from the three isotropy postures in Section 5: the user may reasonably choose to defer any posture commitment until a properly-powered campaign exists, rather than choosing among postures on the current underpowered data.\n");
    }

    /**
     * Reads the committed Phase A artifact's own {@code POOLED_SUMMARY
     * TRANSPORT} row live (never hardcoded) -- {@code [pooledRatio,
     * ciLower, ciUpper, permutationPValue]}.
     */
    private static double[] readPhaseATransportPooledSummary() {
        try {
            List<String> lines = Files.readAllLines(Paths.get(PHASE_A_RELATIVE_PATH));
            for (String line : lines) {
                if (line.startsWith("POOLED_SUMMARY\tTRANSPORT\t")) {
                    String[] parts = line.split("\t");
                    return new double[] { Double.parseDouble(parts[2]),
                                           Double.parseDouble(parts[3]),
                                           Double.parseDouble(parts[4]),
                                           Double.parseDouble(parts[5]) };
                }
            }
            throw new IllegalStateException("Phase A artifact at " + PHASE_A_RELATIVE_PATH
                                             + " has no POOLED_SUMMARY TRANSPORT row");
        } catch (IOException e) {
            throw new IllegalStateException("failed to read Phase A artifact at "
                                             + PHASE_A_RELATIVE_PATH
                                             + " for the CI-vs-permutation cross-campaign reconciliation",
                                             e);
        }
    }

    private static void appendPooledDirectionRows(StringBuilder sb,
                                                     String estimator,
                                                     PooledResult pooled) {
        for (var d : StructureFactor.Direction.values()) {
            var stats = pooled.perDirection().get(d);
            sb.append("LGA_POOLED_DIRECTION\t").append(estimator).append('\t')
              .append(d).append('\t').append(formatPrecise(stats.mean()))
              .append('\t').append(formatPrecise(stats.ciLower())).append('\t')
              .append(formatPrecise(stats.ciUpper())).append('\t')
              .append(stats.nSeeds()).append('\n');
        }
    }

    private static void appendPooledSummaryRow(StringBuilder sb,
                                                  String estimator,
                                                  PooledResult pooled) {
        sb.append("LGA_POOLED_SUMMARY\t").append(estimator).append('\t')
          .append(formatPrecise(pooled.pooledRatio().orElse(Double.NaN)))
          .append('\t').append(formatPrecise(pooled.pooledRatioCiLower()))
          .append('\t').append(formatPrecise(pooled.pooledRatioCiUpper()))
          .append('\t').append(formatPrecise(pooled.permutationPValue()))
          .append('\t').append(formatPrecise(pooled.permutationNull95()))
          .append('\t').append(pooled.permutationCount()).append('\n');
    }

    private static void appendSpectralSection(StringBuilder sb, Report report) {
        sb.append("# === SECTION 2: SPECTRAL BROADENING, three-way K0 -> HYBRID -> LGA ===\n");
        sb.append("# columns(SPECTRAL_SUMMARY)=recordType\tsubstrate\textent\tfftLength\tstride\tnMembersChecked\tmeanPeakFraction\tminPeakFraction\tmaxPeakFraction\tmeanEntropy\tmeanAbsoluteLinewidthRadPerTick\tmeanFractionalLinewidth\tnMembersFractionalDefined\n");
        sb.append("# linewidthDefinition=half-power (FWHM) full-width around each member's dominant power-spectrum bin; ABSOLUTE = radians/tick (SpectralCadence#omegaRadPerTick); FRACTIONAL = absolute/|center frequency rad-per-tick| (dimensionless, undefined/excluded for a DC-peak member) -- user decision D3, BOTH reported, neither promoted to sole headline\n");
        appendSpectralRow(sb, report.k0Spectral());
        appendSpectralRow(sb, report.hybridSpectral());
        appendSpectralRow(sb, report.lgaSpectral());
        boolean hybridBroader = report.hybridSpectral().meanPeakFraction() < report.k0Spectral()
                                                                                     .meanPeakFraction();
        boolean lgaBroader = report.lgaSpectral().meanPeakFraction() < report.k0Spectral()
                                                                               .meanPeakFraction();
        double lgaDrop = report.k0Spectral().meanPeakFraction()
                          - report.lgaSpectral().meanPeakFraction();
        sb.append("# spectralProgressionVerdict=")
          .append(hybridBroader && lgaBroader ? "PROGRESSION PRESENT"
                                               : "PROGRESSION ABSENT OR PARTIAL")
          .append(" (K0 meanPeakFraction=").append(formatPrecise(report.k0Spectral()
                                                                          .meanPeakFraction()))
          .append(", HYBRID=").append(formatPrecise(report.hybridSpectral()
                                                             .meanPeakFraction()))
          .append(", LGA=").append(formatPrecise(report.lgaSpectral()
                                                          .meanPeakFraction()))
          .append(", LGA effect size (K0-LGA drop)=").append(formatPrecise(lgaDrop))
          .append(")\n");
        sb.append("# spectralMethodologyCaveat=K0 (extent 2,2,2) and HYBRID/LGA (extent 4,4,4) have DIFFERENT member counts/seeding draws -- this is an AGGREGATE distributional comparison (mean/min/max peakFraction, mean entropy across active members), not a per-member matched comparison; cadence (stride=225, phaseResolution=3600) IS matched, which is the load-bearing apples-to-apples property per the inviscid-ckn substantive-critique concern. K0's nMembersChecked excludes members seeded with exactly zero quanta (permanently non-rotating, a trivial DC case -- BaselineSpectrumHarness's own convention); HYBRID/LGA's nMembersChecked does NOT apply that filter (a member seeded at zero quanta can still acquire quanta via collision during the run, so a priori exclusion by initial quanta would not mean the same thing there) -- all active members are included for both collision-bearing substrates\n");
        sb.append('\n');
    }

    private static void appendSpectralRow(StringBuilder sb, SpectralSummary s) {
        sb.append("SPECTRAL_SUMMARY\t").append(s.substrate()).append('\t')
          .append(s.extent().x).append(',').append(s.extent().y).append(',')
          .append(s.extent().z).append('\t').append(s.fftLength()).append('\t')
          .append(s.stride()).append('\t').append(s.nMembersChecked())
          .append('\t').append(formatPrecise(s.meanPeakFraction())).append('\t')
          .append(formatPrecise(s.minPeakFraction())).append('\t')
          .append(formatPrecise(s.maxPeakFraction())).append('\t')
          .append(formatPrecise(s.meanEntropy())).append('\t')
          .append(formatPrecise(s.meanAbsoluteLinewidthRadPerTick())).append('\t')
          .append(formatPrecise(s.meanFractionalLinewidth())).append('\t')
          .append(s.nMembersFractionalDefined()).append('\n');
    }

    private static void appendCollisionSection(StringBuilder sb, Report report) {
        sb.append("# === SECTION 3: COLLISION-STATISTICS COMPARABILITY + GAP MECHANISM / BOOT-TRANSIENT ===\n");
        sb.append("# columns(COLLISION_WINDOW)=recordType\twindowStart\twindowEnd\thybridRatePerTick\thybridEffectiveFraction\tlgaRatePerTick\tlgaEffectiveFraction\n");
        List<Double> relativeGaps = new ArrayList<>();
        for (CollisionWindow w : report.windows()) {
            sb.append("COLLISION_WINDOW\t").append(w.windowStart()).append('\t')
              .append(w.windowEnd()).append('\t')
              .append(formatPrecise(w.hybridRatePerTick())).append('\t')
              .append(formatPrecise(w.hybridEffectiveFraction())).append('\t')
              .append(formatPrecise(w.lgaRatePerTick())).append('\t')
              .append(formatPrecise(w.lgaEffectiveFraction())).append('\n');
            relativeGaps.add(Math.abs(w.hybridRatePerTick() - w.lgaRatePerTick())
                              / w.hybridRatePerTick());
        }
        // Reviewer Important #1 / critic Significant #2 (item 3, relay
        // batch): report the FULL series, not first-vs-last -- an
        // endpoint-only comparison hides interior extrema (this series
        // peaks at window 1, not the last window).
        int gapPeakWindow = argmax(relativeGaps);
        sb.append("# gapMechanismVerdict=relativeRateGapByWindow=")
          .append(formatSeries(relativeGaps))
          .append(", peak at window ").append(gapPeakWindow).append(" (")
          .append(formatPrecise(relativeGaps.get(gapPeakWindow)))
          .append(")")
          .append(gapPeakWindow == 0 || gapPeakWindow == relativeGaps.size() - 1
                  ? " -- monotonic-consistent shape"
                  : " -- NON-MONOTONIC: spikes at an INTERIOR window then partially recovers, more consistent with a genuine boot transient than a stable, persistent step-granularity/quantisation-floor effect (which would predict a roughly flat or monotonic series) -- reported as data, not further mechanistically explained")
          .append('\n');
        sb.append("# columns(ANISOTROPY_WITHIN_WINDOW)=recordType\tquartileIndex\ttickStart\ttickEndExclusive\ttotalCollisionsAcrossSeeds (TOTAL, not effective-only -- see class Javadoc)\n");
        for (TickQuartileHomogeneity q : report.anisotropyWithinWindowHomogeneity()) {
            sb.append("ANISOTROPY_WITHIN_WINDOW\t").append(q.quartileIndex())
              .append('\t').append(q.tickStart()).append('\t')
              .append(q.tickEndExclusive()).append('\t')
              .append(q.totalCollisionsAcrossSeeds()).append('\n');
        }
        sb.append("# anisotropyWithinWindowHomogeneityVerdict=")
          .append(anisotropyHomogeneityVerdict(report.anisotropyWithinWindowHomogeneity()))
          .append('\n');
        sb.append("# columns(COLLISION_FIELD)=recordType\tfield\thybridValue\tlgaValue\tverifiedComparable\tnote\n");
        for (FieldComparison fc : report.fieldComparisons()) {
            sb.append("COLLISION_FIELD\t").append(fc.field()).append('\t')
              .append(fc.hybridValue()).append('\t').append(fc.lgaValue())
              .append('\t').append(fc.verifiedComparable()).append('\t')
              .append(fc.note()).append('\n');
        }
        sb.append('\n');
    }

    /**
     * Critic Significant #3 / relay item 5: reports the anisotropy
     * campaign's own within-window quartile homogeneity as DATA, framed
     * as what it can and cannot show (TOTAL collisions, a coarse proxy --
     * see {@link TickQuartileHomogeneity}'s Javadoc) rather than a firm
     * mechanistic conclusion.
     */
    private static String anisotropyHomogeneityVerdict(List<TickQuartileHomogeneity> quartiles) {
        List<Long> counts = new ArrayList<>();
        for (TickQuartileHomogeneity q : quartiles) {
            counts.add(q.totalCollisionsAcrossSeeds());
        }
        long first = counts.get(0);
        long last = counts.get(counts.size() - 1);
        String shape = "quartileCounts=" + counts;
        if (first == 0) {
            return shape + " -- first quartile has zero recorded collisions, ratio-based trend undefined; reported as raw counts only";
        }
        double relChange = (double) (last - first) / first;
        return shape + ", first->last relative change=" + formatPrecise(relChange)
               + (Math.abs(relChange) < 0.1
                  ? " -- roughly HOMOGENEOUS across the 128-tick anisotropy window (TOTAL-collision proxy; does not itself verify direction-resolved homogeneity, only the pooled rate) -- no strong evidence the OLS transport-rate fit's time-homogeneity assumption is violated at this coarse resolution"
                  : " -- NOT homogeneous (TOTAL-collision proxy); the OLS transport-rate fit's implicit time-homogeneity assumption is NOT verified for this campaign -- carried forward as an open question, not resolved here (critic Significant #3)");
    }

    private static int argmax(List<Double> values) {
        int best = 0;
        for (int i = 1; i < values.size(); i++) {
            if (values.get(i) > values.get(best)) {
                best = i;
            }
        }
        return best;
    }

    private static String formatSeries(List<Double> values) {
        StringBuilder s = new StringBuilder("[");
        for (int i = 0; i < values.size(); i++) {
            if (i > 0) {
                s.append(',');
            }
            s.append(formatPrecise(values.get(i)));
        }
        return s.append(']').toString();
    }

    private static void appendEquilibriumSection(StringBuilder sb, Report report) {
        sb.append("# === SECTION 4: EQUILIBRIUM QUANTA-DISTRIBUTION CHARACTERIZATION (SDB caveat, user decision 6) ===\n");
        sb.append("# equilibriumStatisticsAssumeSDB=false (explicit -- CollisionTable's collision rule is an ACCEPTED non-bijection, see class Javadoc)\n");
        appendHistogramRow(sb, report.quantaBefore());
        appendHistogramRow(sb, report.quantaAfter());
        sb.append("# sdbClosedForm=window=").append(report.sdb().window())
          .append(" statesChecked=").append(report.sdb().statesChecked())
          .append(" preimageCountByOutputDiff=").append(report.sdb().preimageCountByOutputDiff())
          .append(" balanced=").append(report.sdb().balanced())
          .append(" (formal parity-floor fact: 3 preimages at diff=0, 2 at |diff|=1, 1 at |diff|>=2 -- odd-difference pairs can never reach exact equality, only oscillate at |diff|=1; even-difference pairs converge to and then permanently self-loop at diff=0)\n");
        // Reviewer Important #1 / critic Significant #2 (item 3, relay
        // batch): full 4-window series, not first-vs-last -- and the
        // ONE-WINDOW LAG against the gapMechanismVerdict series (which
        // peaks at window 1) stated explicitly rather than implying one
        // mechanism.
        List<Double> noOpFractions = new ArrayList<>();
        for (CollisionWindow w : report.windows()) {
            noOpFractions.add(1.0 - w.lgaEffectiveFraction());
        }
        int noOpPeakWindow = argmax(noOpFractions);
        List<Double> relativeGaps = new ArrayList<>();
        for (CollisionWindow w : report.windows()) {
            relativeGaps.add(Math.abs(w.hybridRatePerTick() - w.lgaRatePerTick())
                              / w.hybridRatePerTick());
        }
        int gapPeakWindow = argmax(relativeGaps);
        sb.append("# equalityConcentrationEmpiricalTrend=lgaNoOpFractionByWindow=")
          .append(formatSeries(noOpFractions))
          .append(", peak at window ").append(noOpPeakWindow).append(" (")
          .append(formatPrecise(noOpFractions.get(noOpPeakWindow)))
          .append(")")
          .append(" -- OVERALL RISING (first->last: ").append(formatPrecise(noOpFractions.get(0)))
          .append(" -> ").append(formatPrecise(noOpFractions.get(noOpFractions.size() - 1)))
          .append("), dynamical corroboration of the SDB closed-form's diff=0 preimage concentration (pairs settling into the permanent equality self-loop), but peaks at an INTERIOR window (")
          .append(noOpPeakWindow)
          .append(") then dips slightly rather than rising monotonically. NOTE: Section 3's relativeRateGapByWindow peaks at window ")
          .append(gapPeakWindow).append(", i.e. ")
          .append(noOpPeakWindow == gapPeakWindow + 1
                  ? "exactly ONE WINDOW LATER than this trend's own peak"
                  : "at a different window than this trend's own peak (window " + noOpPeakWindow + ")")
          .append(" -- the two trends are NOT the same mechanism showing up identically; stated explicitly rather than implying a single coupled cause.\n");
        sb.append('\n');
    }

    private static void appendHistogramRow(StringBuilder sb,
                                              QuantaHistogramSummary h) {
        sb.append("QUANTA_HISTOGRAM\t").append(h.label()).append('\t')
          .append(h.n()).append('\t').append(formatPrecise(h.mean())).append('\t')
          .append(formatPrecise(h.variance())).append('\t').append(h.min())
          .append('\t').append(h.max()).append('\n');
    }

    private static void appendPostureSection(StringBuilder sb, Report report) {
        sb.append("# === SECTION 5: ISOTROPY POSTURES -- DATA ONLY, NO POSTURE SELECTED (escalated to user/orchestrator) ===\n");
        var transport = report.lgaAnisotropy().pooledTransport();
        var spectral = report.lgaAnisotropy().pooledSpectral();
        sb.append("# posture(i)_acceptAndCharacterize_evidence=LGA pooled TRANSPORT ratio=")
          .append(formatPrecise(transport.pooledRatio().orElse(Double.NaN)))
          .append(" CI=[").append(formatPrecise(transport.pooledRatioCiLower()))
          .append(',').append(formatPrecise(transport.pooledRatioCiUpper()))
          .append("] permutationP=").append(formatPrecise(transport.permutationPValue()))
          .append("; pooled SPECTRAL ratio=")
          .append(formatPrecise(spectral.pooledRatio().orElse(Double.NaN)))
          .append(" -- this posture accepts whatever these numbers say at face value, characterizing rather than correcting for anisotropy\n");
        sb.append("# posture(ii)_fchcProjection_evidence=NOT constructed in this bead -- the current lattice remains the 12-neighbor FCC arrangement, not a 4D FCHC projection (design memo's own named risk); the measured ratio/CI above is what an FCHC-style projection would need to demonstrably reduce if adopted, but no FCHC variant has been built or measured to compare against\n");
        sb.append("# posture(iii)_orientationalStateRestoresIsotropy_evidence=the pooled TRANSPORT permutation p-value (")
          .append(formatPrecise(transport.permutationPValue()))
          .append(") is ")
          .append(transport.permutationPValue() > 0.05 ? "NOT significant at p<0.05"
                                                         : "significant at p<0.05")
          .append(", i.e. the measured anisotropy is ")
          .append(transport.permutationPValue() > 0.05
                  ? "statistically indistinguishable from an isotropic null at this campaign's power -- CONSISTENT WITH but NOT PROOF OF this hypothesis (small-N caveat above applies: absence of significance here may reflect insufficient statistical power, not genuine isotropy)"
                  : "statistically distinguishable from an isotropic null -- evidence AGAINST this hypothesis as measured, though the small-N caveat above still bounds confidence in this finding")
          .append('\n');
        sb.append("# ESCALATION=this report presents measured anisotropy ratios+CIs from BOTH estimators and evidence for all three postures; NO posture is selected by this bead -- selecting one is out of scope per the locked design and would be exactly the silent-scope-reduction failure the review gates exist to catch. Escalate to the user/orchestrator with these numbers for the posture decision.\n");
    }

    private static String formatPrecise(double v) {
        return String.format(Locale.ROOT, "%.9e", v);
    }

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
