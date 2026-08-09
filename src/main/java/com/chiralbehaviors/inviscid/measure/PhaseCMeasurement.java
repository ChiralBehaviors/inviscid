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
    // Post-.23-gate fix (Important #4): the single named source of truth
    // for "both substrates share phaseResolution=3600 under cadence 2A" --
    // replaces the bare literal 3600 previously scattered across seeding
    // helpers below (a future atlas-header change to the LGA side would
    // otherwise desynchronise silently; newHybrid/newLga now assert
    // against this constant instead).
    public static final int     PHASE_RESOLUTION      = Necronomata.PHASE_RESOLUTION;

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

    /**
     * Quanta-value population summary (distribution-shape data point).
     * {@code bins} (fix-round item 5, S5) is the ACTUAL per-value
     * histogram -- exact population count for every integer value in
     * {@code [min,max]} inclusive, zero-filled for unobserved values --
     * so this record finally carries the histogram its name claims; the
     * scalar {@code mean}/{@code variance} fields remain the distributional
     * SUMMARY, not a replacement for the bins; {@code
     * PhaseCMeasurementTest} reconciles the bins against this summary row
     * (n, mean, variance, and cross-run total-quanta conservation).
     *
     * <p>Round-3 correction (critic S-NEW-3): an earlier version of this
     * Javadoc claimed {@link QuantaExchangeRule}'s pairwise single-unit
     * transfer composes into a discrete maximum principle over the WHOLE
     * system -- that claim is FALSE and has been removed. {@link
     * LatticeGasAutomaton#tick} accumulates every contact's delta into a
     * per-member buffer across an entire tick and applies it once at the
     * end, with each contact's lookup reading the STALE pre-tick quanta --
     * so a member touched by several contacts in one tick can net-move by
     * more than the pairwise +/-1 step (measured: up to 4 contacts on one
     * member in a single tick, net change up to 3). The REAL bound on this
     * path is {@code LatticeGasAutomaton}'s own {@code
     * checkExactnessCeiling} guard ({@code LatticeGasAutomaton.java:
     * 481-488}), which throws rather than let a member's quanta
     * random-walk past {@code CollisionSweep.QUANTA_EXACTNESS_SAFETY_MARGIN}
     * -- that guard, not a maximum principle, is why {@code histogramOf}'s
     * {@code [min,max]} zero-fill loop stays small in practice.
     */
    public record QuantaHistogramSummary(String label, long n, double mean,
                                          double variance, long min, long max,
                                          java.util.SortedMap<Long, Long> bins) {
        public QuantaHistogramSummary {
            bins = java.util.Collections.unmodifiableSortedMap(new java.util.TreeMap<>(bins));
        }
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

    /**
     * Fix-round Critical 1: one anisotropy-campaign seed's excitation/
     * no-op census, characterizing whether the SDB closed form's diff=0
     * preimage concentration is merely an equilibrium TREND (as Section
     * 4's long-run characterization shows) or, on THIS campaign, the
     * absorbing INITIAL CONDITION -- every background member is seeded
     * quanta=0, and {@code phase[i] = floorMod(phase[i]+quanta[i], P)}
     * (bead's own accepted rule) means a zero-quanta member never rotates,
     * so a {@code (0,0)} contact is a PERMANENT no-op, not a transient
     * one. {@code excitedMembersInitial}/{@code cellsTouchedInitial} are
     * measured (not assumed) from the SAME seeded state {@link
     * #lgaSubstrateFactory}'s additive census sink captures at
     * construction, before any tick runs; {@code *Final} are measured
     * from the SAME {@link LatticeGasAutomaton} instance after
     * {@link AnisotropyProbe} has driven every tick (a live reference,
     * not a re-simulation). {@code collisionsAtFirstTick}/{@code
     * collisionsAtLastTick} read {@link CollisionStatistics#collisionsPerTick()}
     * at tick indices 0 and {@code ticks-2} (the last {@code tick()} call
     * argument -- see {@code computeAnisotropyWithinWindowHomogeneity}'s
     * {@code actualTickCount} convention) -- a flat pair is the dynamical
     * signature of a static contact graph.
     */
    public record AnisotropyCampaignEquality(long seed, long totalCollisions,
                                              long effectiveCollisions,
                                              long excitedMembersInitial,
                                              int cellsTouchedInitial,
                                              long excitedMembersFinal,
                                              int cellsTouchedFinal,
                                              long collisionsAtFirstTick,
                                              long collisionsAtLastTick) {
        public double noOpFraction() {
            return totalCollisions == 0L ? Double.NaN
                                          : (double) (totalCollisions
                                                       - effectiveCollisions)
                                            / totalCollisions;
        }
    }

    /**
     * Fix-round Critical 2: one scope's (a seed, or a pooled/long-run
     * substrate) raw per-CONTACT-direction collision census -- {@link
     * FccNeighborhood}'s positive directions {@code 1..6} (the six
     * canonical FCC {@code <110>}-type offsets; {@link
     * CollisionStatistics}'s own Javadoc: negative directions are never
     * populated by a real caller). DISTINCT from {@code LGA_DIRECTION}'s
     * {@link StructureFactor.Direction} ({@code X100}/{@code X110}/
     * {@code X111}) k-space probe directions -- two unrelated direction
     * spaces sharing the word "direction", never to be confused.
     */
    public record ContactDirectionCensus(String scope,
                                           java.util.Map<Integer, Long> perDirection) {
        public ContactDirectionCensus {
            perDirection = java.util.Map.copyOf(perDirection);
        }

        public long total() {
            long sum = 0L;
            for (long v : perDirection.values()) {
                sum += v;
            }
            return sum;
        }
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
                          ContactAtlas.Header atlasHeader,
                          List<AnisotropyCampaignEquality> anisotropyCampaignEquality,
                          List<ContactDirectionCensus> anisotropyContactDirectionPerSeed,
                          List<ContactDirectionCensus> anisotropyEffectiveContactDirectionPerSeed,
                          ContactDirectionCensus longRunHybridContactDirection,
                          ContactDirectionCensus longRunLgaContactDirection) {
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
     * @param perSeedAuditSink OPTIONAL (nullable) side-channel, same
     *                         additive-instrumentation pattern as {@code
     *                         perSeedStatsSink}: if non-null, every {@link
     *                         ConservationAudit} this factory constructs
     *                         is ALSO stashed here keyed by seed, purely
     *                         so a caller can sum {@link
     *                         ConservationAudit#ledger()} sizes AFTER the
     *                         campaign for an exact {@code
     *                         conservationTicksAudited} total (C2 fix:
     *                         read from the audit's own ledger, never
     *                         hand-derived from an assumption about this
     *                         factory's or {@code AnisotropyProbe}'s loop
     *                         shape).
     * @param perSeedInitialCensusSink OPTIONAL (nullable), fix-round
     *                         Critical 1: if non-null, every seed's
     *                         excited-member/cells-touched census (see
     *                         {@link #activeMemberCensus}) is captured
     *                         HERE -- synchronously, immediately after
     *                         {@code seedPacket}, before the {@link
     *                         SubstrateFactory.Substrate} is returned and
     *                         before any tick runs -- so a caller can read
     *                         the campaign's OWN t=0 state without
     *                         re-simulating it. Purely observational: a
     *                         {@code forEachCell}/{@code quantaAt} read,
     *                         no RNG draw, no mutation.
     * @param perSeedLgaSink   OPTIONAL (nullable), same pattern: if
     *                         non-null, the constructed {@link
     *                         LatticeGasAutomaton} instance itself is
     *                         stashed here keyed by seed -- {@link
     *                         AnisotropyProbe} drives this SAME instance's
     *                         ticks externally, so a caller reading it
     *                         back AFTER the campaign completes sees the
     *                         campaign's own final state, not a
     *                         re-simulation.
     */
    private static SubstrateFactory lgaSubstrateFactory(ContactAtlas atlas,
                                                          CollisionTable collisions,
                                                          java.util.Map<Long, CollisionStatistics> perSeedStatsSink,
                                                          java.util.Map<Long, ConservationAudit> perSeedAuditSink,
                                                          java.util.Map<Long, ActiveCensus> perSeedInitialCensusSink,
                                                          java.util.Map<Long, LatticeGasAutomaton> perSeedLgaSink) {
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
            if (perSeedInitialCensusSink != null) {
                perSeedInitialCensusSink.put(seed, activeMemberCensus(lga));
            }
            if (perSeedLgaSink != null) {
                perSeedLgaSink.put(seed, lga);
            }
            ConservationAudit audit = new ConservationAudit(lga, true);
            if (perSeedAuditSink != null) {
                perSeedAuditSink.put(seed, audit);
            }
            AuditedRun run = new AuditedRun(lga, audit);
            return new SubstrateFactory.Substrate(lga, run, lga.statistics());
        };
    }

    /**
     * Fix-round Critical 1: one snapshot's excited-member/cells-touched
     * census over EVERY active (even-parity-cell) member slot -- purely
     * observational ({@link QuantaField#forEachCell}/{@link
     * QuantaField#quantaAt}, no mutation, no RNG draw), reused for both
     * the t=0 and post-campaign final census (see {@link
     * #lgaSubstrateFactory}'s Javadoc).
     */
    record ActiveCensus(long excitedMembers, int cellsTouched) {
    }

    static ActiveCensus activeMemberCensus(QuantaField field) {
        long[] excited = { 0L };
        int[] cells = { 0 };
        field.forEachCell(cell -> {
            int base = field.indexOfCell(cell);
            boolean touched = false;
            for (int local = 0; local < 30; local++) {
                if (field.quantaAt(base + local) != 0L) {
                    excited[0]++;
                    touched = true;
                }
            }
            if (touched) {
                cells[0]++;
            }
        });
        return new ActiveCensus(excited[0], cells[0]);
    }

    /**
     * The even-parity sublattice's cell count for {@code extent} --
     * {@code x*y*z/2}, exact whenever every axis is even ({@link
     * FccNeighborhood}'s own precondition): for each fixed {@code (i,j)},
     * exactly half of the {@code z} values give an even {@code i+j+k}
     * (an even {@code z} alternates parity {@code z/2} times each way).
     * Computed, never hardcoded, so a future extent change follows
     * automatically.
     */
    static long evenParityCellCount(Point3i extent) {
        return ((long) extent.x * extent.y * extent.z) / 2L;
    }

    private static void seedRandomPhases(LatticeGasAutomaton lga, Point3i extent,
                                          long seed) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        // Reads the ALREADY-CONSTRUCTED lga's own phaseResolution() rather
        // than a bare literal (Important #4) -- if a future atlas header
        // ever desynchronises the LGA's computed phaseResolution from the
        // shared cadence-2A value, this seeding draw range follows it
        // automatically instead of silently going stale.
        int phaseResolution = lga.phaseResolution();
        int[] phases = new int[length];
        for (int i = 0; i < length; i++) {
            phases[i] = random.nextInt(phaseResolution);
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
            // PHASE_RESOLUTION (Important #4), not a bare 3600 -- this
            // condition is shared by BOTH substrates (newHybrid/newLga
            // below assert their own phaseResolution() matches it).
            phase[i] = random.nextInt(PHASE_RESOLUTION);
            quanta[i] = random.nextInt(range) + minQuanta;
        }
        return new SharedInitialCondition(phase, quanta);
    }

    private static HybridAutomaton newHybrid(ContactAtlas atlas,
                                              SharedInitialCondition ic,
                                              CollisionStatistics statistics,
                                              Point3i extent) {
        Necronomata automaton = new Necronomata(extent);
        if (automaton.phaseResolution() != PHASE_RESOLUTION) {
            throw new IllegalStateException("hybrid substrate phaseResolution="
                                             + automaton.phaseResolution()
                                             + " does not match the shared cadence-2A PHASE_RESOLUTION="
                                             + PHASE_RESOLUTION
                                             + " the shared initial condition was seeded against (Important #4)");
        }
        int length = ic.phase().length;
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < length; i++) {
                angle[i] = (float) (2.0 * Math.PI * ic.phase()[i] / PHASE_RESOLUTION);
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
        if (lga.phaseResolution() != PHASE_RESOLUTION) {
            throw new IllegalStateException("LGA substrate phaseResolution="
                                             + lga.phaseResolution()
                                             + " does not match the shared cadence-2A PHASE_RESOLUTION="
                                             + PHASE_RESOLUTION
                                             + " the shared initial condition was seeded against (Important #4)");
        }
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

    /**
     * {@link #summarizeK0()}'s result plus the exact number of ticks its
     * OWN strict {@link ConservationAudit} audited -- read from {@link
     * ConservationAudit#ledger()}'s size (post-C2-fix: never hand-derived
     * from an assumption about this method's own loop shape, let alone
     * another driver's).
     */
    private record K0Result(SpectralSummary summary, long conservationTicksAudited) {
    }

    private static K0Result summarizeK0() {
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
        // C2 fix: read the audit's OWN ledger size, never hand-derive it
        // from an assumption about this loop's shape.
        long k0TicksAudited = audit.ledger().size();

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

        SpectralSummary summary = new SpectralSummary("K0", extent, BaselineSpectrumHarness.FFT_LENGTH,
                                    BaselineSpectrumHarness.STRIDE, agg.n(),
                                    agg.sumFraction() / agg.n(), agg.minFraction(),
                                    agg.maxFraction(), agg.sumEntropy() / agg.n(),
                                    agg.sumAbsoluteLinewidth() / agg.n(),
                                    agg.nFractionalDefined() == 0 ? Double.NaN
                                                                   : agg.sumFractionalLinewidth()
                                                                     / agg.nFractionalDefined(),
                                    agg.nFractionalDefined());
        return new K0Result(summary, k0TicksAudited);
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

    /**
     * Fix-round S5: package-private (was {@code private}) for direct unit
     * testing, exactly like this file's other narrative/data pure helpers.
     * {@code bins} is the EXACT per-value population count over
     * {@code [min,max]} inclusive -- see {@link QuantaHistogramSummary}'s
     * Javadoc for why this range is always small (a discrete maximum
     * principle, not an assumption).
     */
    static QuantaHistogramSummary histogramOf(String label, long[] quanta) {
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
        java.util.SortedMap<Long, Long> bins = new java.util.TreeMap<>();
        for (long v = min; v <= max; v++) {
            bins.put(v, 0L);
        }
        for (long q : quanta) {
            bins.merge(q, 1L, Long::sum);
        }
        return new QuantaHistogramSummary(label, n, mean, variance, min, max,
                                           bins);
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
        java.util.Map<Long, ConservationAudit> perSeedAuditSink = new java.util.LinkedHashMap<>();
        // Fix-round Critical 1/2 additive sinks -- see lgaSubstrateFactory's
        // Javadoc for exactly what each captures and when.
        java.util.Map<Long, ActiveCensus> perSeedInitialCensusSink = new java.util.LinkedHashMap<>();
        java.util.Map<Long, LatticeGasAutomaton> perSeedLgaSink = new java.util.LinkedHashMap<>();
        AnisotropyProbe.Report lgaAnisotropy = AnisotropyProbe.runCampaign(AnisotropyProbe.DEFAULT_EXTENT,
                                                                             AnisotropyProbe.DEFAULT_SEEDS,
                                                                             AnisotropyProbe.DEFAULT_TICKS,
                                                                             AnisotropyProbe.DEFAULT_PACKET_QUANTA,
                                                                             lgaSubstrateFactory(atlas,
                                                                                                  collisions,
                                                                                                  perSeedStatsSink,
                                                                                                  perSeedAuditSink,
                                                                                                  perSeedInitialCensusSink,
                                                                                                  perSeedLgaSink));
        // C2 fix: sum every per-seed audit's OWN ledger size -- never a
        // hand-derived assumption about AnisotropyProbe.runOneSeed's loop
        // shape (the STRICT header claim ("every driver, every
        // sub-measurement") is only true if this reads the audits
        // themselves).
        long anisotropyTicksAudited = 0L;
        for (ConservationAudit audit : perSeedAuditSink.values()) {
            anisotropyTicksAudited += audit.ledger().size();
        }
        conservationTicks += anisotropyTicksAudited;
        List<TickQuartileHomogeneity> anisotropyHomogeneity = computeAnisotropyWithinWindowHomogeneity(perSeedStatsSink,
                                                                                                          AnisotropyProbe.DEFAULT_TICKS);

        // Fix-round Critical 1/2: the anisotropy campaign's own
        // excitation/no-op census and raw per-CONTACT-direction census,
        // built ENTIRELY from the additive sinks above -- the campaign
        // itself (AnisotropyProbe) is untouched.
        int anisotropyActualTickCount = AnisotropyProbe.DEFAULT_TICKS - 1;
        List<AnisotropyCampaignEquality> anisotropyCampaignEquality = new ArrayList<>();
        List<ContactDirectionCensus> anisotropyContactDirectionPerSeed = new ArrayList<>();
        List<ContactDirectionCensus> anisotropyEffectiveContactDirectionPerSeed = new ArrayList<>();
        for (SeedResult sr : lgaAnisotropy.perSeed()) {
            long seed = sr.seed();
            ActiveCensus initial = perSeedInitialCensusSink.get(seed);
            ActiveCensus finalCensus = activeMemberCensus(perSeedLgaSink.get(seed));
            CollisionStatistics stats = perSeedStatsSink.get(seed);
            long firstTick = stats.collisionsPerTick().getOrDefault(0, 0L);
            long lastTick = stats.collisionsPerTick()
                                  .getOrDefault(anisotropyActualTickCount - 1,
                                                0L);
            anisotropyCampaignEquality.add(new AnisotropyCampaignEquality(seed,
                                                                             sr.totalCollisions(),
                                                                             sr.effectiveCollisions(),
                                                                             initial.excitedMembers(),
                                                                             initial.cellsTouched(),
                                                                             finalCensus.excitedMembers(),
                                                                             finalCensus.cellsTouched(),
                                                                             firstTick,
                                                                             lastTick));
            anisotropyContactDirectionPerSeed.add(new ContactDirectionCensus(Long.toString(seed),
                                                                                positiveDirectionCounts(stats)));
            anisotropyEffectiveContactDirectionPerSeed.add(new ContactDirectionCensus(Long.toString(seed),
                                                                                          positiveEffectiveDirectionCounts(stats)));
        }

        // --- Sub-measurement 2/4: spectral broadening, three-way ---
        K0Result k0Result = summarizeK0();
        SpectralSummary k0 = k0Result.summary();
        conservationTicks += k0Result.conservationTicksAudited();

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
        // C2 fix: the audit's OWN ledger size, not the tick counter (which
        // happens to agree today but is not itself an audit read).
        conservationTicks += hybridSpectralRun.audit().ledger().size();

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
        // C2 fix: the audit's OWN ledger size (see hybrid counterpart above).
        conservationTicks += lgaSpectralRun.audit().ledger().size();

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
        // C2 fix: sum both long-run audits' OWN ledger sizes -- never a
        // hand-derived "2 substrates x LONG_RUN_TICKS" assumption.
        conservationTicks += hybridLongRun.audit().ledger().size()
                              + lgaLongRun.audit().ledger().size();

        long[] quantaAfter = snapshotQuanta(lgaLong);

        List<FieldComparison> fieldComparisons = new ArrayList<>();
        fieldComparisons.add(new FieldComparison("totalCollisions",
                                                   Long.toString(hybridLongStats.totalCollisions()),
                                                   Long.toString(lgaLongStats.totalCollisions()),
                                                   true,
                                                   "verified comparable: HybridVsLgaConsistencyTest#aggregateStatisticsAgreeBeyondDivergence "
                                                   + "(test 8) asserts the total-collision RATE gap directly (fix-round S1 correction -- "
                                                   + "the prior note cited test 8 for a claim it did not actually assert; test 8 asserted "
                                                   + "only totalCollisions()>0 and the effective-ratio gap), tolerance 0.15, same TOLERANCE "
                                                   + "constant as the effective-ratio assertion below"));
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

        // Fix-round Critical 2: the long run's own per-CONTACT-direction
        // census for both substrates, as first-class DATA rows (previously
        // only a stringified, "NOT independently verified" note on
        // collisionsPerDirection above -- the raw counts are unchanged,
        // only how they are reported).
        ContactDirectionCensus longRunHybridContactDirection = new ContactDirectionCensus("HYBRID",
                                                                                             positiveDirectionCounts(hybridLongStats));
        ContactDirectionCensus longRunLgaContactDirection = new ContactDirectionCensus("LGA",
                                                                                          positiveDirectionCounts(lgaLongStats));

        // conservationViolations stays exactly 0: every driver above is
        // wrapped in a STRICT ConservationAudit, which throws immediately
        // on any violation -- reaching this line is itself the proof.
        return new Report(lgaAnisotropy, k0, hybridSummary, lgaSummary, windows,
                           anisotropyHomogeneity, fieldComparisons,
                           histogramOf("before", quantaBefore),
                           histogramOf("after", quantaAfter), sdb,
                           conservationTicks, conservationViolations,
                           atlas.header(), anisotropyCampaignEquality,
                           anisotropyContactDirectionPerSeed,
                           anisotropyEffectiveContactDirectionPerSeed,
                           longRunHybridContactDirection,
                           longRunLgaContactDirection);
    }

    /**
     * The six canonical FCC contact directions {@code 1..6} only -- {@link
     * CollisionStatistics}'s own class Javadoc: a real caller ({@code
     * CollisionSweep}) only ever supplies {@code Contact.direction()},
     * which never carries a negative direction, so the negative half of
     * {@link CollisionStatistics#collisionsPerDirection()}'s 12-entry map
     * is always exactly zero -- filtered out here rather than carried as
     * dead always-zero columns through every census row.
     */
    static java.util.Map<Integer, Long> positiveDirectionCounts(CollisionStatistics stats) {
        java.util.Map<Integer, Long> positive = new java.util.TreeMap<>();
        for (var e : stats.collisionsPerDirection().entrySet()) {
            if (e.getKey() > 0) {
                positive.put(e.getKey(), e.getValue());
            }
        }
        return positive;
    }

    /**
     * Fix-round item 2c (round 3): the same positive-direction filter as
     * {@link #positiveDirectionCounts}, but over {@link
     * CollisionStatistics#effectiveCollisionsPerDirection()} -- the
     * EFFECTIVE-transfer-only counts (the population that actually moves
     * quanta and drives the hopping-tensor/{@code D_hat} estimate) rather
     * than the raw recorded-contact counts.
     */
    static java.util.Map<Integer, Long> positiveEffectiveDirectionCounts(CollisionStatistics stats) {
        java.util.Map<Integer, Long> positive = new java.util.TreeMap<>();
        for (var e : stats.effectiveCollisionsPerDirection().entrySet()) {
            if (e.getKey() > 0) {
                positive.put(e.getKey(), e.getValue());
            }
        }
        return positive;
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
        sb.append("# phaseResolution=").append(PHASE_RESOLUTION)
          .append(" (both substrates, cadence 2A -- see LatticeGasAutomaton class Javadoc)\n");
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
        // Reads SpectralCadence#nyquistQuantaBound rather than
        // re-deriving the same arithmetic inline a second time.
        int nyquistQuantaBound = new SpectralCadence(PHASE_RESOLUTION,
                                                       SPECTRAL_STRIDE).nyquistQuantaBound();
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
        appendContactDirectionCensusSection(sb, report);
        appendSpectralSection(sb, report);
        appendCollisionSection(sb, report);
        appendEquilibriumSection(sb, report);
        appendAnisotropyCampaignEqualitySection(sb, report);
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
        boolean anyNaNRidge = false;
        for (var e : lga.pooledSpectral().perDirection().entrySet()) {
            double mean = e.getValue().mean();
            if (Double.isNaN(mean)) {
                anyNaNRidge = true;
                continue;
            }
            if (Math.abs(mean) > 1e-12) {
                anyNonzeroRidge = true;
            }
        }
        sb.append(ridgeVerdictText(anyNaNRidge, anyNonzeroRidge, smallN)).append('\n');
        sb.append('\n');
    }

    /**
     * Important #7 pure helper (NaN pooled-SPECTRAL fallthrough fix): a
     * NaN pooled-SPECTRAL mean is an INSTRUMENT ANOMALY, not silently
     * treated as "every direction's ridge slope is exactly zero" (the
     * pre-fix bug -- {@code Math.abs(NaN) > 1e-12} is {@code false}, so
     * the old unconditional two-way branch fell through to RIDGE ABSENT
     * even when the data was undefined). Package-private for direct unit
     * testing of all three branches with synthetic booleans.
     *
     * <p>Fix-round S3: the RIDGE ABSENT branch no longer asserts "purely
     * diffusive dynamics" as the sole explanation for an all-exact-0.0
     * result -- {@code AnisotropyProbe}'s own non-fabrication contract
     * documents that the SAME all-zero signature is produced by the K=0
     * (collision-free, no-signal) baseline BY CONSTRUCTION, so "ridge
     * absent" (measured) and "purely diffusive" (a dynamics CLAIM the
     * measurement alone cannot support) are kept distinct. {@code
     * transportSmallNUsedAsSpectralProxy} branches the wording: when
     * true, both explanations are stated side by side and no dynamics
     * claim is made; when false, the diffusive-signature reading is
     * offered as the more likely account without foreclosing the
     * alternative.
     *
     * <p>Naming honesty (critic MINOR-4, round 3): despite this verdict
     * being about the SPECTRAL ridge, the flag passed in is {@code
     * lgaSmallNEarlyTimeFlag} -- a TRANSPORT-estimator effective-collision
     * statistic, not a spectral-side signal-sufficiency measure (spectral
     * sample sizes are per-direction and independently 8/8/4, per {@code
     * LGA_DIRECTION}). No dedicated spectral low-signal metric exists to
     * substitute cheaply without inventing a new threshold, which would be
     * scope creep for this round -- the parameter is named for what it
     * ACTUALLY is (the transport flag reused as a proxy) rather than what
     * it is used for, so a reader of this signature is not misled into
     * thinking a spectral-specific signal check backs this branch.
     */
    static String ridgeVerdictText(boolean anyNaNRidge, boolean anyNonzeroRidge,
                                    boolean transportSmallNUsedAsSpectralProxy) {
        if (anyNaNRidge) {
            return "INSTRUMENT ANOMALY -- at least one pooled SPECTRAL direction's mean ridge slope is NaN, "
                   + "which is neither a valid RIDGE PRESENT nor RIDGE ABSENT determination; reported as an "
                   + "anomaly rather than silently defaulting to RIDGE ABSENT";
        }
        if (anyNonzeroRidge) {
            return "RIDGE PRESENT -- at least one pooled SPECTRAL direction has a nonzero ridge slope";
        }
        String measured = "RIDGE ABSENT -- every pooled SPECTRAL direction's ridge slope is exactly zero";
        if (transportSmallNUsedAsSpectralProxy) {
            return measured
                   + " -- MEASURED ONLY: this observable alone does NOT distinguish EITHER purely diffusive "
                   + "dynamics (no propagating branch, omega~i*D*k^2, AnisotropyProbe's own documented "
                   + "signature) OR an estimator with insufficient signal to fit (this campaign's own "
                   + "lgaSmallNEarlyTimeFlag above is set, and the campaign's excited-member/no-op "
                   + "characterization elsewhere in this report documents a near-static contact graph, the "
                   + "same signature AnisotropyProbe's own Javadoc names for the K=0 collision-free baseline "
                   + "\"by construction\") -- no dynamics claim is made from this line alone";
        }
        return measured
               + " -- consistent with the purely-diffusive-dynamics signature AnisotropyProbe's own Javadoc "
               + "documents (no propagating branch, omega~i*D*k^2); this campaign is comfortably above the "
               + "small-N threshold, so an insufficient-signal instrument state is a less likely alternative "
               + "reading here, though this line alone still does not itself rule it out";
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
        double phaseARatio = phaseA[0];
        double phaseALower = phaseA[1];
        double phaseAUpper = phaseA[2];
        double phaseAP = phaseA[3];
        double pcLower = pcTransport.pooledRatioCiLower();
        double pcUpper = pcTransport.pooledRatioCiUpper();
        double pcP = pcTransport.permutationPValue();
        boolean phaseAExcludes = ciExcludesOne(phaseALower, phaseAUpper);
        boolean phaseASignificant = isSignificant(phaseAP);
        boolean phaseCExcludes = ciExcludesOne(pcLower, pcUpper);
        boolean phaseCSignificant = isSignificant(pcP);
        sb.append("# ciVsPermutationReconciliation=the pooled TRANSPORT ratio's resample-then-aggregate bootstrap CI (bounded below by 1.0 by construction -- a max/min-of-3 order statistic, upward-biased by seed noise, AnisotropyProbe's own class Javadoc \"STACKED-REVIEW CORRECTION\") is NOT the significance statistic; permutationPValue IS, per that same locked convention. TWO-CAMPAIGN PATTERN: Phase A pooled TRANSPORT ratio=")
          .append(formatPrecise(phaseARatio)).append(" CI=[").append(formatPrecise(phaseALower))
          .append(',').append(formatPrecise(phaseAUpper)).append("] ")
          .append(ciExclusionClause(phaseAExcludes)).append(" permutationP=")
          .append(formatPrecise(phaseAP)).append(' ')
          .append(significanceClause(phaseASignificant))
          .append("; Phase C pooled TRANSPORT ratio=")
          .append(formatPrecise(pcTransport.pooledRatio().orElse(Double.NaN)))
          .append(" CI=[").append(formatPrecise(pcLower))
          .append(',').append(formatPrecise(pcUpper)).append("] ")
          .append(ciExclusionClause(phaseCExcludes)).append(" permutationP=")
          .append(formatPrecise(pcP)).append(' ')
          .append(significanceClause(phaseCSignificant))
          .append(" -- ")
          .append(ciVsPermutationConclusion(phaseAExcludes, phaseASignificant,
                                             phaseCExcludes, phaseCSignificant))
          .append('\n');
    }

    /**
     * C3 pure helper: whether a pooled ratio's CI excludes 1.0 (the
     * no-anisotropy null value) -- DERIVED from the actual CI bounds,
     * never assumed.
     */
    static boolean ciExcludesOne(double ciLower, double ciUpper) {
        return ciLower > 1.0 || ciUpper < 1.0;
    }

    /** C3 pure helper: the locked significance threshold, p&lt;=0.05, matching posture(iii)'s own convention. */
    static boolean isSignificant(double permutationP) {
        return permutationP <= 0.05;
    }

    /** C3 pure helper: label text for {@link #ciExcludesOne}, both branches. */
    static String ciExclusionClause(boolean excludesOne) {
        return excludesOne ? "(excludes 1.0)" : "(includes 1.0)";
    }

    /** C3 pure helper: label text for {@link #isSignificant}, both branches. */
    static String significanceClause(boolean significant) {
        return significant ? "(significant)" : "(NOT significant)";
    }

    /**
     * C3 pure helper (Critical fix): the two-campaign concluding sentence,
     * DERIVED from all four booleans rather than a hardcoded "in BOTH
     * campaigns ..." literal that would silently go stale (and, per the
     * report's own {@code powerRecommendationForGate}, the recommended
     * follow-up campaign is precisely the kind of re-run that could flip
     * this). Package-private for direct unit testing of every branch
     * combination with synthetic booleans.
     */
    static String ciVsPermutationConclusion(boolean phaseAExcludes, boolean phaseASignificant,
                                              boolean phaseCExcludes, boolean phaseCSignificant) {
        boolean phaseAAntiConservative = phaseAExcludes && !phaseASignificant;
        boolean phaseCAntiConservative = phaseCExcludes && !phaseCSignificant;
        if (phaseAAntiConservative && phaseCAntiConservative) {
            return "in BOTH campaigns the CI-excludes-1.0 signal does NOT correspond to permutation significance, "
                   + "consistent with this CI being anti-conservative at this campaign scale, not evidence of a "
                   + "real per-seed direction effect. A reader must not treat \"CI excludes 1.0\" as significance "
                   + "evidence here.";
        }
        if (phaseAAntiConservative != phaseCAntiConservative) {
            return "the two campaigns DISAGREE on the CI-excludes-1.0-without-permutation-significance pattern ("
                   + (phaseAAntiConservative ? "Phase A shows it, Phase C does not"
                                              : "Phase C shows it, Phase A does not")
                   + ") -- this reconciliation must be re-examined against the current numbers, not assumed to "
                   + "still hold uniformly.";
        }
        return "NEITHER campaign shows the CI-excludes-1.0-without-permutation-significance pattern this note "
               + "originally described -- re-examine whether this reconciliation is still applicable rather than "
               + "assuming it.";
    }

    /**
     * Critic relay item 4: a concrete, non-decisional recommendation for
     * .23 -- framed as an input to the gate, not a posture choice.
     *
     * <p>Fix-round Critical 1 correction: seed/tick scaling alone CANNOT
     * remove the absorbing-IC contamination this report's ANISOTROPY
     * CAMPAIGN EQUALITY characterization documents -- every background
     * member is seeded quanta=0 (a permanent no-op preimage under the
     * accepted diff=0 self-loop), so scaling seeds/ticks multiplies MORE
     * copies of the SAME frozen-contact regime, not a differently-powered
     * one. The seeds/ticks arithmetic below is kept (it is still the
     * correct power calculation for the CURRENT seeding), but is
     * explicitly SUBORDINATED to a seeding-design caveat: the
     * remediation this report can identify but not decide is a
     * non-absorbing background seeding (e.g. the long run's own random
     * nonzero quanta, which DO rotate) -- framed here as a USER design
     * decision for the .23 gate, not a selection this class makes.
     *
     * <p>Round-3 wording correction (reviewer R2-4): the seeding caveat
     * previously labelled itself "SUBORDINATE to the arithmetic above",
     * which inverts the finding -- the caveat GOVERNS the arithmetic (the
     * seeds/ticks scaling is valid only once the seeding contamination is
     * removed), not the reverse.
     */
    private static void appendPowerRecommendation(StringBuilder sb,
                                                     double meanEffectivePerSeed) {
        double phaseAMean = readPhaseAMeanEffectiveCollisionsPerSeed();
        double threshold = AnisotropyProbe.SMALL_N_EFFECTIVE_COLLISIONS_THRESHOLD;
        boolean phaseABelow = phaseAMean < threshold;
        boolean phaseCBelow = meanEffectivePerSeed < threshold;
        sb.append("# powerRecommendationForGate=Phase A (mean ")
          .append(formatPrecise(phaseAMean))
          .append(" effective collisions/seed) and Phase C (mean ")
          .append(formatPrecise(meanEffectivePerSeed))
          .append(") ").append(belowThresholdClause(phaseABelow, phaseCBelow))
          .append(" the informational small-N threshold (")
          .append(threshold)
          .append("). SEEDS/TICKS POWER ARITHMETIC (an input to the .23 gate, NOT a decision made here): "
                  + "a properly-powered follow-up anisotropy campaign, seeds 8->24 and ticks 128->400-500. "
                  + "SEEDING CAVEAT (fix-round Critical 1, GOVERNS the arithmetic above: the seeds/ticks "
                  + "scaling is valid only once this seeding contamination is removed, not optional "
                  + "context): this campaign's own background is seeded all-zero-quanta, an "
                  + "ABSORBING initial condition under the accepted diff=0 self-loop (see this report's "
                  + "ANISOTROPY CAMPAIGN EQUALITY characterization) -- seed/tick scaling alone CANNOT remove "
                  + "this contamination: scaling seeds/ticks multiplies the SAME frozen-contact regime, it "
                  + "does NOT power up a differently-behaved one. A "
                  + "non-absorbing background seeding (e.g. random nonzero quanta, as the long run already "
                  + "uses) is the remediation this report can identify but not decide -- that choice, and "
                  + "whether to combine it with the seeds/ticks scaling above, is a USER DESIGN DECISION for "
                  + "the .23 gate, not a selection made here. This whole recommendation is a DISTINCT "
                  + "decision point from the three isotropy postures in Section 5: the user may reasonably "
                  + "choose to defer any posture commitment until a properly-powered, non-absorbing-seeding "
                  + "campaign exists, rather than choosing among postures on the current underpowered, "
                  + "contaminated data.\n");
    }

    /**
     * C3 pure helper: which campaign(s) actually fall below the small-N
     * threshold, DERIVED from the two booleans rather than a hardcoded
     * "both ... fall below" literal. Package-private for direct unit
     * testing of all four combinations with synthetic booleans.
     */
    static String belowThresholdClause(boolean phaseABelow, boolean phaseCBelow) {
        if (phaseABelow && phaseCBelow) {
            return "BOTH fall below";
        }
        if (phaseABelow) {
            return "Phase A falls below (but Phase C does NOT fall below)";
        }
        if (phaseCBelow) {
            return "Phase C falls below (but Phase A does NOT fall below)";
        }
        return "NEITHER falls below";
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

    /**
     * C2/C3 fix: reads the committed Phase A artifact's OWN {@code
     * smallNEarlyTimeFlag} header line live (never hardcoded, replacing
     * the pre-fix literal {@code 29.25}) -- that header's prose carries
     * {@code mean effective collisions/seed=<value>} in scientific
     * notation, parsed here rather than re-measured.
     */
    private static double readPhaseAMeanEffectiveCollisionsPerSeed() {
        try {
            List<String> lines = Files.readAllLines(Paths.get(PHASE_A_RELATIVE_PATH));
            java.util.regex.Pattern pattern = java.util.regex.Pattern.compile("mean effective collisions/seed=([0-9.eE+-]+)");
            for (String line : lines) {
                if (line.startsWith("# smallNEarlyTimeFlag=")) {
                    java.util.regex.Matcher matcher = pattern.matcher(line);
                    if (matcher.find()) {
                        return Double.parseDouble(matcher.group(1));
                    }
                }
            }
            throw new IllegalStateException("Phase A artifact at " + PHASE_A_RELATIVE_PATH
                                             + " has no smallNEarlyTimeFlag header with a "
                                             + "\"mean effective collisions/seed=\" value");
        } catch (IOException e) {
            throw new IllegalStateException("failed to read Phase A artifact at "
                                             + PHASE_A_RELATIVE_PATH
                                             + " for the power recommendation",
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

    /**
     * Fix-round Critical 2: the raw per-CONTACT-direction collision
     * census -- previously collected (via {@link #lgaSubstrateFactory}'s
     * {@code perSeedStatsSink}) but only ever consumed for tick-quartile
     * TOTALS, never surfaced per-direction. Emits per-seed and pooled
     * census rows for the anisotropy campaign, a chi-squared-against-
     * uniform test on the pooled counts, and the long run's own
     * per-direction census for BOTH substrates plus the derived
     * hopping-rate tensor. EVIDENCE ROWS ONLY -- see {@code
     * contactDirectionCensusNote} below for the hard framing constraint.
     *
     * <p>Round-3 additions (critic S-NEW-2): the chi-squared's UNIFORM-null
     * premise (all six directions {@code <110>}-type / O_h-equivalent with
     * identical per-tick exposure) is now stated in-artifact rather than
     * left implicit in code Javadoc only; a double-underflow pooled
     * p-value is annotated rather than left as a bare, misreadable {@code
     * 0.000000000e+00}; and EFFECTIVE-transfer-only per-direction counts
     * (the population that actually moves quanta) are emitted alongside
     * the raw recorded-contact census.
     */
    private static void appendContactDirectionCensusSection(StringBuilder sb,
                                                                Report report) {
        sb.append("# === SECTION 1B: FCC CONTACT-DIRECTION CENSUS (evidence rows only, no posture -- see contactDirectionCensusNote) ===\n");
        sb.append("# columns(LGA_CONTACT_DIRECTION_SEED)=recordType\tseed\td1\td2\td3\td4\td5\td6\ttotal (raw per-CONTACT-direction collision counts, FccNeighborhood's positive directions 1..6 -- DISTINCT from LGA_DIRECTION's X100/X110/X111 k-space probe directions, an unrelated direction space)\n");
        sb.append("# columns(LGA_CONTACT_DIRECTION_POOLED)=recordType\td1\td2\td3\td4\td5\td6\ttotal\tchiSquare\tdf\tpValue (pooled across all 8 anisotropy-campaign seeds; standard Pearson chi-squared goodness-of-fit against a UNIFORM null over the 6 directions, computed live -- see chiSquareCaveat)\n");
        for (ContactDirectionCensus c : report.anisotropyContactDirectionPerSeed()) {
            appendContactDirectionRow(sb, "LGA_CONTACT_DIRECTION_SEED", c.scope(),
                                       c.perDirection());
        }
        java.util.Map<Integer, Long> pooled = pooledDirectionCounts(report.anisotropyContactDirectionPerSeed());
        long[] pooledArray = new long[6];
        for (int d = 1; d <= 6; d++) {
            pooledArray[d - 1] = pooled.get(d);
        }
        double chiSquare = chiSquareStatistic(pooledArray);
        int df = pooledArray.length - 1;
        double pValue = chiSquarePValue(chiSquare, df);
        sb.append("LGA_CONTACT_DIRECTION_POOLED");
        for (int d = 1; d <= 6; d++) {
            sb.append('\t').append(pooled.get(d));
        }
        long pooledTotal = 0;
        for (long v : pooledArray) {
            pooledTotal += v;
        }
        sb.append('\t').append(pooledTotal).append('\t')
          .append(formatPrecise(chiSquare)).append('\t').append(df).append('\t')
          .append(formatPrecise(pValue)).append('\n');
        sb.append("# uniformNullPremise=the chi-squared above tests the pooled counts against a UNIFORM "
                  + "expected count per direction; this is physically motivated, not arbitrary, because all six "
                  + "canonical FCC contact directions (FccNeighborhood's positive directions 1..6) are "
                  + "<110>-type / O_h-equivalent offsets (each a two-nonzero-coordinate unit vector related to "
                  + "the others by cubic point-group symmetry), and every even-parity member queries all six "
                  + "exactly once per tick (LatticeGasAutomaton's per-cell contact loop, direction=1..6) -- equal "
                  + "symmetry class and equal per-tick exposure are what make a uniform expected count the "
                  + "correct null under isotropy, rather than a default chosen for convenience\n");
        if (pValue == 0.0) {
            sb.append("# pValueUnderflowNote=the pValue column above (")
              .append(formatPrecise(pValue))
              .append(") is a DOUBLE UNDERFLOW, not a computed exact zero -- a chi-squared this large (chi2=")
              .append(formatPrecise(chiSquare)).append(", df=").append(df)
              .append(") drives the upper-incomplete-gamma evaluation below double's representable range, so "
                      + "the true value is p < 1.0e-300; reported here so the printed 0.0 is not misread as an "
                      + "exact computed zero\n");
        }
        sb.append("# chiSquareCaveat=this chi-squared tests uniformity of the RAW RECORDED-COLLISION EVENT counts, "
                  + "NOT independence-adjusted -- per this report's ANISOTROPY CAMPAIGN EQUALITY characterization, "
                  + "the anisotropy campaign's contact graph is near-static (same handful of frozen contacts fire "
                  + "every tick), so these ")
          .append(pooledTotal)
          .append(" recorded events are NOT independent draws; read this p-value alongside that characterization, "
                  + "not in isolation -- no independence-adjustment (e.g. an effective-sample-size deflation) is "
                  + "attempted here, since choosing one would itself be a judgment call this instrument does not make\n");

        sb.append("# columns(LGA_CONTACT_DIRECTION_EFFECTIVE_SEED)=recordType\tseed\td1\td2\td3\td4\td5\td6\ttotal "
                  + "(fix-round item 2c: per-seed EFFECTIVE-transfer-only counts, transferMagnitude>0 -- the narrower "
                  + "population that actually moves quanta and drives the hopping-tensor/D_hat estimate, additive-only "
                  + "from the same CollisionStatistics sinks as LGA_CONTACT_DIRECTION_SEED above)\n");
        sb.append("# columns(LGA_CONTACT_DIRECTION_EFFECTIVE_POOLED)=recordType\td1\td2\td3\td4\td5\td6\ttotal "
                  + "(pooled across all 8 anisotropy-campaign seeds; this total must equal ANISOTROPY_CAMPAIGN_EQUALITY_POOLED's "
                  + "effectiveCollisions column -- same underlying per-seed CollisionStatistics#effectiveCollisions(), just "
                  + "broken out by direction here)\n");
        for (ContactDirectionCensus c : report.anisotropyEffectiveContactDirectionPerSeed()) {
            appendContactDirectionRow(sb, "LGA_CONTACT_DIRECTION_EFFECTIVE_SEED", c.scope(),
                                       c.perDirection());
        }
        java.util.Map<Integer, Long> pooledEffective = pooledDirectionCounts(report.anisotropyEffectiveContactDirectionPerSeed());
        sb.append("LGA_CONTACT_DIRECTION_EFFECTIVE_POOLED");
        long pooledEffectiveTotal = 0;
        for (int d = 1; d <= 6; d++) {
            long v = pooledEffective.getOrDefault(d, 0L);
            sb.append('\t').append(v);
            pooledEffectiveTotal += v;
        }
        sb.append('\t').append(pooledEffectiveTotal).append('\n');

        sb.append("# columns(LONG_RUN_CONTACT_DIRECTION)=recordType\tsubstrate\td1\td2\td3\td4\td5\td6\ttotal (the SEPARATE Section 3/4 long run's own per-CONTACT-direction census, both substrates)\n");
        sb.append("# columns(HOPPING_TENSOR)=recordType\tsubstrate\tDxx\tDyy\tDzz\tDxz\tDyz\tDxy (derived from LONG_RUN_CONTACT_DIRECTION via D = sum_d w_d * outer(e_d,e_d), e_d = FccNeighborhood's unit offset for positive direction d)\n");
        appendContactDirectionRow(sb, "LONG_RUN_CONTACT_DIRECTION",
                                   report.longRunHybridContactDirection().scope(),
                                   report.longRunHybridContactDirection().perDirection());
        appendContactDirectionRow(sb, "LONG_RUN_CONTACT_DIRECTION",
                                   report.longRunLgaContactDirection().scope(),
                                   report.longRunLgaContactDirection().perDirection());
        appendHoppingTensorRow(sb, report.longRunHybridContactDirection());
        appendHoppingTensorRow(sb, report.longRunLgaContactDirection());

        sb.append("# contactDirectionCensusNote=this is a FIRST-ORDER directional-structure observable (which of "
                  + "the six FCC contact directions collisions actually occur in), MECHANISM-DISTINCT from the "
                  + "rank-4/FCHC risk the design memo names (a lattice can be exact-cubic-symmetric at 4th order "
                  + "and still show first-order structure here, or vice versa) -- EVIDENCE ONLY, no posture is "
                  + "selected or advocated from this census; interpretation is left to Section 5. The anisotropy "
                  + "campaign in this bead runs the LGA substrate ONLY (see lgaSubstrateFactory) -- no hybrid "
                  + "anisotropy run exists in this bead to compare the LGA_CONTACT_DIRECTION rows against; the "
                  + "LONG_RUN_CONTACT_DIRECTION rows above are the only both-substrate comparison this report can "
                  + "make, and are a DIFFERENT campaign (Section 3/4's shared long run, not the anisotropy packet "
                  + "campaign).\n");
        sb.append('\n');
    }

    private static void appendContactDirectionRow(StringBuilder sb, String recordType,
                                                     String scope,
                                                     java.util.Map<Integer, Long> perDirection) {
        sb.append(recordType).append('\t').append(scope);
        long total = 0;
        for (int d = 1; d <= 6; d++) {
            long v = perDirection.getOrDefault(d, 0L);
            sb.append('\t').append(v);
            total += v;
        }
        sb.append('\t').append(total).append('\n');
    }

    private static void appendHoppingTensorRow(StringBuilder sb,
                                                  ContactDirectionCensus census) {
        double[] tensor = hoppingTensor(census.perDirection());
        sb.append("HOPPING_TENSOR\t").append(census.scope());
        for (double component : tensor) {
            sb.append('\t').append(formatPrecise(component));
        }
        sb.append('\n');
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

    // ------------------------------------------------------------------
    // Fix-round Critical 2: pooled per-CONTACT-direction census + a
    // standard Pearson chi-squared goodness-of-fit test against a uniform
    // null over the six directions. FRAMING (hard constraint, relay item
    // 2): these are EVIDENCE rows/facts, not a posture -- the pooled
    // counts are RAW recorded-collision events, not independent draws
    // (see ANISOTROPY_CAMPAIGN_EQUALITY above: the same handful of frozen
    // contacts fire every tick), so this chi-squared answers "are the raw
    // recorded EVENTS uniform across direction", not "is the underlying
    // physical process isotropic" -- the two questions are NOT the same
    // when observations repeat a near-static contact graph, and this
    // instrument does not conflate them.
    // ------------------------------------------------------------------

    /** Pure helper: sums a list of per-seed censuses into one pooled census, directions 1..6 only. */
    static java.util.Map<Integer, Long> pooledDirectionCounts(List<ContactDirectionCensus> perSeed) {
        java.util.Map<Integer, Long> pooled = new java.util.TreeMap<>();
        for (int d = 1; d <= 6; d++) {
            pooled.put(d, 0L);
        }
        for (ContactDirectionCensus c : perSeed) {
            for (var e : c.perDirection().entrySet()) {
                pooled.merge(e.getKey(), e.getValue(), Long::sum);
            }
        }
        return pooled;
    }

    /**
     * Pearson chi-squared statistic against a UNIFORM null: {@code
     * sum((observed-expected)^2/expected)}, {@code expected = total/n}.
     * Pure, deterministic; {@code observed} must be non-empty with a
     * positive total (both true by construction for a 6-direction
     * pooled census with any recorded collisions).
     */
    static double chiSquareStatistic(long[] observed) {
        double total = 0;
        for (long o : observed) {
            total += o;
        }
        double expected = total / observed.length;
        double stat = 0;
        for (long o : observed) {
            double diff = o - expected;
            stat += diff * diff / expected;
        }
        return stat;
    }

    /**
     * Upper-tail p-value of a chi-squared statistic with {@code df}
     * degrees of freedom -- {@code Q(df/2, chiSquare/2)}, the regularized
     * upper incomplete gamma function (standard chi-squared survival
     * function). Implemented locally (Numerical-Recipes-style series +
     * continued-fraction evaluation of the incomplete gamma function)
     * rather than pulling in a statistics dependency for one function --
     * see class Javadoc "Spartan Design" convention. Deterministic, no
     * RNG. Verified against textbook chi-squared critical values in
     * {@code PhaseCMeasurementTest} (df=1 at x=3.841459 -> p~0.05, df=5
     * at x=11.0705 -> p~0.05).
     */
    static double chiSquarePValue(double chiSquare, int df) {
        if (chiSquare < 0) {
            throw new IllegalArgumentException("chiSquare must be non-negative, was: "
                                                + chiSquare);
        }
        if (df <= 0) {
            throw new IllegalArgumentException("df must be positive, was: " + df);
        }
        if (chiSquare == 0) {
            return 1.0;
        }
        return regularizedUpperIncompleteGamma(df / 2.0, chiSquare / 2.0);
    }

    /** {@code Q(a,x) = 1 - P(a,x)}, via series ({@code x < a+1}) or continued fraction ({@code x >= a+1}). */
    private static double regularizedUpperIncompleteGamma(double a, double x) {
        if (x < a + 1.0) {
            return 1.0 - lowerIncompleteGammaSeries(a, x);
        }
        return upperIncompleteGammaContinuedFraction(a, x);
    }

    private static double lowerIncompleteGammaSeries(double a, double x) {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int n = 1; n <= 200; n++) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (Math.abs(del) < Math.abs(sum) * 1e-15) {
                break;
            }
        }
        return sum * Math.exp(-x + a * Math.log(x) - logGamma(a));
    }

    private static double upperIncompleteGammaContinuedFraction(double a, double x) {
        double tiny = 1e-300;
        double b = x + 1.0 - a;
        double c = 1.0 / tiny;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= 200; i++) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (Math.abs(d) < tiny) {
                d = tiny;
            }
            c = b + an / c;
            if (Math.abs(c) < tiny) {
                c = tiny;
            }
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) < 1e-15) {
                break;
            }
        }
        return Math.exp(-x + a * Math.log(x) - logGamma(a)) * h;
    }

    private static final double[] LANCZOS_COEFFICIENTS = { 76.18009172947146,
                                                             -86.50532032941677,
                                                             24.01409824083091,
                                                             -1.231739572450155,
                                                             0.1208650973866179e-2,
                                                             -0.5395239384953e-5 };

    /** Lanczos-approximation log-gamma, standard textbook coefficients. */
    private static double logGamma(double x) {
        double y = x;
        double tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        double ser = 1.000000000190015;
        for (double c : LANCZOS_COEFFICIENTS) {
            y += 1.0;
            ser += c / y;
        }
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }

    /**
     * The hopping-rate tensor {@code D = sum_d w_d * outer(e_d,e_d)} over
     * the six positive FCC contact directions, {@code e_d} = {@link
     * FccNeighborhood#offsetOf}'s unit vector, {@code w_d} = that
     * direction's collision-census weight. Returns {@code [Dxx, Dyy, Dzz,
     * Dxz, Dyz, Dxy]}. A first-order directional-structure observable --
     * see the artifact's own {@code contactDirectionCensusNote} for the
     * framing constraint (evidence, not a posture).
     */
    static double[] hoppingTensor(java.util.Map<Integer, Long> perDirection) {
        double dxx = 0, dyy = 0, dzz = 0, dxz = 0, dyz = 0, dxy = 0;
        for (int direction = 1; direction <= 6; direction++) {
            Point3i offset = FccNeighborhood.offsetOf(direction);
            double norm = Math.sqrt((double) offset.x * offset.x
                                     + (double) offset.y * offset.y
                                     + (double) offset.z * offset.z);
            double ex = offset.x / norm;
            double ey = offset.y / norm;
            double ez = offset.z / norm;
            long w = perDirection.getOrDefault(direction, 0L);
            dxx += w * ex * ex;
            dyy += w * ey * ey;
            dzz += w * ez * ez;
            dxz += w * ex * ez;
            dyz += w * ey * ez;
            dxy += w * ex * ey;
        }
        return new double[] { dxx, dyy, dzz, dxz, dyz, dxy };
    }

    private static void appendEquilibriumSection(StringBuilder sb, Report report) {
        sb.append("# === SECTION 4: EQUILIBRIUM QUANTA-DISTRIBUTION CHARACTERIZATION (SDB caveat, user decision 6) ===\n");
        sb.append("# equilibriumStatisticsAssumeSDB=false (explicit -- CollisionTable's collision rule is an ACCEPTED non-bijection, see class Javadoc)\n");
        sb.append("# columns(QUANTA_HISTOGRAM)=recordType\tlabel\tn\tmean\tvariance\tmin\tmax (distributional SUMMARY only -- fix-round S5: this row does NOT itself carry the histogram despite its name; see QUANTA_HISTOGRAM_BIN for the actual per-value population counts)\n");
        sb.append("# columns(QUANTA_HISTOGRAM_BIN)=recordType\tlabel\tvalue\tcount (exact per-value population count, one row per integer value in [min,max] inclusive, zero-filled for unobserved values -- localizes WHERE QUANTA_HISTOGRAM's variance drop concentrates, e.g. at value 0 for the SDB diff=0 self-loop, rather than leaving that inferred from variance alone)\n");
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
        double noOpFirst = noOpFractions.get(0);
        double noOpLast = noOpFractions.get(noOpFractions.size() - 1);
        boolean noOpPeakIsInterior = noOpPeakWindow != 0
                                      && noOpPeakWindow != noOpFractions.size() - 1;
        sb.append("# equalityConcentrationEmpiricalTrend=lgaNoOpFractionByWindow=")
          .append(formatSeries(noOpFractions))
          .append(", peak at window ").append(noOpPeakWindow).append(" (")
          .append(formatPrecise(noOpFractions.get(noOpPeakWindow)))
          .append(")")
          .append(" -- ").append(trendShapeStatement(noOpFirst, noOpLast))
          .append(" (first->last: ").append(formatPrecise(noOpFirst))
          .append(" -> ").append(formatPrecise(noOpLast))
          .append(")").append(sdbCorroborationClause(noOpFirst, noOpLast))
          .append(noOpPeakIsInterior
                  ? ", but peaks at an INTERIOR window (" + noOpPeakWindow
                    + ") rather than at the trend's own endpoint"
                  : ", consistent with a monotonic trend across the full window series")
          .append(". NOTE: ")
          .append(lagNote(noOpPeakWindow, gapPeakWindow))
          .append(" -- the two trends are NOT the same mechanism showing up identically; stated explicitly rather than implying a single coupled cause.\n");
        sb.append('\n');
    }

    /**
     * Post-C1-fix pure helper (Critical): states whether {@code
     * noOpPeakWindow} (the no-op-fraction trend's OWN peak -- the
     * grammatical SUBJECT) occurs ONE WINDOW LATER than {@code
     * gapPeakWindow} (Section 3's relativeRateGapByWindow peak), matching
     * the guard's exact semantics ({@code noOpPeakWindow ==
     * gapPeakWindow + 1} means the NO-OP trend peaks LATER). The pre-fix
     * bug made the GAP series the subject instead and inverted the
     * reported direction. Package-private (not {@code private}) so it can
     * be unit-tested directly with synthetic window indices covering both
     * branches -- exactly the coverage gap (Important #2) that hid the
     * original inversion from every existing test.
     */
    static String lagNote(int noOpPeakWindow, int gapPeakWindow) {
        return "this trend's own peak (window " + noOpPeakWindow + ") occurs "
               + (noOpPeakWindow == gapPeakWindow + 1
                  ? "exactly ONE WINDOW LATER than"
                  : "at a different window than")
               + " Section 3's relativeRateGapByWindow peak (window " + gapPeakWindow + ")";
    }

    /**
     * The SDB-corroboration clause must track the same first-vs-last
     * comparison as {@link #trendShapeStatement}: the closed form
     * predicts a RISING no-op fraction, so only a rising trend
     * corroborates it -- a falling or flat trend must say so instead
     * of silently claiming corroboration.
     */
    static String sdbCorroborationClause(double first, double last) {
        if (last > first) {
            return ", dynamical corroboration of the SDB closed-form's diff=0 preimage concentration (pairs settling into the permanent equality self-loop)";
        }
        return ", NOT corroborating the SDB closed-form's diff=0 preimage-concentration prediction (the closed form predicts a RISING no-op fraction; this series is not rising first->last)";
    }

    /**
     * C3 pure helper: the trend-shape headline for a first-vs-last
     * comparison, DERIVED from the data rather than the unconditional
     * "OVERALL RISING" literal the pre-fix code always emitted regardless
     * of the actual series.
     */
    static String trendShapeStatement(double first, double last) {
        if (last > first) {
            return "OVERALL RISING";
        } else if (last < first) {
            return "OVERALL FALLING";
        } else {
            return "OVERALL FLAT (first == last)";
        }
    }

    private static void appendHistogramRow(StringBuilder sb,
                                              QuantaHistogramSummary h) {
        sb.append("QUANTA_HISTOGRAM\t").append(h.label()).append('\t')
          .append(h.n()).append('\t').append(formatPrecise(h.mean())).append('\t')
          .append(formatPrecise(h.variance())).append('\t').append(h.min())
          .append('\t').append(h.max()).append('\n');
        for (var e : h.bins().entrySet()) {
            sb.append("QUANTA_HISTOGRAM_BIN\t").append(h.label()).append('\t')
              .append(e.getKey()).append('\t').append(e.getValue()).append('\n');
        }
    }

    /**
     * Fix-round Critical 1: characterizes the SDB equality-concentration
     * ON THE ANISOTROPY CAMPAIGN ITSELF (extent 8,8,8 / 8 seeds / 127
     * ticks / packet-in-vacuum IC) -- Section 4 above characterizes it
     * only on the SEPARATE long run (extent 4,4,4 / random background),
     * where the effect is mild; this campaign's all-zero-quanta
     * background makes it the campaign's own INITIAL CONDITION and, under
     * the accepted diff=0 self-loop, ABSORBING -- a zero-quanta member's
     * phase never advances, so a {@code (0,0)} contact is a PERMANENT,
     * not merely likely, no-op. Every field here is measured from the
     * additive census sinks {@link #lgaSubstrateFactory} populates (see
     * {@link AnisotropyCampaignEquality}'s Javadoc) -- nothing here is
     * assumed or re-derived from the long run's different campaign.
     *
     * <p>Round-3 corrections: (R2-1) {@code activeMemberSlotsTotal} is now
     * scoped to the WHOLE pooled campaign ({@code
     * AnisotropyProbe.DEFAULT_SEEDS.length} lattices' worth of even-parity
     * member slots), not one lattice, matching the pooled {@code
     * excitedMembers*Total} numerator it is divided against in the
     * narrative below -- the prior single-lattice denominator understated
     * the campaign's contamination 8x. (R2-2) the collisionsPerTick
     * flatness claim is now derived from the measured per-seed
     * first-vs-last delta ({@link #tickFlatnessStatement}) rather than an
     * unconditional "FLAT" literal (seed 47 goes 23-&gt;21, which the old
     * literal did not reflect).
     */
    private static void appendAnisotropyCampaignEqualitySection(StringBuilder sb,
                                                                    Report report) {
        sb.append("# === SECTION 4B: ANISOTROPY CAMPAIGN EQUALITY CHARACTERIZATION (fix-round Critical 1 -- the SDB gate item, ON THE CAMPAIGN THE ISOTROPY STATISTIC ACTUALLY COMES FROM) ===\n");
        sb.append("# columns(ANISOTROPY_CAMPAIGN_EQUALITY_SEED)=recordType\tseed\ttotalCollisions\teffectiveCollisions\tnoOpFraction\texcitedMembersInitial\tcellsTouchedInitial\texcitedMembersFinal\tcellsTouchedFinal\tcollisionsAtFirstTick\tcollisionsAtLastTick\n");
        sb.append("# columns(ANISOTROPY_CAMPAIGN_EQUALITY_POOLED)=recordType\ttotalCollisions\teffectiveCollisions\tnoOpFraction\tactiveMemberSlotsTotal\tevenParityCellsTotal\texcitedMembersInitialTotal\texcitedMembersFinalTotal\n");
        long pooledTotal = 0;
        long pooledEffective = 0;
        long pooledExcitedInitial = 0;
        long pooledExcitedFinal = 0;
        long[] firstTicks = new long[report.anisotropyCampaignEquality().size()];
        long[] lastTicks = new long[report.anisotropyCampaignEquality().size()];
        int seedIdx = 0;
        for (AnisotropyCampaignEquality e : report.anisotropyCampaignEquality()) {
            sb.append("ANISOTROPY_CAMPAIGN_EQUALITY_SEED\t").append(e.seed())
              .append('\t').append(e.totalCollisions()).append('\t')
              .append(e.effectiveCollisions()).append('\t')
              .append(formatPrecise(e.noOpFraction())).append('\t')
              .append(e.excitedMembersInitial()).append('\t')
              .append(e.cellsTouchedInitial()).append('\t')
              .append(e.excitedMembersFinal()).append('\t')
              .append(e.cellsTouchedFinal()).append('\t')
              .append(e.collisionsAtFirstTick()).append('\t')
              .append(e.collisionsAtLastTick()).append('\n');
            pooledTotal += e.totalCollisions();
            pooledEffective += e.effectiveCollisions();
            pooledExcitedInitial += e.excitedMembersInitial();
            pooledExcitedFinal += e.excitedMembersFinal();
            firstTicks[seedIdx] = e.collisionsAtFirstTick();
            lastTicks[seedIdx] = e.collisionsAtLastTick();
            seedIdx++;
        }
        double pooledNoOpFraction = pooledTotal == 0L ? Double.NaN
                                                        : (double) (pooledTotal
                                                                     - pooledEffective)
                                                          / pooledTotal;
        // R2-1: the denominator must be scoped to the WHOLE pooled
        // campaign (every seed's own lattice of even-parity member slots),
        // not one lattice -- the numerator (pooledExcitedInitial/Final) is
        // already summed across all AnisotropyProbe.DEFAULT_SEEDS.length
        // seeds.
        long activeMemberSlotsTotal = evenParityCellCount(AnisotropyProbe.DEFAULT_EXTENT)
                                       * 30L
                                       * AnisotropyProbe.DEFAULT_SEEDS.length;
        long evenParityCellsTotal = evenParityCellCount(AnisotropyProbe.DEFAULT_EXTENT);
        sb.append("ANISOTROPY_CAMPAIGN_EQUALITY_POOLED\t").append(pooledTotal)
          .append('\t').append(pooledEffective).append('\t')
          .append(formatPrecise(pooledNoOpFraction)).append('\t')
          .append(activeMemberSlotsTotal).append('\t')
          .append(evenParityCellsTotal).append('\t')
          .append(pooledExcitedInitial).append('\t')
          .append(pooledExcitedFinal).append('\n');
        sb.append("# anisotropyCampaignEqualityCharacterization=pooled no-op fraction=")
          .append(formatPrecise(pooledNoOpFraction)).append(" (")
          .append(pooledTotal - pooledEffective).append('/').append(pooledTotal)
          .append("), excited members ").append(pooledExcitedInitial)
          .append("/").append(activeMemberSlotsTotal).append(" at t=0 -> ")
          .append(pooledExcitedFinal).append("/").append(activeMemberSlotsTotal)
          .append(" after the campaign's full 127-tick run, collisionsPerTick is ")
          .append(tickFlatnessStatement(firstTicks, lastTicks))
          .append(" (first-tick vs last-tick counts in ANISOTROPY_CAMPAIGN_EQUALITY_SEED above) -- "
                  + "EQUALITY IS THIS CAMPAIGN'S OWN INITIAL CONDITION AND IS ABSORBING under the accepted "
                  + "diff=0 self-loop (CollisionTable's non-bijective rule, Section 4's sdbClosedForm above): "
                  + "every background member is seeded quanta=0, and phase[i]=floorMod(phase[i]+quanta[i],P) "
                  + "means a zero-quanta member NEVER rotates, so a (0,0) contact is a PERMANENT no-op, not a "
                  + "transient equilibrium approach -- this is a DIFFERENT and MORE SEVERE characterization "
                  + "than Section 4's long-run trend (variance contraction from a WIDE random background), and "
                  + "it is THIS campaign's data the isotropy statistic in Section 1/5 is computed from, not the "
                  + "long run's.\n");
        sb.append('\n');
    }

    /**
     * Reviewer R2-2 fix: "collisionsPerTick is FLAT per seed" was an
     * unconditional hardcoded literal even though seed 47's
     * collisionsAtFirstTick/collisionsAtLastTick goes 23-&gt;21 (a C3-class
     * defect -- a hardcoded conclusion not derived from the data -- at a
     * new site). Branches on the MEASURED max {@code |first-last|} delta
     * across seeds; the only non-arbitrary "flat" cutoff that does not
     * itself require an invented tolerance is EXACT equality across every
     * seed, so any nonzero delta is reported as the measured spread
     * instead of a pass/fail judgment call.
     */
    static String tickFlatnessStatement(long[] firstTicks, long[] lastTicks) {
        if (firstTicks.length != lastTicks.length) {
            throw new IllegalArgumentException("firstTicks/lastTicks length mismatch: "
                                                + firstTicks.length + " vs "
                                                + lastTicks.length);
        }
        long maxAbsDelta = 0L;
        for (int i = 0; i < firstTicks.length; i++) {
            maxAbsDelta = Math.max(maxAbsDelta,
                                    Math.abs(lastTicks[i] - firstTicks[i]));
        }
        if (maxAbsDelta == 0L) {
            return "EXACTLY FLAT per seed (first-tick == last-tick for every seed)";
        }
        return "NOT EXACTLY FLAT per seed (max |first-tick - last-tick| delta across seeds = "
               + maxAbsDelta + ")";
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
          .append(" -- CONTAMINATION CAVEAT (fix-round Critical 1): this evidence carries its own "
                  + "contamination, not merely the small-N caveat above -- see this report's ANISOTROPY "
                  + "CAMPAIGN EQUALITY characterization (the campaign's all-zero-quanta background is an "
                  + "ABSORBING initial condition under the accepted diff=0 self-loop, ~99% no-op fraction, "
                  + "<1% of members ever excited), which the permutation null is calibrated FROM, not "
                  + "independently of -- a null calibrated on a contaminated realization cannot itself "
                  + "certify power, so \"NOT significant\" here should be read alongside that "
                  + "characterization, not in isolation")
          .append('\n');
        sb.append("# section1BCrossReference=see Section 1B rows LGA_CONTACT_DIRECTION_*/HOPPING_TENSOR and their "
                  + "premise note (uniformNullPremise/contactDirectionCensusNote) for a first-order FCC "
                  + "contact-direction observable, MECHANISM-DISTINCT from the postures above -- evidence only, "
                  + "no posture is selected or advocated by this pointer\n");
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
