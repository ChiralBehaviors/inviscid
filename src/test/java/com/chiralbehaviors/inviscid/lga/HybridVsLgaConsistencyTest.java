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

package com.chiralbehaviors.inviscid.lga;

import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.BeforeClass;
import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.measure.CollisionStatistics;

/**
 * Tests 7-8 for bead inviscid-0nx.21: the hybrid automaton (Phase A,
 * continuous geometric contact detection) versus {@link
 * LatticeGasAutomaton} (Phase C, table-driven), from a SHARED initial
 * condition, both built against the SAME committed atlas.
 *
 * <h2>Cadence (user's 2A decision, FINAL)</h2>
 * One LGA tick == one hybrid tick; both substrates report {@code
 * phaseResolution() == 3600} (pin {@code
 * QuantaFieldSeamTest.bothSubstratesReportPhaseResolution3600UnderCadence2A}).
 * T2 design-ckn-lattice-seam.md §10 R10 confirms this test is "safe as
 * written" under 2A -- no respecification needed.
 *
 * <h2>Comparison resolution</h2>
 * Compared at the CONTACT-TABLE's own N_lga (24) bin resolution -- per
 * the bead's own wording, "trajectories agree for as long as
 * quantisation error stays below one bin" -- using the SAME {@link
 * PhaseQuantizer} the atlas was transcribed with, not the finer 3600
 * fine-phase resolution. This was empirically necessary, not a stylistic
 * choice: an earlier draft of this test compared fine phase directly
 * (exact 3600-resolution equality) and found "divergence" at tick 1 for
 * every seed tried -- not a real physics signal, but {@code
 * QUANTUM_RATE == (float)(2*pi/3600)}'s float32 rounding occasionally
 * landing a hybrid member's angle within noise of a fine-phase unit
 * boundary, which is essentially guaranteed to happen SOMEWHERE among
 * 1920 slots within one tick and tells you nothing about contact-table
 * fidelity. A bin is 150 fine-phase units wide; float32 rounding noise
 * (~1e-4 fine-phase units) is ~6 orders of magnitude too small to flip a
 * bin boundary on its own, so bin-level comparison isolates the signal
 * this test actually exists to measure: whether the LGA's table-driven
 * collision decisions track the hybrid's live geometric ones.
 *
 * <h2>Expected divergence mechanism (stated, not assumed)</h2>
 * The committed atlas's {@code overlapFraction > 0} ANY-OVERLAP
 * transcription (bead inviscid-gyt) fires a bin the moment ANY fraction
 * of its fine geometric sub-sweep contacts -- a strict superset of the
 * continuous predicate's true contact region. The LGA can therefore fire
 * a collision at a bin the hybrid's live, continuous {@code
 * ContactPredicate} evaluation would not (yet) fire at, the moment
 * either member's angle enters such a bin. Once one collision decision
 * differs, the two trajectories' quanta -- and hence all subsequent
 * phase advancement -- diverge permanently. This is a property of the
 * transcription's coarseness, not a bug in either substrate.
 *
 * @author halhildebrand
 */
public class HybridVsLgaConsistencyTest {

    private static final String RESOURCE_PATH = "lga/contact-atlas-v2.tsv";

    private static ContactAtlas   ATLAS;
    private static CollisionTable COLLISIONS;
    private static Point3i        EXTENT;

    @BeforeClass
    public static void loadFixtures() throws IOException {
        URL resource = HybridVsLgaConsistencyTest.class.getClassLoader()
                                                          .getResource(RESOURCE_PATH);
        Path path = Paths.get(resource.getPath());
        ATLAS = ContactAtlas.read(path);
        COLLISIONS = CollisionTable.buildFromPhaseARule(new QuantaExchangeRule());
        EXTENT = ATLAS.header().extent();
    }

    /** One shared initial condition, applied identically to both substrates. */
    private record SharedInitialCondition(int[] phase, long[] quanta) {
    }

    private static SharedInitialCondition sharedInitialCondition(long seed,
                                                                   int quantaBound) {
        Random random = new Random(seed);
        int length = 30 * EXTENT.x * EXTENT.y * EXTENT.z;
        int[] phase = new int[length];
        long[] quanta = new long[length];
        for (int i = 0; i < length; i++) {
            phase[i] = random.nextInt(3600);
            quanta[i] = random.nextInt(2 * quantaBound + 1) - quantaBound;
        }
        return new SharedInitialCondition(phase, quanta);
    }

    private static HybridAutomaton newHybrid(SharedInitialCondition ic,
                                              CollisionStatistics statistics) {
        Necronomata automaton = new Necronomata(EXTENT);
        int length = ic.phase().length;
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < length; i++) {
                angle[i] = (float) (2.0 * Math.PI * ic.phase()[i] / 3600.0);
                frequency[i] = ic.quanta()[i];
            }
        });
        FccNeighborhood neighborhood = new FccNeighborhood(EXTENT);
        ContactPredicate predicate = new ContactPredicate(new MemberGeometry(ATLAS.header()
                                                                                    .geometryResolution(),
                                                                               ATLAS.header()
                                                                                    .memberRadius()));
        ContactScan scan = new ContactScan(automaton, neighborhood, predicate);
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        return new HybridAutomaton(automaton, sweep);
    }

    private static LatticeGasAutomaton newLga(SharedInitialCondition ic,
                                               CollisionStatistics statistics) {
        LatticeGasAutomaton lga = new LatticeGasAutomaton(EXTENT, ATLAS,
                                                            COLLISIONS,
                                                            statistics);
        lga.process((phase, quanta) -> {
            System.arraycopy(ic.phase(), 0, phase, 0, ic.phase().length);
            System.arraycopy(ic.quanta(), 0, quanta, 0, ic.quanta().length);
        });
        return lga;
    }

    private static int hybridBin(Necronomata automaton, PhaseQuantizer quantizer,
                                  int slot) {
        float[][] box = new float[1][];
        automaton.process((angle, frequency, deltaA, deltaF) -> box[0] = angle);
        return quantizer.bin(box[0][slot]);
    }

    private static int lgaBin(LatticeGasAutomaton lga, int slot) {
        int[][] box = new int[1][];
        lga.process((phase, quanta) -> box[0] = phase);
        return box[0][slot] / lga.subBinSteps();
    }

    /**
     * @return the first tick (1-indexed count of ticks RUN) at which any
     *         member's N_lga BIN disagrees between the two substrates -
     *         see class Javadoc, "Comparison resolution" - or {@code
     *         maxTicks} if they agreed for the entire window (never
     *         expected to be reached here - see test assertions).
     */
    private static int firstDivergenceTick(HybridAutomaton hybrid,
                                            LatticeGasAutomaton lga,
                                            PhaseQuantizer quantizer,
                                            int maxTicks) {
        int length = lga.slotCount();
        for (int tick = 0; tick < maxTicks; tick++) {
            hybrid.tick(tick);
            lga.tick(tick);
            for (int slot = 0; slot < length; slot++) {
                int hybridBin = hybridBin(hybrid.automaton(), quantizer, slot);
                int lgaBin = lgaBin(lga, slot);
                if (hybridBin != lgaBin) {
                    return tick + 1;
                }
            }
        }
        return maxTicks;
    }

    /**
     * Test 7. Divergence is asserted to actually OCCUR within the test's
     * window (proving the two substrates are not literally the same
     * model, not a vacuously-loose threshold that always "passes" by
     * never being exercised).
     *
     * <p><b>No lower-bound assertion (final-review Significant D fix, T2
     * critique-final-0nx21-automaton-arc.md [21949]).</b> An earlier
     * version of this test additionally asserted {@code divergenceTick
     * >= 1}. That assertion was STRUCTURALLY VACUOUS, not merely weak:
     * {@link #firstDivergenceTick}'s own control flow starts at {@code
     * tick=0} and returns {@code tick+1} on the first mismatch (or {@code
     * maxTicks} on none), so the return value can never be less than 1
     * for ANY implementation, correct or badly broken -- the assertion
     * would pass identically either way. Dropped rather than kept as
     * dead ceremony.
     *
     * <p><b>The empirical finding itself is still true and still
     * reported, just not asserted as a numeric floor.</b> Measured
     * divergence at tick 1 for this seed/config both before and after
     * the fine-step contact-firing revision (bead inviscid-0nx.21's
     * fine-step fix improved test 8's AGGREGATE agreement dramatically,
     * but did not change how quickly individual TRAJECTORIES first
     * disagree at N_lga=24 bin resolution -- with 1920 independently
     * seeded slots, some pair's contact decision differs from the live
     * hybrid's within the first tick essentially always). This is
     * exactly why test 8 -- not this test -- carries the claim that
     * actually matters for Phase C: "same physics, not the same
     * trajectory."
     */
    @Test
    public void lgaAndHybridAgreeOnASharedInitialConditionForTheFirstNTicks() {
        SharedInitialCondition ic = sharedInitialCondition(42L, 6);
        CollisionStatistics hybridStats = new CollisionStatistics();
        CollisionStatistics lgaStats = new CollisionStatistics();
        HybridAutomaton hybrid = newHybrid(ic, hybridStats);
        LatticeGasAutomaton lga = newLga(ic, lgaStats);
        PhaseQuantizer quantizer = PhaseQuantizer.of(ATLAS.header());

        int maxTicks = 64;
        int divergenceTick = firstDivergenceTick(hybrid, lga, quantizer,
                                                  maxTicks);

        assertTrue("expected divergence to actually occur within the "
                   + maxTicks
                   + "-tick window (a threshold that never fires is vacuous) - "
                   + "observed divergence at tick " + divergenceTick,
                   divergenceTick < maxTicks);
    }

    /**
     * Effective/total collision-ratio agreement tolerance for test 8 -
     * UNCHANGED across the bin-level-firing -> fine-step-firing fix (see
     * that test's Javadoc, "HISTORY" / "CURRENT"): the fix closed the gap,
     * not the tolerance.
     */
    private static final double TOLERANCE = 0.15;

    /**
     * Test 8. The claim that actually matters: the LGA is the same
     * PHYSICS as the hybrid, not the same trajectory. Aggregate collision
     * statistics (effective/total ratio -- T2
     * analysis-73v-spectral-conversion-and-cadence.md §6's recommended
     * cadence-correctness discriminator, since {@code phaseResolution()}
     * equality can no longer distinguish the substrates under 2A) agree
     * within a documented tolerance long past the trajectory divergence
     * point test 7 found.
     *
     * <p><b>HISTORY (bin-level firing, superseded -- kept as the
     * quantified motivation for fine-step semantics, not dead text).</b>
     * When {@link LatticeGasAutomaton} fired contacts directly off {@link
     * ContactTable}'s bin-level ANY-OVERLAP bitset, the SAME 2000-tick run
     * (extent (4,4,4), 1920 slots, quanta in [-6,6], seed 42L) measured:
     * hybrid collision rate = 26.9/tick, effective ratio 17.5%; LGA
     * collision rate = 178.1/tick, effective ratio 68.6% -- roughly 6.6x
     * and 3.9x respectively, not "within CI" at any tolerance that would
     * still mean something. Root-caused (not an LGA implementation bug):
     * {@link ContactTable#of(ContactAtlas)} correctly consumed the atlas's
     * ANY-OVERLAP {@code overlapFraction > 0} signal (bead inviscid-gyt,
     * deliberately widened from bin-centre-only {@code contact} to fix a
     * ~12% true-contact reproduction rate) -- a bin fires the moment ANY
     * fraction of its fine geometric sub-sweep contacts, a strict superset
     * of {@code ContactPredicate}'s live point evaluation. This finding
     * was escalated (not silently toleranced) and led directly to the
     * USER DECISION (2026-08-08, FINAL) that produced {@link
     * FineStepContactTable}: contact firing moved to the atlas's fine
     * 360-step geometry grid, replicating point-evaluation semantics.
     *
     * <p><b>CURRENT (fine-step firing) -- GREEN, legitimately.</b> Same
     * 2000-tick run, same seed: hybrid collision rate = 26.9/tick,
     * effective ratio 17.5% (unchanged, as expected -- the hybrid itself
     * was never touched); LGA collision rate = 23.9/tick, effective ratio
     * 18.6% -- rate within ~11% (was 6.6x over), effective-ratio gap 0.011
     * (was 0.51, ~46x tighter). No tolerance was widened to reach this;
     * {@link #TOLERANCE} (0.15) is unchanged from the original, pre-fix
     * design.
     *
     * <p><b>Fix-round S1 correction (bead inviscid-0nx.23).</b> The
     * .23-gate critic found that {@code PhaseCMeasurement}'s
     * {@code COLLISION_FIELD totalCollisions} row cited THIS test for a
     * claim it did not actually assert -- the prose above documents the
     * ~11% rate gap, but until this fix-round the only executable
     * assertion below was on the EFFECTIVE-RATIO gap, not the totals
     * themselves. The assertion on {@code relativeRateGap} directly below
     * closes that gap: an LGA regression that doubled the collision rate
     * at a constant effective ratio would have passed this test before
     * the fix, and now does not.
     */
    @Test
    public void aggregateStatisticsAgreeBeyondDivergence() {
        SharedInitialCondition ic = sharedInitialCondition(42L, 6);
        CollisionStatistics hybridStats = new CollisionStatistics();
        CollisionStatistics lgaStats = new CollisionStatistics();
        HybridAutomaton hybrid = newHybrid(ic, hybridStats);
        LatticeGasAutomaton lga = newLga(ic, lgaStats);

        int ticks = 2000;
        for (int t = 0; t < ticks; t++) {
            hybrid.tick(t);
            lga.tick(t);
        }

        assertTrue("expected the hybrid to have collided over " + ticks
                   + " ticks", hybridStats.totalCollisions() > 0);
        assertTrue("expected the LGA to have collided over " + ticks
                   + " ticks", lgaStats.totalCollisions() > 0);

        double hybridRatio = (double) hybridStats.effectiveCollisions()
                              / hybridStats.totalCollisions();
        double lgaRatio = (double) lgaStats.effectiveCollisions()
                           / lgaStats.totalCollisions();

        assertTrue("expected effective/total collision ratios to agree within "
                   + TOLERANCE + " (hybrid=" + hybridRatio + " ["
                   + hybridStats.totalCollisions() + " total/" + ticks
                   + " ticks], lga=" + lgaRatio + " ["
                   + lgaStats.totalCollisions() + " total/" + ticks
                   + " ticks]) -- same PHYSICS, not the same trajectory "
                   + "(test 7 already shows trajectories diverge)",
                   Math.abs(hybridRatio - lgaRatio) < TOLERANCE);

        // Fix-round S1: the total-collision RATE gap itself, not merely
        // the effective/total ratio -- the claim PhaseCMeasurement's
        // COLLISION_FIELD row actually needs verified. ticks is the same
        // for both substrates, so the rate ratio is exactly proportional
        // to the totals ratio.
        double hybridRate = (double) hybridStats.totalCollisions() / ticks;
        double lgaRate = (double) lgaStats.totalCollisions() / ticks;
        double relativeRateGap = Math.abs(hybridRate - lgaRate) / hybridRate;
        assertTrue("expected total-collision RATE gap to agree within " + TOLERANCE
                   + " (hybridRate=" + hybridRate + " [" + hybridStats.totalCollisions()
                   + " total/" + ticks + " ticks], lgaRate=" + lgaRate + " ["
                   + lgaStats.totalCollisions() + " total/" + ticks
                   + " ticks], relativeRateGap=" + relativeRateGap
                   + ") -- this is the assertion PhaseCMeasurement's COLLISION_FIELD "
                   + "totalCollisions row cites (fix-round S1)",
                   relativeRateGap < TOLERANCE);
    }
}
