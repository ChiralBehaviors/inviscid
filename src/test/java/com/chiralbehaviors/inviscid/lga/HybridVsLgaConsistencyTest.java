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

import static com.chiralbehaviors.inviscid.Constants.TWO_PI;
import static org.junit.Assert.assertEquals;
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
 * quantisation error stays below one bin" -- using {@link #binOf(float,
 * int)} at the {@code phaseResolutionNLga} the atlas itself was
 * transcribed with ({@link ContactAtlas.Header}), not the finer 3600
 * fine-phase resolution. This was empirically necessary, not a stylistic
 * choice: an earlier draft of this test compared fine phase directly
 * (exact 3600-resolution equality) and found "divergence" at tick 1 for
 * every seed tried -- not a real physics signal, but {@code
 * QUANTUM_RATE == (float)(2*pi/3600)}'s float32 rounding occasionally
 * landing a hybrid member's angle within noise of a fine-phase unit
 * boundary, which is essentially guaranteed to happen SOMEWHERE among
 * 1920 slots within one tick and tells you nothing about contact-table
 * fidelity. A bin is 150 fine-phase units wide, so bin-level comparison
 * suppresses most of that noise and brings the signal this test exists to
 * measure -- whether the LGA's table-driven collision decisions track the
 * hybrid's live geometric ones -- above it.
 *
 * <p><b>It suppresses most of that noise, not all of it (measured, bead
 * inviscid-0nx.27).</b> An earlier version of this paragraph asserted that
 * float32 rounding noise (~1e-4 fine-phase units) is "~6 orders of
 * magnitude too small to flip a bin boundary on its own". That is false,
 * and it matters for reading {@link #firstDivergenceIsPinnedForSeed42}.
 * Measured on the seed-42 run: after one tick, 5 of the 1920 slots
 * disagree at BIN level, and every one of the 5 sits between 2.6e-5 and
 * 1.1e-4 fine-phase units BELOW a fine phase the LGA holds EXACTLY on a
 * bin boundary (300, 900, 1800, 3150, 3300), so each lands one bin low.
 * Noise that small does flip a bin, whenever the angle lands that close to
 * a boundary. What bin-level comparison buys is the SIZE of the residue,
 * not its absence: a handful of boundary-straddling slots (5 at tick 1;
 * between 1 and 12 across all 64 ticks) rather than the
 * essentially-any-slot noise of exact fine-phase comparison. And at no
 * tick in the 64-tick window does any slot's fine-phase separation exceed
 * one bin -- measured: no tick's maximum exceeds 150, and at tick 64 the
 * largest across slots is 34.
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
 * <p><b>Scope of that mechanism (measured, bead inviscid-0nx.27).</b> It
 * describes what is expected to drive genuine TRAJECTORY divergence. It is
 * NOT the mechanism behind the tick-1 disagreement pinned by {@link
 * #firstDivergenceIsPinnedForSeed42}, and reading it as such is a mistake
 * this Javadoc previously invited. Measured on the seed-42 run: the two
 * substrates' fine phases stay within 5.0e-4 fine-phase units through tick
 * 3, first exceed 1 fine-phase unit at tick 4 (2 slots), and by tick 64
 * the largest separation across slots is 34 fine-phase units. The bin
 * disagreement pinned at tick 1 therefore predates any measurable
 * trajectory separation between the substrates. Whether the any-overlap
 * transcription is what produces the tick-4-onward separation was NOT
 * measured on this bead; only the separation itself was.
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

    /**
     * The hybrid side's continuous-angle -> {@code N_lga} bin conversion.
     *
     * <p><b>Provenance (bead inviscid-0nx.27, E.0).</b> This is a
     * deliberately LITERAL transcription of the {@code bin(float)} method
     * of the phase-quantizer class that bead retired (the original is at
     * git {@code 3af3a47}). Compared against that revision line by line:
     * the same static import of {@code Constants.TWO_PI} (float, never
     * {@code 2 * Math.PI}), the same {@code TWO_PI / nLga} bin width (a
     * {@code private final float} field there, a {@code float} local
     * here), the same {@code % TWO_PI} normalisation and negative
     * fix-up, and the same top-of-circle CLAMP - the same float
     * operations in the same order.
     *
     * <p>An exhaustive differential sweep over every one of the ~2.17e9
     * float bit patterns in {@code [-2*pi, +2*pi]} (the entire domain
     * reachable after {@code % TWO_PI} normalisation) found 0 mismatches
     * between the two while both existed. Given the textual identity
     * above, <b>that sweep confirmed an identity the source already
     * shows, and was never load-bearing</b> - the same float operations
     * on the same input cannot disagree. It is recorded in T2 as history,
     * not as evidence this method needs, and not (as an earlier version
     * of this Javadoc framed it) as a regrettably unrepeatable
     * measurement. What actually defends this method going forward is
     * {@link #binOfClampsTopOfCircleRatherThanWrapping}, which IS
     * runnable.
     *
     * <p>Do not "simplify" this to {@code % nLga}, to {@code double}
     * arithmetic, or to a {@code stepOf}-routed rederivation. All three
     * change results, and {@link
     * #binOfClampsTopOfCircleRatherThanWrapping} asserts against each one
     * by name - see its Javadoc for the measured cases. Note in particular
     * that a {@code double} rewrite is NOT caught by the top-of-circle
     * cases (it reproduces them exactly); case 3 of that test is what
     * covers it, and deleting case 3 silently reopens that axis. Case 3's
     * coverage has a measured limit of its own - see its Javadoc for
     * exactly which {@code double} rewrites it does and does not
     * separate.
     *
     * @param angle
     *            continuous rotation angle, radians, any sign or magnitude
     * @param nLga
     *            the contact table's own bin count, sourced from the atlas
     *            header ({@code phaseResolutionNLga}) - never hardcoded, and
     *            never {@code Necronomata.PHASE_RESOLUTION}
     * @return the phase bin (0..{@code nLga}-1) {@code angle} falls into
     */
    static int binOf(float angle, int nLga) {
        float binWidth = TWO_PI / nLga;
        float normalized = angle % TWO_PI;
        if (normalized < 0) {
            normalized += TWO_PI;
        }
        int idx = (int) (normalized / binWidth);
        return idx >= nLga ? nLga - 1 : idx;
    }

    private static int hybridBin(Necronomata automaton, int nLga, int slot) {
        float[][] box = new float[1][];
        automaton.process((angle, frequency, deltaA, deltaF) -> box[0] = angle);
        return binOf(box[0][slot], nLga);
    }

    private static int lgaBin(LatticeGasAutomaton lga, int slot) {
        int[][] box = new int[1][];
        lga.process((phase, quanta) -> box[0] = phase);
        return box[0][slot] / lga.subBinSteps();
    }

    /**
     * The first point at which the two substrates' N_lga bins disagree.
     *
     * @param tick
     *            1-indexed count of ticks RUN when the disagreement was
     *            found, or {@code maxTicks} if they agreed throughout
     * @param slot
     *            the member slot that disagreed, or -1 if none did
     */
    private record Divergence(int tick, int slot) {
    }

    /**
     * @return the first tick (1-indexed count of ticks RUN) at which any
     *         member's N_lga BIN disagrees between the two substrates -
     *         see class Javadoc, "Comparison resolution" - together with
     *         the slot that disagreed, or {@code (maxTicks, -1)} if they
     *         agreed for the entire window (never expected to be reached
     *         here - see test assertions).
     */
    private static Divergence firstDivergence(HybridAutomaton hybrid,
                                                LatticeGasAutomaton lga,
                                                int nLga, int maxTicks) {
        int length = lga.slotCount();
        for (int tick = 0; tick < maxTicks; tick++) {
            hybrid.tick(tick);
            lga.tick(tick);
            for (int slot = 0; slot < length; slot++) {
                int hybridBin = hybridBin(hybrid.automaton(), nLga, slot);
                int lgaBin = lgaBin(lga, slot);
                if (hybridBin != lgaBin) {
                    return new Divergence(tick + 1, slot);
                }
            }
        }
        return new Divergence(maxTicks, -1);
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
     * {@link #firstDivergence}'s own control flow starts at {@code
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

        int maxTicks = 64;
        int divergenceTick = firstDivergence(hybrid, lga,
                                              ATLAS.header()
                                                   .phaseResolutionNLga(),
                                              maxTicks).tick();

        assertTrue("expected divergence to actually occur within the "
                   + maxTicks
                   + "-tick window (a threshold that never fires is vacuous) - "
                   + "observed divergence at tick " + divergenceTick,
                   divergenceTick < maxTicks);
    }

    /**
     * FIRST-BIN-DISAGREEMENT PIN (bead inviscid-0nx.27, E.0; renamed from
     * "substrate-agreement characterization pin" in that bead's second fix
     * round, because the measurements below do not support the stronger
     * name). The exact first point - tick AND slot - at which the hybrid
     * and LGA substrates' N_lga bins part company on the seed-42 /
     * 64-tick run. It pins the slot as well as the tick because the tick
     * alone is a weak pin: tick 1 is the minimum {@link #firstDivergence}
     * can return, and the class Javadoc's "Comparison resolution" section
     * has the measured reason a bin-level disagreement at tick 1 is
     * expected (5 slots straddle a bin boundary there).
     *
     * <p><b>This is NOT a guard for the binning swap that bead performed,
     * and must not be described as one.</b> Measured, not assumed: the
     * pinned {@code (tick 1, slot 161)} is INVARIANT under nine rejected
     * alternatives to {@link #binOf(float, int)} - a {@code % nLga} wrap;
     * four {@code double} variants (true-{@code double} divisor, float
     * divisor widened, fully-{@code double} normalisation, and {@code 2 *
     * Math.PI}); and four {@code stepOf}-routed rederivations. Every one
     * leaves this test GREEN. The binning semantics are guarded by {@link
     * #binOfClampsTopOfCircleRatherThanWrapping} and by that test ALONE;
     * if it is ever weakened, nothing else in the suite notices.
     *
     * <p><b>The invariance is EMPIRICAL; its mechanism is NOT
     * established.</b> An earlier version of this Javadoc explained it as
     * a full-bin substrate disagreement that swamps the boundary-level
     * effects those variants change. That explanation was measured and is
     * false. At tick 1 the two substrates agree: their fine phases match
     * on all 1920 slots to within 2.6e-4 fine-phase units, and slot 161 is
     * hybrid fine phase 899.99996 against LGA fine phase exactly 900 - a
     * gap of 4.3e-5 against a bin 150 fine-phase units wide, far SMALLER
     * than a bin rather than larger. What this test pins is a float
     * boundary artifact of the comparison INSTRUMENT: the hybrid's float
     * angle {@code 1.5707963f} gives quotient {@code 5.9999995} and hence
     * bin 5, while the LGA's exact integer {@code 900 / 150} gives bin 6.
     * Nor is the invariance universal - a two-stage FLOAT regroup (angle
     * -> 360 steps -> bin) returns 6 at slot 161, agrees with the LGA
     * there, and MOVES the pin to slot 758. So the invariance is a fact
     * about which specific angles those nine variants happen to move, and
     * no more than that has been measured.
     *
     * <p>What this pin IS for: it records where this whole configuration
     * - the two substrates, the atlas, the shared initial condition, AND
     * the {@link #binOf} instrument - first disagrees. Do not read it as a
     * characterization of the substrate pair alone: measured above, the
     * tick-1 disagreement is instrument-dominated, and at tick 1 the
     * substrates have not separated above float noise - 2.6e-4 fine-phase
     * units at maximum, measured above (the class Javadoc's
     * "Scope of that mechanism" paragraph has the by-tick numbers). A
     * later bead that changes the substrates and expects THIS pin to move
     * as its evidence is relying on something not established here. These
     * literals are OBSERVATIONS, not derived quantities - they
     * legitimately change whenever the substrates, the atlas, the shared
     * initial condition, or the binning instrument change. A failure here
     * means "something moved - go find out what", not "the code is wrong".
     */
    @Test
    public void firstDivergenceIsPinnedForSeed42() {
        SharedInitialCondition ic = sharedInitialCondition(42L, 6);
        HybridAutomaton hybrid = newHybrid(ic, new CollisionStatistics());
        LatticeGasAutomaton lga = newLga(ic, new CollisionStatistics());

        Divergence divergence = firstDivergence(hybrid, lga,
                                                  ATLAS.header()
                                                       .phaseResolutionNLga(),
                                                  64);

        assertEquals("first divergence TICK moved - the binning path or a "
                     + "substrate changed", 1, divergence.tick());
        assertEquals("first divergence SLOT moved - the binning path or a "
                     + "substrate changed", 161, divergence.slot());
    }

    /**
     * MIGRATED REGRESSION GUARD (bead inviscid-0nx.27, E.0). This test
     * inherits the duty of the deleted phase-quantizer test's {@code
     * topOfCircleMapsToLastBinNotFirst}: {@link #binOf(float, int)} must
     * CLAMP a top-of-circle division overshoot to the last bin, never WRAP
     * it to bin 0. Retiring that class removed the only executable
     * assertion of this semantics, so it moves here rather than dying with
     * the class.
     *
     * <p><b>Why CLAMP is the physically right answer here, not merely the
     * inherited one.</b> The provenance argument ("the retired class did
     * it this way") is not a justification and would leave this guard
     * rationale-less once that class is out of living memory. The real
     * reason is the shape of the thing being compared against: the LGA
     * side derives its bin as {@code phase[slot] / subBinSteps} from a
     * WRAPPED INTEGER accumulator confined to {@code [0, phaseResolution)}
     * - it can never produce an out-of-range index, and an angle sitting
     * just below {@code 2*pi} is genuinely its LAST bin. A wrapping binner
     * on the hybrid side would map that same angle to bin 0 and hand
     * {@link #firstDivergence} a FULL-BIN disagreement that exists purely
     * as an artifact of the comparison instrument. Clamping is what makes
     * the hybrid side's index range match the LGA side's by construction;
     * that is what makes this guard a statement about the instrument
     * rather than a test agreeing with itself.
     *
     * <p><b>All three rejected alternatives are asserted against
     * explicitly</b>, so this guard is non-vacuous - it names the wrong
     * answers, not just the right one. Measured on this bead (values below
     * are empirical, not assumed):
     * <ul>
     * <li><b>Case 1, float-division overshoot near {@code 2*pi}.</b> At
     * {@code nLga=360}, {@code nextDown(TWO_PI)/binWidth} rounds up to
     * exactly {@code 360.0f}: clamp gives 359, a {@code % nLga} wrap gives
     * 0. This is the trigger the inviscid-0nx.18 code review found. It does
     * NOT fire at production {@code nLga=24} (measured: clamp and wrap both
     * give 23), which is precisely why the guard cannot be written at 24
     * alone.</li>
     * <li><b>Case 2, tiny-negative-angle absorption.</b> For {@code angle =
     * -eps} with {@code eps} below about half a ULP of {@code TWO_PI}
     * (~2.4e-7), the float add {@code normalized += TWO_PI} rounds back to
     * exactly {@code TWO_PI}, so the division yields exactly {@code nLga}.
     * This one DOES fire at production {@code nLga=24}: clamp gives 23,
     * while BOTH a {@code % nLga} wrap AND a {@code
     * ContactAtlasGenerator.stepOf}-routed rederivation give 0. Measured
     * boundary: {@code -1e-7} is absorbed, {@code -1e-6} is not. The
     * stepOf-routed 0 is asserted below rather than merely claimed here
     * (it was prose-only until the second fix round of this bead);
     * {@code MemberGeometry.stepOf}'s body gives the same 0, measured,
     * but is not asserted because the two bodies are duplicates and one
     * assertion pins the shape.</li>
     * <li><b>Case 3, float-vs-double index division at a BIN BOUNDARY.</b>
     * The two cases above are top-of-circle float rounding artifacts, and
     * a {@code double}-arithmetic {@code binOf} that RETAINED the clamp
     * would pass both of them verbatim (measured: 359 and 23) - so they do
     * not cover the {@code double} axis at all. This case covers PART of
     * that axis, and the part was measured exactly. CAUGHT (each returns
     * 359 / 23 / <b>8</b>, so only this case separates it): a {@code
     * double} index division with a true-{@code double} divisor {@code
     * (double) TWO_PI / nLga}; one with the FLOAT divisor widened, {@code
     * (double) (TWO_PI / nLga)}; and a fully-{@code double} normalisation
     * that still uses the float {@code TWO_PI} constant. NOT CAUGHT -
     * measured, each returns 359 / 23 / <b>9</b> and passes all three
     * cases here: a {@code double} rewrite that ALSO swaps the constant to
     * {@code 2 * Math.PI}, live at 5 of the 122880 sampled angles (as is a
     * {@code stepOf}-shaped variant using the same constant); and a
     * two-stage FLOAT regroup (angle -> 360 steps -> bin, clamp
     * retained), live at <b>9</b> of the same 122880 - the highest
     * liveness of the three escapers. This case therefore guards the
     * WIDENING, not the CONSTANT, and the constant is not the only escape
     * route: the two-stage regroup differs on STAGE COUNT, changing
     * neither constant nor width. {@link
     * #binOf(float, int)}'s Javadoc names the constant in prose ("never
     * {@code 2 * Math.PI}") and nothing in this suite asserts it. At the
     * float nearest {@code 3*pi/4} ({@code 2.3561945f}) and production
     * {@code nLga=24}, the float quotient is EXACTLY {@code 9.0} while the
     * double quotient is {@code 8.999999772327}: float gives bin 9, double
     * gives bin 8. This class is not exotic and not top-of-circle - there
     * are 10 such angles across {@code [0, 2*pi]} at {@code nLga=24} (328
     * at {@code nLga=360}), they are bin boundaries including {@code
     * 3*pi/4} and {@code 3*pi/2}, and unlike cases 1 and 2 this class is
     * LIVE: the seed-42 / 64-tick run of {@link
     * #firstDivergenceIsPinnedForSeed42} produces float/double bin
     * disagreements on 10 of its 122880 sampled angles, {@code
     * 2.3561945f} among them. Which answer is "more correct" is beside the
     * point - float is what the atlas and every pinned artifact were
     * transcribed against, so {@code binOf} must reproduce the float
     * rounding, not improve on it.</li>
     * </ul>
     *
     * <p><b>Why {@code binOf} is a literal transcription and not a call to
     * the {@code stepOf} siblings.</b> {@code MemberGeometry.stepOf} and
     * {@code ContactAtlasGenerator.stepOf} WRAP where this CLAMPS, compute
     * the index in {@code double} where this uses {@code float}, and reach
     * the bin in two stages (angle -> 360 steps -> bin) rather than one.
     * Those three differences do not cancel: at {@code nLga=100} a
     * stepOf-routed rederivation disagrees on 40/360 bin centres (measured
     * on bead inviscid-0nx.18, whose substantive critique first reported
     * the number; reproduced on this bead), agreeing at today's {@code
     * nLga=24} only because {@code 24 | 360} (measured: 0/360 bin centres
     * at {@code nLga=24}). Even at {@code nLga=24} the substitution is
     * live rather than latent: <b>3</b> of the same 122880 sampled angles,
     * on 2 distinct angles ({@code 2.6179938} and {@code 4.9741883}), bin
     * differently under a stepOf-routed rederivation.
     *
     * <p><b>That 3 depends on one detail, named here because a
     * measurement on this very bead got it wrong.</b> Both real {@code
     * stepOf} bodies compute {@code angularResolution} as {@code TWO_PI /
     * resolution} - a FLOAT division that is only then widened to {@code
     * double} ({@code 0.017453292384743690}) - NOT a true {@code double}
     * division ({@code (double) TWO_PI / resolution} =
     * {@code 0.017453293005625408}). Scoring the true-{@code double}
     * divisor instead reports 83 on the same sample, ~28x too many, and
     * 83 is what this Javadoc claimed until the second fix round. The
     * qualitative conclusion is the same under either divisor: the count
     * is {@literal >} 0, so the axis is LIVE. See bead inviscid-ann for
     * the still-open question of whether the {@code stepOf} siblings
     * should adopt this clamp.
     */
    @Test
    public void binOfClampsTopOfCircleRatherThanWrapping() {
        float justUnderTwoPi = Math.nextDown(TWO_PI);
        assertEquals("nLga=360 top-of-circle must clamp to the LAST bin; 0 "
                     + "means a '% nLga' wrap was reintroduced", 359,
                     binOf(justUnderTwoPi, 360));

        float tinyNegative = -1e-7f;
        assertEquals("a tiny negative angle normalises to exactly TWO_PI in "
                     + "float and must clamp to the LAST bin; 0 means either "
                     + "a '% nLga' wrap or a stepOf-routed rederivation was "
                     + "substituted", 23, binOf(tinyNegative, 24));

        // Case 3: the FLOAT-vs-DOUBLE axis. Cases 1 and 2 are top-of-circle
        // artifacts that a double-arithmetic binOf reproduces exactly (it
        // returns 359 and 23 for them), so they leave that axis unguarded.
        // This angle - the float nearest 3*pi/4, and one that the seed-42
        // run actually reaches - separates them at production nLga=24.
        float threeQuarterPi = 2.3561945f;
        assertEquals("float-precision index division must be preserved; 8 "
                     + "means the division was widened to double", 9,
                     binOf(threeQuarterPi, 24));

        // Non-vacuity: each case must actually be capable of discriminating
        // the variant it exists to reject. If the arithmetic ever changed
        // such that these angles stopped landing on their boundaries, the
        // assertions above would still pass while guarding nothing.
        assertEquals("case 1 must sit exactly on the overshoot boundary", 360,
                     (int) (normalise(justUnderTwoPi) / (TWO_PI / 360)));
        assertEquals("case 2 must normalise to exactly TWO_PI", TWO_PI,
                     normalise(tinyNegative), 0.0f);
        assertEquals("case 3's float quotient must be exactly 9.0 - if it is "
                     + "not, this angle no longer separates float from double",
                     9.0f, normalise(threeQuarterPi) / (TWO_PI / 24), 0.0f);
        assertEquals("case 3's double quotient must floor BELOW the float "
                     + "one, or the double variant is not being rejected", 8,
                     (int) ((double) normalise(threeQuarterPi)
                            / ((double) TWO_PI / 24)));
        int geometryResolution = ATLAS.header().geometryResolution();
        assertEquals("case 2 must actually reject a stepOf-routed "
                     + "rederivation - 23 here means EITHER that "
                     + "ContactAtlasGenerator.stepOf stopped wrapping this "
                     + "angle to step 0, OR that the atlas geometryResolution "
                     + "moved off {180, 360, 720} (measured: 24, 40, 72, 120, "
                     + "1800 and 3600 all give 23, the quotient flooring "
                     + "below geometryResolution; 8 and 16 give 21 and 22, "
                     + "the most their step ranges can reach); either way "
                     + "case 2 "
                     + "no longer discriminates that substitution (see bead "
                     + "inviscid-ann)", 0,
                     ContactAtlasGenerator.binOfStep(ContactAtlasGenerator.stepOf(tinyNegative,
                                                                                   geometryResolution),
                                                      24, geometryResolution));
    }

    private static float normalise(float angle) {
        float normalized = angle % TWO_PI;
        return normalized < 0 ? normalized + TWO_PI : normalized;
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
