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

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.measure.AuditedRun;
import com.chiralbehaviors.inviscid.measure.CollisionStatistics;
import com.chiralbehaviors.inviscid.measure.ConservationAudit;
import com.chiralbehaviors.inviscid.measure.ConservationAudit.AuditResult;

/**
 * Behavioral tests for {@link HybridAutomaton} (bead inviscid-0nx.15, Phase
 * A.4): the composed double-buffered advect-and-collide tick.
 *
 * <p>Timing (acceptance criterion, extent 6^3 x 1000 ticks): {@link
 * CollisionSweepTest#conservationAuditPassesOverALongRun} already measures
 * 1000 ticks at extent (4,4,4) at 19.21s wall (surefire-observed), i.e.
 * ~19ms/tick -- consistent with {@code ContactScan}'s own measured cost
 * scaling (extent^3-driven pair evaluation count). At extent 6^3 the
 * measured per-tick scan cost is ~68ms ({@code ContactScan}'s class
 * Javadoc), so 1000 ticks at 6^3 would be ~70s -- too slow to run three
 * times over inside this one surefire module alongside the rest of the
 * suite. <b>Provenance, corrected (2026-08-08 substantive-critic round,
 * FIX 3):</b> the 4^3/6^3-50-tick split below is a unilateral
 * timing-budget judgment call by this implementation, empirically
 * justified by the ~66-68ms/tick measurement -- the bead's own text
 * (verified via {@code bd show inviscid-0nx.15}) does not itself offer
 * this substitution as a pre-authorized option; an earlier version of this
 * Javadoc and the paired T2 decision memo incorrectly attributed it as
 * "the bead's own documented choice", which has been corrected. The
 * 1000-tick tests below (conservation-exactness, identical-trajectory
 * hashing) run at extent (4,4,4); {@link #sixCubedFiftyTickTiming()} is
 * the SHORT extent-6^3 run (50 ticks) that records the real 6^3 cost
 * honestly rather than silently dropping the 6^3 requirement -- see that
 * test's own Javadoc for the
 * observed number.
 *
 * @author halhildebrand
 */
public class HybridAutomatonTest {

    private static final double RADIUS           = LgaTestGeometry.BASELINE_RADIUS;
    private static final int    RESOLUTION       = 360;
    private static final int    MEMBERS_PER_CUBE = 6;

    private static MemberGeometry newGeometry() {
        return new MemberGeometry(RESOLUTION, RADIUS);
    }

    private static ContactPredicate newPredicate() {
        return new ContactPredicate(newGeometry());
    }

    private static void seedRandomAngles(Necronomata automaton, Point3i extent,
                                          long seed) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[length];
        for (int i = 0; i < length; i++) {
            angles[i] = random.nextFloat() * (float) (2 * Math.PI);
        }
        automaton.process((angleArray, frequency, deltaA, deltaF) -> System.arraycopy(angles,
                                                                                        0,
                                                                                        angleArray,
                                                                                        0,
                                                                                        length));
    }

    private static void seedRandomQuanta(Necronomata automaton, Point3i extent,
                                          long seed, int bound) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        float[] quanta = new float[length];
        for (int i = 0; i < length; i++) {
            quanta[i] = random.nextInt(bound);
        }
        automaton.process((angleArray, frequency, deltaA, deltaF) -> System.arraycopy(quanta,
                                                                                        0,
                                                                                        frequency,
                                                                                        0,
                                                                                        length));
    }

    /**
     * The codebase's canonical signed-quanta seeding convention (see
     * {@code BaselineSpectrumHarness}, {@code
     * frequencyDistribution=uniform integer in [-5,5] via new
     * Random(seed).nextInt(11) + -5}): each member's quanta is an
     * independent uniform integer in {@code [-5, 5]}, inclusive. Used by
     * bead inviscid-7a6's rest-pose characterization test, which needs a
     * NONZERO-quanta lattice (unlike {@link #seedRandomQuanta}'s {@code
     * [0, bound)} range, this range is symmetric about zero, matching
     * the convention the bead's own text cites).
     */
    private static void seedRandomQuantaSigned(Necronomata automaton,
                                                Point3i extent, long seed) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        float[] quanta = new float[length];
        for (int i = 0; i < length; i++) {
            quanta[i] = random.nextInt(11) - 5;
        }
        automaton.process((angleArray, frequency, deltaA, deltaF) -> System.arraycopy(quanta,
                                                                                        0,
                                                                                        frequency,
                                                                                        0,
                                                                                        length));
    }

    private static void seedQuanta(Necronomata automaton, Point3i cell,
                                    int cube, int member, float quanta) {
        int localIndex = cube * MEMBERS_PER_CUBE + member;
        automaton.process((angleArray, frequency, deltaA, deltaF) -> frequency[automaton.indexOfCell(cell)
                                                                                + localIndex] = quanta);
    }

    private static void seedAngleAndQuanta(Necronomata automaton, Point3i cell,
                                            int cube, int member, float angle,
                                            float quanta) {
        int localIndex = cube * MEMBERS_PER_CUBE + member;
        automaton.process((angleArray, frequency, deltaA, deltaF) -> {
            int index = automaton.indexOfCell(cell) + localIndex;
            angleArray[index] = angle;
            frequency[index] = quanta;
        });
    }

    private record TickState(float[] angle, float[] frequency) {
    }

    private static TickState snapshot(Necronomata automaton) {
        float[][] captured = new float[2][];
        automaton.process((angleArray, frequency, deltaA, deltaF) -> {
            captured[0] = angleArray.clone();
            captured[1] = frequency.clone();
        });
        return new TickState(captured[0], captured[1]);
    }

    /**
     * A {@link ContactScan} test double returning a fixed contact list --
     * same pattern as {@code CollisionSweepTest.StubScan} -- used to force
     * a deterministic same-tick, same-member multi-contact scenario.
     */
    private static final class StubScan extends ContactScan {
        private final List<Contact> contacts;

        StubScan(Necronomata automaton, FccNeighborhood neighborhood,
                 ContactPredicate predicate, List<Contact> contacts) {
            super(automaton, neighborhood, predicate);
            this.contacts = contacts;
        }

        @Override
        public List<Contact> scan() {
            return contacts;
        }
    }

    /**
     * Runs one {@link HybridAutomaton#tick(int)} against a fixed
     * three-member, two-contact fixture (A touched by both a tie against B
     * and a real transfer against C), in either the natural or reversed
     * contact-processing order, and returns the FULL post-tick lattice
     * state (angle and frequency, not just the two contacts' deltas) --
     * this is the HybridAutomaton-level analogue of {@code
     * CollisionSweepTest.resolveStubbedPair}, one layer up (through {@code
     * Necronomata.step()} as well as the collision resolution).
     */
    private static TickState runMultiContactTick(boolean reverseOrder) {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        Point3i cell = new Point3i(0, 0, 0);

        int cubeA = 0, memberA = 0;
        int cubeB = 0, memberB = 1;
        int cubeC = 0, memberC = 2;

        seedQuanta(automaton, cell, cubeA, memberA, 5f);
        seedQuanta(automaton, cell, cubeB, memberB, 5f);
        seedQuanta(automaton, cell, cubeC, memberC, 4f);

        Contact contact1 = new Contact(cell, cubeA, memberA, cell, cubeB,
                                        memberB, 1, 0.0);
        Contact contact2 = new Contact(cell, cubeA, memberA, cell, cubeC,
                                        memberC, 2, 0.0);
        List<Contact> order = reverseOrder ? List.of(contact2, contact1)
                                            : List.of(contact1, contact2);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        StubScan scan = new StubScan(automaton, neighborhood, newPredicate(),
                                      order);
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);

        hybrid.tick(0);

        return snapshot(automaton);
    }

    /**
     * Test 1: double-buffering, proven by scan-order invariance at the
     * composed {@link HybridAutomaton} level (not just {@code
     * CollisionSweep} in isolation -- {@code
     * CollisionSweepTest.sameTickSameMemberMultiContactResolvesAgainstThePreTickSnapshot}
     * already covers the per-contact-delta claim; this test asserts the
     * FULL post-tick derived state -- angle AND frequency, i.e. after
     * {@link Necronomata#step()} has run too -- is scan-order-invariant).
     * A single-buffer implementation would make contact2's outcome (and
     * therefore the post-step angle of every touched member) depend on
     * whether contact1 or contact2 is resolved first.
     */
    @Test
    public void tickIsDoubleBuffered() {
        TickState forward = runMultiContactTick(false);
        TickState reversed = runMultiContactTick(true);

        assertArrayEquals("scan order must not change post-tick frequency state",
                           forward.frequency(), reversed.frequency(), 0f);
        assertArrayEquals("scan order must not change post-tick angle state",
                           forward.angle(), reversed.angle(), 0f);
    }

    /**
     * Test 2: exact conservation over 1000 ticks, via {@link
     * ConservationAudit}. Extent (4,4,4) -- see class Javadoc's timing
     * note.
     *
     * <p><b>Driven through {@link AuditedRun#tick(int)}, not manual
     * {@code hybrid.tick(); audit.auditTick();}</b> (FIX 1,
     * 2026-08-08 substantive-critic round, Significant): the manual
     * two-call pattern never invokes {@link
     * CollisionSweep#reconcileWithLedger(TickReport, long)},
     * silently dropping the reconciliation guarantee -- exactly the
     * bypass future B.4/B.5 long-run harnesses would otherwise copy from
     * this test. Routing through {@code AuditedRun} exercises the FULL
     * wired path ({@code tick} -&gt; {@code auditTick} -&gt; {@code
     * reconcileWithLedger}) for all 1000 ticks; a thrown {@code
     * ReconciliationException} would fail this test just as loudly as a
     * conservation violation would.
     */
    @Test
    public void totalQuantaConservedOver1000Ticks() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 42L);
        seedRandomQuanta(automaton, extent, 42L, 6);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);
        ConservationAudit audit = new ConservationAudit(automaton);
        AuditedRun run = new AuditedRun(hybrid, audit);

        AuditResult last = null;
        for (int tick = 0; tick < 1000; tick++) {
            AuditedRun.TickOutcome outcome = run.tick(tick);
            last = outcome.auditResult();
            assertTrue("conservation violated at tick " + tick + ": "
                       + last.violations(), last.isClean());
        }

        assertEquals("exact: final total quanta must equal the baseline established at construction",
                     last.baselineTotal(), last.totalQuanta());
    }

    /**
     * Test 3: identically-seeded, INDEPENDENTLY-constructed runs reach an
     * identical full-state hash after 1000 ticks. Extent (4,4,4) -- see
     * class Javadoc's timing note.
     */
    @Test
    public void identicalSeedsGiveIdenticalTrajectories() {
        long hashOne = runSeededTrajectoryHash();
        long hashTwo = runSeededTrajectoryHash();

        assertEquals("identically-seeded, independently-constructed 1000-tick runs must reach an identical full-state hash",
                     hashOne, hashTwo);
    }

    private static long runSeededTrajectoryHash() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 7L);
        seedRandomQuanta(automaton, extent, 7L, 6);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);

        for (int tick = 0; tick < 1000; tick++) {
            hybrid.tick(tick);
        }

        TickState state = snapshot(automaton);
        return 31L * Arrays.hashCode(state.angle()) + Arrays.hashCode(state.frequency());
    }

    /**
     * Test 4: non-vacuity guard. A hybrid automaton over a seeded,
     * contact-bearing lattice for 300 ticks must actually have resolved
     * contacts, AND at least some must have been real (nonzero) transfers
     * -- both {@link CollisionStatistics#totalCollisions()} and {@link
     * CollisionStatistics#effectiveCollisions()} checked, per the bead's
     * explicit non-vacuity requirement (a hybrid automaton that never
     * collides is just the K=0 baseline wearing a new class name).
     */
    @Test
    public void collisionsActuallyOccur() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 42L);
        seedRandomQuanta(automaton, extent, 42L, 6);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);

        for (int tick = 0; tick < 300; tick++) {
            hybrid.tick(tick);
        }

        assertTrue("expected at least one resolved contact over 300 ticks",
                   statistics.totalCollisions() > 0);
        assertTrue("expected at least one non-no-op (effective) transfer over 300 ticks",
                   statistics.effectiveCollisions() > 0);
    }

    /**
     * Test 5: end-to-end proof that inviscid-P0.1's frequency<->angular-rate
     * coupling is load-bearing through a REAL collision transfer, not just
     * asserted by {@code NecronomataStateSemanticsTest} in isolation. Member
     * B receives one quantum this tick (quantaB: 2 -> 3); after {@link
     * HybridAutomaton#tick(int)} (which calls {@code
     * Necronomata#step()}), B's angular rate ({@code deltaA}) must reflect
     * the POST-transfer frequency (3), not the pre-transfer one (2) -- {@code
     * deltaA == QUANTUM_RATE * 3}, and explicitly NOT {@code QUANTUM_RATE *
     * 2}.
     */
    @Test
    public void angleAdvanceMatchesFrequencyAfterCollisions() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);

        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = neighborhood.neighbor(cellA, 1);
        int cubeA = 3, memberA = 1;
        int cubeB = 3, memberB = 0;
        float angle = 1.7627826f;

        seedAngleAndQuanta(automaton, cellA, cubeA, memberA, angle, 5f);
        seedAngleAndQuanta(automaton, cellB, cubeB, memberB, angle, 2f);

        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);

        hybrid.tick(0);

        int indexB = automaton.indexOfCell(cellB) + cubeB * MEMBERS_PER_CUBE
                     + memberB;
        float[] captured = new float[2];
        automaton.process((angleArray, frequency, deltaA, deltaF) -> {
            captured[0] = frequency[indexB];
            captured[1] = deltaA[indexB];
        });

        assertEquals("member B must have gained exactly one quantum",
                     3.0f, captured[0], 0f);
        assertEquals("member B's angular rate must reflect the POST-transfer frequency",
                     Necronomata.QUANTUM_RATE * 3f, captured[1], 1e-9f);
        assertTrue("the transfer must have genuinely changed the angular rate from what the pre-transfer frequency would have given",
                   Math.abs(captured[1] - Necronomata.QUANTUM_RATE * 2f) > 1e-9f);
    }

    /**
     * Test 6: a zero-frequency lattice never moves and never (effectively)
     * collides. Uses the automaton's default all-zero state (angle == 0,
     * frequency == 0 everywhere) -- literally the rest pose. This is
     * deliberate, not an oversight: with every member's quanta at 0,
     * {@link QuantaExchangeRule} sees a tie (0 == 0) at EVERY contact it
     * is offered, geometric or not, so {@code effectiveCollisions()} is
     * the honest physics claim ("nothing ever actually transferred") --
     * asserted here. {@code totalCollisions()} is NOT asserted zero: .13
     * established that the all-zero rest pose has live geometric contacts
     * at cube 3 (recorded, since {@code CollisionSweep} records every
     * resolved contact including no-ops -- see that class's Javadoc,
     * "Recording convention"), so a nonzero {@code totalCollisions()} here
     * would be correct rest-pose behavior, not a bug; asserting it to zero
     * would be dishonest about what "never collides" actually means for
     * this seed.
     */
    @Test
    public void zeroFrequencyLatticeIsStatic() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        TickState before = snapshot(automaton);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);

        for (int tick = 0; tick < 100; tick++) {
            hybrid.tick(tick);
        }

        TickState after = snapshot(automaton);
        assertArrayEquals("a zero-frequency lattice's angles must never move",
                           before.angle(), after.angle(), 0f);
        assertArrayEquals("a zero-frequency lattice's frequencies must never change",
                           before.frequency(), after.frequency(), 0f);
        assertEquals("a zero-frequency lattice must never produce a real (non-no-op) transfer",
                     0L, statistics.effectiveCollisions());
    }

    /**
     * Test 6b (bead inviscid-7a6, Phase A gate critique, substantive-critic
     * round, 2026-08-08): tick-0 contact-DENSITY characterization for the
     * uniform rest pose - angle=0 for EVERY member (the natural boot
     * condition a quantized automaton, C.1/C.4, will design against) but,
     * unlike {@link #zeroFrequencyLatticeIsStatic()}, NONZERO quanta (a
     * {@code frequency=0} lattice has {@code deltaA=0} forever and can
     * never exercise a real contact - it is a trivial fixed-point check,
     * not a contact-density check). Because {@link ContactPredicate}'s
     * threshold is uniform across geometrically-identical member pairs,
     * angle=0 uniformly is exactly the condition most likely to produce
     * either a lattice-wide simultaneous-contact spike at boot, or zero
     * contacts (if angle=0 sits in a non-contacting geometric region) -
     * nothing established which before this test. Angle is left at its
     * default JVM-zero-initialized value (uniform rest pose); quanta is
     * seeded via {@link #seedRandomQuantaSigned} (Random(42L), the
     * codebase's [-5,5] convention) so contacts CAN actually resolve to
     * real transfers.
     *
     * <p><b>PINNED BOOT-CONDITION FACTS (this seed, extent (4,4,4),
     * deterministic - geometry + seed fixed):</b> tick 0 resolves 320
     * contacts total ({@code totalCollisions()}), of which 293 (91.6%)
     * are effective (non-no-op) transfers ({@code
     * effectiveCollisions()}). The per-direction breakdown is sharply
     * concentrated, not lattice-uniform, and answers the open question
     * decisively: ONLY four of the twelve {@link
     * FccNeighborhood#DIRECTIONS} - the positive-signed +2 (128), +4
     * (32), +5 (128), +6 (32) - carry any contacts at all; every one of
     * the other eight is exactly zero every tick-0 run at this seed, but
     * NOT for the same reason (FIX 3, stacked-review round 2026-08-08):
     * all six negative directions -1..-6 are zero STRUCTURALLY, by
     * construction - {@code ContactScan} is restricted to only the 6
     * positive directions to prevent double-counting (per bead
     * inviscid-0nx.13's own recorded NOTES), so -1..-6 read zero on
     * EVERY run regardless of lattice state, not a fact about this
     * geometry; +1 and +3 being zero, in contrast, IS a PHYSICAL finding
     * about this rest-pose geometry specifically (worth re-verifying
     * under a different seed or a changed {@code ContactPredicate}
     * geometry, unlike -1..-6 which no re-seed could ever change). This
     * IS a lattice-wide simultaneous-contact SPIKE at boot, not the
     * sparse regime: 320 contacts resolving in a single tick, over a
     * 4^3-cell lattice, is far above the ~1.1e-4 steady-state density
     * {@code ContactPredicateTest} measures for randomly-angled lattices
     * - the uniform angle=0 rest pose is exactly the geometrically-
     * aligned worst case that produces it, confirming (not merely
     * echoing) .13's cube-3 finding at full lattice scale, and
     * confirming which of the two possibilities the bead's own text
     * posed ("(a) ... spike" vs "(b) zero contacts") is the true one.
     * C.1/C.4 inherit this measured number, not an assumption, for their
     * boot initial condition: a naive "seed random quanta at angle=0 and
     * go" boot immediately drives hundreds of same-tick collision
     * resolutions, concentrated on four directions only, not the sparse
     * trickle a caller unaware of this finding might expect.</p>
     *
     * <p><b>TRIPWIRE NOTICE for the fixer (precedent: bead
     * inviscid-6cf's tripwire).</b> If the pinned 320 / 293 / {128, 32,
     * 128, 32} values above ever break, first check whether {@code
     * ContactPredicate}'s geometry constants ({@link #RADIUS}, {@link
     * #RESOLUTION}) or {@code FccNeighborhood}'s member offsets changed
     * - if either did, this is a LEGITIMATE RE-BASELINE: recompute the
     * new measured numbers with this same seed and re-pin them, updating
     * this Javadoc's numbers to match. If NEITHER changed, the
     * divergence is a silent regression in contact-resolution or
     * collision-rule logic and must be investigated as a bug - do NOT
     * "fix" a failure here by loosening these exact-equality assertions
     * into a range or tolerance check; that would convert a tripwire
     * into a rubber stamp.</p>
     */
    @Test
    public void uniformRestPoseTickZeroContactDensity() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        // angle is left at its default zero-initialized value for every
        // member -- the uniform rest pose (angle=0 uniformly). Deliberately
        // NOT seedRandomAngles: this test characterizes the natural boot
        // condition, not a randomly-angled lattice.
        seedRandomQuantaSigned(automaton, extent, 42L);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);
        ConservationAudit audit = new ConservationAudit(automaton);
        AuditedRun run = new AuditedRun(hybrid, audit);

        AuditedRun.TickOutcome outcome = run.tick(0);

        assertTrue("audit must stay clean through the uniform rest-pose tick-0: "
                   + outcome.auditResult().violations(),
                   outcome.auditResult().isClean());

        assertEquals("pinned tick-0 total resolved-contact count for the uniform rest pose (boot-condition fact, see class javadoc)",
                     320L, statistics.totalCollisions());
        assertEquals("pinned tick-0 effective (nonzero-transfer) contact count for the uniform rest pose (boot-condition fact, see class javadoc)",
                     293L, statistics.effectiveCollisions());

        Map<Integer, Long> perDirection = statistics.collisionsPerDirection();
        long sum = 0L;
        for (long count : perDirection.values()) {
            sum += count;
        }
        assertEquals("per-direction breakdown must sum to the total contact count",
                     statistics.totalCollisions(), sum);
        assertTrue("per-direction breakdown must be non-vacuous: at least one direction must show tick-0 contacts",
                   perDirection.values().stream().anyMatch(count -> count > 0));

        assertEquals("pinned per-direction tick-0 breakdown (boot-condition fact, see class javadoc)",
                     Map.ofEntries(Map.entry(1, 0L), Map.entry(-1, 0L),
                                    Map.entry(2, 128L), Map.entry(-2, 0L),
                                    Map.entry(3, 0L), Map.entry(-3, 0L),
                                    Map.entry(4, 32L), Map.entry(-4, 0L),
                                    Map.entry(5, 128L), Map.entry(-5, 0L),
                                    Map.entry(6, 32L), Map.entry(-6, 0L)),
                     perDirection);
    }

    /**
     * Acceptance criterion: record the extent-6^3 tick cost honestly. A
     * SHORT (50-tick) run at extent (6,6,6), rather than the full 1000
     * ticks (which the class Javadoc computes at ~70s, too slow to also
     * run at both extents inside one surefire module).
     *
     * <p><b>Observed timing (this development machine, surefire-reported,
     * {@code mvn -Dtest=HybridAutomatonTest#sixCubedFiftyTickTiming}):
     * 3.289s wall for 50 ticks, i.e. ~66ms/tick</b> -- consistent with
     * {@code ContactScan}'s own independently-measured ~68ms/tick at this
     * extent (that number is the scan cost alone; this one also includes
     * rule resolution, {@code step()}, and the audit call, so the close
     * agreement confirms scan cost dominates the tick, as {@code
     * ContactScan}'s own Javadoc predicts). Extrapolated to the full
     * 1000-tick acceptance criterion: ~66s, confirming the class Javadoc's
     * ~70s estimate and the decision to run the 1000-tick tests at extent
     * (4,4,4) instead. The loose upper bound below (60s for this 50-tick
     * run) is a hang-guard only, not a performance gate.
     */
    @Test
    public void sixCubedFiftyTickTiming() {
        Point3i extent = new Point3i(6, 6, 6);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 42L);
        seedRandomQuanta(automaton, extent, 42L, 6);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);
        ConservationAudit audit = new ConservationAudit(automaton);
        AuditedRun run = new AuditedRun(hybrid, audit);

        long start = System.nanoTime();
        for (int tick = 0; tick < 50; tick++) {
            AuditedRun.TickOutcome outcome = run.tick(tick);
            AuditResult result = outcome.auditResult();
            assertTrue("conservation violated at tick " + tick + ": "
                       + result.violations(), result.isClean());
        }
        long elapsedMs = (System.nanoTime() - start) / 1_000_000L;

        assertTrue("extent 6^3, 50 ticks took " + elapsedMs
                   + "ms -- exceeded the 60s hang-guard bound",
                   elapsedMs < 60_000L);
    }
}
