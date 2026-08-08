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
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.lga.CollisionSweep.AppliedCollision;
import com.chiralbehaviors.inviscid.lga.CollisionSweep.ReconciliationException;
import com.chiralbehaviors.inviscid.lga.CollisionSweep.TickResult;
import com.chiralbehaviors.inviscid.measure.CollisionStatistics;
import com.chiralbehaviors.inviscid.measure.ConservationAudit;
import com.chiralbehaviors.inviscid.measure.ConservationAudit.AuditResult;
import com.chiralbehaviors.inviscid.measure.ConservationAudit.LedgerEntry;
import com.chiralbehaviors.inviscid.measure.ConservationAudit.Violation;

/**
 * Behavioral tests for {@link CollisionSweep} (bead inviscid-0nx.14's rule
 * loop, and inviscid-1yk's direction-attribution seam).
 *
 * @author halhildebrand
 */
public class CollisionSweepTest {

    private static final double RADIUS     = 0.015;
    private static final int    RESOLUTION = 360;
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
     * A single cell pair, seeded so a contact fires (same fixture
     * {@code ContactScanTest} uses) and with unequal quanta so the rule
     * makes a real (non-no-op) decision.
     */
    private static final class Fixture {
        final Point3i          extent = new Point3i(4, 4, 4);
        final Necronomata      automaton = new Necronomata(extent);
        final FccNeighborhood  neighborhood = new FccNeighborhood(extent);
        final ContactScan      scan = new ContactScan(automaton, neighborhood,
                                                        newPredicate());
        final CollisionStatistics statistics = new CollisionStatistics();
        final CollisionRule    rule = new QuantaExchangeRule();

        Fixture() {
            Point3i cellA = new Point3i(0, 0, 0);
            Point3i cellB = neighborhood.neighbor(cellA, 1);
            seedMember(cellA, 3, 1, 1.7627826f, 5f);
            seedMember(cellB, 3, 0, 1.7627826f, 2f);
        }

        void seedMember(Point3i cell, int cube, int member, float angle,
                        float quanta) {
            int localIndex = cube * MEMBERS_PER_CUBE + member;
            automaton.process((angleArray, frequency, deltaA, deltaF) -> {
                int index = automaton.indexOfCell(cell) + localIndex;
                angleArray[index] = angle;
                frequency[index] = quanta;
            });
        }

        long quantaAt(Point3i cell, int cube, int member) {
            int localIndex = cube * MEMBERS_PER_CUBE + member;
            float[] value = new float[1];
            automaton.process((angleArray, frequency, deltaA, deltaF) -> value[0] = frequency[automaton.indexOfCell(cell)
                                                                                                + localIndex]);
            return Math.round((double) value[0]);
        }
    }

    @Test
    public void resolvesFixtureContactAndAppliesViaDeltaF() {
        Fixture fixture = new Fixture();
        CollisionSweep sweep = new CollisionSweep(fixture.automaton,
                                                    fixture.scan, fixture.rule,
                                                    fixture.statistics);

        TickResult result = sweep.tick(0);

        assertTrue("expected the seeded fixture contact to be found",
                   !result.applied().isEmpty());
        fixture.automaton.step();

        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = fixture.neighborhood.neighbor(cellA, 1);
        assertEquals("higher-quanta member A must have given up one quantum",
                     4L, fixture.quantaAt(cellA, 3, 1));
        assertEquals("lower-quanta member B must have gained one quantum",
                     3L, fixture.quantaAt(cellB, 3, 0));
    }

    @Test
    public void everyResolvedContactIsRecordedIncludingNoOps() {
        Fixture fixture = new Fixture();
        CollisionSweep sweep = new CollisionSweep(fixture.automaton,
                                                    fixture.scan, fixture.rule,
                                                    fixture.statistics);

        TickResult result = sweep.tick(0);

        assertEquals("every resolved contact must be recorded (including no-ops)",
                     result.applied().size(),
                     fixture.statistics.totalCollisions());
    }

    /**
     * Determinism of {@link CollisionSweep#tick(int)} itself (code-review
     * Important finding on the .14 stacked review: the prior version of
     * this test only re-scanned the same {@code ContactScan} twice and
     * never called {@code tick()} at all). Two INDEPENDENTLY-constructed,
     * identically-seeded fixtures must resolve to identical {@link
     * TickResult#applied()} lists and reach identical post-tick {@code
     * frequency} state - proving determinism end to end, not just that
     * {@code ContactScan} alone is deterministic (already covered by
     * {@code ContactScanTest.scanIsDeterministic}).
     */
    @Test
    public void tickIsDeterministicAcrossIndependentlySeededFixtures() {
        Fixture fixtureOne = new Fixture();
        Fixture fixtureTwo = new Fixture();
        CollisionSweep sweepOne = new CollisionSweep(fixtureOne.automaton,
                                                       fixtureOne.scan,
                                                       fixtureOne.rule,
                                                       fixtureOne.statistics);
        CollisionSweep sweepTwo = new CollisionSweep(fixtureTwo.automaton,
                                                       fixtureTwo.scan,
                                                       fixtureTwo.rule,
                                                       fixtureTwo.statistics);

        TickResult resultOne = sweepOne.tick(0);
        TickResult resultTwo = sweepTwo.tick(0);
        assertEquals("identically-seeded, independently-constructed fixtures must resolve identical applied-collision lists",
                     resultOne.applied(), resultTwo.applied());

        fixtureOne.automaton.step();
        fixtureTwo.automaton.step();
        assertArrayEquals("identically-seeded, independently-constructed fixtures must reach identical post-tick frequency state",
                           frequencySnapshot(fixtureOne.automaton),
                           frequencySnapshot(fixtureTwo.automaton), 0f);
    }

    private static float[] frequencySnapshot(Necronomata automaton) {
        float[][] captured = new float[1][];
        automaton.process((angleArray, frequency, deltaA, deltaF) -> captured[0] = frequency.clone());
        return captured[0];
    }

    @Test
    public void signedTransferTotalIsAlwaysZero() {
        Fixture fixture = new Fixture();
        CollisionSweep sweep = new CollisionSweep(fixture.automaton,
                                                    fixture.scan, fixture.rule,
                                                    fixture.statistics);

        TickResult result = sweep.tick(0);
        assertEquals(0L, result.signedTransferTotal());
    }

    @Test
    public void reconcileWithLedgerPassesWhenLedgerDeltaIsZero() {
        Fixture fixture = new Fixture();
        CollisionSweep sweep = new CollisionSweep(fixture.automaton,
                                                    fixture.scan, fixture.rule,
                                                    fixture.statistics);

        TickResult result = sweep.tick(0);
        // Should not throw: a conservative rule's signed transfer total is
        // always zero, and a correctly-applied tick's real ledger delta is
        // also zero.
        CollisionSweep.reconcileWithLedger(result, 0L);
    }

    @Test
    public void reconcileWithLedgerThrowsOnMismatch() {
        Fixture fixture = new Fixture();
        CollisionSweep sweep = new CollisionSweep(fixture.automaton,
                                                    fixture.scan, fixture.rule,
                                                    fixture.statistics);

        TickResult result = sweep.tick(0);
        assertThrows(ReconciliationException.class,
                     () -> CollisionSweep.reconcileWithLedger(result, 7L));
    }

    /**
     * ce3 reconciliation, negative control: a {@link CollisionSweep}
     * subclass that lies about the magnitude it hands to {@link
     * CollisionStatistics#recordCollision} - the "deliberately
     * mis-recorded transfer" case - must be caught by {@link
     * CollisionSweep#tick(int)}'s recording-integrity cross-check.
     */
    @Test
    public void deliberatelyMisRecordedTransferIsCaught() {
        Fixture fixture = new Fixture();
        CollisionSweep lying = new CollisionSweep(fixture.automaton,
                                                    fixture.scan, fixture.rule,
                                                    fixture.statistics) {
            @Override
            protected long magnitudeToRecord(CollisionRule.Delta delta) {
                return Math.abs(delta.deltaA()) + 1;
            }
        };

        assertThrows(ReconciliationException.class, () -> lying.tick(0));
    }

    @Test
    public void directionsTouchingFindsTheAppliedDirection() {
        Fixture fixture = new Fixture();
        CollisionSweep sweep = new CollisionSweep(fixture.automaton,
                                                    fixture.scan, fixture.rule,
                                                    fixture.statistics);

        TickResult result = sweep.tick(0);
        Point3i cellA = new Point3i(0, 0, 0);

        List<Integer> directions = CollisionSweep.directionsTouching(result,
                                                                       cellA,
                                                                       3, 1);
        assertTrue("expected direction +1 to be attributed to the touched member",
                   directions.contains(1));
    }

    @Test
    public void directionsTouchingIsEmptyForAnUntouchedMember() {
        Fixture fixture = new Fixture();
        CollisionSweep sweep = new CollisionSweep(fixture.automaton,
                                                    fixture.scan, fixture.rule,
                                                    fixture.statistics);

        TickResult result = sweep.tick(0);
        Point3i farCell = new Point3i(2, 2, 2);

        List<Integer> directions = CollisionSweep.directionsTouching(result,
                                                                       farCell,
                                                                       0, 0);
        assertTrue(directions.isEmpty());
    }

    /**
     * Test 6: 1000 ticks on a seeded lattice, full loop (scan -> rule ->
     * step), {@code ConservationAudit} clean and EXACT every tick.
     * Angles and quanta are deliberately seeded away from the all-zero
     * rest pose (cube 3's rest-position contacts, propagated from bead
     * inviscid-0nx.13) so tick-0 collisions are an accounted-for
     * consequence of the seed, not a surprise.
     */
    @Test
    public void conservationAuditPassesOverALongRun() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 42L);
        seedRandomQuanta(automaton, extent, 42L, 6);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionRule rule = new QuantaExchangeRule();
        CollisionSweep sweep = new CollisionSweep(automaton, scan, rule,
                                                    statistics);
        ConservationAudit audit = new ConservationAudit(automaton);

        for (int tick = 0; tick < 1000; tick++) {
            TickResult result = sweep.tick(tick);
            automaton.step();
            AuditResult auditResult = audit.auditTick(tick);
            assertTrue("conservation violated at tick " + tick + ": "
                       + auditResult.violations(), auditResult.isClean());

            List<LedgerEntry> ledger = audit.ledger();
            LedgerEntry entry = ledger.get(ledger.size() - 1);
            CollisionSweep.reconcileWithLedger(result,
                                                entry.totalAfter()
                                                - entry.totalBefore());
        }

        // Non-vacuity (code-review Important 2 / critic Significant 2):
        // over 1000 ticks on a seeded, contact-bearing lattice, the loop
        // must actually have resolved contacts, AND at least some of them
        // must have been real (nonzero) transfers -- not merely a clean
        // audit over a run where nothing ever collided or every collision
        // happened to tie.
        assertTrue("expected at least one resolved contact over 1000 ticks",
                   statistics.totalCollisions() > 0);
        assertTrue("expected at least one non-no-op (effective) transfer over 1000 ticks",
                   statistics.effectiveCollisions() > 0);
    }

    /**
     * Same-tick, same-member multi-contact regression (bead inviscid-72s,
     * Critical finding on the .14 stacked review): every contact this
     * tick must resolve against the FROZEN pre-tick snapshot, never
     * against a same-tick partially-accumulated {@code deltaF} - matching
     * bead inviscid-0nx.15's explicit double-buffered mandate ("all
     * contacts are detected against the pre-tick state... a single-buffer
     * implementation makes the result depend on scan order") and Phase
     * C's synchronous collision-table update. A rejected "sequential"
     * alternative (resolve against {@code frequency[i] + deltaF[i]},
     * accumulated in scan order) was tried first and is exactly what this
     * test rules out - it makes the outcome depend on which contact is
     * processed first, which is empirically common (measured ~18.8% of
     * ticks have a same-member multi-contact in exploratory runs) and
     * risks masquerading as spurious anisotropy once B.5's directional
     * statistics are compared. The overdraw concern that might otherwise
     * motivate sequential resolution does not apply: quanta are signed
     * {@code long}s with no floor, so "overdrawing" a member below zero
     * is legal, not a bug to guard against.
     *
     * <p>Fixture: member A (quanta 5) is touched by two contacts this
     * tick - contact1 against B (quanta 5, a tie under the pre-tick
     * snapshot) and contact2 against C (quanta 4, A &gt; C under the
     * pre-tick snapshot). Under snapshot semantics both contacts read the
     * SAME pre-tick values regardless of resolution order, so contact1 is
     * always a no-op and contact2 always transfers exactly one quantum
     * from A to C - the per-contact deltas, not just the net lattice
     * state, are asserted identical for both processing orders.
     */
    @Test
    public void sameTickSameMemberMultiContactResolvesAgainstThePreTickSnapshot() {
        CollisionRule.Delta[] forward = resolveStubbedPair(false);
        CollisionRule.Delta[] reversed = resolveStubbedPair(true);

        CollisionRule.Delta expectedContact1 = CollisionRule.Delta.noop();
        CollisionRule.Delta expectedContact2 = new CollisionRule.Delta(-1L,
                                                                         1L);

        assertEquals("contact1 (A vs B, tied at 5/5 pre-tick) must be a no-op regardless of order",
                     expectedContact1, forward[0]);
        assertEquals("contact2 (A vs C, 5>4 pre-tick) must transfer one quantum regardless of order",
                     expectedContact2, forward[1]);
        assertEquals("reversing resolution order must not change contact1's outcome",
                     expectedContact1, reversed[0]);
        assertEquals("reversing resolution order must not change contact2's outcome",
                     expectedContact2, reversed[1]);
    }

    /**
     * Runs a fresh fixture (member A=5, B=5, C=4, all in one cell, distinct
     * cube/member slots) through {@link CollisionSweep#tick(int)} with a
     * stubbed {@link ContactScan} returning exactly two contacts - A-vs-B
     * and A-vs-C - in either forward or reversed order, and returns the
     * two contacts' resolved deltas indexed [contact1(A-vs-B),
     * contact2(A-vs-C)] regardless of which order they were actually
     * processed in (so the caller can compare forward vs reversed
     * apples-to-apples).
     */
    private static CollisionRule.Delta[] resolveStubbedPair(boolean reverseOrder) {
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

        TickResult result = sweep.tick(0);

        CollisionRule.Delta deltaContact1 = null;
        CollisionRule.Delta deltaContact2 = null;
        for (AppliedCollision applied : result.applied()) {
            if (applied.contact().equals(contact1)) {
                deltaContact1 = applied.delta();
            } else if (applied.contact().equals(contact2)) {
                deltaContact2 = applied.delta();
            }
        }
        return new CollisionRule.Delta[] { deltaContact1, deltaContact2 };
    }

    private static void seedQuanta(Necronomata automaton, Point3i cell,
                                    int cube, int member, float quanta) {
        int localIndex = cube * MEMBERS_PER_CUBE + member;
        automaton.process((angleArray, frequency, deltaA, deltaF) -> frequency[automaton.indexOfCell(cell)
                                                                                + localIndex] = quanta);
    }

    /**
     * A {@link ContactScan} test double that returns a fixed contact list
     * instead of running real geometry - used to force a deterministic
     * same-tick multi-contact scenario ({@link
     * #sameTickSameMemberMultiContactResolvesAgainstThePreTickSnapshot}).
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
     * Closes bead inviscid-1yk end to end: a REAL {@link
     * ConservationAudit.Violation} (not a hand-built one) is captured
     * from {@link ConservationAudit#auditTick(int)} and fed into {@link
     * CollisionSweep#directionsTouching(TickResult, Point3i, int, int)},
     * asserting the returned direction set is exactly the seeded
     * fixture's known contact direction.
     *
     * <p>Conservation-only auditing has a documented blind spot ({@code
     * ConservationAudit}'s own "Scope" Javadoc): a correctly-applied,
     * zero-sum collision alone never produces a {@code Violation} - the
     * whole point of the rule being zero-sum. To get a REAL violation
     * report at all, this test deliberately corrupts the SAME member the
     * seeded contact already touches, via a direct {@code frequency}
     * write through {@code process()} - bypassing the sanctioned {@code
     * deltaF} path on purpose, simulating an external conservation bug
     * unrelated to the collision rule itself.
     */
    @Test
    public void directionsTouchingMatchesARealConservationAuditViolation() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        // Verified non-contact baseline (ContactScanTest): the default
        // all-zero rest pose has spurious cube-3 contacts, so start from
        // the known-clean uniform angle instead.
        fillUniformAngle(automaton, extent, 1.0f);

        Point3i cellA = new Point3i(0, 0, 0);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        Point3i cellB = neighborhood.neighbor(cellA, 1);
        int fixtureCubeA = 3;
        int fixtureMemberA = 1;
        int fixtureCubeB = 3;
        int fixtureMemberB = 0;
        float fixtureAngle = 1.7627826f;

        seedAngleAndQuanta(automaton, cellA, fixtureCubeA, fixtureMemberA,
                            fixtureAngle, 5f);
        seedAngleAndQuanta(automaton, cellB, fixtureCubeB, fixtureMemberB,
                            fixtureAngle, 2f);

        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        ConservationAudit audit = new ConservationAudit(automaton);

        TickResult result = sweep.tick(0);
        assertEquals("expected exactly the one seeded contact",
                     1, result.applied().size());

        // Deliberate corruption: direct frequency write, bypassing deltaF.
        automaton.process((angleArray, frequency, deltaA, deltaF) -> {
            int index = automaton.indexOfCell(cellA)
                        + fixtureCubeA * MEMBERS_PER_CUBE + fixtureMemberA;
            frequency[index] = frequency[index] + 5f;
        });
        automaton.step();

        AuditResult auditResult = audit.auditTick(0);
        assertTrue("expected the deliberate corruption to be reported as a violation",
                   !auditResult.isClean());

        Violation violation = auditResult.violations()
                                          .stream()
                                          .filter(v -> v.cell().equals(cellA)
                                                       && v.cube() == fixtureCubeA
                                                       && v.member() == fixtureMemberA)
                                          .findFirst()
                                          .orElseThrow(() -> new AssertionError("expected a violation at the corrupted member; got: "
                                                                                 + auditResult.violations()));

        List<Integer> directions = CollisionSweep.directionsTouching(result,
                                                                       violation.cell(),
                                                                       violation.cube(),
                                                                       violation.member());
        assertEquals("expected the corrupted member's violation to be attributed to exactly the seeded contact's direction",
                     List.of(1), directions);
    }

    private static void fillUniformAngle(Necronomata automaton, Point3i extent,
                                          float angle) {
        int length = 30 * extent.x * extent.y * extent.z;
        automaton.process((angleArray, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < length; i++) {
                angleArray[i] = angle;
            }
        });
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
}
