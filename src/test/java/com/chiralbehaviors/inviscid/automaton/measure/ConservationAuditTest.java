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

package com.chiralbehaviors.inviscid.automaton.measure;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.List;
import java.util.Map;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.automaton.Necronomata;
import com.chiralbehaviors.inviscid.automaton.measure.ConservationAudit.Violation;

/**
 * The acceptance instrument for the (not-yet-written) collision rules,
 * bead inviscid-0nx.14 / .20: built first so the audit cannot be fitted to
 * the rules it is meant to judge.
 *
 * @author halhildebrand
 */
public class ConservationAuditTest {

    /**
     * A 2x2x2 extent has exactly one even-parity cell: (0,0,0). That's
     * enough surface area for these tests (30 members: 5 cubes x 6
     * members) without dragging in wrap/neighbor concerns that belong to
     * FccNeighborhoodTest.
     */
    private Necronomata freshAutomaton() {
        return new Necronomata(new Point3i(2, 2, 2));
    }

    @Test
    public void detectsInjectedQuantum() {
        Necronomata automaton = freshAutomaton();
        ConservationAudit audit = new ConservationAudit(automaton);

        // Plant one extra quantum directly into member 3 of cell (0,0,0),
        // bypassing the deltaF collision-transfer path entirely -- exactly
        // the kind of silent injection the audit must catch.
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 3] = frequency[base + 3] + 1f;
        });

        ConservationAudit.AuditResult result = audit.auditTick(1);

        assertFalse("planted injection must be detected", result.isClean());
        List<Violation> violations = result.violations();
        assertFalse("violation report must not be empty", violations.isEmpty());
        Violation v = violations.get(0);
        assertEquals(new Point3i(0, 0, 0), v.cell());
        assertEquals(0, v.cube());
        assertEquals(3, v.member());
        assertEquals(1, v.tick());
    }

    @Test
    public void detectsDeletedQuantum() {
        Necronomata automaton = freshAutomaton();

        // Seed a quantum first so there is something to delete.
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 7] = 2f;
        });

        ConservationAudit audit = new ConservationAudit(automaton);

        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 7] = frequency[base + 7] - 1f;
        });

        ConservationAudit.AuditResult result = audit.auditTick(1);

        assertFalse("planted deletion must be detected", result.isClean());
        Violation v = result.violations().get(0);
        assertEquals(new Point3i(0, 0, 0), v.cell());
        assertEquals(1, v.cube());
        assertEquals(1, v.member());
        assertEquals(1, v.tick());
    }

    @Test
    public void passesOnConservativeSyntheticExchange() {
        Necronomata automaton = freshAutomaton();
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 0] = 3f;
            frequency[base + 1] = 0f;
        });

        ConservationAudit audit = new ConservationAudit(automaton);

        // Hand-written two-member exchange: member 0 gives one quantum to
        // member 1. Total is conserved.
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 0] = frequency[base + 0] - 1f;
            frequency[base + 1] = frequency[base + 1] + 1f;
        });

        ConservationAudit.AuditResult result = audit.auditTick(1);

        assertTrue("conservative exchange must audit clean: "
                   + result.violations(), result.isClean());
        assertTrue(result.violations().isEmpty());
    }

    @Test
    public void perDirectionCountsSumToTotal() {
        CollisionStatistics stats = new CollisionStatistics();
        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = new Point3i(1, 1, 0);

        stats.recordCollision(cellA, 0, 0, cellB, 0, 1, 1, 1L, 1);
        stats.recordCollision(cellA, 0, 2, cellB, 0, 3, -1, 1L, 1);
        stats.recordCollision(cellA, 0, 4, cellB, 0, 5, 4, 2L, 2);
        stats.recordCollision(cellA, 1, 0, cellB, 1, 1, 4, 1L, 3);
        stats.recordCollision(cellA, 1, 2, cellB, 1, 3, -6, 3L, 3);

        long total = stats.totalCollisions();
        long sum = 0L;
        for (long count : stats.collisionsPerDirection().values()) {
            sum += count;
        }

        assertEquals(total, sum);
        assertEquals(5L, total);
        assertEquals(2L, stats.collisionsInDirection(4));
    }

    @Test
    public void auditIsExactNotApproximate() {
        Necronomata automaton = freshAutomaton();
        ConservationAudit audit = new ConservationAudit(automaton);

        // Corrupt the representation: a frequency slot that is no longer
        // an integer value in its float slot. This must be flagged as a
        // violation -- never silently absorbed by an epsilon tolerance.
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 10] = 0.0001f;
        });

        ConservationAudit.AuditResult result = audit.auditTick(1);

        assertFalse("non-integer quanta representation must be flagged, "
                    + "never rounded away by a tolerance", result.isClean());

        // A separate, clean lattice with no corruption must audit clean.
        Necronomata clean = freshAutomaton();
        ConservationAudit cleanAudit = new ConservationAudit(clean);
        ConservationAudit.AuditResult cleanResult = cleanAudit.auditTick(1);
        assertTrue(cleanResult.isClean());

        // Sanity: total is computed via exact long arithmetic, not float
        // summation (which would itself be an epsilon trap over many
        // entries). If this ever needs a tolerance, the state
        // representation is wrong -- that is exactly what this test
        // exists to catch.
        try {
            long total = cleanAudit.currentTotalQuanta();
            assertEquals(0L, total);
        } catch (Exception e) {
            fail("currentTotalQuanta() must not throw on a clean lattice: "
                 + e);
        }
    }

    /**
     * Regression for a bug caught by code review: {@code auditTick} was
     * diffing only against the PREVIOUS tick's snapshot, so a violation
     * that stabilizes (inject once, never touch again) fell out of the
     * per-slot diff on the very next tick and {@code isClean()} silently
     * flipped back to {@code true} even though the lattice-wide total
     * never returned to baseline. A leak-once collision rule must never
     * pass silently forever.
     */
    @Test
    public void auditDetectsPersistentDriftAcrossTicks() {
        Necronomata automaton = freshAutomaton();
        ConservationAudit audit = new ConservationAudit(automaton);

        // Plant one extra quantum, then never touch the lattice again.
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 3] = frequency[base + 3] + 1f;
        });

        ConservationAudit.AuditResult first = audit.auditTick(1);
        ConservationAudit.AuditResult second = audit.auditTick(2);

        assertFalse("tick 1 must report the freshly localized violation",
                    first.isClean());
        assertFalse("tick 2 must still report non-clean -- the drift never "
                    + "healed, it just stopped producing fresh per-slot "
                    + "deltas", second.isClean());
        assertEquals(second.baselineTotal() + 1, second.totalQuanta());
        assertEquals(Violation.Kind.RESIDUAL_DRIFT,
                     second.violations().get(0).kind());
    }

    /**
     * Regression for a false positive caught by the critic on the round-1
     * fix: once a drift is unhealed, {@code auditTick} was reporting
     * EVERY subsequent tick's per-slot changes as
     * {@code CONSERVATION_VIOLATION} even when that tick's own delta
     * summed to zero -- a legitimate balanced transfer occurring after an
     * earlier leak was mislabeled as two fresh violations. A tick whose
     * own net delta is zero must never emit {@code CONSERVATION_VIOLATION}
     * entries, regardless of a pre-existing unhealed drift; the persisting
     * divergence is still reported, but only via
     * {@code RESIDUAL_DRIFT}.
     */
    @Test
    public void balancedTransferAfterUnhealedDriftIsNotMislabeled() {
        Necronomata automaton = freshAutomaton();
        ConservationAudit audit = new ConservationAudit(automaton);

        // Tick 1: plant an unhealed leak.
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 3] = frequency[base + 3] + 1f;
        });
        ConservationAudit.AuditResult tick1 = audit.auditTick(1);
        assertFalse("tick 1 must report the freshly localized violation",
                    tick1.isClean());
        assertTrue("tick 1's violation must be a real localized violation",
                   tick1.violations().stream()
                        .anyMatch(v -> v.kind() == Violation.Kind.CONSERVATION_VIOLATION));

        // Tick 2: a legitimate balanced transfer (member0 -> member1) on
        // top of the still-unhealed tick-1 leak. Net delta for THIS tick
        // is zero.
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 0] = frequency[base + 0] - 1f;
            frequency[base + 1] = frequency[base + 1] + 1f;
        });
        ConservationAudit.AuditResult tick2 = audit.auditTick(2);

        assertFalse("tick 2 must still report non-clean: the tick-1 leak "
                    + "never healed", tick2.isClean());
        assertTrue("a balanced transfer must never be mislabeled as "
                   + "CONSERVATION_VIOLATION: " + tick2.violations(),
                   tick2.violations().stream()
                        .noneMatch(v -> v.kind() == Violation.Kind.CONSERVATION_VIOLATION));
        assertEquals(1, tick2.violations().size());
        assertEquals(Violation.Kind.RESIDUAL_DRIFT,
                     tick2.violations().get(0).kind());

        // Tick 3: no activity at all. Same expectation: RESIDUAL_DRIFT
        // only.
        ConservationAudit.AuditResult tick3 = audit.auditTick(3);
        assertFalse(tick3.isClean());
        assertEquals(1, tick3.violations().size());
        assertEquals(Violation.Kind.RESIDUAL_DRIFT,
                     tick3.violations().get(0).kind());
    }

    @Test
    public void strictModeThrowsOnPlantedViolation() {
        Necronomata nonStrictAutomaton = freshAutomaton();
        ConservationAudit nonStrictAudit = new ConservationAudit(nonStrictAutomaton);
        nonStrictAutomaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = nonStrictAutomaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 3] = frequency[base + 3] + 1f;
        });
        ConservationAudit.AuditResult expected = nonStrictAudit.auditTick(1);
        assertFalse(expected.isClean());

        Necronomata strictAutomaton = freshAutomaton();
        ConservationAudit strictAudit = new ConservationAudit(strictAutomaton,
                                                                true);
        assertTrue(strictAudit.isStrict());
        strictAutomaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = strictAutomaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 3] = frequency[base + 3] + 1f;
        });

        try {
            strictAudit.auditTick(1);
            fail("strict mode must throw on a planted violation");
        } catch (ConservationAudit.ConservationViolationException e) {
            assertFalse(e.result().isClean());
            assertEquals("strict-mode exception must carry the same result "
                         + "a non-strict call would have returned", expected,
                         e.result());
        }
    }

    @Test
    public void checkpointDoesNotResetBaselineInvariant() {
        Necronomata automaton = freshAutomaton();
        ConservationAudit audit = new ConservationAudit(automaton);

        // Plant a permanent drift.
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 5] = frequency[base + 5] + 1f;
        });
        ConservationAudit.AuditResult firstTick = audit.auditTick(1);
        assertFalse(firstTick.isClean());

        // Re-synchronize the localization reference to the current
        // (drifted) state -- checkpoint() must NOT erase the fact the
        // baseline invariant is still violated.
        audit.checkpoint();

        ConservationAudit.AuditResult afterCheckpoint = audit.auditTick(2);
        assertFalse("checkpoint() must not mask a persisting baseline "
                    + "divergence", afterCheckpoint.isClean());
        assertEquals(1, afterCheckpoint.violations().size());
        assertEquals(Violation.Kind.RESIDUAL_DRIFT,
                     afterCheckpoint.violations().get(0).kind());
    }

    /**
     * Routes a conservative transfer through the REAL {@code deltaF} +
     * {@code step()} path (not a direct {@code frequency} write), so a
     * future change to {@code step()}'s apply-then-zero contract is
     * caught here too, not only by NecronomataStateSemanticsTest.
     */
    @Test
    public void passesOnConservativeExchangeThroughStepPath() {
        Necronomata automaton = freshAutomaton();
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            frequency[base + 0] = 3f;
            frequency[base + 1] = 0f;
        });

        ConservationAudit audit = new ConservationAudit(automaton);

        automaton.process((angle, frequency, deltaA, deltaF) -> {
            int base = automaton.indexOfCell(new Point3i(0, 0, 0));
            deltaF[base + 0] = -1f;
            deltaF[base + 1] = 1f;
        });
        automaton.step();

        ConservationAudit.AuditResult result = audit.auditTick(1);

        assertTrue("conservative exchange via deltaF+step() must audit "
                   + "clean: " + result.violations(), result.isClean());
    }

    @Test
    public void collisionsForMemberPairIsSymmetric() {
        CollisionStatistics stats = new CollisionStatistics();
        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = new Point3i(1, 1, 0);

        // Cross-cube pair: (cube 0, member 2) <-> (cube 1, member 3),
        // recorded once in each side order -- symmetry must hold
        // regardless of which side of the collision each slot was passed
        // as.
        stats.recordCollision(cellA, 0, 2, cellB, 1, 3, 1, 1L, 1);
        stats.recordCollision(cellA, 1, 3, cellB, 0, 2, 1, 1L, 2);

        assertEquals(2L, stats.collisionsForMemberPair(0, 2, 1, 3));
        assertEquals(2L, stats.collisionsForMemberPair(1, 3, 0, 2));
        assertEquals(0L, stats.collisionsForMemberPair(0, 2, 0, 3));
    }

    /**
     * The bead's named proof test (inviscid-xew): a within-cube member
     * index alone is not a unique key -- two collisions sharing the same
     * within-cube {@code (memberA, memberB)} indices but occurring in
     * different cubes must land in different
     * {@code collisionsForMemberPair} buckets. Losing the cube axis would
     * silently conflate distinct face-type pairs.
     */
    @Test
    public void memberPairKeyDistinguishesCubes() {
        CollisionStatistics stats = new CollisionStatistics();
        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = new Point3i(1, 1, 0);

        // Same within-cube member indices (3, 4) in both collisions, but
        // different cubes (0 vs 1).
        stats.recordCollision(cellA, 0, 3, cellB, 0, 4, 1, 1L, 1);
        stats.recordCollision(cellA, 1, 3, cellB, 1, 4, 1, 1L, 2);

        assertEquals("cube-0 pair must be recorded independently of the "
                     + "cube-1 pair sharing the same member indices", 1L,
                     stats.collisionsForMemberPair(0, 3, 0, 4));
        assertEquals("cube-1 pair must be recorded independently of the "
                     + "cube-0 pair sharing the same member indices", 1L,
                     stats.collisionsForMemberPair(1, 3, 1, 4));
        assertEquals("a cross-cube combination that was never recorded "
                     + "must not be conflated with either same-cube bucket",
                     0L, stats.collisionsForMemberPair(0, 3, 1, 4));
        assertEquals(2L, stats.totalCollisions());
    }

    @Test
    public void recordCollisionValidatesCubeAndMemberRanges() {
        CollisionStatistics stats = new CollisionStatistics();
        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = new Point3i(1, 1, 0);

        assertRejectsWithArgumentName(() -> stats.recordCollision(cellA, -1,
                                                                    0, cellB,
                                                                    0, 0, 1,
                                                                    1L, 1),
                                       "cubeA");
        assertRejectsWithArgumentName(() -> stats.recordCollision(cellA, 5,
                                                                    0, cellB,
                                                                    0, 0, 1,
                                                                    1L, 1),
                                       "cubeA");
        assertRejectsWithArgumentName(() -> stats.recordCollision(cellA, 0,
                                                                    -1, cellB,
                                                                    0, 0, 1,
                                                                    1L, 1),
                                       "memberA");
        assertRejectsWithArgumentName(() -> stats.recordCollision(cellA, 0,
                                                                    6, cellB,
                                                                    0, 0, 1,
                                                                    1L, 1),
                                       "memberA");
        assertRejectsWithArgumentName(() -> stats.recordCollision(cellA, 0,
                                                                    0, cellB,
                                                                    -1, 0, 1,
                                                                    1L, 1),
                                       "cubeB");
        assertRejectsWithArgumentName(() -> stats.recordCollision(cellA, 0,
                                                                    0, cellB,
                                                                    5, 0, 1,
                                                                    1L, 1),
                                       "cubeB");
        assertRejectsWithArgumentName(() -> stats.recordCollision(cellA, 0,
                                                                    0, cellB,
                                                                    0, -1, 1,
                                                                    1L, 1),
                                       "memberB");
        assertRejectsWithArgumentName(() -> stats.recordCollision(cellA, 0,
                                                                    0, cellB,
                                                                    0, 6, 1,
                                                                    1L, 1),
                                       "memberB");
        assertEquals("none of the rejected calls may have been recorded",
                     0L, stats.totalCollisions());
    }

    private static void assertRejectsWithArgumentName(Runnable call,
                                                        String argName) {
        try {
            call.run();
            fail("expected IllegalArgumentException naming " + argName);
        } catch (IllegalArgumentException expected) {
            assertTrue("exception message must name the rejected argument ("
                       + argName + "): " + expected.getMessage(),
                       expected.getMessage().contains(argName));
        }
    }

    @Test
    public void transferMagnitudeHistogramExactCounts() {
        CollisionStatistics stats = new CollisionStatistics();
        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = new Point3i(1, 1, 0);

        long[] magnitudes = { 1L, 1L, 2L, 3L, 3L, 3L };
        int tick = 1;
        for (long magnitude : magnitudes) {
            stats.recordCollision(cellA, 0, 0, cellB, 0, 1, 1, magnitude,
                                   tick++);
        }

        Map<Long, Long> histogram = stats.transferMagnitudeHistogram();
        assertEquals(3, histogram.size());
        assertEquals(Long.valueOf(2L), histogram.get(1L));
        assertEquals(Long.valueOf(1L), histogram.get(2L));
        assertEquals(Long.valueOf(3L), histogram.get(3L));
    }

    @Test
    public void meanFreePathProxyDegenerateAndExact() {
        CollisionStatistics empty = new CollisionStatistics();
        assertTrue("zero collisions must yield NaN, never a divide-by-zero "
                   + "artifact", Double.isNaN(empty.meanFreePathProxy()));

        CollisionStatistics stats = new CollisionStatistics();
        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = new Point3i(1, 1, 0);
        // 3 collisions spanning ticks 1..5 inclusive -> span 5, proxy 5/3.
        stats.recordCollision(cellA, 0, 0, cellB, 0, 1, 1, 1L, 1);
        stats.recordCollision(cellA, 0, 0, cellB, 0, 1, 1, 1L, 2);
        stats.recordCollision(cellA, 0, 0, cellB, 0, 1, 1, 1L, 5);

        assertEquals(5.0 / 3.0, stats.meanFreePathProxy(), 0.0);
    }

    @Test
    public void recordCollisionRejectsNegativeTransferMagnitude() {
        CollisionStatistics stats = new CollisionStatistics();
        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = new Point3i(1, 1, 0);
        try {
            stats.recordCollision(cellA, 0, 0, cellB, 0, 1, 1, -1L, 1);
            fail("negative transferMagnitude must be rejected");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }
}
