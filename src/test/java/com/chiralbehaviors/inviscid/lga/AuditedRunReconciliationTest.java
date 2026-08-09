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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertThrows;

import java.util.List;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.lga.CollisionSweep.ReconciliationException;
import com.chiralbehaviors.inviscid.lga.CollisionSweep.TickResult;
import com.chiralbehaviors.inviscid.measure.AuditedRun;
import com.chiralbehaviors.inviscid.measure.CollisionStatistics;
import com.chiralbehaviors.inviscid.measure.ConservationAudit;
import com.chiralbehaviors.inviscid.measure.ConservationAudit.AuditResult;
import com.chiralbehaviors.inviscid.measure.ConservationAudit.LedgerEntry;

/**
 * The reconciliation NEGATIVE CONTROLS for bead inviscid-0nx.15's mandatory
 * pre-close requirement (that bead's NOTES, from the .14 review): {@link
 * CollisionSweep}'s TWO independent reconciliation cross-checks (see that
 * class's own Javadoc, "Two independent cross-checks") must each be
 * demonstrated to actually catch its own named failure class when driven
 * through the wiring that invokes {@link
 * CollisionSweep#reconcileWithLedger(TickReport, long)} per
 * tick alongside {@code ConservationAudit.auditTick} -- {@link AuditedRun}.
 * <ul>
 * <li>{@link #deliberatelyMisRecordedTransferIsCaughtByTheWiredReconciliation()}
 * -- the RECORDING-integrity check (a lie about what magnitude was
 * recorded to {@code CollisionStatistics}).</li>
 * <li>{@link #directFrequencyCorruptionIsCaughtByTheLedgerComparison()} --
 * the LEDGER-comparison check specifically (a stray write that bypasses
 * {@code deltaF} entirely, corrupting the real lattice total that {@code
 * reconcileWithLedger} observes, while the resolved contacts' own {@code
 * signedTransferTotal} stays provably zero).</li>
 * </ul>
 *
 * <p>Lives in the {@code lga} package -- not alongside {@link AuditedRun}
 * in {@code measure}, where {@code com.chiralbehaviors.inviscid.measure.
 * AuditedRunTest} covers the positive path -- because the first test's
 * seam, {@link CollisionSweep#magnitudeToRecord(CollisionRule.Delta)}, is
 * deliberately package-private to {@code lga} (see that method's own
 * Javadoc, "Visibility contract": "ONLY so the same-package reconciliation
 * negative-control test... can override it"). A test class in {@code
 * measure} cannot override a package-private {@code lga} method at all;
 * {@link AuditedRun} itself is fully public and lives in {@code measure}
 * regardless of where this test is compiled from. The second test does not
 * strictly need this package, but is kept alongside the first since both
 * are the same class of negative control (the two halves of {@code
 * CollisionSweep}'s reconciliation contract).
 *
 * @author halhildebrand
 */
public class AuditedRunReconciliationTest {

    private static final double RADIUS           = LgaTestGeometry.BASELINE_RADIUS;
    private static final int    RESOLUTION       = 360;
    private static final int    MEMBERS_PER_CUBE = 6;

    private static MemberGeometry newGeometry() {
        return new MemberGeometry(RESOLUTION, RADIUS);
    }

    private static ContactPredicate newPredicate() {
        return new ContactPredicate(newGeometry());
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

    /**
     * A deliberately mis-recorded transfer -- via the same {@code
     * CollisionSweep.magnitudeToRecord} seam {@code
     * CollisionSweepTest.deliberatelyMisRecordedTransferIsCaught} exercises
     * at the {@code CollisionSweep} unit level -- must be caught when
     * driven through the FULL {@link AuditedRun} per-tick wiring ({@code
     * HybridAutomaton.tick()} -&gt; {@code CollisionSweep.tick()} -&gt;
     * {@code ConservationAudit.auditTick} -&gt; {@code
     * CollisionSweep.reconcileWithLedger}), not just in isolation.
     *
     * <p><b>Honesty note on which check fires.</b> {@code
     * CollisionSweep.tick(int)}'s own recording-integrity cross-check
     * (independently-computed applied-magnitude total vs. what was
     * actually recorded) fires BEFORE this run ever reaches {@code
     * ConservationAudit.auditTick} or {@code reconcileWithLedger} -- the
     * lie is in what gets recorded to {@code CollisionStatistics}, not in
     * the real {@code deltaF}/{@code frequency} values the ledger
     * observes, so the ledger-delta comparison structurally cannot be the
     * one that catches THIS particular seam (see {@code CollisionSweep}'s
     * own class Javadoc, "Two independent cross-checks"). What this test
     * proves is that {@link AuditedRun}'s per-tick wiring does not mask or
     * swallow that exception -- it propagates cleanly out of {@code
     * run.tick(...)}, all the way through both classes' composition.
     */
    @Test
    public void deliberatelyMisRecordedTransferIsCaughtByTheWiredReconciliation() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);

        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = neighborhood.neighbor(cellA, 1);
        float angle = 1.7627826f;
        seedAngleAndQuanta(automaton, cellA, 3, 1, angle, 5f);
        seedAngleAndQuanta(automaton, cellB, 3, 0, angle, 2f);

        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep lyingSweep = new CollisionSweep(automaton, scan,
                                                         new QuantaExchangeRule(),
                                                         statistics) {
            @Override
            long magnitudeToRecord(CollisionRule.Delta delta) {
                return Math.abs(delta.deltaA()) + 1;
            }
        };
        HybridAutomaton hybrid = new HybridAutomaton(automaton, lyingSweep);
        ConservationAudit audit = new ConservationAudit(automaton);
        AuditedRun run = new AuditedRun(hybrid, audit);

        assertThrows(ReconciliationException.class, () -> run.tick(0));
    }

    /**
     * FIX 2 (2026-08-08 substantive-critic round, Significant): proves the
     * LEDGER-comparison half of {@code CollisionSweep}'s reconciliation
     * contract actually fires for its own named failure class -- {@code
     * reconcileWithLedger}'s own Javadoc, "Failure mode this guards
     * against": "a stray write bypassing {@code deltaF}" -- using a REAL
     * in-situ corruption, not a fabricated {@code long} handed straight to
     * the static method (that narrower claim is already covered by {@code
     * CollisionSweepTest.reconcileWithLedgerThrowsOnMismatch}, which does
     * pass a literal mismatched value).
     *
     * <p><b>Structure.</b> A real, non-lying {@link CollisionSweep} resolves
     * one genuine contact via {@link HybridAutomaton#tick(int)} -- its
     * {@link TickResult#signedTransferTotal()} is therefore exactly zero by
     * {@link CollisionRule.Delta}'s own construction-time invariant, and
     * nothing about that resolution is tampered with. ONLY AFTER that tick
     * completes, a frequency slot NOT involved in the resolved contact is
     * corrupted directly via {@code Necronomata.process(Processor)} --
     * bypassing {@code deltaF} entirely, exactly the failure mode {@code
     * reconcileWithLedger}'s Javadoc names. A NON-STRICT {@link
     * ConservationAudit} (the default constructor) is used so {@code
     * auditTick} reports the corruption as a violation via {@code
     * AuditResult.isClean()} rather than throwing {@code
     * ConservationViolationException} itself -- that exception must not
     * pre-empt reaching the ledger-comparison assertion below (bead's
     * explicit instruction).
     *
     * <p>The real, audit-computed ledger delta for this tick (from {@code
     * ConservationAudit.ledger()}, not a hand-picked literal) is then fed
     * into {@link CollisionSweep#reconcileWithLedger(TickReport, long)}
     * directly -- the resolved tick's zero {@code signedTransferTotal}
     * against the corrupted lattice's genuinely nonzero observed delta --
     * and asserted to throw. Sanity assertions along the way (zero signed
     * transfer total, nonzero real ledger delta, a genuinely non-clean
     * audit result) pin down that this is authentic in-situ corruption
     * reaching a real mismatch, not a vacuously-passing assertion.
     *
     * <p><b>Non-vacuity, verified by deletion (not left in the file):</b>
     * removing the corruption step (or asserting {@code
     * reconcileWithLedger} does NOT throw) makes this test fail -- the
     * uncorrupted path has {@code ledgerDelta == 0 == signedTransferTotal}
     * and {@code reconcileWithLedger} does not throw, so {@code
     * assertThrows} itself fails with no exception thrown. This was
     * verified manually during development (temporarily commenting out the
     * corruption block reproduces exactly that failure) rather than kept
     * as a permanent mutation in this file.
     */
    @Test
    public void directFrequencyCorruptionIsCaughtByTheLedgerComparison() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);

        Point3i cellA = new Point3i(0, 0, 0);
        Point3i cellB = neighborhood.neighbor(cellA, 1);
        float angle = 1.7627826f;
        seedAngleAndQuanta(automaton, cellA, 3, 1, angle, 5f);
        seedAngleAndQuanta(automaton, cellB, 3, 0, angle, 2f);

        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);
        // Non-strict (default): auditTick() must report, not throw, so the
        // ledger-comparison assertion below is actually reached.
        ConservationAudit audit = new ConservationAudit(automaton);

        TickResult collisionResult = hybrid.tick(0);
        assertEquals("the resolved tick's own signed transfer total must stay provably zero -- nothing about the collision resolution itself is tampered with",
                     0L, collisionResult.signedTransferTotal());

        // Deliberate corruption AFTER the tick completes: a direct
        // frequency write bypassing deltaF entirely, at a member NOT part
        // of the resolved contact -- the exact failure mode
        // reconcileWithLedger's Javadoc names ("a stray write bypassing
        // deltaF").
        Point3i untouchedCell = new Point3i(2, 2, 2);
        automaton.process((angleArray, frequency, deltaA, deltaF) -> {
            int index = automaton.indexOfCell(untouchedCell);
            frequency[index] = frequency[index] + 3f;
        });

        AuditResult auditResult = audit.auditTick(0);
        assertFalse("the corruption must be a real, detectable conservation violation (non-strict audit: reported, not thrown)",
                    auditResult.isClean());

        List<LedgerEntry> ledger = audit.ledger();
        LedgerEntry entry = ledger.get(ledger.size() - 1);
        long ledgerDelta = entry.totalAfter() - entry.totalBefore();
        assertNotEquals("the real, audit-observed ledger delta for this tick must be genuinely nonzero -- authentic corruption, not a fabricated mismatch",
                         0L, ledgerDelta);

        assertThrows("the ledger-comparison check must fire for real in-situ corruption",
                     ReconciliationException.class,
                     () -> CollisionSweep.reconcileWithLedger(collisionResult,
                                                               ledgerDelta));
    }
}
