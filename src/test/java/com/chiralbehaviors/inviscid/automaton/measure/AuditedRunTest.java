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

import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.automaton.Necronomata;
import com.chiralbehaviors.inviscid.automaton.QuantaField;
import com.chiralbehaviors.inviscid.automaton.lga.CollisionSweep;
import com.chiralbehaviors.inviscid.automaton.lga.ContactPredicate;
import com.chiralbehaviors.inviscid.automaton.lga.ContactScan;
import com.chiralbehaviors.inviscid.automaton.lga.FccNeighborhood;
import com.chiralbehaviors.inviscid.automaton.lga.HybridAutomaton;
import com.chiralbehaviors.inviscid.automaton.lga.LgaTestGeometry;
import com.chiralbehaviors.inviscid.automaton.lga.MemberGeometry;
import com.chiralbehaviors.inviscid.automaton.lga.QuantaExchangeRule;
import com.chiralbehaviors.inviscid.automaton.lga.TickDriver;
import com.chiralbehaviors.inviscid.automaton.lga.TickReport;
import com.chiralbehaviors.inviscid.automaton.measure.AuditedRun.TickOutcome;

/**
 * Behavioral tests for {@link AuditedRun} (bead inviscid-0nx.15's
 * pre-close requirement, per that bead's NOTES): {@link
 * CollisionSweep#reconcileWithLedger(TickReport, long)}
 * MUST be invoked per tick alongside {@link
 * ConservationAudit#auditTick(int)}.
 *
 * <p>The {@code CollisionSweep}-SPECIFIC reconciliation negative control
 * (a deliberately mis-recorded transfer, caught by the wiring) lives in
 * {@code com.chiralbehaviors.inviscid.automaton.lga.AuditedRunReconciliationTest},
 * not here: {@code CollisionSweep.magnitudeToRecord} -- the seam that
 * test exercises -- is deliberately package-private to {@code lga} (see
 * that method's own Javadoc, "Visibility contract"), so a subclass
 * overriding it can only be written from a test class in the {@code lga}
 * package, not from this ({@code measure}) package -- even though
 * {@link AuditedRun} itself lives in {@code measure} and is fully
 * public. The GENERIC reconciliation negative control -- proving the
 * WIDENED {@code reconcileWithLedger(TickReport, long)} contract checks
 * any {@link com.chiralbehaviors.inviscid.automaton.lga.TickDriver}, not just a
 * {@code CollisionSweep}-backed one -- lives here, see {@link
 * #genericLyingTickReportTripsReconciliationRegardlessOfCollisionSweep()}.
 *
 * @author halhildebrand
 */
public class AuditedRunTest {

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
     * Positive path: {@link AuditedRun#tick(int)} wires {@code
     * ConservationAudit.auditTick} and {@code
     * CollisionSweep.reconcileWithLedger} together, per tick, over a
     * seeded, contact-bearing 200-tick run, and neither check ever fires
     * (a correctly-behaving hybrid automaton is, by construction, both
     * conservation-clean and ledger-reconciled every tick). Non-vacuity:
     * collisions actually occurred (same reasoning as {@code
     * HybridAutomatonTest.collisionsActuallyOccur} -- a run with zero
     * collisions would pass this test vacuously).
     */
    @Test
    public void reconciliationWiredPerTickOverALongRun() {
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

        for (int tick = 0; tick < 200; tick++) {
            TickOutcome outcome = run.tick(tick);
            assertTrue("conservation violated at tick " + tick + ": "
                       + outcome.auditResult().violations(),
                       outcome.auditResult().isClean());
        }

        assertTrue("expected at least one resolved contact over 200 ticks",
                   statistics.totalCollisions() > 0);
        assertTrue("expected at least one non-no-op (effective) transfer over 200 ticks",
                   statistics.effectiveCollisions() > 0);
    }

    /**
     * A {@link TickReport} that LIES: it always claims a nonzero
     * {@code signedTransferTotal}, regardless of what actually happened
     * to the lattice this tick.
     */
    private static final class LyingTickReport implements TickReport {
        private final int tick;

        LyingTickReport(int tick) {
            this.tick = tick;
        }

        @Override
        public int tick() {
            return tick;
        }

        @Override
        public long signedTransferTotal() {
            return 1L;
        }
    }

    /**
     * A {@link TickDriver}, deliberately NOT backed by {@code
     * CollisionSweep}/{@code HybridAutomaton}, that never mutates
     * {@code field} and always reports a {@link LyingTickReport}.
     */
    private static final class StaticLyingTickDriver implements TickDriver {
        private final QuantaField field;

        StaticLyingTickDriver(QuantaField field) {
            this.field = field;
        }

        @Override
        public TickReport tick(int tickNumber) {
            return new LyingTickReport(tickNumber);
        }

        @Override
        public QuantaField field() {
            return field;
        }
    }

    /**
     * GENERIC reconciliation negative control (bead inviscid-ckn /
     * inviscid-0nx.21, critic finding on T2
     * critique-checkpoint-0nx21-steps0-7.md [21922]): the {@code
     * CollisionSweep}-specific negative control in {@code
     * AuditedRunReconciliationTest} exercises only a real, physically
     * zero-sum {@code TickResult} gone wrong -- it never proves the
     * WIDENED {@code reconcileWithLedger(TickReport, long)} contract
     * (bead inviscid-ckn) actually checks an arbitrary {@link
     * TickDriver}/{@link TickReport} pair with no relationship to
     * {@code CollisionSweep} at all. This test drives exactly that: a
     * driver that never touches the lattice (so {@code
     * ConservationAudit}'s real observed ledger delta is always
     * {@code 0}) paired with a report that unconditionally claims
     * {@code signedTransferTotal() == 1}. The mismatch MUST be caught.
     */
    @Test
    public void genericLyingTickReportTripsReconciliationRegardlessOfCollisionSweep() {
        Point3i extent = new Point3i(2, 2, 2);
        Necronomata field = new Necronomata(extent);
        StaticLyingTickDriver driver = new StaticLyingTickDriver(field);
        ConservationAudit audit = new ConservationAudit(field);
        AuditedRun run = new AuditedRun(driver, audit);

        assertThrows(CollisionSweep.ReconciliationException.class,
                     () -> run.tick(0));
    }
}
