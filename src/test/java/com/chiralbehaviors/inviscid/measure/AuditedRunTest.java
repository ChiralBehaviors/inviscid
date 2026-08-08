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

import static org.junit.Assert.assertTrue;

import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.lga.CollisionSweep;
import com.chiralbehaviors.inviscid.lga.ContactPredicate;
import com.chiralbehaviors.inviscid.lga.ContactScan;
import com.chiralbehaviors.inviscid.lga.FccNeighborhood;
import com.chiralbehaviors.inviscid.lga.HybridAutomaton;
import com.chiralbehaviors.inviscid.lga.MemberGeometry;
import com.chiralbehaviors.inviscid.lga.QuantaExchangeRule;
import com.chiralbehaviors.inviscid.measure.AuditedRun.TickOutcome;

/**
 * Behavioral tests for {@link AuditedRun} (bead inviscid-0nx.15's
 * pre-close requirement, per that bead's NOTES): {@link
 * CollisionSweep#reconcileWithLedger(CollisionSweep.TickResult, long)}
 * MUST be invoked per tick alongside {@link
 * ConservationAudit#auditTick(int)}.
 *
 * <p>The reconciliation NEGATIVE CONTROL (a deliberately mis-recorded
 * transfer, caught by the wiring) lives in {@code
 * com.chiralbehaviors.inviscid.lga.AuditedRunReconciliationTest}, not
 * here: {@code CollisionSweep.magnitudeToRecord} -- the seam that test
 * exercises -- is deliberately package-private to {@code lga} (see that
 * method's own Javadoc, "Visibility contract"), so a subclass overriding
 * it can only be written from a test class in the {@code lga} package,
 * not from this ({@code measure}) package -- even though {@link
 * AuditedRun} itself lives in {@code measure} and is fully public.
 *
 * @author halhildebrand
 */
public class AuditedRunTest {

    private static final double RADIUS           = 0.015;
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
}
