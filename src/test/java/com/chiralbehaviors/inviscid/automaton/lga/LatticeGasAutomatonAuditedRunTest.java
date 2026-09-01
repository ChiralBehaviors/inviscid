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

package com.chiralbehaviors.inviscid.automaton.lga;

import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.BeforeClass;
import org.junit.Test;

import com.chiralbehaviors.inviscid.automaton.measure.AuditedRun;
import com.chiralbehaviors.inviscid.automaton.measure.AuditedRun.TickOutcome;
import com.chiralbehaviors.inviscid.automaton.measure.CollisionStatistics;
import com.chiralbehaviors.inviscid.automaton.measure.ConservationAudit;

/**
 * Per the locked design of record: {@code AuditedRun} is THE way to run
 * an audited simulation, and bead inviscid-0nx.21's LGA driver MUST flow
 * through it, with {@code reconcileWithLedger} invoked per tick -- not
 * optional. This class proves that wiring for {@link LatticeGasAutomaton}
 * specifically: a positive path (R6, empirically: {@link
 * CollisionStatistics} populates identically to the hybrid path, and no
 * tick is ever flagged unclean or unreconciled), and the LGA-SPECIFIC
 * reconciliation negative control (distinct from {@code
 * AuditedRunTest}'s GENERIC fake-driver negative control, T2
 * critique-checkpoint-0nx21-steps0-7.md [21922] finding).
 *
 * @author halhildebrand
 */
public class LatticeGasAutomatonAuditedRunTest {

    private static final String RESOURCE_PATH = "lga/contact-atlas-v2.tsv";

    private static ContactAtlas   ATLAS;
    private static CollisionTable COLLISIONS;
    private static Point3i        EXTENT;

    /** Overrides the reconciliation seam to lie -- see {@link LatticeGasAutomaton#signedTransferTotalToReport}. */
    private static final class LyingLatticeGasAutomaton extends LatticeGasAutomaton {
        LyingLatticeGasAutomaton(Point3i extent, ContactAtlas atlas,
                                  CollisionTable collisions,
                                  CollisionStatistics statistics) {
            super(extent, atlas, collisions, statistics);
        }

        @Override
        long signedTransferTotalToReport(long computed) {
            return computed + 1L;
        }
    }

    @BeforeClass
    public static void loadFixtures() throws IOException {
        URL resource = LatticeGasAutomatonAuditedRunTest.class.getClassLoader()
                                                                .getResource(RESOURCE_PATH);
        Path path = Paths.get(resource.getPath());
        ATLAS = ContactAtlas.read(path);
        COLLISIONS = CollisionTable.buildFromPhaseARule(new QuantaExchangeRule());
        EXTENT = ATLAS.header().extent();
    }

    private static void seedDeterministic(LatticeGasAutomaton automaton,
                                           long seed, int quantaBound) {
        Random random = new Random(seed);
        int length = 30 * EXTENT.x * EXTENT.y * EXTENT.z;
        int[] phases = new int[length];
        long[] quanta = new long[length];
        for (int i = 0; i < length; i++) {
            phases[i] = random.nextInt(3600);
            quanta[i] = random.nextInt(2 * quantaBound + 1) - quantaBound;
        }
        automaton.process((phase, q) -> {
            System.arraycopy(phases, 0, phase, 0, length);
            System.arraycopy(quanta, 0, q, 0, length);
        });
    }

    /**
     * Positive path (mirrors {@code AuditedRunTest.reconciliationWiredPerTickOverALongRun}
     * for the LGA driver specifically): {@link ConservationAudit} stays
     * clean and {@code reconcileWithLedger} never fires over a real,
     * contact-bearing 200-tick audited run. R6, empirically: {@link
     * CollisionStatistics} populates via this driver exactly like it
     * does via {@code CollisionSweep} -- non-vacuity below.
     */
    @Test
    public void auditedRunStaysCleanOverALongLgaRun() {
        LatticeGasAutomaton lga = new LatticeGasAutomaton(EXTENT, ATLAS,
                                                            COLLISIONS,
                                                            new CollisionStatistics());
        seedDeterministic(lga, 42L, 6);
        ConservationAudit audit = new ConservationAudit(lga);
        AuditedRun run = new AuditedRun(lga, audit);

        for (int tick = 0; tick < 200; tick++) {
            TickOutcome outcome = run.tick(tick);
            assertTrue("conservation violated at tick " + tick + ": "
                       + outcome.auditResult().violations(),
                       outcome.auditResult().isClean());
        }

        assertTrue("R6: expected CollisionStatistics to be populated by the LGA driver, "
                   + "same as the hybrid path (empirical, not assumed)",
                   lga.statistics().totalCollisions() > 0);
        assertTrue("R6: expected at least one effective (non-no-op) transfer",
                   lga.statistics().effectiveCollisions() > 0);
    }

    /**
     * LGA-SPECIFIC negative control: a driver that lies about its
     * {@code signedTransferTotal} by exactly 1 must trip {@code
     * reconcileWithLedger}'s cross-check, proving the reconciliation
     * wiring is load-bearing for THIS driver, not merely inherited by
     * assumption from the generic {@code TickReport} contract test.
     */
    @Test
    public void lyingLgaSignedTransferTotalTripsReconciliation() {
        LatticeGasAutomaton lga = new LyingLatticeGasAutomaton(EXTENT, ATLAS,
                                                                  COLLISIONS,
                                                                  new CollisionStatistics());
        seedDeterministic(lga, 42L, 6);
        ConservationAudit audit = new ConservationAudit(lga);
        AuditedRun run = new AuditedRun(lga, audit);

        assertThrows(CollisionSweep.ReconciliationException.class,
                     () -> run.tick(0));
    }
}
