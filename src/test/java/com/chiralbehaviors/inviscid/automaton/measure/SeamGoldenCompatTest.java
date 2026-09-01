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

import java.util.Map;
import java.util.Random;
import java.util.function.Consumer;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.automaton.Necronomata;
import com.chiralbehaviors.inviscid.automaton.QuantaField;
import com.chiralbehaviors.inviscid.automaton.lga.CollisionSweep;
import com.chiralbehaviors.inviscid.automaton.lga.ContactPredicate;
import com.chiralbehaviors.inviscid.automaton.lga.ContactScan;
import com.chiralbehaviors.inviscid.automaton.lga.FccNeighborhood;
import com.chiralbehaviors.inviscid.automaton.lga.HybridAutomaton;
import com.chiralbehaviors.inviscid.automaton.lga.MemberGeometry;
import com.chiralbehaviors.inviscid.automaton.lga.QuantaExchangeRule;
import com.chiralbehaviors.inviscid.automaton.measure.AnisotropyProbe.EstimatorResult;
import com.chiralbehaviors.inviscid.automaton.measure.AnisotropyProbe.SeedResult;
import com.chiralbehaviors.inviscid.automaton.measure.AuditedRun.TickOutcome;

/**
 * Golden-compatibility pins for the {@link QuantaField} / {@link
 * com.chiralbehaviors.inviscid.automaton.lga.TickReport} seam (bead inviscid-ckn /
 * inviscid-0nx.21, T2 design-ckn-lattice-seam.md §6.1-§6.3). Tests
 * numbered per the design memo's test list (6-10; 1-5 live in
 * {@link QuantaFieldSeamTest}).
 *
 * <p><b>Pin provenance (§6.3's ordering constraint).</b> Every literal
 * pinned below was captured against commit {@code c1b6e92} -- HEAD at
 * the moment inviscid-0nx.20 closed CLEAN (full {@code mvn test} green,
 * clean working tree) and BEFORE any inviscid-0nx.21 seam edit landed --
 * via a temporary, since-deleted harness that called the pre-seam
 * {@code coarseGrainedField(Necronomata)}, {@code
 * ConservationAudit(Necronomata)}, and {@code AnisotropyProbe.runOneSeed}
 * APIs directly. Pins captured after the seam lands would prove only
 * that the code equals itself (§6.3); these were not.
 *
 * <p><b>Exactness argument (§6.1).</b> Zero tolerance ({@code 0.0} /
 * exact {@code long} equality) is legitimate, not a fudge: every
 * legitimate {@code frequency} slot is an exact integer, so {@code
 * (double)Math.round((double)v) == (double)v} bit-for-bit, and Java is
 * strictfp-by-default (JEP 306, since Java 17; this project targets
 * Java 20) so summation order cannot be legally reassociated. Together
 * these prove the seam is output-identical for the quanta path; the
 * pins below are the empirical confirmation of that proof.
 *
 * @author halhildebrand
 */
public class SeamGoldenCompatTest {

    private static final double RADIUS     = 0.015;
    private static final int    RESOLUTION = 360;

    // --- Test 7 pin: coarseGrainedField, extent (4,4,4), deterministic
    // frequency[i] = (i % 7) - 3 pattern, captured pre-seam at c1b6e92.
    private static final double[] COARSE_GRAINED_FIELD_PIN = {
        -5.0, 0.0, 3.0, 0.0, 0.0, 1.0, 0.0, -5.0, -1.0, 0.0, 0.0, 0.0, 0.0,
        5.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 5.0, 0.0, -1.0, 0.0, 0.0, -3.0,
        0.0, 5.0, -5.0, 0.0, 3.0, 0.0, -3.0, 0.0, 5.0, 0.0, 0.0, 3.0, 0.0,
        -3.0, 1.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -5.0, 0.0, 3.0,
        0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -3.0, 0.0, 5.0, 0.0 };

    // --- Test 8 pin: ConservationAudit ledger, extent (4,4,4), seed 42L
    // (angles + quanta bound 6), 10 ticks of a real HybridAutomaton run,
    // captured pre-seam at c1b6e92. Every tick was clean: before==after
    // ==4696, drift==0 (a correctly-behaving hybrid run is conservation-
    // clean and ledger-reconciled every tick by construction, matching
    // AuditedRunTest.reconciliationWiredPerTickOverALongRun's finding).
    private static final long   LEDGER_TOTAL             = 4696L;
    private static final long   LEDGER_TOTAL_COLLISIONS   = 40L;
    private static final long   LEDGER_EFFECTIVE_COLLISIONS = 26L;

    // --- Test 9 pin: AnisotropyProbe.runOneSeed, extent (6,6,6), seed
    // 42L, 16 ticks, packetQuanta 100, captured pre-seam at c1b6e92.
    // Per T2 design-ckn-lattice-seam.md §6.2 test 9 / bead .22's plan-
    // audit F4 correction: pinned as a String.format("%.9e", v) literal,
    // compared as STRINGS -- never a raw-double byte/tolerance compare,
    // since a 9-decimal-digit-formatted literal is not bit-identical to
    // the double it was rendered from.
    private static final String TRANSPORT_RATIO_PIN         = "1.831884058e+00";
    private static final String TRANSPORT_X111_PIN          = "1.032679739e-04";
    private static final String TRANSPORT_X100_PIN          = "5.637254902e-05";
    private static final String TRANSPORT_X110_PIN          = "6.029411765e-05";
    private static final long   RUN_ONE_SEED_TOTAL_COLLISIONS = 246L;
    private static final long   RUN_ONE_SEED_EFFECTIVE_COLLISIONS = 6L;

    private static MemberGeometry newGeometry() {
        return new MemberGeometry(RESOLUTION, RADIUS);
    }

    private static void seedRandomAngles(Necronomata automaton, Point3i extent,
                                          long seed) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[length];
        for (int i = 0; i < length; i++) {
            angles[i] = random.nextFloat() * (float) (2 * Math.PI);
        }
        automaton.process((angleArray, frequency, deltaA,
                            deltaF) -> System.arraycopy(angles, 0,
                                                         angleArray, 0,
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
        automaton.process((angleArray, frequency, deltaA,
                            deltaF) -> System.arraycopy(quanta, 0, frequency,
                                                         0, length));
    }

    /** Test 6 -- pins the §6.1 delta (the single behavioural difference). */
    @Test
    public void corruptSlotIsRoundedNotSummedFractionallyThroughTheSeam() {
        // A 2x2x2 extent has exactly one even-parity cell: (0,0,0).
        Point3i extent = new Point3i(2, 2, 2);
        Necronomata automaton = new Necronomata(extent);
        automaton.process((angle, frequency, deltaA, deltaF) -> frequency[3] = 1.5f);

        double[] field = StructureFactor.coarseGrainedField(automaton);

        // Pre-seam behaviour summed the raw float (1.5) fractionally.
        // Through the seam (quantaAt -> long), the corrupt slot is
        // ROUNDED first: Math.round(1.5) == 2. This is the documented,
        // deliberate delta -- not an accident.
        assertEquals("corrupt slot must be rounded, not summed fractionally",
                     2.0, field[0], 0.0);
    }

    /** Test 7. */
    @Test
    public void coarseGrainedFieldThroughTheSeamMatchesPinnedPreSeamSnapshots() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < frequency.length; i++) {
                frequency[i] = (i % 7) - 3;
            }
        });

        double[] field = StructureFactor.coarseGrainedField(automaton);
        assertArrayEqualsExact(COARSE_GRAINED_FIELD_PIN, field);
    }

    /** Test 8. */
    @Test
    public void conservationLedgerThroughTheSeamMatchesThePinnedPreSeamLedger() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 42L);
        seedRandomQuanta(automaton, extent, 42L, 6);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            new ContactPredicate(newGeometry()));
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);
        ConservationAudit audit = new ConservationAudit(automaton);
        AuditedRun run = new AuditedRun(hybrid, audit);

        for (int tick = 0; tick < 10; tick++) {
            TickOutcome outcome = run.tick(tick);
            ConservationAudit.LedgerEntry entry = audit.ledger()
                                                         .get(audit.ledger()
                                                                   .size()
                                                              - 1);
            assertEquals("tick " + tick, LEDGER_TOTAL, entry.totalBefore());
            assertEquals("tick " + tick, LEDGER_TOTAL, entry.totalAfter());
            assertEquals("tick " + tick, 0L, entry.cumulativeDrift());
            assertTrue("tick " + tick + " must be clean",
                       outcome.auditResult().isClean());
        }
        assertEquals(LEDGER_TOTAL_COLLISIONS, statistics.totalCollisions());
        assertEquals(LEDGER_EFFECTIVE_COLLISIONS,
                     statistics.effectiveCollisions());
    }

    /**
     * Test 9. Also the pin for §4's RNG draw-order rule: {@code
     * runOneSeed} still calls {@code seedRandomAngles} then {@code
     * seedPacket}, exactly as at capture time.
     */
    @Test
    public void runOneSeedThroughTheSeamMatchesPinnedPhaseANumerics() {
        Point3i extent = new Point3i(6, 6, 6);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);
        SeedResult result = AnisotropyProbe.runOneSeed(extent, 42L, 16, 100,
                                                         origin);

        EstimatorResult transport = result.transport();
        assertTrue(transport.ratio().isPresent());
        assertEquals(TRANSPORT_RATIO_PIN,
                     String.format("%.9e", transport.ratio().getAsDouble()));
        Map<StructureFactor.Direction, AnisotropyProbe.DirectionMagnitude> perDirection = transport.perDirection();
        assertEquals(TRANSPORT_X111_PIN,
                     String.format("%.9e",
                                   perDirection.get(StructureFactor.Direction.X111)
                                               .magnitude()));
        assertEquals(TRANSPORT_X100_PIN,
                     String.format("%.9e",
                                   perDirection.get(StructureFactor.Direction.X100)
                                               .magnitude()));
        assertEquals(TRANSPORT_X110_PIN,
                     String.format("%.9e",
                                   perDirection.get(StructureFactor.Direction.X110)
                                               .magnitude()));

        assertEquals(RUN_ONE_SEED_TOTAL_COLLISIONS, result.totalCollisions());
        assertEquals(RUN_ONE_SEED_EFFECTIVE_COLLISIONS,
                     result.effectiveCollisions());
    }

    /**
     * Test 10 -- non-vacuity for tests 7-9: an injected faulty
     * {@link QuantaField} that perturbs one slot must turn the
     * compatibility pin red, proving the exact-equality assertions above
     * are not passing vacuously (e.g. via a broken/empty comparison).
     */
    @Test
    public void deliberatelyOffByOneQuantaFieldBreaksTheCompatPin() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < frequency.length; i++) {
                frequency[i] = (i % 7) - 3;
            }
        });

        QuantaField offByOne = new QuantaField() {
            @Override
            public Point3i extent() {
                return automaton.extent();
            }

            @Override
            public int slotCount() {
                return automaton.slotCount();
            }

            @Override
            public long quantaAt(int slot) {
                // Deliberately wrong: perturb exactly one slot by +1.
                long v = automaton.quantaAt(slot);
                return slot == 0 ? v + 1 : v;
            }

            @Override
            public boolean isExactAt(int slot) {
                return automaton.isExactAt(slot);
            }

            @Override
            public float phaseAt(int slot) {
                return automaton.phaseAt(slot);
            }

            @Override
            public int phaseResolution() {
                return automaton.phaseResolution();
            }

            @Override
            public void forEachCell(Consumer<? super Point3i> action) {
                automaton.forEachCell(action);
            }

            @Override
            public int indexOfCell(Point3i cell) {
                return automaton.indexOfCell(cell);
            }
        };

        double[] field = StructureFactor.coarseGrainedField(offByOne);
        assertFalse("a +1-perturbed QuantaField must NOT reproduce the pinned snapshot",
                    exactArrayEquals(COARSE_GRAINED_FIELD_PIN, field));
    }

    private static void assertArrayEqualsExact(double[] expected,
                                                 double[] actual) {
        assertEquals("array length", expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals("index " + i, expected[i], actual[i], 0.0);
        }
    }

    private static boolean exactArrayEquals(double[] a, double[] b) {
        if (a.length != b.length) {
            return false;
        }
        for (int i = 0; i < a.length; i++) {
            if (a[i] != b[i]) {
                return false;
            }
        }
        return true;
    }
}
