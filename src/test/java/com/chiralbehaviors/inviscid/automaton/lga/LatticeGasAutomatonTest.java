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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import javax.vecmath.Point3i;

import org.junit.BeforeClass;
import org.junit.Test;

import com.chiralbehaviors.inviscid.automaton.measure.CollisionStatistics;

/**
 * The 6 named failing tests for bead inviscid-0nx.21's {@link
 * LatticeGasAutomaton} (C.4: synchronous even-parity update). Uses the
 * COMMITTED atlas ({@code contact-atlas-v2.tsv}, N_lga=24, extent
 * (4,4,4)) - the real production artifact, loaded instantly rather than
 * regenerated - paired with the frozen conservation-exact {@link
 * CollisionTable} (bead inviscid-0nx.20).
 *
 * @author halhildebrand
 */
public class LatticeGasAutomatonTest {

    private static final String RESOURCE_PATH = "lga/contact-atlas-v2.tsv";

    private static ContactAtlas   ATLAS;
    private static CollisionTable COLLISIONS;
    private static Point3i        EXTENT;

    @BeforeClass
    public static void loadFixtures() throws IOException {
        URL resource = LatticeGasAutomatonTest.class.getClassLoader()
                                                      .getResource(RESOURCE_PATH);
        Path path = Paths.get(resource.getPath());
        ATLAS = ContactAtlas.read(path);
        COLLISIONS = CollisionTable.buildFromPhaseARule(new QuantaExchangeRule());
        EXTENT = ATLAS.header().extent();
    }

    private static LatticeGasAutomaton newAutomaton() {
        return new LatticeGasAutomaton(EXTENT, ATLAS, COLLISIONS,
                                        new CollisionStatistics());
    }

    /**
     * A {@link LatticeGasAutomaton} whose scan phase visits cells in the
     * REVERSE of the natural order - see {@link
     * LatticeGasAutomaton#cellVisitOrder()}'s Javadoc.
     */
    private static final class ReversedOrderAutomaton extends LatticeGasAutomaton {
        ReversedOrderAutomaton(Point3i extent, ContactAtlas atlas,
                                CollisionTable collisions,
                                CollisionStatistics statistics) {
            super(extent, atlas, collisions, statistics);
        }

        @Override
        List<Point3i> cellVisitOrder() {
            List<Point3i> natural = super.cellVisitOrder();
            List<Point3i> reversed = new java.util.ArrayList<>(natural);
            Collections.reverse(reversed);
            return reversed;
        }
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

    private static long totalQuanta(LatticeGasAutomaton automaton) {
        long total = 0L;
        for (int i = 0; i < automaton.slotCount(); i++) {
            total += automaton.quantaAt(i);
        }
        return total;
    }

    private static long[] quantaSnapshot(LatticeGasAutomaton automaton) {
        long[] snapshot = new long[automaton.slotCount()];
        for (int i = 0; i < snapshot.length; i++) {
            snapshot[i] = automaton.quantaAt(i);
        }
        return snapshot;
    }

    private static float[] phaseSnapshot(LatticeGasAutomaton automaton) {
        float[] snapshot = new float[automaton.slotCount()];
        for (int i = 0; i < snapshot.length; i++) {
            snapshot[i] = automaton.phaseAt(i);
        }
        return snapshot;
    }

    /**
     * Test 1: no cell sees a partially-updated lattice - proven by
     * running the SAME seeded initial condition through two automatons
     * that visit cells in opposite orders during the scan phase, and
     * asserting the post-tick state is bit-for-bit identical at EVERY
     * tick over a window long enough to contain a genuine same-tick
     * multi-touch (a member touched by more than one resolved collision
     * within a single tick) - see "Non-vacuity" below. If the algorithm
     * read partially-updated state mid-sweep, reversing traversal order
     * would change what a later-visited, repeat-touched member sees.
     *
     * <p><b>Non-vacuity (final-review Critical A fix).</b> A window with
     * ZERO same-tick multi-touches cannot distinguish a synchronous
     * update from a naive sequential one: with no member touched twice
     * in the same tick, there is nothing for scan order to disagree
     * about, and natural/reversed state would match even under a
     * regression that writes {@code quanta} live during the scan (the
     * critic's mutation-confirmed finding: tick(0) alone on this exact
     * seed/quantaBound has zero repeat-touches, so the single-tick
     * version of this test passed under that regression). Repeat-touches
     * were measured to begin around tick ~29 for this seed/quantaBound
     * and occur in the large majority of ticks thereafter, so this test
     * runs a 50-tick window (comfortable margin past the onset) and
     * explicitly asserts at least one same-tick multi-touch was observed
     * - instrumented via {@link LatticeGasAutomaton#lastTouchedIndices()}
     * - so a future regression that shrinks the window back to a
     * multi-touch-free range cannot silently pass.
     */
    @Test
    public void updateIsSynchronous() {
        LatticeGasAutomaton natural = newAutomaton();
        LatticeGasAutomaton reversed = new ReversedOrderAutomaton(EXTENT,
                                                                    ATLAS,
                                                                    COLLISIONS,
                                                                    new CollisionStatistics());
        seedDeterministic(natural, 42L, 6);
        seedDeterministic(reversed, 42L, 6);

        int ticksToRun = 50;
        boolean sawSameTickMultiTouch = false;
        for (int t = 0; t < ticksToRun; t++) {
            natural.tick(t);
            reversed.tick(t);

            assertArrayEqualsLong("tick " + t + " quanta",
                                   quantaSnapshot(natural),
                                   quantaSnapshot(reversed));
            assertArrayEqualsFloat("tick " + t + " phase",
                                    phaseSnapshot(natural),
                                    phaseSnapshot(reversed));

            if (hasRepeatIndex(natural.lastTouchedIndices())) {
                sawSameTickMultiTouch = true;
            }
        }

        assertTrue("non-vacuity: expected at least one same-tick multi-touch "
                   + "(a member touched by more than one resolved collision "
                   + "within a single tick) somewhere in the " + ticksToRun
                   + "-tick window - without one, this test cannot "
                   + "distinguish a synchronous update from a naive "
                   + "sequential one (see Javadoc)", sawSameTickMultiTouch);
    }

    private static boolean hasRepeatIndex(List<Integer> touched) {
        Set<Integer> seen = new HashSet<>();
        for (int index : touched) {
            if (!seen.add(index)) {
                return true;
            }
        }
        return false;
    }

    /** Test 2. */
    @Test
    public void totalQuantaConservedExactlyOver10000Ticks() {
        LatticeGasAutomaton automaton = newAutomaton();
        seedDeterministic(automaton, 7L, 6);

        long before = totalQuanta(automaton);
        for (int t = 0; t < 10000; t++) {
            automaton.tick(t);
        }
        long after = totalQuanta(automaton);

        assertEquals("total quanta must be conserved exactly over 10000 ticks",
                     before, after);
    }

    /** Test 3: full-state equality, checked every tick, not just the end. */
    @Test
    public void identicalSeedsGiveIdenticalTrajectories() {
        LatticeGasAutomaton a = newAutomaton();
        LatticeGasAutomaton b = newAutomaton();
        seedDeterministic(a, 1729L, 5);
        seedDeterministic(b, 1729L, 5);

        for (int t = 0; t < 300; t++) {
            a.tick(t);
            b.tick(t);
            assertArrayEqualsLong("tick " + t + " quanta",
                                   quantaSnapshot(a), quantaSnapshot(b));
            assertArrayEqualsFloat("tick " + t + " phase", phaseSnapshot(a),
                                    phaseSnapshot(b));
        }
    }

    /** Test 4: non-vacuity. */
    @Test
    public void collisionsOccur() {
        LatticeGasAutomaton automaton = newAutomaton();
        seedDeterministic(automaton, 42L, 6);

        for (int t = 0; t < 200; t++) {
            automaton.tick(t);
        }

        assertTrue("expected at least one resolved contact over 200 ticks",
                   automaton.statistics().totalCollisions() > 0);
        assertTrue("expected at least one non-no-op (effective) transfer over 200 ticks",
                   automaton.statistics().effectiveCollisions() > 0);
    }

    /**
     * Test 5: an all-zero-quanta lattice is static - no rotation
     * (phase advances by quanta, which is zero) AND no collisions ever
     * fire a nonzero transfer (equal-quanta contacts are the
     * QuantaExchangeRule's own defined no-op), regardless of phase.
     */
    @Test
    public void zeroQuantaLatticeIsStatic() {
        LatticeGasAutomaton automaton = newAutomaton();
        Random random = new Random(3L);
        int length = 30 * EXTENT.x * EXTENT.y * EXTENT.z;
        int[] phases = new int[length];
        for (int i = 0; i < length; i++) {
            phases[i] = random.nextInt(3600);
        }
        automaton.process((phase, quanta) -> System.arraycopy(phases, 0,
                                                                phase, 0,
                                                                length));

        long[] quantaBefore = quantaSnapshot(automaton);
        float[] phaseBefore = phaseSnapshot(automaton);

        for (int t = 0; t < 100; t++) {
            automaton.tick(t);
        }

        assertArrayEqualsLong("quanta", quantaBefore,
                               quantaSnapshot(automaton));
        assertArrayEqualsFloat("phase", phaseBefore,
                                phaseSnapshot(automaton));
    }

    /**
     * Test 6: no floating point in the tick path - verified structurally
     * (bead's own words: "assert by inspection/structure"), two ways:
     * (a) reflection confirms the state-carrying fields are {@code int[]}
     * / {@code long[]}, structurally incapable of holding a fractional
     * value; (b) recorded here, in the test itself, as the human-
     * readable verification record the bead asks for: {@link
     * LatticeGasAutomaton#tick(int)}'s source was read end-to-end and
     * contains zero {@code float}/{@code double} operators or locals -
     * the only floating-point arithmetic anywhere in the class is
     * {@link LatticeGasAutomaton#phaseAt(int)}, a measurement read
     * entirely outside {@code tick(int)}'s call graph (it is never
     * called BY {@code tick}).
     */
    @Test
    public void noFloatingPointInTheTickPath() throws NoSuchFieldException {
        Field phaseField = LatticeGasAutomaton.class.getDeclaredField("phase");
        Field quantaField = LatticeGasAutomaton.class.getDeclaredField("quanta");

        assertEquals("phase state must be int[]", int[].class,
                     phaseField.getType());
        assertEquals("quanta state must be long[]", long[].class,
                     quantaField.getType());
        assertTrue("phase field must not be a floating-point type",
                   !phaseField.getType().equals(float[].class)
                   && !phaseField.getType().equals(double[].class));
        assertTrue("quanta field must not be a floating-point type",
                   !quantaField.getType().equals(float[].class)
                   && !quantaField.getType().equals(double[].class));
        // Non-vacuity: the class must actually declare these fields as
        // instance state (private, non-static), not e.g. constants that
        // would make the above checks trivially pass on an empty class.
        assertTrue("phase must be an instance field",
                   !Modifier.isStatic(phaseField.getModifiers()));
        assertTrue("quanta must be an instance field",
                   !Modifier.isStatic(quantaField.getModifiers()));
    }

    /**
     * Header-pairing guard (bead inviscid-0nx.19's contract, T2
     * design-ckn-lattice-seam.md §9 step 10): a {@link ContactTable}
     * paired with a DIFFERENT atlas's header must be refused at
     * construction, not silently accepted.
     */
    @Test
    public void rejectsAContactTableHeaderMismatchedWithTheAtlasHeader() throws IOException {
        ContactAtlas otherAtlas = ContactAtlasGenerator.generate(12,
                                                                   new Point3i(4,
                                                                               4,
                                                                               4),
                                                                   42L, 50);
        ContactTable table = ContactTable.of(ATLAS);

        IllegalStateException thrown = assertThrows(IllegalStateException.class,
                                                       () -> new LatticeGasAutomaton(EXTENT,
                                                                                      table,
                                                                                      otherAtlas.header(),
                                                                                      COLLISIONS,
                                                                                      new CollisionStatistics()));
        assertTrue("expected the pairing violation named in the message: "
                   + thrown.getMessage(),
                   thrown.getMessage().contains("header pairing violated"));

        // Positive control: the SAME header pairs cleanly.
        new LatticeGasAutomaton(EXTENT, table, ATLAS.header(), COLLISIONS,
                                 new CollisionStatistics());
    }

    /**
     * Fine-step divisibility guard (USER DECISION 2026-08-08, FINAL:
     * {@code accumulator/10} must be exact integer arithmetic since
     * {@code phaseResolution (3600) == geometryResolution (360) * 10}
     * exactly - "a loud guard tying the two resolutions, divisibility-
     * style"). {@code geometryResolution=96} is chosen deliberately:
     * {@code 3600 % 96 != 0} violates THIS guard, while 96 satisfies every
     * OTHER divisibility constraint that could fire first and mask it. The
     * one that is still LIVE is {@link MemberGeometry}'s constructor
     * requirement that {@code resolution % 8 == 0} - {@code 96 = 8 * 12}
     * satisfies it. ({@code 24 | 96} was the original reason for this
     * value: it dodged the {@code nLga}-divides-{@code geometryResolution}
     * guard of the phase-quantizer class that bead inviscid-0nx.27
     * retired. That guard is gone and nothing replaced it, so the {@code
     * 24 | 96} property is now incidental rather than load-bearing; the
     * value is retained because a fixture that violates exactly one guard
     * at a time is the point.)
     */
    @Test
    public void rejectsAPhaseResolutionNotDivisibleByGeometryResolution() {
        ContactAtlas.Header badHeader = new ContactAtlas.Header(ContactAtlas.ATLAS_VERSION,
                                                                   "test",
                                                                   "test-commit",
                                                                   ContactAtlasGenerator.RADIUS,
                                                                   96,
                                                                   1.0,
                                                                   24,
                                                                   "Cubes[0]",
                                                                   EXTENT, 42L,
                                                                   0, 150);
        ContactAtlas badAtlas = new ContactAtlas(badHeader, java.util.List.of());
        ContactTable badTable = ContactTable.of(badAtlas);

        IllegalStateException thrown = assertThrows(IllegalStateException.class,
                                                       () -> new LatticeGasAutomaton(EXTENT,
                                                                                      badTable,
                                                                                      badHeader,
                                                                                      COLLISIONS,
                                                                                      new CollisionStatistics()));
        assertTrue("expected the divisibility violation named in the message: "
                   + thrown.getMessage(),
                   thrown.getMessage().contains("evenly divisible"));
    }

    private static void assertArrayEqualsLong(String label, long[] expected,
                                                long[] actual) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(label + " index " + i, expected[i], actual[i]);
        }
    }

    private static void assertArrayEqualsFloat(String label,
                                                 float[] expected,
                                                 float[] actual) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(label + " index " + i, expected[i], actual[i], 0.0f);
        }
    }
}
