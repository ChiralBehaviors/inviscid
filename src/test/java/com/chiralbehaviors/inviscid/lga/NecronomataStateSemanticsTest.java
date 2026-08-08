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

import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;

/**
 * Pins the field contract of {@link Necronomata}'s 30-float-per-cell state:
 * {@code frequency} is the conserved integer quanta count, {@code deltaA} is
 * derived from it every tick via {@link Necronomata#QUANTUM_RATE}, and
 * {@code deltaF} is a transient collision accumulator that never survives a
 * tick.
 *
 * @author halhildebrand
 */
public class NecronomataStateSemanticsTest {

    /**
     * A member with nonzero frequency must advance its angle at
     * {@code QUANTUM_RATE * frequency} per step; members with zero frequency
     * must not move at all, regardless of their initial angle.
     */
    @Test
    public void frequencyDrivesAngularRate() {
        Point3i extent = new Point3i(1, 1, 2);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];
        for (int i = 0; i < len; i++) {
            angles[i] = i * 0.1f;
        }
        int m = 5;
        frequency[m] = 3f;
        float[] initialAngles = angles.clone();

        Necronomata automata = new Necronomata(angles, extent, frequency);

        int steps = 4;
        for (int s = 0; s < steps; s++) {
            automata.step();
        }

        float expected = initialAngles[m] + 3f * Necronomata.QUANTUM_RATE * steps;
        assertEquals("driven member must advance at frequency * QUANTUM_RATE per step", expected, angles[m], 1e-5f);

        for (int i = 0; i < len; i++) {
            if (i != m) {
                assertEquals("member " + i + " has zero frequency; its angle must not move", initialAngles[i],
                             angles[i], 1e-5f);
            }
        }
    }

    /**
     * The constructor must not seed rates from the initial angle positions:
     * with zero frequency everywhere, stepping must never advance any
     * angle by any rate derived from its position, no matter what the
     * initial angles were - deliberately including values beyond one
     * revolution (up to {@code 30*0.37 ~ 11.1} radians) to prove that.
     *
     * <p>Updated for inviscid-vb9: {@link Necronomata#step()} now
     * unconditionally wraps {@code angle} into {@code [0, 2*pi)} every
     * tick (the canonicalization is part of step()'s contract, not
     * conditioned on a nonzero rate), so a zero-frequency step no longer
     * leaves an out-of-range seed value byte-for-byte unchanged - it
     * normalizes it to its {@code mod 2*pi} residue. That normalization
     * is exactly what proves the point this test makes: the residue is a
     * pure function of the seed position, with no contribution from any
     * position-derived rate, so the assertion below compares against
     * each seed's own {@code mod 2*pi} value rather than the raw seed.
     */
    @Test
    public void constructorDoesNotSeedRatesFromPositions() {
        Point3i extent = new Point3i(1, 1, 1);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];
        for (int i = 0; i < len; i++) {
            angles[i] = (i + 1) * 0.37f;
        }
        float[] initialAngles = angles.clone();
        double twoPi = 2.0 * Math.PI;
        float[] expectedAngles = new float[len];
        for (int i = 0; i < len; i++) {
            double raw = initialAngles[i];
            expectedAngles[i] = (float) (raw - twoPi * Math.floor(raw / twoPi));
        }

        Necronomata automata = new Necronomata(angles, extent, frequency);
        automata.step();

        assertArrayEquals("with zero frequency everywhere, each angle after step() must equal its own seed value's "
                           + "mod-2*PI residue - no rate derived from position may perturb it further",
                           expectedAngles, angles, 1e-5f);
    }

    /**
     * deltaF is a transient collision-delta accumulator only: whatever a
     * collision writes into it is applied to frequency by the next
     * step(), and deltaF is zeroed so it never leaks into a subsequent
     * tick.
     */
    @Test
    public void deltaFIsATransientCollisionAccumulator() {
        Point3i extent = new Point3i(1, 1, 1);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];
        int m = 2;
        frequency[m] = 5f;

        Necronomata automata = new Necronomata(angles, extent, frequency);

        float collisionDelta = -2f;
        automata.process((angle, freq, deltaA, deltaF) -> deltaF[m] = collisionDelta);

        automata.step();

        assertEquals("frequency must be advanced by exactly the collision delta", 5f + collisionDelta, frequency[m],
                     1e-6f);
        automata.process((angle, freq, deltaA, deltaF) -> assertEquals("deltaF must be zeroed after being applied",
                                                                        0f, deltaF[m], 1e-6f));
    }

    /**
     * Total quanta (the sum of frequency across every member) must be
     * exactly conserved by step() absent any collision write to deltaF -
     * step() itself must never perturb frequency's total.
     */
    @Test
    public void totalQuantaConservedByStep() {
        Point3i extent = new Point3i(2, 2, 2);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];

        Random random = new Random(42L);
        long expectedSum = 0;
        for (int i = 0; i < len; i++) {
            int quanta = random.nextInt(21) - 10;
            frequency[i] = quanta;
            expectedSum += quanta;
        }

        Necronomata automata = new Necronomata(angles, extent, frequency);

        for (int s = 0; s < 100; s++) {
            automata.step();
        }

        long actualSum = 0;
        for (int i = 0; i < len; i++) {
            actualSum += (long) frequency[i];
        }

        assertEquals("total quanta must be exactly conserved by step()", expectedSum, actualSum);
    }

    /**
     * step() applies deltaF to frequency BEFORE recomputing deltaA from it:
     * a quantum absorbed via a collision this tick must move its member's
     * angle on this same step(), not the next one. The derived-rate
     * invariant {@code deltaA == QUANTUM_RATE * frequency} must hold
     * immediately after step() returns.
     */
    @Test
    public void collisionAffectsRateSameTick() {
        Point3i extent = new Point3i(1, 1, 1);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];
        int m = 4;
        float[] initialAngles = angles.clone();

        Necronomata automata = new Necronomata(angles, extent, frequency);

        automata.process((angle, freq, deltaA, deltaF) -> deltaF[m] = 2f);
        automata.step();

        float expectedAfterFirstStep = initialAngles[m] + 2f * Necronomata.QUANTUM_RATE;
        assertEquals("a quantum absorbed this tick must move its member on this same tick", expectedAfterFirstStep,
                     angles[m], 1e-5f);

        automata.process((angle, freq, deltaA, deltaF) -> assertEquals(
                                                                         "deltaA must equal QUANTUM_RATE * frequency immediately after step()",
                                                                         Necronomata.QUANTUM_RATE * freq[m],
                                                                         deltaA[m], 1e-6f));

        automata.step();

        float expectedAfterSecondStep = expectedAfterFirstStep + 2f * Necronomata.QUANTUM_RATE;
        assertEquals("frequency persists as a conserved rate; the member keeps moving at the same rate next tick",
                     expectedAfterSecondStep, angles[m], 1e-5f);
    }

    /**
     * A more realistic conservation scenario than {@link
     * #totalQuantaConservedByStep()}: every tick, quanta actually move
     * between a random pair of members via deltaF (not left untouched at
     * zero). Total quanta must still be exactly conserved, and every
     * frequency slot must remain integer-valued.
     */
    @Test
    public void totalQuantaConservedAcrossTransfers() {
        Point3i extent = new Point3i(2, 2, 2);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];

        Random random = new Random(42L);
        long expectedSum = 0;
        for (int i = 0; i < len; i++) {
            int quanta = random.nextInt(21) - 10;
            frequency[i] = quanta;
            expectedSum += quanta;
        }

        Necronomata automata = new Necronomata(angles, extent, frequency);

        for (int tick = 0; tick < 100; tick++) {
            int a = random.nextInt(len);
            int b = random.nextInt(len);
            int k = random.nextInt(5) + 1;
            automata.process((angle, freq, deltaA, deltaF) -> {
                deltaF[a] -= k;
                deltaF[b] += k;
            });
            automata.step();
        }

        long actualSum = 0;
        for (int i = 0; i < len; i++) {
            actualSum += (long) frequency[i];
            assertEquals("frequency must remain integer-valued at index " + i, Math.rint(frequency[i]),
                         frequency[i], 0f);
        }

        assertEquals("total quanta must be exactly conserved across transfers", expectedSum, actualSum);
    }

    /**
     * inviscid-vb9: angle is a float32 accumulator that must be wrapped
     * into [0, 2*PI) after every step() to eliminate the silent
     * phase-accumulation precision cliff (unbounded accumulation loses
     * rate accuracy well before 2^23 ticks; see Necronomata.step()
     * javadoc). Drives enough steps (1500 at frequency=5) that the
     * unwrapped angle would exceed 2*PI several times over, and checks
     * both that every intermediate angle stays in-range and that the
     * final wrapped value matches the mathematically expected angle mod
     * 2*PI.
     */
    @Test
    public void angleStaysWithinOneRevolution() {
        Point3i extent = new Point3i(1, 1, 1);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];
        int m = 7;
        int quanta = 5;
        frequency[m] = quanta;

        Necronomata automata = new Necronomata(angles, extent, frequency);

        int steps = 1500;
        double twoPi = 2.0 * Math.PI;
        for (int s = 0; s < steps; s++) {
            automata.step();
            assertTrue("angle must stay within [0, 2*PI) after every step, was "
                       + angles[m], angles[m] >= 0f && angles[m] < (float) twoPi);
        }

        double expected = expectedWrappedAngle(steps, quanta, twoPi);
        assertEquals("wrapped angle must match the mathematically expected angle mod 2*PI",
                     expected, angles[m], 5e-3);
    }

    /**
     * inviscid-vb9: negative frequency (negative quanta) must wrap into
     * [0, 2*PI) using floor-mod semantics, never landing on a negative
     * angle value.
     */
    @Test
    public void angleWrapsCorrectlyForNegativeFrequency() {
        Point3i extent = new Point3i(1, 1, 1);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];
        int m = 3;
        int quanta = -5;
        frequency[m] = quanta;

        Necronomata automata = new Necronomata(angles, extent, frequency);

        int steps = 1500;
        double twoPi = 2.0 * Math.PI;
        for (int s = 0; s < steps; s++) {
            automata.step();
            assertTrue("angle must never go negative, was " + angles[m],
                       angles[m] >= 0f);
            assertTrue("angle must stay within [0, 2*PI) after every step, was "
                       + angles[m], angles[m] < (float) twoPi);
        }

        double expected = expectedWrappedAngle(steps, quanta, twoPi);
        assertEquals("wrapped angle must match the mathematically expected angle mod 2*PI (floor-mod, non-negative)",
                     expected, angles[m], 5e-3);
    }

    private static double expectedWrappedAngle(int steps, int quanta, double twoPi) {
        double raw = steps * (double) quanta * Necronomata.QUANTUM_RATE;
        return raw - twoPi * Math.floor(raw / twoPi);
    }

    /**
     * inviscid-5sk: deltaA is DERIVED, never independent - step()
     * unconditionally recomputes it as {@code QUANTUM_RATE * frequency}
     * every tick, so a Processor that writes deltaA directly (outside the
     * documented contract) is self-healing: the stray value does not
     * survive the next step(). This is the weaker of the two raw-array
     * writes the escape hatch permits; see
     * {@link #processorWritingAngleIsVisibleDynamics()} for the one that
     * is NOT self-healing.
     */
    @Test
    public void processorWritingDeltaADirectlyIsOverwrittenByNextStep() {
        Point3i extent = new Point3i(1, 1, 1);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];
        int m = 6;
        frequency[m] = 2f;

        Necronomata automata = new Necronomata(angles, extent, frequency);

        automata.process((angle, freq, deltaA, deltaF) -> deltaA[m] = 999f);
        automata.process((angle, freq, deltaA, deltaF) -> assertEquals(
                          "stray deltaA write must be visible until the next step()",
                          999f, deltaA[m], 0f));

        automata.step();

        automata.process((angle, freq, deltaA, deltaF) -> assertEquals(
                          "deltaA must be fully recomputed from frequency by step(), "
                          + "overwriting any stray write - it is self-healing",
                          Necronomata.QUANTUM_RATE * freq[m], deltaA[m], 1e-6f));
    }

    /**
     * inviscid-5sk negative control: unlike deltaA, a Processor that
     * writes angle directly is NOT self-healed by step() - step() reads
     * angle and adds deltaA to it, so a stray angle write permanently
     * perturbs the trajectory. This test documents that accepted risk of
     * the raw-array escape hatch (Necronomata.process(Processor)); it
     * deliberately asserts the dynamics-changing behavior occurs, rather
     * than guarding against it, per the inviscid-5sk decision to accept
     * this exposure for Phase A ergonomics rather than add a runtime
     * guard that could only check the (already self-healing) deltaA path.
     */
    @Test
    public void processorWritingAngleIsVisibleDynamics() {
        Point3i extent = new Point3i(1, 1, 1);
        int len = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[len];
        float[] frequency = new float[len];
        int m = 1;
        frequency[m] = 0f;

        Necronomata automata = new Necronomata(angles, extent, frequency);

        automata.process((angle, freq, deltaA, deltaF) -> angle[m] = 1.23f);

        assertEquals("a Processor writing angle directly changes it immediately - "
                     + "outside the sanctioned deltaF/frequency contract, and not "
                     + "self-healed the way a stray deltaA write is",
                     1.23f, angles[m], 0f);

        automata.step();

        assertEquals("with zero frequency the member should not have moved on its "
                     + "own, but the stray angle write from before step() persists "
                     + "as the new baseline - proof the write permanently perturbed "
                     + "the trajectory rather than being reset",
                     1.23f, angles[m], 1e-6f);
    }
}
