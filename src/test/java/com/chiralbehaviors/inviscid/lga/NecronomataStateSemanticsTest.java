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
     * with zero frequency everywhere, stepping must never move any angle,
     * no matter what the initial angles were.
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

        Necronomata automata = new Necronomata(angles, extent, frequency);
        automata.step();

        assertArrayEquals("angles must be unchanged when frequency is zero everywhere", initialAngles, angles,
                           1e-6f);
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
}
