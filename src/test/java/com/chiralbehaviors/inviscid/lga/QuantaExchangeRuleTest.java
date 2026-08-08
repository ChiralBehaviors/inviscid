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
import static org.junit.Assert.assertNotEquals;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.lga.CollisionRule.Delta;

/**
 * Behavioral tests for {@link QuantaExchangeRule} (bead inviscid-0nx.14,
 * v1 collision rule, USER-decided 2026-08-08).
 *
 * @author halhildebrand
 */
public class QuantaExchangeRuleTest {

    private static final long  QUANTA_RANGE = 8L;
    private static final Point3i CELL_A     = new Point3i(0, 0, 0);
    private static final Point3i CELL_B     = new Point3i(1, 0, -1);

    private final QuantaExchangeRule rule = new QuantaExchangeRule();

    private static Contact fixtureContact() {
        return new Contact(CELL_A, 3, 1, CELL_B, 3, 0, 1, 0.0);
    }

    private static Contact swapped(Contact contact) {
        return new Contact(contact.cellB(), contact.cubeB(),
                            contact.memberB(), contact.cellA(),
                            contact.cubeA(), contact.memberA(),
                            FccNeighborhood.opposite(contact.direction()),
                            contact.minDistance());
    }

    @Test
    public void deltasSumToZeroForEveryInputPair() {
        Contact contact = fixtureContact();
        for (long quantaA = -QUANTA_RANGE; quantaA <= QUANTA_RANGE; quantaA++) {
            for (long quantaB = -QUANTA_RANGE; quantaB <= QUANTA_RANGE; quantaB++) {
                Delta delta = rule.resolve(contact, 0f, quantaA, 0f, quantaB);
                assertEquals("deltas must sum to exactly zero for quantaA="
                             + quantaA + " quantaB=" + quantaB, 0L,
                             delta.deltaA() + delta.deltaB());
            }
        }
    }

    @Test
    public void ruleIsDeterministic() {
        Contact contact = fixtureContact();
        for (long quantaA = -QUANTA_RANGE; quantaA <= QUANTA_RANGE; quantaA++) {
            for (long quantaB = -QUANTA_RANGE; quantaB <= QUANTA_RANGE; quantaB++) {
                Delta first = rule.resolve(contact, 0.3f, quantaA, 1.1f,
                                            quantaB);
                Delta second = rule.resolve(contact, 0.3f, quantaA, 1.1f,
                                             quantaB);
                assertEquals(first, second);
            }
        }
    }

    /**
     * Swapping which member is "A" (and mirroring the contact's own A/B
     * fields, since a real caller would present the endpoints consistently)
     * must produce mirrored deltas. Asymmetry here would inject a
     * systematic drift indistinguishable from real transport.
     */
    @Test
    public void ruleIsSymmetricUnderParticipantSwap() {
        Contact contact = fixtureContact();
        Contact mirrored = swapped(contact);
        for (long quantaA = -QUANTA_RANGE; quantaA <= QUANTA_RANGE; quantaA++) {
            for (long quantaB = -QUANTA_RANGE; quantaB <= QUANTA_RANGE; quantaB++) {
                Delta ab = rule.resolve(contact, 0.5f, quantaA, 1.5f, quantaB);
                Delta ba = rule.resolve(mirrored, 1.5f, quantaB, 0.5f,
                                         quantaA);
                assertEquals("swapping which member is A must mirror deltaA<->deltaB (quantaA="
                             + quantaA + " quantaB=" + quantaB + ")",
                             ab.deltaA(), ba.deltaB());
                assertEquals("swapping which member is A must mirror deltaA<->deltaB (quantaA="
                             + quantaA + " quantaB=" + quantaB + ")",
                             ab.deltaB(), ba.deltaA());
            }
        }
    }

    /**
     * Non-vacuity: an identity rule conserves everything perfectly and
     * does nothing, which would trivially pass every other test in this
     * class.
     */
    @Test
    public void ruleIsNotTheIdentity() {
        Contact contact = fixtureContact();
        Delta delta = rule.resolve(contact, 0f, 5L, 0f, 2L);
        assertNotEquals("expected a nonzero transfer when quanta differ", 0L,
                         delta.deltaA());
        assertNotEquals(0L, delta.deltaB());
    }

    /**
     * Necronomata's frequency slots are float32 (30 floats/cell); the
     * largest exactly-representable integer in a 32-bit float is 2^24. A
     * transfer near that boundary must not silently push either resulting
     * value off the exact-integer grid, and - separately - the rule must
     * not overflow computing a decision even at long extremes (it never
     * subtracts its two inputs; see the class Javadoc).
     */
    @Test
    public void quantaStayWithinRepresentableRange() {
        Contact contact = fixtureContact();

        long floatExactLimit = 1L << 24;
        Delta near = rule.resolve(contact, 0f, floatExactLimit, 0f,
                                   floatExactLimit - 2L);
        assertEquals(-1L, near.deltaA());
        assertEquals(1L, near.deltaB());
        long resultA = floatExactLimit + near.deltaA();
        long resultB = (floatExactLimit - 2L) + near.deltaB();
        assertEquals("post-transfer quanta must remain exactly representable as float32",
                     resultA, (long) (float) resultA);
        assertEquals("post-transfer quanta must remain exactly representable as float32",
                     resultB, (long) (float) resultB);

        Delta extreme = rule.resolve(contact, 0f, Long.MAX_VALUE, 0f,
                                      Long.MIN_VALUE);
        assertEquals(-1L, extreme.deltaA());
        assertEquals(1L, extreme.deltaB());

        Delta extremeMirrored = rule.resolve(contact, 0f, Long.MIN_VALUE, 0f,
                                              Long.MAX_VALUE);
        assertEquals(1L, extremeMirrored.deltaA());
        assertEquals(-1L, extremeMirrored.deltaB());
    }

    @Test
    public void equalQuantaIsANoOp() {
        Contact contact = fixtureContact();
        Delta delta = rule.resolve(contact, 0f, 3L, 0f, 3L);
        assertEquals(0L, delta.deltaA());
        assertEquals(0L, delta.deltaB());
    }
}
