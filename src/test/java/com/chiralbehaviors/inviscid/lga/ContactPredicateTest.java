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
import static org.junit.Assert.assertTrue;

import java.util.Random;

import org.junit.Test;

/**
 * Behavioral tests for {@link ContactPredicate} (inviscid-0nx.12).
 * <p>
 * The fixtures used by {@link #contactSetIsNonEmptyOverTheAngleSweep()}
 * and {@link #contactSetIsNotUniversal()} are not guessed - they were
 * found by an exhaustive-search probe over every {@code (direction,
 * cubeA, memberA, cubeB, memberB)} combination at {@code RESOLUTION}
 * against the real {@code MemberGeometry}/{@code FccNeighborhood}
 * geometry (radius 0.015, matching {@code NecronomataVisualization}'s
 * call site), reported back with the bead hand-off:
 * <ul>
 * <li>Full coarse sweep (24 angle steps, all 900 member-pair combos, all
 * 12 directions - 6,220,800 samples): min distance 0.0, max distance
 * ~12.60, mean ~5.87, contact fraction ~1.13e-4.</li>
 * <li>Fine sweep restricted to the combo used below (360x360 = 129,600
 * samples): min 0.0 (at angleA=angleB&asymp;1.7628 rad), max exactly
 * {@code Cubes[0].getEdgeLength()} (~5.236, at angleA=angleB=0), mean
 * ~2.39, contact fraction ~0.236% (306 of 129,600 pairs).</li>
 * </ul>
 * Both the never-fires and always-fires degenerate outcomes are
 * therefore ruled out empirically, not assumed.
 *
 * @author halhildebrand
 */
public class ContactPredicateTest {

    private static final double  RADIUS       = LgaTestGeometry.BASELINE_RADIUS;
    private static final int     RESOLUTION   = 360;

    // The empirically-found contact-bearing combo (see class Javadoc).
    private static final int     CUBE_A       = 3;
    private static final int     MEMBER_A     = 1;
    private static final int     CUBE_B       = 3;
    private static final int     MEMBER_B     = 0;
    private static final int     DIRECTION    = 1;

    private static float angleOf(int step) {
        return (float) (step * 2 * Math.PI / RESOLUTION);
    }

    private static ContactPredicate newPredicate() {
        return new ContactPredicate(new MemberGeometry(RESOLUTION, RADIUS));
    }

    /**
     * A member never contacts itself, or a member of the same physical
     * strut, evaluated against itself across every one of the 12 FCC
     * directions - interaction is strictly inter-cell (locked design).
     * If {@link ContactPredicate#minDistance} ever forgot to apply the
     * inter-cell offset, this would fail immediately (identical segments
     * compared to themselves give distance 0, well inside contact
     * range).
     */
    @Test
    public void selfContactIsNeverReported() {
        ContactPredicate predicate = newPredicate();
        for (int cube = 0; cube < 5; cube++) {
            for (int member = 0; member < 6; member++) {
                for (int step = 0; step < RESOLUTION; step += 7) {
                    float angle = angleOf(step);
                    for (int direction : FccNeighborhood.DIRECTIONS) {
                        assertFalse("cube " + cube + " member " + member
                                   + " angle " + angle + " direction "
                                   + direction
                                   + ": self-contact must never be reported",
                                   predicate.contacts(cube, member, angle,
                                                      cube, member, angle,
                                                      direction));
                    }
                }
            }
        }
    }

    /**
     * {@code contacts(A in C, B in C+d) == contacts(B in C+d, A in C via
     * opposite(d))} - verified here against 2,000 seeded-random
     * {@code (cube, member, angle, direction)} combinations. (A prior
     * exploratory probe against 200,000 combinations found zero
     * mismatches; 2,000 is kept here to keep the suite fast while still
     * exercising the property broadly.) If this fails, the offset
     * bookkeeping in {@link ContactPredicate#minDistance} is wrong.
     */
    @Test
    public void contactIsSymmetricUnderDirectionReversal() {
        ContactPredicate predicate = newPredicate();
        Random random = new Random(42L);

        for (int i = 0; i < 2000; i++) {
            int cubeA = random.nextInt(5);
            int memberA = random.nextInt(6);
            float angleA = random.nextFloat() * (float) (2 * Math.PI);
            int cubeB = random.nextInt(5);
            int memberB = random.nextInt(6);
            float angleB = random.nextFloat() * (float) (2 * Math.PI);
            int direction = FccNeighborhood.DIRECTIONS.get(random.nextInt(FccNeighborhood.DIRECTIONS.size()));
            int opposite = FccNeighborhood.opposite(direction);

            boolean forward = predicate.contacts(cubeA, memberA, angleA,
                                                 cubeB, memberB, angleB,
                                                 direction);
            boolean reversed = predicate.contacts(cubeB, memberB, angleB,
                                                   cubeA, memberA, angleA,
                                                   opposite);
            assertEquals("iteration " + i + ": cubeA=" + cubeA + " memberA="
                        + memberA + " angleA=" + angleA + " cubeB=" + cubeB
                        + " memberB=" + memberB + " angleB=" + angleB
                        + " direction=" + direction, forward, reversed);
        }
    }

    /**
     * NON-VACUITY: sweeping both angles over the empirically-found
     * contact-bearing combo must yield at least one contact. A predicate
     * that never fires would produce an automaton with no dynamics and
     * would sail through every other test in this class.
     */
    @Test
    public void contactSetIsNonEmptyOverTheAngleSweep() {
        ContactPredicate predicate = newPredicate();
        int contactCount = 0;
        for (int ia = 0; ia < RESOLUTION; ia++) {
            float angleA = angleOf(ia);
            for (int ib = 0; ib < RESOLUTION; ib++) {
                float angleB = angleOf(ib);
                if (predicate.contacts(CUBE_A, MEMBER_A, angleA, CUBE_B,
                                       MEMBER_B, angleB, DIRECTION)) {
                    contactCount++;
                }
            }
        }
        assertTrue("expected at least one contact over the " + RESOLUTION
                   + "x" + RESOLUTION + " angle sweep, found none",
                   contactCount > 0);
    }

    /**
     * Guards the opposite failure from {@link
     * #contactSetIsNonEmptyOverTheAngleSweep()}: a predicate that always
     * fires is equally useless. angle (0,0) on the same combo is a known
     * non-contact point (separation exactly {@code
     * PhiCoordinates.Cubes[0].getEdgeLength()}, ~5.236 - far outside
     * {@code 2*RADIUS}).
     */
    @Test
    public void contactSetIsNotUniversal() {
        ContactPredicate predicate = newPredicate();
        assertFalse(predicate.contacts(CUBE_A, MEMBER_A, 0f, CUBE_B, MEMBER_B,
                                       0f, DIRECTION));

        int nonContactCount = 0;
        for (int ia = 0; ia < RESOLUTION; ia++) {
            float angleA = angleOf(ia);
            for (int ib = 0; ib < RESOLUTION; ib++) {
                float angleB = angleOf(ib);
                if (!predicate.contacts(CUBE_A, MEMBER_A, angleA, CUBE_B,
                                        MEMBER_B, angleB, DIRECTION)) {
                    nonContactCount++;
                }
            }
        }
        assertTrue("expected at least one non-contact over the "
                   + RESOLUTION + "x" + RESOLUTION + " angle sweep, found none (predicate fires universally)",
                   nonContactCount > 0);
    }

    /**
     * Pure function of its arguments: repeated calls with identical
     * inputs return identical results, and interleaving unrelated calls
     * does not perturb a previously-computed result (no shared mutable
     * state).
     */
    @Test
    public void predicateIsPureAndDeterministic() {
        ContactPredicate predicate = newPredicate();

        double first = predicate.minDistance(CUBE_A, MEMBER_A, 1.7627826f,
                                             CUBE_B, MEMBER_B, 1.7627826f,
                                             DIRECTION);
        boolean firstContact = predicate.contacts(CUBE_A, MEMBER_A,
                                                   1.7627826f, CUBE_B,
                                                   MEMBER_B, 1.7627826f,
                                                   DIRECTION);

        // Interleave unrelated calls across every direction/member/cube.
        for (int cube = 0; cube < 5; cube++) {
            for (int member = 0; member < 6; member++) {
                for (int direction : FccNeighborhood.DIRECTIONS) {
                    predicate.minDistance(cube, member, 3.3f, cube, member,
                                          0.7f, direction);
                }
            }
        }

        double second = predicate.minDistance(CUBE_A, MEMBER_A, 1.7627826f,
                                              CUBE_B, MEMBER_B, 1.7627826f,
                                              DIRECTION);
        boolean secondContact = predicate.contacts(CUBE_A, MEMBER_A,
                                                    1.7627826f, CUBE_B,
                                                    MEMBER_B, 1.7627826f,
                                                    DIRECTION);

        assertEquals(first, second, 0.0);
        assertEquals(firstContact, secondContact);
    }
}
