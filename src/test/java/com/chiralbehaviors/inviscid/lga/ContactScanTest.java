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

import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;

/**
 * Behavioral tests for {@link ContactScan} (inviscid-0nx.13).
 *
 * <h2>Member-direction correspondence branch (test 6)</h2>
 * {@code FccNeighborhood}'s Javadoc records the correspondence between
 * direction index and member index as UNVERIFIED (bead inviscid-0nx.3).
 * {@link #memberDirectionCorrespondenceHolds()} therefore takes the
 * UNVERIFIED branch: it asserts {@link ContactScan} makes NO
 * correspondence assumption and evaluates the full 30x30 member x member
 * sweep for every scanned cell pair, not a diagonal subset.
 *
 * @author halhildebrand
 */
public class ContactScanTest {

    private static final double RADIUS     = 0.015;
    private static final int    RESOLUTION = 360;

    // The empirically-found contact-bearing (cube,member,angle) combo, same
    // fixture ContactPredicateTest uses (see that class's Javadoc for the
    // exhaustive-search provenance): at angleA == angleB == FIXTURE_ANGLE,
    // minDistance is exactly 0.0, well inside 2*RADIUS.
    private static final int    FIXTURE_CUBE_A    = 3;
    private static final int    FIXTURE_MEMBER_A  = 1;
    private static final int    FIXTURE_CUBE_B    = 3;
    private static final int    FIXTURE_MEMBER_B  = 0;
    private static final int    FIXTURE_DIRECTION = 1;
    private static final float  FIXTURE_ANGLE     = 1.7627826f;

    private static MemberGeometry newGeometry() {
        return new MemberGeometry(RESOLUTION, RADIUS);
    }

    private static ContactPredicate newPredicate() {
        return new ContactPredicate(newGeometry());
    }

    private static void seed(Necronomata automaton, Point3i cell, int cube,
                              int member, float angle) {
        int localIndex = cube * 6 + member;
        automaton.process((angleArray, frequency, deltaA, deltaF) -> {
            angleArray[automaton.indexOfCell(cell) + localIndex] = angle;
        });
    }

    /**
     * Same lattice state scanned twice must yield an identical ordered
     * list - no HashSet/HashMap iteration-order dependence, no
     * incidental nondeterminism from e.g. parallel streams.
     */
    @Test
    public void scanIsDeterministic() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 42L);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());

        List<Contact> first = scan.scan();
        List<Contact> second = scan.scan();

        assertEquals("repeated scans of the same state must be identical",
                     first, second);
    }

    /**
     * Canonicalizing every contact (as the unordered pair of
     * (cell,cube,member) endpoints, ignoring direction sign) must never
     * find a duplicate - double-counting a contact would double-apply its
     * collision rule and silently violate conservation.
     * <p>
     * This dedup check is meaningful only because {@link ContactScan}'s
     * own "exactly once" guarantee depends on every canonical direction
     * reaching a distinct neighbor cell. An axis extent of 2 would break
     * that (direction pairs {@code (+1,+3)}, {@code (+2,+5)}, {@code
     * (+4,+6)} alias to the same wrapped neighbor - bead inviscid-cb7),
     * but is unreachable here: {@link FccNeighborhood}'s constructor
     * rejects any axis {@code < 4} (see {@code
     * FccNeighborhoodTest.axisExtentTwoIsRejected} and that class's
     * "Minimum extent 4 per axis" Javadoc), so every {@code
     * FccNeighborhood} this test could construct already satisfies the
     * precondition this check relies on.
     */
    @Test
    public void eachCellPairVisitedExactlyOnce() {
        Point3i extent = new Point3i(6, 6, 6);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 7L);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        List<Contact> contacts = scan.scan();

        // Non-vacuity for this test's own instrumentation: random seeding
        // at this extent is known (see ContactPredicateTest's exhaustive
        // sweep: contact fraction ~1.13e-4) to produce a handful of
        // contacts, so the dedup check below is exercised against a
        // non-empty set, not vacuously true over an empty list.
        assertTrue("expected at least one contact from random seeding at 6^3 "
                   + "to make the dedup check meaningful, found none",
                   !contacts.isEmpty());

        Set<String> canonicalKeys = new HashSet<>();
        for (Contact contact : contacts) {
            String endpointA = endpointKey(contact.cellA(), contact.cubeA(),
                                            contact.memberA());
            String endpointB = endpointKey(contact.cellB(), contact.cubeB(),
                                            contact.memberB());
            // Canonical form: sort the two endpoint keys lexicographically
            // so a contact recorded as A->B is the same key as one that
            // would have been recorded as B->A.
            String canonical = (endpointA.compareTo(endpointB) <= 0)
                                ? endpointA + "|" + endpointB
                                : endpointB + "|" + endpointA;
            assertTrue("cell pair visited more than once: " + canonical,
                       canonicalKeys.add(canonical));
        }
    }

    private static String endpointKey(Point3i cell, int cube, int member) {
        return cell.x + "," + cell.y + "," + cell.z + ":" + cube + ":"
               + member;
    }

    @Test
    public void onlyEvenParityCellsAreScanned() {
        Point3i extent = new Point3i(6, 6, 6);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 99L);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        List<Contact> contacts = scan.scan();

        assertTrue("expected at least one contact to make this check meaningful",
                   !contacts.isEmpty());

        for (Contact contact : contacts) {
            assertEquals("cellA " + contact.cellA()
                         + " must be even-parity", 0,
                         Math.floorMod(contact.cellA().x + contact.cellA().y
                                       + contact.cellA().z, 2));
            assertEquals("cellB " + contact.cellB()
                         + " must be even-parity", 0,
                         Math.floorMod(contact.cellB().x + contact.cellB().y
                                       + contact.cellB().z, 2));
        }
    }

    /**
     * A contact whose only path from cellA to cellB crosses the periodic
     * boundary must still be found. cellA=(3,0,1) (even parity, sum=4) in
     * direction +1 (offset (1,0,-1)) on extent (4,4,4) has a raw target of
     * (4,0,0) - out of range on x - so the neighbor is only reachable via
     * wrap, landing on (0,0,0).
     */
    @Test
    public void scanRespectsPeriodicWrap() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);

        Point3i cellA = new Point3i(3, 0, 1);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        Point3i cellB = neighborhood.neighbor(cellA, FIXTURE_DIRECTION);
        assertEquals("wrap target must be (0,0,0) for this fixture",
                     new Point3i(0, 0, 0), cellB);

        seed(automaton, cellA, FIXTURE_CUBE_A, FIXTURE_MEMBER_A,
             FIXTURE_ANGLE);
        seed(automaton, cellB, FIXTURE_CUBE_B, FIXTURE_MEMBER_B,
             FIXTURE_ANGLE);

        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        List<Contact> contacts = scan.scan();

        boolean found = contacts.stream()
                                 .anyMatch(c -> c.cellA().equals(cellA)
                                               && c.cubeA() == FIXTURE_CUBE_A
                                               && c.memberA() == FIXTURE_MEMBER_A
                                               && c.cellB().equals(cellB)
                                               && c.cubeB() == FIXTURE_CUBE_B
                                               && c.memberB() == FIXTURE_MEMBER_B
                                               && c.direction() == FIXTURE_DIRECTION);
        assertTrue("expected the wrap-only contact between " + cellA + " and "
                   + cellB + " to be found; contacts were: " + contacts,
                   found);
    }

    /**
     * Paired non-vacuity guard (with
     * {@link #scanFindsSomethingOnASeededConfiguration()}): a lattice with
     * every member at a verified non-contact angle must scan clean.
     * <p>
     * <b>NOT the default all-zero lattice.</b> An exploratory probe run
     * while writing this test found that the default zero-angle
     * configuration is NOT contact-free: {@code cubeA=3} has several
     * member combinations (against {@code cubeB=3} of the neighboring
     * cell) whose separation at angle {@code (0,0)} is exactly {@code 0.0}
     * or within double-precision noise of it, for directions {@code 2, 4,
     * 5, 6} - a genuine rest-position geometric coincidence of cube index
     * 3 (the same cube {@code ContactPredicateTest}'s own fixture combo
     * uses), not a bug in either {@code ContactPredicate} or
     * {@code ContactScan}. Since every even-parity cell in any lattice has
     * distinct neighbors in those directions, the all-zero lattice always
     * has contacts and cannot serve as this test's empty fixture.
     * <p>
     * Uniform angle {@code 1.0f} (verified by an exhaustive probe over all
     * 6 canonical directions x 900 member-pair combinations = 5,400
     * evaluations, immediately before writing this assertion: zero
     * contacts) is used instead.
     */
    @Test
    public void scanFindsNothingOnAnEmptyConfiguration() {
        float noContactAngle = 1.0f;
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        fillUniformAngle(automaton, extent, noContactAngle);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());

        assertTrue("expected no contacts on the verified non-contact uniform-angle configuration",
                   scan.scan().isEmpty());
    }

    private static void fillUniformAngle(Necronomata automaton, Point3i extent,
                                         float angle) {
        int length = 30 * extent.x * extent.y * extent.z;
        automaton.process((angleArray, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < length; i++) {
                angleArray[i] = angle;
            }
        });
    }

    /**
     * Paired with {@link #scanFindsNothingOnAnEmptyConfiguration()}: a
     * lattice seeded with the known contact-bearing combo (same fixture
     * ContactPredicateTest.contactSetIsNonEmptyOverTheAngleSweep found
     * empirically, not guessed) must produce at least one contact.
     */
    @Test
    public void scanFindsSomethingOnASeededConfiguration() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        Point3i cellA = new Point3i(0, 0, 0);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        Point3i cellB = neighborhood.neighbor(cellA, FIXTURE_DIRECTION);

        seed(automaton, cellA, FIXTURE_CUBE_A, FIXTURE_MEMBER_A,
             FIXTURE_ANGLE);
        seed(automaton, cellB, FIXTURE_CUBE_B, FIXTURE_MEMBER_B,
             FIXTURE_ANGLE);

        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            newPredicate());
        List<Contact> contacts = scan.scan();

        assertTrue("expected at least one contact from the seeded fixture, "
                   + "found none", !contacts.isEmpty());
    }

    /**
     * THE BRANCH RULING: FccNeighborhood's direction<->member
     * correspondence is UNVERIFIED (inviscid-0nx.3), so ContactScan must
     * make no such assumption - it evaluates every one of the 30x30
     * member combinations for every scanned cell pair, not a diagonal
     * subset. Instrumented with a counting ContactPredicate subclass
     * rather than inspecting scan output, since a diagonal-only scan could
     * coincidentally still "find" the same contacts an unlucky fixture
     * would produce; only a call count distinguishes "evaluated all 900"
     * from "evaluated 30 and got lucky".
     */
    @Test
    public void memberDirectionCorrespondenceHolds() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, 13L);

        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        AtomicLong evaluationCount = new AtomicLong();
        ContactPredicate countingPredicate = new ContactPredicate(newGeometry()) {
            @Override
            public boolean contacts(int cubeA, int memberA, float angleA,
                                    int cubeB, int memberB, float angleB,
                                    int direction) {
                evaluationCount.incrementAndGet();
                return super.contacts(cubeA, memberA, angleA, cubeB, memberB,
                                      angleB, direction);
            }
        };

        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            countingPredicate);
        scan.scan();

        int evenCellCount = automaton.cellCount();
        long expected = (long) evenCellCount * 6L * 900L;
        assertEquals("expected the full 30x30 member sweep (900 combos) for "
                     + "every scanned cell pair (evenCells=" + evenCellCount
                     + " x 6 canonical directions x 900), i.e. no "
                     + "member<->direction correspondence assumption",
                     expected, evaluationCount.get());
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
}
