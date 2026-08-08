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

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.Necronomata;

/**
 * Deterministic per-tick enumeration of every member-member contact in the
 * lattice (inviscid-0nx.13): for each even-parity cell, for each of the 6
 * "positive" FCC directions, for every {@code (cube, member)} pair on the
 * near side against every {@code (cube, member)} pair on the far side,
 * evaluate {@link ContactPredicate#contacts}.
 *
 * <h2>Ordering is load-bearing</h2>
 * Rule application order must be fixed and reproducible or the automaton is
 * not deterministic (bead spec). {@link #scan()} therefore never uses a
 * hash-ordered collection anywhere in its iteration path - only
 * {@link Necronomata#forEach} (row-major {@code i,j,k}, even-parity cells
 * only), a fixed direction list, and plain nested {@code int} loops over
 * {@code cube}/{@code member}. Two scans of the same lattice state always
 * produce bitwise-identical {@link Contact} lists in the same order.
 *
 * <h2>Each unordered cell pair visited exactly once</h2>
 * {@link FccNeighborhood#DIRECTIONS} lists the 12 directions as 6
 * "positive" indices ({@code +1..+6}) followed by their 6 negations. This
 * class restricts iteration to that positive prefix ({@link
 * #CANONICAL_DIRECTIONS}) - the canonicalization the bead spec calls for -
 * so a pair {@code (cellA, cellB)} is only ever visited from whichever side
 * reaches the other via a positive-index offset, never additionally from
 * the far side via the corresponding negative direction. This rules out
 * double-counting by construction (no contact is ever discovered twice and
 * then deduplicated after the fact); double-counting a contact would
 * double-apply its downstream collision rule and silently violate
 * conservation (see {@code ConservationAudit}).
 * <p>
 * <b>Precondition, enforced by the {@code neighborhood} constructor
 * argument, not here.</b> "Exactly once" depends on the 6 canonical
 * directions actually reaching 6 distinct neighbor cells. {@link
 * FccNeighborhood}'s constructor rejects any extent axis {@code < 4} for
 * exactly this reason: at axis extent 2, direction pairs that differ by
 * exactly 2 on that axis ({@code (+1,+3)}, {@code (+2,+5)}, {@code
 * (+4,+6)}) alias to the same wrapped neighbor cell, which both
 * duplicates a {@code Contact} (reached via two different canonical
 * directions from the same {@code cellA}) and corrupts the losing
 * direction's geometry ({@code ContactPredicate.physicalOffset} computes
 * the raw, un-wrapped displacement for whichever direction is asked, so
 * an aliased-but-different direction pairs the right wrapped cell with
 * the wrong physical offset). See {@code FccNeighborhoodTest}'s
 * axis-2-rejection tests and bead inviscid-cb7 (caught by the stacked
 * code-review-expert / substantive-critic gate on this class, both
 * independently, before this precondition existed). This class does not
 * re-validate the extent itself - any {@code FccNeighborhood} instance
 * that could be constructed already satisfies it.
 *
 * <h2>Member-direction correspondence: UNVERIFIED branch taken</h2>
 * {@link FccNeighborhood}'s Javadoc records the correspondence between
 * direction index and member index as UNVERIFIED (inviscid-0nx.3). This
 * class therefore makes NO such assumption: for every scanned cell pair it
 * evaluates the full {@code 5*6 x 5*6 == 900} {@code (cubeA, memberA)} x
 * {@code (cubeB, memberB)} combinations, not a diagonal subset keyed off
 * the direction index. If a future bead verifies the correspondence, this
 * is the place to narrow the sweep - and {@code
 * ContactScanTest.memberDirectionCorrespondenceHolds} the test that would
 * need to change with it.
 *
 * <h2>Cost, measured at extent 6^3</h2>
 * At extent {@code (6,6,6)} there are 108 even-parity cells; with 6
 * canonical directions and 900 member-pair combinations per direction,
 * that is {@code 108 * 6 * 900 == 583,200} {@link ContactPredicate#contacts}
 * evaluations per tick - confirmed by direct instrumentation (a counting
 * {@code ContactPredicate} subclass), not just arithmetic. Measured on the
 * development machine, {@code Random(42L)}-seeded angles, JIT-warmed (3
 * untimed warmup calls, best-of-5 timed thereafter): ~68ms wall time per
 * {@link #scan()} call, finding 17 contacts on that seed (consistent with
 * {@code ContactPredicateTest}'s independently-measured ~1.13e-4 contact
 * fraction: {@code 583,200 * 1.13e-4 ≈ 66} order of magnitude, and the
 * specific count varies with the angle seed since contact is a sparse,
 * angle-dependent predicate). This is the number Phase C (bead
 * inviscid-0nx.16 / the LUT precomputation) should judge its speedup
 * against.
 *
 * @author halhildebrand
 */
public class ContactScan {

    private static final int CUBES_PER_CELL   = 5;
    private static final int MEMBERS_PER_CUBE = 6;

    /**
     * The 6 "positive" FCC directions - {@link FccNeighborhood#DIRECTIONS}'
     * own prefix, not a separately hardcoded list that could drift out of
     * sync with it. See the class Javadoc's "each unordered cell pair
     * visited exactly once" section.
     */
    private static final List<Integer> CANONICAL_DIRECTIONS = FccNeighborhood.DIRECTIONS.subList(0,
                                                                                                    6);

    private final Necronomata      automaton;
    private final FccNeighborhood  neighborhood;
    private final ContactPredicate predicate;

    public ContactScan(Necronomata automaton, FccNeighborhood neighborhood,
                        ContactPredicate predicate) {
        this.automaton = automaton;
        this.neighborhood = neighborhood;
        this.predicate = predicate;
    }

    /**
     * @return every contact present in the automaton's current state, in a
     *         fixed, reproducible order (see class Javadoc).
     */
    public List<Contact> scan() {
        List<Contact> contacts = new ArrayList<>();
        automaton.forEach(cellA -> scanCell(cellA, contacts));
        return contacts;
    }

    private void scanCell(Point3i cellA, List<Contact> contacts) {
        float[] anglesA = automaton.anglesOf(cellA);
        for (int direction : CANONICAL_DIRECTIONS) {
            Point3i cellB = neighborhood.neighbor(cellA, direction);
            float[] anglesB = automaton.anglesOf(cellB);
            scanDirection(cellA, anglesA, cellB, anglesB, direction, contacts);
        }
    }

    private void scanDirection(Point3i cellA, float[] anglesA, Point3i cellB,
                                float[] anglesB, int direction,
                                List<Contact> contacts) {
        for (int cubeA = 0; cubeA < CUBES_PER_CELL; cubeA++) {
            for (int memberA = 0; memberA < MEMBERS_PER_CUBE; memberA++) {
                float angleA = anglesA[cubeA * MEMBERS_PER_CUBE + memberA];
                for (int cubeB = 0; cubeB < CUBES_PER_CELL; cubeB++) {
                    for (int memberB = 0; memberB < MEMBERS_PER_CUBE; memberB++) {
                        float angleB = anglesB[cubeB * MEMBERS_PER_CUBE
                                               + memberB];
                        if (predicate.contacts(cubeA, memberA, angleA, cubeB,
                                               memberB, angleB, direction)) {
                            double minDistance = predicate.minDistance(cubeA,
                                                                       memberA,
                                                                       angleA,
                                                                       cubeB,
                                                                       memberB,
                                                                       angleB,
                                                                       direction);
                            contacts.add(new Contact(cellA, cubeA, memberA,
                                                     cellB, cubeB, memberB,
                                                     direction, minDistance));
                        }
                    }
                }
            }
        }
    }
}
