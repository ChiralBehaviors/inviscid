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

import javax.vecmath.Point3i;

/**
 * One member-member contact found by {@link ContactScan}: member
 * {@code (cubeA, memberA)} of {@code cellA} touches member {@code (cubeB,
 * memberB)} of {@code cellB}, which is {@code cellA}'s neighbor in
 * {@code direction} (per {@link FccNeighborhood#neighbor(Point3i, int)}).
 * {@code minDistance} is the separation {@link ContactPredicate#minDistance}
 * computed for this pair - carried here rather than recomputed by every
 * consumer, per that method's own Javadoc rationale (logging contact
 * margins, not just fire/no-fire).
 *
 * <p>{@code direction} is always one of {@link FccNeighborhood}'s 6
 * "positive" directions ({@code +1..+6}) - {@link ContactScan} visits each
 * unordered cell pair exactly once by restricting to that half of the
 * direction set (see that class's Javadoc), so a {@code Contact} never
 * carries a negative direction and the pair {@code (cellB, cubeB, memberB)}
 * to {@code (cellA, cubeA, memberA)} via {@code opposite(direction)} is
 * never separately reported.
 *
 * <p>{@code cellA} and {@code cellB} are defensively copied in the
 * canonical constructor - {@link javax.vecmath.Point3i} is mutable, and a
 * caller must not be able to corrupt an already-recorded contact by
 * mutating a {@code Point3i} it happens to still hold a reference to.
 *
 * @author halhildebrand
 */
public record Contact(Point3i cellA, int cubeA, int memberA, Point3i cellB,
                       int cubeB, int memberB, int direction,
                       double minDistance) {

    public Contact {
        cellA = new Point3i(cellA);
        cellB = new Point3i(cellB);
    }

    /**
     * @return a defensive copy of {@code cellA}.
     */
    @Override
    public Point3i cellA() {
        return new Point3i(cellA);
    }

    /**
     * @return a defensive copy of {@code cellB}.
     */
    @Override
    public Point3i cellB() {
        return new Point3i(cellB);
    }
}
