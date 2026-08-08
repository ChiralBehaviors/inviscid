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

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Point3i;

/**
 * The 12 nearest-neighbor offsets of the FCC even-parity sublattice used by
 * {@code Necronomata}, plus wrap-aware neighbor lookup. Cells live on the
 * cubic lattice but only the even-parity sublattice ({@code (i+j+k) % 2 ==
 * 0}) is populated (see {@code Necronomata#forEach} /
 * {@code Necronomata#iterator}); every one of these 12 offsets preserves
 * {@code (i+j+k)} parity, so applying one from an even-parity cell lands on
 * another even-parity cell, unbounded (no wrap).
 *
 * <p>Originally transcribed verbatim from the offset comment in
 * {@code Necronomata#process(Point3i)}, direction indices {@code +1..+6}
 * and {@code -1..-6}:
 * <pre>
 * [+1] = { i+1, j  , k-1 }
 * [-1] = { i-1, j  , k+1 }
 * [+2] = { i  , j-1, k+1 }
 * [-2] = { i  , j+1, k-1 }
 * [+3] = { i+1, j  , k+1 }
 * [-3] = { i-1, j  , k-1 }
 * [+4] = { i+1, j+1, k   }
 * [-4] = { i-1, j-1, k   }
 * [+5] = { i  , j+1, k+1 }
 * [-5] = { i  , j-1, k-1 }
 * [+6] = { i-1, j+1, k   }
 * [-6] = { i+1, j-1, k   }
 * </pre>
 *
 * <p><b>Periodic wrap and parity closure.</b> Under periodic boundary
 * conditions on an axis of extent {@code X}, an offset of {@code +1} on
 * that axis wraps to an effective delta of {@code (+1 - X)}, which shifts
 * {@code (i+j+k) mod 2} by {@code X mod 2}. The even-parity sublattice is
 * closed under wrap iff every extent axis is even ({@code X mod 2 == 0}),
 * since subtracting an even number never changes parity. This class
 * enforces that precondition at construction: an odd extent axis throws
 * {@link IllegalArgumentException} naming the offending axis rather than
 * silently manufacturing spurious anisotropy in isotropy measurements.
 *
 * <p><b>Minimum extent 4 per axis.</b> Even alone is not enough: at
 * extent exactly 2 on some axis, two of the canonical direction offsets
 * that differ by exactly 2 in that one axis (the pairs {@code (+1,+3)},
 * {@code (+2,+5)}, {@code (+4,+6)} - see the offset table above) wrap
 * {@code mod 2} to the identical delta on that axis, so both directions
 * land on the SAME neighbor cell instead of two distinct ones. That
 * breaks this class's own "12 distinct neighbors" guarantee (see {@code
 * FccNeighborhoodTest.twelveDistinctOffsetsAllParityPreserving} - true of
 * the offset table in isolation, but not of the wrapped neighbor set at
 * extent 2) for every consumer, not merely {@code ContactScan}: a
 * consumer relying on {@code neighbor(cell, d1) != neighbor(cell, d2)}
 * for {@code d1 != d2} silently gets duplicate/aliased neighbors instead.
 * (Discovered via {@code ContactScan}'s double-count guard firing at
 * extent {@code (4,2,4)}: the y-axis-2 aliasing produced duplicate
 * {@code Contact} entries reached via two different canonical directions,
 * and the losing direction's contact geometry was independently wrong -
 * {@code ContactPredicate.physicalOffset} computes the RAW un-wrapped
 * offset for whichever direction is asked, so an aliased direction is
 * evaluated against the correct wrapped cell but the wrong physical
 * displacement. Bead inviscid-cb7.) The constructor therefore rejects any
 * axis {@code < 4}, not just odd or non-positive axes.
 *
 * <p><b>UNVERIFIED:</b> the correspondence between these 12 direction
 * indices and per-cell member indices (which of the 30 floats-per-cell in
 * {@code Necronomata} a given direction's collision partner reads/writes)
 * is NOT established here. {@code MemberGeometry} (bead .2) does not exist
 * yet; this class only carries lattice-site adjacency, not member
 * geometry. Do not assume direction index {@code d} maps to member index
 * {@code d} (or any other specific mapping) without that future
 * verification.
 *
 * @author halhildebrand
 */
public class FccNeighborhood {

    private static final Map<Integer, Point3i> OFFSET_BY_DIRECTION;

    static {
        Map<Integer, Point3i> offsets = new LinkedHashMap<>();
        offsets.put(1, new Point3i(1, 0, -1));
        offsets.put(-1, new Point3i(-1, 0, 1));
        offsets.put(2, new Point3i(0, -1, 1));
        offsets.put(-2, new Point3i(0, 1, -1));
        offsets.put(3, new Point3i(1, 0, 1));
        offsets.put(-3, new Point3i(-1, 0, -1));
        offsets.put(4, new Point3i(1, 1, 0));
        offsets.put(-4, new Point3i(-1, -1, 0));
        offsets.put(5, new Point3i(0, 1, 1));
        offsets.put(-5, new Point3i(0, -1, -1));
        offsets.put(6, new Point3i(-1, 1, 0));
        offsets.put(-6, new Point3i(1, -1, 0));
        OFFSET_BY_DIRECTION = Collections.unmodifiableMap(offsets);
    }

    /**
     * All 12 direction indices, {@code +1..+6} and {@code -1..-6}.
     * Immutable ({@link List#of}) — unlike an array, a caller cannot
     * corrupt shared traversal order in place (e.g. a contact scan that
     * canonicalizes iteration order and is tempted to reorder in place).
     */
    public static final List<Integer> DIRECTIONS = List.of(1, 2, 3, 4, 5, 6,
                                                             -1, -2, -3, -4,
                                                             -5, -6);

    /**
     * The opposite of {@code direction}: {@code -direction}, validated
     * against the known direction set. Not a lookup table — every offset
     * table entry happens to be the literal negation of its opposite's
     * entry, which is why negation alone suffices here.
     */
    public static int opposite(int direction) {
        validateDirection(direction);
        return -direction;
    }

    /**
     * The lattice offset for {@code direction}, as a defensive copy.
     */
    public static Point3i offsetOf(int direction) {
        validateDirection(direction);
        return new Point3i(OFFSET_BY_DIRECTION.get(direction));
    }

    private static void validateDirection(int direction) {
        if (!OFFSET_BY_DIRECTION.containsKey(direction)) {
            throw new IllegalArgumentException("Unknown FCC direction: "
                                                + direction
                                                + ", must be one of +/-1..+/-6");
        }
    }

    private static void requireEvenAxisAtLeastFour(int value, String axis) {
        if (value <= 0) {
            throw new IllegalArgumentException("extent." + axis
                                                + " must be positive, was: "
                                                + value);
        }
        if (value % 2 != 0) {
            throw new IllegalArgumentException("extent." + axis
                                                + " must be even for the FCC even-parity sublattice to be closed under periodic wrap, was: "
                                                + value);
        }
        if (value < 4) {
            throw new IllegalArgumentException("extent." + axis
                                                + " must be at least 4: at extent 2, canonical FCC direction pairs that differ by exactly 2 in this axis (e.g. +1/+3, +2/+5, +4/+6) alias to the SAME wrapped neighbor cell, breaking the 12-distinct-neighbors guarantee this class exists to provide; was: "
                                                + value);
        }
    }

    private final Point3i extent;

    /**
     * @param extent the periodic-wrap extent; every axis must be positive,
     *               even, and at least 4.
     * @throws IllegalArgumentException if any axis of {@code extent} is
     *                                  zero, negative, odd, or 2; the
     *                                  message names the offending axis
     *                                  (and, for the axis-2 case, the
     *                                  aliasing reason - see class
     *                                  Javadoc).
     */
    public FccNeighborhood(Point3i extent) {
        requireEvenAxisAtLeastFour(extent.x, "x");
        requireEvenAxisAtLeastFour(extent.y, "y");
        requireEvenAxisAtLeastFour(extent.z, "z");
        this.extent = new Point3i(extent);
    }

    /**
     * @return a defensive copy of the periodic-wrap extent.
     */
    public Point3i getExtent() {
        return new Point3i(extent);
    }

    /**
     * The neighbor of {@code cell} in {@code direction}, wrapped
     * periodically on every axis.
     */
    public Point3i neighbor(Point3i cell, int direction) {
        Point3i offset = offsetOf(direction);
        int x = Math.floorMod(cell.x + offset.x, extent.x);
        int y = Math.floorMod(cell.y + offset.y, extent.y);
        int z = Math.floorMod(cell.z + offset.z, extent.z);
        return new Point3i(x, y, z);
    }
}
