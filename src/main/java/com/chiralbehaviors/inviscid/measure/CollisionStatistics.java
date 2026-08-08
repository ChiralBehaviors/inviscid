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

package com.chiralbehaviors.inviscid.measure;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.lga.FccNeighborhood;

/**
 * Collision-recording harness for the (not-yet-written) discrete
 * conservation-exact collision rules (bead inviscid-0nx.14 / .20). Callers
 * that iterate contacts and apply discrete quanta transfers call
 * {@link #recordCollision(Point3i, int, Point3i, int, int, long, int)} once
 * per resolved collision; this class only accumulates statistics, it never
 * decides whether a collision occurs or what it transfers.
 *
 * <p>Member&harr;direction correspondence is UNVERIFIED (see
 * {@link FccNeighborhood}'s class javadoc) -- direction is recorded exactly
 * as the caller reports it, never derived from a member index here.
 *
 * <p>All internal accumulators use deterministic-iteration-order
 * structures ({@link TreeMap}, a fixed-order {@link LinkedHashMap} seeded
 * from {@link FccNeighborhood#DIRECTIONS}, and plain arrays) -- never a
 * bare {@link java.util.HashMap} -- so repeated runs over the same inputs
 * report statistics in the same order.
 *
 * @author halhildebrand
 */
public class CollisionStatistics {

    private final Map<Integer, Long> perDirection;
    private final Map<Integer, Long> perTick = new TreeMap<>();
    private final Map<Long, Long> magnitudeHistogram = new TreeMap<>();
    private final Map<Long, Long> perMemberPair = new TreeMap<>();

    private long totalCollisions = 0L;
    private int minTick = Integer.MAX_VALUE;
    private int maxTick = Integer.MIN_VALUE;

    public CollisionStatistics() {
        Map<Integer, Long> directions = new LinkedHashMap<>();
        for (int direction : FccNeighborhood.DIRECTIONS) {
            directions.put(direction, 0L);
        }
        this.perDirection = directions;
    }

    /**
     * @return the total number of collisions across every direction; must
     *         equal the sum of {@link #collisionsPerDirection()}'s values
     *         by construction (both are updated together).
     */
    public long totalCollisions() {
        return totalCollisions;
    }

    /**
     * @return the per-tick collision count, keyed by tick, in ascending
     *         tick order.
     */
    public Map<Integer, Long> collisionsPerTick() {
        return Collections.unmodifiableMap(perTick);
    }

    /**
     * @return the per-direction collision count, one entry for each of the
     *         12 {@link FccNeighborhood#DIRECTIONS}, in that fixed order
     *         (present with a zero count even if never observed).
     */
    public Map<Integer, Long> collisionsPerDirection() {
        return Collections.unmodifiableMap(perDirection);
    }

    public long collisionsInDirection(int direction) {
        Long count = perDirection.get(direction);
        if (count == null) {
            throw new IllegalArgumentException("Unknown FCC direction: "
                                                + direction
                                                + ", must be one of +/-1..+/-6");
        }
        return count;
    }

    /**
     * @return the collision count for the unordered member-index pair
     *         {@code (memberA, memberB)}; canonicalized so
     *         {@code collisionsForMemberPair(a, b) ==
     *         collisionsForMemberPair(b, a)}.
     */
    public long collisionsForMemberPair(int memberA, int memberB) {
        Long count = perMemberPair.get(memberPairKey(memberA, memberB));
        return count == null ? 0L : count;
    }

    /**
     * @return the histogram of transfer magnitudes, keyed by magnitude, in
     *         ascending magnitude order.
     */
    public Map<Long, Long> transferMagnitudeHistogram() {
        return Collections.unmodifiableMap(magnitudeHistogram);
    }

    /**
     * A cheap proxy for mean free path: the tick-span covered by recorded
     * collisions divided by the number of collisions. Not a physical mean
     * free path (no spatial distance is tracked here) -- a coarse "how
     * often does a collision happen" signal only.
     *
     * @return the proxy value, or {@link Double#NaN} if no collisions have
     *         been recorded.
     */
    public double meanFreePathProxy() {
        if (totalCollisions == 0L) {
            return Double.NaN;
        }
        long span = (long) maxTick - (long) minTick + 1L;
        return (double) span / (double) totalCollisions;
    }

    /**
     * Records one resolved collision. Caller-facing recording API for
     * A.3/A.4's collision rules: they iterate contacts, apply a discrete
     * quanta transfer, and call this once per resolved contact.
     *
     * @param cellA             the first cell involved
     * @param memberA           the member index (within {@code cellA})
     *                          involved
     * @param cellB             the second cell involved
     * @param memberB           the member index (within {@code cellB})
     *                          involved
     * @param direction         the FCC direction from {@code cellA} to
     *                          {@code cellB}, one of
     *                          {@link FccNeighborhood#DIRECTIONS}
     * @param transferMagnitude the magnitude of the quanta transfer this
     *                          collision performed; MUST be non-negative
     *                          (the sign/direction of the transfer is a
     *                          collision-rule concern, not a statistics
     *                          concern) -- a negative value is rejected
     *                          rather than silently recorded
     * @param tick              the tick this collision occurred on
     * @throws IllegalArgumentException if {@code direction} is not one of
     *                                   {@link FccNeighborhood#DIRECTIONS}
     *                                   or {@code transferMagnitude} is
     *                                   negative
     */
    public void recordCollision(Point3i cellA, int memberA, Point3i cellB,
                                 int memberB, int direction,
                                 long transferMagnitude, int tick) {
        if (!perDirection.containsKey(direction)) {
            throw new IllegalArgumentException("Unknown FCC direction: "
                                                + direction
                                                + ", must be one of +/-1..+/-6");
        }
        if (transferMagnitude < 0) {
            throw new IllegalArgumentException("transferMagnitude must be non-negative, was: "
                                                + transferMagnitude);
        }
        totalCollisions++;
        perDirection.merge(direction, 1L, Long::sum);
        perTick.merge(tick, 1L, Long::sum);
        perMemberPair.merge(memberPairKey(memberA, memberB), 1L, Long::sum);
        magnitudeHistogram.merge(transferMagnitude, 1L, Long::sum);
        if (tick < minTick) {
            minTick = tick;
        }
        if (tick > maxTick) {
            maxTick = tick;
        }
    }

    private static long memberPairKey(int memberA, int memberB) {
        int lo = Math.min(memberA, memberB);
        int hi = Math.max(memberA, memberB);
        return ((long) lo << 32) | (hi & 0xFFFFFFFFL);
    }
}
