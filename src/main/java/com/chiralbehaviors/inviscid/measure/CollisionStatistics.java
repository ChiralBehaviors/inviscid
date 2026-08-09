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
 * {@link #recordCollision(Point3i, int, int, Point3i, int, int, int, long, int)}
 * once per resolved collision; this class only accumulates statistics, it
 * never decides whether a collision occurs or what it transfers.
 *
 * <p>Member&harr;direction correspondence is UNVERIFIED (see
 * {@link FccNeighborhood}'s class javadoc) -- direction is recorded exactly
 * as the caller reports it, never derived from a member index here.
 *
 * <p><b>Cube/member addressing (inviscid-xew).</b> Necronomata's per-cell
 * layout is {@value #CUBES_PER_CELL} cubes &times;
 * {@value #MEMBERS_PER_CUBE} members (30 floats/cell). A slot on either
 * side of a collision is therefore addressed by the pair
 * {@code (cube, member)} -- {@code cube} in
 * {@code [0, CUBES_PER_CELL)}, {@code member} in
 * {@code [0, MEMBERS_PER_CUBE)} -- exactly the addressing
 * {@link com.chiralbehaviors.inviscid.measure.ConservationAudit.Violation}
 * uses ({@code cell}, {@code cube}, {@code member} as separate axes) and
 * {@code com.chiralbehaviors.inviscid.lga.Contact} carries
 * ({@code cellA}, {@code cubeA}, {@code memberA}, {@code cellB},
 * {@code cubeB}, {@code memberB}): the first seven positional fields of a
 * {@code Contact} line up with {@code recordCollision}'s first seven
 * parameters. {@code Contact.minDistance()} has no parameter here, and
 * {@code transferMagnitude}/{@code tick} are supplied by the collision
 * rule, not read from the {@code Contact}. Member-pair keying
 * ({@link #collisionsForMemberPair(int, int, int, int)}) canonicalizes the
 * pair of {@code (cube, member)} tuples, not a bare within-cube member
 * index -- two collisions sharing a within-cube member index but occurring
 * in different cubes are distinct pairs. The coarser face-type-pair view
 * (which within-cube member indices collided, regardless of cube -- a
 * 6x6 aggregate a B.5 anisotropy consumer may want) is fully recoverable
 * by summing {@link #collisionsForMemberPair(int, int, int, int)} over
 * all cube combinations for a fixed {@code (memberA, memberB)}; nothing
 * is lost by the finer keying.
 *
 * <p>All internal accumulators use deterministic-iteration-order
 * structures ({@link TreeMap}, a fixed-order {@link LinkedHashMap} seeded
 * from {@link FccNeighborhood#DIRECTIONS}, and plain arrays) -- never a
 * bare {@link java.util.HashMap} -- so repeated runs over the same inputs
 * report statistics in the same order.
 *
 * <p><b>What counts as a "collision" here.</b> Every resolved contact
 * {@code CollisionSweep} passes to {@link #recordCollision} counts toward
 * {@link #totalCollisions()} and {@link #collisionsPerTick()} -- INCLUDING
 * contacts where the collision rule decided a no-op ({@code
 * transferMagnitude == 0}); see {@code CollisionSweep}'s class Javadoc,
 * "Recording convention", for why no-ops are recorded at all. {@link
 * #effectiveCollisions()} is the narrower, transfer-only count ({@code
 * transferMagnitude > 0}). In practice a large fraction of recorded
 * collisions are no-ops (empirically ~60% for {@code QuantaExchangeRule}
 * against typical seed densities, since a tie is a no-op and quanta
 * values cluster) -- a caller wanting a "how much actually moved" signal
 * should read {@link #effectiveCollisions()}, not {@link
 * #totalCollisions()}.
 *
 * <p>Recorded {@code direction} is always one of the 6 canonical
 * "positive" directions ({@code +1..+6}) -- the caller ({@code
 * CollisionSweep}) only ever supplies {@code Contact.direction()}, which
 * by that record's own convention never carries a negative direction.
 *
 * @author halhildebrand
 */
public class CollisionStatistics {

    /** 5 cubes per cell (see Necronomata javadoc). */
    private static final int CUBES_PER_CELL = 5;

    /** 6 members per cube (see Necronomata javadoc). */
    private static final int MEMBERS_PER_CUBE = 6;

    private final Map<Integer, Long> perDirection;
    private final Map<Integer, Long> effectivePerDirection;
    private final Map<Integer, Long> perTick = new TreeMap<>();
    private final Map<Long, Long> magnitudeHistogram = new TreeMap<>();
    private final Map<Long, Long> perMemberPair = new TreeMap<>();

    private long totalCollisions = 0L;
    private long effectiveCollisions = 0L;
    private int minTick = Integer.MAX_VALUE;
    private int maxTick = Integer.MIN_VALUE;

    public CollisionStatistics() {
        Map<Integer, Long> directions = new LinkedHashMap<>();
        Map<Integer, Long> effectiveDirections = new LinkedHashMap<>();
        for (int direction : FccNeighborhood.DIRECTIONS) {
            directions.put(direction, 0L);
            effectiveDirections.put(direction, 0L);
        }
        this.perDirection = directions;
        this.effectivePerDirection = effectiveDirections;
    }

    /**
     * @return the total number of RESOLVED contacts recorded across every
     *         direction -- including zero-transfer no-ops; must equal the
     *         sum of {@link #collisionsPerDirection()}'s values by
     *         construction (both are updated together). See class
     *         Javadoc, "What counts as a collision"; use {@link
     *         #effectiveCollisions()} for the transfer-only count.
     */
    public long totalCollisions() {
        return totalCollisions;
    }

    /**
     * @return the number of recorded collisions with a nonzero {@code
     *         transferMagnitude} -- the narrower, transfer-only count
     *         complementing {@link #totalCollisions()} (see class
     *         Javadoc, "What counts as a collision"). Always {@code <=
     *         totalCollisions()}.
     */
    public long effectiveCollisions() {
        return effectiveCollisions;
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

    /**
     * Fix-round item 2c (bead inviscid-0nx.23, round 3): the narrower,
     * EFFECTIVE-transfer-only counterpart to {@link
     * #collisionsPerDirection()} -- one entry for each of the 12 {@link
     * FccNeighborhood#DIRECTIONS} (present with a zero count even if never
     * observed), counting only recorded collisions with {@code
     * transferMagnitude > 0} (same narrowing {@link #effectiveCollisions()}
     * applies to {@link #totalCollisions()}, see class Javadoc "What counts
     * as a collision"). Additive instrumentation only: accumulated
     * alongside {@link #perDirection} in {@link #recordCollision}, changes
     * nothing about what is recorded or how.
     *
     * @return the per-direction EFFECTIVE collision count, in {@link
     *         FccNeighborhood#DIRECTIONS} order
     */
    public Map<Integer, Long> effectiveCollisionsPerDirection() {
        return Collections.unmodifiableMap(effectivePerDirection);
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
     * @return the collision count for the unordered pair of
     *         {@code (cube, member)} slots {@code (cubeA, memberA)} and
     *         {@code (cubeB, memberB)}; canonicalized so
     *         {@code collisionsForMemberPair(cubeA, memberA, cubeB, memberB)
     *         == collisionsForMemberPair(cubeB, memberB, cubeA, memberA)}.
     *         Two slots sharing a within-cube member index but belonging
     *         to different cubes are distinct pairs -- {@code cube} is
     *         part of the key, never discarded.
     */
    public long collisionsForMemberPair(int cubeA, int memberA, int cubeB,
                                         int memberB) {
        Long count = perMemberPair.get(memberPairKey(slotIndex(cubeA,
                                                                memberA),
                                                       slotIndex(cubeB,
                                                                memberB)));
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
     * often does a collision happen" signal only. Uses {@link
     * #totalCollisions()} (includes no-ops -- see class Javadoc, "What
     * counts as a collision"), so this is a "how often is a member near
     * another" proxy, coarser than a transfer-only mean free path would
     * be; divide the span by {@link #effectiveCollisions()} instead if a
     * transfer-only proxy is wanted.
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
     * quanta transfer, and call this once per resolved contact -- the
     * parameter order mirrors
     * {@code com.chiralbehaviors.inviscid.lga.Contact}'s fields
     * ({@code cellA}, {@code cubeA}, {@code memberA}, {@code cellB},
     * {@code cubeB}, {@code memberB}, {@code direction}) so those seven
     * {@code Contact} fields feed the first seven parameters positionally;
     * {@code Contact.minDistance()} has no parameter here, and
     * {@code transferMagnitude}/{@code tick} are the rule's own outputs.
     *
     * @param cellA             the first cell involved
     * @param cubeA             the cube index (within {@code cellA}),
     *                          must be in {@code [0, CUBES_PER_CELL)}
     *                          i.e. {@code [0, 4]}
     * @param memberA           the member index (within {@code cubeA} of
     *                          {@code cellA}), must be in
     *                          {@code [0, MEMBERS_PER_CUBE)} i.e.
     *                          {@code [0, 5]}
     * @param cellB             the second cell involved
     * @param cubeB             the cube index (within {@code cellB}),
     *                          must be in {@code [0, CUBES_PER_CELL)}
     *                          i.e. {@code [0, 4]}
     * @param memberB           the member index (within {@code cubeB} of
     *                          {@code cellB}), must be in
     *                          {@code [0, MEMBERS_PER_CUBE)} i.e.
     *                          {@code [0, 5]}
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
     * @throws IllegalArgumentException if {@code cubeA}/{@code cubeB} is
     *                                   outside {@code [0, CUBES_PER_CELL)},
     *                                   {@code memberA}/{@code memberB} is
     *                                   outside
     *                                   {@code [0, MEMBERS_PER_CUBE)},
     *                                   {@code direction} is not one of
     *                                   {@link FccNeighborhood#DIRECTIONS},
     *                                   or {@code transferMagnitude} is
     *                                   negative
     */
    public void recordCollision(Point3i cellA, int cubeA, int memberA,
                                 Point3i cellB, int cubeB, int memberB,
                                 int direction, long transferMagnitude,
                                 int tick) {
        validateCube(cubeA, "cubeA");
        validateMember(memberA, "memberA");
        validateCube(cubeB, "cubeB");
        validateMember(memberB, "memberB");
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
        if (transferMagnitude > 0) {
            effectiveCollisions++;
            effectivePerDirection.merge(direction, 1L, Long::sum);
        }
        perDirection.merge(direction, 1L, Long::sum);
        perTick.merge(tick, 1L, Long::sum);
        perMemberPair.merge(memberPairKey(slotIndex(cubeA, memberA),
                                           slotIndex(cubeB, memberB)), 1L,
                             Long::sum);
        magnitudeHistogram.merge(transferMagnitude, 1L, Long::sum);
        if (tick < minTick) {
            minTick = tick;
        }
        if (tick > maxTick) {
            maxTick = tick;
        }
    }

    private static void validateCube(int cube, String argName) {
        if (cube < 0 || cube >= CUBES_PER_CELL) {
            throw new IllegalArgumentException(argName + " must be in [0, "
                                                + (CUBES_PER_CELL - 1)
                                                + "], was: " + cube);
        }
    }

    private static void validateMember(int member, String argName) {
        if (member < 0 || member >= MEMBERS_PER_CUBE) {
            throw new IllegalArgumentException(argName + " must be in [0, "
                                                + (MEMBERS_PER_CUBE - 1)
                                                + "], was: " + member);
        }
    }

    /**
     * @return the flat {@code [0, CUBES_PER_CELL * MEMBERS_PER_CUBE)} slot
     *         index for {@code (cube, member)} -- the same
     *         {@code cube * MEMBERS_PER_CUBE + member} addressing
     *         {@code ConservationAudit} uses internally, kept here only as
     *         a member-pair key component (never exposed).
     */
    private static int slotIndex(int cube, int member) {
        return cube * MEMBERS_PER_CUBE + member;
    }

    private static long memberPairKey(int slotA, int slotB) {
        int lo = Math.min(slotA, slotB);
        int hi = Math.max(slotA, slotB);
        return ((long) lo << 32) | (hi & 0xFFFFFFFFL);
    }
}
