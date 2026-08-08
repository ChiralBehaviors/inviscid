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

import java.util.Random;

import javax.vecmath.Vector3d;

import org.junit.Test;

/**
 * Hand-computed geometry tests for {@link SegmentDistance}, the standard
 * closest-point-between-segments routine with parallel / collinear /
 * degenerate cases handled explicitly.
 *
 * @author halhildebrand
 */
public class SegmentDistanceTest {

    private static final double DELTA = 1e-9;

    /**
     * Two segments crossing in the same plane (z = 0), intersecting exactly
     * at the origin: {@code (-1,-1,0)-(1,1,0)} and {@code (-1,1,0)-(1,-1,0)}.
     * The generic (non-parallel) closest-point solver must land its
     * unclamped {@code (s,t)} exactly at the intersection, giving distance
     * 0.
     */
    @Test
    public void crossingSegmentsHaveZeroDistance() {
        Vector3d p1 = new Vector3d(-1, -1, 0);
        Vector3d q1 = new Vector3d(1, 1, 0);
        Vector3d p2 = new Vector3d(-1, 1, 0);
        Vector3d q2 = new Vector3d(1, -1, 0);

        assertEquals(0.0, SegmentDistance.distance(p1, q1, p2, q2), DELTA);
    }

    /**
     * The classic naive-implementation failure: two parallel, non-
     * overlapping segments where the true closest points are each an
     * <em>endpoint</em> of one segment projected (and clamped) onto the
     * other - not the endpoints the naive Ericson-style fixed-parameter
     * fallback picks.
     * <p>
     * segment 1: {@code (0,0,0)-(1,0,0)} (along x, y=z=0). segment 2:
     * {@code (5,1,0)-(10,1,0)} (along x, offset by y=1, x in [5,10]). The
     * true closest pair is segment 1's {@code (1,0,0)} endpoint against
     * segment 2's {@code (5,1,0)} endpoint: distance
     * {@code sqrt((5-1)^2 + (1-0)^2) == sqrt(17)}. A naive implementation
     * that fixes {@code s=0} for the parallel branch (i.e. always uses
     * segment 1's {@code p1} rather than checking both endpoints) computes
     * distance from {@code (0,0,0)} to the clamped projection of itself
     * onto segment 2, which lands on segment 2's near endpoint
     * {@code (5,1,0)} too, but the DISTANCE is computed from the wrong
     * point on segment 1 ({@code p1} instead of {@code q1}), giving
     * {@code sqrt(5^2+1^2) == sqrt(26)} - larger than the true minimum.
     */
    @Test
    public void parallelSegmentsUseEndpointProjection() {
        Vector3d p1 = new Vector3d(0, 0, 0);
        Vector3d q1 = new Vector3d(1, 0, 0);
        Vector3d p2 = new Vector3d(5, 1, 0);
        Vector3d q2 = new Vector3d(10, 1, 0);

        double expected = Math.sqrt(17.0);
        assertEquals(expected, SegmentDistance.distance(p1, q1, p2, q2), DELTA);
    }

    /**
     * Collinear segments whose parameter ranges overlap:
     * {@code (0,0,0)-(5,0,0)} and {@code (3,0,0)-(10,0,0)} overlap on
     * {@code x in [3,5]} - distance must be exactly 0.
     */
    @Test
    public void collinearOverlappingSegments() {
        Vector3d p1 = new Vector3d(0, 0, 0);
        Vector3d q1 = new Vector3d(5, 0, 0);
        Vector3d p2 = new Vector3d(3, 0, 0);
        Vector3d q2 = new Vector3d(10, 0, 0);

        assertEquals(0.0, SegmentDistance.distance(p1, q1, p2, q2), DELTA);
    }

    /**
     * Genuinely skew segments (not parallel, not intersecting) whose
     * closest points fall strictly inside both parameter ranges - so no
     * clamping is exercised, and the segment distance must match the
     * closed-form skew-line distance {@code |r . (d1 x d2)| / |d1 x d2|}.
     * <p>
     * segment 1: {@code (-1,0,0)-(1,0,0)} (x in [-1,1], y=z=0). segment 2:
     * {@code (0,-1,1)-(0,1,1)} (y in [-1,1], x=0, z=1). The infinite-line
     * closest points are {@code (0,0,0)} (s=0.5, interior) and
     * {@code (0,0,1)} (t=0.5, interior); distance 1. Closed form:
     * {@code d1=(2,0,0)}, {@code d2=(0,2,0)}, {@code d1 x d2 = (0,0,4)},
     * {@code r=p1-p2=(-1,1,-1)}, {@code r . (d1 x d2) = -4},
     * {@code |-4| / |(0,0,4)| == 1}.
     */
    @Test
    public void skewSegmentsMatchClosedForm() {
        Vector3d p1 = new Vector3d(-1, 0, 0);
        Vector3d q1 = new Vector3d(1, 0, 0);
        Vector3d p2 = new Vector3d(0, -1, 1);
        Vector3d q2 = new Vector3d(0, 1, 1);

        assertEquals(1.0, SegmentDistance.distance(p1, q1, p2, q2), DELTA);
    }

    /**
     * {@code distance(a,b) == distance(b,a)} EXACTLY (bit-for-bit, not just
     * within tolerance) for seeded random segment pairs. Guards against
     * implementation asymmetries (e.g. an algorithm that clamps one
     * segment's parameter before the other) that would otherwise silently
     * make {@code ContactPredicate}'s direction-reversal symmetry
     * (inviscid-0nx.12 test 7) depend on argument order.
     */
    @Test
    public void distanceIsSymmetric() {
        Random random = new Random(42L);
        for (int i = 0; i < 500; i++) {
            Vector3d p1 = randomVector(random);
            Vector3d q1 = randomVector(random);
            Vector3d p2 = randomVector(random);
            Vector3d q2 = randomVector(random);

            double forward = SegmentDistance.distance(p1, q1, p2, q2);
            double backward = SegmentDistance.distance(p2, q2, p1, q1);
            assertEquals("iteration " + i, forward, backward, 0.0);

            // Swapping each segment's own endpoints must also be exact -
            // the segment is undirected.
            double forwardFlippedFirst = SegmentDistance.distance(q1, p1, p2,
                                                                   q2);
            assertEquals("iteration " + i + " (flip first)", forward,
                        forwardFlippedFirst, 0.0);
        }
    }

    private static Vector3d randomVector(Random random) {
        return new Vector3d(random.nextDouble() * 20 - 10,
                            random.nextDouble() * 20 - 10,
                            random.nextDouble() * 20 - 10);
    }
}
