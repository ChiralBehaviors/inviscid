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

import javax.vecmath.Vector3d;

/**
 * The minimum Euclidean distance between two line segments in 3-space -
 * the standard closest-point-between-segments routine, with the
 * parallel / collinear / degenerate cases handled EXPLICITLY rather than
 * falling through a single fixed-parameter formula.
 * <p>
 * <b>Why parallel needs its own branch.</b> The textbook closest-point
 * algorithm (Ericson, <i>Real-Time Collision Detection</i> &sect;5.1.9)
 * has a well-known trap: when the two segments are parallel the
 * {@code a*e - b*b} denominator underlying its general-case formula is
 * zero, and the naive fallback fixes one segment's parameter at an
 * arbitrary endpoint (typically {@code s=0}) before solving for the
 * other. That produces a distance that is a valid upper bound but not
 * always the true minimum - see {@code
 * SegmentDistanceTest.parallelSegmentsUseEndpointProjection} for a
 * concrete case where it is wrong. This class instead detects the
 * parallel condition explicitly (via the squared cross-product magnitude
 * of the two direction vectors, which is exactly the Ericson denominator
 * by the Lagrange identity) and falls back to full endpoint projection:
 * project each of the four endpoints onto the opposite segment (clamped
 * to {@code [0,1]}) and take the minimum of the four candidate
 * distances. That is provably correct for parallel segments - overlapping
 * or not, collinear or merely coplanar-parallel - because the closest
 * pair between two parallel segments is always achieved either within a
 * mutual overlap (distance 0, found because the overlapping endpoint's
 * clamped projection lands inside the other segment) or at one
 * segment's nearest endpoint projected onto the other (found directly by
 * the same clamped projection).
 * <p>
 * <b>Exact symmetry.</b> {@link #distance(Vector3d, Vector3d, Vector3d,
 * Vector3d)} canonicalizes its four input points - both which segment is
 * "first" and each segment's own endpoint order - into a single fixed
 * argument order before doing any arithmetic. A segment is a set of two
 * points, not a direction, so every one of the 8 argument permutations
 * describing the same unordered pair of segments resolves to the
 * identical internal call and therefore the identical floating-point
 * result, bit for bit - not merely equal within tolerance. This is what
 * makes {@code SegmentDistanceTest.distanceIsSymmetric} an exact-equality
 * assertion rather than a delta comparison, and is load-bearing for
 * {@code ContactPredicate}'s direction-reversal symmetry
 * (inviscid-0nx.12 test 7): that predicate's boolean outcome must not
 * depend on which cell/member is labeled "A" and which is "B".
 *
 * @author halhildebrand
 */
public final class SegmentDistance {

    /**
     * Below this squared length, a segment is treated as degenerate (a
     * point) rather than risking division by a near-zero direction
     * length. {@code 1e-12} corresponds to a segment shorter than 1e-6 in
     * length - far below any physically meaningful member/strut scale in
     * this codebase (cell spacing and member half-lengths are O(1..10)).
     */
    private static final double DEGENERATE_LENGTH_SQUARED = 1e-12;

    /**
     * Two direction vectors are treated as parallel when {@code |d1 x
     * d2|^2 <= PARALLEL_EPSILON * |d1|^2 * |d2|^2}. By the Lagrange
     * identity {@code |d1 x d2|^2 == |d1|^2|d2|^2 - (d1.d2)^2}, dividing
     * both sides by {@code |d1|^2|d2|^2} shows this threshold is exactly
     * {@code sin^2(angle between d1,d2) <= PARALLEL_EPSILON} - a
     * scale-independent angular tolerance, not an absolute one.
     */
    private static final double PARALLEL_EPSILON = 1e-9;

    /**
     * The minimum Euclidean distance between segment {@code [p1,q1]} and
     * segment {@code [p2,q2]}.
     */
    public static double distance(Vector3d p1, Vector3d q1, Vector3d p2,
                                  Vector3d q2) {
        Vector3d a1 = p1;
        Vector3d b1 = q1;
        if (compareVec(a1, b1) > 0) {
            Vector3d t = a1;
            a1 = b1;
            b1 = t;
        }
        Vector3d a2 = p2;
        Vector3d b2 = q2;
        if (compareVec(a2, b2) > 0) {
            Vector3d t = a2;
            a2 = b2;
            b2 = t;
        }
        int cmp = compareVec(a1, a2);
        if (cmp == 0) {
            cmp = compareVec(b1, b2);
        }
        if (cmp <= 0) {
            return distanceOrdered(a1, b1, a2, b2);
        }
        return distanceOrdered(a2, b2, a1, b1);
    }

    /**
     * Convenience overload for {@link Segment} endpoints.
     */
    public static double distance(Segment a, Segment b) {
        return distance(a.getA(), a.getB(), b.getA(), b.getB());
    }

    private static Vector3d addScaled(Vector3d base, Vector3d dir, double t) {
        Vector3d r = new Vector3d(dir);
        r.scale(t);
        r.add(base);
        return r;
    }

    private static double clamp01(double v) {
        if (v < 0) {
            return 0;
        }
        if (v > 1) {
            return 1;
        }
        return v;
    }

    private static int compareVec(Vector3d a, Vector3d b) {
        int c = Double.compare(a.x, b.x);
        if (c != 0) {
            return c;
        }
        c = Double.compare(a.y, b.y);
        if (c != 0) {
            return c;
        }
        return Double.compare(a.z, b.z);
    }

    private static Vector3d cross(Vector3d a, Vector3d b) {
        Vector3d r = new Vector3d();
        r.cross(a, b);
        return r;
    }

    private static double dist(Vector3d a, Vector3d b) {
        Vector3d d = new Vector3d(a);
        d.sub(b);
        return d.length();
    }

    /**
     * @return the distance from {@code point} to the closest point on
     *         segment {@code [segStart,segEnd]}, clamped to the segment.
     */
    private static double distancePointToSegment(Vector3d point,
                                                  Vector3d segStart,
                                                  Vector3d segEnd) {
        Vector3d d = new Vector3d(segEnd);
        d.sub(segStart);
        double lenSq = d.dot(d);
        double t;
        if (lenSq <= DEGENERATE_LENGTH_SQUARED) {
            t = 0;
        } else {
            Vector3d diff = new Vector3d(point);
            diff.sub(segStart);
            t = clamp01(diff.dot(d) / lenSq);
        }
        return dist(addScaled(segStart, d, t), point);
    }

    /**
     * The core, order-fixed computation. Callers must pass a canonical
     * argument order (see {@link #distance(Vector3d, Vector3d, Vector3d,
     * Vector3d)}) for the exact-symmetry guarantee to hold.
     */
    private static double distanceOrdered(Vector3d p1, Vector3d q1,
                                          Vector3d p2, Vector3d q2) {
        Vector3d d1 = new Vector3d(q1);
        d1.sub(p1);
        Vector3d d2 = new Vector3d(q2);
        d2.sub(p2);
        double a = d1.dot(d1);
        double e = d2.dot(d2);

        if (a <= DEGENERATE_LENGTH_SQUARED
            && e <= DEGENERATE_LENGTH_SQUARED) {
            return dist(p1, p2);
        }
        if (a <= DEGENERATE_LENGTH_SQUARED) {
            return distancePointToSegment(p1, p2, q2);
        }
        if (e <= DEGENERATE_LENGTH_SQUARED) {
            return distancePointToSegment(p2, p1, q1);
        }

        Vector3d crossD = cross(d1, d2);
        double crossLenSq = crossD.dot(crossD);
        if (crossLenSq <= PARALLEL_EPSILON * a * e) {
            return distanceViaEndpointProjection(p1, q1, p2, q2);
        }
        return distanceViaClosestPoints(p1, d1, a, p2, d2, e, crossLenSq);
    }

    /**
     * The general (non-degenerate, non-parallel) closest-point-between-
     * segments solve: minimize {@code |(p1+s*d1) - (p2+t*d2)|} over
     * {@code s,t in [0,1]}, clamping onto whichever segment boundary the
     * unconstrained solution falls outside of.
     */
    private static double distanceViaClosestPoints(Vector3d p1, Vector3d d1,
                                                    double a, Vector3d p2,
                                                    Vector3d d2, double e,
                                                    double denom) {
        Vector3d r = new Vector3d(p1);
        r.sub(p2);
        double b = d1.dot(d2);
        double c = d1.dot(r);
        double f = d2.dot(r);

        double s = clamp01((b * f - c * e) / denom);
        double t = (b * s + f) / e;
        if (t < 0) {
            t = 0;
            s = clamp01(-c / a);
        } else if (t > 1) {
            t = 1;
            s = clamp01((b - c) / a);
        }
        return dist(addScaled(p1, d1, s), addScaled(p2, d2, t));
    }

    /**
     * Parallel-segment fallback: the minimum of the four endpoint-to-
     * opposite-segment projected distances. See the class Javadoc for why
     * this is correct where the naive fixed-parameter formula is not.
     */
    private static double distanceViaEndpointProjection(Vector3d p1,
                                                         Vector3d q1,
                                                         Vector3d p2,
                                                         Vector3d q2) {
        double d = distancePointToSegment(p1, p2, q2);
        d = Math.min(d, distancePointToSegment(q1, p2, q2));
        d = Math.min(d, distancePointToSegment(p2, p1, q1));
        d = Math.min(d, distancePointToSegment(q2, p1, q1));
        return d;
    }

    private SegmentDistance() {
    }
}
