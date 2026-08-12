/**
 * Copyright (c) 2026 Chiral Behaviors, LLC, all rights reserved.
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

package com.chiralbehaviors.inviscid.jitterbug;

import java.util.ArrayList;
import java.util.List;

/**
 * Self-intersection of the eight triangles, counted as edge-through-face events.
 *
 * <p>
 * Moller-Trumbore ray/triangle intersection over every ordered pair of distinct
 * triangles, with a strict interior test: an edge must pierce the far triangle
 * strictly inside its boundary and strictly between its own endpoints.
 *
 * <p>
 * <b>Hinge contacts are excluded.</b> Two triangles that share a hinge vertex
 * touch there by construction; that is legal, so any pierce point that lands on
 * a vertex the two triangles share is discarded.
 *
 * <p>
 * On the symmetric family the exclusion turns out to be <b>inert</b>: switching
 * it off changes the count at no angle over the whole circle, measured in
 * {@code InterpenetrationTest.theHingeExclusionChangesNothingOnTheSymmetricPath}.
 * That strengthens the finding below rather than weakening it — the contact
 * boundary is not an artefact of the contact rule, since it survives with the
 * rule removed. The rule is kept because it is correct for off-path
 * configurations, where a shared vertex can fall strictly inside the far
 * triangle.
 *
 * <p>
 * <b>What this found.</b> Nobody was looking for it: over the whole circle the
 * count takes only the values 0 and 24, and the collision set is two intervals.
 * In the {@code eps -> 0} limit they are {@code (60, 120)} and
 * {@code (240, 300)}. At the {@link #DEFAULT_EPS} margin actually used the
 * measured transition sits slightly inside — 60.001109528 at the lower end, and
 * symmetrically inside 120 at the upper — because <b>onset is tangential</b>, so
 * the observed onset angle is margin-limited and scales as {@code sqrt(eps)};
 * both the onset angle and the scaling are measured in
 * {@code InterpenetrationTest.onsetIsAtTheOctahedronAndTheResidualGapIsToleranceLimited}.
 * "Onset precisely at the octahedron" is therefore the limit statement, not the
 * measurement, and the two must not be quoted interchangeably.
 *
 * <p>
 * That still makes {@code a = 60} the <b>contact boundary</b> — a real
 * distinction, after its previously claimed distinction as a rank singularity
 * was refuted and the phenomenon was nearly discarded along with the mechanism.
 * The admissible symmetric path is therefore the interval [-60, 60], not a
 * circle, and there is no 2*pi circuit.
 *
 * <p>
 * No new dependencies: plain arithmetic only.
 *
 * @author halhildebrand
 */
public final class Interpenetration {

    /**
     * Default barycentric / parametric margin. Matches the Python reference.
     * Because onset is tangential, the observed onset angle scales as
     * {@code sqrt(eps)} — see the test that measures the scaling.
     */
    public static final double DEFAULT_EPS   = 1e-9;

    /** How close a pierce point must be to a shared vertex to count as contact. */
    public static final double CONTACT_TOL   = 1e-7;

    /** Below this, a triangle or a ray is treated as degenerate and skipped. */
    private static final double DEGENERATE   = 1e-12;

    private Interpenetration() {
    }

    /** Edge-through-face count at the default margin. */
    public static int count(double[][][] x) {
        return count(x, DEFAULT_EPS);
    }

    /**
     * Edge-through-face count over all 56 ordered pairs of distinct triangles.
     *
     * @param eps strict-interior margin in barycentric and ray parameters
     */
    public static int count(double[][][] x, double eps) {
        int total = 0;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                if (i == j) {
                    continue;
                }
                total += edgesThrough(x[i], x[j], sharedCorners(x, i, j), eps);
            }
        }
        return total;
    }

    /**
     * Count edges of triangle {@code p} that pierce the interior of triangle
     * {@code q}, ignoring pierce points at any of the {@code shared} vertices.
     */
    static int edgesThrough(double[][] p, double[][] q, List<double[]> shared,
                            double eps) {
        double[] e1 = sub(q[1], q[0]);
        double[] e2 = sub(q[2], q[0]);
        if (norm(cross(e1, e2)) < DEGENERATE) {
            return 0;
        }
        int hits = 0;
        for (int k = 0; k < 3; k++) {
            double[] o = p[k];
            double[] d = sub(p[(k + 1) % 3], p[k]);

            double[] pv = cross(d, e2);
            double det = dot(e1, pv);
            if (Math.abs(det) < DEGENERATE) {
                continue; // parallel; a coplanar overlap is not an edge pierce
            }
            double inv = 1.0 / det;
            double[] tv = sub(o, q[0]);
            double u = dot(tv, pv) * inv;
            if (u < eps || u > 1 - eps) {
                continue;
            }
            double[] qv = cross(tv, e1);
            double v = dot(d, qv) * inv;
            if (v < eps || u + v > 1 - eps) {
                continue;
            }
            double t = dot(e2, qv) * inv;
            if (t < eps || t > 1 - eps) {
                continue;
            }
            double[] hit = { o[0] + t * d[0], o[1] + t * d[1],
                             o[2] + t * d[2] };
            boolean atHinge = false;
            for (double[] s : shared) {
                if (JitterbugGeometry.distance(hit, s) < CONTACT_TOL) {
                    atHinge = true;
                    break;
                }
            }
            if (!atHinge) {
                hits++;
            }
        }
        return hits;
    }

    /**
     * The corners of face {@code i} that are hinged to face {@code j}, as
     * positions.
     */
    static List<double[]> sharedCorners(double[][][] x, int i, int j) {
        List<double[]> out = new ArrayList<>(3);
        for (int[] pair : JitterbugLinkage.hingePairs()) {
            if (pair[0] == i && pair[2] == j) {
                out.add(x[i][pair[1]]);
            } else if (pair[2] == i && pair[0] == j) {
                out.add(x[i][pair[3]]);
            }
        }
        return out;
    }

    private static double[] cross(double[] a, double[] b) {
        return new double[] { a[1] * b[2] - a[2] * b[1],
                              a[2] * b[0] - a[0] * b[2],
                              a[0] * b[1] - a[1] * b[0] };
    }

    private static double dot(double[] a, double[] b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    private static double norm(double[] a) {
        return Math.sqrt(dot(a, a));
    }

    private static double[] sub(double[] a, double[] b) {
        return new double[] { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
    }
}
