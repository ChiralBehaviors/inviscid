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

import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.A_ICOSAHEDRON_DEG;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.A_OCTAHEDRON_DEG;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.L_EDGE;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.R_CIRC;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.Z;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.cluster;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.corners;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.strutLengths;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

import org.junit.Test;

/**
 * The symmetric jitterbug family, pinned.
 *
 * <p>
 * Every number asserted here was measured independently by the Python reference
 * harness in {@code analysis/jitterbug-variety/} and is reproduced here by this
 * Java implementation. Where the two disagree one of them is wrong; that is the
 * point of keeping both.
 *
 * <p>
 * <b>A check that cannot fail is not a check.</b> The reference
 * {@code jb_a_family.py} verified the icosahedron's "exactly 30 equal short
 * edges" by slicing the sorted distance list to its 30 smallest entries and then
 * counting how many of <em>those</em> matched the minimum — an assertion with no
 * failing branch. {@link #icosahedronHasExactlyThirtyEqualShortEdges()} counts
 * over all 66 vertex pairs instead, so a 29 or a 31 would fail it.
 *
 * @author halhildebrand
 */
public class JitterbugGeometryTest {

    private static final double MERGE_TOL = 1e-7;

    /**
     * Item 1. At a=60 the parameterisation returns the octahedron vertices
     * exactly: six points, one per signed coordinate axis, at circumradius R.
     */
    @Test
    public void octahedronAtSixtyDegrees() {
        JitterbugGeometry.Clustering c = cluster(corners(A_OCTAHEDRON_DEG),
                                                 MERGE_TOL);
        assertEquals("a=60 must collapse to the 6 octahedron vertices", 6,
                     c.representatives().length);

        double[][] expected = { { R_CIRC, 0, 0 }, { -R_CIRC, 0, 0 },
                                { 0, R_CIRC, 0 }, { 0, -R_CIRC, 0 },
                                { 0, 0, R_CIRC }, { 0, 0, -R_CIRC } };
        for (double[] want : expected) {
            double best = Double.MAX_VALUE;
            for (double[] got : c.representatives()) {
                best = Math.min(best, distance(want, got));
            }
            assertEquals("octahedron vertex " + Arrays.toString(want)
                         + " not produced at a=60", 0.0, best, 1e-12);
        }
    }

    /**
     * Item 2, part 1. Twelve distinct corners of multiplicity 2 at a generic a.
     */
    @Test
    public void genericAngleGivesTwelveDoubledCorners() {
        for (double a : new double[] { 0, 15, 22.24, 30, 45, 59, 61, 90, 137,
                                       180, 200, 271, 333 }) {
            JitterbugGeometry.Clustering c = cluster(corners(a), MERGE_TOL);
            assertEquals("distinct corners at a=" + a, 12,
                         c.representatives().length);
            for (int m : c.multiplicities()) {
                assertEquals("multiplicity at a=" + a, 2, m);
            }
        }
    }

    /**
     * Item 2, part 2. The 12 to 6 merge happens at a = 60, 120, 240 and 300 —
     * <em>not</em> at 60 alone, which is what the superseded memo recorded. A
     * 0.5 degree sweep of the whole circle finds those four angles and no
     * others.
     */
    @Test
    public void sixfoldMergeHappensAtFourAnglesNotOne() {
        List<Double> merged = new ArrayList<>();
        for (int k = 0; k < 720; k++) {
            double a = k * 0.5;
            JitterbugGeometry.Clustering c = cluster(corners(a), MERGE_TOL);
            if (c.representatives().length != 12) {
                assertEquals("merge at a=" + a + " must be 6x4", 6,
                             c.representatives().length);
                for (int m : c.multiplicities()) {
                    assertEquals("multiplicity at a=" + a, 4, m);
                }
                merged.add(a);
            }
        }
        assertEquals("merge angles on a 0.5 degree sweep of [0,360)",
                     Arrays.asList(60.0, 120.0, 240.0, 300.0), merged);
    }

    /**
     * Item 2, part 3. Chirality is forced. With sigma = +1 on every face the
     * sharing breaks: 24 distinct corners of multiplicity 1 at a generic a. The
     * correct sigma = sx*sy*sz is not a convention, it is the only choice that
     * makes the linkage a linkage.
     */
    @Test
    public void uniformSigmaBreaksVertexSharing() {
        for (double a : new double[] { 5, 15, 30, 45, 59, 61, 75 }) {
            JitterbugGeometry.Clustering c = cluster(corners(a, true),
                                                     MERGE_TOL);
            assertEquals("uniform-sigma distinct corners at a=" + a, 24,
                         c.representatives().length);
            for (int m : c.multiplicities()) {
                assertEquals(1, m);
            }
        }
        // The angles where the two sigma choices happen to agree, and why
        // sampling only these would have hidden the difference entirely.
        // Measured, not assumed: a=90 was originally written here as a
        // sharing-broken angle and the implementation disagreed; the reference
        // harness confirmed the implementation.
        assertEquals("a=0: no rotation applied yet, both conventions coincide",
                     12, cluster(corners(0, true), MERGE_TOL)
                                                             .representatives().length);
        assertEquals("a=60: zero rotation angle, both give the octahedron", 6,
                     cluster(corners(60, true), MERGE_TOL)
                                                          .representatives().length);
        assertEquals("a=120: 6x4 for uniform sigma too", 6,
                     cluster(corners(120, true), MERGE_TOL)
                                                           .representatives().length);
        for (double a : new double[] { 90, 180, 270 }) {
            assertEquals("uniform sigma recovers 12x2 at a=" + a, 12,
                         cluster(corners(a, true), MERGE_TOL)
                                                             .representatives().length);
        }
    }

    /**
     * Item 3. The icosahedron sits at a = 22.238756093 degrees, where the open
     * square diagonals reach strut length and there are exactly 30 equal short
     * edges.
     *
     * <p>
     * The angle is located here by bisection on
     * {@link JitterbugGeometry#shortestNonStrutDistance(double)} rather than
     * read from the constant, so {@code A_ICOSAHEDRON_DEG} is itself under
     * test.
     */
    @Test
    public void icosahedronAngleIsLocatedByBisection() {
        double lo = 1.0;
        double hi = 40.0;
        double flo = JitterbugGeometry.shortestNonStrutDistance(lo) - L_EDGE;
        double fhi = JitterbugGeometry.shortestNonStrutDistance(hi) - L_EDGE;
        assertTrue("root must be bracketed", flo * fhi < 0);
        for (int i = 0; i < 200; i++) {
            double mid = 0.5 * (lo + hi);
            double fm = JitterbugGeometry.shortestNonStrutDistance(mid) - L_EDGE;
            if (flo * fm <= 0) {
                hi = mid;
                fhi = fm;
            } else {
                lo = mid;
                flo = fm;
            }
        }
        double measured = 0.5 * (lo + hi);
        assertEquals("icosahedron angle", A_ICOSAHEDRON_DEG, measured, 1e-9);
        assertEquals("icosahedron angle to the recorded 9 decimals",
                     22.238756093, measured, 5e-10);
    }

    /**
     * Item 3, the falsifiable half. Count <em>all</em> 66 vertex pairs within
     * tolerance of the minimum distance and require exactly 30. The reference
     * implementation could not fail this; this one can.
     */
    @Test
    public void icosahedronHasExactlyThirtyEqualShortEdges() {
        JitterbugGeometry.Clustering c = cluster(corners(A_ICOSAHEDRON_DEG),
                                                 MERGE_TOL);
        double[][] v = c.representatives();
        assertEquals(12, v.length);

        List<Double> all = new ArrayList<>();
        for (int i = 0; i < v.length; i++) {
            for (int j = i + 1; j < v.length; j++) {
                all.add(distance(v[i], v[j]));
            }
        }
        assertEquals("an icosahedron's 12 vertices have 66 pairs", 66,
                     all.size());
        double min = all.stream().mapToDouble(Double::doubleValue).min()
                        .getAsDouble();
        assertEquals("the short edge is the octahedron strut length", L_EDGE,
                     min, 1e-9);

        // The count must be stable across four decades of tolerance: if it were
        // an artefact of a lucky epsilon, widening would sweep in the next
        // distance shell.
        for (double tol : new double[] { 1e-9, 1e-7, 1e-6, 1e-4 }) {
            long n = all.stream().filter(d -> d < min + tol).count();
            assertEquals("short edges within " + tol + " of the minimum", 30, n);
        }

        // ...and there really is a next shell, so 30 is a boundary, not a cap.
        double nextShell = all.stream().mapToDouble(Double::doubleValue)
                              .filter(d -> d > min + 1e-4).min().getAsDouble();
        assertTrue("the 31st distance must be well separated (measured 2.288...)",
                   nextShell > min + 0.5);
        assertEquals(2.2882456, nextShell, 1e-6);
    }

    /**
     * Item 2 support. The vector equilibrium at a=0 is the cuboctahedron:
     * circumradius equals edge length.
     */
    @Test
    public void vectorEquilibriumHasRadiusEqualToEdge() {
        JitterbugGeometry.Clustering c = cluster(corners(0.0), MERGE_TOL);
        assertEquals(12, c.representatives().length);
        for (double[] p : c.representatives()) {
            assertEquals("VE circumradius == edge", L_EDGE,
                         Math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]),
                         1e-12);
        }
    }

    /**
     * NOT A VERIFICATION — A GUARD. Equal struts at every a is a tautology of a
     * rigid-body parameterisation: the corners of each triangle are moved by a
     * single rotation and translation, so their mutual distances cannot change.
     * This test can only fail if the implementation stops being a rigid-body
     * parameterisation (a per-corner scale, a wrong pivot, a lost translation).
     * It is kept for that, and must not be cited as evidence that the geometry
     * is right.
     */
    @Test
    public void strutLengthsAreRigid_guardNotVerification() {
        for (double a : new double[] { 0, 22.24, 30, 60, 90, 180, 271 }) {
            double[] e = strutLengths(a);
            assertEquals(24, e.length);
            for (double len : e) {
                assertEquals("strut length at a=" + a, L_EDGE, len, 1e-12);
            }
        }
    }

    /**
     * The constants the whole harness rests on.
     */
    @Test
    public void constants() {
        assertEquals(1.0, R_CIRC, 0.0);
        assertEquals(Math.sqrt(2.0), L_EDGE, 0.0);
        assertEquals("Z == 2R/sqrt(3), matching Jitterbug.java's translate()",
                     2.0 / Math.sqrt(3.0), Z, 1e-15);
        assertEquals(1.154701, Z, 1e-6);
    }

    /**
     * The sign convention is immaterial and must not be "fixed".
     *
     * <p>
     * These classes use sigma = +(sx*sy*sz), following the Python harness.
     * {@code PhiCoordinates.JITTERBUG_INVERSES} drives the JavaFX
     * {@code Jitterbug} with sigma = -(sx*sy*sz). The two families are the same
     * set of configurations traversed in opposite senses: the negated-sigma
     * family at a is the positive-sigma family at -a, to machine precision.
     */
    @Test
    public void negatedSigmaIsTheSameFamilyTraversedBackwards() {
        double worst = 0;
        for (double a : new double[] { 0, 10, 22.238756093, 30, 45, 60, 75, 90,
                                       120, 180 }) {
            worst = Math.max(worst, hausdorff(corners(a, false, -1.0),
                                              corners(-a, false, 1.0)));
        }
        assertEquals("set-Hausdorff between the two sign conventions", 0.0,
                     worst, 1e-14);

        // Control: at the SAME a they are genuinely different configurations,
        // so the agreement above is not vacuous.
        assertTrue("same-a control must not agree",
                   hausdorff(corners(30, false, -1.0),
                             corners(30, false, 1.0)) > 0.8);
    }

    /**
     * The 24 struts are geometrically distinct segments at a generic a and
     * coincide in pairs at a = 60, 90 and 120. a=90 joins the merge angles here
     * even though its vertices do not merge: cos(90)=0 puts every triangle's
     * centroid at the origin, so opposite triangles superpose.
     */
    @Test
    public void strutSegmentsCoincideInPairsAtSixtyNinetyAndOneTwenty() {
        assertEquals(24, distinctSegments(0.0));
        assertEquals(24, distinctSegments(30.0));
        assertEquals(12, distinctSegments(60.0));
        assertEquals(12, distinctSegments(90.0));
        assertEquals(12, distinctSegments(120.0));
        assertEquals(24, distinctSegments(180.0));
    }

    private static int distinctSegments(double a) {
        double[][][] x = corners(a);
        List<double[]> mids = new ArrayList<>();
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                double[] p = x[i][j];
                double[] q = x[i][(j + 1) % 3];
                mids.add(new double[] { 0.5 * (p[0] + q[0]),
                                        0.5 * (p[1] + q[1]),
                                        0.5 * (p[2] + q[2]) });
            }
        }
        List<double[]> reps = new ArrayList<>();
        for (double[] m : mids) {
            boolean seen = false;
            for (double[] r : reps) {
                if (distance(m, r) < MERGE_TOL) {
                    seen = true;
                    break;
                }
            }
            if (!seen) {
                reps.add(m);
            }
        }
        return reps.size();
    }

    private static double distance(double[] p, double[] q) {
        double dx = p[0] - q[0];
        double dy = p[1] - q[1];
        double dz = p[2] - q[2];
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    private static double hausdorff(double[][][] a, double[][][] b) {
        return Math.max(directed(a, b), directed(b, a));
    }

    private static double directed(double[][][] a, double[][][] b) {
        double worst = 0;
        for (double[][] fa : a) {
            for (double[] p : fa) {
                double best = Double.MAX_VALUE;
                for (double[][] fb : b) {
                    for (double[] q : fb) {
                        best = Math.min(best, distance(p, q));
                    }
                }
                worst = Math.max(worst, best);
            }
        }
        return worst;
    }

    /**
     * Sanity on the clustering primitive itself, which every other assertion in
     * this class routes through. A TreeSet keeps the assertion order-free.
     */
    @Test
    public void clusteringLabelsPartitionTheCorners() {
        JitterbugGeometry.Clustering c = cluster(corners(30.0), MERGE_TOL);
        TreeSet<Integer> seen = new TreeSet<>();
        int total = 0;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                seen.add(c.labels()[i][j]);
                total++;
            }
        }
        assertEquals(24, total);
        assertEquals(12, seen.size());
        assertEquals(Integer.valueOf(0), seen.first());
        assertEquals(Integer.valueOf(11), seen.last());
        int sum = 0;
        for (int m : c.multiplicities()) {
            sum += m;
        }
        assertEquals(24, sum);
    }
}
