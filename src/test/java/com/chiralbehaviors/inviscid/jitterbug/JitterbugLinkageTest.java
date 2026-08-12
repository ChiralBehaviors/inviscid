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

import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.corners;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugLinkage.hingePairs;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugLinkage.jacobian;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugLinkage.residual;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

/**
 * The jitterbug as a LINKAGE: eight rigid triangles (48 body DOF) tied together
 * at twelve shared vertices (36 scalar constraints).
 *
 * <p>
 * This class pins the results that overturned the recorded "1 degree of freedom,
 * bifurcating at a=60" model. Measured here and independently in
 * {@code analysis/jitterbug-variety/jb_b_variety.py} and {@code jb_c_branches.py}:
 * the internal DOF count is <b>6</b>, and a=60 is <b>not</b> a rank singularity —
 * the only rank drops on the whole circle are at a=90 and a=270.
 *
 * @author halhildebrand
 */
public class JitterbugLinkageTest {

    private static final double RANK_TOL = 1e-9;

    /**
     * Item 6, and the harness's own guard. The 36x48 Jacobian must be the exact
     * derivative of the nonlinear hinge residual with respect to the body
     * motions. This finite-difference cross-check is what caught a sign error in
     * the {@code omega} block of the reference implementation, so it is ported
     * as a test rather than left as a one-off script.
     */
    @Test
    public void jacobianAgreesWithFiniteDifferencesOfTheExactResidual() {
        double worst = 0;
        for (double a : new double[] { 0, 17, 30, 60, 90, 137, 271 }) {
            worst = Math.max(worst, maxJacobianDeviation(corners(a),
                                                         jacobian(corners(a))));
        }
        assertTrue("finite-difference agreement (measured max 1.1e-10), got "
                   + worst, worst < 1e-9);
    }

    /**
     * The finite-difference check above must be able to fail. Negating the
     * angular-velocity block — the exact class of error it was written to catch
     * — has to blow it up by nine orders of magnitude.
     */
    @Test
    public void theFiniteDifferenceCheckIsFalsifiable() {
        double[][][] x = corners(30.0);
        double[][] corrupted = jacobian(x);
        for (double[] row : corrupted) {
            for (int c = 0; c < 24; c++) {
                row[c] = -row[c];
            }
        }
        double deviation = maxJacobianDeviation(x, corrupted);
        assertTrue("a sign-flipped omega block must be rejected, deviation was "
                   + deviation, deviation > 0.1);
    }

    /**
     * Item 7. Rank is 36 through a=60 — sigma_36 stays at 0.5987572, three
     * ulp-scale steps either side — so a=60 is a smooth point of the variety and
     * cannot be a branch point. Constant rank on a neighbourhood means locally a
     * manifold.
     */
    @Test
    public void aSixtyIsNotARankSingularity() {
        for (double a : new double[] { 59.999999, 60.0, 60.000001 }) {
            double[] s = Linear.singularValues(jacobian(corners(a)));
            assertEquals("rank at a=" + a, 36,
                         Linear.rank(jacobian(corners(a)), RANK_TOL));
            assertEquals("sigma_36 at a=" + a, 0.5987572, s[35], 1e-7);
        }
    }

    /**
     * Item 7, the whole circle. A 0.1 degree sweep of [0,360) finds exactly two
     * rank drops, at a=90 and a=270, where sigma_36 collapses to ~1e-16. Nothing
     * else on the circle drops. Generic internal DOF = nullity - 6 global = 6.
     */
    @Test
    public void theOnlyRankDropsOnTheCircleAreAtNinetyAndTwoSeventy() {
        List<Double> drops = new ArrayList<>();
        for (int k = 0; k < 3600; k++) {
            double a = k * 0.1;
            double[][] j = jacobian(corners(a));
            int rank = Linear.rank(j, RANK_TOL);
            if (rank != 36) {
                assertEquals("rank drop at a=" + a + " must be to 35", 35, rank);
                drops.add(Math.round(a * 10.0) / 10.0);
            } else {
                assertEquals("internal DOF at a=" + a, 6, 48 - rank - 6);
            }
        }
        assertEquals("rank drops on a 0.1 degree sweep of [0,360)",
                     List.of(90.0, 270.0), drops);

        for (double a : new double[] { 90.0, 270.0 }) {
            double[] s = Linear.singularValues(jacobian(corners(a)));
            assertTrue("sigma_36 must collapse at a=" + a + ", was " + s[35],
                       s[35] < 1e-14);
            assertTrue("sigma_35 must NOT collapse at a=" + a + " (the drop is "
                       + "by exactly one), was " + s[34], s[34] > 0.1);
            assertEquals("local dimension at a=" + a, 7, 48 - 35 - 6);
        }
    }

    /**
     * Item 8. The six global rigid motions — three translations and three
     * rotations about the origin, the latter carrying
     * {@code dt_i = omega x c_i} — must lie in the Jacobian's null space
     * exactly. This is a real guard: any sign or transpose error in the
     * {@code omega} block breaks it loudly.
     */
    @Test
    public void globalRigidMotionsAreInTheNullSpace() {
        for (double a : new double[] { 0, 30, 60, 90, 137, 270 }) {
            double[][][] x = corners(a);
            double[][] j = jacobian(x);
            double[][] g = JitterbugLinkage.globalRigidMotions(x);
            assertEquals(6, g.length);
            for (int k = 0; k < 6; k++) {
                double n = Linear.norm(Linear.apply(j, g[k]));
                assertTrue("global motion " + k + " at a=" + a + " gave ||J g||="
                           + n, n < 1e-9);
            }
        }
    }

    /**
     * The companion falsifiability check for the one above: the rotation rows
     * only satisfy J because of the {@code omega x c_i} translation coupling.
     * Drop that coupling and they must stop being null directions, otherwise the
     * assertion above is passing for the wrong reason.
     */
    @Test
    public void globalRotationsNeedTheirTranslationCoupling() {
        double[][][] x = corners(30.0);
        double[][] j = jacobian(x);
        double[][] g = JitterbugLinkage.globalRigidMotions(x);
        double worst = 0;
        for (int k = 3; k < 6; k++) {
            double[] stripped = g[k].clone();
            for (int c = 24; c < 48; c++) {
                stripped[c] = 0;
            }
            worst = Math.max(worst, Linear.norm(Linear.apply(j, stripped)));
        }
        assertTrue("stripping omega x c must break the null-space property",
                   worst > 0.1);
    }

    /**
     * The linkage's combinatorics are read off once at a generic angle and then
     * held fixed. At a=60 the twelve shared vertices merge to six, but the hinge
     * incidences do not change — the merge is a coincidence of positions, not a
     * re-wiring. Re-deriving the pairing at a=60 would find a 6x4 clustering and
     * silently invent a different linkage.
     */
    @Test
    public void hingePairingIsFixedAndSurvivesTheMerge() {
        int[][] pairs = hingePairs();
        assertEquals(12, pairs.length);

        Set<Integer> covered = new HashSet<>();
        for (int[] p : pairs) {
            assertEquals(4, p.length);
            assertTrue("a hinge must join two DIFFERENT faces", p[0] != p[2]);
            covered.add(p[0] * 3 + p[1]);
            covered.add(p[2] * 3 + p[3]);
        }
        assertEquals("the 12 hinges use 24 distinct (face,corner) slots", 24,
                     covered.size());

        // The pairing satisfies the constraints identically along the whole
        // symmetric path, including through the merge and past it.
        for (double a : new double[] { 0, 22.238756093, 45, 59.9, 60, 60.1, 90,
                                       120, 180, 240, 300 }) {
            assertEquals("hinge residual at a=" + a, 0.0,
                         Linear.norm(residual(corners(a))), 1e-14);
        }
    }

    /**
     * Item 11. The hinge incidence graph on the eight triangles is the cube
     * graph Q3: 12 edges, 3-regular, triangle-free, and — the structural
     * statement, not just the numerology — two faces are adjacent exactly when
     * their sign triples differ in one coordinate.
     */
    @Test
    public void hingeIncidenceGraphIsTheCubeGraphQ3() {
        Set<Long> edges = new HashSet<>();
        for (int[] p : hingePairs()) {
            edges.add((long) Math.min(p[0], p[2]) * 8 + Math.max(p[0], p[2]));
        }
        assertEquals("Q3 has 12 edges", 12, edges.size());

        int[] degree = new int[8];
        for (long e : edges) {
            degree[(int) (e / 8)]++;
            degree[(int) (e % 8)]++;
        }
        for (int i = 0; i < 8; i++) {
            assertEquals("Q3 is 3-regular; face " + i, 3, degree[i]);
        }

        int triangles = 0;
        for (int a = 0; a < 8; a++) {
            for (int b = a + 1; b < 8; b++) {
                for (int c = b + 1; c < 8; c++) {
                    if (edges.contains((long) a * 8 + b)
                        && edges.contains((long) b * 8 + c)
                        && edges.contains((long) a * 8 + c)) {
                        triangles++;
                    }
                }
            }
        }
        assertEquals("Q3 is triangle-free", 0, triangles);

        for (long e : edges) {
            assertEquals("adjacent faces differ in exactly one sign", 1,
                         hamming((int) (e / 8), (int) (e % 8)));
        }
        // ...and the converse, so this really is Q3 and not a subgraph of it.
        for (int a = 0; a < 8; a++) {
            for (int b = a + 1; b < 8; b++) {
                boolean adjacent = edges.contains((long) a * 8 + b);
                assertEquals("Hamming-1 pair (" + a + "," + b + ") adjacency",
                             hamming(a, b) == 1, adjacent);
            }
        }
    }

    /**
     * The body-motion parameterisation the Jacobian's columns refer to: face
     * {@code i} gets a rotation about its <em>own current centroid</em> plus a
     * translation. Zero motion must be the identity, and the motions must be
     * genuinely rigid.
     */
    @Test
    public void bodyMotionsAreRigidAndIdentityAtZero() {
        double[][][] x0 = corners(30.0);
        double[][][] same = JitterbugLinkage.applyBodyMotions(x0, new double[48]);
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                // Tolerance 1e-15, not 0.0. Zero motion computes (x - c) + c,
                // which round-trips bit-exactly at a=30 by luck rather than by
                // any property of applyBodyMotions: 24 of the 72 one-ULP
                // coordinate perturbations break exact equality, by 2.5e-32 to
                // 4.9e-32. Asserting 0.0 pins the luck and would flake on other
                // hardware. 1e-15 is ~17 orders above the observed drift and
                // still ~13 below the strut length, so a real identity failure
                // cannot hide under it.
                assertEquals(0.0, JitterbugGeometry.distance(x0[i][j],
                                                             same[i][j]),
                             1e-15);
            }
        }

        java.util.Random rng = new java.util.Random(20260811L);
        double[] z = new double[48];
        for (int k = 0; k < 48; k++) {
            z[k] = 0.7 * (rng.nextDouble() - 0.5);
        }
        double[][][] moved = JitterbugLinkage.applyBodyMotions(x0, z);
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals("body " + i + " must move rigidly",
                             JitterbugGeometry.distance(x0[i][j],
                                                        x0[i][(j + 1) % 3]),
                             JitterbugGeometry.distance(moved[i][j],
                                                        moved[i][(j + 1) % 3]),
                             1e-13);
            }
        }
        // A generic body motion tears the hinges apart; if it did not, the
        // constraint set would be everything and the rank results meaningless.
        assertFalse("a random body motion must violate the hinges",
                    Linear.norm(residual(moved)) < 1.0);
    }

    /**
     * The constraint Jacobian must be the hinge-pair row difference of the
     * position Jacobian. The two are written independently, so this catches a
     * transcription slip in either without either one validating itself.
     */
    @Test
    public void constraintJacobianIsTheRowDifferenceOfThePositionJacobian() {
        for (double a : new double[] { 0, 30, 60, 90 }) {
            double[][][] x = corners(a);
            double[][] j = jacobian(x);
            double[][] p = JitterbugLinkage.positionJacobian(x);
            int[][] pairs = hingePairs();
            for (int c = 0; c < 12; c++) {
                for (int d = 0; d < 3; d++) {
                    int ra = 9 * pairs[c][0] + 3 * pairs[c][1] + d;
                    int rb = 9 * pairs[c][2] + 3 * pairs[c][3] + d;
                    for (int col = 0; col < 48; col++) {
                        assertEquals("row " + (3 * c + d) + " col " + col
                                     + " at a=" + a, p[ra][col] - p[rb][col],
                                     j[3 * c + d][col], 1e-15);
                    }
                }
            }
        }
    }

    private static int hamming(int faceA, int faceB) {
        double[] a = JitterbugGeometry.faceAxis(faceA);
        double[] b = JitterbugGeometry.faceAxis(faceB);
        int d = 0;
        for (int k = 0; k < 3; k++) {
            if (Math.signum(a[k]) != Math.signum(b[k])) {
                d++;
            }
        }
        return d;
    }

    /**
     * Central-difference the exact nonlinear residual of the body-motion
     * parameterisation at {@code z = 0} and report the worst disagreement with
     * {@code claimed}.
     */
    private static double maxJacobianDeviation(double[][][] x0,
                                               double[][] claimed) {
        double h = 1e-6;
        double worst = 0;
        for (int c = 0; c < 48; c++) {
            double[] zp = new double[48];
            double[] zm = new double[48];
            zp[c] = h;
            zm[c] = -h;
            double[] rp = residual(JitterbugLinkage.applyBodyMotions(x0, zp));
            double[] rm = residual(JitterbugLinkage.applyBodyMotions(x0, zm));
            for (int r = 0; r < 36; r++) {
                double fd = (rp[r] - rm[r]) / (2 * h);
                worst = Math.max(worst, Math.abs(claimed[r][c] - fd));
            }
        }
        return worst;
    }
}
