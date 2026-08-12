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

import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.L_EDGE;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Reachability of the tetrahedron from the vector equilibrium — <b>through
 * interpenetrating configurations</b>.
 *
 * <p>
 * The qualifier is not decoration. Every demonstrated path here passes through
 * configurations in which the triangles pass through each other, at every step
 * refinement tried. Whether a collision-free path exists is <b>open</b>, and
 * nothing in this class may be cited as showing one does. "The tetrahedron is
 * reachable" without the qualifier is a claim this suite does not support.
 *
 * <p>
 * <b>Warrant — and the part of it that was withdrawn.</b> This class previously
 * claimed refinement invariance, "the same verdict at 40, 80, 160 and 200
 * waypoints", as what made a constructed path evidence. That claim is false and
 * is retracted: the per-target verdicts are <em>not</em> invariant under
 * refinement, and at least one target flips inside the very set of step counts
 * that was cited as proving invariance. Nothing may be cited from here as
 * showing refinement invariance.
 *
 * <p>
 * What is warranted is narrower and does survive: at every refinement run here
 * at least one seating is reached, and every waypoint carries a post-projection
 * hinge residual at machine precision, which is what puts the path <em>on</em>
 * the variety rather than near it. The residual alone would only cover the
 * waypoints, so it is necessary and not sufficient.
 *
 * <p>
 * <b>What is a property of the variety and what is a property of the search.</b>
 * "At least one tetrahedron is reachable" is the former. <em>Which</em> targets
 * are reached, and <em>how many</em>, is the latter — a joint function of the
 * solver and of the step count, and this suite therefore states no such count
 * anywhere. This class marches with an incremental analytic-Jacobian
 * Gauss-Newton; the Python reference drives scipy least-squares over a global
 * body-motion vector and reaches a different subset of the same eight. Neither
 * subset decomposes the variety into connected components, and an unreached
 * target is not evidence of a separate component — only of a march that got
 * stuck.
 *
 * <p>
 * <b>Why no count is pinned.</b> The march amplifies last-bit differences: in
 * the review recorded at T2 {@code inviscid/jitterbug-java-harness-review-findings.md},
 * perturbing a single corner coordinate by one ulp changed which targets were
 * reached. Any exact count out of this solver is reproducible only while every
 * arithmetic operation is bit-identical, which makes it a fact about a build
 * rather than about the linkage. The step counts here stop at 200 and nothing is
 * claimed beyond them.
 *
 * @author halhildebrand
 */
public class JitterbugReachabilityTest {

    private static final int    WAYPOINTS = 80;

    private static List<int[][]> TETS;
    private static double[][][]  VE;

    @BeforeClass
    public static void setUp() {
        TETS = TetrahedronSeatings.tetrahedra();
        VE = JitterbugGeometry.corners(0.0);
    }

    /**
     * The substantive claim, and the whole of it: at least one tetrahedron
     * seating is reachable from the vector equilibrium along a path that stays
     * on the constraint variety. How many of the eight sampled seatings this
     * solver happens to reach is deliberately not asserted — see the class
     * javadoc.
     */
    @Test
    public void theTetrahedronIsReachableFromTheVectorEquilibrium() {
        int reached = 0;
        for (int t = 0; t < 8; t++) {
            VarietyWalk.Result r = walkTo(t, WAYPOINTS);
            assertTrue("waypoint " + t + " left the variety: worst residual "
                       + r.worstResidual(), r.worstResidual() < 1e-10);
            if (r.reached()) {
                reached++;
                assertRegularTetrahedron(r.configuration(), "target " + t);
            } else {
                // A labelling-free check, because the RMS metric scores against
                // a FIXED corner correspondence: a walk landing on a genuine
                // tetrahedron with a different labelling would score 1/sqrt(3)
                // = 0.5774 and read as a failure. These are not that — they end
                // on 12-corner configurations, so they are real misses.
                assertEquals("target " + t + " scored a miss and must really be"
                             + " one, not a relabelled tetrahedron", 12,
                             JitterbugGeometry.cluster(r.configuration(), 1e-7)
                                              .representatives().length);
            }
        }
        assertTrue("at least one tetrahedron must be reachable", reached >= 1);
    }

    /**
     * What refinement actually buys, now that invariance has been withdrawn: the
     * existence claim holds at every step count tried, and every waypoint of
     * every march stays on the variety.
     *
     * <p>
     * Deliberately <b>not</b> asserted, because it is false or fitted: that a
     * given target keeps its verdict across refinements, that the number of
     * reached targets is stable, or that any particular endpoint RMS recurs. The
     * previous version of this test pinned one reached target and one unreached
     * target and had to truncate the unreached one's refinement list to keep its
     * RMS window true — the truncation was load-bearing, which is the tell.
     */
    @Test
    public void theExistenceClaimSurvivesEveryRefinement() {
        for (int n : new int[] { 40, 80, 160, 200 }) {
            int reached = 0;
            for (int t = 0; t < 8; t++) {
                VarietyWalk.Result r = walkTo(t, n);
                assertTrue("target " + t + " at n=" + n + " left the variety: "
                           + "worst residual " + r.worstResidual(),
                           r.worstResidual() < 1e-10);
                if (r.reached()) {
                    reached++;
                    assertRegularTetrahedron(r.configuration(),
                                             "target " + t + " at n=" + n);
                }
            }
            assertTrue("at least one tetrahedron must be reachable at n=" + n,
                       reached >= 1);
        }
    }

    /**
     * The demonstrated paths are <b>not collision-free</b>. Every refinement
     * interpenetrates on a large fraction of its waypoints, so the reachability
     * result must always be stated with that qualifier.
     *
     * <p>
     * Only the two robust facts are asserted: the count is never zero, and it is
     * never a negligible sliver of the path. The exact counts are not pinned and
     * neither is their trend with refinement — both move by up to a factor of two
     * under a last-bit perturbation of the corner coordinates, so a trend read
     * off this solver is not a property of the path.
     */
    @Test
    public void everyDemonstratedPathPassesThroughInterpenetration() {
        assertInterpenetrating(walkTo(1, 40));
        assertInterpenetrating(walkTo(1, 80));
        assertInterpenetrating(walkTo(1, 160));
        assertInterpenetrating(walkTo(1, 200));
    }

    /**
     * The Newton projection must actually project. Perturbing a configuration
     * off the variety and projecting has to bring the hinge residual back to
     * machine zero, and the projection must be a small correction rather than a
     * jump to somewhere unrelated.
     */
    @Test
    public void newtonProjectionReturnsToTheVariety() {
        java.util.Random rng = new java.util.Random(20260811L);
        double[] z = new double[48];
        for (int k = 0; k < 48; k++) {
            z[k] = 0.02 * (rng.nextDouble() - 0.5);
        }
        double[][][] off = JitterbugLinkage.applyBodyMotions(VE, z);
        double before = Linear.norm(JitterbugLinkage.residual(off));
        assertTrue("the perturbation must leave the variety, was " + before,
                   before > 1e-4);

        double[][][] on = VarietyWalk.project(off);
        assertTrue("projection must reach machine zero, got "
                   + Linear.norm(JitterbugLinkage.residual(on)),
                   Linear.norm(JitterbugLinkage.residual(on)) < 1e-13);

        double shift = 0;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                shift = Math.max(shift,
                                 JitterbugGeometry.distance(off[i][j],
                                                            on[i][j]));
            }
        }
        assertTrue("minimum-norm projection must be a small correction, was "
                   + shift, shift < 10 * before);
    }

    /**
     * Kabsch alignment must be a rigid motion and must not reflect: a mirrored
     * fit would let a walk "reach" the mirror image of its target and score it
     * as success.
     */
    @Test
    public void alignmentIsRigidAndNeverReflects() {
        double[][][] a = JitterbugGeometry.corners(0.0);
        double[][][] b = JitterbugGeometry.corners(37.0);
        double[][][] moved = VarietyWalk.align(a, b);
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals("alignment must preserve strut lengths", L_EDGE,
                             JitterbugGeometry.distance(moved[i][j],
                                                        moved[i][(j + 1) % 3]),
                             1e-12);
            }
        }
        // A mirror image cannot be aligned onto the original: the best rigid fit
        // leaves a residual, which is exactly what excluding reflections buys.
        double[][][] mirror = new double[8][3][3];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                mirror[i][j] = new double[] { -a[i][j][0], a[i][j][1],
                                              a[i][j][2] };
            }
        }
        double worst = 0;
        double[][][] fit = VarietyWalk.align(mirror, a);
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                worst = Math.max(worst,
                                 JitterbugGeometry.distance(fit[i][j],
                                                            a[i][j]));
            }
        }
        assertTrue("a reflection must not be absorbed by the alignment, worst "
                   + worst, worst > 0.1);
    }

    private static void assertInterpenetrating(VarietyWalk.Result r) {
        assertTrue("the path must interpenetrate somewhere at n="
                   + r.waypoints() + ", got " + r.interpenetrating(),
                   r.interpenetrating() > 0);
        assertTrue("a substantial fraction must interpenetrate at n="
                   + r.waypoints() + ", got " + r.interpenetrating(),
                   r.interpenetrating() > 0.3 * r.waypoints());
    }

    private static void assertRegularTetrahedron(double[][][] x, String what) {
        JitterbugGeometry.Clustering c = JitterbugGeometry.cluster(x, 1e-7);
        assertEquals(what + " must have 4 distinct corners", 4,
                     c.representatives().length);
        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                assertEquals(what + " must be regular with edge L", L_EDGE,
                             JitterbugGeometry.distance(c.representatives()[i],
                                                        c.representatives()[j]),
                             1e-9);
            }
        }
        assertEquals(what + " must satisfy the hinges", 0.0,
                     Linear.norm(JitterbugLinkage.residual(x)), 1e-12);
    }

    private static VarietyWalk.Result walkTo(int target, int waypoints) {
        return VarietyWalk.walk(VE, TetrahedronSeatings.build(TETS.get(target)),
                                waypoints);
    }
}
