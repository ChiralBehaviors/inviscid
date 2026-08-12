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

import static com.chiralbehaviors.inviscid.jitterbug.Interpenetration.count;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.corners;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

/**
 * The collision boundary.
 *
 * <p>
 * This is the finding nobody was looking for, and it reframes {@code a = 60}
 * after that angle's claimed status as a rank singularity was refuted. The
 * model self-intersects on {@code (60, 120)} and on {@code (240, 300)} in the
 * {@code eps -> 0} limit, and nowhere else.
 *
 * <p>
 * These tests measure the interference; they do not rule on it. Under USER
 * DECISION 16 (2026-08-12, bead {@code inviscid-qvf.8}) interference is
 * <b>permitted</b>, so the symmetric path continues through both intervals and
 * the configuration space is the whole variety. An earlier version of this
 * javadoc concluded that "the admissible symmetric path is the interval, not
 * the circle" and that "there is no 2*pi circuit for a periodicity argument to
 * run on"; that conclusion is <b>withdrawn</b>. It was a verdict on the
 * measurement rather than the measurement, and every count asserted below is
 * untouched by its withdrawal.
 *
 * <p>
 * (The hbar periodicity argument stays dead regardless — it fails on three
 * further independent counts that have nothing to do with interference.)
 *
 * <p>
 * At the finite margin the counter actually runs at, the endpoints sit a little
 * inside those intervals, and the tests below probe there rather than sampling
 * only where the limit statement happens to hold.
 *
 * @author halhildebrand
 */
public class InterpenetrationTest {

    /**
     * Item 5, with the endpoints probed rather than stepped over. The count is
     * 0 or 24 and nothing else; it is 0 well below 60 and well above 120, and 24
     * across the interior.
     *
     * <p>
     * <b>The universal "24 throughout {@code 60 < a < 120}" is not asserted,
     * because it is false at the margin this counter runs at.</b> The angles
     * 60.0005 and 119.9989 lie strictly inside the limiting interval and count
     * <em>zero</em>: onset is tangential, so the finite margin pushes both
     * endpoints inward by about a thousandth of a degree. The earlier version of
     * this test sampled 59.9 and then jumped to 61, which probed the universal
     * only where it holds — structurally the same slice-then-count defect this
     * port was written to correct, so the inward angles are now sampled
     * explicitly and on both ends.
     */
    @Test
    public void theCountIsZeroOrTwentyFourAndTheEndpointsAreMarginLimited() {
        for (double a : new double[] { 0, 15, 30, 45, 55, 59, 59.9, 60 }) {
            assertEquals("a=" + a + " must be collision-free", 0,
                         count(corners(a)));
        }
        // Strictly inside (60,120) and strictly inside (240,300), and still 0:
        // this is what refutes the universal.
        for (double a : new double[] { 60.0005, 60.001, 119.9989, 240.0005,
                                       299.9989 }) {
            assertEquals("a=" + a + " is inside the limiting interval but the"
                         + " finite margin has not caught up yet", 0,
                         count(corners(a)));
        }
        for (double a : new double[] { 60.002, 60.01, 60.1, 60.5, 61, 70, 80, 90,
                                       100, 110, 119, 119.9, 119.99 }) {
            assertEquals("a=" + a + " must show 24 interpenetrations", 24,
                         count(corners(a)));
        }
        for (double a : new double[] { 120, 120.1, 130, 150, 180 }) {
            assertEquals("a=" + a + " must be collision-free again", 0,
                         count(corners(a)));
        }
    }

    /**
     * Item 5, over the whole circle. A 1 degree sweep of [0,360) sees only the
     * values 0 and 24, and on that grid the collision set is the two open
     * intervals (60,120) and (240,300).
     *
     * <p>
     * A 1 degree grid cannot resolve the margin-limited endpoints — the nearest
     * sample above 60 is 61, a thousand times further out than the transition.
     * That is exactly why the sweep is not the whole story and
     * {@link #theCountIsZeroOrTwentyFourAndTheEndpointsAreMarginLimited()}
     * probes inside the interval. What the sweep does establish, and it is worth
     * having, is that no third value occurs anywhere.
     */
    @Test
    public void theCollisionSetIsExactlyTwoOpenIntervals() {
        List<String> transitions = new ArrayList<>();
        int previous = -1;
        for (int k = 0; k < 360; k++) {
            int n = count(corners(k));
            assertTrue("only 0 and 24 may occur; a=" + k + " gave " + n,
                       n == 0 || n == 24);
            boolean inside = (k > 60 && k < 120) || (k > 240 && k < 300);
            assertEquals("collision membership at a=" + k, inside ? 24 : 0, n);
            if (n != previous) {
                transitions.add(k + ":" + n);
                previous = n;
            }
        }
        assertEquals("[0:0, 61:24, 120:0, 241:24, 300:0]",
                     transitions.toString());
    }

    /**
     * Item 5, onset. The reference recorded the onset as "precisely at the
     * octahedron" while measuring it at 60.0011 with its default margin. The
     * first is the {@code eps -> 0} limit and the second is the measurement;
     * they are not interchangeable, and the second explains the first. Onset is
     * <em>tangential</em>, so the observed onset angle is limited by the
     * strict-interior margin and scales as sqrt(eps). Tightening eps by 10^3
     * moves the onset 10^1.5 closer to 60 — measured here across three decades,
     * which is what turns "onset at 60" from an assertion into an extrapolation
     * with a measured rate.
     */
    @Test
    public void onsetIsAtTheOctahedronAndTheResidualGapIsToleranceLimited() {
        double coarse = bisectOnset(1e-9);
        double medium = bisectOnset(1e-12);
        double fine = bisectOnset(1e-14);

        assertEquals("onset at eps=1e-9", 60.001109528, coarse, 1e-7);
        assertEquals("onset at eps=1e-12", 60.000035090, medium, 1e-8);
        assertTrue("onset at eps=1e-14 is within 1e-5 of 60, was " + fine,
                   Math.abs(fine - 60.0) < 1e-5);

        assertTrue("onset must approach 60 monotonically from above",
                   60.0 < fine && fine < medium && medium < coarse);

        // sqrt scaling: three decades tighter in eps buys 10^1.5 in the gap.
        double ratio = (coarse - 60.0) / (medium - 60.0);
        assertEquals("gap scales as sqrt(eps)", Math.sqrt(1000.0), ratio, 2.0);
    }

    /**
     * The hinge exclusion is <b>inert on the symmetric family</b> — and that is
     * worth knowing, not hiding.
     *
     * <p>
     * Removing it changes the count at no angle on the whole circle. The reason
     * is that a shared hinge vertex is always an <em>endpoint</em> of the
     * piercing edge, and the strict ray test already requires
     * {@code t} strictly inside {@code (0,1)}; the exclusion would only bite if
     * a shared vertex fell in the far triangle's interior, which never happens
     * here. So the collision boundary at {@code a = 60} is not an artefact of
     * the contact rule: it survives with the rule switched off entirely. The
     * rule is still correct to keep for off-path configurations, and
     * {@link #theHingeExclusionSuppressesAContactAtASharedVertex()} tests the
     * mechanism directly.
     */
    @Test
    public void theHingeExclusionChangesNothingOnTheSymmetricPath() {
        for (int k = 0; k < 360; k++) {
            double[][][] x = corners(k);
            int naive = 0;
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 8; j++) {
                    if (i != j) {
                        naive += Interpenetration.edgesThrough(x[i], x[j],
                                                               List.of(),
                                                               Interpenetration.DEFAULT_EPS);
                    }
                }
            }
            assertEquals("hinge exclusion must not matter at a=" + k, count(x),
                         naive);
        }
    }

    /**
     * The exclusion mechanism itself, exercised on the only geometry that can
     * reach it: an edge piercing a triangle's interior exactly at a vertex the
     * two bodies share.
     */
    @Test
    public void theHingeExclusionSuppressesAContactAtASharedVertex() {
        double[][] face = { { -1, -1, 0 }, { 1, -1, 0 }, { 0, 1, 0 } };
        double[][] spike = { { 0, 0, -1 }, { 0, 0, 1 }, { 0.8, 0.8, 0.5 } };
        double[] origin = { 0, 0, 0 };

        assertEquals("the spike pierces the face at the origin", 1,
                     Interpenetration.edgesThrough(spike, face, List.of(),
                                                   Interpenetration.DEFAULT_EPS));
        assertEquals("declaring the origin a shared hinge vertex suppresses it",
                     0, Interpenetration.edgesThrough(spike, face,
                                                      List.of(origin),
                                                      Interpenetration.DEFAULT_EPS));
        assertEquals("an unrelated shared vertex must not suppress anything", 1,
                     Interpenetration.edgesThrough(spike, face,
                                                   List.of(new double[] { 0.9,
                                                                          -0.9,
                                                                          0 }),
                                                   Interpenetration.DEFAULT_EPS));
    }

    /**
     * The counter must be able to see an intersection at all. Two triangles
     * pushed through each other have to register, otherwise the zeros above mean
     * nothing.
     */
    @Test
    public void theCounterDetectsAnActualCrossing() {
        double[][] a = { { -1, -1, 0 }, { 1, -1, 0 }, { 0, 1, 0 } };
        double[][] b = { { 0, 0, -1 }, { 0.5, 0, 1 }, { -0.5, 0, 1 } };
        assertEquals("b's edges must pierce a", 2,
                     Interpenetration.edgesThrough(b, a, List.of(),
                                                   Interpenetration.DEFAULT_EPS));
        assertEquals("a's edges do not pierce b here", 0,
                     Interpenetration.edgesThrough(a, b, List.of(),
                                                   Interpenetration.DEFAULT_EPS));

        // Separated: nothing.
        double[][] far = { { 0, 0, 9 }, { 0.5, 0, 11 }, { -0.5, 0, 11 } };
        assertEquals(0, Interpenetration.edgesThrough(far, a, List.of(),
                                                      Interpenetration.DEFAULT_EPS));
    }

    /**
     * The shared-corner lookup must find exactly the hinge that joins two
     * adjacent faces, and nothing for non-adjacent ones. The interpenetration
     * count is only as trustworthy as this.
     */
    @Test
    public void sharedCornerLookupMatchesTheHingeGraph() {
        double[][][] x = corners(30.0);
        int adjacent = 0;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                if (i == j) {
                    continue;
                }
                int n = Interpenetration.sharedCorners(x, i, j).size();
                assertTrue("faces share at most one corner, " + i + "/" + j
                           + " shared " + n, n <= 1);
                if (n == 1) {
                    adjacent++;
                }
            }
        }
        assertEquals("12 hinges seen from both ends", 24, adjacent);
    }

    /**
     * Bisect the first angle above 60 at which the count becomes nonzero.
     */
    private static double bisectOnset(double eps) {
        double lo = 60.0;
        double hi = 60.5;
        for (int i = 0; i < 50; i++) {
            double mid = 0.5 * (lo + hi);
            if (count(corners(mid), eps) == 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return hi;
    }
}
