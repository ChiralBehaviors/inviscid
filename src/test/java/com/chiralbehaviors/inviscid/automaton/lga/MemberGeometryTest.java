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

import static com.chiralbehaviors.inviscid.Constants.ROOT_2;
import static com.chiralbehaviors.inviscid.Constants.TWO_PI;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import javax.vecmath.Vector3d;

import org.junit.Test;

import com.chiralbehaviors.inviscid.LengthTable;
import com.chiralbehaviors.inviscid.PhiCoordinates;

/**
 * Headless behavioral tests for {@link MemberGeometry}, independent of the
 * JavaFX composition it ports (see {@link MemberGeometryJavaFxParityTest}
 * for the port-verification anchor).
 *
 * @author halhildebrand
 */
public class MemberGeometryTest {

    private static final double DELTA      = 1e-9;
    private static final int    RESOLUTION = 360;

    /**
     * Same convention {@code MemberGeometry.memberTip} uses: even members
     * expose the segment's "b" endpoint as tip, odd members expose "a".
     */
    private static Vector3d tipOf(Segment segment, int member) {
        return (member % 2 == 0) ? segment.getB() : segment.getA();
    }

    @Test
    public void geometryIsPureFunctionOfAngle() {
        MemberGeometry geometry = new MemberGeometry(RESOLUTION, LgaTestGeometry.BASELINE_RADIUS);

        for (int cube = 0; cube < 5; cube++) {
            for (int member = 0; member < 6; member++) {
                for (float angle = 0f; angle < TWO_PI; angle += 0.37f) {
                    Segment first = geometry.memberSegment(cube, member, angle);
                    Segment second = geometry.memberSegment(cube, member, angle);
                    assertEquals("a.x differs across repeated calls",
                                first.getA().x, second.getA().x, 0.0);
                    assertEquals("a.y differs across repeated calls",
                                first.getA().y, second.getA().y, 0.0);
                    assertEquals("a.z differs across repeated calls",
                                first.getA().z, second.getA().z, 0.0);
                    assertEquals("b.x differs across repeated calls",
                                first.getB().x, second.getB().x, 0.0);
                    assertEquals("b.y differs across repeated calls",
                                first.getB().y, second.getB().y, 0.0);
                    assertEquals("b.z differs across repeated calls",
                                first.getB().z, second.getB().z, 0.0);
                }
            }
        }

        // No shared mutable state: interleaving calls for different
        // (cube, member, angle) must not perturb an already-computed result.
        Segment baseline = geometry.memberSegment(1, 3, 1.234f);
        for (int cube = 0; cube < 5; cube++) {
            for (int member = 0; member < 6; member++) {
                geometry.memberSegment(cube, member, 5.5f);
            }
        }
        Segment afterInterleaving = geometry.memberSegment(1, 3, 1.234f);
        assertEquals(baseline.getA().x, afterInterleaving.getA().x, 0.0);
        assertEquals(baseline.getA().y, afterInterleaving.getA().y, 0.0);
        assertEquals(baseline.getA().z, afterInterleaving.getA().z, 0.0);
        assertEquals(baseline.getB().x, afterInterleaving.getB().x, 0.0);
        assertEquals(baseline.getB().y, afterInterleaving.getB().y, 0.0);
        assertEquals(baseline.getB().z, afterInterleaving.getB().z, 0.0);
    }

    @Test
    public void lengthModulationMatchesLengthTable() {
        MemberGeometry geometry = new MemberGeometry(RESOLUTION, LgaTestGeometry.BASELINE_RADIUS);
        LengthTable table = new LengthTable(RESOLUTION);
        double edgeLength = PhiCoordinates.Cubes[0].getEdgeLength();
        double angularResolution = TWO_PI / RESOLUTION;
        // full segment span = 2 * halfSegmentLength * length
        //                   = edgeLength * ROOT_2 * length
        double scale = edgeLength * ROOT_2;

        for (int cube = 0; cube < 5; cube++) {
            for (int member = 0; member < 6; member++) {
                for (int step = 0; step < RESOLUTION; step++) {
                    // sample mid-step to stay clear of boundary rounding
                    float angle = (float) ((step + 0.5) * angularResolution);
                    Segment segment = geometry.memberSegment(cube, member, angle);
                    double expected = table.lengthAt(step) * scale;
                    assertEquals("cube " + cube + " member " + member
                                + " step " + step, expected,
                                segment.length(), 1e-6);
                }
            }
        }
    }

    @Test
    public void oppositeMembersAreAntipodalAtZeroAngle() {
        MemberGeometry geometry = new MemberGeometry(RESOLUTION, LgaTestGeometry.BASELINE_RADIUS);
        int[][] pairs = { { 0, 1 }, { 2, 3 }, { 4, 5 } };

        for (int cube = 0; cube < 5; cube++) {
            for (int[] pair : pairs) {
                Segment segmentEven = geometry.memberSegment(cube, pair[0], 0f);
                Segment segmentOdd = geometry.memberSegment(cube, pair[1], 0f);
                Vector3d tipEven = tipOf(segmentEven, pair[0]);
                Vector3d tipOdd = tipOf(segmentOdd, pair[1]);

                assertEquals("cube " + cube + " pair " + pair[0] + "/"
                            + pair[1] + " x", -tipEven.x, tipOdd.x, DELTA);
                assertEquals("cube " + cube + " pair " + pair[0] + "/"
                            + pair[1] + " y", -tipEven.y, tipOdd.y, DELTA);
                assertEquals("cube " + cube + " pair " + pair[0] + "/"
                            + pair[1] + " z", -tipEven.z, tipOdd.z, DELTA);
            }
        }
    }

    /**
     * Verifies the tip structure the 3 member pairs actually produce at
     * angle 0.
     * <p>
     * Bead-spec correction (code review, tracked in inviscid-18i and bead
     * notes on inviscid-0nx.2): the bead's original test name -
     * {@code sixMembersSpanThreeOrthogonalAxesAtZero} - claimed all three
     * representative tips are mutually orthogonal. That is FALSE for this
     * geometry; asserting it would be exactly the vacuous/wrong-assertion
     * failure mode this bead exists to prevent, so this test was renamed to
     * match what is actually, verifiably true instead of what the bead
     * name assumed. What genuinely, algebraically holds (derived from the
     * construction, not curve-fit to observed numbers):
     * <p>
     * At step 0, {@code LengthTable.lengthAt(0) == 1/ROOT_2} exactly (the
     * angle-0 ray hits the fundamental octant's bounding line at
     * {@code x = halfEdge = 1/ROOT_2}, see {@code LengthTable}'s
     * constructor). Since {@code halfSegmentLength == edgeLength*ROOT_2/2}
     * and {@code halfInterval == edgeLength/2}, that makes
     * {@code halfSegmentLength * lengthAt(0) == halfInterval} exactly - the
     * member's scaled half-arm length exactly equals the cell's translate
     * offset at angle 0. That identity is what forces:
     * <ul>
     * <li>member2's tip (rotation about X, offset along X) and member4's
     * tip (the same offset, but rotation pre/post-composed with a fixed
     * 90-degree Z rotation) to be exactly perpendicular - the X-axis
     * rotation family and its Z-reoriented sibling never cancel and never
     * reinforce along the same component.</li>
     * <li>member0's tip (rotation about Z, offset along Z - a different
     * pair entirely) to sit at 60 degrees (dot / (|a||b|) == 0.5) from
     * both member2's and member4's tips, by the same magnitude, rather
     * than at 90 degrees.</li>
     * <li>all three representative tips to be equal in magnitude, and
     * linearly independent (nonzero scalar triple product) - the six
     * members genuinely span 3 dimensions, they don't collapse onto a
     * plane or a line.</li>
     * </ul>
     */
    @Test
    public void memberAxesFormVerifiedStructureAtZero() {
        MemberGeometry geometry = new MemberGeometry(RESOLUTION, LgaTestGeometry.BASELINE_RADIUS);

        for (int cube = 0; cube < 5; cube++) {
            Vector3d t0 = geometry.memberTip(cube, 0, 0f);
            Vector3d t2 = geometry.memberTip(cube, 2, 0f);
            Vector3d t4 = geometry.memberTip(cube, 4, 0f);

            double len0 = t0.length();
            double len2 = t2.length();
            double len4 = t4.length();
            assertEquals("cube " + cube + " |t0| != |t2|", len0, len2, DELTA);
            assertEquals("cube " + cube + " |t0| != |t4|", len0, len4, DELTA);

            // member2 <-> member4: exactly orthogonal.
            assertEquals("cube " + cube + " t2 not orthogonal to t4", 0.0,
                        t2.dot(t4), DELTA);

            // member0 <-> member2 / member4: 60 degrees, equal on both.
            double cos0v2 = t0.dot(t2) / (len0 * len2);
            double cos0v4 = t0.dot(t4) / (len0 * len4);
            assertEquals("cube " + cube + " angle(t0,t2) != 60deg", 0.5,
                        cos0v2, 1e-6);
            assertEquals("cube " + cube + " angle(t0,t4) != 60deg", 0.5,
                        cos0v4, 1e-6);

            // Linearly independent: nonzero scalar triple product, i.e.
            // the six members' tips genuinely span three dimensions.
            Vector3d cross = new Vector3d();
            cross.cross(t2, t4);
            double tripleProduct = t0.dot(cross);
            assertTrue("cube " + cube + " t0,t2,t4 are coplanar/degenerate: "
                      + tripleProduct, Math.abs(tripleProduct) > 1e-3);
        }
    }

    @Test
    public void memberRadiusIsConstructorProvidedValue() {
        MemberGeometry geometry = new MemberGeometry(RESOLUTION, LgaTestGeometry.BASELINE_RADIUS);
        assertEquals(LgaTestGeometry.BASELINE_RADIUS, geometry.memberRadius(), 0.0);

        MemberGeometry other = new MemberGeometry(RESOLUTION, 0.5);
        assertEquals(0.5, other.memberRadius(), 0.0);
    }

    /**
     * Code review Important (inviscid-0nx.2 FIX 4b): {@code stepOf}'s
     * exact-boundary rounding ({@code (int) (normalized / angularResolution)}
     * at a point where the true mathematical quotient is an integer) is
     * currently only safe because {@code Constants.TWO_PI} and {@code
     * Constants.ROOT_2} are {@code float}, not {@code double} - a boundary
     * angle built from those floats and re-divided lands reliably on one
     * side of the integer rather than jittering across it from rounding
     * noise. That was an implicit precision coincidence, not a tested
     * contract. This test pins the two boundary behaviors that matter -
     * exact step boundaries, and the {@code angle == 2*PI} wrap - so a
     * future change to that arithmetic (e.g. switching to double angular
     * resolution, or reordering the modulo) shows up as a failing test
     * instead of a silent off-by-one at render/contact-detection
     * boundaries.
     */
    @Test
    public void stepBoundaryAndWrapAreLockedContract() {
        MemberGeometry geometry = new MemberGeometry(RESOLUTION, LgaTestGeometry.BASELINE_RADIUS);
        LengthTable table = new LengthTable(RESOLUTION);
        double edgeLength = PhiCoordinates.Cubes[0].getEdgeLength();
        double angularResolution = TWO_PI / RESOLUTION;
        double scale = edgeLength * ROOT_2;

        // angle == 2*PI must wrap to exactly the same geometry as angle ==
        // 0 - float modulo of an exact multiple of the modulus is 0, not
        // the modulus itself, so this exercises the "normalized < 0"
        // branch never firing here, landing on step 0, not step
        // `resolution`.
        Segment atZero = geometry.memberSegment(0, 0, 0f);
        Segment atTwoPi = geometry.memberSegment(0, 0, TWO_PI);
        assertEquals(atZero.getA().x, atTwoPi.getA().x, 0.0);
        assertEquals(atZero.getA().y, atTwoPi.getA().y, 0.0);
        assertEquals(atZero.getA().z, atTwoPi.getA().z, 0.0);
        assertEquals(atZero.getB().x, atTwoPi.getB().x, 0.0);
        assertEquals(atZero.getB().y, atTwoPi.getB().y, 0.0);
        assertEquals(atZero.getB().z, atTwoPi.getB().z, 0.0);

        // Exact step boundaries: angle == step * angularResolution,
        // computed with the same float/double arithmetic MemberGeometry's
        // stepOf() uses internally (angle is a float parameter; TWO_PI is
        // float; angularResolution is TWO_PI/resolution widened to
        // double). This pins the currently-observed rounding outcome as a
        // regression contract rather than an untested assumption.
        for (int step = 0; step < RESOLUTION; step += 37) {
            float boundaryAngle = (float) (step * angularResolution);
            int expectedStep = expectedStepOf(boundaryAngle,
                                              angularResolution);
            Segment actual = geometry.memberSegment(0, 0, boundaryAngle);
            double expectedLength = table.lengthAt(expectedStep) * scale;
            assertEquals("step " + step + " boundary angle " + boundaryAngle,
                        expectedLength, actual.length(), 1e-6);
        }
    }

    /**
     * Locked-regression duplicate of {@code MemberGeometry.stepOf}'s exact
     * arithmetic (that private method can't be called directly from here).
     * NOT an independent derivation - a latent rounding bug in stepOf would
     * be reproduced here, not caught. Its value is pinning the currently
     * observed boundary behavior against future change, in {@link
     * #stepBoundaryAndWrapAreLockedContract()} - see that test's Javadoc.
     */
    private static int expectedStepOf(float angle, double angularResolution) {
        double normalized = angle % TWO_PI;
        if (normalized < 0) {
            normalized += TWO_PI;
        }
        return ((int) (normalized / angularResolution)) % RESOLUTION;
    }

    @Test
    public void rejectsMalformedResolution() {
        try {
            new MemberGeometry(7, LgaTestGeometry.BASELINE_RADIUS);
            assertTrue("expected IllegalArgumentException for resolution not divisible by 8",
                      false);
        } catch (IllegalArgumentException expected) {
            // expected
        }
        try {
            new MemberGeometry(0, LgaTestGeometry.BASELINE_RADIUS);
            assertTrue("expected IllegalArgumentException for non-positive resolution",
                      false);
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }
}
