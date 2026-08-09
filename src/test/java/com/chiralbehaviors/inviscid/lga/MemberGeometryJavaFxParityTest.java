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

import static com.chiralbehaviors.inviscid.Constants.ROOT_2;
import static com.chiralbehaviors.inviscid.Constants.TWO_PI;
import static com.chiralbehaviors.inviscid.CubicGrid.yAxis;
import static org.junit.Assert.assertEquals;

import javax.vecmath.Vector3d;

import org.junit.Test;

import com.chiralbehaviors.inviscid.LengthTable;
import com.chiralbehaviors.inviscid.PhiCoordinates;

import javafx.geometry.Point3D;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Scale;
import javafx.scene.transform.Transform;
import javafx.scene.transform.Translate;

/**
 * The anchor test for {@link MemberGeometry}: verifies its pure-vecmath tip
 * / segment geometry against an independently-built JAVAFX {@code
 * Transform} composition (real {@code javafx.scene.transform.Rotate} /
 * {@code Translate} / {@code Scale} / {@code Transform.createConcatenation}
 * objects doing the actual matrix math), following the exact formula
 * {@code NecronomataVisualization} documents and builds from:
 *
 * <pre>
 *   base(yAxis(Cubes[cube])) * rotations[member][step]
 *          * translate(+/-halfInterval) * scale(1, length, 1)
 * </pre>
 *
 * (world-space {@code position} is omitted on both sides - see {@link
 * MemberGeometry}'s class Javadoc). {@code javafx.scene.transform.*}
 * classes are plain value/matrix types - they do not require
 * {@code Platform.startup} / toolkit initialization to construct or apply,
 * so no fallback was needed here; this runs directly under surefire.
 *
 * <h2>The buildLengths divergence — resolved (inviscid-6cf)</h2>
 * {@code NecronomataVisualization.buildLengths} used to override its
 * rendered per-step length at indices 1, 2, 3 (forced to index 0's value)
 * and index 5 (forced to index 4's value), regardless of what {@code
 * LengthTable.lengthAt(step)} computed for those steps — a genuine,
 * resolution-independent discrepancy from the physical truth (nx_plan_audit
 * F2 ruling; inviscid-t33). That override has been removed: {@code
 * buildLengths} now renders {@code LengthTable.lengthAt(step)} directly for
 * every step, matching {@link MemberGeometry}'s physics-authoritative
 * geometry. This test therefore asserts plain parity ({@link #assertClose})
 * at every step, with no exceptions.
 */
public class MemberGeometryJavaFxParityTest {

    private static final Point3D CANONICAL_Y_AXIS = new Point3D(0, 1, 0);
    private static final double  DELTA            = 1e-5;

    private static final int     RESOLUTION       = 360;

    /**
     * Verbatim port of {@code NecronomataVisualization.base(Point3D)}
     * (private there) - real JavaFX {@code Rotate} built from a real
     * {@code Point3D} cross product / dot product, so this exercises
     * JavaFX's own axis-angle math, not a hand-rolled substitute.
     */
    private static Rotate base(Point3D yAxis) {
        Point3D axisOfRotation = yAxis.crossProduct(CANONICAL_Y_AXIS);
        double angle = Math.acos(yAxis.normalize()
                                      .dotProduct(CANONICAL_Y_AXIS));
        return new Rotate(-Math.toDegrees(angle), axisOfRotation);
    }

    /**
     * Verbatim port of the per-member rotation construction in {@code
     * NecronomataVisualization.buildRotations}, evaluated ad-hoc for a
     * given step instead of precomputed into an array.
     */
    private static Transform rotationFor(int member, int step,
                                         double angularResolutionRad) {
        double degrees = Math.toDegrees(step * angularResolutionRad);
        switch (member) {
        case 0:
            return new Rotate(degrees, new Point3D(0, 0, 1));
        case 1:
            return new Rotate(-degrees, new Point3D(0, 0, 1));
        case 2:
            return new Rotate(degrees, new Point3D(1, 0, 0));
        case 3:
            return new Rotate(-degrees, new Point3D(1, 0, 0));
        case 4:
            return new Rotate(90,
                              new Point3D(0, 0,
                                         1)).createConcatenation(new Rotate(degrees,
                                                                            new Point3D(1,
                                                                                       0,
                                                                                       0)));
        case 5:
            return new Rotate(90,
                              new Point3D(0, 0,
                                         1)).createConcatenation(new Rotate(-degrees,
                                                                            new Point3D(1,
                                                                                       0,
                                                                                       0)));
        default:
            throw new IllegalArgumentException("member: " + member);
        }
    }

    /**
     * Verbatim port of the per-member translate assignment in {@code
     * NecronomataVisualization.createCellAnimators} (the xPos / xNeg /
     * yPos / yNeg locals there).
     */
    private static Translate translateFor(int member, double halfInterval) {
        switch (member) {
        case 0:
            return new Translate(0, 0, halfInterval);
        case 1:
            return new Translate(0, 0, -halfInterval);
        case 2:
        case 4:
            return new Translate(halfInterval, 0, 0);
        case 3:
        case 5:
            return new Translate(-halfInterval, 0, 0);
        default:
            throw new IllegalArgumentException("member: " + member);
        }
    }

    private static Vector3d toVector(Point3D p) {
        return new Vector3d(p.getX(), p.getY(), p.getZ());
    }

    @Test
    public void matchesJavaFxTransformComposition() {
        double angularResolutionRad = TWO_PI / RESOLUTION;
        double edgeLength = PhiCoordinates.Cubes[0].getEdgeLength();
        double halfInterval = edgeLength / 2.0;
        double halfSegmentLength = edgeLength * ROOT_2 / 2.0;
        LengthTable table = new LengthTable(RESOLUTION);

        MemberGeometry geometry = new MemberGeometry(RESOLUTION, LgaTestGeometry.BASELINE_RADIUS);

        int compared = 0;
        for (int cube = 0; cube < 5; cube++) {
            Rotate baseRotation = base(yAxis(PhiCoordinates.Cubes[cube]));
            for (int member = 0; member < 6; member++) {
                Translate translate = translateFor(member, halfInterval);
                for (int step = 0; step < RESOLUTION; step++) {
                    float angle = (float) ((step + 0.5) * angularResolutionRad);
                    double length = table.lengthAt(step);
                    Transform rotate = rotationFor(member, step,
                                                   angularResolutionRad);
                    Scale scale = new Scale(1.0, length, 1.0);

                    Transform combined = baseRotation.createConcatenation(rotate)
                                                     .createConcatenation(translate)
                                                     .createConcatenation(scale);

                    Point3D expectedA = combined.transform(new Point3D(0,
                                                                       -halfSegmentLength,
                                                                       0));
                    Point3D expectedB = combined.transform(new Point3D(0,
                                                                       halfSegmentLength,
                                                                       0));

                    Segment actual = geometry.memberSegment(cube, member,
                                                            angle);

                    assertClose("cube " + cube + " member " + member
                               + " step " + step + " endpoint a",
                               toVector(expectedA), actual.getA());
                    assertClose("cube " + cube + " member " + member
                               + " step " + step + " endpoint b",
                               toVector(expectedB), actual.getB());
                    compared++;
                }
            }
        }

        // 5 cubes * 6 members * 360 steps = 10800: plain parity, no steps
        // excluded (inviscid-6cf removed the buildLengths override that
        // used to force a divergence check at steps 1, 2, 3, 5) -
        // comfortably above the bead's ">= 64 angles" floor.
        assertEquals(5 * 6 * RESOLUTION, compared);
    }

    private void assertClose(String message, Vector3d expected,
                             Vector3d actual) {
        assertEquals(message + " x", expected.x, actual.x, DELTA);
        assertEquals(message + " y", expected.y, actual.y, DELTA);
        assertEquals(message + " z", expected.z, actual.z, DELTA);
    }
}
