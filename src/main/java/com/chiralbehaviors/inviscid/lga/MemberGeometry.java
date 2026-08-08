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

import javax.vecmath.Matrix3d;
import javax.vecmath.Vector3d;

import com.chiralbehaviors.inviscid.LengthTable;
import com.chiralbehaviors.inviscid.PhiCoordinates;

/**
 * Headless (pure javax.vecmath) reconstruction of the member tip / strut
 * geometry that {@code NecronomataVisualization} composes as a JavaFX
 * {@code Transform} chain:
 *
 * <pre>
 *   position * base(yAxis(Cubes[cube])) * rotations[member][step]
 *            * translate(+/-halfInterval) * scale(1, LengthTable.lengthAt(step), 1)
 * </pre>
 *
 * This class reproduces every factor except {@code position} (world-space
 * cell placement is orthogonal to member geometry and is left to callers),
 * so member tip positions - and therefore member-to-member contact
 * distances - can be computed without a JavaFX toolkit. Two separate
 * guarantees, not one: repeated calls with identical inputs are bitwise
 * identical (no shared mutable state, see {@code
 * MemberGeometryTest.geometryIsPureFunctionOfAngle}), and the geometry
 * itself is verified to agree with the actual JavaFX composition to 1e-5
 * (floating-point tolerance, not bitwise) by {@code
 * MemberGeometryJavaFxParityTest}.
 *
 * <h2>Angular quantization</h2>
 * Member rotation AND length modulation are driven off a discrete
 * {@code step} derived from the continuous {@code angle}, exactly as
 * {@code NecronomataVisualization.setState} does: {@code step =
 * floor(normalize(angle) / (2*PI/resolution)) % resolution}. This is not
 * cosmetic - both the visualization and this class key length and rotation
 * off the same LUT index, so any {@code resolution}-scale discretization is
 * shared between the two.
 *
 * <h2>Length modulation: {@code LengthTable} is authoritative</h2>
 * {@code NecronomataVisualization.buildLengths} overrides four of its
 * rendered per-step lengths - indices 1, 2, 3 (forced to match index 0) and
 * index 5 (forced to match index 4) - regardless of what
 * {@code LengthTable.lengthAt(step)} computes for those steps:
 * {@code lengths[1] = lengths[2] = lengths[3] = lengths[0]; lengths[5] =
 * lengths[4];}. This is NOT a general small-resolution phenomenon - an
 * exhaustive sweep of every resolution divisible by 8 from 8 to 400
 * (verified directly against {@code LengthTable.lengthAt}, not inferred)
 * found exactly three resolutions where any of the four overridden indices
 * happens to already equal its override target:
 * <ul>
 * <li>{@code resolution == 8}: degenerate - {@code LengthTable} holds only
 * one distinct sample ({@code resolution/8 == 1}), so every step trivially
 * equals every other step, including all four override targets.</li>
 * <li>{@code resolution == 16}: step 3 only ({@code lengthAt(3) ==
 * lengthAt(0)} coincidentally, via {@code LengthTable}'s own reflection
 * mapping step 3 to its index-0 sample at that specific resolution).</li>
 * <li>{@code resolution == 40}: step 5 only ({@code lengthAt(5) ==
 * lengthAt(4)}, same kind of isolated coincidence).</li>
 * </ul>
 * Steps 1 and 2 never coincide with step 0 at any of the 49 non-degenerate
 * resolutions checked. This codebase's established precedent resolution is
 * 3600 (see the design doc's "3600-step LUT precedent"), at which none of
 * the four indices coincide - so at the resolution this code is actually
 * exercised with, steps 1, 2, 3 and 5 sit deep inside {@code
 * LengthTable}'s first, unreflected octant sample range (the first fold is
 * at step {@code resolution/8 == 450}) and the override is a genuine,
 * independent, hardcoded flattening of the first few rendered steps -
 * unrelated to the reflection scheme, and not something that happens to
 * "wash out" at realistic resolutions. Per the design doc,
 * member length modulation is "load-bearing for contact detection, not
 * cosmetic", so this class follows {@code LengthTable.lengthAt(step)}
 * directly - the physical truth - and does NOT reproduce the
 * visualization's four-step override. {@code
 * MemberGeometryJavaFxParityTest} therefore excludes steps 1, 2, 3 and 5
 * from its sweep (same reasoning recorded at the exclusion site); the
 * visualization-side discrepancy is reported back as a separate bug bead
 * rather than silently carried into this, the physics-authoritative
 * geometry source.
 *
 * @author halhildebrand
 */
public class MemberGeometry {

    private static final Vector3d CANONICAL_Y_AXIS = new Vector3d(0, 1, 0);
    private static final Vector3d X_AXIS           = new Vector3d(1, 0, 0);
    private static final Vector3d Z_AXIS           = new Vector3d(0, 0, 1);

    /**
     * Rotation matrix about {@code axis} (need not be unit length; this
     * normalizes it, mirroring what {@code javafx.scene.transform.Rotate}
     * does internally before applying the standard axis-angle / Rodrigues
     * formula - required so {@link #computeBase(int)}'s non-unit
     * cross-product axis produces the same rotation JavaFX would build from
     * it) by {@code angleRad} radians.
     */
    private static Matrix3d rotation(double angleRad, Vector3d axis) {
        Matrix3d m = new Matrix3d();
        double mag = axis.length();
        if (mag == 0.0) {
            m.setIdentity();
            return m;
        }
        double x = axis.x / mag;
        double y = axis.y / mag;
        double z = axis.z / mag;
        double sin = Math.sin(angleRad);
        double cos = Math.cos(angleRad);
        double t = 1 - cos;

        m.m00 = t * x * x + cos;
        m.m01 = t * x * y - sin * z;
        m.m02 = t * x * z + sin * y;
        m.m10 = t * x * y + sin * z;
        m.m11 = t * y * y + cos;
        m.m12 = t * y * z - sin * x;
        m.m20 = t * x * z - sin * y;
        m.m21 = t * y * z + sin * x;
        m.m22 = t * z * z + cos;
        return m;
    }

    private static Vector3d transform(Matrix3d m, Vector3d in) {
        Vector3d out = new Vector3d();
        m.transform(in, out);
        return out;
    }

    private final double      angularResolution;
    /**
     * {@code base(cube)} for each of the 5 sub-cubes, computed once at
     * construction. {@code PhiCoordinates.Cubes[cube].getFaces().get(2)
     * .centroid()} triangulates the face on every call, and {@code base()}
     * chains that with an {@code acos}/cross-product/Rodrigues-matrix
     * build - all of it invariant per cube. Caching turns a per-{@code
     * memberSegment()}-call cost into a one-time, 5-entry cost; the
     * contact predicate (inviscid-0nx.12) and phase-quantized LUT
     * (inviscid-0nx.18) both call this in cells x members x timesteps hot
     * loops. Never mutated after construction (rotation() always returns a
     * freshly allocated Matrix3d, and {@link #transform(Matrix3d, Vector3d)}
     * only reads its matrix argument) and never exposed outside this
     * class, so entries are used directly rather than defensively copied
     * per call - copying would reintroduce the allocation this cache
     * exists to avoid.
     */
    private final Matrix3d[]  baseByCube;
    private final double      halfInterval;
    private final double      halfSegmentLength;
    private final double      radius;
    private final int         resolution;
    private final LengthTable table;

    /**
     * @param resolution
     *            angular quantization step count; must be > 0 and divisible
     *            by 8 (mirrors {@code NecronomataVisualization}'s
     *            constructor assertion and {@link LengthTable}'s own
     *            constructor requirement).
     * @param radius
     *            member (strut) radius. This is a model parameter, not a
     *            derived geometric quantity - {@code
     *            NecronomataVisualization} passes {@code 0.015f} at its call
     *            site (an arbitrary visual-thickness choice at cube-edge
     *            scale, with no bearing on centerline / tip positions).
     *            Callers doing contact-distance math should pass whatever
     *            radius the automaton model defines the member's physical
     *            thickness to be; 0.015 is only the value the existing
     *            visualization happens to render with.
     */
    public MemberGeometry(int resolution, double radius) {
        if (resolution <= 0 || resolution % 8 != 0) {
            throw new IllegalArgumentException("resolution must be > 0 and divisible by 8: "
                                                + resolution);
        }
        this.resolution = resolution;
        this.radius = radius;
        this.angularResolution = TWO_PI / resolution;
        this.table = new LengthTable(resolution);
        // NecronomataVisualization always derives the exemplar mesh length
        // and the cell-placement half-interval from Cubes[0]'s edge length,
        // regardless of which of the 5 cube indices a given strut belongs
        // to - Cubes[cube] is used ONLY for orientation, via base(). Match
        // that exactly, or geometry for cube != 0 would silently diverge
        // from the JavaFX path.
        double edgeLength = PhiCoordinates.Cubes[0].getEdgeLength();
        this.halfSegmentLength = edgeLength * ROOT_2 / 2.0;
        this.halfInterval = edgeLength / 2.0;

        this.baseByCube = new Matrix3d[PhiCoordinates.Cubes.length];
        for (int cube = 0; cube < baseByCube.length; cube++) {
            baseByCube[cube] = computeBase(cube);
        }
    }

    /**
     * @return the rotation that reorients the canonical (world Y-axis)
     *         member layout onto the given cube's own Y axis, i.e.
     *         {@code NecronomataVisualization.base(yAxis(Cubes[cube]))}.
     *         Called once per cube at construction; see {@link
     *         #baseByCube}.
     */
    private static Matrix3d computeBase(int cube) {
        Vector3d yAxis = PhiCoordinates.Cubes[cube].getFaces()
                                                    .get(2)
                                                    .centroid();
        Vector3d axisOfRotation = new Vector3d();
        axisOfRotation.cross(yAxis, CANONICAL_Y_AXIS);
        Vector3d normalized = new Vector3d(yAxis);
        normalized.normalize();
        double angle = Math.acos(normalized.dot(CANONICAL_Y_AXIS));
        return rotation(-angle, axisOfRotation);
    }

    /**
     * @return this instance's member radius, as supplied at construction.
     *         See the constructor Javadoc for provenance.
     */
    public double memberRadius() {
        return radius;
    }

    /**
     * @param cube
     *            which of the 5 sub-cubes (0..4) the member belongs to
     * @param member
     *            which of the 6 members of the cube (0..5)
     * @param angle
     *            continuous rotation angle, radians; quantized internally to
     *            the same {@code step} the JavaFX visualization would use
     * @return the two endpoints of the member at the given angle, in
     *         cell-local coordinates (everything in the JavaFX composition
     *         except the world-space {@code position} translate).
     */
    public Segment memberSegment(int cube, int member, float angle) {
        validateCube(cube);
        validateMember(member);

        int step = stepOf(angle);
        double length = table.lengthAt(step);
        double stepAngleRad = step * angularResolution;

        Vector3d offset = translateOffset(member);

        Vector3d rawNeg = new Vector3d(0, -halfSegmentLength * length, 0);
        rawNeg.add(offset);
        Vector3d rawPos = new Vector3d(0, halfSegmentLength * length, 0);
        rawPos.add(offset);

        Matrix3d rotate = rotationForMember(member, stepAngleRad);
        Matrix3d baseRotation = baseByCube[cube];

        Vector3d a = transform(baseRotation, transform(rotate, rawNeg));
        Vector3d b = transform(baseRotation, transform(rotate, rawPos));

        return new Segment(a, b);
    }

    /**
     * @return the member's tip: the endpoint of
     *         {@link #memberSegment(int, int, float)} corresponding to the
     *         pre-transform local point {@code (0, +halfSegmentLength, 0)}
     *         for even members (0, 2, 4), or {@code (0, -halfSegmentLength,
     *         0)} for odd members (1, 3, 5) - i.e. {@code Segment.getB()}
     *         for even, {@code Segment.getA()} for odd.
     *         <p>
     *         This is a convention, not a rotation-invariant geometric
     *         distinction: because every member's translate offset is
     *         perpendicular to its own pre-rotation local Y axis, both
     *         segment endpoints are always exactly equidistant from the
     *         cell-local origin (offset o and half-arm h0y satisfy
     *         o . h0y == 0, so |o + h0y| == |o - h0y|) - "farthest from
     *         origin" cannot distinguish them. The even/odd split above is
     *         the one that makes each mirror pair (0,1), (2,3), (4,5) -
     *         built with opposite rotation sign AND opposite translate
     *         sign - land on true geometric negatives of each other at
     *         angle 0 (see {@code MemberGeometryTest
     *         .oppositeMembersAreAntipodalAtZeroAngle}), which is the
     *         behaviorally meaningful property for a contact predicate.
     */
    public Vector3d memberTip(int cube, int member, float angle) {
        Segment segment = memberSegment(cube, member, angle);
        return (member % 2 == 0) ? segment.getB() : segment.getA();
    }

    private Matrix3d rotationForMember(int member, double stepAngleRad) {
        switch (member) {
        case 0:
            return rotation(stepAngleRad, Z_AXIS);
        case 1:
            return rotation(-stepAngleRad, Z_AXIS);
        case 2:
            return rotation(stepAngleRad, X_AXIS);
        case 3:
            return rotation(-stepAngleRad, X_AXIS);
        case 4:
            return concatenate(rotation(Math.PI / 2.0, Z_AXIS),
                               rotation(stepAngleRad, X_AXIS));
        case 5:
            return concatenate(rotation(Math.PI / 2.0, Z_AXIS),
                               rotation(-stepAngleRad, X_AXIS));
        default:
            throw new IllegalArgumentException("member: " + member);
        }
    }

    /**
     * @return the matrix equivalent of JavaFX's {@code outer.createConcatenation(inner)}
     *         - i.e. apply {@code inner} first, then {@code outer}.
     */
    private Matrix3d concatenate(Matrix3d outer, Matrix3d inner) {
        Matrix3d m = new Matrix3d();
        m.mul(outer, inner);
        return m;
    }

    private int stepOf(float angle) {
        double normalized = angle % TWO_PI;
        if (normalized < 0) {
            normalized += TWO_PI;
        }
        return ((int) (normalized / angularResolution)) % resolution;
    }

    private Vector3d translateOffset(int member) {
        switch (member) {
        case 0:
            return new Vector3d(0, 0, halfInterval);
        case 1:
            return new Vector3d(0, 0, -halfInterval);
        case 2:
        case 4:
            return new Vector3d(halfInterval, 0, 0);
        case 3:
        case 5:
            return new Vector3d(-halfInterval, 0, 0);
        default:
            throw new IllegalArgumentException("member: " + member);
        }
    }

    private void validateCube(int cube) {
        if (cube < 0 || cube >= PhiCoordinates.Cubes.length) {
            throw new IllegalArgumentException("cube: " + cube);
        }
    }

    private void validateMember(int member) {
        if (member < 0 || member > 5) {
            throw new IllegalArgumentException("member: " + member);
        }
    }
}
