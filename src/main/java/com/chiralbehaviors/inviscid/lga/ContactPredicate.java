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

import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import com.chiralbehaviors.inviscid.PhiCoordinates;

/**
 * The geometric contact predicate (inviscid-0nx.12): do member
 * {@code (cubeA, memberA)} of cell {@code C} and member
 * {@code (cubeB, memberB)} of the neighboring cell {@code C + offset
 * (direction)} touch, given their angles? Members are capsules (a
 * {@link Segment} plus {@link MemberGeometry#memberRadius()}); contact
 * holds iff the segment-segment distance is less than
 * {@code 2 * memberRadius}.
 * <p>
 * This class takes explicit {@code (cube, member)} pairs and a
 * {@code direction} rather than cell coordinates - per {@code
 * FccNeighborhood}'s own Javadoc, the correspondence between direction
 * index and member index is UNVERIFIED, so nothing here assumes member
 * {@code m} faces direction {@code m}; callers supply both.
 *
 * <h2>Inter-cell spacing derivation</h2>
 * {@link MemberGeometry#memberSegment(int, int, float)} returns member
 * endpoints in CELL-LOCAL coordinates - the entire JavaFX transform chain
 * except the world-space {@code position} translate (see that class's
 * Javadoc). To compare a member of cell {@code C} against a member of a
 * neighboring cell, the neighbor's cell-local segment must be translated
 * by the physical world-space offset between the two cells' {@code
 * position} translates.
 * <p>
 * That offset is derived from how {@code NecronomataVisualization}
 * actually places cells, not assumed:
 * <pre>
 * CubicGrid grid = new CubicGrid(Neighborhood.SIX, PhiCoordinates.Cubes[3],
 *                                automata.getExtent().x);
 * ...
 * Transform position = grid.positionTransform(location.x - ..., location.y - ..., location.z - ...);
 * </pre>
 * {@code CubicGrid.positionTransform(i,j,k)} computes {@code i *
 * intervalX * xAxis + j * intervalY * yAxis + k * intervalZ * zAxis},
 * where - for the {@code CubicGrid(Neighborhood, Cube, int)} constructor
 * used above - {@code intervalX == intervalY == intervalZ ==
 * cube.getEdgeLength()} and {@code xAxis / yAxis / zAxis} are the
 * (normalized) face centroids of {@code Cubes[3]} at face indices 1, 2,
 * 0 respectively ({@code CubicGrid.xAxis/yAxis/zAxis} static helpers).
 * {@code (i,j,k)} here are exactly the same cell lattice coordinates
 * {@code FccNeighborhood}'s {@code Point3i} offsets are expressed in
 * (both are transcribed from / consumed by {@code Necronomata}'s
 * {@code Point3i} cell addressing - see {@code FccNeighborhood}'s class
 * Javadoc). So the physical world-space displacement between cell
 * {@code C} and {@code C + offsetOf(direction)} is:
 * <pre>
 * offset.x * spacing * xAxis + offset.y * spacing * yAxis + offset.z * spacing * zAxis
 * </pre>
 * where {@code offset = FccNeighborhood.offsetOf(direction)} and
 * {@code spacing} is the cell edge length.
 * <p>
 * This class recomputes {@code xAxis/yAxis/zAxis} from {@code
 * PhiCoordinates.Cubes[3]} directly (mirroring {@code
 * CubicGrid.xAxis/yAxis/zAxis}) rather than hardcoding the numeric
 * result, so it stays correct if the underlying coordinate table ever
 * changes. Numerically (verified by direct computation against the
 * actual {@code PhiCoordinates}/{@code CubicGrid} classes, not asserted
 * blind) those three axes come out to be the exact signed permutation
 * {@code (0,-1,0)}, {@code (0,0,-1)}, {@code (1,0,0)} of the canonical
 * basis - i.e. genuinely orthonormal and axis-aligned in the same
 * ambient {@code Vector3d} frame {@code MemberGeometry}'s segments live
 * in, not merely orthogonal. There was no ambiguity to resolve: all 5
 * {@code PhiCoordinates.Cubes[i]} share the same edge length to within
 * floating-point precision (~1e-7 relative, from {@code float}-precision
 * {@code PHI} arithmetic upstream), so {@code spacing} is taken from
 * {@code Cubes[0].getEdgeLength()} - the exact same source {@code
 * MemberGeometry} itself uses for its cell-local coordinate scale
 * (halfInterval / halfSegmentLength) - rather than {@code Cubes[3]},
 * keeping every length in this class traceable to one source instead of
 * two numerically-almost-but-not-exactly-equal ones.
 *
 * @author halhildebrand
 */
public class ContactPredicate {

    /**
     * {@code CubicGrid.xAxis(Cubes[3])}: face index 1's centroid,
     * normalized. See the class Javadoc's spacing derivation.
     */
    private static final Vector3d GRID_AXIS_X = unit(PhiCoordinates.Cubes[3].getFaces()
                                                                             .get(1)
                                                                             .centroid());

    /**
     * {@code CubicGrid.yAxis(Cubes[3])}: face index 2's centroid,
     * normalized.
     */
    private static final Vector3d GRID_AXIS_Y = unit(PhiCoordinates.Cubes[3].getFaces()
                                                                             .get(2)
                                                                             .centroid());

    /**
     * {@code CubicGrid.zAxis(Cubes[3])}: face index 0's centroid,
     * normalized.
     */
    private static final Vector3d GRID_AXIS_Z = unit(PhiCoordinates.Cubes[3].getFaces()
                                                                             .get(0)
                                                                             .centroid());

    /**
     * The cell-to-cell lattice spacing, in the same length units {@code
     * MemberGeometry} uses for member geometry. Sourced from {@code
     * Cubes[0]}, matching {@code MemberGeometry}'s own source - see the
     * class Javadoc.
     */
    private static final double CELL_SPACING = PhiCoordinates.Cubes[0].getEdgeLength();

    private static Vector3d unit(Vector3d v) {
        Vector3d copy = new Vector3d(v);
        copy.normalize();
        return copy;
    }

    private final MemberGeometry geometry;

    public ContactPredicate(MemberGeometry geometry) {
        this.geometry = geometry;
    }

    /**
     * @return {@code true} iff member {@code (cubeA, memberA)} of the
     *         local cell and member {@code (cubeB, memberB)} of the
     *         neighboring cell in {@code direction} are within {@code
     *         2 * memberRadius} of each other - i.e. their capsules
     *         touch.
     */
    public boolean contacts(int cubeA, int memberA, float angleA, int cubeB,
                            int memberB, float angleB, int direction) {
        double minDistance = minDistance(cubeA, memberA, angleA, cubeB,
                                         memberB, angleB, direction);
        return minDistance < 2 * geometry.memberRadius();
    }

    /**
     * @return the minimum Euclidean separation between member
     *         {@code (cubeA, memberA)} of the local cell and member
     *         {@code (cubeB, memberB)} of the neighboring cell in
     *         {@code direction} - exposed (rather than only the boolean
     *         {@link #contacts}) so inviscid-0nx.13/.15 can log contact
     *         margins, not just fire/no-fire.
     */
    public double minDistance(int cubeA, int memberA, float angleA,
                              int cubeB, int memberB, float angleB,
                              int direction) {
        Segment segmentA = geometry.memberSegment(cubeA, memberA, angleA);
        Segment localB = geometry.memberSegment(cubeB, memberB, angleB);
        Vector3d offset = physicalOffset(direction);

        Vector3d bA = new Vector3d(localB.getA());
        bA.add(offset);
        Vector3d bB = new Vector3d(localB.getB());
        bB.add(offset);

        return SegmentDistance.distance(segmentA.getA(), segmentA.getB(),
                                        bA, bB);
    }

    /**
     * @return the world-space displacement between a cell and its
     *         neighbor in {@code direction} - see the class Javadoc's
     *         spacing derivation. Package-private (rather than private)
     *         so {@code ContactPredicateGridParityTest} can cross-check
     *         it directly against {@code CubicGrid.positionTransform}
     *         (inviscid-egm) - a natural inspection point for a future
     *         A.5 atlas as well.
     */
    Vector3d physicalOffset(int direction) {
        Point3i offset = FccNeighborhood.offsetOf(direction);
        Vector3d displacement = new Vector3d();
        displacement.scaleAdd(offset.x, GRID_AXIS_X, displacement);
        displacement.scaleAdd(offset.y, GRID_AXIS_Y, displacement);
        displacement.scaleAdd(offset.z, GRID_AXIS_Z, displacement);
        displacement.scale(CELL_SPACING);
        return displacement;
    }
}
