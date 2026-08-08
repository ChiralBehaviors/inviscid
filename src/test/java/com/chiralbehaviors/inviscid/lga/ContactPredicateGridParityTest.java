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

import static org.junit.Assert.assertEquals;

import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import org.junit.Test;

import com.chiralbehaviors.inviscid.CubicGrid;
import com.chiralbehaviors.inviscid.CubicGrid.Neighborhood;
import com.chiralbehaviors.inviscid.PhiCoordinates;

import javafx.geometry.Point3D;
import javafx.scene.transform.Transform;

/**
 * Regression test for {@link ContactPredicate}'s physical cross-cell offset
 * derivation (inviscid-egm, raised in substantive-critique of
 * inviscid-0nx.12). {@link ContactPredicate#physicalOffset(int)} recomputes
 * {@code GRID_AXIS_X/Y/Z} and {@code CELL_SPACING} independently from
 * {@code PhiCoordinates.Cubes[3]} rather than delegating to {@code
 * CubicGrid} - see that class's Javadoc for the full derivation. That
 * independence is exactly the risk: an axis swap, a sign flip, or a
 * {@code PhiCoordinates}/{@code CubicGrid} refactor could silently desync
 * the two, producing a geometrically-plausible but WRONG-neighbor
 * predicate that every other {@code ContactPredicate} test (fixture-driven
 * from the predicate's own output, or invariant under coherent axis
 * relabeling - see {@code ContactPredicateTest}'s class Javadoc) would
 * still pass.
 * <p>
 * This test builds the REAL {@code CubicGrid} the same way {@code
 * NecronomataVisualization} does -
 * {@code new CubicGrid(Neighborhood.SIX, PhiCoordinates.Cubes[3], extent)}
 * (see that class's constructor call site) - and cross-checks, for all 12
 * {@code FccNeighborhood} directions, that {@code
 * grid.positionTransform(offset.x, offset.y, offset.z)} (applied to the
 * origin) equals {@code ContactPredicate.physicalOffset(direction)}. The
 * {@code extent} argument does not affect {@code positionTransform} (only
 * {@code xAxis/yAxis/zAxis} and {@code intervalX/Y/Z} do, which are
 * extent-independent), so any positive value mirrors the visualization's
 * call shape faithfully; 4 is used as the smallest extent {@code
 * FccNeighborhood} itself accepts.
 *
 * <h2>Tolerance</h2>
 * {@code ContactPredicate}'s class Javadoc documents a ~7e-9 relative
 * discrepancy between {@code Cubes[0]} (its {@code CELL_SPACING} source)
 * and {@code Cubes[3]} (its axis source) edge lengths - at this class's
 * {@code CELL_SPACING} magnitude (~5.236, per {@code
 * ContactPredicateTest}'s class Javadoc), that is on the order of 4e-8
 * absolute. {@link #DELTA} (1e-6) comfortably clears that known,
 * benign source of noise while remaining many orders of magnitude
 * tighter than the O(1) (~5.236 magnitude) displacement an axis
 * permutation or sign error would introduce - perturbation evidence
 * recorded as a comment on bead inviscid-egm (an axis swap in
 * physicalOffset was observed to fail this test immediately).
 *
 * @author halhildebrand
 */
public class ContactPredicateGridParityTest {

    private static final double DELTA  = 1e-6;
    private static final int    EXTENT = 4;

    @Test
    public void physicalOffsetMatchesCubicGridPositionTransformForAllTwelveDirections() {
        CubicGrid grid = new CubicGrid(Neighborhood.SIX,
                                       PhiCoordinates.Cubes[3], EXTENT);
        ContactPredicate predicate = new ContactPredicate(new MemberGeometry(360,
                                                                              0.015));

        for (int direction : FccNeighborhood.DIRECTIONS) {
            Point3i offset = FccNeighborhood.offsetOf(direction);

            Transform positionTransform = grid.positionTransform(offset.x,
                                                                  offset.y,
                                                                  offset.z);
            Point3D gridOffset = positionTransform.transform(Point3D.ZERO);

            Vector3d predicateOffset = predicate.physicalOffset(direction);

            assertEquals("direction " + direction + " x: CubicGrid "
                        + gridOffset + " vs ContactPredicate "
                        + predicateOffset, gridOffset.getX(),
                        predicateOffset.x, DELTA);
            assertEquals("direction " + direction + " y: CubicGrid "
                        + gridOffset + " vs ContactPredicate "
                        + predicateOffset, gridOffset.getY(),
                        predicateOffset.y, DELTA);
            assertEquals("direction " + direction + " z: CubicGrid "
                        + gridOffset + " vs ContactPredicate "
                        + predicateOffset, gridOffset.getZ(),
                        predicateOffset.z, DELTA);
        }
    }
}
