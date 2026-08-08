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
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.HashSet;
import java.util.Set;

import javax.vecmath.Point3i;

import org.junit.Test;

/**
 * @author halhildebrand
 *
 */
public class FccNeighborhoodTest {

    @Test
    public void twelveDistinctOffsetsAllParityPreserving() {
        assertEquals(12, FccNeighborhood.DIRECTIONS.size());
        Set<Point3i> offsets = new HashSet<>();
        for (int d : FccNeighborhood.DIRECTIONS) {
            Point3i offset = FccNeighborhood.offsetOf(d);
            int sum = offset.x + offset.y + offset.z;
            assertEquals("offset for direction " + d + " (" + offset
                         + ") must preserve (i+j+k) parity, sum was " + sum,
                         0, Math.floorMod(sum, 2));
            offsets.add(offset);
        }
        assertEquals("all 12 offsets must be distinct", 12, offsets.size());
    }

    @Test
    public void offsetsFormSixOppositePairs() {
        for (int d : FccNeighborhood.DIRECTIONS) {
            int opposite = FccNeighborhood.opposite(d);
            Point3i offset = FccNeighborhood.offsetOf(d);
            Point3i oppositeOffset = FccNeighborhood.offsetOf(opposite);

            assertEquals("direction " + d + " and its opposite " + opposite
                         + " must sum to the zero vector", 0,
                         offset.x + oppositeOffset.x);
            assertEquals("direction " + d + " and its opposite " + opposite
                         + " must sum to the zero vector", 0,
                         offset.y + oppositeOffset.y);
            assertEquals("direction " + d + " and its opposite " + opposite
                         + " must sum to the zero vector", 0,
                         offset.z + oppositeOffset.z);

            assertEquals("opposite(opposite(d)) must equal d", d,
                         FccNeighborhood.opposite(opposite));
        }
    }

    @Test
    public void allTwelveNeighborsAreEquidistantInLatticeMetric() {
        Integer expectedLengthSquared = null;
        for (int d : FccNeighborhood.DIRECTIONS) {
            Point3i offset = FccNeighborhood.offsetOf(d);
            int lengthSquared = offset.x * offset.x + offset.y * offset.y
                                 + offset.z * offset.z;
            if (expectedLengthSquared == null) {
                expectedLengthSquared = lengthSquared;
            } else {
                assertEquals("direction " + d
                             + " is not equidistant with the other FCC neighbors",
                             expectedLengthSquared.intValue(), lengthSquared);
            }
        }
    }

    @Test
    public void periodicWrapPreservesEvenParityOnEvenExtents() {
        Point3i extent = new Point3i(6, 6, 6);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);

        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    if ((i + j + k) % 2 != 0) {
                        continue;
                    }
                    Point3i cell = new Point3i(i, j, k);
                    for (int d : FccNeighborhood.DIRECTIONS) {
                        Point3i neighbor = neighborhood.neighbor(cell, d);
                        int sum = neighbor.x + neighbor.y + neighbor.z;
                        assertEquals("wrapped neighbor of " + cell
                                     + " in direction " + d + " (-> "
                                     + neighbor
                                     + ") must remain even-parity", 0,
                                     Math.floorMod(sum, 2));
                    }
                }
            }
        }
    }

    @Test
    public void oddExtentIsRejected() {
        assertAxisRejected(new Point3i(5, 6, 6), "x");
        assertAxisRejected(new Point3i(6, 5, 6), "y");
        assertAxisRejected(new Point3i(6, 6, 5), "z");
    }

    @Test
    public void zeroExtentIsRejected() {
        assertAxisRejected(new Point3i(0, 6, 6), "x");
        assertAxisRejected(new Point3i(6, 0, 6), "y");
        assertAxisRejected(new Point3i(6, 6, 0), "z");
    }

    /**
     * Extent 2 is even and positive, but rejected anyway: canonical
     * direction pairs that differ by exactly 2 on that axis - {@code
     * (+1,+3)}, {@code (+2,+5)}, {@code (+4,+6)} - alias to the same
     * wrapped neighbor cell at axis extent 2, breaking the "12 distinct
     * neighbors" guarantee (bead inviscid-cb7, caught by stacked review on
     * {@code ContactScan}: empirically confirmed via 32 duplicate {@code
     * Contact} entries at extent {@code (4,2,4)}). See {@code
     * FccNeighborhood}'s class Javadoc "Minimum extent 4 per axis"
     * section.
     */
    @Test
    public void axisExtentTwoIsRejected() {
        assertAxisRejected(new Point3i(2, 6, 6), "x");
        assertAxisRejected(new Point3i(6, 2, 6), "y");
        assertAxisRejected(new Point3i(6, 6, 2), "z");
    }

    @Test
    public void negativeExtentIsRejected() {
        assertAxisRejected(new Point3i(-6, 6, 6), "x");
        assertAxisRejected(new Point3i(6, -6, 6), "y");
        assertAxisRejected(new Point3i(6, 6, -6), "z");
    }

    @Test
    public void directionsIsImmutable() {
        assertThrows(UnsupportedOperationException.class,
                     () -> FccNeighborhood.DIRECTIONS.set(0, 99));
    }

    @Test
    public void neighborRelationIsSymmetric() {
        Point3i extent = new Point3i(6, 6, 6);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);

        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    Point3i cell = new Point3i(i, j, k);
                    for (int d : FccNeighborhood.DIRECTIONS) {
                        Point3i forward = neighborhood.neighbor(cell, d);
                        Point3i back = neighborhood.neighbor(forward,
                                                              FccNeighborhood.opposite(d));
                        assertEquals("neighbor(neighbor(c,d), opposite(d)) must equal c for cell "
                                     + cell + " direction " + d, cell, back);
                    }
                }
            }
        }
    }

    private void assertAxisRejected(Point3i extent, String axis) {
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                                                   () -> new FccNeighborhood(extent));
        assertTrue("exception message must name the offending axis '" + axis
                   + "' but was: " + e.getMessage(),
                   e.getMessage() != null && e.getMessage().contains(axis));
    }
}
