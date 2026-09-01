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

package com.chiralbehaviors.inviscid.automaton;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.HashSet;
import java.util.Set;

import javax.vecmath.Point3i;

import org.junit.Test;

/**
 * Regression tests for {@link Necronomata#indexOfCell(int, int, int)}.
 *
 * @author halhildebrand
 */
public class NecronomataIndexTest {

    /**
     * Non-cubic extent: every cell's indexOfCell must be a multiple of 30,
     * distinct per (i,j,k), and in bounds.
     */
    @Test
    public void indexOfCellIsRowMajorWithThirtyFloatStride() {
        int x = 3;
        int y = 4;
        int z = 5;
        Necronomata automata = new Necronomata(x, y, z);

        Set<Integer> seen = new HashSet<>();
        int maxIndex = 30 * (x * y * z - 1);

        for (int i = 0; i < x; i++) {
            for (int j = 0; j < y; j++) {
                for (int k = 0; k < z; k++) {
                    int idx = automata.indexOfCell(i, j, k);

                    assertEquals("index must be a multiple of 30 for cell (" + i + "," + j + "," + k + ")", 0,
                                 idx % 30);
                    assertTrue("index must be >= 0 for cell (" + i + "," + j + "," + k + ")", idx >= 0);
                    assertTrue("index must be <= " + maxIndex + " for cell (" + i + "," + j + "," + k + ")",
                               idx <= maxIndex);
                    assertTrue("index must be distinct per cell, duplicate at (" + i + "," + j + "," + k + ")",
                               seen.add(idx));

                    int expected = 30 * ((i * y + j) * z + k);
                    assertEquals("index must match row-major formula for cell (" + i + "," + j + "," + k + ")",
                                 expected, idx);
                }
            }
        }
    }

    /**
     * Round-trip: write cell-identifying values into the angles array at the
     * expected row-major offsets, then verify anglesOf(cell) returns exactly
     * that 30-float slice.
     */
    @Test
    public void anglesOfReturnsExactCellSliceAtRowMajorOffset() {
        int x = 3;
        int y = 4;
        int z = 5;

        float[] angles = new float[30 * x * y * z];
        float[] frequency = new float[30 * x * y * z];

        for (int i = 0; i < x; i++) {
            for (int j = 0; j < y; j++) {
                for (int k = 0; k < z; k++) {
                    int cellIndex = (i * y + j) * z + k;
                    for (int s = 0; s < 30; s++) {
                        angles[30 * cellIndex + s] = cellIndex * 100 + s;
                    }
                }
            }
        }

        Necronomata automata = new Necronomata(angles, new Point3i(x, y, z), frequency);

        Point3i[] cellsToCheck = { new Point3i(0, 0, 0), new Point3i(x - 1, y - 1, z - 1), new Point3i(1, 2, 3),
                                    new Point3i(2, 0, 4), new Point3i(0, 3, 0) };

        for (Point3i cell : cellsToCheck) {
            int cellIndex = (cell.x * y + cell.y) * z + cell.z;
            float[] expected = new float[30];
            for (int s = 0; s < 30; s++) {
                expected[s] = cellIndex * 100 + s;
            }

            float[] actual = automata.anglesOf(cell);
            assertArrayEquals("anglesOf must return exact 30-float slice for cell " + cell, expected, actual, 0.0f);
        }
    }
}
