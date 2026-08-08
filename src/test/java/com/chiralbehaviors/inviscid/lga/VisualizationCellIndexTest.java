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
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;

/**
 * Headless-testable pin for inviscid-0nx.24: {@code
 * NecronomataVisualization.setState} computed its flat-array offset as
 * {@code cell * STRUTS_PER_CELL} where {@code cell} is the ordinal position
 * in the even-parity-only cell list built by {@code Necronomata.forEach}.
 * {@code Necronomata.indexOfCell} is a full row-major index over ALL cells
 * (even and odd parity), so the two mappings agree only at cell 0 — every
 * later visualized cell read a wrong 30-float slice.
 *
 * <p>{@link VisualizationCellIndex#cellOffsets(Necronomata)} is the
 * extracted, headless-testable fix: it replays the exact even-parity
 * visitation order {@code NecronomataVisualization.createCellAnimators}
 * uses (via {@code Necronomata.forEach}) and records each visited cell's
 * true {@code indexOfCell} offset, so {@code setState} can index by
 * {@code cellOffsets[cell]} instead of the broken {@code cell *
 * STRUTS_PER_CELL} arithmetic.
 */
public class VisualizationCellIndexTest {

    private static final int STRUTS_PER_CELL = 30;

    @Test
    public void visualizationCellOrderMapsToIndexOfCell() {
        Necronomata automata = new Necronomata(new Point3i(4, 6, 8));

        List<Point3i> cellList = new ArrayList<>();
        automata.forEach(cellList::add);

        int[] offsets = VisualizationCellIndex.cellOffsets(automata);

        assertEquals(cellList.size(), offsets.length);
        for (int n = 0; n < cellList.size(); n++) {
            assertEquals("cell " + n, automata.indexOfCell(cellList.get(n)),
                         offsets[n]);
        }
    }

    @Test
    public void everyVisualizedCellReadsADistinctSlice() {
        Necronomata automata = new Necronomata(new Point3i(4, 6, 8));

        int[] offsets = VisualizationCellIndex.cellOffsets(automata);

        Set<Integer> seen = new HashSet<>();
        for (int offset : offsets) {
            assertTrue("duplicate offset: " + offset, seen.add(offset));
        }
        assertEquals(automata.cellCount(), seen.size());
    }

    @Test
    public void sliceOffsetsAreInBounds() {
        Point3i extent = new Point3i(4, 6, 8);
        Necronomata automata = new Necronomata(extent);
        int arrayLength = STRUTS_PER_CELL * extent.x * extent.y * extent.z;

        int[] offsets = VisualizationCellIndex.cellOffsets(automata);

        assertTrue(offsets.length > 0);
        for (int offset : offsets) {
            assertTrue("offset " + offset + " out of bounds",
                      offset + STRUTS_PER_CELL <= arrayLength);
        }
    }
}
