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

import java.util.ArrayList;
import java.util.List;

import com.chiralbehaviors.inviscid.automaton.Necronomata;

/**
 * Headless-testable fix for inviscid-0nx.24: maps a visualized-cell ordinal
 * (position in the even-parity-only cell list {@code
 * NecronomataVisualization.createCellAnimators} builds via {@link
 * Necronomata#forEach}) to the cell's true flat-array offset ({@link
 * Necronomata#indexOfCell(javax.vecmath.Point3i)}).
 *
 * <p>{@code NecronomataVisualization} previously computed this offset as
 * {@code cell * STRUTS_PER_CELL}, where {@code cell} was simply the
 * ordinal position in that list. {@code indexOfCell} is a row-major index
 * over ALL cells (even and odd parity), so the two mappings coincide only
 * for cell 0 — every subsequent visualized cell read a wrong 30-float
 * slice of the automaton's state array. This class replays the same
 * visitation order and records each cell's real offset, so callers index
 * by {@code cellOffsets(automata)[cell]} instead.
 */
public class VisualizationCellIndex {

    private VisualizationCellIndex() {
    }

    /**
     * @return the flat-array offset (30 floats/cell) for each visualized
     *         cell, in {@link Necronomata#forEach} visitation order —
     *         i.e. {@code offsets[n]} is {@code
     *         automata.indexOfCell(nth even-parity cell visited)}.
     */
    public static int[] cellOffsets(Necronomata automata) {
        List<Integer> offsets = new ArrayList<>();
        automata.forEach(cell -> offsets.add(automata.indexOfCell(cell)));
        int[] result = new int[offsets.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = offsets.get(i);
        }
        return result;
    }
}
