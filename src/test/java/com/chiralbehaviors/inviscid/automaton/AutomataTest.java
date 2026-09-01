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

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3i;

import org.junit.Test;

/**
 * @author halhildebrand
 *
 */
public class AutomataTest {
    @Test
    public void testIteration() {
        // (7, 6, 6): odd x extent. This is a non-PBC configuration —
        // Necronomata itself accepts any extent, but
        // com.chiralbehaviors.inviscid.automaton.lga.FccNeighborhood's periodic-wrap
        // neighbor lookup requires all-even extents to keep the
        // even-parity sublattice closed under wrap (inviscid-0nx.3). This
        // test only exercises forEach()/iterator() parity traversal, not
        // wrap-around neighbor lookup, so the odd extent is fine as-is.
        Necronomata automata = new Necronomata(7, 6, 6);
        List<Point3i> loop = new ArrayList<>();

        automata.forEach(c -> loop.add(c));

        List<Point3i> collected = new ArrayList<>();
        for (Point3i c : automata) {
            collected.add(c);
        }

        assertEquals(loop, collected);

    }
}
