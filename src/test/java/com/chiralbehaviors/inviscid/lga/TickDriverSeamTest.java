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
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.QuantaField;
import com.chiralbehaviors.inviscid.measure.CollisionStatistics;

/**
 * Conformance tests for the {@link TickReport} / {@link TickDriver} tick
 * seam (bead inviscid-ckn / inviscid-0nx.21, T2
 * design-ckn-lattice-seam.md §3). {@link CollisionSweep.TickResult}
 * adopts {@link TickReport} with a zero-body change (it already declares
 * both accessors under those exact names); {@link HybridAutomaton} adopts
 * {@link TickDriver} by adding {@link HybridAutomaton#field()}.
 *
 * @author halhildebrand
 */
public class TickDriverSeamTest {

    private static final double RADIUS     = 0.015;
    private static final int    RESOLUTION = 360;

    @Test
    public void collisionSweepTickResultIsATickReport() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            new ContactPredicate(new MemberGeometry(RESOLUTION,
                                                                                     RADIUS)));
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    new CollisionStatistics());

        CollisionSweep.TickResult result = sweep.tick(0);
        TickReport report = result;
        assertEquals(result.tick(), report.tick());
        assertEquals(result.signedTransferTotal(),
                     report.signedTransferTotal());
    }

    @Test
    public void hybridAutomatonIsATickDriverWhoseFieldIsItsOwnAutomaton() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        FccNeighborhood neighborhood = new FccNeighborhood(extent);
        ContactScan scan = new ContactScan(automaton, neighborhood,
                                            new ContactPredicate(new MemberGeometry(RESOLUTION,
                                                                                     RADIUS)));
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                    new QuantaExchangeRule(),
                                                    new CollisionStatistics());
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);

        TickDriver driver = hybrid;
        QuantaField field = driver.field();
        assertSame("field() must be the exact automaton instance HybridAutomaton drives",
                   automaton, field);

        TickReport report = driver.tick(0);
        assertEquals(0, report.tick());
        assertTrue("a provably zero-sum tick must report signedTransferTotal == 0",
                   report.signedTransferTotal() == 0L);
    }
}
