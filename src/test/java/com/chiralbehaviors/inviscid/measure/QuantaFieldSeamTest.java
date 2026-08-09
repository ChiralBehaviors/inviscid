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

package com.chiralbehaviors.inviscid.measure;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.HashSet;
import java.util.Set;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.QuantaField;

/**
 * Conformance tests for the {@link QuantaField} read seam (bead
 * inviscid-ckn / inviscid-0nx.21, T2 {@code
 * inviscid/design-ckn-lattice-seam.md} §6.2). Tests numbered per the
 * design memo's test list.
 *
 * <p>SCOPE (checkpoint after seam checklist steps 0-7): only
 * {@link Necronomata}'s {@link QuantaField} implementation is exercised
 * here. Tests 4a ({@code phaseAtIsTheFineAccumulatorNotTheContactBinCentre})
 * and 4b ({@code bothSubstratesReportPhaseResolution3600UnderCadence2A})
 * are LGA-side pins and are deferred to when {@code LatticeGasAutomaton}
 * (checklist step 10) lands — they cannot be written against a substrate
 * that does not exist yet. This is a deliberate, reported scope
 * narrowing, not a silent omission.
 *
 * @author halhildebrand
 */
public class QuantaFieldSeamTest {

    private static final int MEMBERS_PER_CELL = 30;

    private Necronomata freshAutomaton(Point3i extent) {
        return new Necronomata(extent);
    }

    /** Test 1. */
    @Test
    public void slotCountIsThirtyPerCellForBothSubstrates() {
        Point3i extent = new Point3i(4, 4, 4);
        QuantaField field = freshAutomaton(extent);
        assertEquals(MEMBERS_PER_CELL * extent.x * extent.y * extent.z,
                     field.slotCount());
    }

    /** Test 2. */
    @Test
    public void indexOfCellRoundTripsForEveryEvenParityCell() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = freshAutomaton(extent);
        QuantaField field = automaton;
        for (Point3i cell : automaton) {
            assertEquals("indexOfCell must agree between Necronomata and the seam for "
                         + cell, automaton.indexOfCell(cell),
                         field.indexOfCell(cell));
        }
    }

    /** Test 3. */
    @Test
    public void forEachCellVisitsExactlyTheEvenParitySublatticeAndNothingElse() {
        Point3i extent = new Point3i(4, 4, 4);
        QuantaField field = freshAutomaton(extent);
        Set<Point3i> visited = new HashSet<>();
        field.forEachCell(visited::add);

        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    Point3i cell = new Point3i(i, j, k);
                    boolean evenParity = (i + j + k) % 2 == 0;
                    assertEquals("cell " + cell
                                 + " inclusion must match even-parity predicate",
                                 evenParity, visited.contains(cell));
                }
            }
        }
    }

    /** Test 4. */
    @Test
    public void quantaAtAgreesWithTheRawFrequencyArrayForEverySlot() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = freshAutomaton(extent);
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < frequency.length; i++) {
                frequency[i] = (i % 5) - 2;
            }
        });
        float[][] oracle = new float[1][];
        automaton.process((angle, frequency, deltaA,
                            deltaF) -> oracle[0] = frequency);

        QuantaField field = automaton;
        for (int i = 0; i < field.slotCount(); i++) {
            assertEquals("slot " + i, Math.round((double) oracle[0][i]),
                         field.quantaAt(i));
        }
    }

    /** Test 5 — negative control. */
    @Test
    public void isExactAtFlagsADeliberatelyCorruptedSlot() {
        Point3i extent = new Point3i(2, 2, 2);
        Necronomata automaton = freshAutomaton(extent);
        automaton.process((angle, frequency, deltaA, deltaF) -> frequency[3] = 1.5f);

        QuantaField field = automaton;
        assertFalse("slot 3 was deliberately corrupted to 1.5",
                    field.isExactAt(3));
        assertTrue("slot 0 was never touched and must still read exact",
                   field.isExactAt(0));
    }
}
