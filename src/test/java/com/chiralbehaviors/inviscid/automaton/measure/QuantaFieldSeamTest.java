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

package com.chiralbehaviors.inviscid.automaton.measure;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Set;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.automaton.Necronomata;
import com.chiralbehaviors.inviscid.automaton.QuantaField;
import com.chiralbehaviors.inviscid.automaton.lga.CollisionTable;
import com.chiralbehaviors.inviscid.automaton.lga.ContactAtlas;
import com.chiralbehaviors.inviscid.automaton.lga.LatticeGasAutomaton;
import com.chiralbehaviors.inviscid.automaton.lga.QuantaExchangeRule;

/**
 * Conformance tests for the {@link QuantaField} read seam (bead
 * inviscid-ckn / inviscid-0nx.21, T2 {@code
 * inviscid/design-ckn-lattice-seam.md} §6.2). Tests numbered per the
 * design memo's test list.
 *
 * <p>Tests 1-5 exercise {@link Necronomata}'s implementation; tests 4a/4b
 * exercise {@link LatticeGasAutomaton}'s, now that it exists (checklist
 * step 10 landed).
 *
 * @author halhildebrand
 */
public class QuantaFieldSeamTest {

    private static final String RESOURCE_PATH = "lga/contact-atlas-v2.tsv";

    private static final int MEMBERS_PER_CELL = 30;

    private Necronomata freshAutomaton(Point3i extent) {
        return new Necronomata(extent);
    }

    private static LatticeGasAutomaton freshLga() throws IOException {
        URL resource = QuantaFieldSeamTest.class.getClassLoader()
                                                  .getResource(RESOURCE_PATH);
        Path path = Paths.get(resource.getPath());
        ContactAtlas atlas = ContactAtlas.read(path);
        CollisionTable collisions = CollisionTable.buildFromPhaseARule(new QuantaExchangeRule());
        return new LatticeGasAutomaton(atlas.header().extent(), atlas,
                                        collisions, new CollisionStatistics());
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

    /**
     * Test 3-lga (final-review Significant G fix, T2
     * critique-final-0nx21-automaton-arc.md [21949]): {@link
     * LatticeGasAutomaton#forEachCell} is an INDEPENDENT re-implementation
     * of the {@code (i+j+k)%2==0} even-parity predicate, not delegated to
     * {@code Necronomata} or shared code - test 3 above exercised only
     * {@code Necronomata}'s copy, leaving the LGA's own copy untested
     * despite this exact file being touched to add LGA-specific tests
     * 4a/4b. A transcription slip in either copy (off-by-one, wrong
     * modulus, swapped axis) would go undetected without this twin.
     */
    @Test
    public void forEachCellVisitsExactlyTheEvenParitySublatticeAndNothingElseForLga() throws IOException {
        LatticeGasAutomaton lga = freshLga();
        Point3i extent = lga.extent();
        Set<Point3i> visited = new HashSet<>();
        lga.forEachCell(visited::add);

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

    /**
     * Test 4a (rev 3, T2 design-ckn-lattice-seam.md §6.2): under the
     * committed atlas's M=150 sub-bin steps, driving the accumulator
     * through several sub-bin steps WITHIN one contact bin (bin =
     * phase/150, so phase in {0, 10, ..., 140} all stay in bin 0) must
     * produce a CHANGING {@link QuantaField#phaseAt(int)} reading. Under
     * the REJECTED spec -- reporting the coarse contact-bin CENTRE angle,
     * {@code (bin + 0.5) * 2*pi/N_lga}, instead of the fine accumulator --
     * this would be constant across all 150 sub-steps -- see T2
     * analysis-73v-spectral-conversion-and-cadence.md §5 for the
     * quantisation-sideband artifact that regression would (re)introduce.
     */
    @Test
    public void phaseAtIsTheFineAccumulatorNotTheContactBinCentre() throws IOException {
        LatticeGasAutomaton lga = freshLga();
        int slot = 3;
        int subBinSteps = lga.subBinSteps();
        assertTrue("fixture assumption: need at least a few sub-bin steps to drive through",
                   subBinSteps >= 10);

        float previous = Float.NaN;
        int previousBin = -1;
        for (int phase = 0; phase < subBinSteps; phase += 10) {
            int fixedPhase = phase;
            lga.process((p, quanta) -> p[slot] = fixedPhase);
            int bin = phase / subBinSteps;
            float current = lga.phaseAt(slot);

            if (previousBin == bin && !Float.isNaN(previous)) {
                assertNotEquals("phaseAt must change between sub-bin steps "
                                 + (phase - 10) + " and " + phase
                                 + " within the SAME bin " + bin
                                 + " -- a constant reading here means "
                                 + "phaseAt regressed to the REJECTED spec: "
                                 + "reporting the coarse contact-bin CENTRE "
                                 + "angle ((bin + 0.5) * 2*pi/N_lga), which is "
                                 + "invariant within a bin, instead of the fine "
                                 + "accumulator phase",
                                 previous, current, 0.0f);
            }
            previous = current;
            previousBin = bin;
        }
    }

    /**
     * Test 4b (rev 3): under the user's 2A cadence decision, BOTH
     * substrates report {@code phaseResolution() == 3600} -- {@code
     * N_lga} (24, for the committed atlas) keys contact-table lookups
     * ONLY and is never returned here. This is what keeps a per-substrate
     * axis-rescaling need out of any angle-spectrum instrument (T2
     * analysis-73v-spectral-conversion-and-cadence.md §2.1).
     */
    @Test
    public void bothSubstratesReportPhaseResolution3600UnderCadence2A() throws IOException {
        QuantaField necronomata = freshAutomaton(new Point3i(4, 4, 4));
        QuantaField lga = freshLga();

        assertEquals(3600, necronomata.phaseResolution());
        assertEquals(3600, lga.phaseResolution());
    }
}
