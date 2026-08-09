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

package com.chiralbehaviors.inviscid;

import java.util.function.Consumer;

import javax.vecmath.Point3i;

/**
 * The measurement read seam over a jitterbug lattice's conserved quanta
 * (bead inviscid-ckn / inviscid-0nx.21). Derived from a call-site census
 * of {@code ConservationAudit}, {@code StructureFactor.coarseGrainedField},
 * and {@code SpectrumAnalyzer.recordAngleSeries} — not speculation. Exposes
 * exact {@code long} quanta per slot, never a {@code float[]}, so a
 * table-driven substrate (the formal LGA) never has to materialise floats
 * it does not have.
 *
 * <p>LAYOUT CONTRACT (both implementations MUST honour, and the
 * conformance test pins): slots are 30 per cell (5 cubes x 6 members),
 * linearised as {@code 30 * ((i*extent.y + j)*extent.z + k)}; within a
 * cell, {@code cube = local/6}, {@code member = local%6}; cells are the
 * even-parity sublattice ({@code (i+j+k)%2==0}). This contract — not a
 * method — is what makes an instrument's output identical across
 * substrates.
 *
 * <p>Deliberately NOT split into a separate phase-only interface: there
 * are exactly two implementations, both of which implement every method
 * here. See T2 {@code inviscid/design-ckn-lattice-seam.md} §2 for the
 * full rationale, including why {@link #isExactAt(int)} earns its place
 * despite having no float-storage meaning for an integer-backed
 * substrate (it must report "structurally exact", never present a
 * tautology as evidence).
 *
 * @author halhildebrand
 */
public interface QuantaField {

    /** Defensive copy of the lattice extent. */
    Point3i extent();

    /** Total slot count == 30 * extent.x * extent.y * extent.z. */
    int slotCount();

    /**
     * Exact signed quanta in {@code slot}. Integer by definition of the
     * conserved quantity; never a rounded view of a fractional value.
     */
    long quantaAt(int slot);

    /**
     * Whether {@code slot}'s NATIVE representation holds an exact
     * integer. Necronomata: {@code Math.rint(frequency[slot]) ==
     * frequency[slot]}. Integer-backed substrates: constant {@code true}
     * — and callers MUST report that as "structurally exact", never as
     * "corruption check passed" (T2 seam memo §2, risk R1).
     */
    boolean isExactAt(int slot);

    /**
     * Member phase in radians, wrapped to {@code [0, 2*pi)}, at the
     * FINEST resolution the substrate holds. Necronomata: the raw
     * {@code angle[slot]}. A future formal-LGA implementation under the
     * user's 2A cadence decision (T2 {@code
     * inviscid/analysis-73v-spectral-conversion-and-cadence.md} §5) MUST
     * return its fine sub-bin phase accumulator here, NEVER a
     * contact-bin centre — a bin-centre read injects a deterministic
     * quantisation-sideband artifact into any angle-spectrum instrument
     * that samples this accessor.
     */
    float phaseAt(int slot);

    /**
     * Phase resolution in steps per revolution, matching
     * {@link #phaseAt}'s quantisation. Necronomata: 3600
     * ({@code Necronomata.PHASE_RESOLUTION}). Any bin&lt;-&gt;quanta
     * conversion MUST read this accessor, never a hardwired constant —
     * see {@code SpectrumAnalyzer}'s {@code expectedBinForFrequency} /
     * {@code frequencyForBin}, which is inviscid-73v's hardening target.
     */
    int phaseResolution();

    /** Visits exactly the even-parity cells, in Necronomata.forEach order. */
    void forEachCell(Consumer<? super Point3i> action);

    /** Flat base index of {@code cell}'s 30-slot block. */
    int indexOfCell(Point3i cell);
}
