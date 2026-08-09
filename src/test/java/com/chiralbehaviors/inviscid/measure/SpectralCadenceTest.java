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
import static org.junit.Assert.assertThrows;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;

/**
 * Conformance tests for {@link SpectralCadence} (bead inviscid-73v,
 * option 1B, T2 analysis-73v-spectral-conversion-and-cadence.md §2).
 *
 * @author halhildebrand
 */
public class SpectralCadenceTest {

    private static Necronomata freshNecronomata() {
        return new Necronomata(new Point3i(4, 4, 4));
    }

    @Test
    public void perTickIsStrideOneAtTheFieldsPhaseResolution() {
        SpectralCadence cadence = SpectralCadence.perTick(freshNecronomata());
        assertEquals(3600, cadence.phaseResolution());
        assertEquals(1, cadence.stride());
    }

    /**
     * Reproduces the committed {@code BaselineSpectrumHarness} constants
     * exactly (73v §2.2's independent-validation claim): {@code
     * STRIDE=225}, Nyquist bound 8, at {@code P=3600}.
     */
    @Test
    public void alignedReproducesTheCommittedPhaseAConstantsAt3600() {
        SpectralCadence cadence = SpectralCadence.aligned(freshNecronomata());
        assertEquals(3600, cadence.phaseResolution());
        assertEquals(225, cadence.stride());
        assertEquals(225, cadence.alignmentStride());
        assertEquals(8, cadence.nyquistQuantaBound());
    }

    /**
     * Arithmetic-identity pin (memo: "arithmetically identical"):
     * {@code new SpectralCadence(3600, 1).binFor(f, n) ==
     * SpectrumAnalyzer.expectedBinForFrequency(f, n)} for every {@code f}
     * in a representative range, at several power-of-two {@code n}.
     */
    @Test
    public void binForAtStrideOneMatchesExpectedBinForFrequencyExactly() {
        SpectralCadence cadence = new SpectralCadence(3600, 1);
        for (int n : new int[] { 4, 8, 16, 32 }) {
            for (float f = -20f; f <= 20f; f += 1f) {
                assertEquals("f=" + f + " n=" + n,
                             SpectrumAnalyzer.expectedBinForFrequency(f, n),
                             cadence.binFor(f, n));
            }
        }
    }

    /** Arithmetic-identity pin, inverse direction. */
    @Test
    public void quantaForAtStrideOneMatchesFrequencyForBinExactly() {
        SpectralCadence cadence = new SpectralCadence(3600, 1);
        for (int n : new int[] { 4, 8, 16, 32 }) {
            for (int bin = 0; bin < n; bin++) {
                assertEquals("bin=" + bin + " n=" + n,
                             SpectrumAnalyzer.frequencyForBin(bin, n),
                             cadence.quantaFor(bin, n), 0.0);
            }
        }
    }

    @Test
    public void cyclesPerTickIsQuantaOverPhaseResolution() {
        SpectralCadence cadence = new SpectralCadence(3600, 1);
        assertEquals(0.25, cadence.cyclesPerTick(900), 0.0);
    }

    @Test
    public void omegaRadPerTickDividesBySampleStride() {
        SpectralCadence cadence = new SpectralCadence(3600, 225);
        assertEquals(1.0, cadence.omegaRadPerTick(225.0), 0.0);
    }

    @Test
    public void rejectsNonPositivePhaseResolutionOrStride() {
        assertThrows(IllegalArgumentException.class,
                     () -> new SpectralCadence(0, 1));
        assertThrows(IllegalArgumentException.class,
                     () -> new SpectralCadence(3600, 0));
    }
}
