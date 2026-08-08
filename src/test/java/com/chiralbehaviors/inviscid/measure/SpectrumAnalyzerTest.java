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
import static org.junit.Assert.assertTrue;

import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.measure.SpectrumAnalyzer.WindowFunction;

/**
 * @author halhildebrand
 *
 */
public class SpectrumAnalyzerTest {

    @Test
    public void constantRateRotorGivesSingleSpectralLine() {
        int n = 256;
        // bin-aligned frequencies (Necronomata.frequency quanta
        // convention): k1 = 4, k2 = 8 cycles over the n-sample window, so
        // omega = 2*pi*k/n lands exactly on a bin with no spectral
        // leakage from a fractional-bin rate. Rectangular window: this is
        // the exactly-periodic K=0 case (bead inviscid-0nx.7 / B.2's
        // >=0.95 peak-concentration requirement), not the
        // collision-perturbed case Hann is for.
        float frequency1 = 56.25f;
        float frequency2 = 112.5f;
        int expectedBin1 = SpectrumAnalyzer.expectedBinForFrequency(frequency1,
                                                                      n);
        int expectedBin2 = SpectrumAnalyzer.expectedBinForFrequency(frequency2,
                                                                      n);

        double[] spectrum1 = SpectrumAnalyzer.powerSpectrum(
        syntheticRotorSeries(n, frequency1), WindowFunction.RECTANGULAR);
        double[] spectrum2 = SpectrumAnalyzer.powerSpectrum(
        syntheticRotorSeries(n, frequency2), WindowFunction.RECTANGULAR);

        int peak1 = SpectrumAnalyzer.peakBin(spectrum1);
        int peak2 = SpectrumAnalyzer.peakBin(spectrum2);

        assertEquals("peak bin for frequency1", expectedBin1, peak1);
        assertEquals("peak bin for frequency2", expectedBin2, peak2);
        assertEquals("peak bin scales linearly with rate (ratio 2)", 2,
                     peak2 / peak1);

        assertTrue("bin " + peak1 + " should hold >=99% of power",
                   dominantFraction(spectrum1, peak1) >= 0.99);
        assertTrue("bin " + peak2 + " should hold >=99% of power",
                   dominantFraction(spectrum2, peak2) >= 0.99);
    }

    /**
     * SIGNIFICANT 1 (critique inviscid/critique-B1-fft-spectrum-inviscid-0nx.6):
     * the whole justification for choosing {@code exp(i*angle)} over
     * {@code sin(angle)} is that it discriminates rotation direction
     * instead of producing a symmetric mirror-bin pair - asserted in
     * javadoc, previously never tested. A negative-frequency rotor must
     * peak at bin {@code (N-|k|) mod N}, and its positive-frequency
     * mirror bin {@code k} must NOT also dominate (which is exactly what
     * a real-valued sin/cos analysis would produce instead).
     */
    @Test
    public void constantRateRotorNegativeFrequencyGivesMirroredSpectralLine() {
        int n = 256;
        float frequency = -56.25f;
        int expectedBin = SpectrumAnalyzer.expectedBinForFrequency(frequency,
                                                                     n);
        int mirrorBin = SpectrumAnalyzer.expectedBinForFrequency(-frequency,
                                                                   n);

        double[] spectrum = SpectrumAnalyzer.powerSpectrum(
        syntheticRotorSeries(n, frequency), WindowFunction.RECTANGULAR);

        assertEquals("negative-frequency rotor wraps to N-|k|", n - 4,
                     expectedBin);
        int peak = SpectrumAnalyzer.peakBin(spectrum);
        assertEquals("peak bin for negative frequency", expectedBin, peak);

        assertTrue("expected bin " + expectedBin + " should hold >=99% of power",
                   dominantFraction(spectrum, expectedBin) >= 0.99);
        assertTrue("positive-frequency mirror bin " + mirrorBin
                   + " must NOT dominate (would indicate no direction discrimination)",
                   spectrum[mirrorBin] < 0.01 * spectrum[expectedBin]);
    }

    /**
     * CRITICAL 1 (same critique): pins Hann's proven 2/3 peak-bin
     * concentration ceiling as documented, tested behavior, alongside
     * rectangular's ~100% for the same bin-aligned tone - so neither the
     * ceiling nor the fix regresses silently.
     */
    @Test
    public void windowFunctionControlsPeakBinConcentration() {
        int n = 256;
        float frequency = 56.25f; // bin 4
        float[] series = syntheticRotorSeries(n, frequency);
        int bin = SpectrumAnalyzer.expectedBinForFrequency(frequency, n);

        double[] rectangular = SpectrumAnalyzer.powerSpectrum(series,
                                                                WindowFunction.RECTANGULAR);
        double[] hann = SpectrumAnalyzer.powerSpectrum(series,
                                                         WindowFunction.HANN);

        double rectangularFraction = dominantFraction(rectangular, bin);
        double hannFraction = dominantFraction(hann, bin);

        assertTrue("rectangular window should concentrate >=99% of power in the home bin, was "
                   + rectangularFraction, rectangularFraction >= 0.99);
        assertTrue("Hann window caps concentration near 2/3 (was "
                   + hannFraction + ")",
                   Math.abs(hannFraction - 2.0 / 3.0) < 0.01);
    }

    /**
     * SIGNIFICANT 2 (same critique): the bin<->frequency conversion
     * utility, checked directly against known exact values (independent
     * of the FFT machinery) and, in the two tests above, against actually
     * observed rotor peaks for both positive and negative frequencies.
     */
    @Test
    public void expectedBinForFrequencyRoundTripsWithFrequencyForBin() {
        int n = 256;

        assertEquals(4, SpectrumAnalyzer.expectedBinForFrequency(56.25f, n));
        assertEquals(252,
                     SpectrumAnalyzer.expectedBinForFrequency(-56.25f, n));
        assertEquals(0, SpectrumAnalyzer.expectedBinForFrequency(0f, n));

        assertEquals(56.25, SpectrumAnalyzer.frequencyForBin(4, n), 1e-9);
        assertEquals(-56.25, SpectrumAnalyzer.frequencyForBin(252, n), 1e-9);
        assertEquals(0.0, SpectrumAnalyzer.frequencyForBin(0, n), 1e-9);
    }

    @Test
    public void spectrumIsDeterministicAcrossRuns() {
        int n = 128;
        float[] series = syntheticRotorSeries(n, 4.7f);

        double[] first = SpectrumAnalyzer.powerSpectrum(series.clone());
        double[] second = SpectrumAnalyzer.powerSpectrum(series.clone());

        assertEquals(first.length, second.length);
        for (int i = 0; i < first.length; i++) {
            assertEquals("bin " + i, first[i], second[i], 0.0);
        }
    }

    /**
     * Supplementary coverage (not one of the seven named failing tests, but
     * exercises the "records a per-member angle time series from a
     * Necronomata" half of the deliverable): a member seeded with a fixed
     * frequency free-rotates under repeated {@link Necronomata#step()}
     * calls, and {@link SpectrumAnalyzer#recordAngleSeries} must capture
     * exactly the angle trajectory {@code step()} produces - checked
     * against a manual re-derivation of the same recursion, not against
     * SpectrumAnalyzer's own output.
     */
    @Test
    public void recordAngleSeriesTracksNecronomataStep() {
        int steps = 32;
        int globalMemberIndex = 3;
        Necronomata automaton = new Necronomata(new Point3i(1, 1, 1));
        automaton.process((angle, frequency, deltaA, deltaF) -> frequency[globalMemberIndex] = 5f);

        float[] recorded = SpectrumAnalyzer.recordAngleSeries(automaton,
                                                                globalMemberIndex,
                                                                steps);

        float expectedAngle = 0f;
        float expectedFrequency = 5f;
        for (int t = 0; t < steps; t++) {
            assertEquals("angle at step " + t, expectedAngle, recorded[t],
                         1e-6f);
            expectedAngle += Necronomata.QUANTUM_RATE * expectedFrequency;
        }
    }

    private static double dominantFraction(double[] spectrum, int bin) {
        double total = 0;
        for (double v : spectrum) {
            total += v;
        }
        return spectrum[bin] / total;
    }

    /**
     * theta(n) = theta0 + (QUANTUM_RATE*frequency)*n: the linear phase
     * ramp a free rotor with {@code Necronomata.frequency == frequency}
     * produces under repeated {@code Necronomata.step()} - unbounded,
     * never wrapped mod 2*pi, exactly like the production angle array.
     * Generated fresh per sample (not accumulated) - deliberately does
     * NOT reproduce the float32 accumulation drift documented in {@link
     * SpectrumAnalyzer}'s class javadoc.
     */
    private static float[] syntheticRotorSeries(int n, float frequency) {
        Random random = new Random(11L);
        double theta0 = random.nextDouble() * 2 * Math.PI;
        double rate = Necronomata.QUANTUM_RATE * frequency;
        float[] series = new float[n];
        for (int i = 0; i < n; i++) {
            series[i] = (float) (theta0 + rate * i);
        }
        return series;
    }
}
