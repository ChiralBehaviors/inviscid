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
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Random;

import org.junit.Test;

/**
 * @author halhildebrand
 *
 */
public class FftTest {

    /**
     * THE ORACLE. The hand-rolled radix-2 FFT is never trusted on its own
     * output alone: every other test in this file exercises structural
     * properties (flatness, single-bin concentration, invertibility), but
     * none of them would catch a subtly wrong twiddle-factor sign or
     * butterfly wiring. This test cross-checks the FFT against a
     * completely independent O(N^2) direct-summation DFT on seeded random
     * input, which has no shared code path with {@link Fft} and so cannot
     * share its bugs.
     */
    @Test
    public void matchesNaiveDftOnSeededRandomInput() {
        int n = 256;
        Random random = new Random(42L);
        double[] re = new double[n];
        double[] im = new double[n];
        for (int i = 0; i < n; i++) {
            re[i] = random.nextDouble() * 2 - 1;
            im[i] = random.nextDouble() * 2 - 1;
        }

        double[] expectedRe = new double[n];
        double[] expectedIm = new double[n];
        naiveDft(re, im, expectedRe, expectedIm);

        double[] actualRe = re.clone();
        double[] actualIm = im.clone();
        Fft.fft(actualRe, actualIm);

        for (int k = 0; k < n; k++) {
            assertEquals("re[" + k + "]", expectedRe[k], actualRe[k], 1e-9);
            assertEquals("im[" + k + "]", expectedIm[k], actualIm[k], 1e-9);
        }
    }

    @Test
    public void impulseGivesFlatSpectrum() {
        int n = 64;
        double[] re = new double[n];
        double[] im = new double[n];
        re[0] = 1.0;

        Fft.fft(re, im);

        for (int k = 0; k < n; k++) {
            double power = re[k] * re[k] + im[k] * im[k];
            assertEquals("power at bin " + k, 1.0, power, 1e-9);
        }
    }

    @Test
    public void pureToneGivesSingleBin() {
        int n = 64;
        int bin = 5;
        double[] re = new double[n];
        double[] im = new double[n];
        for (int i = 0; i < n; i++) {
            double theta = 2 * Math.PI * bin * i / n;
            re[i] = Math.cos(theta);
            im[i] = Math.sin(theta);
        }

        Fft.fft(re, im);

        double[] power = new double[n];
        double total = 0;
        for (int k = 0; k < n; k++) {
            power[k] = re[k] * re[k] + im[k] * im[k];
            total += power[k];
        }

        assertTrue("bin " + bin + " should hold >99% of total power",
                   power[bin] > 0.99 * total);
    }

    @Test
    public void rejectsNonPowerOfTwoLength() {
        double[] re = new double[100];
        double[] im = new double[100];
        try {
            Fft.fft(re, im);
            fail("expected IllegalArgumentException for non-power-of-two length");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    @Test
    public void roundTripsThroughInverse() {
        int n = 128;
        Random random = new Random(7L);
        double[] originalRe = new double[n];
        double[] originalIm = new double[n];
        for (int i = 0; i < n; i++) {
            originalRe[i] = random.nextDouble() * 2 - 1;
            originalIm[i] = random.nextDouble() * 2 - 1;
        }

        double[] re = originalRe.clone();
        double[] im = originalIm.clone();
        Fft.fft(re, im);
        Fft.ifft(re, im);

        for (int i = 0; i < n; i++) {
            assertEquals("re[" + i + "]", originalRe[i], re[i], 1e-9);
            assertEquals("im[" + i + "]", originalIm[i], im[i], 1e-9);
        }
    }

    /**
     * Independent O(N^2) direct-summation DFT, sharing no code with
     * {@link Fft}. X[k] = sum_n x[n] * exp(-2*pi*i*k*n/N).
     */
    private static void naiveDft(double[] re, double[] im, double[] outRe,
                                  double[] outIm) {
        int n = re.length;
        for (int k = 0; k < n; k++) {
            double sumRe = 0;
            double sumIm = 0;
            for (int t = 0; t < n; t++) {
                double angle = -2 * Math.PI * k * t / n;
                double c = Math.cos(angle);
                double s = Math.sin(angle);
                sumRe += re[t] * c - im[t] * s;
                sumIm += re[t] * s + im[t] * c;
            }
            outRe[k] = sumRe;
            outIm[k] = sumIm;
        }
    }
}
