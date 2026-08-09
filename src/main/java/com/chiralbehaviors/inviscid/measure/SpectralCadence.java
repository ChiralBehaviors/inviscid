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

import com.chiralbehaviors.inviscid.QuantaField;

/**
 * A substrate's phase-resolution / sampling-cadence pair (bead
 * inviscid-73v, option 1B, T2
 * analysis-73v-spectral-conversion-and-cadence.md §2.3). Parameterises
 * the bin&lt;-&gt;quanta conversion that {@code SpectrumAnalyzer}'s
 * {@code expectedBinForFrequency}/{@code frequencyForBin} previously
 * hardwired to {@code Necronomata.PHASE_RESOLUTION} (3600) -- any
 * {@link QuantaField} implementation's {@link QuantaField#phaseResolution()}
 * now drives the conversion instead, so a future substrate with a
 * different {@code P} (or a different sampling {@code stride}) cannot
 * silently mis-convert bins.
 *
 * <p><b>Derivation (73v §2.1-2.2, not a free design choice).</b> For a
 * substrate with phase resolution {@code P} sampled every {@code stride}
 * ticks, a member with conserved quanta {@code q} has:
 * <pre>
 *   cycles per tick    = q / P
 *   cycles per sample  = q * stride / P
 *   FFT bin (length n) = round(q * stride * n / P), wrapped mod n
 *   inverse            = q = signedBin * P / (stride * n)
 * </pre>
 * Writing {@code P = 2^a * odd(P)}: exact bin alignment for every
 * integer {@code q} requires {@code stride == odd(P)} (the minimal
 * stride), and the Nyquist quanta bound at that stride is
 * {@code P/(2*stride) == 2^(a-1)}. At {@code P=3600} this reproduces the
 * committed {@code BaselineSpectrumHarness} constants ({@code
 * STRIDE=225}, {@code |quanta| < 8}) exactly -- independent validation
 * of the derivation, not a coincidence.
 *
 * @param phaseResolution steps per revolution (e.g. 3600 for both
 *                        Necronomata and the formal LGA under the 2A
 *                        cadence decision)
 * @param stride          ticks per recorded sample (1 for a per-tick
 *                        series; {@link #alignmentStride()} for exact
 *                        FFT-bin alignment)
 * @author halhildebrand
 */
public record SpectralCadence(int phaseResolution, int stride) {

    public SpectralCadence {
        if (phaseResolution <= 0) {
            throw new IllegalArgumentException("phaseResolution must be positive, was "
                                                + phaseResolution);
        }
        if (stride <= 0) {
            throw new IllegalArgumentException("stride must be positive, was "
                                                + stride);
        }
    }

    /** Per-tick cadence (stride 1) for {@code field}'s phase resolution. */
    public static SpectralCadence perTick(QuantaField field) {
        return new SpectralCadence(field.phaseResolution(), 1);
    }

    /**
     * Alignment cadence: {@code stride == oddPart(phaseResolution)}, the
     * minimal stride at which every integer quanta count binds to an
     * exact FFT bin with no rounding (73v §2.2).
     */
    public static SpectralCadence aligned(QuantaField field) {
        int p = field.phaseResolution();
        return new SpectralCadence(p, oddPart(p));
    }

    /** {@code oddPart(phaseResolution)} -- the minimal exact-alignment stride. */
    public int alignmentStride() {
        return oddPart(phaseResolution);
    }

    /**
     * The largest quanta magnitude that does not alias at this cadence:
     * {@code phaseResolution / (2*stride)}.
     */
    public int nyquistQuantaBound() {
        return phaseResolution / (2 * stride);
    }

    /** Cycles of rotation per tick for a member with {@code quanta}. */
    public double cyclesPerTick(double quanta) {
        return quanta / phaseResolution;
    }

    /** Converts an angular rate measured per SAMPLE to radians per TICK. */
    public double omegaRadPerTick(double omegaPerSample) {
        return omegaPerSample / stride;
    }

    /**
     * The expected FFT bin (0..n-1, wrapped) for a member advancing
     * {@code quanta} phase-steps per sample (== {@code quanta} raw
     * conserved quanta when {@code stride == 1}; the caller pre-folds
     * {@code stride} into {@code quanta} for a coarser sampling cadence,
     * matching {@code BaselineSpectrumHarness}'s existing convention).
     */
    public int binFor(double quanta, int n) {
        if (n <= 0) {
            throw new IllegalArgumentException("n must be positive, was " + n);
        }
        double k = quanta * stride * n / phaseResolution;
        long rounded = Math.round(k);
        long wrapped = ((rounded % n) + n) % n;
        return (int) wrapped;
    }

    /** Inverse of {@link #binFor(double, int)} (exact when no rounding occurred). */
    public double quantaFor(int bin, int n) {
        if (n <= 0) {
            throw new IllegalArgumentException("n must be positive, was " + n);
        }
        int normalized = ((bin % n) + n) % n;
        int signed = normalized > n / 2 ? normalized - n : normalized;
        return (double) signed * phaseResolution / ((double) stride * n);
    }

    private static int oddPart(int p) {
        while (p % 2 == 0) {
            p /= 2;
        }
        return p;
    }
}
