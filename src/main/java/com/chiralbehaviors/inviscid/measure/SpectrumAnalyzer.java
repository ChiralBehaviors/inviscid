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

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.QuantaField;

/**
 * Records a per-member angle time series from a {@link Necronomata} and
 * turns it into a deterministic power spectrum.
 *
 * <h2>The ramp-vs-sinusoid choice</h2>
 * A free rotor's angle is a linear phase ramp - {@code
 * Necronomata.step()} advances {@code angle[m] += QUANTUM_RATE *
 * frequency[m]} every tick (wrapped into {@code [0, 2*pi)} since
 * inviscid-vb9) - so a K=0 baseline member's recorded series is
 * {@code theta(n) = (theta0 + omega*n) mod 2*pi}: a sawtooth, the
 * wrapped image of a ramp. Three ways to turn that into something
 * spectrally meaningful were considered; {@code exp(i*angle)}, the
 * complex analytic-signal mapping, was chosen (the wrap strengthens the
 * case: a raw wrapped series has periodic sawtooth discontinuities that
 * {@code exp(i*angle)}, being exactly mod-2*pi invariant, never sees):
 *
 * <ul>
 * <li><b>Raw angle.</b> Rejected outright: a ramp is not periodic within a
 * finite analysis window, so its DFT smears broadband across many bins
 * (looks like a chirp) even though the physical rotor has one well-defined
 * frequency. This is "the ramp problem" this class exists to solve.</li>
 *
 * <li><b>{@code sin(angle)} / {@code cos(angle)}.</b> Periodic and
 * tempting, but real-valued: a real signal's DFT is conjugate-symmetric
 * ({@code X[N-k] == conj(X[k])}), so one physical rotation frequency
 * produces <i>two</i> equal-magnitude bins (k and N-k), not the single
 * dominant line the K=0 baseline should show. Making "one line" true would
 * require an extra half-spectrum convention, and multiple simultaneous
 * rotors would still waste half the spectrum on mirror images that can
 * collide with each other under folding.</li>
 *
 * <li><b>Differencing</b> ({@code angle[n] - angle[n-1]}). For a pure free
 * rotor the difference is <i>constant</i> ({@code == QUANTUM_RATE *
 * frequency}), i.e. differencing collapses exactly the quantity a spectrum
 * is supposed to reveal - the rotation frequency - into DC (bin 0),
 * regardless of what that frequency actually is. Useful later as a
 * collision-event / residual detector once rotors are no longer perfectly
 * free, wrong as a baseline spectral estimator.</li>
 *
 * <li><b>{@code exp(i*angle) = cos(angle) + i*sin(angle)}, the choice made
 * here.</b> Treats the angle as the instantaneous phase of a complex
 * analytic signal on the unit circle. For the linear ramp this becomes
 * {@code exp(i*theta0) * exp(i*omega*n)} - an <i>exact</i> pure complex
 * exponential at digital frequency {@code omega}, landing on exactly one
 * non-mirrored bin with no detrending step and no phantom mirror image,
 * and - unlike {@code sin}/{@code cos} - it discriminates rotation
 * <i>direction</i>: a negative-rate rotor lands at bin {@code (N-|k|) mod
 * N}, not at the same bin as its positive-rate mirror (verified by
 * {@code SpectrumAnalyzerTest#constantRateRotorNegativeFrequencyGivesMirroredSpectralLine}).
 * This is the standard analytic-signal / instantaneous-phase trick, it
 * requires no assumption that the ramp stay perfectly linear (a
 * collision-perturbed phase trajectory still maps to a well-defined
 * instantaneous complex phase), and it is what makes "one dominant
 * spectral line per constant-rate rotor" true by construction rather than
 * by convention.</li>
 * </ul>
 *
 * <h2>Window choice</h2>
 * {@link #powerSpectrum(float[], WindowFunction)} windows the {@code
 * exp(i*angle)} series before the FFT. Two windows are exposed:
 *
 * <ul>
 * <li>{@link WindowFunction#RECTANGULAR} - no windowing (multiply by 1).
 * Correct for the K=0 baseline: a bin-aligned pure tone is already exactly
 * periodic in the analysis window, so there is no edge discontinuity to
 * suppress, and rectangular is the only one of the two windows that
 * concentrates (near) all of a bin-aligned tone's power into its single
 * bin (empirically &gt;99%; see {@code
 * SpectrumAnalyzerTest#windowFunctionControlsPeakBinConcentration}). Use
 * this for any K=0 acceptance check that thresholds peak-bin
 * concentration (e.g. bead inviscid-0nx.7 / B.2's &gt;=0.95 requirement) -
 * {@link WindowFunction#HANN} cannot reach that threshold, see below.</li>
 *
 * <li>{@link WindowFunction#HANN} ({@code
 * 0.5*(1-cos(2*pi*n/(N-1)))}). <b>Caps peak-bin concentration at exactly
 * 2/3 for any bin-aligned tone</b> - this is a provable property of Hann's
 * three-point (0.25, 0.5, 0.25) frequency-domain smearing across bins
 * k-1, k, k+1, not an artifact of this implementation (confirmed
 * empirically at N=256 for k=4 and k=8: 0.664-0.667). An earlier version
 * of this javadoc claimed Hann was "close to a no-op" for the K=0 case;
 * that was wrong - it leaks a third of the tone's power out of its home
 * bin unconditionally. Hann still earns its place once collision physics
 * (bead inviscid-0nx.14) perturb the phase trajectory: the series is then
 * no longer a perfect ramp between window edges, and Hann's low sidelobes
 * keep that leakage from swamping genuine collision-broadened lines in a
 * way rectangular's high sidelobes would not. Choose it for that
 * non-periodic case (Phase A / bead inviscid-0nx.9's S(k,omega)), not for
 * K=0 peak-concentration checks.</li>
 * </ul>
 *
 * {@link #powerSpectrum(float[])} - the no-window-argument convenience
 * overload - defaults to {@link WindowFunction#HANN} for that same
 * forward-looking reason; callers that need the K=0
 * &gt;=0.95-concentration behavior must call {@link
 * #powerSpectrum(float[], WindowFunction)} with {@link
 * WindowFunction#RECTANGULAR} explicitly.
 *
 * <h2>Float32 angle-accumulation precision ceiling (RETIRED by inviscid-vb9)</h2>
 * {@code Necronomata.angle} is {@code float} (32-bit). As of inviscid-vb9,
 * {@code step()} wraps {@code angle} into {@code [0, 2*pi)} (floor-mod, in
 * double precision) after every tick, so the value is always bounded to
 * one revolution - it can no longer walk arbitrarily far from zero, and
 * the per-tick rounding error is now just one bounded mod-reduction plus
 * one add on a value under {@code 2*pi}, which does <b>not</b> grow with
 * warm-up tick count the way unbounded accumulation used to. The residual
 * error at any tick count is on the order of a single float32 ULP near
 * {@code 2*pi} (~2.4e-7 rad), not an accumulating drift.
 *
 * <p><b>Historical note</b> (true of the code before inviscid-vb9, kept
 * for context - do not use these numbers to reason about current
 * behavior): before the wrap, {@code angle} accumulated unbounded and
 * every tick after the first pushed the value further from zero, so
 * float32's fixed 24-bit mantissa meant the ULP grew with magnitude and a
 * fixed-size {@code deltaA} increment eventually rounded away partially,
 * then completely. Measured (actual Java {@code float} accumulation) as
 * the relative error of the locally-observed rotation rate within a
 * 4096-sample analysis window, as a function of prior warm-up ticks:
 *
 * <pre>
 *     10,000 prior ticks:        0.016% rate error
 *     1,000,000 prior ticks:     0.72%  rate error
 *     2^23 (~8,388,608) ticks:  -10.5%  rate error   (universal float32 ULP-crossing threshold)
 *     50,000,000 ticks:         -100%   rate error   (angle stops advancing entirely)
 * </pre>
 *
 * That {@code 2^23}-tick cliff was frequency-independent and silent: it
 * would present in a spectrum as line broadening or a spurious frequency
 * shift, indistinguishable from genuine collision broadening (the
 * "instrument contaminates the measurement" trap) - which is exactly why
 * inviscid-vb9 fixed it at the source ({@code Necronomata.step()}) rather
 * than working around it here. {@code SpectrumAnalyzerTest}'s synthetic
 * rotor series are generated fresh per sample (not accumulated), so they
 * never exercised this failure mode either way.
 *
 * @author halhildebrand
 */
public final class SpectrumAnalyzer {

    /**
     * Window applied to the {@code exp(i*angle)} series before the FFT.
     * See the class javadoc's "Window choice" section for the tradeoff.
     */
    public enum WindowFunction {
        /** No windowing (multiply by 1 everywhere). */
        RECTANGULAR,
        /** {@code 0.5*(1-cos(2*pi*n/(N-1)))}. Caps peak-bin concentration
         * at 2/3 for a bin-aligned tone; see class javadoc. */
        HANN
    }

    /**
     * A recorded angle series carrying the {@link SpectralCadence} it was
     * sampled at (bead inviscid-ckn / inviscid-0nx.21, 73v option 1D):
     * unlike a bare {@code float[]}, a {@code PhaseSeries} makes it
     * structurally impossible to analyse a series against the wrong
     * phase resolution -- the analyzer reads {@code cadence} instead of
     * assuming {@code Necronomata.PHASE_RESOLUTION}.
     */
    public record PhaseSeries(float[] angles, SpectralCadence cadence) {
    }

    private SpectrumAnalyzer() {
    }

    /**
     * QuantaField-typed, read-only phase sampler (bead inviscid-ckn /
     * inviscid-0nx.21, 73v option 1D): records {@code steps} consecutive
     * {@link QuantaField#phaseAt(int)} readings for slot {@code
     * globalMemberIndex}, advancing the substrate between samples via
     * {@code advanceOneTick} -- a callback the CALLER's own tick loop
     * supplies (e.g. {@code Necronomata::step} for a free rotor, or
     * {@code driver::tick} for an audited run), never a {@code step()}
     * this method owns itself. This is what "read-only" means here: the
     * method only ever calls {@link QuantaField#phaseAt(int)}; advancing
     * is entirely the caller's concern, matching the design memo's "do
     * not put step() on the seam" rule. The first recorded value is the
     * member's phase <i>before</i> any of these {@code steps} advances
     * run, matching {@link #recordAngleSeries}'s existing convention.
     */
    public static PhaseSeries sampleSeries(QuantaField field,
                                            int globalMemberIndex, int steps,
                                            Runnable advanceOneTick) {
        if (steps <= 0) {
            throw new IllegalArgumentException("steps must be positive, was "
                                                + steps);
        }
        float[] series = new float[steps];
        for (int t = 0; t < steps; t++) {
            series[t] = field.phaseAt(globalMemberIndex);
            advanceOneTick.run();
        }
        return new PhaseSeries(series, SpectralCadence.perTick(field));
    }

    /**
     * The expected FFT bin (0..n-1, wrapped) for a member rotating at
     * {@code frequency} quanta (the same signed quantity {@code
     * Necronomata.frequency} holds - {@code deltaA == QUANTUM_RATE *
     * frequency} radians/step). Negative frequencies wrap to the upper
     * half of the bin range ({@code (N-|k|) mod N}), matching where
     * {@link #powerSpectrum} actually places a negative-rate rotor's
     * single line (see the ramp-vs-sinusoid section on direction
     * discrimination). Inverse of {@link #frequencyForBin(int, int)}.
     *
     * @deprecated hardwires {@code Necronomata.PHASE_RESOLUTION} (bead
     *             inviscid-73v). Arithmetically identical delegation to
     *             {@code new SpectralCadence(Necronomata.PHASE_RESOLUTION, 1).binFor(...)}
     *             -- prefer {@link SpectralCadence#binFor(double, int)}
     *             directly, reading {@code phaseResolution} from a
     *             {@link QuantaField} instead of assuming 3600.
     */
    @Deprecated
    public static int expectedBinForFrequency(float frequency, int n) {
        return new SpectralCadence(Necronomata.PHASE_RESOLUTION, 1).binFor(frequency,
                                                                             n);
    }

    /**
     * The signed frequency (quanta, {@code Necronomata.frequency}
     * convention) whose rotor would peak at bin {@code bin} of an
     * {@code n}-sample spectrum. Bins in the upper half ({@code bin >
     * n/2}) are interpreted as negative frequencies, per the standard FFT
     * negative-frequency convention. Approximate inverse of {@link
     * #expectedBinForFrequency(float, int)} (exact when the forward
     * mapping did not need rounding).
     *
     * @deprecated hardwires {@code Necronomata.PHASE_RESOLUTION} (bead
     *             inviscid-73v). Arithmetically identical delegation to
     *             {@code new SpectralCadence(Necronomata.PHASE_RESOLUTION, 1).quantaFor(...)}
     *             -- prefer {@link SpectralCadence#quantaFor(int, int)}
     *             directly.
     */
    @Deprecated
    public static double frequencyForBin(int bin, int n) {
        return new SpectralCadence(Necronomata.PHASE_RESOLUTION,
                                    1).quantaFor(bin, n);
    }

    /**
     * Peak bin: the index of the largest value in {@code powerSpectrum}.
     * Ties resolve to the lowest index.
     */
    public static int peakBin(double[] powerSpectrum) {
        int peak = 0;
        for (int i = 1; i < powerSpectrum.length; i++) {
            if (powerSpectrum[i] > powerSpectrum[peak]) {
                peak = i;
            }
        }
        return peak;
    }

    /**
     * {@link #powerSpectrum(float[], WindowFunction)} with {@link
     * WindowFunction#HANN} - see the class javadoc's "Window choice"
     * section for why that is the default and when callers must instead
     * pass {@link WindowFunction#RECTANGULAR} explicitly.
     */
    public static double[] powerSpectrum(float[] angleSeries) {
        return powerSpectrum(angleSeries, WindowFunction.HANN);
    }

    /**
     * The deterministic power spectrum of {@code angleSeries}: {@code
     * exp(i*angle[n])} windowed with {@code window}, then FFT'd; {@code
     * power[k] = re[k]^2 + im[k]^2}. Length must be a power of two (see
     * {@link Fft}).
     */
    public static double[] powerSpectrum(float[] angleSeries,
                                          WindowFunction window) {
        int n = angleSeries.length;
        double[] w = window(window, n);
        double[] re = new double[n];
        double[] im = new double[n];
        for (int i = 0; i < n; i++) {
            double theta = angleSeries[i];
            re[i] = Math.cos(theta) * w[i];
            im[i] = Math.sin(theta) * w[i];
        }

        Fft.fft(re, im);

        double[] power = new double[n];
        for (int k = 0; k < n; k++) {
            power[k] = re[k] * re[k] + im[k] * im[k];
        }
        return power;
    }

    /**
     * Records {@code steps} consecutive angle values for member {@code
     * globalMemberIndex} (the flat 0..(30*cellCount-1) index {@link
     * Necronomata#indexOfCell} addresses into), advancing {@code
     * automaton} one tick between each sample via {@link
     * Necronomata#step()}. The first recorded value is the member's angle
     * <i>before</i> any of these {@code steps} ticks run.
     *
     * <p>Uses {@link Necronomata#process(Necronomata.Processor)} strictly
     * read-only, per that method's javadoc contract - it never writes
     * {@code angle}/{@code deltaA} (or any other array).
     *
     * <p>See the class javadoc's "Float32 angle-accumulation precision
     * ceiling" section: since inviscid-vb9, {@code automaton}'s angle is
     * wrapped into {@code [0, 2*pi)} every tick, so its per-tick rounding
     * error stays bounded regardless of prior warm-up length (that
     * section's error table is retired, historical-only behavior).
     *
     * <p>Necronomata-typed delegation (bead inviscid-ckn /
     * inviscid-0nx.21) to {@link #sampleSeries(QuantaField, int, int,
     * Runnable)}, kept byte-identical so Phase B stays untouched -- see
     * {@link #sampleSeries} for the substrate-agnostic, read-only
     * replacement.
     */
    public static float[] recordAngleSeries(Necronomata automaton,
                                             int globalMemberIndex,
                                             int steps) {
        return sampleSeries(automaton, globalMemberIndex, steps,
                             automaton::step).angles();
    }

    private static double[] hannWindow(int n) {
        double[] w = new double[n];
        if (n == 1) {
            w[0] = 1.0;
            return w;
        }
        for (int i = 0; i < n; i++) {
            w[i] = 0.5 * (1 - Math.cos(2 * Math.PI * i / (n - 1)));
        }
        return w;
    }

    private static double[] rectangularWindow(int n) {
        double[] w = new double[n];
        java.util.Arrays.fill(w, 1.0);
        return w;
    }

    private static double[] window(WindowFunction fn, int n) {
        switch (fn) {
        case RECTANGULAR:
            return rectangularWindow(n);
        case HANN:
            return hannWindow(n);
        default:
            throw new IllegalArgumentException("unhandled window: " + fn);
        }
    }
}
