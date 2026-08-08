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

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.OptionalDouble;
import java.util.Random;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.measure.SpectrumAnalyzer.WindowFunction;

/**
 * Generates the K=0 (collision-free) baseline spectrum golden artifact for
 * bead inviscid-0nx.7 (B.2). Builds a headless {@link Necronomata} with a
 * seeded integer frequency (quanta) field, runs it with NO collision
 * handling ({@link Necronomata#process(Point3i)} is still the no-op it was
 * left as by inviscid-0nx.14's prerequisite), and records a per-member
 * power spectrum via {@link SpectrumAnalyzer}. Deliberately sequenced
 * ahead of any collision rule: a baseline measured after collisions land
 * would not be a baseline, and this ordering is what makes the later
 * "collision-broadened lines" claim falsifiable.
 *
 * <h2>The bin-alignment problem and its fix (STRIDE)</h2>
 * A free rotor with {@code Necronomata.frequency == quanta} advances at
 * {@code deltaA = QUANTUM_RATE * quanta} radians/tick, i.e. digital
 * frequency {@code quanta / PHASE_RESOLUTION} (= quanta/3600) cycles per
 * tick. Sampling once per tick over an {@code N}-sample window (N a power
 * of two, required by {@link Fft}) would land the tone at bin {@code
 * quanta*N/3600}. Since {@code 3600 = 2^4 * 3^2 * 5^2 = 16*225} and a
 * power-of-two {@code N} can never supply the factor of {@code 225}, that
 * bin is a NON-INTEGER for essentially every small integer quanta value
 * (e.g. quanta=1) - the tone leaks across bins (sinc leakage under a
 * rectangular window) and the required &gt;=0.95 peak-bin concentration is
 * mathematically unreachable that way, regardless of N.
 *
 * <p>Fix: sample every {@link #STRIDE} ticks instead of every tick (advance
 * the automaton {@code STRIDE} steps between recorded samples). The
 * recorded series' digital frequency becomes {@code quanta*STRIDE/3600}
 * cycles/sample, landing at bin {@code quanta*STRIDE*N/3600}. Choosing
 * {@code STRIDE = 225} (the full odd part of 3600) and {@code N} any
 * power of two {@code >= 16} makes {@code STRIDE*N/3600 == N/16} an EXACT
 * INTEGER, so the bin {@code quanta*(N/16)} is bin-aligned for EVERY
 * integer quanta value, not just specially chosen ones. With {@link
 * #FFT_LENGTH} = 256, bin spacing per unit quanta is {@code 256/16 = 16} -
 * well separated across the seeded range (empirically confirmed: 0/110
 * bin mismatches across every nonzero-quanta member for the default
 * parameters).
 *
 * <p>Consequence for bin<->frequency conversion: because the recorded
 * series' effective rate is {@code quanta*STRIDE}, not {@code quanta}, the
 * correct call into the shared {@link SpectrumAnalyzer#expectedBinForFrequency}
 * convention is {@code expectedBinForFrequency(quanta * STRIDE, N)}.
 *
 * <h2>Float32 precision budget - RETIRED by inviscid-vb9, kept for
 * historical context</h2>
 * As of inviscid-vb9, {@code Necronomata.step()} wraps {@code angle} into
 * {@code [0, 2*pi)} every tick, so the per-tick rounding error stays
 * bounded regardless of total tick count - the sweep below (run against
 * the pre-wrap unbounded accumulation) no longer describes current
 * behavior; peak-bin concentration at the default parameters is now
 * &gt;=0.9999 for every nonzero-quanta member (see the regenerated golden
 * artifact), not the 0.9914 worst case this section originally measured.
 * {@link #FFT_LENGTH} = 256 remains the chosen value (no need to grow it
 * now that the precision ceiling that motivated keeping it small is
 * gone), but the reasoning below is retained only as a record of why it
 * was originally picked.
 *
 * <p>The bin-alignment fix's own total tick count is {@code
 * (N-1)*STRIDE}, which GROWS with {@code N} (since the minimal {@code
 * STRIDE} is a fixed 225 for any {@code N >= 16}) - before inviscid-vb9,
 * picking a large FFT length "for better resolution" directly fought the
 * float32 angle-accumulation ceiling documented in {@link
 * SpectrumAnalyzer}'s class javadoc. That javadoc's retired error table
 * was measured for a short {@code recordAngleSeries} window taken AFTER a
 * long prior warm-up (only the window's own few thousand ticks
 * accumulated further rounding against an already-large angle offset);
 * THIS harness had no separate warm-up phase - the recorded window itself
 * spanned the entire accumulation from angle&asymp;0 to
 * angle&asymp;{@code quanta*QUANTUM_RATE*(N-1)*STRIDE}, so the
 * instantaneous rate drifted (a float32-rounding-induced chirp, not a
 * stationary tone) across the SAME window being analyzed. Empirically
 * swept (not guessed) with {@code STRIDE=225}, extent (2,2,2), quanta in
 * [-5,5], against the PRE-WRAP code:
 *
 * <pre>
 *     N=128  (  28,575 ticks): worst nonzero-quanta peak fraction 0.9996
 *     N=256  (  57,375 ticks): worst nonzero-quanta peak fraction 0.9914  &lt;- chosen
 *     N=512  ( 114,975 ticks): worst nonzero-quanta peak fraction 0.9306  (FAILS 0.95)
 *     N=1024 ( 230,175 ticks): worst nonzero-quanta peak fraction 0.9208  (FAILS 0.95)
 *     N=2048 ( 460,575 ticks): worst nonzero-quanta peak fraction 0.4228  (FAILS badly)
 * </pre>
 *
 * Bin alignment itself (peakBin == expectedBin) held exactly at every N
 * swept - only peak-bin CONCENTRATION degraded with tick count, confirming
 * this was chirp-style leakage from accumulation rounding, not a
 * bin-arithmetic error. Concentration was NOT monotonic in {@code
 * |quanta|}: at the default parameters the worst case was quanta=-4
 * (0.9914), not the range extremes quanta=&plusmn;5 - swept empirically,
 * not guessed. This non-monotonicity was itself a symptom of the
 * pre-wrap rounding-chirp mechanism and should not be assumed to still
 * hold post-wrap.
 *
 * <h2>Nyquist / aliasing bound on seeded quanta</h2>
 * The bin-alignment math (previous section) maps quanta to a digital
 * frequency of {@code quanta*STRIDE/PHASE_RESOLUTION} cycles/sample, which
 * is only meaningful modulo 1 cycle - so two quanta values {@code
 * PHASE_RESOLUTION/STRIDE} apart are indistinguishable (alias to the same
 * bin), and the usual Nyquist folding halves that: {@code minQuanta} and
 * {@code maxQuanta} must satisfy {@code |quanta| <
 * PHASE_RESOLUTION/(2*STRIDE)}. At the default {@code STRIDE=225} that
 * bound is {@code 3600/450 = 8}, comfortably above the default {@code
 * [-5,5]} range. Concretely verified: quanta=10 and quanta=-6 (16 apart -
 * exactly one full aliasing period) both land on bin 160 of a 256-sample
 * spectrum at stride 225 - an aliased quanta value and its wrap-equivalent
 * are numerically indistinguishable, so {@code
 * k0LineFrequencyEqualsQuantaTimesQuantumRate} would keep passing even
 * though {@code expectedBin} was computed from the wrong provenance
 * (whichever of the two quanta values the caller actually meant). {@link
 * #run(Point3i, long, int, int, int, int)} enforces this bound (and
 * {@code fftLength} power-of-two, {@code stride > 0}, {@code minQuanta <=
 * maxQuanta}) up front so a violation fails at the call site, not
 * mysteriously inside {@link Fft}.
 *
 * <p><b>Re-sweep before changing any parameter.</b> The 0.95 concentration
 * margin documented above is verified ONLY for the default parameters
 * (extent (2,2,2), seed 42L, quanta range [-5,5], stride 225, N=256).
 * Changing the seed, quanta range, extent, stride, or FFT length
 * invalidates that sweep - re-run it (see this class's sweep methodology
 * in the section above) before regenerating the golden artifact with new
 * parameters; do not assume the margin transfers.
 *
 * <h2>Provenance / golden-artifact convention</h2>
 * The committed TSV's {@code gitCommit} header field is always the literal
 * {@code UNCOMMITTED} - the commit that actually introduces the artifact
 * cannot be known before it is made (a chicken-and-egg problem), so
 * embedding any other value would just be wrong by the time it lands. Data
 * rows are per-member SUMMARIES (memberIndex, quanta, expectedBin,
 * peakBin, peakFraction, spectralEntropy) printed at fixed precision
 * ({@code %.9e}), not full spectra - keeps the artifact small and
 * diffable. {@link BaselineK0SpectrumTest#goldenArtifactMatchesRegeneration}
 * compares the printed forms numerically (tolerance {@value #TOLERANCE}),
 * never byte-for-byte, per this bead's plan-audit correction: strict JVM
 * floating point is not guaranteed identical across platforms/JIT for a
 * transcendental-heavy pipeline like this one.
 *
 * <p><b>Row-count reconciliation.</b> The golden artifact has 120 data
 * rows (one per active member: 4 even-parity cells * 30 members at the
 * default extent), but the pure-tone ({@code k0SpectrumIsPureTones}) and
 * line-frequency ({@code k0LineFrequencyEqualsQuantaTimesQuantumRate})
 * tests only check 110 of them - the 10 members whose seeded quanta is
 * exactly 0 are deliberately skipped by both tests (a non-rotating
 * member's spectrum is a trivial all-DC degenerate case, not "a single
 * dominant line" in the sense those tests are checking), NOT a bug in
 * either the artifact or the test filtering. Both tests assert a nonzero
 * checked-count precisely to guard against this filter silently reducing
 * to zero.
 *
 * @author halhildebrand
 */
public final class BaselineSpectrumHarness {

    /** memberIndex, quanta, expectedBin, peakBin, peakFraction, spectralEntropy */
    public static final class MemberSpectrum {
        public final int    memberIndex;
        public final float  quanta;
        public final int    expectedBin;
        public final int    peakBin;
        public final double peakFraction;
        public final double spectralEntropy;

        MemberSpectrum(int memberIndex, float quanta, int expectedBin,
                       int peakBin, double peakFraction,
                       double spectralEntropy) {
            this.memberIndex = memberIndex;
            this.quanta = quanta;
            this.expectedBin = expectedBin;
            this.peakBin = peakBin;
            this.peakFraction = peakFraction;
            this.spectralEntropy = spectralEntropy;
        }
    }

    public static final class Result {
        public final Point3i             extent;
        public final long                seed;
        public final int                 fftLength;
        public final int                 stride;
        public final int                 minQuanta;
        public final int                 maxQuanta;
        public final List<MemberSpectrum> members;
        public final float[]             frequencyBefore;
        public final float[]             frequencyAfter;
        public final boolean             frequencyFieldUnchanged;

        Result(Point3i extent, long seed, int fftLength, int stride,
               int minQuanta, int maxQuanta, List<MemberSpectrum> members,
               float[] frequencyBefore, float[] frequencyAfter) {
            this.extent = extent;
            this.seed = seed;
            this.fftLength = fftLength;
            this.stride = stride;
            this.minQuanta = minQuanta;
            this.maxQuanta = maxQuanta;
            this.members = members;
            this.frequencyBefore = frequencyBefore;
            this.frequencyAfter = frequencyAfter;
            this.frequencyFieldUnchanged = Arrays.equals(frequencyBefore,
                                                           frequencyAfter);
        }
    }

    public static final Point3i DEFAULT_EXTENT     = new Point3i(2, 2, 2);
    public static final long    DEFAULT_SEED        = 42L;
    /** See class javadoc "Float32 precision budget" for the empirical sweep behind this value. */
    public static final int     FFT_LENGTH          = 256;
    /** See class javadoc "The bin-alignment problem and its fix". */
    public static final int     STRIDE              = 225;
    public static final int     MIN_QUANTA          = -5;
    public static final int     MAX_QUANTA          = 5;

    static final double         TOLERANCE           = 1e-6;

    private static final String GOLDEN_RELATIVE_PATH = "src/test/resources/lga/baseline-k0-spectrum.tsv";

    private BaselineSpectrumHarness() {
    }

    /**
     * Regenerates the K=0 baseline with the default (documented) parameters
     * and overwrites the committed golden artifact. Run manually
     * (IDE/classpath invocation - no exec plugin is configured in this
     * project) whenever the baseline parameters intentionally change; the
     * regenerated file must then be reviewed and committed by hand.
     */
    public static void main(String[] args) throws IOException {
        Result result = run();
        String tsv = toTsv(result);
        Path path = Paths.get(GOLDEN_RELATIVE_PATH);
        Files.createDirectories(path.getParent());
        Files.write(path, tsv.getBytes(StandardCharsets.UTF_8));
        System.out.println("Wrote " + path.toAbsolutePath() + " ("
                            + result.members.size() + " member rows)");
    }

    /**
     * The upper bound on Shannon spectral entropy (nats) consistent with a
     * peak-bin power fraction of at least {@code minConcentration} over
     * {@code n} bins: the entropy of the worst case (maximally spread)
     * remainder, {@code -p*ln(p) - (1-p)*ln((1-p)/(n-1))}. For {@code p >
     * 0.5} this bound is monotonically DEcreasing in {@code p}, so the
     * bound evaluated at exactly {@code minConcentration} is a valid upper
     * bound for any measured distribution whose peak fraction is at least
     * {@code minConcentration} - a derived threshold, not an arbitrary
     * one.
     */
    public static double maxSpectralEntropyForConcentration(double minConcentration,
                                                              int n) {
        if (minConcentration <= 0.0 || minConcentration >= 1.0) {
            throw new IllegalArgumentException(
            "minConcentration must be in (0,1), was " + minConcentration);
        }
        if (n <= 1) {
            throw new IllegalArgumentException("n must be > 1, was " + n);
        }
        double rest = 1.0 - minConcentration;
        double perBin = rest / (n - 1);
        double h = -minConcentration * Math.log(minConcentration);
        h -= rest * Math.log(perBin);
        return h;
    }

    /**
     * A coarse-grained transport statistic for the quanta (frequency)
     * field: the ratio of its across-member variance after the run to
     * before. Guarded: if the field did not change at all (zero collision
     * events - exactly the K=0 case, since {@link
     * Necronomata#process(Point3i)} never writes {@code deltaF}), or if
     * the before-variance is zero, the statistic is UNDEFINED and this
     * returns {@link OptionalDouble#empty()} rather than a numeric ratio -
     * a naive computation on identical arrays would otherwise silently
     * return exactly 1.0, which reads as "measured, normal diffusion" when
     * in fact no transport process occurred at all to measure.
     */
    public static OptionalDouble quantaSpreadRatio(float[] frequencyBefore,
                                                     float[] frequencyAfter) {
        if (Arrays.equals(frequencyBefore, frequencyAfter)) {
            return OptionalDouble.empty();
        }
        double varianceBefore = variance(frequencyBefore);
        if (varianceBefore == 0.0) {
            return OptionalDouble.empty();
        }
        return OptionalDouble.of(variance(frequencyAfter) / varianceBefore);
    }

    /**
     * Runs the baseline with the default (documented) parameters: extent
     * {@link #DEFAULT_EXTENT}, seed {@link #DEFAULT_SEED}, {@link
     * #FFT_LENGTH}-sample spectra strided by {@link #STRIDE}, quanta
     * seeded uniformly in {@code [}{@link #MIN_QUANTA}{@code ,}{@link
     * #MAX_QUANTA}{@code ]}.
     */
    public static Result run() {
        return run(DEFAULT_EXTENT, DEFAULT_SEED, FFT_LENGTH, STRIDE,
                    MIN_QUANTA, MAX_QUANTA);
    }

    public static Result run(Point3i extent, long seed, int fftLength,
                              int stride, int minQuanta, int maxQuanta) {
        validateParameters(fftLength, stride, minQuanta, maxQuanta);
        Necronomata automaton = new Necronomata(extent);
        int[] memberIndices = collectActiveMemberIndices(automaton);
        float[] quanta = seedFrequencies(automaton, memberIndices, seed,
                                          minQuanta, maxQuanta);
        float[] frequencyBefore = snapshotFrequency(automaton);

        float[][] series = recordStridedAngleSeries(automaton,
                                                      memberIndices,
                                                      fftLength, stride);

        float[] frequencyAfter = snapshotFrequency(automaton);

        List<MemberSpectrum> summaries = new ArrayList<>(memberIndices.length);
        for (int m = 0; m < memberIndices.length; m++) {
            double[] power = SpectrumAnalyzer.powerSpectrum(series[m],
                                                              WindowFunction.RECTANGULAR);
            int peakBin = SpectrumAnalyzer.peakBin(power);
            double total = sum(power);
            double peakFraction = power[peakBin] / total;
            double entropy = spectralEntropy(power, total);
            int expectedBin = SpectrumAnalyzer.expectedBinForFrequency(quanta[m]
                                                                        * stride,
                                                                        fftLength);
            summaries.add(new MemberSpectrum(memberIndices[m], quanta[m],
                                              expectedBin, peakBin,
                                              peakFraction, entropy));
        }

        return new Result(extent, seed, fftLength, stride, minQuanta,
                           maxQuanta, summaries, frequencyBefore,
                           frequencyAfter);
    }

    /** Data rows (excludes provenance header) as printed for the TSV: [memberIndex, quanta, expectedBin, peakBin, peakFraction, spectralEntropy]. */
    public static List<String[]> toDataRows(Result result) {
        List<String[]> rows = new ArrayList<>(result.members.size());
        for (MemberSpectrum m : result.members) {
            rows.add(new String[] { Integer.toString(m.memberIndex),
                                     Integer.toString((int) m.quanta),
                                     Integer.toString(m.expectedBin),
                                     Integer.toString(m.peakBin),
                                     formatPrecise(m.peakFraction),
                                     formatPrecise(m.spectralEntropy) });
        }
        return rows;
    }

    public static String toTsv(Result result) {
        StringBuilder sb = new StringBuilder();
        sb.append("# BaselineSpectrumHarness golden artifact - K=0 (collision-free) baseline spectrum\n");
        sb.append("# bead=inviscid-0nx.7\n");
        sb.append("# extent=").append(result.extent.x).append(',')
          .append(result.extent.y).append(',').append(result.extent.z)
          .append('\n');
        sb.append("# seed=").append(result.seed).append('\n');
        sb.append("# fftLength=").append(result.fftLength).append('\n');
        sb.append("# stride=").append(result.stride).append('\n');
        sb.append("# quantumRate=").append(Necronomata.QUANTUM_RATE)
          .append('\n');
        sb.append("# phaseResolution=").append(Necronomata.PHASE_RESOLUTION)
          .append('\n');
        sb.append("# frequencyDistribution=uniform integer in [")
          .append(result.minQuanta).append(',').append(result.maxQuanta)
          .append("] via new Random(seed).nextInt(").append(result.maxQuanta
                                                              - result.minQuanta
                                                              + 1)
          .append(") + ").append(result.minQuanta).append('\n');
        sb.append("# window=RECTANGULAR\n");
        sb.append("# collisionRule=NONE (Necronomata.process(Point3i) is a no-op - K=0 baseline, predates inviscid-0nx.14)\n");
        sb.append("# generator=com.chiralbehaviors.inviscid.measure.BaselineSpectrumHarness\n");
        sb.append("# gitCommit=UNCOMMITTED (see class javadoc \"Provenance / golden-artifact convention\")\n");
        sb.append("# comparisonTolerance=").append(TOLERANCE)
          .append(" (numeric, not byte-exact - see BaselineK0SpectrumTest#goldenArtifactMatchesRegeneration)\n");
        sb.append("# precision=%.9e\n");
        sb.append("# columns=memberIndex\tquanta\texpectedBin\tpeakBin\tpeakFraction\tspectralEntropy\n");
        for (String[] row : toDataRows(result)) {
            sb.append(String.join("\t", row)).append('\n');
        }
        return sb.toString();
    }

    private static int[] collectActiveMemberIndices(Necronomata automaton) {
        List<Integer> indices = new ArrayList<>();
        automaton.forEach(cell -> {
            int base = automaton.indexOfCell(cell);
            for (int local = 0; local < 30; local++) {
                indices.add(base + local);
            }
        });
        int[] result = new int[indices.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = indices.get(i);
        }
        return result;
    }

    private static String formatPrecise(double v) {
        return String.format(Locale.ROOT, "%.9e", v);
    }

    /**
     * Records {@code samples} angle values for every member in {@code
     * memberIndices} in a single synchronized pass - one {@link
     * Necronomata#step()} call advances every member at once regardless of
     * how many are being tracked, so recording all members together (as
     * opposed to replaying the run once per member) keeps the total step
     * count independent of member count. Sample {@code s} is captured
     * BEFORE {@code stride} more ticks run (sample 0 is the pre-run
     * state), matching {@link SpectrumAnalyzer#recordAngleSeries}'s
     * "first recorded value is before any ticks" contract.
     */
    private static float[][] recordStridedAngleSeries(Necronomata automaton,
                                                        int[] memberIndices,
                                                        int samples,
                                                        int stride) {
        float[][] series = new float[memberIndices.length][samples];
        float[][] angleBox = new float[1][];
        for (int s = 0; s < samples; s++) {
            automaton.process((angle, frequency, deltaA,
                                deltaF) -> angleBox[0] = angle);
            float[] angle = angleBox[0];
            for (int m = 0; m < memberIndices.length; m++) {
                series[m][s] = angle[memberIndices[m]];
            }
            if (s < samples - 1) {
                for (int t = 0; t < stride; t++) {
                    automaton.step();
                }
            }
        }
        return series;
    }

    private static float[] seedFrequencies(Necronomata automaton,
                                            int[] memberIndices, long seed,
                                            int minQuanta, int maxQuanta) {
        Random rnd = new Random(seed);
        int range = maxQuanta - minQuanta + 1;
        float[] quanta = new float[memberIndices.length];
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int m = 0; m < memberIndices.length; m++) {
                int q = rnd.nextInt(range) + minQuanta;
                frequency[memberIndices[m]] = q;
                quanta[m] = q;
            }
        });
        return quanta;
    }

    private static float[] snapshotFrequency(Necronomata automaton) {
        float[][] box = new float[1][];
        automaton.process((angle, frequency, deltaA,
                            deltaF) -> box[0] = frequency.clone());
        return box[0];
    }

    private static double spectralEntropy(double[] power, double total) {
        double h = 0.0;
        for (double p : power) {
            if (p <= 0.0) {
                continue;
            }
            double pi = p / total;
            h -= pi * Math.log(pi);
        }
        return h;
    }

    /**
     * Fails at the call site (near the mistake) rather than inside {@link
     * Fft} or as a mysteriously-low concentration/wrong-provenance result.
     * See the class javadoc's "Nyquist / aliasing bound on seeded quanta"
     * section for why the {@code |quanta| < PHASE_RESOLUTION/(2*stride)}
     * bound is required.
     */
    private static void validateParameters(int fftLength, int stride,
                                             int minQuanta, int maxQuanta) {
        if (fftLength <= 0 || (fftLength & (fftLength - 1)) != 0) {
            throw new IllegalArgumentException(
            "fftLength must be a power of two, was " + fftLength);
        }
        if (stride <= 0) {
            throw new IllegalArgumentException("stride must be positive, was "
                                                + stride);
        }
        if (minQuanta > maxQuanta) {
            throw new IllegalArgumentException("minQuanta (" + minQuanta
                                                + ") must be <= maxQuanta ("
                                                + maxQuanta + ")");
        }
        double bound = Necronomata.PHASE_RESOLUTION / (2.0 * stride);
        if (Math.abs((double) minQuanta) >= bound
            || Math.abs((double) maxQuanta) >= bound) {
            throw new IllegalArgumentException(
            "seeded quanta range [" + minQuanta + "," + maxQuanta
            + "] violates the Nyquist/aliasing bound for stride=" + stride
            + ": |quanta| must be strictly less than PHASE_RESOLUTION/(2*stride) = "
            + bound
            + " (quanta at or beyond this bound alias onto another quanta's"
            + " bin - e.g. quanta=10 and quanta=-6 both land on bin 160 at"
            + " stride=225, N=256 - making the line-frequency test pass with"
            + " the wrong provenance)");
        }
    }

    private static double sum(double[] xs) {
        double total = 0.0;
        for (double x : xs) {
            total += x;
        }
        return total;
    }

    private static double variance(float[] xs) {
        double mean = 0.0;
        for (float x : xs) {
            mean += x;
        }
        mean /= xs.length;
        double variance = 0.0;
        for (float x : xs) {
            double d = x - mean;
            variance += d * d;
        }
        return variance / xs.length;
    }
}
