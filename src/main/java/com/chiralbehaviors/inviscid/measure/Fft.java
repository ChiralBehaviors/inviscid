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

/**
 * In-place radix-2 decimation-in-time complex FFT. No external dependency
 * (Commons Math etc deliberately not pulled in for ~60 lines of textbook
 * Cooley-Tukey) - correctness instead rests on {@code FftTest}'s
 * independent O(N^2) direct-summation DFT oracle, not on this
 * implementation being trusted on its own output.
 *
 * <p>{@code re}/{@code im} are mutated in place; length must be a power of
 * two (including 1, the trivial no-op case) or {@link #fft(double[],
 * double[])} throws.
 *
 * @author halhildebrand
 */
public final class Fft {

    private Fft() {
    }

    /**
     * Forward transform: {@code X[k] = sum_n x[n] * exp(-2*pi*i*k*n/N)}, in
     * place.
     *
     * @throws IllegalArgumentException if {@code re.length != im.length} or
     *             the length is not a power of two.
     */
    public static void fft(double[] re, double[] im) {
        transform(re, im, -1.0);
    }

    /**
     * Inverse transform: {@code x[n] = (1/N) * sum_k X[k] *
     * exp(+2*pi*i*k*n/N)}, in place. Forward then inverse recovers the
     * original signal (to floating-point precision).
     *
     * @throws IllegalArgumentException if {@code re.length != im.length} or
     *             the length is not a power of two.
     */
    public static void ifft(double[] re, double[] im) {
        transform(re, im, 1.0);
        int n = re.length;
        for (int i = 0; i < n; i++) {
            re[i] /= n;
            im[i] /= n;
        }
    }

    /**
     * Bit-reversal permutation, the standard iterative-FFT precondition:
     * after this, element {@code i} holds what will become the butterfly
     * input at the bit-reversed index of {@code i}.
     */
    private static void bitReverse(double[] re, double[] im) {
        int n = re.length;
        for (int i = 1, j = 0; i < n; i++) {
            int bit = n >> 1;
            for (; (j & bit) != 0; bit >>= 1) {
                j ^= bit;
            }
            j ^= bit;
            if (i < j) {
                double tr = re[i];
                re[i] = re[j];
                re[j] = tr;
                double ti = im[i];
                im[i] = im[j];
                im[j] = ti;
            }
        }
    }

    /**
     * Shared forward/inverse butterfly network. {@code sign} is -1 for the
     * forward transform's {@code exp(-i*theta)} twiddle, +1 for the
     * inverse's {@code exp(+i*theta)}; the inverse's {@code 1/N} scaling is
     * applied by the caller, not here.
     */
    private static void transform(double[] re, double[] im, double sign) {
        if (re.length != im.length) {
            throw new IllegalArgumentException(
            "re and im must be the same length: " + re.length + " != "
            + im.length);
        }
        int n = re.length;
        if (n == 0 || (n & (n - 1)) != 0) {
            throw new IllegalArgumentException(
            "length must be a power of two, was " + n);
        }

        bitReverse(re, im);

        for (int size = 2; size <= n; size <<= 1) {
            int half = size / 2;
            double theta = sign * 2 * Math.PI / size;
            double wRe = Math.cos(theta);
            double wIm = Math.sin(theta);
            for (int start = 0; start < n; start += size) {
                double curRe = 1.0;
                double curIm = 0.0;
                for (int j = 0; j < half; j++) {
                    int evenIdx = start + j;
                    int oddIdx = evenIdx + half;
                    double tRe = re[oddIdx] * curRe - im[oddIdx] * curIm;
                    double tIm = re[oddIdx] * curIm + im[oddIdx] * curRe;
                    re[oddIdx] = re[evenIdx] - tRe;
                    im[oddIdx] = im[evenIdx] - tIm;
                    re[evenIdx] += tRe;
                    im[evenIdx] += tIm;
                    double nextRe = curRe * wRe - curIm * wIm;
                    double nextIm = curRe * wIm + curIm * wRe;
                    curRe = nextRe;
                    curIm = nextIm;
                }
            }
        }
    }
}
