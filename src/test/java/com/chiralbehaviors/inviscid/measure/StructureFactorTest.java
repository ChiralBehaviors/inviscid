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

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.List;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.Necronomata;

/**
 * B.4 (bead inviscid-0nx.9): the dynamic structure factor S(k,omega)
 * instrument that turns "collective branches" from a hope into a
 * measurement.
 *
 * @author halhildebrand
 */
public class StructureFactorTest {

    /**
     * THE DEFINITIVE GUARD (stacked-review Critical 1). {@code
     * DispersionPoint.k} must report the Euclidean magnitude
     * {@code |k_vec|}, not the per-axis reciprocal-space component: for
     * X110/X111 those differ by {@code sqrt(2)}/{@code sqrt(3)}. A
     * per-axis-component bug is otherwise invisible - each single-direction
     * oracle test (X100/X110/X111) independently self-consistently checks
     * its OWN {@code k}, so nothing catches the cross-direction magnitude
     * error except a test that measures the SAME physical dispersion
     * relation along more than one direction and compares the results.
     *
     * <p>Constructs three INDEPENDENT synthetic fields, one per direction,
     * each obeying the SAME isotropic dispersion relation {@code omega =
     * c*|k_vec|} (a physically trivial case - sound in a uniform medium,
     * say) for several {@code m} along that direction; fits a ridge to
     * each; and asserts the three fitted slopes (each an estimate of the
     * same {@code c}) agree with each other well within the tolerance that
     * separates them from a {@code sqrt(2)}/{@code sqrt(3)} per-axis-k bug
     * (a ~41%/~73% deviation) - this test is EXPECTED, and was verified
     * (bead inviscid-0nx.9 stacked review), to FAIL against the
     * per-axis-component code with exactly that magnitude of cross-slope
     * disagreement.
     */
    @Test
    public void isotropicDispersionGivesConsistentSlopeAcrossDirections() {
        Point3i extent = new Point3i(8, 8, 8);
        StructureFactor sf = new StructureFactor(extent);
        int T = 256;
        double c = 0.25;
        int[] ms = { 1, 2, 3 };

        double slope100 = fittedIsotropicSlope(sf, extent, T, c,
                                                StructureFactor.Direction.X100,
                                                1, ms);
        double slope110 = fittedIsotropicSlope(sf, extent, T, c,
                                                StructureFactor.Direction.X110,
                                                2, ms);
        double slope111 = fittedIsotropicSlope(sf, extent, T, c,
                                                StructureFactor.Direction.X111,
                                                3, ms);

        double tolerance = 0.15 * c;
        assertEquals("X100 vs X110 fitted slope for the SAME isotropic omega=c|k| relation "
                    + "must agree within " + tolerance + " - a per-axis-k bug would show "
                    + "a spurious sqrt(2)x (~41%) discrepancy here",
                    slope100, slope110, tolerance);
        assertEquals("X100 vs X111 fitted slope for the SAME isotropic omega=c|k| relation "
                    + "must agree within " + tolerance + " - a per-axis-k bug would show "
                    + "a spurious sqrt(3)x (~73%) discrepancy here",
                    slope100, slope111, tolerance);
        assertEquals("X110 vs X111 fitted slope for the SAME isotropic omega=c|k| relation "
                    + "must agree within " + tolerance,
                    slope110, slope111, tolerance);
    }

    private static double fittedIsotropicSlope(StructureFactor sf,
                                                 Point3i extent, int T,
                                                 double c,
                                                 StructureFactor.Direction direction,
                                                 int nd, int[] ms) {
        int n = extent.x * extent.y * extent.z;
        double[][] re = new double[T][n];
        double[][] im = new double[T][n];
        double sqrtNd = Math.sqrt(nd);
        for (int t = 0; t < T; t++) {
            for (int i = 0; i < extent.x; i++) {
                for (int j = 0; j < extent.y; j++) {
                    for (int k = 0; k < extent.z; k++) {
                        if (((i + j + k) & 1) != 0) {
                            continue;
                        }
                        int idx = (i * extent.y + j) * extent.z + k;
                        double sumOverM = 0;
                        double sumImOverM = 0;
                        for (int m : ms) {
                            double kComponent = 2 * Math.PI * m / extent.x;
                            double kMag = kComponent * sqrtNd;
                            double omega = c * kMag;
                            double dotX;
                            switch (direction) {
                            case X100:
                                dotX = kComponent * i;
                                break;
                            case X110:
                                dotX = kComponent * (i + j);
                                break;
                            case X111:
                                dotX = kComponent * (i + j + k);
                                break;
                            default:
                                throw new IllegalStateException();
                            }
                            double theta = dotX - omega * t;
                            sumOverM += Math.cos(theta);
                            sumImOverM += Math.sin(theta);
                        }
                        re[t][idx] = sumOverM;
                        im[t][idx] = sumImOverM;
                    }
                }
            }
        }
        List<StructureFactor.DispersionPoint> points = sf.spectrum(re, im,
                                                                     direction);
        List<StructureFactor.DispersionPoint> onRidge = new java.util.ArrayList<>();
        for (int m : ms) {
            onRidge.add(points.get(m));
        }
        return sf.extractRidge(onRidge).slope();
    }

    /**
     * THE MIRROR-REINFORCEMENT GUARD (stacked-review round 2, critic).
     * {@code DispersionPoint.k} must be SIGNED, not folded to a
     * non-negative magnitude via {@code Math.abs()}: a real-valued
     * field's spatial DFT is Hermitian - {@code F(-k,t) ==
     * conj(F(k,t))} for any real {@code f}, exactly, since flipping the
     * sign of every {@code m} negates {@code theta} in
     * {@link StructureFactor#spatialDft} (re unchanged, im flips) - which
     * propagates through the temporal transform to {@code power(-k,omega)
     * == power(k,-omega)}. The "mirror" index ({@code axisExtent-m}) of
     * ANY real single-direction propagating wave therefore lands at the
     * SAME {@code |k|} with the OPPOSITE-sign {@code omega}. An
     * {@code abs()}-folded {@code k} makes that mirror collide with the
     * forward point at an IDENTICAL reported {@code k} - and an
     * unweighted OLS fit over the raw, UNFILTERED {@code spectrum(...)}
     * output (exactly what {@link #extractRidge} does, and exactly what a
     * real caller - bead .10 - does) sees two points at the same k with
     * opposite-sign omega CANCEL rather than reinforce. With a genuinely
     * signed k, the mirror instead lands at {@code -k} with {@code
     * -omega} - exactly on the SAME dispersion line - and reinforces the
     * fit instead.
     *
     * <p>Deliberately uses the PUBLIC real-only
     * {@code spectrum(fieldByTick, direction)} entry point (not the
     * complex-oracle overload the other tests use) over its full,
     * un-restricted {@code m} range, and fits {@link #extractRidge}
     * directly on that UNFILTERED list - no manual pre-selection of a
     * "safe" sub-range, which is exactly the trap a caller unaware of the
     * Hermitian-mirror issue would fall into.
     *
     * <h2>Why the field also injects a DC term and a time-independent
     * Nyquist-patterned term</h2>
     * A field built from ONLY {@code ms={1,2,3}} leaves {@code m=0} and
     * {@code m=extent/2} (Nyquist) with near-zero total power - and an
     * unweighted argmax over near-machine-epsilon noise picks an
     * essentially ARBITRARY bin there (verified empirically during
     * development: peak fractions around 0.1 and wildly wrong omegas, not
     * the naively-expected clean {@code (k,0)}). Those garbage points then
     * get EQUAL WEIGHT to the genuine high-power points in the unweighted
     * OLS, corrupting the fit for reasons unrelated to the sign-convention
     * bug under test. Adding a constant (DC, always {@code omega=0}, no
     * ambiguity) grounds {@code m=0} robustly. The Nyquist point
     * ({@code m=extent/2}) is trickier: it is self-mirrored (its own
     * mirror index), so an OSCILLATING real injection there splits its
     * power EXACTLY 50/50 between {@code +omega} and {@code -omega} with a
     * genuinely ambiguous argmax tie-break (verified to disagree between
     * two independently-correct DFT implementations - not a bug in
     * either, just a real tie). This test instead injects a
     * TIME-INDEPENDENT term there ({@code cos(nyquistComponent . x)}, no
     * {@code t} dependence) - unambiguously {@code omega=0}, deterministic
     * across platforms, but consequently a KNOWN, fixed point OFF the true
     * {@code omega=c*k} line (at {@code (extent/2-component, 0)} rather
     * than the line's {@code (extent/2-component, c*pi)}).
     *
     * <p>That one deterministic off-line point biases the fit down from a
     * perfect match: empirically verified (bead inviscid-0nx.9 stacked
     * review round 2, both directions) - {@code Math.abs()}-folded k gives
     * slope EXACTLY {@code 0.0} (total cancellation: the six genuine
     * mirror-paired points cancel each other exactly, leaving only the
     * already-zero DC and Nyquist contributions); signed k gives slope
     * {@code ~0.667*c} (the six mirror-paired points reinforce correctly;
     * the fit is pulled down from a perfect {@code 1.0*c} only by the one
     * deliberately-off-line Nyquist point). The two outcomes - EXACT zero
     * vs a robust, substantial {@code ~2/3} of {@code c} - are separated
     * by far more than either implementation's numerical noise.
     */
    @Test
    public void realFieldRidgeReinforcesAcrossMirrorPointsX100() {
        Point3i extent = new Point3i(8, 8, 8);
        StructureFactor sf = new StructureFactor(extent);
        double c = 0.25;
        double slope = realOneDirectionalRidgeSlope(sf, extent, 32, c,
                                                      StructureFactor.Direction.X100,
                                                      1, new int[] { 1, 2, 3 });
        assertTrue("real-field, full-UNFILTERED-range ridge slope must be a substantial "
                   + "fraction of the true dispersion constant c=" + c + ", was " + slope
                   + " - a value at or near zero (the abs()-folded-k signature: the six "
                   + "mirror-paired points cancel EXACTLY to slope=0.0) means the mirror "
                   + "points at (-k,-omega) are CANCELING instead of reinforcing the fit",
                   slope > 0.5 * c);
        assertTrue("fitted slope must not exceed the true dispersion constant by more "
                   + "than a modest margin, was " + slope, slope < 0.85 * c);
    }

    /**
     * Same guard, X110 (the coordinator's explicit "and X110" ask) - the
     * mirror-collision mechanism is direction-agnostic (it is a property
     * of {@code k}-folding, not of how many components a probe direction
     * has), so the same reinforcement-not-cancellation property must hold
     * there too. Uses a larger {@code T} (256 vs X100's 32) since X110's
     * {@code sqrt(2)} magnitude factor is irrational, so its injected
     * {@code omega} values (other than the DC/Nyquist terms) are not
     * exactly bin-aligned; a larger {@code T} keeps the resulting
     * quantization error small.
     */
    @Test
    public void realFieldRidgeReinforcesAcrossMirrorPointsX110() {
        Point3i extent = new Point3i(8, 8, 8);
        StructureFactor sf = new StructureFactor(extent);
        double c = 0.25;
        double slope = realOneDirectionalRidgeSlope(sf, extent, 256, c,
                                                      StructureFactor.Direction.X110,
                                                      2, new int[] { 1, 2, 3 });
        assertTrue("X110 real-field, full-UNFILTERED-range ridge slope must be a "
                   + "substantial fraction of c=" + c + ", was " + slope + " - near zero "
                   + "means mirror cancellation (the abs()-folded-k bug)", slope > 0.5 * c);
        assertTrue("fitted slope must not exceed the true dispersion constant by more "
                   + "than a modest margin, was " + slope, slope < 0.85 * c);
    }

    private static double realOneDirectionalRidgeSlope(StructureFactor sf,
                                                         Point3i extent, int T,
                                                         double c,
                                                         StructureFactor.Direction direction,
                                                         int nd, int[] ms) {
        int n = extent.x * extent.y * extent.z;
        double[][] field = new double[T][n];
        double sqrtNd = Math.sqrt(nd);
        int nyquistM = extent.x / 2;
        double nyquistKComponent = 2 * Math.PI * nyquistM / extent.x;
        // DC term (m=0) grounds the otherwise near-zero-power/argmax-noise
        // point at a robust, unambiguous omega=0 - see the X100 test's
        // javadoc. The self-mirrored Nyquist point (m=extent.x/2) gets a
        // TIME-INDEPENDENT spatial term (omega=0 by construction, not by
        // an oscillation that would split its power 50/50 between +omega
        // and -omega with a genuinely ambiguous, floating-point-path-
        // dependent argmax tie-break - verified to disagree between two
        // independently correct DFT implementations when tried). This
        // deterministically pins the Nyquist point OFF the true dispersion
        // line at a KNOWN, reproducible (extent.x/2-component, 0) rather
        // than leaving it either noisy or platform-dependent.
        for (int t = 0; t < T; t++) {
            for (int i = 0; i < extent.x; i++) {
                for (int j = 0; j < extent.y; j++) {
                    for (int k = 0; k < extent.z; k++) {
                        if (((i + j + k) & 1) != 0) {
                            continue;
                        }
                        int idx = (i * extent.y + j) * extent.z + k;
                        double nyquistDotX;
                        switch (direction) {
                        case X100:
                            nyquistDotX = nyquistKComponent * i;
                            break;
                        case X110:
                            nyquistDotX = nyquistKComponent * (i + j);
                            break;
                        case X111:
                            nyquistDotX = nyquistKComponent * (i + j + k);
                            break;
                        default:
                            throw new IllegalStateException();
                        }
                        double sum = 1.0 + Math.cos(nyquistDotX);
                        for (int m : ms) {
                            double kComponent = 2 * Math.PI * m / extent.x;
                            double omega = c * kComponent * sqrtNd;
                            double dotX;
                            switch (direction) {
                            case X100:
                                dotX = kComponent * i;
                                break;
                            case X110:
                                dotX = kComponent * (i + j);
                                break;
                            case X111:
                                dotX = kComponent * (i + j + k);
                                break;
                            default:
                                throw new IllegalStateException();
                            }
                            sum += Math.cos(dotX - omega * t);
                        }
                        field[t][idx] = sum;
                    }
                }
            }
        }
        // REAL field: the public real-only entry point, not the
        // complex-oracle overload - and the FULL unfiltered points list,
        // deliberately not pre-selected to ms.
        List<StructureFactor.DispersionPoint> points = sf.spectrum(field,
                                                                     direction);
        return sf.extractRidge(points).slope();
    }

    /**
     * THE ORACLE. A synthetic complex plane wave sampled on the
     * even-parity sublattice must produce a single, essentially
     * uncontaminated (k0, omega0) spectral peak - the instrument
     * validated independently of the automaton.
     *
     * <p>Uses the complex analytic form {@code exp(i*(k0.x - omega0*t))},
     * not the bead-description's literal {@code cos(...)}: this mirrors
     * {@code SpectrumAnalyzer}'s already-locked precedent (see that
     * class's "ramp-vs-sinusoid" javadoc section) - a real cos/sin series
     * splits power 50/50 across a mirror-image (-k0,-omega0) bin, which
     * cannot clear a "single dominant peak, &gt;90% concentration" bar
     * except in a degenerate case. The physical production field
     * ({@link StructureFactor#coarseGrainedField(Necronomata)}) is always
     * real; this complex-field capability is exposed only via the
     * package-private {@code spectrum(re, im, direction)} overload used
     * here and in {@link #evenParitySublatticeIndexingIsCorrect()}.
     */
    @Test
    public void planeWaveGivesSingleKOmegaPeak() {
        Point3i extent = new Point3i(4, 4, 4);
        StructureFactor sf = new StructureFactor(extent);
        int T = 8;
        int m0 = 1;
        int p0 = 1;
        double k0 = 2 * Math.PI * m0 / extent.x;
        double omega0 = 2 * Math.PI * p0 / T;
        int n = extent.x * extent.y * extent.z;

        double[][] re = new double[T][n];
        double[][] im = new double[T][n];
        for (int t = 0; t < T; t++) {
            for (int i = 0; i < extent.x; i++) {
                for (int j = 0; j < extent.y; j++) {
                    for (int k = 0; k < extent.z; k++) {
                        if (((i + j + k) & 1) != 0) {
                            continue;
                        }
                        int idx = (i * extent.y + j) * extent.z + k;
                        double theta = k0 * i - omega0 * t;
                        re[t][idx] = Math.cos(theta);
                        im[t][idx] = Math.sin(theta);
                    }
                }
            }
        }

        List<StructureFactor.DispersionPoint> points = sf.spectrum(re, im,
                                                                     StructureFactor.Direction.X100);

        // Grand total power across the WHOLE (m, omega) grid, not just the
        // per-k peak, so "single peak" is falsifiable against leakage into
        // OTHER k's too, not merely against other omega bins at k0.
        double grandTotal = 0;
        for (int m = 0; m < extent.x; m++) {
            double[] power = sf.powerSpectrumAt(re, im, m, 0, 0);
            for (double p : power) {
                grandTotal += p;
            }
        }

        StructureFactor.DispersionPoint hit = points.get(m0);
        assertEquals(m0, hit.m());
        assertEquals(k0, hit.k(), 1e-9);
        assertEquals(omega0, hit.omega(), 1e-9);
        assertTrue("grand total power must be positive", grandTotal > 0);
        assertTrue("peak power fraction across the FULL (k,omega) grid must exceed 0.90, was "
                   + (hit.power() / grandTotal),
                   hit.power() / grandTotal >= 0.90);
        assertTrue("per-k peak fraction at k0 must exceed 0.90, was "
                   + hit.peakFraction(), hit.peakFraction() >= 0.90);
    }

    /**
     * THE ORACLE, X110 (stacked-review Critical 2 - X110 had zero
     * oracle coverage). Mirrors {@link #planeWaveGivesSingleKOmegaPeak()}
     * but along the [110] ray: {@code k.x == kComponent*(i+j)} for
     * {@code k = (kComponent, kComponent, 0)}. Directly exercises the
     * Critical-1 magnitude fix - the expected {@code hit.k()} is
     * {@code kComponent*sqrt(2)}, not the bare component.
     */
    @Test
    public void planeWaveAlongX110GivesSingleKOmegaPeak() {
        Point3i extent = new Point3i(4, 4, 4);
        StructureFactor sf = new StructureFactor(extent);
        int T = 8;
        int m0 = 1;
        int p0 = 1;
        double kComponent = 2 * Math.PI * m0 / extent.x;
        double expectedK = kComponent * Math.sqrt(2);
        double omega0 = 2 * Math.PI * p0 / T;
        int n = extent.x * extent.y * extent.z;

        double[][] re = new double[T][n];
        double[][] im = new double[T][n];
        for (int t = 0; t < T; t++) {
            for (int i = 0; i < extent.x; i++) {
                for (int j = 0; j < extent.y; j++) {
                    for (int k = 0; k < extent.z; k++) {
                        if (((i + j + k) & 1) != 0) {
                            continue;
                        }
                        int idx = (i * extent.y + j) * extent.z + k;
                        double theta = kComponent * (i + j) - omega0 * t;
                        re[t][idx] = Math.cos(theta);
                        im[t][idx] = Math.sin(theta);
                    }
                }
            }
        }

        List<StructureFactor.DispersionPoint> points = sf.spectrum(re, im,
                                                                     StructureFactor.Direction.X110);
        double grandTotal = 0;
        for (int m = 0; m < extent.x; m++) {
            double[] power = sf.powerSpectrumAt(re, im, m, m, 0);
            for (double p : power) {
                grandTotal += p;
            }
        }

        StructureFactor.DispersionPoint hit = points.get(m0);
        assertEquals(m0, hit.m());
        assertEquals(expectedK, hit.k(), 1e-9);
        assertEquals(omega0, hit.omega(), 1e-9);
        assertTrue("grand total power must be positive", grandTotal > 0);
        assertTrue("peak power fraction across the FULL (k,omega) grid must exceed 0.90, was "
                   + (hit.power() / grandTotal),
                   hit.power() / grandTotal >= 0.90);
        assertTrue("per-k peak fraction at k0 must exceed 0.90, was "
                   + hit.peakFraction(), hit.peakFraction() >= 0.90);
    }

    /**
     * The bead's Deliverable-named "full-shell option" (stacked-review
     * Important 3): the enumerated grid must contain no reciprocal-lattice
     * ghost-duplicate pairs (see class javadoc's {@code F(k+G0)==F(k)}
     * identity), and each probe direction's ray, restricted to the range
     * where its OWN tuple representation is the canonical one
     * ({@code m < extent.x/2} - the same range {@link
     * StructureFactor.Direction#X111} is always restricted to), must be a
     * subset of the shell's reported triples.
     */
    @Test
    public void fullShellExcludesGhostDuplicatesAndContainsProbeDirections() {
        Point3i extent = new Point3i(4, 4, 4);
        StructureFactor sf = new StructureFactor(extent);
        int T = 8;
        int n = extent.x * extent.y * extent.z;
        double[] snapshot = new double[n];
        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    if (((i + j + k) & 1) != 0) {
                        continue;
                    }
                    int idx = (i * extent.y + j) * extent.z + k;
                    snapshot[idx] = i + 2 * j + 3 * k + 1;
                }
            }
        }
        double[][] field = new double[T][];
        for (int t = 0; t < T; t++) {
            field[t] = snapshot;
        }

        List<StructureFactor.ShellPoint> shell = sf.fullShell(field, 10);
        assertTrue("shell must not be empty", !shell.isEmpty());

        java.util.Set<List<Integer>> reported = new java.util.HashSet<>();
        for (StructureFactor.ShellPoint p : shell) {
            reported.add(List.of(p.mx(), p.my(), p.mz()));
        }
        for (StructureFactor.ShellPoint p : shell) {
            int partnerMx = (p.mx() + extent.x / 2) % extent.x;
            int partnerMy = (p.my() + extent.y / 2) % extent.y;
            int partnerMz = (p.mz() + extent.z / 2) % extent.z;
            assertTrue("ghost partner of (" + p.mx() + "," + p.my() + ","
                       + p.mz() + ") must not also be reported",
                       !reported.contains(List.of(partnerMx, partnerMy,
                                                   partnerMz)));
        }

        for (int m = 0; m < extent.x / 2; m++) {
            assertTrue("X100 triple (" + m + ",0,0) must be in the shell",
                       reported.contains(List.of(m, 0, 0)));
            assertTrue("X110 triple (" + m + "," + m + ",0) must be in the shell",
                       reported.contains(List.of(m, m, 0)));
            assertTrue("X111 triple (" + m + "," + m + "," + m
                       + ") must be in the shell",
                       reported.contains(List.of(m, m, m)));
        }
    }

    /**
     * A time-independent (but spatially varying) field must show all its
     * power on the omega=0 ridge, for every k with any power at all.
     */
    @Test
    public void staticFieldGivesZeroFrequencyRidgeOnly() {
        Point3i extent = new Point3i(4, 4, 4);
        StructureFactor sf = new StructureFactor(extent);
        int T = 8;
        int n = extent.x * extent.y * extent.z;
        double[] snapshot = new double[n];
        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    if (((i + j + k) & 1) != 0) {
                        continue;
                    }
                    int idx = (i * extent.y + j) * extent.z + k;
                    snapshot[idx] = i + 2 * j + 3 * k + 1;
                }
            }
        }
        double[][] field = new double[T][];
        for (int t = 0; t < T; t++) {
            field[t] = snapshot;
        }

        List<StructureFactor.DispersionPoint> points = sf.spectrum(field,
                                                                     StructureFactor.Direction.X100);
        int checked = 0;
        for (StructureFactor.DispersionPoint p : points) {
            if (p.power() <= 0) {
                continue;
            }
            assertEquals("k index " + p.m()
                        + " should peak at omega=0 for a time-static field",
                        0.0, p.omega(), 1e-9);
            assertTrue("k index " + p.m()
                       + " peak fraction should be ~1.0 for a time-static field, was "
                       + p.peakFraction(), p.peakFraction() >= 0.999);
            checked++;
        }
        assertTrue("no nonzero-power k indices were checked - test is vacuous",
                   checked > 0);

        StructureFactor.Ridge ridge = sf.extractRidge(points);
        assertEquals("ridge slope should be ~0 for a purely static field",
                     0.0, ridge.slope(), 1e-6);
    }

    /**
     * The subtle one. The even-parity sublattice is not a simple cubic
     * grid - it is (in the conventional-cell description standard in
     * solid-state physics) an FCC Bravais lattice. See
     * {@link StructureFactor}'s class javadoc for the full
     * reciprocal-lattice derivation this test exercises directly:
     *
     * <p>(a) the [111] probe direction is the one direction whose
     * reciprocal-lattice ghost partner (a shift of {@code extent/2} on
     * ALL three axes simultaneously) maps the probed ray back onto
     * itself, so {@code spectrum(..., X111)} must report only the first
     * half of the m range;
     *
     * <p>(b) that restriction discards no independent information - the
     * excluded "ghost" half is PROVABLY (not just presumed) an exact
     * numeric duplicate of its retained partner, per the reciprocal-space
     * identity {@code F(k+G0) == F(k)} derived in the class javadoc;
     *
     * <p>(c) the non-diagonal [100] probe direction never self-intersects
     * that same ghost (its partner always leaves the probed ray), so it
     * reports the FULL, un-restricted m range with no contamination
     * either way.
     */
    @Test
    public void evenParitySublatticeIndexingIsCorrect() {
        Point3i extent = new Point3i(4, 4, 4);
        StructureFactor sf = new StructureFactor(extent);
        int T = 8;
        int m0 = 1;
        int p0 = 1;
        double k0 = 2 * Math.PI * m0 / extent.x;
        double omega0 = 2 * Math.PI * p0 / T;
        int n = extent.x * extent.y * extent.z;

        double[][] re = new double[T][n];
        double[][] im = new double[T][n];
        for (int t = 0; t < T; t++) {
            for (int i = 0; i < extent.x; i++) {
                for (int j = 0; j < extent.y; j++) {
                    for (int k = 0; k < extent.z; k++) {
                        if (((i + j + k) & 1) != 0) {
                            continue;
                        }
                        int idx = (i * extent.y + j) * extent.z + k;
                        // plane wave along [111]: k . x == k0 * (i + j + k)
                        double theta = k0 * (i + j + k) - omega0 * t;
                        re[t][idx] = Math.cos(theta);
                        im[t][idx] = Math.sin(theta);
                    }
                }
            }
        }

        List<StructureFactor.DispersionPoint> x111 = sf.spectrum(re, im,
                                                                   StructureFactor.Direction.X111);
        assertEquals("X111 probe must report only the first half of the m range "
                    + "(Brillouin-zone restriction, see class javadoc)",
                    extent.x / 2, x111.size());

        StructureFactor.DispersionPoint hit = x111.get(m0);
        // hit.k() reports the Euclidean |k_vec|, not the bare per-axis
        // component k0 used to build the phase above - X111 has 3 equal
        // nonzero components, so |k_vec| == k0*sqrt(3).
        assertEquals(k0 * Math.sqrt(3), hit.k(), 1e-9);
        assertEquals(omega0, hit.omega(), 1e-9);
        assertTrue("per-k peak fraction at k0 must exceed 0.90, was "
                   + hit.peakFraction(), hit.peakFraction() >= 0.90);

        // The excluded ghost partner is a PROVABLE exact duplicate.
        double[] retained = sf.powerSpectrumAt(re, im, m0, m0, m0);
        int ghostM = m0 + extent.x / 2;
        double[] ghost = sf.powerSpectrumAt(re, im, ghostM, ghostM, ghostM);
        assertArrayEquals("the discarded X111 ghost partner must be numerically "
                         + "identical to the retained entry it duplicates - proving "
                         + "it is not lost information, and not a distorted aliasing "
                         + "artifact (the option-(a) failure mode this class avoids)",
                           retained, ghost, 1e-9);

        // The [100] probe never self-intersects that ghost: full, un-restricted
        // m range, no contamination.
        List<StructureFactor.DispersionPoint> x100 = sf.spectrum(re, im,
                                                                   StructureFactor.Direction.X100);
        assertEquals(extent.x, x100.size());
    }

    /**
     * Constructing against an odd, or sub-4, extent axis must throw -
     * mirrors {@code FccNeighborhood}'s own precondition (reused here,
     * not re-implemented) rather than silently manufacturing spurious
     * anisotropy in the reported spectrum.
     */
    @Test
    public void rejectsNonPeriodicConfiguration() {
        try {
            new StructureFactor(new Point3i(3, 4, 4));
            fail("odd x extent must be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("even"));
        }
        try {
            new StructureFactor(new Point3i(2, 4, 4));
            fail("extent below the 4-per-axis floor must be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("4"));
        }
    }

    /**
     * X110 requires {@code extent.x == extent.y} (a diagonal probe needs a
     * matching discrete k grid on both axes); a legal (even, &gt;=4 per
     * axis) but MISMATCHED extent must still be rejected at probe time,
     * not silently miscomputed.
     */
    @Test
    public void rejectsMismatchedExtentForX110Probe() {
        Point3i extent = new Point3i(4, 6, 4);
        StructureFactor sf = new StructureFactor(extent);
        double[][] field = new double[8][extent.x * extent.y * extent.z];
        try {
            sf.spectrum(field, StructureFactor.Direction.X110);
            fail("X110 probe with extent.x != extent.y must be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("X110"));
        }
    }

    /**
     * X111 requires all three axes equal; a legal but unequal extent must
     * be rejected at probe time.
     */
    @Test
    public void rejectsMismatchedExtentForX111Probe() {
        Point3i extent = new Point3i(4, 4, 6);
        StructureFactor sf = new StructureFactor(extent);
        double[][] field = new double[8][extent.x * extent.y * extent.z];
        try {
            sf.spectrum(field, StructureFactor.Direction.X111);
            fail("X111 probe with unequal axes must be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(),
                       expected.getMessage().contains("X111"));
        }
    }

    /**
     * Coverage gap (stacked-review Important 4a+4b): every prior oracle
     * used a cubic, all-power-of-two-friendly extent. This exercises the
     * X100 probe - which never requires matching axes - on a NON-cubic
     * extent whose PROBED axis (extent.x = 6) is a legal FCC extent (even,
     * &gt;=4) that is NOT a power of two, directly exercising the class
     * javadoc's "direct DFT removes the power-of-two constraint" claim on
     * the axis that actually matters (the one the spatial DFT sums over).
     */
    @Test
    public void planeWaveOracleOnNonCubicNonPowerOfTwoExtent() {
        Point3i extent = new Point3i(6, 4, 8);
        StructureFactor sf = new StructureFactor(extent);
        int T = 8;
        int m0 = 1;
        int p0 = 1;
        double k0 = 2 * Math.PI * m0 / extent.x;
        double omega0 = 2 * Math.PI * p0 / T;
        int n = extent.x * extent.y * extent.z;

        double[][] re = new double[T][n];
        double[][] im = new double[T][n];
        for (int t = 0; t < T; t++) {
            for (int i = 0; i < extent.x; i++) {
                for (int j = 0; j < extent.y; j++) {
                    for (int k = 0; k < extent.z; k++) {
                        if (((i + j + k) & 1) != 0) {
                            continue;
                        }
                        int idx = (i * extent.y + j) * extent.z + k;
                        double theta = k0 * i - omega0 * t;
                        re[t][idx] = Math.cos(theta);
                        im[t][idx] = Math.sin(theta);
                    }
                }
            }
        }

        List<StructureFactor.DispersionPoint> points = sf.spectrum(re, im,
                                                                     StructureFactor.Direction.X100);
        assertEquals("X100 probe must report the full un-restricted extent.x range",
                     extent.x, points.size());

        double grandTotal = 0;
        for (int m = 0; m < extent.x; m++) {
            double[] power = sf.powerSpectrumAt(re, im, m, 0, 0);
            for (double p : power) {
                grandTotal += p;
            }
        }

        StructureFactor.DispersionPoint hit = points.get(m0);
        assertEquals(m0, hit.m());
        assertEquals(k0, hit.k(), 1e-9);
        assertEquals(omega0, hit.omega(), 1e-9);
        assertTrue("grand total power must be positive", grandTotal > 0);
        assertTrue("peak power fraction across the FULL (k,omega) grid must exceed 0.90, was "
                   + (hit.power() / grandTotal),
                   hit.power() / grandTotal >= 0.90);
    }

    /**
     * Non-vacuity guard (K=0 baseline, bead inviscid-0nx.7 / B.2's
     * pattern): a free-rotor field has a CONSTANT frequency (quanta)
     * field - {@code Necronomata.step()} never fires a collision, so
     * {@code deltaF} is always zero and the coarse-grained field is
     * bit-identical every tick. The instrument must report exactly that:
     * every k peaks at omega=0, and the extracted ridge has ~zero slope.
     * A propagating branch here would mean the instrument is finding
     * physics that cannot exist at K=0.
     */
    @Test
    public void k0BaselineProducesNoPropagatingRidge() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        // Seed exactly ONE member per even cell (position-dependent,
        // nonzero) rather than every member uniformly: coarseGrainedField
        // sums 30 members per cell, and 30 is a multiple of 5, so a
        // seed pattern like (globalIndex % 5) applied to all 30 members
        // sums to exactly zero every cell - a degenerate, vacuous field.
        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < extent.x; i++) {
                for (int j = 0; j < extent.y; j++) {
                    for (int k = 0; k < extent.z; k++) {
                        if (((i + j + k) & 1) != 0) {
                            continue;
                        }
                        int base = automaton.indexOfCell(i, j, k);
                        frequency[base] = i + 2 * j + 3 * k + 1;
                    }
                }
            }
        });

        StructureFactor sf = new StructureFactor(extent);
        int T = 8;
        double[][] field = new double[T][];
        for (int t = 0; t < T; t++) {
            field[t] = StructureFactor.coarseGrainedField(automaton);
            automaton.step();
        }

        List<StructureFactor.DispersionPoint> points = sf.spectrum(field,
                                                                     StructureFactor.Direction.X100);
        int checked = 0;
        for (StructureFactor.DispersionPoint p : points) {
            if (p.power() <= 0) {
                continue;
            }
            assertEquals("K=0 (collision-free) baseline: the coarse-grained quanta "
                        + "field is provably constant tick-to-tick (no deltaF ever "
                        + "fires), so every k must peak at omega=0",
                        0.0, p.omega(), 1e-9);
            assertTrue("k index " + p.m() + " peak fraction should be ~1.0, was "
                       + p.peakFraction(), p.peakFraction() >= 0.999);
            checked++;
        }
        assertTrue("no nonzero-power k indices were checked - test is vacuous",
                   checked > 0);

        StructureFactor.Ridge ridge = sf.extractRidge(points);
        assertEquals("no propagating branch: ridge slope must be ~0",
                     0.0, ridge.slope(), 1e-6);
    }

    /**
     * Positive-control bracket for {@link #k0BaselineProducesNoPropagatingRidge()}
     * (stacked-review Significant 5): that test only ever exercises a
     * FROZEN field (K=0, nothing changes), so the {@code step()}-between-
     * snapshots wiring itself (write {@code deltaF} via the sanctioned
     * collision-rule escape hatch -&gt; {@code step()} folds it into
     * {@code frequency} -&gt; next {@code coarseGrainedField} sees the new
     * value) was never actually exercised end to end. This test
     * hand-injects a REAL, alternating quanta transfer into a single
     * member every tick and asserts the instrument reports genuine
     * nonzero-omega structure - proof the wiring is live, not merely
     * proof the guard doesn't false-positive on a static input.
     *
     * <p>The injected sequence alternates the transferred member's
     * frequency between {@code +3} and {@code -3} (zero mean, period 2) so
     * ALL of its power lands on the Nyquist bin with nothing at DC -
     * deliberately avoiding a same-magnitude DC/Nyquist tie that a naive
     * 0/+6 alternation would hit (verified by hand during development: a
     * [0,6,0,6,...] sequence has EQUAL power at bin 0 and the Nyquist bin,
     * and the argmax tie-break silently resolves to bin 0 - the wrong,
     * vacuous answer for this test's purpose).
     */
    @Test
    public void dynamicQuantaTransferProducesNonzeroOmegaStructure() {
        Point3i extent = new Point3i(4, 4, 4);
        Necronomata automaton = new Necronomata(extent);
        int base = automaton.indexOfCell(0, 0, 0);
        automaton.process((angle, frequency, deltaA,
                            deltaF) -> frequency[base] = 3f);

        StructureFactor sf = new StructureFactor(extent);
        int T = 8;
        double[][] field = new double[T][];
        for (int t = 0; t < T; t++) {
            field[t] = StructureFactor.coarseGrainedField(automaton);
            int delta = (t % 2 == 0) ? -6 : 6;
            automaton.process((angle, frequency, deltaA,
                                deltaF) -> deltaF[base] += delta);
            automaton.step();
        }

        List<StructureFactor.DispersionPoint> points = sf.spectrum(field,
                                                                     StructureFactor.Direction.X100);
        StructureFactor.DispersionPoint k0 = points.get(0);
        assertEquals("a genuine, hand-injected quanta-transfer sequence (deltaF writes "
                    + "between snapshots via step(), not a frozen field) must show its "
                    + "energy at the Nyquist ridge, not omega=0 - bracketing the K=0 "
                    + "guard test with a real-dynamics positive control",
                    Math.PI, Math.abs(k0.omega()), 1e-9);
        assertTrue("peak fraction should be ~1.0 for a clean period-2 alternation, was "
                   + k0.peakFraction(), k0.peakFraction() >= 0.999);
    }
}
