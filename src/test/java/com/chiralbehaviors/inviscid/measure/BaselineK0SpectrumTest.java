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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.OptionalDouble;

import javax.vecmath.Point3i;

import org.junit.BeforeClass;
import org.junit.Test;

import com.chiralbehaviors.inviscid.measure.BaselineSpectrumHarness.MemberSpectrum;
import com.chiralbehaviors.inviscid.measure.BaselineSpectrumHarness.Result;

/**
 * B.2 (bead inviscid-0nx.7): the K=0 (collision-free) baseline spectrum,
 * captured BEFORE any collision rule exists (bead inviscid-0nx.14).
 * Deliberately sequenced ahead of that bead so the "collision-broadened
 * lines" claim it will make is falsifiable against a golden artifact that
 * predates any collision code.
 *
 * @author halhildebrand
 */
public class BaselineK0SpectrumTest {

    private static final double MIN_CONCENTRATION = 0.95;

    private static Result       result;

    @BeforeClass
    public static void generate() {
        result = BaselineSpectrumHarness.run();
    }

    /**
     * NON-VACUITY (hard requirement): at K=0 there is no collision rule
     * (Necronomata.process(Point3i) is still a no-op), so deltaF is always
     * zero and the frequency (quanta) field is bit-identical before and
     * after the entire run - by construction, not by coincidence. A naive
     * transport statistic (e.g. a before/after variance ratio) would
     * silently compute exactly 1.0 on identical arrays, which reads as "a
     * clean isotropic diffusion coefficient" - a false positive, since no
     * transport process exists at all. BaselineSpectrumHarness.quantaSpreadRatio
     * must guard against this and report DEGENERATE (empty), never 1.0.
     */
    @Test
    public void k0HasZeroTransport() {
        assertTrue("frequency (quanta) field changed during a K=0 (collision-free) run - "
                   + "transport should be impossible with no collision rule",
                   result.frequencyFieldUnchanged);

        OptionalDouble ratio = BaselineSpectrumHarness.quantaSpreadRatio(result.frequencyBefore,
                                                                          result.frequencyAfter);
        assertFalse("transport statistic returned a value despite zero collisions - "
                   + "a clean number here (e.g. a naive ratio of identical arrays would be "
                   + "exactly 1.0) is a false positive that would misreport 'no transport' as "
                   + "'normal isotropic transport'; it must be reported as DEGENERATE/undefined",
                   ratio.isPresent());
    }

    /**
     * Ties the measured spectral peak back to the seeded integer quanta
     * count via SpectrumAnalyzer.expectedBinForFrequency - the shared
     * P0.1-aligned bin<->frequency convention (not re-derived here). See
     * BaselineSpectrumHarness's class javadoc for why the expected bin
     * uses {@code quanta * STRIDE}, not {@code quanta} directly.
     */
    @Test
    public void k0LineFrequencyEqualsQuantaTimesQuantumRate() {
        int checked = 0;
        for (MemberSpectrum m : result.members) {
            if (m.quanta == 0f) {
                continue;
            }
            assertEquals("member " + m.memberIndex + " (quanta=" + m.quanta
                        + ") peak bin does not match the seeded quanta count",
                        m.expectedBin, m.peakBin);
            checked++;
        }
        assertTrue("no nonzero-frequency members were checked - test is vacuous",
                   checked > 0);
    }

    /**
     * Every nonzero-frequency member's spectrum must be a single dominant
     * line: peak-bin power fraction >= 0.95 (numeric, asserted directly -
     * not "a peak exists"), and spectral entropy below a threshold DERIVED
     * from that same 0.95 concentration bound (see
     * BaselineSpectrumHarness#maxSpectralEntropyForConcentration), not an
     * arbitrary magic number.
     */
    @Test
    public void k0SpectrumIsPureTones() {
        double entropyBound = BaselineSpectrumHarness.maxSpectralEntropyForConcentration(MIN_CONCENTRATION,
                                                                                           result.fftLength);
        int checked = 0;
        for (MemberSpectrum m : result.members) {
            if (m.quanta == 0f) {
                continue;
            }
            assertTrue("member " + m.memberIndex + " (quanta=" + m.quanta
                       + ") peak-bin power fraction " + m.peakFraction
                       + " is below the required " + MIN_CONCENTRATION,
                       m.peakFraction >= MIN_CONCENTRATION);
            assertTrue("member " + m.memberIndex + " (quanta=" + m.quanta
                       + ") spectral entropy " + m.spectralEntropy
                       + " exceeds the derived bound " + entropyBound
                       + " for >=" + MIN_CONCENTRATION + " concentration",
                       m.spectralEntropy <= entropyBound);
            checked++;
        }
        assertTrue("no nonzero-frequency members were checked - test is vacuous",
                   checked > 0);
    }

    /**
     * Regenerates the baseline in-test and compares it, field by field,
     * against the committed golden TSV within the documented numeric
     * tolerance (byte/string-exact float comparison is deliberately
     * avoided - see the artifact's provenance header and the
     * inviscid-0nx.7 plan-audit correction on this bead).
     */
    @Test
    public void goldenArtifactMatchesRegeneration() throws IOException {
        List<String[]> golden = readGoldenDataRows();
        List<String[]> regenerated = BaselineSpectrumHarness.toDataRows(result);

        assertEquals("regenerate with BaselineSpectrumHarness and review the diff: "
                    + "golden artifact has a different number of data rows than a fresh run",
                    golden.size(), regenerated.size());

        for (int i = 0; i < golden.size(); i++) {
            String[] g = golden.get(i);
            String[] r = regenerated.get(i);
            String rowContext = "row " + i + " (member " + g[0] + ")";

            assertEquals("regenerate with BaselineSpectrumHarness and review the diff: "
                        + rowContext + " memberIndex mismatch", g[0], r[0]);
            assertEquals("regenerate with BaselineSpectrumHarness and review the diff: "
                        + rowContext + " quanta mismatch", g[1], r[1]);
            assertEquals("regenerate with BaselineSpectrumHarness and review the diff: "
                        + rowContext + " expectedBin mismatch", g[2], r[2]);
            assertEquals("regenerate with BaselineSpectrumHarness and review the diff: "
                        + rowContext + " peakBin mismatch", g[3], r[3]);

            double gPeakFraction = Double.parseDouble(g[4]);
            double rPeakFraction = Double.parseDouble(r[4]);
            if (Math.abs(gPeakFraction - rPeakFraction) > BaselineSpectrumHarness.TOLERANCE) {
                fail("regenerate with BaselineSpectrumHarness and review the diff: "
                   + rowContext + " peakFraction differs by more than " + BaselineSpectrumHarness.TOLERANCE
                   + " (golden=" + gPeakFraction + ", regenerated=" + rPeakFraction + ")");
            }

            double gEntropy = Double.parseDouble(g[5]);
            double rEntropy = Double.parseDouble(r[5]);
            if (Math.abs(gEntropy - rEntropy) > BaselineSpectrumHarness.TOLERANCE) {
                fail("regenerate with BaselineSpectrumHarness and review the diff: "
                   + rowContext + " spectralEntropy differs by more than " + BaselineSpectrumHarness.TOLERANCE
                   + " (golden=" + gEntropy + ", regenerated=" + rEntropy + ")");
            }
        }
    }

    /**
     * FIX 1 (stacked review, code-review-expert Important #1): {@code run}
     * must validate its parameters up front and fail near the mistake,
     * not mysteriously inside {@link Fft} or via a silently-wrong
     * peak-bin result. In particular quanta at or beyond the Nyquist/
     * aliasing bound {@code PHASE_RESOLUTION/(2*stride)} (8 at the
     * default stride=225 - see BaselineSpectrumHarness's class javadoc
     * "Nyquist / aliasing bound on seeded quanta") must be rejected, and
     * the failure message must name the bound.
     */
    @Test
    public void run_rejectsParametersThatViolateItsPreconditions() {
        assertRejectedWithBound("maxQuanta at the Nyquist bound (8) must be rejected",
                                 2, 2, 2, 42L, 256, 225, -5, 8);
        assertRejectedWithBound("minQuanta at the negative Nyquist bound (-8) must be rejected",
                                 2, 2, 2, 42L, 256, 225, -8, 5);

        try {
            BaselineSpectrumHarness.run(new Point3i(2, 2, 2), 42L, 300, 225,
                                         -5, 5);
            fail("fftLength=300 (not a power of two) must be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue("message should mention fftLength: " + expected.getMessage(),
                       expected.getMessage().contains("fftLength"));
        }

        try {
            BaselineSpectrumHarness.run(new Point3i(2, 2, 2), 42L, 256, 0,
                                         -5, 5);
            fail("stride=0 must be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue("message should mention stride: " + expected.getMessage(),
                       expected.getMessage().contains("stride"));
        }

        try {
            BaselineSpectrumHarness.run(new Point3i(2, 2, 2), 42L, 256, 225,
                                         5, -5);
            fail("minQuanta > maxQuanta must be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue("message should mention minQuanta/maxQuanta: "
                       + expected.getMessage(),
                       expected.getMessage().contains("minQuanta"));
        }
    }

    private static void assertRejectedWithBound(String description, int ex,
                                                  int ey, int ez, long seed,
                                                  int fftLength, int stride,
                                                  int minQuanta,
                                                  int maxQuanta) {
        try {
            BaselineSpectrumHarness.run(new Point3i(ex, ey, ez), seed,
                                         fftLength, stride, minQuanta,
                                         maxQuanta);
            fail(description);
        } catch (IllegalArgumentException expected) {
            assertTrue("message should state the Nyquist/aliasing bound (8): "
                       + expected.getMessage(),
                       expected.getMessage().contains("8"));
        }
    }

    private static List<String[]> readGoldenDataRows() throws IOException {
        URL resource = BaselineK0SpectrumTest.class.getClassLoader()
                                                     .getResource("lga/baseline-k0-spectrum.tsv");
        if (resource == null) {
            fail("regenerate with BaselineSpectrumHarness and review the diff: "
               + "golden artifact src/test/resources/lga/baseline-k0-spectrum.tsv is missing");
        }
        Path path = Paths.get(resource.getPath());
        List<String[]> rows = new ArrayList<>();
        for (String line : Files.readAllLines(path)) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            rows.add(line.split("\t"));
        }
        return rows;
    }
}
