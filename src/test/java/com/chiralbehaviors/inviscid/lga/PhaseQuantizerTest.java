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

package com.chiralbehaviors.inviscid.lga;

import static com.chiralbehaviors.inviscid.Constants.TWO_PI;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Failing-tests-first (TDD) coverage for {@link PhaseQuantizer}: phase
 * binning of continuous angles into N_lga discrete bins, and the
 * precomputed member-segment geometry LUT keyed by (cube, member, bin).
 *
 * <p>N_lga, the geometry LUT resolution, and the member radius are all
 * sourced from the single committed atlas header ({@code
 * contact-atlas-v2.tsv}, bead inviscid-0nx.16) - the same idiom {@code
 * CommittedContactAtlasTest} uses - rather than hardcoded a second time
 * here.
 *
 * @author halhildebrand
 */
public class PhaseQuantizerTest {

    private static final String RESOURCE_PATH = "lga/contact-atlas-v2.tsv";
    private static final double DELTA         = 1e-6;

    private static ContactAtlas.Header HEADER;

    @BeforeClass
    public static void loadHeader() throws IOException {
        URL resource = PhaseQuantizerTest.class.getClassLoader()
                                                .getResource(RESOURCE_PATH);
        if (resource == null) {
            fail("committed atlas src/test/resources/lga/contact-atlas-v2.tsv is missing");
        }
        Path path = Paths.get(resource.getPath());
        HEADER = ContactAtlas.read(path).header();
    }

    @Test
    public void binsPartitionTheCircleExactly() {
        PhaseQuantizer quantizer = PhaseQuantizer.of(HEADER);
        int nLga = quantizer.nLga();
        float binWidth = TWO_PI / nLga;

        boolean[] visited = new boolean[nLga];
        int previous = -1;
        for (float angle = 0f; angle < TWO_PI; angle += binWidth / 37f) {
            int bin = quantizer.bin(angle);
            assertTrue("bin out of range: " + bin, bin >= 0 && bin < nLga);
            visited[bin] = true;
            if (previous != -1) {
                assertTrue("bin regressed at angle " + angle + ": " + previous
                          + " -> " + bin, bin >= previous);
            }
            previous = bin;
        }
        // Honestly cover the top of the circle, including the largest
        // float strictly below 2*pi - the exact case the code-review
        // clamp fix (topOfCircleMapsToLastBinNotFirst) addresses.
        int lastBin = quantizer.bin(Math.nextDown(TWO_PI));
        assertEquals("largest float below 2*pi must be the last bin", nLga - 1, lastBin);
        visited[lastBin] = true;
        assertTrue("last bin regressed relative to the sweep", lastBin >= previous);

        assertEquals("last bin never reached", nLga - 1, lastBin);
        for (int i = 0; i < nLga; i++) {
            assertTrue("bin " + i + " never visited", visited[i]);
        }
    }

    @Test
    public void negativeAndOverflowAnglesNormalise() {
        PhaseQuantizer quantizer = PhaseQuantizer.of(HEADER);

        assertEquals("negative angle should land in the same bin as just-under-2*pi",
                    quantizer.bin(TWO_PI - 0.1f), quantizer.bin(-0.1f));
        assertEquals("2*pi-overflowed angle should land in the same bin as its residue",
                    quantizer.bin(0.1f),
                    quantizer.bin((float) (2 * Math.PI) + 0.1f));
    }

    /**
     * Regression for the code-review Important finding on inviscid-0nx.18:
     * {@code bin()}'s final division was entirely float-precision, unlike
     * every sibling ({@code MemberGeometry.stepOf}, {@code
     * NecronomataVisualization.setState}, {@code
     * ContactComboCache#angleOf}) which does the equivalent division in
     * double. For the largest float below {@code TWO_PI}, {@code
     * (int)(normalized/binWidth)} can equal {@code nLga} itself (a
     * division-rounding artifact, not a real 25th bin), and a bare {@code
     * % nLga} then wraps that near-2*pi sliver to bin 0 instead of the
     * true last bin. {@code nLga=360} is a reviewer-verified trigger under
     * strict float32 (this test's own {@code nLga=24} production value
     * happens not to trigger it).
     */
    @Test
    public void topOfCircleMapsToLastBinNotFirst() {
        int triggeringNLga = 360;
        PhaseQuantizer quantizer = new PhaseQuantizer(triggeringNLga,
                                                        HEADER.geometryResolution(),
                                                        HEADER.memberRadius());
        float largestBelowTwoPi = Math.nextDown(TWO_PI);
        assertEquals("largest float below 2*pi must land in the last bin, not wrap to bin 0",
                    triggeringNLga - 1, quantizer.bin(largestBelowTwoPi));
    }

    @Test
    public void centreThenBinIsIdentity() {
        PhaseQuantizer quantizer = PhaseQuantizer.of(HEADER);
        for (int bin = 0; bin < quantizer.nLga(); bin++) {
            assertEquals("bin " + bin, bin, quantizer.bin(quantizer.centre(bin)));
        }
    }

    /**
     * Division of responsibility (substantive-critique Minor finding 2):
     * this test verifies {@code segment()}/{@code buildLut}'s WIRING -
     * that the precomputed LUT entry equals {@code
     * geometry.memberSegment(cube, member, centre)} for the SAME {@code
     * centre = quantizer.centre(bin)} value fed to both sides. It does
     * NOT independently verify that {@code centre()} derives the
     * geometrically-true bin-center angle - both expected and actual are
     * downstream of the identical {@code centre(bin)} call, so a
     * systematic derivation error in {@code centre()} (e.g. left-edge
     * instead of center) would not be caught here. That independent
     * geometric check is {@link #quantisationErrorIsBoundedByHalfBin()},
     * which sweeps angles directly and bounds the error against the
     * half-bin claim without routing through the LUT.
     */
    @Test
    public void lutMatchesMemberGeometryAtBinCentres() {
        PhaseQuantizer quantizer = PhaseQuantizer.of(HEADER);
        MemberGeometry geometry = new MemberGeometry(HEADER.geometryResolution(),
                                                       HEADER.memberRadius());

        for (int cube = 0; cube < 5; cube++) {
            for (int member = 0; member < 6; member++) {
                for (int bin = 0; bin < quantizer.nLga(); bin++) {
                    float centre = quantizer.centre(bin);
                    Segment expected = geometry.memberSegment(cube, member, centre);
                    Segment actual = quantizer.segment(cube, member, bin);
                    String context = "cube " + cube + " member " + member + " bin "
                                    + bin;
                    assertEquals(context + " a.x", expected.getA().x,
                                actual.getA().x, DELTA);
                    assertEquals(context + " a.y", expected.getA().y,
                                actual.getA().y, DELTA);
                    assertEquals(context + " a.z", expected.getA().z,
                                actual.getA().z, DELTA);
                    assertEquals(context + " b.x", expected.getB().x,
                                actual.getB().x, DELTA);
                    assertEquals(context + " b.y", expected.getB().y,
                                actual.getB().y, DELTA);
                    assertEquals(context + " b.z", expected.getB().z,
                                actual.getB().z, DELTA);
                }
            }
        }
    }

    @Test
    public void quantisationErrorIsBoundedByHalfBin() {
        PhaseQuantizer quantizer = PhaseQuantizer.of(HEADER);
        int nLga = quantizer.nLga();
        float binWidth = TWO_PI / nLga;
        float halfBin = binWidth / 2f;

        for (float angle = 0f; angle < TWO_PI; angle += 0.013f) {
            int bin = quantizer.bin(angle);
            float centre = quantizer.centre(bin);
            float diff = Math.abs(angle - centre);
            float circular = Math.min(diff, TWO_PI - diff);
            assertTrue("angle " + angle + " bin " + bin + " centre " + centre
                      + " circular error " + circular + " exceeds half-bin "
                      + halfBin, circular <= halfBin + 1e-4f);
        }
    }

    @Test(expected = IllegalArgumentException.class)
    public void constructorRejectsNonPositiveNLga() {
        new PhaseQuantizer(0, HEADER.geometryResolution(), HEADER.memberRadius());
    }

    /**
     * Regression for substantive-critique Significant finding 1
     * (inviscid-0nx.18, post-code-review): {@code bin(angle)} and {@code
     * ContactAtlasGenerator.binOfStep(step)} are independent derivations
     * that agree only when {@code nLga} divides {@code
     * geometryResolution}. Mutation-verified by the critic: {@code
     * nLga=100, geometryResolution=360} produces 40/360 mismatches (today's
     * production {@code nLga=24} agrees only because 24 happens to divide
     * 360). A future N_lga revision that isn't a divisor would silently
     * misalign live {@code bin()} classification against the
     * inviscid-0nx.19-transcribed contact table, so this must fail loudly
     * at construction rather than silently corrupt collision lookups.
     */
    @Test(expected = IllegalArgumentException.class)
    public void constructorRejectsNLgaThatDoesNotDivideGeometryResolution() {
        new PhaseQuantizer(100, HEADER.geometryResolution(), HEADER.memberRadius());
    }

    @Test(expected = IllegalArgumentException.class)
    public void segmentRejectsNegativeCube() {
        PhaseQuantizer quantizer = PhaseQuantizer.of(HEADER);
        quantizer.segment(-1, 0, 0);
    }

    @Test(expected = IllegalArgumentException.class)
    public void segmentRejectsOutOfRangeMember() {
        PhaseQuantizer quantizer = PhaseQuantizer.of(HEADER);
        quantizer.segment(0, 6, 0);
    }

    /**
     * Pins the seam inviscid-0nx.19 (contact table transcription) will
     * stand on (substantive-critique Significant finding 1): {@code
     * PhaseQuantizer.bin(angle)} (continuous-angle division) and {@code
     * ContactAtlasGenerator.binOfStep(step)} (integer step-index division)
     * are independent derivations of "which N_lga bin does this geometry
     * step belong to" - this asserts they agree for EVERY step at the
     * committed atlas header's values. Passes immediately on current code
     * (nLga=24 divides geometryResolution=360); this is a pin against
     * regression, not a fix for a bug - the divisibility guard on the
     * constructor ({@code
     * constructorRejectsNLgaThatDoesNotDivideGeometryResolution}) is what
     * actually prevents a future non-divisor N_lga from breaking this
     * equivalence.
     */
    @Test
    public void binAgreesWithAtlasBinOfStepAtCommittedValues() {
        PhaseQuantizer quantizer = PhaseQuantizer.of(HEADER);
        int geometryResolution = HEADER.geometryResolution();
        int nLga = quantizer.nLga();

        for (int step = 0; step < geometryResolution; step++) {
            int expected = ContactAtlasGenerator.binOfStep(step, nLga, geometryResolution);
            int actual = quantizer.bin(ContactComboCache.angleOf(step, geometryResolution));
            assertEquals("step " + step, expected, actual);
        }
    }
}
