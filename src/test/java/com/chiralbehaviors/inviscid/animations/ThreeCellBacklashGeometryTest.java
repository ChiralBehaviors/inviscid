/**
 * Copyright (c) 2026 Chiral Behaviors, LLC, all rights reserved.
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

package com.chiralbehaviors.inviscid.animations;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry;

/**
 * The geometric premise of {@link ThreeCellBacklashAnimation}, made falsifiable.
 *
 * <p>
 * An animation is a picture, and a picture that is subtly wrong looks exactly
 * like one that is right — which is the failure this project has paid for
 * repeatedly ("a redundancy check on POSITIONS does not validate IDENTITY", and
 * the mirrored-sigma trap that makes a rendered cell the mirror of the computed
 * one). So the claim the animation rests on is asserted here rather than
 * eyeballed: three cells at FIXED sites, one dowel spacing apart, with phases
 * alternating {@code a, a+60, a}, actually SHARE their triangular faces.
 *
 * <p>
 * If this test fails, the animation is drawing a structure that does not close,
 * and nothing seen in it means anything.
 */
public class ThreeCellBacklashGeometryTest {

    private static final double A = ThreeCellBacklashAnimation.A_REF;

    private static final int[][] PERMS = { { 0, 1, 2 }, { 0, 2, 1 },
                                           { 1, 0, 2 }, { 1, 2, 0 },
                                           { 2, 0, 1 }, { 2, 1, 0 } };

    /** Site spacing: the two radial distances to the shared face, summed. */
    private static double sep0() {
        return ThreeCellBacklashAnimation.radial(Math.toRadians(A))
               + ThreeCellBacklashAnimation.radial(Math.toRadians(A + 60.0));
    }

    /**
     * Worst corner discrepancy across the shared face between a cell at phase
     * {@code gA} sitting at the origin and one at {@code gB} sitting one site
     * along the dowel. Resolved by permutation search, because which corner
     * meets which is a bookkeeping question on an equilateral triangle -- the
     * same resolution jb_hc's {@code face_pairing} uses.
     */
    private static double mismatch(double gA, double gB) {
        double n = Math.sqrt(3.0);
        double[] axis = { 1 / n, 1 / n, 1 / n };
        int fv = ThreeCellBacklashAnimation.faceToward(new double[] { 1, 1, 1 });
        int fo = ThreeCellBacklashAnimation.faceToward(new double[] { -1, -1,
                                                                      -1 });
        double[][][] a = JitterbugGeometry.corners(gA);
        double[][][] b = JitterbugGeometry.corners(gB);
        double best = Double.POSITIVE_INFINITY;
        for (int[] p : PERMS) {
            double worst = 0;
            for (int c = 0; c < 3; c++) {
                double d2 = 0;
                for (int t = 0; t < 3; t++) {
                    double u = a[fv][p[c]][t];
                    double v = b[fo][c][t] + sep0() * axis[t];
                    d2 += (u - v) * (u - v);
                }
                worst = Math.max(worst, Math.sqrt(d2));
            }
            best = Math.min(best, worst);
        }
        return best;
    }

    /**
     * The spacing is not a fitted number: it is the honeycomb's own lattice
     * step along a body diagonal. jb_hc measures {@code lattice(-30)} as the
     * fold half-diagonal 1.154700538, and the offset it places a triangular
     * neighbour at is that times the UNNORMALISED (1,1,1), i.e. exactly 2.
     */
    @Test
    public void siteSpacingIsTheHoneycombsOwnLatticeStep() {
        assertEquals(2.0, sep0(), 1e-12);
    }

    /** The two faces the dowel runs through are genuinely antipodal. */
    @Test
    public void theDowelPassesThroughOppositeFaces() {
        int fv = ThreeCellBacklashAnimation.faceToward(new double[] { 1, 1, 1 });
        int fo = ThreeCellBacklashAnimation.faceToward(new double[] { -1, -1,
                                                                      -1 });
        double[] u = JitterbugGeometry.faceAxis(fv);
        double[] w = JitterbugGeometry.faceAxis(fo);
        double dot = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
        assertEquals("the shared-face axes must be exactly antipodal", -1.0,
                     dot, 1e-12);
    }

    /**
     * THE LOAD-BEARING ONE. At the reference the faces coincide to machine
     * precision, in BOTH joint orientations -- cell 0 to cell 1 is VE-to-hole
     * and cell 1 to cell 2 is hole-to-VE, and a placement that closed only one
     * of them would draw a chain that comes apart at every other joint.
     */
    @Test
    public void theThreeCellChainActuallyClosesItsSharedFaces() {
        assertTrue("cell 0 -> 1 must close", mismatch(A, A + 60.0) < 1e-12);
        assertTrue("cell 1 -> 2 must close", mismatch(A + 60.0, A) < 1e-12);
    }

    /**
     * And the joint is what absorbs a phase difference. This is the animation's
     * whole subject, so it is asserted rather than assumed: the mismatch is
     * FIRST ORDER in the relative phase, which is why play buys a bounded
     * gradient at all (jb_pr, T2 [23654]). Two-sided -- it must be zero at zero
     * and grow proportionally, so a model that absorbed the difference for free
     * would fail here.
     */
    @Test
    public void aPhaseDifferenceCostsClearanceAtFirstOrder() {
        double half = mismatch(A, A + 60.0 + 0.5);
        double one = mismatch(A, A + 60.0 + 1.0);
        double two = mismatch(A, A + 60.0 + 2.0);
        assertTrue("a phase difference must cost something", half > 1e-6);
        assertEquals("mismatch must double when the phase difference doubles",
                     2.0, one / half, 0.02);
        assertEquals(2.0, two / one, 0.02);
    }

    /**
     * The coefficient jb_pr measures as EL/sqrt(2) = 1 per radian, reproduced
     * from the Java geometry. If the Java and Python sides disagreed here, the
     * animation would be showing a different medium from the one the harness
     * measures.
     */
    @Test
    public void theClearanceCostPerRadianMatchesTheHarness() {
        double d = 1e-4;
        double perRadian = mismatch(A, A + 60.0 + d) / Math.toRadians(d);
        assertEquals("EL/sqrt(2) per radian, as jb_pr R1 measures",
                     JitterbugGeometry.L_EDGE / Math.sqrt(2.0), perRadian,
                     1e-5);
    }

    // ---------------------------------------------------------------------
    // The DYNAMICS. Correct geometry with a mistranscribed integrator would
    // still look plausible on screen, so the Java transcription of jb_ct's
    // law is checked against the three things jb_ct itself gates.
    // ---------------------------------------------------------------------

    /** Run the animation's own integrator headlessly. Returns arrival times
     *  (NaN where a cell never moved) and the relative energy drift. */
    private static double[] run(double tmax, double h) {
        int n = ThreeCellBacklashAnimation.CELLS;
        double sep0 = sep0();
        double[] g = new double[n];
        double[] gd = new double[n];
        for (int k = 0; k < n; k++) {
            g[k] = Math.toRadians(A + 60.0 * (k % 2));
        }
        gd[0] = ThreeCellBacklashAnimation.KICK;
        double e0 = ThreeCellBacklashAnimation.energy(g, gd);
        double[] arrive = new double[n];
        java.util.Arrays.fill(arrive, Double.NaN);
        arrive[0] = 0.0;
        double now = 0, worst = 0;
        while (now < tmax) {
            ThreeCellBacklashAnimation.advance(g, gd, sep0, h);
            now += h;
            for (int k = 1; k < n; k++) {
                if (Double.isNaN(arrive[k])
                    && Math.abs(gd[k]) > 0.02 * ThreeCellBacklashAnimation.KICK) {
                    arrive[k] = now;
                }
            }
            for (int k = 0; k + 1 < n; k++) {
                worst = Math.max(worst,
                                 Math.abs(ThreeCellBacklashAnimation.gap(g, k, sep0))
                                 - ThreeCellBacklashAnimation.PLAY);
            }
        }
        double drift = Math.abs(ThreeCellBacklashAnimation.energy(g, gd) - e0) / e0;
        return new double[] { arrive[1], arrive[2], drift, worst };
    }

    /**
     * THE CLAIM THE ANIMATION IS FOR. The kick reaches cell 1 before cell 2,
     * and neither arrives at t = 0 -- which is the whole difference from the
     * rigid model, where the constraint loads every cell before any time
     * passes. Two-sided: an instantaneous onset fails here, and so does a
     * disturbance that never crosses.
     */
    @Test
    public void theKickTakesTimeToCrossAndArrivesInOrder() {
        double[] r = run(6.0, 1e-4);
        assertTrue("cell 1 must be reached", !Double.isNaN(r[0]));
        assertTrue("cell 2 must be reached", !Double.isNaN(r[1]));
        assertTrue("the onset is NOT instantaneous at cell 1", r[0] > 1e-3);
        assertTrue("cell 2 must arrive strictly after cell 1", r[1] > r[0]);
    }

    /**
     * V = 0, so energy is the only audit there is, and a contact chain is
     * exactly where a careless integrator loses it. jb_ct records that nudging
     * a position back onto the stop instead of resolving the impact bleeds 20
     * to 40 percent -- this is the row that would catch that transcription.
     */
    @Test
    public void energySurvivesTheContactsAndTheJointsKeepTheirPlay() {
        double[] r = run(6.0, 1e-4);
        assertTrue("energy drift " + r[2] + " must stay negligible", r[2] < 1e-9);
        assertTrue("no joint may exceed its play, worst overshoot " + r[3],
                   r[3] < 1e-9);
    }
}
