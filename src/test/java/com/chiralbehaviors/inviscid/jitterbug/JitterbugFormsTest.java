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

package com.chiralbehaviors.inviscid.jitterbug;

import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.A_ICOSAHEDRON_DEG;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.L_EDGE;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.cluster;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.corners;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

/**
 * Fuller's volume sequence, measured on the actual parameterisation.
 *
 * <p>
 * Normalised so that the regular tetrahedron of edge L has volume 1, the
 * jitterbug's convex hull runs 20 (vector equilibrium) to 18.512296
 * (icosahedron) to 4 (octahedron) to 1 (tetrahedron).
 *
 * <p>
 * The hull is computed by {@link ConvexHullVolume}, a test-scope brute-force
 * supporting-plane enumeration written because commons-math3 has no 3D hull and
 * {@code src/main} is not allowed a new dependency. The helper is itself pinned
 * against three solids of known volume before it is trusted here.
 *
 * @author halhildebrand
 */
public class JitterbugFormsTest {

    /** Volume of a regular tetrahedron of edge L: {@code L^3 / (6 sqrt 2)}. */
    private static final double V_TET = L_EDGE * L_EDGE * L_EDGE
                                        / (6 * Math.sqrt(2.0));

    private static final double HULL_TOL = 1e-9;

    /**
     * The hull helper must be right before anything is concluded from it. A unit
     * cube, a regular tetrahedron of edge sqrt(2) and a unit-circumradius
     * octahedron, all with volumes known in closed form. The cube also exercises
     * the coplanar-facet path, which is the case an incremental hull would have
     * had to special-case.
     */
    @Test
    public void convexHullVolumeIsCorrectOnKnownSolids() {
        double[][] cube = new double[8][];
        int n = 0;
        for (int x = 0; x < 2; x++) {
            for (int y = 0; y < 2; y++) {
                for (int z = 0; z < 2; z++) {
                    cube[n++] = new double[] { x, y, z };
                }
            }
        }
        assertEquals("unit cube", 1.0,
                     ConvexHullVolume.volume(cube, HULL_TOL), 1e-12);

        assertEquals("regular tetrahedron of edge sqrt(2)", V_TET,
                     ConvexHullVolume.volume(TetrahedronSeatings.TET, HULL_TOL),
                     1e-12);
        assertEquals(1.0 / 3.0, V_TET, 1e-15);

        double[][] octa = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 },
                            { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
        assertEquals("octahedron of circumradius 1", 4.0 / 3.0,
                     ConvexHullVolume.volume(octa, HULL_TOL), 1e-12);

        // An interior point must not change the hull: guards against the
        // enumeration mistaking a non-supporting plane for a facet.
        double[][] cubePlusCentre = new double[9][];
        System.arraycopy(cube, 0, cubePlusCentre, 0, 8);
        cubePlusCentre[8] = new double[] { 0.5, 0.5, 0.5 };
        assertEquals(1.0, ConvexHullVolume.volume(cubePlusCentre, HULL_TOL),
                     1e-12);
    }

    /**
     * Fuller's sequence, reproducing the reference values 20 / 18.512296 / 4 / 1
     * to six decimals.
     */
    @Test
    public void volumeSequenceIsTwentyEighteenPointFiveFourAndOne() {
        assertEquals("vector equilibrium", 20.0, relativeVolume(0.0), 1e-9);
        assertEquals("icosahedron", 18.512296,
                     relativeVolume(A_ICOSAHEDRON_DEG), 1e-6);
        assertEquals("octahedron", 4.0, relativeVolume(60.0), 1e-9);
        assertEquals("tetrahedron",
                     1.0, ConvexHullVolume.volume(TetrahedronSeatings.TET,
                                                  HULL_TOL)
                          / V_TET,
                     1e-12);
    }

    /**
     * The vector equilibrium is a <b>local minimum</b> of hull volume, flanked
     * by two equal maxima.
     *
     * <p>
     * Written first as "volume falls monotonically from the VE to the
     * octahedron" and immediately falsified: V(1 degree) = 20.112946 &gt; V(0) =
     * 20. Ternary search over [0,20] puts a maximum at a = 7.2189519 degrees
     * with V/V_tet = 20.4430915, about 2.2% above Fuller's 20. Reproduced by the
     * Python reference to nine decimals.
     *
     * <p>
     * <b>"The maximum" was the wrong article.</b> V is even about a = 0 — checked
     * below at seven angles, agreeing to better than 1e-12 — so the same maximum
     * occurs at a = -7.2189519, and the VE sits at a critical point <em>between
     * them</em> which is a local minimum, not a peak. Over the full circle there
     * are four such maxima. The searches here cover [-20, 20] and say nothing
     * about the rest of the circle.
     *
     * <p>
     * Not a contradiction of the recorded 20: that is the vector equilibrium's
     * own volume, and it is exact. The measured explanation is
     * {@link #theTwelveVerticesAreCosphericalAtEveryAngle()}: all 12 vertices lie
     * on a common sphere at every angle, and its radius falls monotonically from
     * 1.4142136 at the VE to 1.0 at the octahedron. The cuboctahedron is a poor
     * volume-maximiser among 12 points on a sphere, so for the first few degrees
     * the gain from a better-distributed shape beats the loss from the shrinking
     * radius; past 7.2 degrees the radius wins and V falls the rest of the way.
     *
     * <p>
     * An earlier version of this javadoc explained the excess by saying the hull
     * volume and the jitterbug's own enclosed volume "part company as soon as the
     * squares buckle". That reconciliation is <b>unmeasured and is withdrawn</b>:
     * a probe found the own enclosed volume also exceeds 20 under one natural
     * triangulation of the skew openings and falls monotonically under the other,
     * so "own enclosed volume" is not even well defined here without naming a
     * triangulation.
     */
    @Test
    public void hullVolumePeaksOffTheVectorEquilibrium() {
        assertEquals("V just past the VE exceeds V at the VE", 20.112946,
                     relativeVolume(1.0), 1e-6);

        double peak = ternaryPeak(0.0, 20.0);
        assertEquals("hull volume peaks at 7.2189519 degrees", 7.2189519, peak,
                     1e-5);
        assertEquals("peak hull volume", 20.4430915, relativeVolume(peak),
                     1e-6);

        // Even about a = 0, so the VE is a local MINIMUM between two maxima.
        for (double a : new double[] { 0.001, 0.1, 1.0, 3.0, 7.2189519, 12.0,
                                       20.0 }) {
            assertEquals("V must be even about the VE at a=" + a,
                         relativeVolume(a), relativeVolume(-a), 1e-12);
        }
        double mirrored = ternaryPeak(-20.0, 0.0);
        assertEquals("the mirror maximum", -7.2189519, mirrored, 1e-5);
        assertEquals("at the same height", relativeVolume(peak),
                     relativeVolume(mirrored), 1e-9);
        assertEquals("and the VE is the local minimum between them", 20.0,
                     relativeVolume(0.0), 1e-9);
        org.junit.Assert.assertTrue("strictly a minimum, not a plateau",
                                    relativeVolume(0.001) > 20.0
                                                                 && relativeVolume(-0.001) > 20.0);

        // Strictly rising up to the peak, strictly falling from it to the
        // octahedron. Sampled, so a wobble in the hull routine would show.
        double previous = -1;
        for (double a = 0; a <= 7.0; a += 0.5) {
            double v = relativeVolume(a);
            org.junit.Assert.assertTrue("volume must rise up to the peak; a=" + a,
                                        v > previous);
            previous = v;
        }
        previous = Double.MAX_VALUE;
        for (int k = 8; k <= 60; k++) {
            double v = relativeVolume(k);
            org.junit.Assert.assertTrue("volume must fall past the peak; a=" + k
                                        + " got " + v + " after " + previous,
                                        v < previous);
            previous = v;
        }
        assertEquals(4.0, previous, 1e-9);
    }

    /**
     * The measured reason the hull can exceed the cuboctahedron's own 20: the
     * 12 vertices never leave a common sphere, so the only things trading off
     * are the sphere's radius, which falls monotonically, and how well the 12
     * points are spread on it.
     */
    @Test
    public void theTwelveVerticesAreCosphericalAtEveryAngle() {
        // Seeded from the actual first sample, so the very first comparison in
        // the loop is a real one.
        double previous = radius(0.0);
        assertEquals("circumradius at the vector equilibrium", Math.sqrt(2.0),
                     previous, 1e-12);
        for (int k = 1; k <= 600; k++) {
            double a = k * 0.1;
            double r = radius(a);
            org.junit.Assert.assertTrue("circumradius must fall at a=" + a
                                        + ": " + r + " after " + previous,
                                        r < previous);
            previous = r;
        }
        assertEquals("circumradius at the icosahedron", 1.3449970239,
                     radius(A_ICOSAHEDRON_DEG), 1e-9);
        assertEquals("circumradius at the octahedron", 1.0, radius(60.0),
                     1e-12);
    }

    /**
     * The common distance of all 24 corners from the origin at {@code aDeg},
     * asserting on the way that there is in fact only one such distance. Taken
     * over the corners rather than the clustered representatives so the check
     * still means the same thing at a = 60, where the 12 vertices collapse onto
     * 6.
     */
    private static double radius(double aDeg) {
        double[][][] x = corners(aDeg);
        double lo = Double.MAX_VALUE;
        double hi = 0;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                double[] p = x[i][j];
                double r = Math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                lo = Math.min(lo, r);
                hi = Math.max(hi, r);
            }
        }
        assertEquals("every corner on one sphere at a=" + aDeg, lo, hi, 1e-12);
        return hi;
    }

    /** Ternary search for the interior maximum of {@link #relativeVolume}. */
    private static double ternaryPeak(double from, double to) {
        double lo = from;
        double hi = to;
        for (int i = 0; i < 200; i++) {
            double m1 = lo + (hi - lo) / 3.0;
            double m2 = hi - (hi - lo) / 3.0;
            if (relativeVolume(m1) < relativeVolume(m2)) {
                lo = m1;
            } else {
                hi = m2;
            }
        }
        return 0.5 * (lo + hi);
    }

    private static double relativeVolume(double aDeg) {
        return ConvexHullVolume.volume(cluster(corners(aDeg), 1e-7)
                                                                  .representatives(),
                                       HULL_TOL)
               / V_TET;
    }
}
