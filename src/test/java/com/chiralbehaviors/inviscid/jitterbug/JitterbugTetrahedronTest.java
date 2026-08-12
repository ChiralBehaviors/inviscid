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

import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.L_EDGE;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugLinkage.jacobian;
import static com.chiralbehaviors.inviscid.jitterbug.JitterbugLinkage.residual;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * The tetrahedron inside the linkage.
 *
 * <p>
 * Items 9 and 10: 1344 hinge-consistent ways to seat two triangles on each of
 * the four faces of a regular tetrahedron, and the explicit placement satisfies
 * the hinge constraints exactly with every strut at the octahedron edge length.
 *
 * <p>
 * The tetrahedron is not on the symmetric path — it is a point of the
 * six-dimensional variety that the one-parameter family never visits.
 *
 * @author halhildebrand
 */
public class JitterbugTetrahedronTest {

    private static List<int[][]> ALL;
    private static List<int[][]> TETS;

    @BeforeClass
    public static void enumerate() {
        ALL = TetrahedronSeatings.all();
        TETS = TetrahedronSeatings.tetrahedra();
    }

    /**
     * Item 9. 1344 hinge-consistent assignments seat two triangles on each of
     * the four tetrahedron faces, out of 9216 with no occupancy rule at all.
     *
     * <p>
     * The full occupancy spectrum is asserted as well, not just the headline:
     * the five profiles must sum back to 9216, which is the arithmetic check
     * that no branch of the search was lost.
     */
    @Test
    public void thereAre1344TetrahedronSeatingsOutOf9216() {
        assertEquals("hinge-consistent assignments with no occupancy rule",
                     9216, ALL.size());
        assertEquals("seatings with 2 triangles on each of 4 faces", 1344,
                     TETS.size());

        Map<String, Integer> profiles = new TreeMap<>();
        for (int[][] s : ALL) {
            profiles.merge(Arrays.toString(TetrahedronSeatings
                                                              .occupancyProfile(s)),
                           1, Integer::sum);
        }
        Map<String, Integer> expected = new TreeMap<>();
        expected.put("[8]", 96);
        expected.put("[6, 2]", 1152);
        expected.put("[4, 4]", 864);
        expected.put("[4, 2, 2]", 5760);
        expected.put("[2, 2, 2, 2]", 1344);
        assertEquals("occupancy spectrum", expected, profiles);

        int sum = profiles.values().stream().mapToInt(Integer::intValue).sum();
        assertEquals("the spectrum must account for every assignment", 9216,
                     sum);
    }

    /**
     * Item 10. The explicit placement is exact: hinge residual identically zero,
     * every strut at the octahedron edge length, and the 24 corners sitting on
     * exactly 4 distinct points.
     *
     * <p>
     * Checked over <em>all</em> 1344 seatings rather than the first three the
     * reference sampled — the property is structural and should not depend on
     * which one is picked.
     */
    @Test
    public void everyTetrahedronSeatingSatisfiesTheHingesExactly() {
        for (int t = 0; t < TETS.size(); t++) {
            double[][][] x = TetrahedronSeatings.build(TETS.get(t));
            assertEquals("hinge residual for seating " + t, 0.0,
                         Linear.norm(residual(x)), 0.0);
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 3; j++) {
                    assertEquals("strut of seating " + t, L_EDGE,
                                 JitterbugGeometry.distance(x[i][j],
                                                            x[i][(j + 1) % 3]),
                                 1e-15);
                }
            }
            JitterbugGeometry.Clustering c = JitterbugGeometry.cluster(x, 1e-9);
            assertEquals("seating " + t + " must occupy 4 points", 4,
                         c.representatives().length);
            for (int m : c.multiplicities()) {
                assertEquals("6 corners at each tetrahedron vertex", 6, m);
            }
        }
    }

    /**
     * A correction to the reference. {@code jb_d_tet.py} reported "rank 35 at
     * the tetrahedron" from a single sampled seating and the memo recorded it as
     * a property of the configuration. It is not uniform: across all 1344
     * seatings both 35 and 36 occur.
     *
     * <p>
     * So the tetrahedron is not, as a class, a singular point of the variety —
     * some seatings are singular and some are ordinary. The distribution is
     * <b>576 seatings at rank 35 and 768 at rank 36</b>, and it is pinned below
     * rather than merely promised: the previous version of this test asserted
     * only that both keys occur, which is a strictly weaker statement than the
     * javadoc was making. The correct sentence for the record is "576 of the
     * 1344 tetrahedron seatings are singular points of the variety (rank 35) and
     * the other 768 are ordinary points (rank 36)"; the reference's sampled
     * seating was in the 43% minority.
     *
     * <p>
     * Unlike the reachability counts, this one is safe to pin: it is an
     * exhaustive enumeration over a finite exactly-constructed set, with no
     * iterative solver between the geometry and the answer.
     */
    @Test
    public void tetrahedronRankIsNotUniformAcrossSeatings() {
        Map<Integer, Integer> histogram = new HashMap<>();
        for (int[][] s : TETS) {
            int rank = Linear.rank(jacobian(TetrahedronSeatings.build(s)), 1e-9);
            histogram.merge(rank, 1, Integer::sum);
        }
        assertEquals("no rank outside {35,36}", 2, histogram.size());
        assertEquals("singular seatings (rank 35)", Integer.valueOf(576),
                     histogram.get(35));
        assertEquals("ordinary seatings (rank 36)", Integer.valueOf(768),
                     histogram.get(36));
        assertEquals("every seating accounted for", 1344,
                     histogram.values().stream().mapToInt(Integer::intValue)
                              .sum());

        // The reference's sampled seatings, individually, so the correction is
        // reproducible rather than aggregate.
        assertEquals(35, Linear.rank(jacobian(TetrahedronSeatings
                                                                 .build(TETS.get(0))),
                                     1e-9));
        assertEquals(35, Linear.rank(jacobian(TetrahedronSeatings
                                                                 .build(TETS.get(1))),
                                     1e-9));
        assertEquals("the seating that refutes 'the tetrahedron has rank 35'",
                     36, Linear.rank(jacobian(TetrahedronSeatings
                                                                 .build(TETS.get(2))),
                                     1e-9));
    }

    /**
     * The tetrahedron is genuinely off the symmetric path: no angle produces a
     * four-point configuration, so it cannot be reached by turning the
     * jitterbug's single handle.
     */
    @Test
    public void noAngleOnTheSymmetricPathIsATetrahedron() {
        List<Integer> counts = new ArrayList<>();
        for (int k = 0; k < 3600; k++) {
            int n = JitterbugGeometry
                                     .cluster(JitterbugGeometry.corners(k * 0.1),
                                              1e-7)
                                     .representatives().length;
            if (!counts.contains(n)) {
                counts.add(n);
            }
        }
        counts.sort(null);
        assertEquals("distinct-corner counts over the whole circle",
                     List.of(6, 12), counts);
    }

    /**
     * The enumeration must be doing real work, and this is a <b>magnitude
     * check</b> that it is: with no hinge constraint at all there are 24^8
     * assignments, and the consistency filter cuts that to 9216 — better than
     * seven orders of magnitude. If the filter were inert, 1344 would be an
     * accident of the search shape rather than a count of anything.
     *
     * <p>
     * It is not the corruption experiment an earlier version of this javadoc
     * described; no hinge is corrupted here. The test is nonetheless not
     * vacuous — weakening {@code TetrahedronSeatings.consistent()} raises
     * {@code all()} to 479232 and fails this assertion — but it is a coarse
     * check with real headroom, so read it as a floor and not as a
     * characterisation of the filter.
     */
    @Test
    public void theSeatingSearchActuallyEnforcesTheHinges() {
        // GUARD: a property of java.lang.Math, not of the enumeration. Present
        // so the 10^7 factor below has a stated numerator.
        double unconstrained = Math.pow(24, 8);
        assertEquals(1.10075314176E11, unconstrained, 1.0);

        assertTrue("the hinge check must be doing the work",
                   unconstrained / ALL.size() > 1.0e7);
    }
}
