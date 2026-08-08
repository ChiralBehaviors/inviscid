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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.BeforeClass;
import org.junit.Test;

import com.chiralbehaviors.inviscid.PhiCoordinates;

/**
 * The COMMITTED contact atlas ({@code src/test/resources/lga/contact-atlas-v1.tsv},
 * bead inviscid-0nx.16 Stage 2): validates the actual Phase A -> Phase C
 * handoff artifact at the USER-DECIDED N_lga=24 (recorded on the bead's
 * {@code --design} field, full justification in T2 {@code
 * inviscid/analysis-nlga-candidates.md}) - as distinct from {@link
 * ContactAtlasTest}'s generated-in-test contract tests at a fixed,
 * clearly-illustrative {@code TEST_N_LGA=12}.
 *
 * <h2>Why no full-regeneration comparison, unlike {@code
 * BaselineK0SpectrumTest#goldenArtifactMatchesRegeneration}</h2>
 * That precedent's golden artifact regenerates in well under a second, so
 * comparing it field-by-field against a fresh in-test run is cheap. The
 * committed atlas here was generated at {@code ticksObserved=15000} - a
 * real {@link com.chiralbehaviors.inviscid.measure.AuditedRun} over that
 * many ticks measured ~294s (~4m54s) wall (see the Stage 2 generation
 * report) - regenerating it inside surefire on every build would make the
 * suite unusable. Instead this class:
 * <ul>
 * <li>{@link #committedHeaderMatchesTheUserDecidedParameters()} validates
 * the full header against the known-correct parameters (the "regenerate
 * with ContactAtlasGenerator" recovery path if this ever drifts).</li>
 * <li>{@link #committedAtlasAgreesWithLiveContactPredicate()} spot-checks
 * >= 500 SEEDED rows sampled from the COMMITTED file itself (not a fresh
 * generation) against a live {@link ContactPredicate} evaluated at the
 * committed header's own bin centers.</li>
 * <li>{@link #geometricRowsReproduceUnderRegeneration()} specifically
 * re-derives the GEOMETRIC half of the artifact - a direct {@link
 * ContactPredicate} evaluation at bin centers, exactly what {@code
 * ContactAtlasGenerator#sweepGeometricGroundTruth} computes - for a
 * sample, and compares BOTH the {@code contact} verdict and {@code
 * minDistance} bit-for-bit. This is the closest cheap analogue to
 * {@code goldenArtifactMatchesRegeneration}: it exercises the exact same
 * deterministic computation the real generator used, just without paying
 * for the {@code AuditedRun} dynamic half, which this test does not need
 * in order to confirm the geometric ground truth is reproducible.</li>
 * </ul>
 * Both are cheap (well under a second) because neither drives the dynamic
 * {@code AuditedRun} a full regeneration would require.
 *
 * @author halhildebrand
 */
public class CommittedContactAtlasTest {

    private static final String RESOURCE_PATH = "lga/contact-atlas-v1.tsv";

    private static final int     EXPECTED_N_LGA               = 24;
    private static final int     EXPECTED_GEOMETRY_RESOLUTION = ContactAtlasGenerator.GEOMETRY_RESOLUTION;
    private static final double  EXPECTED_RADIUS              = ContactAtlasGenerator.RADIUS;
    private static final Point3i EXPECTED_EXTENT              = ContactAtlasGenerator.DEFAULT_EXTENT;
    private static final long    EXPECTED_SEED                = ContactAtlasGenerator.DEFAULT_SEED;
    private static final int     EXPECTED_TICKS                = 15000;

    private static ContactAtlas ATLAS;

    @BeforeClass
    public static void loadCommittedAtlas() throws IOException {
        ATLAS = readCommitted();
    }

    private static ContactAtlas readCommitted() throws IOException {
        URL resource = CommittedContactAtlasTest.class.getClassLoader()
                                                        .getResource(RESOURCE_PATH);
        if (resource == null) {
            fail("regenerate with ContactAtlasGenerator and review the diff: "
               + "committed atlas src/test/resources/lga/contact-atlas-v1.tsv is missing");
        }
        Path path = Paths.get(resource.getPath());
        return ContactAtlas.read(path);
    }

    private static ContactPredicate newPredicate() {
        return new ContactPredicate(new MemberGeometry(EXPECTED_GEOMETRY_RESOLUTION,
                                                        EXPECTED_RADIUS));
    }

    /**
     * Non-vacuity precondition every test here relies on.
     */
    @Test
    public void committedAtlasIsNonEmpty() {
        assertFalse("regenerate with ContactAtlasGenerator and review the diff: "
                    + "committed atlas has zero rows",
                    ATLAS.rows().isEmpty());
    }

    /**
     * The full header, checked field-by-field against the user-decided
     * generation parameters (bead acceptance criterion: "N_lga chosen
     * with a written justification... recorded in the header").
     * {@code gitCommit} is checked for presence/non-blankness only, not an
     * exact SHA - it legitimately changes on every regeneration.
     */
    @Test
    public void committedHeaderMatchesTheUserDecidedParameters() {
        ContactAtlas.Header header = ATLAS.header();
        String prefix = "regenerate with ContactAtlasGenerator and review the diff: ";

        assertEquals(prefix + "atlasVersion", ContactAtlas.ATLAS_VERSION,
                    header.atlasVersion());
        assertEquals(prefix + "generatorClass",
                    ContactAtlasGenerator.class.getName(),
                    header.generatorClass());
        assertFalse(prefix + "gitCommit must be populated, not blank",
                    header.gitCommit() == null || header.gitCommit().isBlank());
        assertEquals(prefix + "memberRadius", EXPECTED_RADIUS,
                    header.memberRadius(), 0.0);
        assertEquals(prefix + "geometryResolution",
                    EXPECTED_GEOMETRY_RESOLUTION, header.geometryResolution());
        assertEquals(prefix + "cubeEdgeLength",
                    PhiCoordinates.Cubes[0].getEdgeLength(),
                    header.cubeEdgeLength(), 0.0);
        assertEquals(prefix + "phaseResolutionNLga (the USER-DECIDED N_lga)",
                    EXPECTED_N_LGA, header.phaseResolutionNLga());
        assertEquals(prefix + "phiCoordinatesCubeSet", "Cubes[0]",
                    header.phiCoordinatesCubeSet());
        assertEquals(prefix + "extent", EXPECTED_EXTENT, header.extent());
        assertEquals(prefix + "seed", EXPECTED_SEED, header.seed());
        assertEquals(prefix + "ticksObserved", EXPECTED_TICKS,
                    header.ticksObserved());
    }

    /**
     * For >= 500 seeded sample rows drawn from the COMMITTED atlas's own
     * rows (not a freshly-regenerated one, and not the full bin-grid
     * universe - {@link ContactAtlasTest#atlasAgreesWithLiveContactPredicate()}
     * already covers that at the illustrative test N_lga), every row's
     * {@code contact} verdict agrees with a live {@link ContactPredicate}
     * evaluation at the row's own bin centers.
     */
    @Test
    public void committedAtlasAgreesWithLiveContactPredicate() {
        ContactPredicate predicate = newPredicate();
        List<ContactAtlas.Row> rows = ATLAS.rows();
        Random random = new Random(42L);

        int sampleCount = 600;
        for (int i = 0; i < sampleCount; i++) {
            ContactAtlas.Row row = rows.get(random.nextInt(rows.size()));
            float angleA = (float) ContactAtlasGenerator.binCenter(row.phaseBinA(),
                                                                    EXPECTED_N_LGA);
            float angleB = (float) ContactAtlasGenerator.binCenter(row.phaseBinB(),
                                                                    EXPECTED_N_LGA);
            boolean live = predicate.contacts(row.cubeA(), row.memberA(),
                                              angleA, row.cubeB(),
                                              row.memberB(), angleB,
                                              row.direction());
            assertEquals("regenerate with ContactAtlasGenerator and review the diff: "
                         + "sample " + i + " " + row
                         + " disagrees with a live ContactPredicate evaluation",
                         row.contact(), live);
        }
    }

    /**
     * Paired non-vacuity over the committed data: some rows contact, some
     * do not; and - the bead's core {@code observedCount} requirement -
     * some rows were dynamically reached (a real Phase A run actually hit
     * that bin) and some were not (geometrically possible, dynamically
     * unreached, which the bead's own text calls "legal and informative").
     */
    @Test
    public void committedAtlasDistinguishesGeometricFromDynamicReachability() {
        long contactRows = ATLAS.rows().stream().filter(ContactAtlas.Row::contact)
                                 .count();
        long noContactRows = ATLAS.rows().size() - contactRows;
        long observedRows = ATLAS.rows().stream()
                                  .filter(row -> row.observedCount() > 0)
                                  .count();
        long unobservedRows = ATLAS.rows().size() - observedRows;

        String prefix = "regenerate with ContactAtlasGenerator and review the diff: ";
        assertTrue(prefix + "expected at least one contacting row, found none",
                   contactRows > 0);
        assertTrue(prefix + "expected at least one non-contacting row, found none (atlas looks total, not partial)",
                   noContactRows > 0);
        assertTrue(prefix + "expected at least one dynamically-observed row (observedCount>0), found none",
                   observedRows > 0);
        assertTrue(prefix + "expected at least one geometrically-possible-but-dynamically-unreached row (observedCount==0), found none",
                   unobservedRows > 0);
    }

    /**
     * Re-derives the GEOMETRIC half of a sample of committed rows - the
     * exact computation {@code ContactAtlasGenerator#sweepGeometricGroundTruth}
     * performed at generation time (bin-center {@link ContactPredicate}
     * evaluation) - and compares BOTH {@code contact} and {@code
     * minDistance} bit-for-bit (deterministic floating-point recomputation
     * of the same inputs, not a numeric-tolerance comparison). See class
     * Javadoc for why this stands in for a full {@code
     * goldenArtifactMatchesRegeneration}-style comparison without paying
     * for the 15,000-tick dynamic run.
     */
    @Test
    public void geometricRowsReproduceUnderRegeneration() {
        ContactPredicate predicate = newPredicate();
        List<ContactAtlas.Row> rows = ATLAS.rows();
        Random random = new Random(1729L);

        int sampleCount = 150;
        int contactSeen = 0;
        int noContactSeen = 0;
        for (int i = 0; i < sampleCount; i++) {
            ContactAtlas.Row row = rows.get(random.nextInt(rows.size()));
            float angleA = (float) ContactAtlasGenerator.binCenter(row.phaseBinA(),
                                                                    EXPECTED_N_LGA);
            float angleB = (float) ContactAtlasGenerator.binCenter(row.phaseBinB(),
                                                                    EXPECTED_N_LGA);

            boolean regeneratedContact = predicate.contacts(row.cubeA(),
                                                             row.memberA(),
                                                             angleA,
                                                             row.cubeB(),
                                                             row.memberB(),
                                                             angleB,
                                                             row.direction());
            double regeneratedMinDistance = predicate.minDistance(row.cubeA(),
                                                                   row.memberA(),
                                                                   angleA,
                                                                   row.cubeB(),
                                                                   row.memberB(),
                                                                   angleB,
                                                                   row.direction());

            assertEquals("regenerate with ContactAtlasGenerator and review the diff: "
                         + "sample " + i + " " + row + " contact verdict differs on regeneration",
                         row.contact(), regeneratedContact);
            assertEquals("regenerate with ContactAtlasGenerator and review the diff: "
                         + "sample " + i + " " + row + " minDistance differs on regeneration",
                         row.minDistance(), regeneratedMinDistance, 1e-12);

            if (row.contact()) {
                contactSeen++;
            } else {
                noContactSeen++;
            }
        }
        assertTrue("regenerate with ContactAtlasGenerator and review the diff: "
                  + "sample contained no contacting rows - widen the sample",
                  contactSeen > 0);
        assertTrue("regenerate with ContactAtlasGenerator and review the diff: "
                  + "sample contained no non-contacting rows - widen the sample",
                  noContactSeen > 0);
    }
}
