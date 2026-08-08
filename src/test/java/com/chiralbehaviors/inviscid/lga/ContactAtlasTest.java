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
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import javax.vecmath.Point3i;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Behavioral tests for {@link ContactAtlas} / {@link ContactAtlasGenerator}
 * (bead inviscid-0nx.16, A.5). Run at a FIXED test {@code nLga} ({@link
 * #TEST_N_LGA}, 12) - clearly a TEST parameter, not the epic's user-reserved
 * production N_lga decision (which stays in {@code {8, 12, 16, 24}}, chosen
 * downstream from the measurement campaign this bead's generator also
 * drives - see {@code ContactAtlasGenerator}'s own Javadoc).
 *
 * <p>The atlas under test ({@link #ATLAS}) is generated once ({@link
 * BeforeClass}) and shared read-only across every test method - generation
 * itself (geometric sweep + a real {@link
 * com.chiralbehaviors.inviscid.measure.AuditedRun}) is the expensive part
 * and every test here only inspects the result, so there is no reason to
 * regenerate per test.
 *
 * @author halhildebrand
 */
public class ContactAtlasTest {

    private static final double  RADIUS     = ContactAtlasGenerator.RADIUS;
    private static final int     RESOLUTION = ContactAtlasGenerator.GEOMETRY_RESOLUTION;
    private static final int     TEST_N_LGA = 12;
    private static final Point3i EXTENT     = new Point3i(4, 4, 4);
    private static final long    SEED       = 42L;
    private static final int     TICKS      = 300;

    private static ContactAtlas ATLAS;

    private record Key(int direction, int cubeA, int memberA, int cubeB,
                        int memberB, int phaseBinA, int phaseBinB) {
    }

    @BeforeClass
    public static void generateAtlas() {
        ATLAS = ContactAtlasGenerator.generate(TEST_N_LGA, EXTENT, SEED,
                                                TICKS);
    }

    private static ContactPredicate newPredicate() {
        return new ContactPredicate(new MemberGeometry(RESOLUTION, RADIUS));
    }

    private static Map<Key, ContactAtlas.Row> indexByKey(ContactAtlas atlas) {
        Map<Key, ContactAtlas.Row> byKey = new HashMap<>();
        for (ContactAtlas.Row row : atlas.rows()) {
            byKey.put(new Key(row.direction(), row.cubeA(), row.memberA(),
                              row.cubeB(), row.memberB(), row.phaseBinA(),
                              row.phaseBinB()), row);
        }
        return byKey;
    }

    /**
     * Non-vacuity precondition every test here relies on: a real atlas at
     * the test parameters actually has rows to inspect.
     */
    @Test
    public void generatedAtlasIsNonEmpty() {
        assertFalse("expected at least one contacting combination at nLga="
                    + TEST_N_LGA, ATLAS.rows().isEmpty());
    }

    /**
     * Write, read, deep-equals - including the header (bead's test 1).
     */
    @Test
    public void roundTripsThroughSerialization() throws IOException {
        Path path = Files.createTempFile("contact-atlas-roundtrip", ".tsv");
        try {
            ATLAS.write(path);
            ContactAtlas reloaded = ContactAtlas.read(path);

            assertEquals(ATLAS.header(), reloaded.header());
            assertEquals(ATLAS.rows(), reloaded.rows());
        } finally {
            Files.deleteIfExists(path);
        }
    }

    /**
     * THE key test (bead's own words): a stale-parameter atlas is refused
     * loudly, not silently consumed - both when a required header
     * parameter is entirely ABSENT (a hand-corrupted or truncated file)
     * and when every parameter is present but one VALUE disagrees with
     * what the caller expects (a table generated for a different {@code
     * nLga}, extent, seed, etc).
     */
    @Test
    public void rejectsAtlasWithMissingOrMismatchedHeader() throws IOException {
        Path validPath = Files.createTempFile("contact-atlas-valid", ".tsv");
        Path missingPath = Files.createTempFile("contact-atlas-missing-header",
                                                 ".tsv");
        try {
            ATLAS.write(validPath);

            // Missing: drop the "# seed=..." header line entirely.
            List<String> lines = Files.readAllLines(validPath);
            List<String> corrupted = lines.stream()
                                           .filter(line -> !line.startsWith("# seed="))
                                           .toList();
            Files.write(missingPath, corrupted);

            ContactAtlas.HeaderMismatchException missing = assertThrows(ContactAtlas.HeaderMismatchException.class,
                                                                          () -> ContactAtlas.read(missingPath));
            assertTrue("expected the missing key named in the failure: "
                       + missing.getMessage(),
                       missing.getMessage().contains("seed"));

            // Mismatched: every parameter present, but the caller's
            // expectation disagrees on phaseResolutionNLga.
            ContactAtlas.Header actual = ATLAS.header();
            ContactAtlas.Header wrongExpectation = new ContactAtlas.Header(actual.atlasVersion(),
                                                                            actual.generatorClass(),
                                                                            actual.gitCommit(),
                                                                            actual.memberRadius(),
                                                                            actual.geometryResolution(),
                                                                            actual.cubeEdgeLength(),
                                                                            actual.phaseResolutionNLga()
                                                                            + 1,
                                                                            actual.phiCoordinatesCubeSet(),
                                                                            actual.extent(),
                                                                            actual.seed(),
                                                                            actual.ticksObserved());
            ContactAtlas.HeaderMismatchException mismatched = assertThrows(ContactAtlas.HeaderMismatchException.class,
                                                                             () -> ContactAtlas.readValidated(validPath,
                                                                                                               wrongExpectation));
            assertTrue("expected the mismatched key named in the failure: "
                       + mismatched.getMessage(),
                       mismatched.getMessage()
                                 .contains("phaseResolutionNLga"));

            // Positive control: readValidated against the TRUE header
            // does not throw.
            ContactAtlas.readValidated(validPath, actual);
        } finally {
            Files.deleteIfExists(validPath);
            Files.deleteIfExists(missingPath);
        }
    }

    /**
     * bead inviscid-0nx.16.2: {@code geometryResolution} (the {@code
     * MemberGeometry} angle-quantization LUT resolution - a real tunable
     * of {@code ContactAtlasGenerator}'s 6/7-arg {@code generate}
     * overloads that changes which bins geometrically contact) must be
     * part of the staleness contract exactly like every other header
     * parameter - mirrors {@code rejectsAtlasWithMissingOrMismatchedHeader}'s
     * mismatch half, mechanically, for this specific key.
     */
    @Test
    public void rejectsAtlasWithMismatchedGeometryResolution() throws IOException {
        Path path = Files.createTempFile("contact-atlas-geometry-resolution",
                                          ".tsv");
        try {
            ATLAS.write(path);

            ContactAtlas.Header actual = ATLAS.header();
            ContactAtlas.Header wrongExpectation = new ContactAtlas.Header(actual.atlasVersion(),
                                                                            actual.generatorClass(),
                                                                            actual.gitCommit(),
                                                                            actual.memberRadius(),
                                                                            actual.geometryResolution()
                                                                            + 1,
                                                                            actual.cubeEdgeLength(),
                                                                            actual.phaseResolutionNLga(),
                                                                            actual.phiCoordinatesCubeSet(),
                                                                            actual.extent(), actual.seed(),
                                                                            actual.ticksObserved());

            ContactAtlas.HeaderMismatchException mismatched = assertThrows(ContactAtlas.HeaderMismatchException.class,
                                                                             () -> ContactAtlas.readValidated(path,
                                                                                                               wrongExpectation));
            assertTrue("expected the mismatched key named in the failure: "
                       + mismatched.getMessage(),
                       mismatched.getMessage()
                                 .contains("geometryResolution"));

            // Positive control.
            ContactAtlas.readValidated(path, actual);
        } finally {
            Files.deleteIfExists(path);
        }
    }

    /**
     * For >= 500 seeded sample points drawn from the FULL {@code
     * (direction, cubeA, memberA, cubeB, memberB, phaseBinA, phaseBinB)}
     * universe at {@link #TEST_N_LGA} - not merely existing atlas rows -
     * the atlas's verdict (present with {@code contact=true}, present
     * with {@code contact=false}, or absent) always agrees with a live
     * {@link ContactPredicate} evaluation at the same bin centers (bead's
     * test 3).
     */
    @Test
    public void atlasAgreesWithLiveContactPredicate() {
        ContactPredicate predicate = newPredicate();
        Map<Key, ContactAtlas.Row> byKey = indexByKey(ATLAS);
        Random random = new Random(42L);

        int sampleCount = 600;
        int checked = 0;
        for (int i = 0; i < sampleCount; i++) {
            int direction = FccNeighborhood.DIRECTIONS.get(random.nextInt(FccNeighborhood.DIRECTIONS.size()));
            int cubeA = random.nextInt(5);
            int memberA = random.nextInt(6);
            int cubeB = random.nextInt(5);
            int memberB = random.nextInt(6);
            int binA = random.nextInt(TEST_N_LGA);
            int binB = random.nextInt(TEST_N_LGA);

            float angleA = (float) ContactAtlasGenerator.binCenter(binA,
                                                                    TEST_N_LGA);
            float angleB = (float) ContactAtlasGenerator.binCenter(binB,
                                                                    TEST_N_LGA);
            boolean live = predicate.contacts(cubeA, memberA, angleA, cubeB,
                                              memberB, angleB, direction);

            Key key = new Key(direction, cubeA, memberA, cubeB, memberB,
                              binA, binB);
            ContactAtlas.Row row = byKey.get(key);
            if (row == null) {
                assertFalse("sample " + i + " " + key
                            + ": absent from atlas but live predicate says contact",
                            live);
            } else {
                assertEquals("sample " + i + " " + key
                             + ": atlas row disagrees with live predicate",
                             live, row.contact());
            }
            checked++;
        }
        assertTrue("expected >= 500 checked samples, got " + checked,
                   checked >= 500);
    }

    /**
     * Non-vacuity, paired: some sampled combinations DO contact, some do
     * NOT - neither degenerate extreme (bead's test 4).
     */
    @Test
    public void atlasIsNeitherEmptyNorTotal() {
        assertFalse("atlas must not be empty", ATLAS.rows().isEmpty());

        long totalUniverse = (long) FccNeighborhood.DIRECTIONS.size() * 5 * 6
                              * 5 * 6 * TEST_N_LGA * TEST_N_LGA;
        assertTrue("atlas rows (" + ATLAS.rows().size()
                   + ") must be a proper subset of the full bin-grid universe ("
                   + totalUniverse + ") - not every combination contacts",
                   ATLAS.rows().size() < totalUniverse);
    }

    /**
     * {@code contacts(A in C, B in C+d)} implies {@code contacts(B in
     * C+d, A in C, opposite(d))} - the same property {@code
     * ContactPredicateTest} verifies directly on {@link ContactPredicate},
     * checked here at the ATLAS level: since the geometric sweep
     * enumerates every one of the 12 {@link FccNeighborhood#DIRECTIONS}
     * (not just the canonical 6), a contacting row's direction-reversed
     * mirror must independently appear in the atlas too, not merely hold
     * true if re-evaluated (bead's test 5).
     */
    @Test
    public void atlasIsSymmetricUnderDirectionReversal() {
        ContactPredicate predicate = newPredicate();
        Map<Key, ContactAtlas.Row> byKey = indexByKey(ATLAS);

        List<ContactAtlas.Row> contactingRows = ATLAS.rows()
                                                       .stream()
                                                       .filter(ContactAtlas.Row::contact)
                                                       .toList();
        assertFalse("need at least one contacting row to check symmetry on",
                    contactingRows.isEmpty());

        Random random = new Random(42L);
        int sampleSize = Math.min(100, contactingRows.size());
        int checked = 0;
        for (int i = 0; i < sampleSize; i++) {
            ContactAtlas.Row row = contactingRows.get(random.nextInt(contactingRows.size()));
            int oppositeDirection = FccNeighborhood.opposite(row.direction());

            float angleA = (float) ContactAtlasGenerator.binCenter(row.phaseBinA(),
                                                                    TEST_N_LGA);
            float angleB = (float) ContactAtlasGenerator.binCenter(row.phaseBinB(),
                                                                    TEST_N_LGA);
            assertTrue("ContactPredicate itself must be symmetric under direction reversal for "
                       + row,
                       predicate.contacts(row.cubeB(), row.memberB(), angleB,
                                         row.cubeA(), row.memberA(), angleA,
                                         oppositeDirection));

            Key mirrorKey = new Key(oppositeDirection, row.cubeB(),
                                    row.memberB(), row.cubeA(), row.memberA(),
                                    row.phaseBinB(), row.phaseBinA());
            ContactAtlas.Row mirrorRow = byKey.get(mirrorKey);
            assertTrue("expected the direction-reversed mirror of " + row
                       + " to also appear in the atlas as a contacting row",
                       mirrorRow != null && mirrorRow.contact());
            checked++;
        }
        assertTrue(checked > 0);
    }

    /**
     * The header's {@code phaseResolutionNLga} matches the test's N_lga,
     * and every row's bin indices are honoured (in {@code [0, nLga)}
     * (bead's test 6).
     */
    @Test
    public void phaseResolutionIsRecordedAndHonoured() {
        assertEquals(TEST_N_LGA, ATLAS.header().phaseResolutionNLga());
        for (ContactAtlas.Row row : ATLAS.rows()) {
            assertTrue("phaseBinA out of range: " + row,
                       row.phaseBinA() >= 0 && row.phaseBinA() < TEST_N_LGA);
            assertTrue("phaseBinB out of range: " + row,
                       row.phaseBinB() >= 0 && row.phaseBinB() < TEST_N_LGA);
        }
    }
}
