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

package com.chiralbehaviors.inviscid.automaton.lga;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.junit.BeforeClass;
import org.junit.Test;

import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.sun.management.ThreadMXBean;

/**
 * Failing-tests-first (TDD) coverage for {@link ContactTable} (bead
 * inviscid-0nx.19, C.2): the runtime contact-predicate table transcribed
 * from the committed Phase A atlas ({@code contact-atlas-v2.tsv}) plus
 * {@link MemberGeometry} geometry - see bead inviscid-gyt's ANY-OVERLAP
 * transcription-semantics decision (table fires iff a row's {@code
 * overlapFraction > 0}), which supersedes the bead's original
 * bin-center-only wording.
 *
 * <p>The committed atlas is loaded once ({@link BeforeClass}) and shared
 * read-only across test methods, mirroring {@link CommittedContactAtlasTest}'s
 * own idiom.
 *
 * @author halhildebrand
 */
public class ContactTableTest {

    private static final String RESOURCE_PATH = "lga/contact-atlas-v2.tsv";

    /**
     * The four header fields {@link #rejectsStaleAtlas()} mutates one at a
     * time - exactly the set the bead's own text names (memberRadius,
     * N_lga, geometryResolution, cubeEdgeLength).
     */
    private static final List<String> MUTATED_FIELDS = List.of("memberRadius",
                                                                 "geometryResolution",
                                                                 "phaseResolutionNLga",
                                                                 "cubeEdgeLength");

    private static ContactAtlas ATLAS;
    private static Path         ATLAS_PATH;
    private static ContactTable TABLE;

    @BeforeClass
    public static void loadCommittedAtlasAndTable() throws IOException {
        URL resource = ContactTableTest.class.getClassLoader()
                                              .getResource(RESOURCE_PATH);
        if (resource == null) {
            fail("committed atlas src/test/resources/lga/contact-atlas-v2.tsv is missing");
        }
        ATLAS_PATH = Paths.get(resource.getPath());
        ATLAS = ContactAtlas.read(ATLAS_PATH);
        TABLE = ContactTable.load(ATLAS_PATH, ATLAS.header());
    }

    /**
     * The transcription-fidelity test and the point of the bead
     * (corrected per bead inviscid-gyt): for EVERY one of the committed
     * atlas's 4096 rows, the table's verdict equals {@code overlapFraction
     * > 0} - exhaustive, not a sample, and therefore exceeds the bead
     * acceptance criterion's ">= 5000 sampled combinations" floor together
     * with {@link #packedIndexingRoundTrips()}'s exhaustive full-domain
     * sweep (6,220,800 combinations at the committed N_lga=24).
     */
    @Test
    public void tableTranscribesAtlasOverlapExactly() {
        int mismatches = 0;
        for (ContactAtlas.Row row : ATLAS.rows()) {
            boolean expected = row.overlapFraction() > 0.0;
            boolean actual = TABLE.contacts(row.direction(), row.cubeA(),
                                            row.memberA(), row.phaseBinA(),
                                            row.cubeB(), row.memberB(),
                                            row.phaseBinB());
            if (expected != actual) {
                mismatches++;
            }
        }
        assertEquals("table verdict disagreed with (overlapFraction>0) for "
                     + mismatches + " of " + ATLAS.rows().size()
                     + " committed atlas rows - table transcription is not exact",
                     0, mismatches);
    }

    /**
     * DISCRIMINATION test (code-review Important finding 1, T1 scratch
     * 6325adaa): the committed atlas has ZERO rows with {@code
     * overlapFraction <= 0} ({@code
     * ContactAtlasGenerator.sweepComboOverlap}'s {@code if (hits == 0)
     * continue;} guard means every committed row already satisfies {@code
     * overlapFraction > 0}), so {@link #tableTranscribesAtlasOverlapExactly()}
     * alone cannot distinguish "transcribe iff {@code overlapFraction >
     * 0}" from "transcribe every atlas row unconditionally" - a regression
     * that dropped {@link ContactTable#of(ContactAtlas)}'s {@code
     * row.overlapFraction() > 0.0} guard would stay green against
     * committed data alone. This test builds a SYNTHETIC in-memory atlas
     * (same committed header, two fabricated rows at distinct cells - one
     * {@code overlapFraction=0.0}, one {@code >0}) that CAN discriminate.
     * <p>
     * MUTATION-VERIFIED, not merely asserted: the guard in {@code
     * ContactTable.of(ContactAtlas)} was temporarily weakened to
     * unconditional (every row transcribed regardless of {@code
     * overlapFraction}), {@code mvn test -Dtest=ContactTableTest} was
     * re-run, and this test - ONLY this test, every other {@code
     * ContactTableTest} method stayed green against the weakened guard,
     * confirming committed-data-only coverage cannot see this class of
     * regression - failed with the zero-row assertion (verbatim red
     * output recorded in the implementing agent's final report / T1
     * scratch). The guard was then restored exactly and this suite
     * re-confirmed green.
     */
    @Test
    public void tableDiscriminatesZeroFromPositiveOverlapFraction() {
        ContactAtlas.Header header = ATLAS.header();
        ContactAtlas.Row zeroRow = new ContactAtlas.Row(1, 0, 0, 0, 1, 0, 0,
                                                         false, 0.0, 1.0, 0L);
        ContactAtlas.Row positiveRow = new ContactAtlas.Row(1, 0, 0, 0, 2, 0,
                                                             0, true, 0.5,
                                                             0.001, 5L);
        ContactAtlas fixture = new ContactAtlas(header,
                                                List.of(zeroRow, positiveRow));

        ContactTable table = ContactTable.of(fixture);

        assertFalse("row with overlapFraction=0.0 must NOT be transcribed as contacting",
                   table.contacts(1, 0, 0, 0, 0, 1, 0));
        assertTrue("row with overlapFraction>0 must be transcribed as contacting",
                  table.contacts(1, 0, 0, 0, 0, 2, 0));
    }

    /**
     * Code-review Important finding 2 (T1 scratch 6325adaa):
     * {@link ContactTable#memoryFootprintBytes()} and {@link
     * ContactTable#firedCount()} exist specifically to satisfy the bead's
     * "memory footprint recorded" acceptance criterion, but had zero
     * regression-proof test coverage (the 777,600-byte figure lived only
     * in a T1 scratch note). {@code firedCount}'s expected value is
     * computed by an INDEPENDENT oracle here (a {@link Set} of distinct
     * fired 7-tuple keys, not by calling into {@link
     * ContactTable#of(ContactAtlas)}'s own bit-counting logic) rather than
     * trusted from the atlas's raw row count, so a future atlas with
     * duplicate or additional zero-overlap rows would not silently
     * invalidate this assertion.
     */
    @Test
    public void memoryFootprintAndFiredCountArePinned() {
        assertEquals("measured memory footprint must match the committed N_lga=24 design calculation (97,200 longs * 8 bytes)",
                    777_600L, TABLE.memoryFootprintBytes());

        Set<List<Integer>> firedKeys = new HashSet<>();
        for (ContactAtlas.Row row : ATLAS.rows()) {
            if (row.overlapFraction() > 0.0) {
                firedKeys.add(List.of(row.direction(), row.cubeA(),
                                      row.memberA(), row.cubeB(),
                                      row.memberB(), row.phaseBinA(),
                                      row.phaseBinB()));
            }
        }
        assertEquals("firedCount must equal the independently-counted distinct fired-cell keys from the committed atlas",
                    (long) firedKeys.size(), TABLE.firedCount());
    }

    /**
     * Substantive-critic Significant finding (final review round, T1
     * scratch for inviscid-0nx.19): direction-reversal symmetry - a
     * physical contact between {@code (cubeA, memberA)} in cell {@code C}
     * and {@code (cubeB, memberB)} in {@code C}'s {@code direction}
     * neighbor is the SAME physical event as {@code (cubeB, memberB)}
     * seeing {@code (cubeA, memberA)} via {@code opposite(direction)} -
     * checked here exhaustively over all 4096 committed atlas rows. This
     * is a PIN, not a fix: the critic independently verified the
     * invariant already holds bit-exactly (0 mismatches) against the
     * committed data - {@code ContactAtlasGenerator.sweepOverlapAndCenter}
     * sweeps all 12 directions per combo independently (not merely
     * mirrored), so both orientations of a physical contact are
     * genuinely, independently geometrically computed, and this test
     * confirms {@link ContactTable} preserves that symmetry through
     * transcription rather than introducing a directional bug. {@code
     * opposite(d)} uses {@link FccNeighborhood#opposite}'s convention
     * (raw signed-value negation, {@code -direction}) - the same
     * convention {@code ContactAtlasGenerator
     * .mirrorNegativeDirectionObservedCounts} and {@link
     * ContactTable#directionIndex} both key off.
     */
    @Test
    public void symmetryHoldsAcrossAllAtlasRows() {
        int mismatches = 0;
        for (ContactAtlas.Row row : ATLAS.rows()) {
            boolean forward = TABLE.contacts(row.direction(), row.cubeA(),
                                             row.memberA(), row.phaseBinA(),
                                             row.cubeB(), row.memberB(),
                                             row.phaseBinB());
            int opposite = FccNeighborhood.opposite(row.direction());
            boolean reverse = TABLE.contacts(opposite, row.cubeB(),
                                             row.memberB(), row.phaseBinB(),
                                             row.cubeA(), row.memberA(),
                                             row.phaseBinA());
            if (forward != reverse) {
                mismatches++;
            }
        }
        assertEquals("expected direction-reversal symmetry (forward cell and its "
                     + "opposite-direction/A-B-swapped mirror must agree) across every "
                     + "committed atlas row, found " + mismatches + " mismatch(es) of "
                     + ATLAS.rows().size(),
                     0, mismatches);
    }

    /**
     * DEDUP GUARD test (final review round, T1 scratch for
     * inviscid-0nx.19): the committed atlas has no rows sharing the same
     * 7-tuple key (independently confirmed: 4096 rows, 4096 distinct
     * keys), so no committed-data test can discriminate {@link
     * ContactTable#of(ContactAtlas)}'s dedup guard ({@code if
     * (!getBit(bits, idx))}) from an unconditional {@code firedCount++}.
     * This test builds a SYNTHETIC in-memory atlas with two rows sharing
     * the IDENTICAL 7-tuple key (both {@code overlapFraction > 0}) and
     * asserts the cell is counted ONCE, not twice.
     * <p>
     * MUTATION-VERIFIED, not merely asserted: the dedup guard in {@code
     * ContactTable.of(ContactAtlas)} was temporarily replaced with an
     * unconditional {@code firedCount++} (removing the {@code
     * !getBit(bits, idx)} check), {@code mvn test -Dtest=ContactTableTest}
     * was re-run, and this test - ONLY this test; {@link
     * #memoryFootprintAndFiredCountArePinned()} stayed green because the
     * committed atlas's 4096 keys are already distinct, so the mutation is
     * invisible there too - failed with the duplicate-key assertion
     * (verbatim red output recorded in the implementing agent's final
     * report / T1 scratch). The guard was then restored exactly (diff-
     * verified byte-identical to pre-mutation) and this suite re-confirmed
     * green.
     */
    @Test
    public void firedCountDeduplicatesRepeatedKeys() {
        ContactAtlas.Header header = ATLAS.header();
        ContactAtlas.Row first = new ContactAtlas.Row(1, 0, 0, 0, 3, 0, 0,
                                                       true, 0.3, 0.01, 2L);
        ContactAtlas.Row duplicateKey = new ContactAtlas.Row(1, 0, 0, 0, 3,
                                                              0, 0, true, 0.9,
                                                              0.02, 7L);
        ContactAtlas fixture = new ContactAtlas(header,
                                                List.of(first, duplicateKey));

        ContactTable table = ContactTable.of(fixture);

        assertEquals("two atlas rows sharing the same 7-tuple key must be counted as ONE fired cell, not two",
                    1L, table.firedCount());
        assertTrue("the shared cell must still be transcribed as contacting",
                  table.contacts(1, 0, 0, 0, 0, 3, 0));
    }

    /**
     * One-sided live check (correction #2 on the relay, superseding the
     * bead's original "table verdict == ContactPredicate at the bin
     * centre" wording): bin-centre contact, evaluated LIVE against {@link
     * ContactPredicate} (not read back from the atlas's own {@code
     * contact} column), is always a subset of ANY-OVERLAP - the bin
     * centre is itself one of the fine sweep's own samples (see {@code
     * ContactAtlasGenerator.sweepOverlapAndCenter}'s proof), so bin-centre
     * contact MUST imply table contact. The converse does not hold (that
     * is the entire point of the any-overlap decision) and is not
     * asserted here.
     */
    @Test
    public void binCentreContactImpliesTableContact() {
        ContactPredicate predicate = new ContactPredicate(new MemberGeometry(ATLAS.header()
                                                                                    .geometryResolution(),
                                                                              ATLAS.header()
                                                                                    .memberRadius()));
        int nLga = ATLAS.header().phaseResolutionNLga();
        int checked = 0;
        int centreContactsSeen = 0;
        for (ContactAtlas.Row row : ATLAS.rows()) {
            float angleA = (float) ContactAtlasGenerator.binCenter(row.phaseBinA(),
                                                                     nLga);
            float angleB = (float) ContactAtlasGenerator.binCenter(row.phaseBinB(),
                                                                     nLga);
            boolean liveCentreContact = predicate.contacts(row.cubeA(),
                                                            row.memberA(),
                                                            angleA,
                                                            row.cubeB(),
                                                            row.memberB(),
                                                            angleB,
                                                            row.direction());
            checked++;
            if (liveCentreContact) {
                centreContactsSeen++;
                boolean tableContact = TABLE.contacts(row.direction(),
                                                       row.cubeA(),
                                                       row.memberA(),
                                                       row.phaseBinA(),
                                                       row.cubeB(),
                                                       row.memberB(),
                                                       row.phaseBinB());
                assertTrue("bin-centre contact (live ContactPredicate) at " + row
                          + " must imply table contact - bin-centre is a special case of any-overlap",
                          tableContact);
            }
        }
        assertTrue("expected at least one checked row", checked > 0);
        assertTrue("expected at least one live bin-centre contact to actually exercise the implication - test would be vacuous otherwise",
                  centreContactsSeen > 0);
    }

    /**
     * The .19-side defense (correction #3): {@link ContactTable#load}
     * REFUSES to load an atlas whose on-disk header disagrees with the
     * {@code header} a caller expects (the same header a paired geometry
     * consumer - {@link MemberGeometry}, {@link FineStepContactTable} -
     * would be built from) - for each of the four
     * fields the bead's own text names, a copy of the committed atlas with
     * that single header field mutated must be refused, and the failure
     * message must name the offending field.
     */
    @Test
    public void rejectsStaleAtlas() throws IOException {
        for (String field : MUTATED_FIELDS) {
            Path tmp = writeAtlasWithMutatedHeaderField(field);
            try {
                ContactAtlas.HeaderMismatchException thrown = null;
                try {
                    ContactTable.load(tmp, ATLAS.header());
                    fail("expected refusal to load an atlas with a stale " + field);
                } catch (ContactAtlas.HeaderMismatchException e) {
                    thrown = e;
                }
                assertTrue("exception message for stale " + field
                          + " must name the mismatched field: " + thrown.getMessage(),
                          thrown.getMessage().contains(field));
            } finally {
                Files.deleteIfExists(tmp);
            }
        }
    }

    /**
     * No allocation in the hot path: measured directly via {@code
     * com.sun.management.ThreadMXBean}'s per-thread allocation counter
     * (JIT-warmed first, so the measurement window is not dominated by
     * interpreter/compilation noise unrelated to {@code contacts()}
     * itself). "Constant time" is a structural property of the
     * implementation (one fixed-arithmetic {@code packIndex} computation
     * plus one array-index bit test - no loop, no data-dependent branch
     * beyond bounds checks) rather than something this test times
     * empirically; recorded here, not (re)proven by a timing assertion,
     * per the bead's own "or by inspection recorded in notes" allowance.
     */
    private static final int[] DIRECTIONS = { 1, 2, 3, 4, 5, 6, -1, -2, -3,
                                              -4, -5, -6 };

    @Test
    public void lookupIsConstantTimeAndAllocationFree() {
        ThreadMXBean threadBean = (ThreadMXBean) ManagementFactory.getThreadMXBean();
        assertTrue("JVM must support per-thread allocation counting for this test",
                  threadBean.isThreadAllocatedMemorySupported());
        threadBean.setThreadAllocatedMemoryEnabled(true);
        long threadId = Thread.currentThread().threadId();

        Random random = new Random(7L);
        int nLga = TABLE.nLga();

        // Warm up the JIT well past typical C2 compilation thresholds
        // before measuring, so the measured allocation delta reflects the
        // steady-state compiled contacts() path, not startup noise.
        boolean sink = false;
        for (int i = 0; i < 200_000; i++) {
            sink ^= randomLookup(random, nLga);
        }

        long before = threadBean.getThreadAllocatedBytes(threadId);
        int iterations = 1_000_000;
        for (int i = 0; i < iterations; i++) {
            sink ^= randomLookup(random, nLga);
        }
        long after = threadBean.getThreadAllocatedBytes(threadId);

        // Prevent the JIT from eliding the loop as dead code.
        assertTrue("sink sentinel", sink || !sink);

        assertEquals("contacts() allocated " + (after - before) + " byte(s) over "
                     + iterations + " calls - the hot path must be allocation-free",
                     0L, after - before);
    }

    /**
     * NOTE: {@link #DIRECTIONS} is a hoisted, class-level constant -
     * deliberately NOT a local array literal - so this harness method's
     * own random-argument selection allocates nothing either; an
     * in-loop {@code int[]} literal here would corrupt the measurement in
     * {@link #lookupIsConstantTimeAndAllocationFree()} by attributing the
     * TEST's own per-call array allocation to {@link
     * ContactTable#contacts}, which is exactly what the first (red) run
     * of this test caught.
     */
    private boolean randomLookup(Random random, int nLga) {
        int direction = DIRECTIONS[random.nextInt(DIRECTIONS.length)];
        int cubeA = random.nextInt(PhiCoordinates.Cubes.length);
        int memberA = random.nextInt(6);
        int cubeB = random.nextInt(PhiCoordinates.Cubes.length);
        int memberB = random.nextInt(6);
        int binA = random.nextInt(nLga);
        int binB = random.nextInt(nLga);
        return TABLE.contacts(direction, cubeA, memberA, binA, cubeB, memberB,
                              binB);
    }

    /**
     * Non-vacuity, mirroring inviscid-A.5 test 4: the table has at least
     * one {@code true} cell (every committed atlas row is one, by {@link
     * #tableTranscribesAtlasOverlapExactly()}) AND at least one {@code
     * false} cell drawn from OUTSIDE the atlas's row set - not merely
     * inferred from {@code firedCount < domainSize}, but directly probed,
     * so this test would fail loudly if {@code contacts()} were
     * hard-wired to always return {@code true}.
     */
    @Test
    public void tableIsNeitherEmptyNorTotal() {
        long fired = 0;
        for (ContactAtlas.Row row : ATLAS.rows()) {
            if (TABLE.contacts(row.direction(), row.cubeA(), row.memberA(),
                               row.phaseBinA(), row.cubeB(), row.memberB(),
                               row.phaseBinB())) {
                fired++;
            }
        }
        assertTrue("expected at least one true cell in the table (every committed atlas row)",
                  fired > 0);
        assertTrue("table must not be total: fired count (" + fired
                  + ") should be a small fraction of the full domain ("
                  + TABLE.domainSize() + ")",
                  fired < TABLE.domainSize());

        Set<List<Integer>> atlasKeys = new HashSet<>();
        for (ContactAtlas.Row row : ATLAS.rows()) {
            atlasKeys.add(List.of(row.direction(), row.cubeA(), row.memberA(),
                                  row.cubeB(), row.memberB(), row.phaseBinA(),
                                  row.phaseBinB()));
        }

        int nLga = TABLE.nLga();
        Random random = new Random(99L);
        int falseSamplesFound = 0;
        int attempts = 0;
        while (falseSamplesFound < 20 && attempts < 1_000_000) {
            attempts++;
            int direction = DIRECTIONS[random.nextInt(DIRECTIONS.length)];
            int cubeA = random.nextInt(PhiCoordinates.Cubes.length);
            int memberA = random.nextInt(6);
            int cubeB = random.nextInt(PhiCoordinates.Cubes.length);
            int memberB = random.nextInt(6);
            int binA = random.nextInt(nLga);
            int binB = random.nextInt(nLga);
            List<Integer> key = List.of(direction, cubeA, memberA, cubeB,
                                        memberB, binA, binB);
            if (atlasKeys.contains(key)) {
                continue;
            }
            boolean tableContact = TABLE.contacts(direction, cubeA, memberA,
                                                  binA, cubeB, memberB, binB);
            assertFalse("combo " + key
                       + " is not an atlas row (never fired) but table.contacts() returned true",
                       tableContact);
            falseSamplesFound++;
        }
        assertTrue("expected to find at least one non-atlas combo to probe (found "
                  + falseSamplesFound + " in " + attempts + " attempts)",
                  falseSamplesFound > 0);
    }

    /**
     * Bead acceptance criterion: the (direction, cubeA, memberA, cubeB,
     * memberB, binA, binB) -> packed bit index mapping is injective over
     * its FULL domain (exhaustive, not sampled - a collision here would
     * corrupt the collision table silently, per the bead's own text).
     * Also confirms the mapping is a bijection onto {@code [0,
     * domainSize)}: every domain slot is visited exactly once.
     */
    @Test
    public void packedIndexingRoundTrips() {
        int nLga = ATLAS.header().phaseResolutionNLga();
        int cubeCount = PhiCoordinates.Cubes.length;
        int memberCount = 6;

        long domainSize = (long) DIRECTIONS.length * cubeCount * memberCount
                          * cubeCount * memberCount * nLga * nLga;
        assertTrue("domain size must fit an int-addressable BitSet for this exhaustive test",
                  domainSize <= Integer.MAX_VALUE);
        assertEquals("sanity: domain size formula must match ContactTable's own",
                    TABLE.domainSize(), domainSize);

        BitSet seen = new BitSet((int) domainSize);
        long collisions = 0;
        long checked = 0;
        for (int direction : DIRECTIONS) {
            for (int cubeA = 0; cubeA < cubeCount; cubeA++) {
                for (int memberA = 0; memberA < memberCount; memberA++) {
                    for (int binA = 0; binA < nLga; binA++) {
                        for (int cubeB = 0; cubeB < cubeCount; cubeB++) {
                            for (int memberB = 0; memberB < memberCount; memberB++) {
                                for (int binB = 0; binB < nLga; binB++) {
                                    long idx = ContactTable.packIndex(direction,
                                                                      cubeA,
                                                                      memberA,
                                                                      binA,
                                                                      cubeB,
                                                                      memberB,
                                                                      binB,
                                                                      nLga);
                                    assertTrue("index out of domain range: "
                                              + idx,
                                              idx >= 0 && idx < domainSize);
                                    int i = (int) idx;
                                    if (seen.get(i)) {
                                        collisions++;
                                    } else {
                                        seen.set(i);
                                    }
                                    checked++;
                                }
                            }
                        }
                    }
                }
            }
        }

        assertEquals("expected the packed index to be injective over the full 7-tuple domain, found "
                     + collisions + " collision(s) over " + checked
                     + " checked combinations",
                     0, collisions);
        assertEquals("expected every domain slot to be visited exactly once (checked count must equal domain size)",
                    domainSize, checked);
        assertEquals("expected the packed index to be a bijection onto [0, domainSize)",
                    domainSize, seen.cardinality());
    }

    private static Path writeAtlasWithMutatedHeaderField(String field) throws IOException {
        List<String> lines = Files.readAllLines(ATLAS_PATH);
        List<String> mutated = new ArrayList<>(lines.size());
        String prefix = "# " + field + "=";
        for (String line : lines) {
            if (line.startsWith(prefix)) {
                mutated.add(prefix + mutatedValueFor(field));
            } else {
                mutated.add(line);
            }
        }
        Path tmp = Files.createTempFile("contact-atlas-stale-" + field, ".tsv");
        Files.write(tmp, mutated);
        return tmp;
    }

    private static String mutatedValueFor(String field) {
        return switch (field) {
        case "memberRadius" -> "0.999";
        case "geometryResolution" -> "8";
        case "phaseResolutionNLga" -> "8";
        case "cubeEdgeLength" -> "0.001";
        default -> throw new IllegalStateException("unknown mutated field: "
                                                    + field);
        };
    }
}
