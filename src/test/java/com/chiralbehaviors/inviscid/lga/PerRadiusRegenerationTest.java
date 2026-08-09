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
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.UUID;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.vecmath.Point3i;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Conformance tests for {@link PerRadiusRegeneration} (bead
 * inviscid-0nx.30, E.3), the per-{@code r} artifact regeneration harness
 * that implements {@code design-seeding-radius.md} §D-B ("{@code r} is an
 * explicit parameter of every run and artifact header, never a hardcoded
 * constant in new code") and §Engineering-constraints ("contact-atlas-v2
 * .tsv and all committed artifacts are PINNED; new {@code (r, rho)}
 * physics produces atlas v3 + NEW versioned campaign artifacts alongside
 * old (pin-capture rule) ... Never overwrite v2").
 *
 * <h2>Why this suite runs at a reduced geometry resolution</h2>
 * A regeneration at the campaign's real {@link
 * ContactAtlasGenerator#GEOMETRY_RESOLUTION} of 360 pays a full {@link
 * ContactComboCache#sweepExhaustively} at the new radius (the checked-in
 * cache's header declares {@code geometryResolution=360, memberRadius=
 * 0.015}, so ANY other {@code (resolution, radius)} pair misses it - see
 * that class's cache-staleness contract) - several minutes of wall time,
 * which is a manual-run cost, not a {@code mvn test} cost. These tests
 * therefore drive the harness at a small resolution and tick count where
 * the same code paths run in seconds. The radius, the atlas version, and
 * the output-path policy - everything this suite actually asserts - are
 * resolution-independent. The full-fidelity {@code r=0.05} run at
 * resolution 360 is {@link PerRadiusRegeneration#main}'s job.
 *
 * @author halhildebrand
 */
public class PerRadiusRegenerationTest {

    /**
     * Deliberately NOT {@link LgaTestGeometry#BASELINE_RADIUS}: the whole
     * point of this suite is to prove the harness stamps a radius that is
     * not the pinned baseline. Wiring this to the shared constant would
     * make every assertion below vacuous the moment the baseline moved.
     */
    private static final double SWEEP_RADIUS = 0.05;

    /** Divisible by 8 ({@code MemberGeometry}'s constructor contract). */
    private static final int  TEST_RESOLUTION = 24;
    private static final int  TEST_NLGA       = 6;
    private static final int  TEST_TICKS      = 50;
    private static final long TEST_SEED       = 42L;

    /**
     * Every committed artifact the bead's acceptance criteria declare
     * PINNED, by FULL path (the plan audit's finding 3: bare filenames
     * silently miss files or match the wrong ones).
     */
    private static final List<String> PINNED_ARTIFACTS = List.of("src/test/resources/lga/contact-atlas-v2.tsv",
                                                                  "src/test/resources/lga/anisotropy-report-phaseA.tsv",
                                                                  "src/test/resources/lga/measurement-report-phaseC.tsv",
                                                                  "src/test/resources/lga/baseline-k0-spectrum.tsv",
                                                                  "src/main/resources/lga/discovered-combos-cache.tsv");

    /**
     * The trees those artifacts live in. Snapshotting the TREES rather
     * than the five names is what lets pollution - a new file written into
     * a pinned directory - be detected at all.
     */
    private static final List<String> PINNED_TREES = List.of("src/test/resources/lga",
                                                              "src/main/resources/lga");

    /**
     * Snapshot value for a tree entry that has no bytes to hash - a
     * directory, in practice. Present so directories PARTICIPATE in the
     * byte-identity bracket rather than being filtered out of it.
     */
    private static final String NOT_A_REGULAR_FILE = "<not-a-regular-file>";

    private static Map<String, String> PINNED_TREES_ON_ENTRY;

    /**
     * THE PIN GUARD, bracketing the WHOLE class rather than one test.
     *
     * <h2>Two mutation-proved defects this shape exists to avoid</h2>
     * <ol>
     * <li>Hashing five artifacts BY NAME cannot see a file ADDED to a
     * pinned directory, which is exactly what pollution looks like. With
     * the refusal guard neutered, a regeneration wrote four new files into
     * the pinned trees while a fixed-name hash set stayed green. So the
     * snapshot is the full recursive LISTING plus hashes of both trees:
     * additions, deletions and mutations all show. DIRECTORIES included -
     * an earlier spelling filtered on {@code isRegularFile}, which made
     * the bracket blind to an added EMPTY directory, which is precisely
     * the pollution the un-normalized path guard actually produced (see
     * {@link #refusesToWriteThroughADotDotTailThatEscapesTheTargetDirectory}).</li>
     * <li>Bracketing a SINGLE test is order-dependent, and JUnit 4's
     * {@code MethodSorters.DEFAULT} is hash-ordered. Proved, not assumed:
     * with a stray write into {@code src/test/resources/lga} injected into
     * {@code regenerate}, an in-test before/after pair PASSED, because
     * another test in this class had already called {@code regenerate} and
     * created the stray file before the "before" snapshot was taken. A
     * guard that only works when it happens to run first is not a
     * guard.</li>
     * </ol>
     * Snapshotting in {@code @BeforeClass} and comparing in {@code
     * @AfterClass} covers EVERY regeneration any test here performs, in any
     * order.
     */
    @BeforeClass
    public static void snapshotPinnedTreesBeforeAnyRegeneration() throws IOException {
        PINNED_TREES_ON_ENTRY = snapshotPinnedTrees();
    }

    /**
     * @see #snapshotPinnedTreesBeforeAnyRegeneration
     */
    @AfterClass
    public static void noTestInThisClassTouchedAPinnedArtifact() throws IOException {
        assertEquals("running this suite must leave both pinned resource trees byte-identical, with no file added, changed or removed",
                     PINNED_TREES_ON_ENTRY, snapshotPinnedTrees());
    }

    /**
     * The headline acceptance criterion: a v3 atlas regenerated at {@code
     * r != 0.015} carries that radius and {@code atlasVersion=3} in its
     * header.
     *
     * <p>
     * The pin-safety half of the old version of this test moved to
     * {@link #snapshotPinnedTreesBeforeAnyRegeneration} / {@link
     * #noTestInThisClassTouchedAPinnedArtifact}, which are order-
     * independent; see there for why the in-test spelling was not. What
     * stays here is the ACTIVE half: the refusal is a mechanism, so drive
     * it at the most dangerous destination there is rather than infer its
     * existence from the default path happening to land elsewhere.
     */
    @Test
    public void regeneratesAV3AtlasAtANonBaselineRadiusWithoutTouchingPinnedArtifacts() throws IOException {
        try {
            PerRadiusRegeneration.regenerate(smallSweepParameters(),
                                              Path.of("src/test/resources/lga"));
            fail("the harness must REFUSE to regenerate into the pinned resource tree, not merely decline to by default");
        } catch (IllegalArgumentException expected) {
            assertTrue("refusal must name the offending path, was: "
                       + expected.getMessage(),
                       expected.getMessage()
                                .contains("src/test/resources/lga"));
        }

        PerRadiusRegeneration.Result result = PerRadiusRegeneration.regenerate(smallSweepParameters());

        Map<String, String> atlasHeader = readHeader(result.atlasPath());
        assertEquals("regenerated atlas must declare atlasVersion=3", "3",
                     atlasHeader.get("atlasVersion"));
        assertEquals("memberRadius is the atlas's physics-identity field and must carry the requested radius verbatim",
                     Double.toString(SWEEP_RADIUS),
                     atlasHeader.get("memberRadius"));
        assertEquals(Integer.toString(TEST_RESOLUTION),
                     atlasHeader.get("geometryResolution"));
        assertEquals(Integer.toString(TEST_NLGA),
                     atlasHeader.get("phaseResolutionNLga"));

        Map<String, String> cacheHeader = readHeader(result.combosCachePath());
        assertEquals("regenerated combos cache must carry the requested radius",
                     Double.toString(SWEEP_RADIUS),
                     cacheHeader.get("memberRadius"));
        assertEquals(Integer.toString(TEST_RESOLUTION),
                     cacheHeader.get("geometryResolution"));
        assertEquals("cache header comboCount must agree with the rows actually written",
                     Integer.toString(result.comboCount()),
                     cacheHeader.get("comboCount"));
    }

    /**
     * NON-VACUITY for the test above. If the harness silently ignored the
     * radius it was handed and regenerated at the baseline, every header
     * assertion would still be checking SOMETHING, but nothing about
     * {@code r}. So: the swept radius must actually change the discovered
     * physics. {@code r=0.05} is more than three times the baseline's
     * {@code 0.015} and contact is a monotone overlap predicate on
     * member radius (a pair overlapping at {@code r} still overlaps at
     * any {@code r' > r}), so the {@code r=0.05} combo universe must be a
     * STRICT SUPERSET of the baseline's at the same resolution.
     */
    @Test
    public void aLargerRadiusDiscoversStrictlyMoreCombosAtTheSameResolution() {
        List<ContactComboCache.Combo> baseline = ContactComboCache.combosFor(predicateAt(LgaTestGeometry.BASELINE_RADIUS),
                                                                              TEST_RESOLUTION,
                                                                              LgaTestGeometry.BASELINE_RADIUS);
        List<ContactComboCache.Combo> swept = ContactComboCache.combosFor(predicateAt(SWEEP_RADIUS),
                                                                            TEST_RESOLUTION,
                                                                            SWEEP_RADIUS);

        assertTrue("non-vacuity floor: the baseline sweep must discover at least one combo, was "
                   + baseline.size(), baseline.size() > 0);
        assertTrue("r=" + SWEEP_RADIUS + " must discover strictly more combos than r="
                   + LgaTestGeometry.BASELINE_RADIUS + " (monotone overlap predicate), was "
                   + swept.size() + " vs " + baseline.size(),
                   swept.size() > baseline.size());
        assertTrue("the larger radius's combo set must CONTAIN the smaller's (monotonicity, not merely a different count)",
                   swept.containsAll(baseline));
    }

    /**
     * The output-path policy is a MECHANISM, not a convention. The bead's
     * rule - "never into {@code src/test/resources/lga/} and never over
     * {@code contact-atlas-v2.tsv} or the committed combos cache" - is
     * enforced by a refusal in the harness, so a mistaken caller fails
     * loudly instead of destroying a pinned artifact.
     *
     * <h2>The three cases with {@code ..} inside a MISSING tail</h2>
     * The last three entries are not variations on {@code target/../src}:
     * they are the class of escape that survived the first fix round.
     * {@code canonicalize} resolves the deepest EXISTING ancestor for real
     * and re-appends the missing tail, so {@code ..} inside that tail is a
     * component {@code toRealPath} never sees - and {@link
     * Path#startsWith} is purely lexical. All three were ACCEPTED by the
     * un-normalized guard. {@code target/nonexistent/../../src/main/
     * resources/lga} in particular, driven through {@link
     * PerRadiusRegeneration#regenerate} with no attacker and no race,
     * created a directory INSIDE the pinned tree (see {@link
     * #refusesToWriteThroughADotDotTailThatEscapesTheTargetDirectory});
     * with a concurrent {@code mkdir target/racy} during the sweep the
     * same shape wrote both artifacts into {@code src/main/resources/
     * lga/}.
     *
     * <p>
     * FALSIFIER: removing the {@code normalize()} from {@code
     * canonicalize}'s return - the three {@code ..}-tail cases then throw
     * no {@link IllegalArgumentException} at all.
     */
    @Test
    public void refusesToWriteOutsideTheTargetDirectory() throws IOException {
        for (String forbidden : List.of("src/test/resources/lga",
                                         "src/main/resources/lga", "src",
                                         "/tmp", "target/../src",
                                         "target/../src/main/resources/lga",
                                         "target/nonexistent/../../src/main/resources/lga",
                                         "target/a/b/../../../src/main/resources/lga",
                                         "target/./nope/../../src/test/resources/lga")) {
            try {
                PerRadiusRegeneration.resolveOutputDirectory(Path.of(forbidden));
                fail("expected a refusal for output directory outside target/: "
                     + forbidden);
            } catch (IllegalArgumentException expected) {
                assertTrue("refusal must name the offending path, was: "
                           + expected.getMessage(),
                           expected.getMessage().contains(forbidden));
            }
        }
        assertEquals("an accepted directory must come back CANONICALIZED, so the write path carries no .. or symlink indirection the guard already resolved away",
                     Path.of("target").toAbsolutePath().toRealPath(),
                     PerRadiusRegeneration.resolveOutputDirectory(Path.of("target")));
        assertEquals("a not-yet-existing directory under target/ must still be accepted, canonicalized",
                     Path.of("target").toAbsolutePath().toRealPath()
                          .resolve("not-created-by-this-test"),
                     PerRadiusRegeneration.resolveOutputDirectory(Path.of("target/not-created-by-this-test")));
    }

    /**
     * THE {@code ..}-IN-A-MISSING-TAIL ESCAPE, driven through the FULL
     * {@link PerRadiusRegeneration#regenerate} entry point rather than
     * poked at {@link PerRadiusRegeneration#resolveOutputDirectory}, and
     * asserting on the pinned directory's LISTING - because the two ways
     * this failed were both invisible to a weaker test.
     *
     * <p>
     * Against the un-normalized guard this exact call ACCEPTED the path
     * and {@code ContactComboCache.write}'s {@code Files.createDirectories}
     * created {@code src/main/resources/lga/pwned} inside the PINNED tree,
     * before the file write died with {@code NoSuchFileException} - not
     * the documented {@link IllegalArgumentException}. Two things follow,
     * and both are asserted here rather than assumed:
     * <ul>
     * <li>The refusal must be the DOCUMENTED exception type. A {@code
     * NoSuchFileException} escaping instead means the guard did not fire
     * and the harness merely tripped over the consequences.</li>
     * <li>The check must be on the LISTING, not on file hashes. The
     * pollution was an empty DIRECTORY. {@link #snapshotPinnedTrees} used
     * to filter on {@code isRegularFile} and so was structurally blind to
     * those, which would have let the class-wide byte-identity bracket
     * stay green through it; that filter is gone and the bracket now
     * carries directories too, but this in-test listing assertion stays,
     * as the immediate and locally-attributable check. No pinned FILE was
     * overwritten only because
     * {@code Files.createDirectories} internally normalizes, which is
     * incidental JDK behaviour and not this harness's guard working.</li>
     * </ul>
     *
     * <p>
     * FALSIFIER: removing the {@code normalize()} from {@code
     * PerRadiusRegeneration.canonicalize}'s return.
     */
    @Test
    public void refusesToWriteThroughADotDotTailThatEscapesTheTargetDirectory() throws IOException {
        Path pinnedTree = Path.of("src/main/resources/lga");
        assertTrue("the escape's destination must be a real pinned-artifact directory for this test to mean anything",
                   Files.isDirectory(pinnedTree));
        Set<String> before = listing(pinnedTree);

        Path escape = Path.of("target/nonexistent/../../src/main/resources/lga/pwned");
        try {
            PerRadiusRegeneration.regenerate(smallSweepParameters(), escape);
            fail("expected a refusal: " + escape
                 + " has .. inside a missing tail and resolves into a PINNED artifact directory");
        } catch (IllegalArgumentException expected) {
            assertTrue("refusal must name the offending path, was: "
                       + expected.getMessage(),
                       expected.getMessage().contains(escape.toString()));
            assertTrue("refusal must name the RESOLVED destination so the escape is legible, was: "
                       + expected.getMessage(),
                       expected.getMessage()
                                .contains(pinnedTree.toAbsolutePath()
                                                     .toRealPath().toString()));
        }

        assertEquals("a refused regeneration must not have created ANYTHING in the pinned tree - not a file, and not the empty directory the un-normalized guard left behind",
                     before, listing(pinnedTree));
    }

    /**
     * THE ESCAPE THE LEXICAL GUARD MISSED, as a regression test.
     *
     * <p>
     * {@code Path.normalize()} collapses {@code ..} but is blind to
     * symlinks: {@code target/escape -> src/main/resources/lga} normalizes
     * to itself and passes a {@code startsWith("target")} comparison. This
     * is not a theoretical gap - against the lexical-only guard a real
     * {@code regenerate(params, Path.of("target/escape"))} wrote BOTH
     * artifacts into {@code src/main/resources/lga/}, i.e. the pin-capture
     * rule fully defeated by one {@code ln -s}.
     *
     * <p>
     * So this test does not merely poke {@link
     * PerRadiusRegeneration#resolveOutputDirectory}: it drives the SAME
     * full {@code regenerate} entry point the exploit used, and then
     * asserts the link target's directory LISTING is unchanged - because a
     * refusal that still created a file would be no refusal at all.
     *
     * <p>
     * FALSIFIER: reverting {@code resolveOutputDirectory} to
     * {@code normalize()}-only comparison.
     */
    @Test
    public void refusesToWriteThroughASymlinkThatEscapesTheTargetDirectory() throws IOException {
        Path linkTarget = Path.of("src/main/resources/lga");
        assertTrue("the escape's destination must be a real pinned-artifact directory for this test to mean anything",
                   Files.isDirectory(linkTarget));
        Path link = Path.of("target",
                             "escape-" + UUID.randomUUID() + "-lga");
        Files.createDirectories(Path.of("target"));
        Files.createSymbolicLink(link, linkTarget.toAbsolutePath());
        try {
            Set<String> before = listing(linkTarget);

            PerRadiusRegeneration.Parameters parameters = smallSweepParameters();
            try {
                PerRadiusRegeneration.regenerate(parameters, link);
                fail("expected a refusal: " + link
                     + " is a symlink out of target/ and into a PINNED artifact directory");
            } catch (IllegalArgumentException expected) {
                assertTrue("refusal must name the resolved destination so the escape is legible, was: "
                           + expected.getMessage(),
                           expected.getMessage()
                                    .contains(linkTarget.toAbsolutePath()
                                                         .toRealPath()
                                                         .toString()));
            }

            assertEquals("a refused regeneration must not have created ANY file in the symlink's destination",
                         before, listing(linkTarget));
        } finally {
            Files.deleteIfExists(link);
        }
    }

    private static Set<String> listing(Path directory) throws IOException {
        try (Stream<Path> entries = Files.list(directory)) {
            return entries.map(p -> p.getFileName().toString())
                           .collect(Collectors.toCollection(TreeSet::new));
        }
    }

    private static PerRadiusRegeneration.Parameters smallSweepParameters() {
        return new PerRadiusRegeneration.Parameters(TEST_RESOLUTION,
                                                     SWEEP_RADIUS, TEST_NLGA,
                                                     new Point3i(4, 4, 4),
                                                     TEST_SEED, TEST_TICKS);
    }

    /**
     * Artifact names are r-STAMPED, so two radii in the same {@code
     * target/} never collide (the pin-capture rule's "alongside old", not
     * "over old").
     */
    @Test
    public void artifactNamesAreRadiusStamped() {
        Path atlasA = PerRadiusRegeneration.atlasPathFor(Path.of("target"),
                                                          0.05);
        Path atlasB = PerRadiusRegeneration.atlasPathFor(Path.of("target"),
                                                          0.15);
        assertNotEquals("two radii must not share an atlas filename", atlasA,
                        atlasB);
        assertTrue("atlas name must declare v3 and the radius, was "
                   + atlasA.getFileName(),
                   atlasA.getFileName().toString().contains("v3")
                   && atlasA.getFileName().toString().contains("0.05"));

        Path cacheA = PerRadiusRegeneration.combosCachePathFor(Path.of("target"),
                                                                0.05);
        assertNotEquals(cacheA,
                        PerRadiusRegeneration.combosCachePathFor(Path.of("target"),
                                                                  0.15));
        assertFalse("the r-stamped cache must never collide with the committed cache's name",
                    cacheA.getFileName().toString()
                          .equals("discovered-combos-cache.tsv"));
    }

    /**
     * CODIFIES A KNOWN LIMITATION rather than leaving it as a surprise for
     * the next bead: {@code ContactAtlas.checkVersion} refuses any {@code
     * atlasVersion} other than {@link ContactAtlas#ATLAS_VERSION} (== 2),
     * so a v3 atlas is WRITE-ONLY under the current reader. Widening that
     * check belongs to whichever bead first needs to LOAD a v3 atlas; this
     * assertion exists so that bead discovers the constraint from a red
     * test naming it, not from a runtime failure mid-campaign.
     */
    @Test
    public void v3AtlasIsNotYetReadableByTheV2Reader() throws IOException {
        PerRadiusRegeneration.Result result = PerRadiusRegeneration.regenerate(smallSweepParameters());
        try {
            ContactAtlas.read(result.atlasPath());
            fail("expected the v2 reader to refuse a v3 atlas - if this now passes, ContactAtlas.checkVersion was widened and this test should become a round-trip assertion");
        } catch (ContactAtlas.HeaderMismatchException expected) {
            assertTrue("refusal must name both versions, was: "
                       + expected.getMessage(),
                       expected.getMessage().contains("atlasVersion=3"));
        }
    }

    /**
     * {@code r} is the physics-identity field stamped into both artifact
     * headers, so a non-finite one is not a harmless oddity: it produces a
     * {@code contact-atlas-v3-rNaN.tsv} whose header claims a radius no
     * comparison can ever match. {@code NaN} defeats every ordering
     * comparison, so the obvious {@code radius <= 0} check waves it
     * through - the reason the guard is spelled with {@code isFinite}.
     */
    @Test
    public void rejectsANonFiniteRadiusRatherThanStampingItIntoAHeader() {
        for (double bad : new double[] { Double.NaN,
                                          Double.POSITIVE_INFINITY,
                                          Double.NEGATIVE_INFINITY, 0.0,
                                          -0.015 }) {
            try {
                new PerRadiusRegeneration.Parameters(TEST_RESOLUTION, bad,
                                                      TEST_NLGA,
                                                      new Point3i(4, 4, 4),
                                                      TEST_SEED, TEST_TICKS);
                fail("expected a refusal for memberRadius=" + bad);
            } catch (IllegalArgumentException expected) {
                assertTrue("refusal must name the offending radius, was: "
                           + expected.getMessage(),
                           expected.getMessage()
                                    .contains(Double.toString(bad)));
            }
        }
    }

    private static ContactPredicate predicateAt(double radius) {
        return new ContactPredicate(new MemberGeometry(TEST_RESOLUTION,
                                                        radius));
    }

    /**
     * @return every ENTRY under both pinned resource trees - regular files
     *         mapped to their SHA-256, directories (and anything else) to
     *         {@link #NOT_A_REGULAR_FILE}. Comparing two of these detects
     *         mutation of a known artifact, DELETION of one, and - the
     *         case a fixed filename list structurally cannot see -
     *         ADDITION of a new entry into a pinned tree. Directories are
     *         included deliberately and not incidentally: the pollution
     *         actually observed against the un-normalized path guard was
     *         an EMPTY directory ({@code src/main/resources/lga/pwned}),
     *         which a regular-files-only walk is structurally blind to,
     *         and which would therefore have passed this bracket
     *         unnoticed.
     */
    private static Map<String, String> snapshotPinnedTrees() throws IOException {
        Map<String, String> snapshot = new TreeMap<>();
        for (String tree : PINNED_TREES) {
            try (Stream<Path> walk = Files.walk(Path.of(tree))) {
                for (Path entry : walk.toList()) {
                    snapshot.put(entry.toString(),
                                 Files.isRegularFile(entry) ? sha256(Files.readAllBytes(entry))
                                                             : NOT_A_REGULAR_FILE);
                }
            }
        }
        for (String pinned : PINNED_ARTIFACTS) {
            assertTrue("pinned artifact must exist and be covered by the tree snapshot (a snapshot that stopped seeing these would compare empty to empty): "
                       + pinned, snapshot.containsKey(pinned));
        }
        return snapshot;
    }

    private static String sha256(byte[] bytes) {
        try {
            byte[] digest = MessageDigest.getInstance("SHA-256")
                                          .digest(bytes);
            StringBuilder sb = new StringBuilder(digest.length * 2);
            for (byte b : digest) {
                sb.append(Character.forDigit((b >> 4) & 0xF, 16))
                  .append(Character.forDigit(b & 0xF, 16));
            }
            return sb.toString();
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 unavailable", e);
        }
    }

    private static Map<String, String> readHeader(Path path) throws IOException {
        Map<String, String> kv = new LinkedHashMap<>();
        for (String line : Files.readAllLines(path, StandardCharsets.UTF_8)) {
            if (!line.startsWith("#")) {
                break;
            }
            int eq = line.indexOf('=');
            if (eq > 0) {
                kv.put(line.substring(1, eq).trim(),
                       line.substring(eq + 1).trim());
            }
        }
        return kv;
    }
}
