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
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Conformance tests for {@link FineStepContactTable} (USER DECISION
 * 2026-08-08, FINAL, bead inviscid-0nx.21): the fine-step contact-firing
 * anchor that replaces {@link ContactTable}'s bin-level ANY-OVERLAP
 * transcription as the LGA's contact-firing decision.
 *
 * @author halhildebrand
 */
public class FineStepContactTableTest {

    private static final int    GEOMETRY_RESOLUTION = LgaTestGeometry.BASELINE_GEOMETRY_RESOLUTION;
    private static final double RADIUS              = LgaTestGeometry.BASELINE_RADIUS;
    private static final String RESOURCE_PATH       = "lga/contact-atlas-v2.tsv";

    /**
     * {@code FineStepContactTable}'s identity-index array size: {@code
     * CANONICAL_DIRECTIONS * CUBES_PER_CELL * MEMBERS_PER_CUBE *
     * CUBES_PER_CELL * MEMBERS_PER_CUBE == 6 * 5*6 * 5*6}. Mirrored here
     * (the production fields are private) so this suite can derive the
     * exact expected footprint rather than band-match it.
     */
    private static final int IDENTITY_DOMAIN = 6 * 5 * 6 * 5 * 6;

    /**
     * The geometry resolution every {@link ComboPin#footprintBytes} in
     * {@link #PINNED_BY_RADIUS} was MEASURED at, as a hard literal -
     * deliberately NOT {@link #GEOMETRY_RESOLUTION}, which is an alias of
     * a constant that can move. Footprint scales with {@code
     * ceil(res^2 / 64)}, so a resolution move invalidates every registered
     * literal at once; {@link
     * #pinnedFootprintsAgreeWithTheStructuralIdentity} compares the two
     * and says so in those terms rather than reporting each pin as a typo.
     */
    private static final int PIN_MEASUREMENT_RESOLUTION = 360;

    /**
     * A per-radius combo/footprint expectation, MEASURED at that radius.
     *
     * @param discovered      combos {@link ContactComboCache} discovers at
     *                        this radius (both direction signs)
     * @param canonical       the positive-direction subset {@link
     *                        FineStepContactTable} actually stores
     * @param footprintBytes  {@code memoryFootprintBytes()} at this radius
     *                        and {@link #GEOMETRY_RESOLUTION}
     */
    private record ComboPin(int discovered, int canonical,
                             long footprintBytes) {
    }

    /**
     * PER-RADIUS PINS, keyed by a HARD LITERAL radius (bead
     * inviscid-0nx.30, E.3).
     *
     * <h2>Why a fail-closed registry and not {@code if (RADIUS ==
     * BASELINE_RADIUS)}</h2>
     * The obvious spelling - guard the exact assertions behind a
     * comparison against {@link LgaTestGeometry#BASELINE_RADIUS} - is a
     * SILENT-SKIP TRAP. {@code design-seeding-radius.md} §D-B exists
     * precisely to MOVE that anchor off {@code 0.015}; the first time
     * someone does, the guard goes false, the {@code 446} pin evaporates,
     * and the suite stays green while asserting nothing about combo
     * counts. A test that stops testing when you change the thing it was
     * pinning is worse than no test.
     *
     * <p>
     * So: look the radius up here. An UNKNOWN radius FAILS, naming what
     * must be measured. Changing the anchor therefore turns the suite RED
     * and forces the new radius's counts to be measured and recorded -
     * the bead's "must become per-r or anchored-r assertions" as a
     * mechanism rather than a convention.
     *
     * <p>
     * The key is deliberately the literal {@code 0.015d} and NOT {@code
     * LgaTestGeometry.BASELINE_RADIUS}. If key and lookup moved together,
     * retargeting the anchor to (say) {@code 0.05} would match a pin
     * holding counts measured at {@code 0.015} - a FALSE GREEN, strictly
     * worse than the silent skip this registry replaces.
     *
     * <h2>Measured evidence that these counts really are r-sensitive</h2>
     * The {@code 0.05} entry is not decoration. A full-fidelity {@link
     * PerRadiusRegeneration} run at {@code r=0.05}, {@code
     * geometryResolution=360} (bead inviscid-0nx.30, 2026-08-09) discovered
     * <strong>832</strong> combos, <strong>416</strong> of them canonical -
     * against {@code 446}/{@code 223} at the baseline. The single candidate
     * radius the design memo names first moves BOTH by ~87%, so the
     * baseline numbers below are r-sensitive facts, not exact-looking
     * constants.
     *
     * <h2>On the {@code footprintBytes} of the second entry</h2>
     * It is DERIVED from {@link #footprintOf}, not separately observed, and
     * an earlier draft withheld the whole entry for that reason. That
     * reasoning was self-refuting: {@code 0.015}'s
     * {@code 3_634_200} is equally derivable, and {@link
     * #memoryFootprintIsExactlyTheCanonicalOnlyPackedSize} asserts the
     * derivation as an EXACT identity at every radius. Withholding a second entry on a
     * ground that condemns the first bought nothing and left the registry
     * unable to demonstrate the multi-radius behaviour it exists for.
     *
     * <h2>Exactly how much each field is worth, per radius</h2>
     * Stated field by field so nobody over-reads the registry.
     * <ul>
     * <li>At the radius CURRENTLY UNDER TEST, {@code discovered} and
     * {@code canonical} are checked against a live measurement by {@link
     * #comboCountsArePinnedAtThisRadius}. That is where this registry's
     * falsifying power actually lives.</li>
     * <li>At any OTHER registered radius, {@code discovered} and {@code
     * canonical} are unfalsifiable in-suite: the only assertion that
     * reaches them is {@link
     * #pinnedFootprintsAgreeWithTheStructuralIdentity}'s {@code canonical
     * < discovered}, which any plausible transcription satisfies. The
     * {@code 832}/{@code 416} pair rests on the recorded {@link
     * PerRadiusRegeneration} run and nothing else; a typo in it would sit
     * inert until the anchor moved onto {@code 0.05}, and would then
     * surface as a combo-count mismatch. That is a real gap, and closing
     * it means measuring at that radius in-suite - minutes of sweep - not
     * writing another assertion over the same literals.</li>
     * <li>{@code footprintBytes} adds NO falsifying power at ANY radius:
     * it is derived from {@code canonical} by {@link #footprintOf}, the
     * same function {@link
     * #memoryFootprintIsExactlyTheCanonicalOnlyPackedSize} checks the live
     * table against. Its job is to keep a registered entry internally
     * consistent, which {@link
     * #pinnedFootprintsAgreeWithTheStructuralIdentity} enforces at every
     * registered radius.</li>
     * </ul>
     */
    private static final Map<Double, ComboPin> PINNED_BY_RADIUS = Map.of(0.015d,
                                                                          new ComboPin(446,
                                                                                        223,
                                                                                        3_634_200L),
                                                                          0.05d,
                                                                          new ComboPin(832,
                                                                                        416,
                                                                                        6_760_800L));

    /**
     * {@code FineStepContactTable}'s exact footprint for {@code canonical}
     * stored combos at {@link #GEOMETRY_RESOLUTION}: a packed bitset of
     * {@code ceil(res^2 / 64)} longs per combo, plus the identity index.
     * The single spelling of the identity that {@link
     * #memoryFootprintIsExactlyTheCanonicalOnlyPackedSize} checks the live
     * table against and {@link
     * #pinnedFootprintsAgreeWithTheStructuralIdentity} checks the
     * registry's literals against.
     */
    private static long footprintOf(int canonical) {
        long wordsPerCombo = ((long) GEOMETRY_RESOLUTION * GEOMETRY_RESOLUTION
                              + 63) / 64;
        return (long) canonical * wordsPerCombo * 8L + IDENTITY_DOMAIN * 4L;
    }

    private static ContactPredicate    PREDICATE;
    private static FineStepContactTable FINE;
    private static ContactAtlas        ATLAS;
    private static ContactTable        BIN_TABLE;
    private static List<ContactComboCache.Combo> DISCOVERED;
    private static List<ContactComboCache.Combo> CANONICAL;

    @BeforeClass
    public static void buildFixtures() throws IOException {
        // FAIL FAST on an unregistered radius. THIS is the falsifier for
        // "every radius under test has a recorded pin" - not the
        // identically-named test below, which by construction can never
        // run on the failing path because this aborts the whole class
        // first. The condition and the remediation message were MOVED
        // here, not duplicated for belt-and-braces; the named test is what
        // is left behind, and it is documentation with an @Test annotation.
        //
        // The move is about WHEN the failure is reported. Reaching the
        // named test costs a full resolution-360 build plus an exhaustive
        // sweep - minutes of work whose only possible outcome is the
        // message already known right here.
        if (!PINNED_BY_RADIUS.containsKey(RADIUS)) {
            throw new AssertionError("no ComboPin recorded for memberRadius="
                                      + RADIUS
                                      + " - run PerRadiusRegeneration at this radius, read comboCount from the generated cache header, and add a ComboPin entry to PINNED_BY_RADIUS (see everyRadiusUnderTestHasARecordedComboPin)");
        }
        PREDICATE = new ContactPredicate(new MemberGeometry(GEOMETRY_RESOLUTION,
                                                              RADIUS));
        long start = System.nanoTime();
        FINE = FineStepContactTable.buildFor(PREDICATE, GEOMETRY_RESOLUTION,
                                              RADIUS);
        long elapsedMs = (System.nanoTime() - start) / 1_000_000;
        System.out.println("FineStepContactTable build: " + elapsedMs + "ms");

        DISCOVERED = ContactComboCache.combosFor(PREDICATE,
                                                  GEOMETRY_RESOLUTION, RADIUS);
        CANONICAL = DISCOVERED.stream()
                                .filter(c -> c.direction() >= 1
                                             && c.direction() <= 6)
                                .toList();

        URL resource = FineStepContactTableTest.class.getClassLoader()
                                                        .getResource(RESOURCE_PATH);
        ATLAS = ContactAtlas.read(Paths.get(resource.getPath()));
        BIN_TABLE = ContactTable.of(ATLAS);
    }

    /**
     * FIXTURE-CONSISTENCY GUARD (new with the E.3 literal consolidation).
     * {@link #FINE} is built at {@link #RADIUS} while {@link #ATLAS} and
     * {@link #BIN_TABLE} are read from the COMMITTED {@code
     * contact-atlas-v2.tsv}, whose radius is fixed at {@code 0.015}. Those
     * two were independent literals before this bead; they are one edit
     * apart now. If the shared anchor moves without the atlas being
     * regenerated, {@link #fineFiringImpliesBinLevelAnyOverlapFiring}
     * would compare a fine table at the new radius against a bin table at
     * the old one and fail for an entirely misleading reason. This test
     * exists to say what to do about it.
     *
     * <p>
     * It does NOT run first, and no claim that it does would be true:
     * JUnit 4's {@code MethodSorters.DEFAULT} orders by a hash of the
     * method name, and at {@code r=0.05} the misleading {@link
     * #fineFiringImpliesBinLevelAnyOverlapFiring} failure was observed
     * alongside this one. Read this one's message first when several go
     * red together.
     *
     * <p>
     * FALSIFIER: any radius change not accompanied by a regenerated atlas.
     */
    @Test
    public void committedAtlasMatchesTheRadiusTheFineTableIsBuiltAt() {
        assertEquals("the committed atlas was generated at a DIFFERENT radius than this suite builds the fine table at - regenerate via PerRadiusRegeneration and repoint RESOURCE_PATH, do not just move the constant",
                     RADIUS, ATLAS.header().memberRadius(), 0.0);
    }

    /**
     * The registry's own non-vacuity rule, WRITTEN DOWN WHERE IT IS
     * FINDABLE: a radius with no recorded pin must fail loudly rather than
     * let every count assertion below quietly degrade to structure-only
     * checking.
     *
     * <p>
     * BE HONEST ABOUT WHAT THIS METHOD IS. It cannot fail. {@link
     * #buildFixtures} asserts the identical condition, with the identical
     * remediation message, before this class's fixtures are built - so on
     * the only input that could red this one, the class is already
     * aborted. The falsifier is not lost, it MOVED; what stands here is a
     * named, greppable statement of the rule for someone reading the
     * suite, addressed to a human rather than to the build. It is not a
     * second line of defence, and counting it as one would be counting the
     * same check twice.
     *
     * <p>
     * FALSIFIER: none reachable, by construction - see {@link
     * #buildFixtures}, which carries it.
     */
    @Test
    public void everyRadiusUnderTestHasARecordedComboPin() {
        assertTrue("no ComboPin recorded for memberRadius=" + RADIUS
                   + " - run PerRadiusRegeneration at this radius, read comboCount from the generated cache header, and add a ComboPin entry. Do NOT delete this assertion: the pins are what make the count assertions r-sensitive.",
                   PINNED_BY_RADIUS.containsKey(RADIUS));
    }

    /**
     * Build-vs-artifact decision pin: measured comfortably under the ~10s
     * decision threshold (build() itself already timed and printed in
     * {@link #buildFixtures()}; this re-times a cache hit, which must be
     * effectively instantaneous).
     */
    @Test
    public void cachedRebuildIsFast() {
        long start = System.nanoTime();
        FineStepContactTable again = FineStepContactTable.buildFor(PREDICATE,
                                                                     GEOMETRY_RESOLUTION,
                                                                     RADIUS);
        long elapsedMs = (System.nanoTime() - start) / 1_000_000;
        assertSame("cache must return the identical instance", FINE, again);
        assertTrue("cache hit must be well under 10ms, was " + elapsedMs,
                   elapsedMs < 10);
    }

    /**
     * LAYER 3 - NON-VACUITY FLOOR. Everything below derives from the
     * loaded combo cache, so a cache that were empty, truncated, or
     * missing its negative half would make the derived assertions
     * trivially satisfiable. This fixes the floor first.
     *
     * <p>
     * FALSIFIER (r): a radius small enough that NO pair ever overlaps
     * (as {@code r -> 0} the contact set empties and {@code canonical}
     * goes to 0). Also fires on a corrupt/half-written cache at any
     * radius - {@code discovered > canonical} is what proves BOTH
     * direction signs survived the load, which no other assertion here
     * checks.
     */
    @Test
    public void loadedComboCacheIsNonDegenerate() {
        assertTrue("non-vacuity floor: the cache must discover at least one combo, was "
                   + DISCOVERED.size(), DISCOVERED.size() > 0);
        assertTrue("non-vacuity floor: at least one CANONICAL (positive-direction) combo must exist or the fine table stores nothing, was "
                   + CANONICAL.size(), CANONICAL.size() > 0);
        assertTrue("the discovered set must contain negative-direction combos too (discovered="
                   + DISCOVERED.size() + ", canonical=" + CANONICAL.size()
                   + ") - if these are equal the cache lost its negative half",
                   DISCOVERED.size() > CANONICAL.size());
        assertEquals("the fine table must store exactly the canonical subset",
                     CANONICAL.size(), FINE.comboCount());
    }

    /**
     * LAYER 1 - STRUCTURAL IDENTITY. The footprint is not approximately
     * anything: it is EXACTLY {@code canonicalCombos * ceil(res^2 / 64) *
     * 8 bytes} of packed bitset plus {@code IDENTITY_DOMAIN * 4 bytes} of
     * index. Deriving it from the loaded cache's own canonical count
     * replaces the old {@code 2.5 MiB < x < 5.0 MiB} band, which passed
     * for any packing scheme landing anywhere in a 2x window.
     *
     * <p>
     * FALSIFIER (r): NONE, and that is stated deliberately - this is an
     * algebraic identity of the packing code and holds at EVERY radius.
     * Its r-sensitivity comes from {@link
     * #comboCountsArePinnedAtThisRadius}, not from here. What it DOES
     * falsify: storing all discovered combos instead of the canonical
     * subset (footprint doubles - the "Canonical directions only" claim in
     * {@code FineStepContactTable}'s Javadoc, otherwise untested);
     * switching from bitset to byte-per-cell (8x); changing {@code
     * CANONICAL_DIRECTIONS} away from 6.
     */
    @Test
    public void memoryFootprintIsExactlyTheCanonicalOnlyPackedSize() {
        long wordsPerCombo = ((long) GEOMETRY_RESOLUTION * GEOMETRY_RESOLUTION
                              + 63) / 64;

        assertEquals("footprint must equal canonicalCombos(" + CANONICAL.size()
                     + ") * wordsPerCombo(" + wordsPerCombo
                     + ") * 8 + identityIndex(" + IDENTITY_DOMAIN
                     + ") * 4 - a mismatch means the storage scheme changed, NOT that the radius changed",
                     footprintOf(CANONICAL.size()),
                     FINE.memoryFootprintBytes());

        assertTrue("floor: the footprint must exceed the bare identity index, or no combo bitset was stored at all",
                   FINE.memoryFootprintBytes() > (long) IDENTITY_DOMAIN * 4L);
    }

    /**
     * REGISTRY SELF-CONSISTENCY, at every registered radius - including
     * radii not currently under test, which no other assertion here
     * reaches.
     *
     * <p>
     * {@code footprintBytes} is derived from {@code canonical} by {@link
     * #footprintOf}, so a hand-entered pair can disagree, and a wrong
     * literal would sit inert until the anchor moved onto it and then fail
     * as a confusing footprint mismatch rather than a typo. Checking it at
     * registration time instead turns that into an immediate, local red.
     *
     * <p>
     * FALSIFIER: any registered {@code (canonical, footprintBytes)} pair
     * that does not satisfy the packing identity - e.g. registering
     * {@code r=0.05}'s 416 canonical combos against the baseline's
     * {@code 3_634_200}.
     *
     * <p>
     * THE RESOLUTION GUARD COMES FIRST, and is not decoration. Every
     * registered {@code footprintBytes} was measured at {@link
     * #PIN_MEASUREMENT_RESOLUTION}, while {@link #footprintOf} computes
     * against the CURRENT {@link #GEOMETRY_RESOLUTION}. If the suite's
     * resolution ever moves, every pin in the registry fails at once -
     * and, without the guard, each with a message accusing its literal of
     * being a typo, which is the one thing that would not be true. The
     * guard converts "all your literals are wrong" into "the resolution
     * moved; the literals are stale by construction and must be
     * re-measured", which is what actually happened.
     */
    @Test
    public void pinnedFootprintsAgreeWithTheStructuralIdentity() {
        assertEquals("every footprintBytes in PINNED_BY_RADIUS was MEASURED at geometryResolution="
                     + PIN_MEASUREMENT_RESOLUTION
                     + ", but this suite now builds at " + GEOMETRY_RESOLUTION
                     + ". The registered literals are not typos - they are stale by construction, because footprint scales with ceil(res^2/64). Re-measure each radius via PerRadiusRegeneration at the new resolution and re-register, then update PIN_MEASUREMENT_RESOLUTION.",
                     PIN_MEASUREMENT_RESOLUTION, GEOMETRY_RESOLUTION);
        assertTrue("the registry must hold more than one radius, or it demonstrates nothing about per-r behaviour",
                   PINNED_BY_RADIUS.size() > 1);
        PINNED_BY_RADIUS.forEach((radius, pin) -> {
            assertTrue("a pin's canonical subset must be a strict subset of what was discovered at r="
                       + radius, pin.canonical() < pin.discovered());
            assertEquals("footprintBytes for the pin at r=" + radius
                         + " must equal the packing identity applied to its own canonical count ("
                         + pin.canonical() + ")",
                         footprintOf(pin.canonical()), pin.footprintBytes());
        });
    }

    /**
     * LAYER 2 - THE r-ANCHORED PIN. The exact counts, per radius, from
     * {@link #PINNED_BY_RADIUS}.
     *
     * <p>
     * FALSIFIER (r): ANY radius whose exhaustive sweep discovers a
     * different number of ever-contacting combos - which is to say, any
     * radius that changes the physics this campaign is measuring. Contact
     * is a MONOTONE overlap predicate on member radius (a pair overlapping
     * at {@code r} still overlaps at every {@code r' > r}), so every
     * {@code r > 0.015} that crosses even ONE combo's discovery threshold
     * fails this, rising toward the full 10,800-combo universe as {@code
     * r} approaches the universal-contact wall at {@code ~6.3} (design
     * memo §D-B); every {@code r < 0.015} that drops a combo fails it too.
     * The candidate sweep values in the design memo - {@code 0.05}, {@code
     * 0.15}, {@code 0.5} - all fail this assertion until measured and
     * registered, which is the intent.
     *
     * <p>
     * This is also what makes {@link
     * #memoryFootprintIsExactlyTheCanonicalOnlyPackedSize} and {@link
     * #fineTableAgreesWithLiveContactPredicateAtEveryStepPairExhaustive}
     * r-sensitive: both derive from {@code CANONICAL.size()}, which this
     * test pins exactly.
     */
    @Test
    public void comboCountsArePinnedAtThisRadius() {
        ComboPin pin = PINNED_BY_RADIUS.get(RADIUS);
        assertNotNull("no ComboPin for memberRadius=" + RADIUS
                      + " - see everyRadiusUnderTestHasARecordedComboPin",
                      pin);

        assertEquals("discovered combo count at memberRadius=" + RADIUS,
                     pin.discovered(), DISCOVERED.size());
        assertEquals("canonical (positive-direction) combo count at memberRadius="
                     + RADIUS, pin.canonical(), CANONICAL.size());
        assertEquals("memoryFootprintBytes at memberRadius=" + RADIUS
                     + " and geometryResolution=" + GEOMETRY_RESOLUTION,
                     pin.footprintBytes(), FINE.memoryFootprintBytes());
    }

    @Test
    public void geometryResolutionMatchesTheAtlas() {
        assertEquals(GEOMETRY_RESOLUTION, FINE.geometryResolution());
        assertEquals(ATLAS.header().geometryResolution(), FINE.geometryResolution());
    }

    /**
     * Test 3a (bead's naming): the fine table agrees with a FRESH,
     * independent {@link ContactPredicate} evaluation at EVERY (combo,
     * stepA, stepB) - exhaustive, not sampled (measured feasible: ~4-5s,
     * comfortably fast enough to run exhaustively rather than sample).
     * This is the new transcription-fidelity anchor.
     */
    @Test
    public void fineTableAgreesWithLiveContactPredicateAtEveryStepPairExhaustive() {
        List<ContactComboCache.Combo> combos = DISCOVERED;
        long checked = 0;
        long mismatches = 0;
        for (ContactComboCache.Combo combo : combos) {
            if (combo.direction() < 1 || combo.direction() > 6) {
                continue;
            }
            for (int a = 0; a < GEOMETRY_RESOLUTION; a++) {
                float angleA = ContactComboCache.angleOf(a, GEOMETRY_RESOLUTION);
                for (int b = 0; b < GEOMETRY_RESOLUTION; b++) {
                    float angleB = ContactComboCache.angleOf(b,
                                                              GEOMETRY_RESOLUTION);
                    boolean live = PREDICATE.contacts(combo.cubeA(),
                                                       combo.memberA(), angleA,
                                                       combo.cubeB(),
                                                       combo.memberB(), angleB,
                                                       combo.direction());
                    boolean fine = FINE.contacts(combo.direction(),
                                                  combo.cubeA(),
                                                  combo.memberA(), a,
                                                  combo.cubeB(),
                                                  combo.memberB(), b);
                    if (live != fine) {
                        mismatches++;
                    }
                    checked++;
                }
            }
        }
        // LAYER 1 - STRUCTURAL IDENTITY, replacing the old
        // `checked > 20_000_000L` band. That band tolerated the canonical
        // count sitting anywhere in [155, infinity): a loop nest that
        // silently skipped a third of the combos still passed it. The
        // exact product does not.
        //
        // FALSIFIER (r): none directly - CANONICAL.size() is the r-carrier
        // and comboCountsArePinnedAtThisRadius pins it to 223 at r=0.015,
        // which pins `checked` to 28,900,800 transitively. What THIS line
        // falsifies is a loop-nest bug: any combo skipped, any step pair
        // not visited.
        assertEquals("exhaustive sweep must visit exactly canonicalCombos("
                     + CANONICAL.size() + ") x " + GEOMETRY_RESOLUTION + "^2 pairs",
                     (long) CANONICAL.size() * GEOMETRY_RESOLUTION
                             * GEOMETRY_RESOLUTION,
                     checked);
        assertTrue("non-vacuity floor: the exhaustive sweep must actually check something, was "
                   + checked, checked > 0);
        assertEquals("fine table disagreed with live ContactPredicate on "
                     + mismatches + " of " + checked + " (combo, stepA, stepB) triples",
                     0L, mismatches);
    }

    /**
     * Non-vacuity for a combo NOT in the discovered set: must report no
     * contact everywhere (the "-1 sentinel" path), never throw.
     */
    @Test
    public void unknownComboNeverContacts() {
        // direction=1, cubeA=0, memberA=0, cubeB=0, memberB=0 is
        // extremely unlikely to be a real contacting combo (adjacent-cube
        // same-member pairs at direction 1 are geometrically implausible
        // at this radius) - verified empirically below rather than
        // asserted a priori.
        boolean anyContact = false;
        for (int a = 0; a < GEOMETRY_RESOLUTION && !anyContact; a++) {
            for (int b = 0; b < GEOMETRY_RESOLUTION; b++) {
                if (FINE.contacts(1, 0, 0, a, 0, 0, b)) {
                    anyContact = true;
                    break;
                }
            }
        }
        assertTrue("fixture assumption: (dir=1,0,0,0,0) must be a non-contacting combo for this test to probe the sentinel path",
                   !anyContact);
    }

    /**
     * Test 3b (bead's naming): CONSISTENCY - fine-step firing implies
     * bin-level ANY-OVERLAP firing (fine true => .19 table true for the
     * containing bins). Exhaustive over the fine table's own fired bits
     * (cheap: only ~28,384 of 57.8M pairs fire).
     */
    @Test
    public void fineFiringImpliesBinLevelAnyOverlapFiring() {
        int nLga = BIN_TABLE.nLga();
        List<ContactComboCache.Combo> combos = DISCOVERED;
        long fineFires = 0;
        long subsetViolations = 0;
        for (ContactComboCache.Combo combo : combos) {
            if (combo.direction() < 1 || combo.direction() > 6) {
                continue;
            }
            for (int a = 0; a < GEOMETRY_RESOLUTION; a++) {
                for (int b = 0; b < GEOMETRY_RESOLUTION; b++) {
                    if (!FINE.contacts(combo.direction(), combo.cubeA(),
                                        combo.memberA(), a, combo.cubeB(),
                                        combo.memberB(), b)) {
                        continue;
                    }
                    fineFires++;
                    int binA = ContactAtlasGenerator.binOfStep(a, nLga,
                                                                GEOMETRY_RESOLUTION);
                    int binB = ContactAtlasGenerator.binOfStep(b, nLga,
                                                                GEOMETRY_RESOLUTION);
                    if (!BIN_TABLE.contacts(combo.direction(), combo.cubeA(),
                                             combo.memberA(), binA,
                                             combo.cubeB(), combo.memberB(),
                                             binB)) {
                        subsetViolations++;
                    }
                }
            }
        }
        assertTrue("non-vacuity: expected the fine table to fire somewhere",
                   fineFires > 0);
        assertEquals(subsetViolations + " of " + fineFires
                     + " fine-fired (combo,stepA,stepB) triples did NOT have "
                     + "their containing bin fire in the ANY-OVERLAP bin table "
                     + "(subset relation violated)",
                     0L, subsetViolations);
    }
}
