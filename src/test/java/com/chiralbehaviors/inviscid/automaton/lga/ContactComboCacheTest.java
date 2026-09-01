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

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Assume;
import org.junit.Test;

import com.chiralbehaviors.inviscid.automaton.lga.ContactComboCache.Combo;

/**
 * Behavioral tests for {@link ContactComboCache} (bead inviscid-gyt).
 * {@link ContactComboCache#sweepExhaustively} is exercised at a small,
 * cheap {@code geometryResolution} (8, the minimum {@link MemberGeometry}
 * accepts - divisible by 8) rather than the production {@code 360}: a
 * live, uncached exhaustive sweep at 360 costs several minutes (see that
 * class's Javadoc), which is exactly the cost this class's caching exists
 * to avoid paying inside the ordinary test suite. The PRODUCTION
 * {@code (360, RADIUS)} pair is exercised via the checked-in cache
 * resource, not a live sweep, in {@link #cachedCombosMatchTheResourceHeaderCount()}.
 *
 * @author halhildebrand
 */
public class ContactComboCacheTest {

    private static final int    SMALL_RESOLUTION = 8;
    private static final double RADIUS           = ContactAtlasGenerator.RADIUS;

    private static ContactPredicate newPredicate(int resolution) {
        return new ContactPredicate(new MemberGeometry(resolution, RADIUS));
    }

    /**
     * bead inviscid-gyt Phase A gate rework: {@link
     * ContactComboCache#angleOf}'s step-CENTER, {@code Constants.TWO_PI}
     * (float)-based reconstruction must round-trip through {@link
     * ContactAtlasGenerator#stepOf} for EVERY step - i.e. the angle this
     * class hands to {@link ContactPredicate#contacts} for "step s" is
     * genuinely quantized back to step {@code s} by {@link MemberGeometry}
     * internally, not some neighboring step. An earlier edge-sampled,
     * double-precision-{@code TWO_PI} version of this method failed this
     * property for ~38% of steps at {@code resolution=360} - silently
     * evaluating the WRONG step's geometry for over a third of the fine
     * sweep, which is what originally produced observedCount>0/
     * overlapFraction==0 anomalies in the generated atlas.
     */
    @Test
    public void angleOfRoundTripsThroughStepOfForEveryStep() {
        int resolution = ContactAtlasGenerator.GEOMETRY_RESOLUTION;
        int mismatches = 0;
        for (int step = 0; step < resolution; step++) {
            float angle = ContactComboCache.angleOf(step, resolution);
            int roundTripped = ContactAtlasGenerator.stepOf(angle, resolution);
            if (roundTripped != step) {
                mismatches++;
            }
        }
        assertEquals("expected every step to round-trip through stepOf exactly, found "
                     + mismatches + "/" + resolution + " mismatches", 0,
                     mismatches);
    }

    /**
     * Non-vacuity precondition: a live sweep at a small resolution finds a
     * proper, non-empty, non-total subset of the {@code 12*5*6*5*6==10,800}
     * combo universe - neither degenerate extreme.
     */
    @Test
    public void exhaustiveSweepFindsAProperNonEmptySubset() {
        ContactPredicate predicate = newPredicate(SMALL_RESOLUTION);
        List<Combo> combos = ContactComboCache.sweepExhaustively(predicate,
                                                                   SMALL_RESOLUTION);

        long totalUniverse = (long) FccNeighborhood.DIRECTIONS.size() * 5 * 6
                              * 5 * 6;
        assertFalse("expected at least one ever-contacting combo",
                    combos.isEmpty());
        assertTrue("expected a proper subset of the full combo universe ("
                   + totalUniverse + "), found " + combos.size(),
                   combos.size() < totalUniverse);
    }

    /**
     * Determinism: two independent sweeps at the same resolution over the
     * same (stateless) predicate produce the identical combo set, in the
     * identical order (fixed nested-loop iteration, no hash-ordered
     * collection anywhere in the path).
     */
    @Test
    public void exhaustiveSweepIsDeterministic() {
        ContactPredicate predicateA = newPredicate(SMALL_RESOLUTION);
        ContactPredicate predicateB = newPredicate(SMALL_RESOLUTION);

        List<Combo> first = ContactComboCache.sweepExhaustively(predicateA,
                                                                  SMALL_RESOLUTION);
        List<Combo> second = ContactComboCache.sweepExhaustively(predicateB,
                                                                   SMALL_RESOLUTION);

        assertEquals(first, second);
    }

    /**
     * Every combo the sweep returns genuinely contacts somewhere on the
     * same grid it was discovered from (spot-checked directly against a
     * live {@link ContactPredicate}, not just trusted).
     */
    @Test
    public void everyDiscoveredComboGenuinelyContactsSomewhere() {
        ContactPredicate predicate = newPredicate(SMALL_RESOLUTION);
        List<Combo> combos = ContactComboCache.sweepExhaustively(predicate,
                                                                   SMALL_RESOLUTION);
        assertFalse(combos.isEmpty());

        for (Combo combo : combos) {
            boolean contactsSomewhere = false;
            outer:
            for (int a = 0; a < SMALL_RESOLUTION; a++) {
                float angleA = ContactComboCache.angleOf(a, SMALL_RESOLUTION);
                for (int b = 0; b < SMALL_RESOLUTION; b++) {
                    float angleB = ContactComboCache.angleOf(b,
                                                              SMALL_RESOLUTION);
                    if (predicate.contacts(combo.cubeA(), combo.memberA(),
                                           angleA, combo.cubeB(),
                                           combo.memberB(), angleB,
                                           combo.direction())) {
                        contactsSomewhere = true;
                        break outer;
                    }
                }
            }
            assertTrue("combo " + combo + " was discovered but does not contact anywhere",
                       contactsSomewhere);
        }
    }

    /**
     * {@link ContactComboCache#combosFor} falls back to a live sweep (and
     * does not crash / silently return the wrong geometry's combos) for a
     * {@code (geometryResolution, memberRadius)} pair the checked-in cache
     * resource does not cover - the staleness/miss-handling half of the
     * cache contract. Deliberately cheap ({@link #SMALL_RESOLUTION}, not
     * the production 360) so this stays a fast test.
     */
    @Test
    public void combosForFallsBackToLiveSweepForAnUncachedResolution() {
        ContactPredicate predicate = newPredicate(SMALL_RESOLUTION);
        List<Combo> viaCombosFor = ContactComboCache.combosFor(predicate,
                                                                 SMALL_RESOLUTION,
                                                                 RADIUS);
        List<Combo> viaDirectSweep = ContactComboCache.sweepExhaustively(newPredicate(SMALL_RESOLUTION),
                                                                           SMALL_RESOLUTION);
        assertEquals(viaDirectSweep, viaCombosFor);
    }

    /**
     * {@link ContactComboCache#combosFor} caches per {@code
     * (geometryResolution, memberRadius)} pair - a second call with the
     * same pair returns the identical (== reference-equal via {@code
     * List.copyOf} + the cache map, not merely equal-by-value) list rather
     * than resweeping.
     */
    @Test
    public void combosForCachesRepeatedCallsForTheSamePair() {
        ContactPredicate predicateA = newPredicate(SMALL_RESOLUTION);
        ContactPredicate predicateB = newPredicate(SMALL_RESOLUTION);

        List<Combo> first = ContactComboCache.combosFor(predicateA,
                                                          SMALL_RESOLUTION,
                                                          RADIUS);
        List<Combo> second = ContactComboCache.combosFor(predicateB,
                                                           SMALL_RESOLUTION,
                                                           RADIUS);

        assertTrue("expected the second combosFor call to hit the in-memory cache and return the same list instance",
                   first == second);
    }

    /**
     * The PRODUCTION {@code (geometryResolution=360, memberRadius=RADIUS)}
     * pair every real caller uses: {@link ContactComboCache#combosFor}
     * must resolve via the checked-in cache resource ({@link
     * ContactComboCache#RESOURCE_PATH}), not a live sweep (which would
     * make this test itself take several minutes - defeating the whole
     * point of the cache). A non-empty, plausible-sized ({@code < 10,800},
     * matching the documented ~446 discovery) result is the load-bearing
     * assertion; exact count is pinned more tightly by {@code
     * CommittedContactAtlasTest} against the real committed atlas.
     */
    @Test
    public void cachedCombosMatchTheResourceHeaderCount() {
        ContactPredicate predicate = newPredicate(ContactAtlasGenerator.GEOMETRY_RESOLUTION);
        long start = System.nanoTime();
        List<Combo> combos = ContactComboCache.combosFor(predicate,
                                                           ContactAtlasGenerator.GEOMETRY_RESOLUTION,
                                                           RADIUS);
        long elapsedMs = (System.nanoTime() - start) / 1_000_000;

        assertFalse("expected the production combo cache to be non-empty",
                    combos.isEmpty());
        assertTrue("expected a proper subset of the full 10,800-combo universe, found "
                   + combos.size(), combos.size() < 10_800);
        assertTrue("cached combosFor took " + elapsedMs
                   + "ms - expected a fast (<10s) resource load, not a live multi-minute sweep; "
                   + "is src/main/resources/lga/discovered-combos-cache.tsv missing or stale?",
                   elapsedMs < 10_000);

        Set<Combo> distinct = new HashSet<>(combos);
        assertEquals("expected no duplicate combos in the discovered set",
                     combos.size(), distinct.size());
    }

    /**
     * {@link FccNeighborhood#DIRECTIONS} lists 6 positive directions
     * followed by their 6 negations; {@link ContactPredicate}'s own
     * direction-reversal symmetry (verified directly by {@code
     * ContactAtlasTest.atlasIsSymmetricUnderDirectionReversal}) means every
     * discovered positive-direction combo's A/B-swapped, opposite-direction
     * mirror is independently discoverable too.
     *
     * <p>EXHAUSTIVE over every positive-direction combo (final-review
     * Significant C, T2 critique-final-0nx21-automaton-arc.md [21949]):
     * previously spot-checked a random 50/223 sample. The combo universe
     * is small (223 canonical combos, cheap boolean/set-membership
     * compares - not a 360x360 geometric sweep), so there is no cost
     * reason to sample rather than cover every combo; sampling only left
     * ~173 combos' mirror-closure unverified for no benefit.
     */
    @Test
    public void discoveredCombosAreClosedUnderDirectionReversalMirroring() {
        ContactPredicate predicate = newPredicate(ContactAtlasGenerator.GEOMETRY_RESOLUTION);
        List<Combo> combos = ContactComboCache.combosFor(predicate,
                                                           ContactAtlasGenerator.GEOMETRY_RESOLUTION,
                                                           RADIUS);
        Set<Combo> comboSet = new HashSet<>(combos);
        List<Combo> positive = combos.stream()
                                      .filter(c -> c.direction() > 0)
                                      .toList();
        assertFalse(positive.isEmpty());

        int checked = 0;
        for (Combo combo : positive) {
            Combo mirror = new Combo(FccNeighborhood.opposite(combo.direction()),
                                      combo.cubeB(), combo.memberB(),
                                      combo.cubeA(), combo.memberA());
            assertTrue("expected mirror " + mirror + " of discovered combo "
                       + combo + " to also be discovered",
                       comboSet.contains(mirror));
            checked++;
        }
        assertEquals("expected every positive-direction combo to be checked, not a subset",
                     positive.size(), checked);
        assertTrue(checked > 0);
    }

    /**
     * OPT-IN, slow (several-minutes) live-drift tripwire (code review
     * follow-up on bead inviscid-gyt): {@link ContactComboCache#combosFor}'s
     * cache-validity check is a {@code (geometryResolution, memberRadius)}
     * header match plus an internal {@code comboCount} cross-check only -
     * neither catches a {@link ContactPredicate}/{@link MemberGeometry}/
     * {@link com.chiralbehaviors.inviscid.PhiCoordinates} algorithm change
     * that leaves those two parameters unchanged but silently changes
     * which combos actually contact. This test regenerates the combo
     * universe LIVE at the PRODUCTION resolution (360) via {@link
     * ContactComboCache#sweepExhaustively} and diffs it against the
     * checked-in {@code discovered-combos-cache.tsv} resource - any
     * disagreement means the committed cache has silently gone stale.
     *
     * <p><b>Why this reads the RESOURCE and not {@link
     * ContactComboCache#combosFor}.</b> The subject of this tripwire is
     * THE FILE. It used to reach the file through {@code combosFor},
     * which was safe only while the in-JVM memo could be populated from
     * nowhere else. {@link ContactComboCache#rebuild} now publishes its
     * SWEEP into that memo (bead inviscid-0nx.30, E.3), so a single
     * earlier {@code PerRadiusRegeneration.regenerate} at the committed
     * {@code (360, RADIUS)} pair in the same JVM would have made this
     * assertion sweep-vs-sweep: unconditionally green, and green precisely
     * when the committed cache was most likely to be stale. No caller does
     * that today; E.4/E.5 regenerating a baseline for side-by-side
     * comparison is exactly the code that would. {@link
     * ContactComboCache#loadCommittedCache} bypasses the memo, so this
     * comparison is file-vs-sweep by construction rather than by luck of
     * what ran first.
     *
     * <p>Gated behind the {@code -Dinviscid.slowTests=true} system
     * property (skipped otherwise via {@link Assume#assumeTrue}, so
     * {@code mvn test} without the property reports this test SKIPPED,
     * not run) - a live {@code 360x360} exhaustive sweep over the full
     * 10,800-combo universe costs several minutes wall time (see class
     * Javadoc), too slow for the ordinary fast test suite. Run explicitly
     * with {@code mvn test -Dtest=ContactComboCacheTest
     * -Dinviscid.slowTests=true}.
     */
    @Test
    public void liveSweepMatchesTheCommittedCacheAtProductionResolution() {
        Assume.assumeTrue("skipped unless -Dinviscid.slowTests=true is set (several-minutes live sweep - see javadoc)",
                           Boolean.getBoolean("inviscid.slowTests"));

        int resolution = ContactAtlasGenerator.GEOMETRY_RESOLUTION;

        List<Combo> cached = ContactComboCache.loadCommittedCache(resolution,
                                                                    RADIUS);
        assertFalse("the committed discovered-combos-cache.tsv must be present and its header must declare geometryResolution="
                    + resolution + ", memberRadius=" + RADIUS
                    + " - a missing or foreign-header resource would make this tripwire compare against nothing",
                    cached == null);
        List<Combo> live = ContactComboCache.sweepExhaustively(newPredicate(resolution),
                                                                 resolution);

        Set<Combo> cachedSet = new HashSet<>(cached);
        Set<Combo> liveSet = new HashSet<>(live);

        Set<Combo> missingFromCache = new HashSet<>(liveSet);
        missingFromCache.removeAll(cachedSet);
        Set<Combo> staleInCache = new HashSet<>(cachedSet);
        staleInCache.removeAll(liveSet);

        assertTrue("committed discovered-combos-cache.tsv has drifted from a live sweep at geometryResolution="
                   + resolution + ": " + missingFromCache.size()
                   + " combos newly contact but are missing from the cache, "
                   + staleInCache.size()
                   + " cached combos no longer contact live - regenerate via ContactComboCache.main. "
                   + "missingFromCache=" + missingFromCache + " staleInCache=" + staleInCache,
                   missingFromCache.isEmpty() && staleInCache.isEmpty());
    }
}
