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

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The fine-step ({@code geometryResolution x geometryResolution} per combo)
 * contact structure (USER DECISION 2026-08-08, FINAL, recorded on bead
 * inviscid-0nx.21): replicates the hybrid's {@link ContactPredicate}
 * point-evaluation semantics at the atlas's own {@code geometryResolution}
 * (360), for every one of {@link ContactComboCache}'s exhaustively-
 * discovered ever-contacting combos.
 *
 * <h2>Why this exists (test 8's finding)</h2>
 * The bin-level {@link ContactTable} (bead inviscid-0nx.19) transcribes
 * geometry via ANY-OVERLAP semantics (bead inviscid-gyt): a bin fires if
 * ANY fraction of its fine geometric sub-sweep contacts - a strict superset
 * of a live point evaluation, deliberately widened to fix a ~12% true-
 * contact reproduction rate under bin-centre-only transcription. That
 * widening, measured empirically for the first time once {@code
 * LatticeGasAutomaton} existed to drive it dynamically, produced a ~6.6x
 * collision-rate / ~3.9x effective-ratio gap against the hybrid (recorded
 * in {@code HybridVsLgaConsistencyTest}'s Javadoc as the quantified
 * motivation for this class). This class is the fix: contact FIRING moves
 * to the atlas's fine 360-step geometry grid, replicating point-evaluation
 * at that resolution. gyt and D1 are NOT overturned - {@link ContactTable}
 * and {@link CollisionTable} stand unchanged as the key/outcome layer;
 * any-overlap remains the BIN table's own transcription semantics, just no
 * longer the LGA's contact-firing decision.
 *
 * <h2>Build vs. artifact (measured, not guessed)</h2>
 * Built at construction, NOT committed as an artifact: measured 4.3-4.4s
 * wall, cold and warm alike (446 combos x 360x360 = 57,801,600 {@link
 * ContactPredicate#contacts} evaluations), comfortably under the ~10s
 * decision threshold. Cached per {@code (geometryResolution, memberRadius)}
 * pair for the JVM's lifetime (mirrors {@link ContactComboCache}'s own
 * in-memory cache), so repeated {@link LatticeGasAutomaton} construction
 * across a test suite pays the cost once.
 *
 * <h2>Reuse, not reimplementation</h2>
 * Every contact bit is {@link ContactPredicate#contacts} evaluated at
 * {@link ContactComboCache#angleOf} step angles - the EXACT same geometry
 * call and step-to-angle reconstruction {@link
 * ContactAtlasGenerator#sweepComboOverlap} uses to build the committed
 * atlas, applied to every step pair instead of aggregated into bins. The
 * combo universe itself is {@link ContactComboCache#combosFor} (the
 * checked-in 446-combo cache) - a combo NOT in that set never contacts at
 * any step, by that class's own exhaustiveness proof, so it costs one
 * array lookup (a "not found" sentinel), never a geometry re-check.
 *
 * <h2>Canonical directions only</h2>
 * Only combos whose {@code direction} is in {@code 1..6} are stored -
 * {@link LatticeGasAutomaton}'s scan phase only ever queries canonical
 * (positive) directions (matching {@link ContactScan}'s own
 * canonicalization), so a combo discovered only under a negative direction
 * is dropped: it is never looked up.
 *
 * <h2>No bin-level pre-filter (deliberate, reported)</h2>
 * The bin-level {@link ContactTable} was considered as a coarse pre-filter
 * ahead of the fine lookup. Rejected: this class's lookup is already O(1)
 * (one array index, one bit test, no allocation) - as cheap as a pre-filter
 * check would be - and since ANY-OVERLAP is a strict SUPERSET of the fine
 * predicate, a bin-level "no" already implies a fine "no" but a bin-level
 * "yes" (the common case, given the measured ~6.6x over-firing) almost
 * never lets the pre-filter skip the fine check. Adding it would cost a
 * bin-index computation and a second table lookup for negligible skip
 * benefit.
 *
 * <h2>Memory</h2>
 * One packed bitset per stored combo, {@code geometryResolution^2} bits
 * (129,600 bits = 16,200 bytes at 360). Measured total ~3.5 MiB (~223
 * canonical-direction combos - roughly half of the naive "all 446" ~6.9
 * MiB estimate, since only positive-direction combos are stored) - see
 * {@link #memoryFootprintBytes()} and its test pin.
 *
 * @author halhildebrand
 */
final class FineStepContactTable {

    private static final int CUBES_PER_CELL   = 5;
    private static final int MEMBERS_PER_CUBE = 6;
    private static final int CANONICAL_DIRECTIONS = 6;

    private record CacheKey(int geometryResolution, double memberRadius) {
    }

    private static final Map<CacheKey, FineStepContactTable> CACHE = new ConcurrentHashMap<>();

    private final int      geometryResolution;
    private final int[]    comboIndexByIdentity;
    private final long[][] bitsByCombo;
    private final int      wordsPerCombo;

    /**
     * @return the (cached, built-once-per-geometry) fine-step contact
     *         table for {@code predicate}'s geometry.
     */
    static FineStepContactTable buildFor(ContactPredicate predicate,
                                          int geometryResolution,
                                          double memberRadius) {
        return CACHE.computeIfAbsent(new CacheKey(geometryResolution,
                                                    memberRadius),
                                      key -> build(predicate, key));
    }

    private static FineStepContactTable build(ContactPredicate predicate,
                                                CacheKey key) {
        int resolution = key.geometryResolution();
        List<ContactComboCache.Combo> discovered = ContactComboCache.combosFor(predicate,
                                                                                  resolution,
                                                                                  key.memberRadius());
        List<ContactComboCache.Combo> canonical = discovered.stream()
                                                              .filter(c -> c.direction() >= 1
                                                                           && c.direction() <= CANONICAL_DIRECTIONS)
                                                              .toList();

        int identityDomain = CANONICAL_DIRECTIONS * CUBES_PER_CELL
                              * MEMBERS_PER_CUBE * CUBES_PER_CELL
                              * MEMBERS_PER_CUBE;
        int[] comboIndexByIdentity = new int[identityDomain];
        Arrays.fill(comboIndexByIdentity, -1);

        long domainBits = (long) resolution * resolution;
        int wordsPerCombo = (int) ((domainBits + 63) / 64);
        long[][] bitsByCombo = new long[canonical.size()][wordsPerCombo];

        for (int i = 0; i < canonical.size(); i++) {
            ContactComboCache.Combo combo = canonical.get(i);
            comboIndexByIdentity[identityIndex(combo.direction(),
                                                combo.cubeA(), combo.memberA(),
                                                combo.cubeB(),
                                                combo.memberB())] = i;
            long[] bits = bitsByCombo[i];
            for (int a = 0; a < resolution; a++) {
                float angleA = ContactComboCache.angleOf(a, resolution);
                for (int b = 0; b < resolution; b++) {
                    float angleB = ContactComboCache.angleOf(b, resolution);
                    if (predicate.contacts(combo.cubeA(), combo.memberA(),
                                            angleA, combo.cubeB(),
                                            combo.memberB(), angleB,
                                            combo.direction())) {
                        long idx = (long) a * resolution + b;
                        bits[(int) (idx >>> 6)] |= 1L << (idx & 63);
                    }
                }
            }
        }

        return new FineStepContactTable(resolution, comboIndexByIdentity,
                                         bitsByCombo, wordsPerCombo);
    }

    private FineStepContactTable(int geometryResolution,
                                  int[] comboIndexByIdentity,
                                  long[][] bitsByCombo, int wordsPerCombo) {
        this.geometryResolution = geometryResolution;
        this.comboIndexByIdentity = comboIndexByIdentity;
        this.bitsByCombo = bitsByCombo;
        this.wordsPerCombo = wordsPerCombo;
    }

    /**
     * The hot-path lookup: one array index to resolve the combo, one bit
     * test - no allocation, no geometry, no trig.
     *
     * @param direction canonical (positive) FCC direction, {@code 1..6}
     */
    boolean contacts(int direction, int cubeA, int memberA, int stepA,
                      int cubeB, int memberB, int stepB) {
        int comboIndex = comboIndexByIdentity[identityIndex(direction, cubeA,
                                                              memberA, cubeB,
                                                              memberB)];
        if (comboIndex < 0) {
            return false;
        }
        long idx = (long) stepA * geometryResolution + stepB;
        return (bitsByCombo[comboIndex][(int) (idx >>> 6)] & (1L << (idx
                                                                       & 63))) != 0L;
    }

    int geometryResolution() {
        return geometryResolution;
    }

    int comboCount() {
        return bitsByCombo.length;
    }

    /** Total heap held by the packed bitsets + the identity index array. */
    long memoryFootprintBytes() {
        return (long) bitsByCombo.length * wordsPerCombo * 8L
               + (long) comboIndexByIdentity.length * 4L;
    }

    private static int identityIndex(int direction, int cubeA, int memberA,
                                      int cubeB, int memberB) {
        int d = direction - 1;
        return (((d * CUBES_PER_CELL + cubeA) * MEMBERS_PER_CUBE + memberA)
                * CUBES_PER_CELL + cubeB) * MEMBERS_PER_CUBE + memberB;
    }
}
