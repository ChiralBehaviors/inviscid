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
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

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

    private static final int    GEOMETRY_RESOLUTION = ContactAtlasGenerator.GEOMETRY_RESOLUTION;
    private static final double RADIUS              = ContactAtlasGenerator.RADIUS;
    private static final String RESOURCE_PATH       = "lga/contact-atlas-v2.tsv";

    private static ContactPredicate    PREDICATE;
    private static FineStepContactTable FINE;
    private static ContactAtlas        ATLAS;
    private static ContactTable        BIN_TABLE;

    @BeforeClass
    public static void buildFixtures() throws IOException {
        PREDICATE = new ContactPredicate(new MemberGeometry(GEOMETRY_RESOLUTION,
                                                              RADIUS));
        long start = System.nanoTime();
        FINE = FineStepContactTable.buildFor(PREDICATE, GEOMETRY_RESOLUTION,
                                              RADIUS);
        long elapsedMs = (System.nanoTime() - start) / 1_000_000;
        System.out.println("FineStepContactTable build: " + elapsedMs + "ms");

        URL resource = FineStepContactTableTest.class.getClassLoader()
                                                        .getResource(RESOURCE_PATH);
        ATLAS = ContactAtlas.read(Paths.get(resource.getPath()));
        BIN_TABLE = ContactTable.of(ATLAS);
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
     * Memory footprint pin. Canonical-direction-only storage (class
     * Javadoc, "Canonical directions only") roughly HALVES the naive
     * "all 446 combos" estimate (~6.9 MiB) to ~3.5 MiB, since the 446
     * discovered combos split roughly evenly between positive and
     * negative directions and only the positive half is ever looked up.
     */
    @Test
    public void memoryFootprintIsPinnedInTheExpectedRange() {
        long bytes = FINE.memoryFootprintBytes();
        double mib = bytes / (1024.0 * 1024.0);
        assertTrue("expected ~3.5 MiB (canonical-only storage), was " + mib
                   + " MiB", mib > 2.5 && mib < 5.0);
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
        List<ContactComboCache.Combo> combos = ContactComboCache.combosFor(PREDICATE,
                                                                              GEOMETRY_RESOLUTION,
                                                                              RADIUS);
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
        // ~223 canonical (positive-direction) combos x 360x360 -
        // roughly half of the naive "all 446 combos" 57.8M estimate,
        // since only positive-direction combos are ever stored/queried
        // (class Javadoc, "Canonical directions only").
        assertTrue("expected a substantial exhaustive check count, was "
                   + checked, checked > 20_000_000L);
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
        List<ContactComboCache.Combo> combos = ContactComboCache.combosFor(PREDICATE,
                                                                              GEOMETRY_RESOLUTION,
                                                                              RADIUS);
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
