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

import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.regex.Pattern;

import org.junit.Test;

/**
 * Guards the DOCUMENTED EXCEPTIONS to the {@link LgaTestGeometry}
 * consolidation (bead inviscid-0nx.30, E.3).
 *
 * <p>
 * The consolidation's value is that a reader who sees {@code
 * LgaTestGeometry.BASELINE_RADIUS} knows the value tracks the campaign
 * anchor. Its risk is the mirror image: a reader who sees a BARE {@code
 * 0.015} left behind cannot tell whether it is an oversight from the
 * migration or a deliberate exception, and the tidy-up that "finishes the
 * job" is in one case exactly wrong. {@code measure.SeamGoldenCompatTest}
 * is that case, and it is a PINNED seam this bead may not edit - so the
 * protection has to live here.
 *
 * @author halhildebrand
 */
public class LgaTestGeometryTest {

    private static final Path GOLDEN_SEAM_SOURCE = Path.of("src/test/java/com/chiralbehaviors/inviscid/measure/SeamGoldenCompatTest.java");

    /**
     * {@link LgaTestGeometry#GOLDEN_PHASE_A_RADIUS} must be a FROZEN
     * literal, never an alias. If it were sourced from {@link
     * LgaTestGeometry#BASELINE_RADIUS} (or from {@code
     * ContactAtlasGenerator.RADIUS}, which is the same thing) it would move
     * with the anchor and stop denoting the value Phase A's goldens were
     * captured at, which is the only thing it is for.
     *
     * <p>
     * FALSIFIER: rewriting that constant as an alias, or editing its value.
     */
    @Test
    public void thePhaseAGoldenRadiusIsFrozenAtTheValueTheGoldensWereCapturedAt() {
        assertEquals("Phase A's goldens were captured at r=0.015; this constant denotes THAT, not the current anchor",
                     0.015, LgaTestGeometry.GOLDEN_PHASE_A_RADIUS, 0.0);
    }

    /**
     * THE BEHAVIOURAL HALF of a two-part guard: the seam's {@code RADIUS}
     * field really does hold the frozen Phase A value, read by reflection
     * off the compiled class rather than inferred from source text. It is
     * not the primary check and does not subsume {@link
     * #phaseAGoldenSeamKeepsItsOwnFrozenRadius} - see below for the case
     * only the source check covers.
     *
     * <p>
     * This is what {@link #phaseAGoldenSeamKeepsItsOwnFrozenRadius}'s
     * source regex was standing in for, and it is strictly better at the
     * job: it is immune to formatting, to comments, and to the regex's
     * prefix-looseness (an un-terminated {@code RADIUS\s*=\s*0.015} match
     * accepts {@code 0.0150} and {@code 0.0155} as happily as {@code
     * 0.015}). If the field is ever consolidated onto a moved anchor, the
     * value read here changes and this goes red naming both numbers.
     *
     * <p>
     * It reads the FIELD, so it deliberately does not cover the case where
     * the field is left alone and {@code newGeometry()} is rewired to
     * {@link LgaTestGeometry#BASELINE_RADIUS} directly, leaving {@code
     * RADIUS} unused. That is the residual case the source check still
     * carries; the two are complements, not duplicates.
     *
     * <p>
     * Note that a {@code static final double} initialized from a literal
     * is a compile-time {@code ConstantValue}, so this reads the value the
     * class file records for the field - the declared value - not a
     * runtime-computed one.
     *
     * <p>
     * FALSIFIER: editing that field's value, or sourcing it from anything
     * that can move.
     */
    @Test
    public void phaseAGoldenSeamsRadiusFieldHoldsTheFrozenValue() throws Exception {
        Field radius = Class.forName("com.chiralbehaviors.inviscid.measure.SeamGoldenCompatTest")
                             .getDeclaredField("RADIUS");
        radius.setAccessible(true);
        assertEquals("SeamGoldenCompatTest.RADIUS must still be the FROZEN Phase A golden radius, not the current campaign anchor ("
                     + LgaTestGeometry.BASELINE_RADIUS + ")",
                     LgaTestGeometry.GOLDEN_PHASE_A_RADIUS,
                     radius.getDouble(null), 0.0);
    }

    /**
     * THE EXCEPTION, MECHANIZED. {@code SeamGoldenCompatTest} must keep its
     * OWN radius literal and must not be consolidated onto the shared
     * anchor.
     *
     * <h2>What the failure actually looks like</h2>
     * Not a silent green - an earlier draft of this Javadoc claimed that,
     * and it was wrong. That file reads {@code RADIUS} at exactly ONE site
     * ({@code newGeometry()}), feeding exactly one test, whose
     * radius-sensitive quantities are HARD LITERALS ({@code
     * LEDGER_TOTAL_COLLISIONS = 40L}, {@code LEDGER_EFFECTIVE_COLLISIONS =
     * 26L}). Nothing there is regenerated as a comparison input, so there
     * is no self-comparison to go green: consolidating the file and then
     * retargeting the anchor to {@code 0.05}, ON A CLEAN BUILD, was
     * measured to produce {@code expected:<40> but was:<162>}.
     *
     * <p>
     * The clean build is not incidental, and reproducing this without one
     * is a trap that has already caught a first run: incrementally, the
     * ledger test stays GREEN, because {@link
     * LgaTestGeometry#BASELINE_RADIUS} is a compile-time {@code
     * ConstantValue} still inlined in the unrecompiled consumer. See
     * {@link LgaTestGeometry#GOLDEN_PHASE_A_RADIUS} for the full account -
     * including WHICH of that seam's pins this exception actually governs,
     * which is fewer than all of them.
     *
     * <p>
     * The failure is LOUD but UNATTRIBUTABLE, and that is what earns a
     * source-reading test. It arrives at the moment of the ANCHOR MOVE,
     * not the moment of the edit; it arrives alongside every other
     * consequence of that move; and its message names a collision count,
     * with nothing connecting it to a radius consolidation made weeks
     * earlier. This test converts it into a named red at the point of the
     * edit that caused it. The file is a PINNED seam, so no in-file
     * mechanism THAT RUNS - no assertion, no guard of its own - can be
     * added there; the guard has to live on this side. A COMMENT is a
     * different matter, and the exclusion below is deliberately shaped to
     * leave room for one.
     *
     * <h2>Why the exclusion is comment-blind</h2>
     * The narrow reading is on CODE, not on the raw file. A previous
     * spelling asserted the file does not {@code contains("LgaTestGeometry")}
     * at all, which forbade the single most useful thing anyone could add
     * to it - a comment saying "deliberately does NOT use LgaTestGeometry"
     * - and so foreclosed the call-site-discoverability problem it was
     * meant to sit alongside. Comments are stripped first; only a real
     * code reference fails.
     *
     * <p>
     * The check remains intentionally narrow: the file exists, it still
     * declares its own terminated {@code 0.015} literal, and its code does
     * not reference {@link LgaTestGeometry}. It says nothing about the
     * rest of that file, which remains a pinned seam.
     *
     * <p>
     * FALSIFIER: any edit replacing that file's local literal with a
     * reference to the shared constant, or reaching for that constant
     * anywhere else in its code - i.e. exactly the well-intentioned
     * tidy-up this exception exists to forbid.
     */
    @Test
    public void phaseAGoldenSeamKeepsItsOwnFrozenRadius() throws IOException {
        assertTrue("the pinned golden seam must exist for this guard to mean anything: "
                   + GOLDEN_SEAM_SOURCE, Files.isRegularFile(GOLDEN_SEAM_SOURCE));
        String source = Files.readString(GOLDEN_SEAM_SOURCE,
                                          StandardCharsets.UTF_8);

        assertTrue("SeamGoldenCompatTest must keep its own frozen radius literal ("
                   + LgaTestGeometry.GOLDEN_PHASE_A_RADIUS
                   + "): it pins Phase A goldens by contract and must NOT follow a moving campaign anchor",
                   Pattern.compile("RADIUS\\s*=\\s*"
                                    + Pattern.quote(Double.toString(LgaTestGeometry.GOLDEN_PHASE_A_RADIUS))
                                    + "\\s*;")
                           .matcher(source)
                           .find());

        assertFalse("SeamGoldenCompatTest's CODE must not reference LgaTestGeometry - consolidating it would retarget a pinned Phase A golden onto a moving anchor, surfacing later as an unattributable collision-count mismatch. A COMMENT mentioning LgaTestGeometry is fine, and encouraged.",
                    withoutComments(source).contains("LgaTestGeometry"));
    }

    /**
     * @return {@code source} with block and line comments removed, so a
     *         code-reference assertion cannot be tripped by prose. Good
     *         enough for this one pinned file, and not more: it would also
     *         strip a {@code //} or a comment opener occurring inside a
     *         string literal, of which that file has none. If it ever
     *         acquires one, this method needs a real lexer, not a wider
     *         regex.
     */
    private static String withoutComments(String source) {
        return source.replaceAll("(?s)/\\*.*?\\*/", "")
                      .replaceAll("(?m)//[^\n]*", "");
    }
}
