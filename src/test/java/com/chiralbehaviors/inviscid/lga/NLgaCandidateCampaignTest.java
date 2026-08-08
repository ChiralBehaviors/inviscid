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

import org.junit.Test;

/**
 * Fast, deterministic unit tests for {@link NLgaCandidateCampaign#circularWidth}
 * - the pure-math seam the full campaign's per-combo width measurement is
 * built on. The campaign itself (bead inviscid-0nx.16's measurement
 * campaign - real geometric sweeps plus a real dynamic run per candidate)
 * is a main-able harness run outside surefire (see that class's own
 * Javadoc); this test only pins down the width-extraction arithmetic so a
 * regression there is caught fast, without paying the campaign's own cost
 * in the ordinary build.
 *
 * @author halhildebrand
 */
public class NLgaCandidateCampaignTest {

    private static final double TOLERANCE = 1e-9;

    @Test
    public void noMarksIsZeroWidth() {
        boolean[] marks = new boolean[360];
        assertEquals(0.0,
                     NLgaCandidateCampaign.circularWidth(marks, 360),
                     TOLERANCE);
    }

    @Test
    public void everyMarkIsFullCircle() {
        boolean[] marks = new boolean[360];
        java.util.Arrays.fill(marks, true);
        assertEquals(2 * Math.PI,
                     NLgaCandidateCampaign.circularWidth(marks, 360),
                     TOLERANCE);
    }

    @Test
    public void singleMarkIsZeroWidth() {
        boolean[] marks = new boolean[360];
        marks[42] = true;
        assertEquals(0.0,
                     NLgaCandidateCampaign.circularWidth(marks, 360),
                     TOLERANCE);
    }

    /**
     * A contiguous, non-wrapping block of marks: indices [10, 20)
     * (inclusive of both endpoints 10 and 19) span 9 steps at
     * {@code 2*pi/360} radians/step.
     */
    @Test
    public void contiguousBlockMeasuresSpanBetweenExtremes() {
        boolean[] marks = new boolean[360];
        for (int i = 10; i < 20; i++) {
            marks[i] = true;
        }
        double expected = 9 * (2 * Math.PI / 360);
        assertEquals(expected,
                     NLgaCandidateCampaign.circularWidth(marks, 360),
                     TOLERANCE);
    }

    /**
     * A block that wraps across the 0/2*pi boundary (indices 355..359 and
     * 0..4) must be measured as the SHORT arc through the wrap (9 steps),
     * not the long way around - the entire reason {@link
     * NLgaCandidateCampaign#circularWidth} uses the largest-gap
     * complement instead of a naive max-minus-min.
     */
    @Test
    public void wrappingBlockMeasuresShortArcThroughTheBoundary() {
        boolean[] marks = new boolean[360];
        for (int i = 355; i < 360; i++) {
            marks[i] = true;
        }
        for (int i = 0; i <= 4; i++) {
            marks[i] = true;
        }
        double expected = 9 * (2 * Math.PI / 360);
        assertEquals(expected,
                     NLgaCandidateCampaign.circularWidth(marks, 360),
                     TOLERANCE);
    }

    /**
     * Two well-separated single-step blobs: the largest gap is between
     * them (not the wrap gap), so the measured "width" is the smaller
     * span from the first blob through the second - exercising the
     * general (non-contiguous) case the largest-gap-complement formula
     * must also handle.
     */
    @Test
    public void twoBlobsMeasuresSpanExcludingTheLargestGap() {
        boolean[] marks = new boolean[360];
        marks[0] = true;
        marks[90] = true;
        double expected = 90 * (2 * Math.PI / 360);
        assertEquals(expected,
                     NLgaCandidateCampaign.circularWidth(marks, 360),
                     TOLERANCE);
    }
}
