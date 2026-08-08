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

package com.chiralbehaviors.inviscid.animations;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.LengthTable;
import com.chiralbehaviors.inviscid.Necronomata;

import javafx.scene.transform.Scale;
import javafx.scene.transform.Transform;

/**
 * Regression test for inviscid-6cf: {@code
 * NecronomataVisualization.buildLengths} must render {@link
 * LengthTable#lengthAt(int)} directly for every step, with no per-step
 * override. Constructs the real {@code NecronomataVisualization} — no
 * {@code Platform.startup} / toolkit initialization needed, since {@code
 * javafx.scene.transform.*} and {@code javafx.scene.shape.MeshView} /
 * {@code javafx.scene.Group} construction (as opposed to display) does
 * not require a running JavaFX application thread; this runs directly
 * under surefire, exactly like {@code MemberGeometryJavaFxParityTest}.
 *
 * <p>Filed as inviscid-ghp (substantive-critic finding on the inviscid-6cf
 * fix): the 6cf change removed {@code buildLengths}'s override but shipped
 * with no test touching the production method itself — {@code
 * MemberGeometryJavaFxParityTest} only exercises the headless {@code
 * MemberGeometry} mirror, never {@code NecronomataVisualization}. This
 * test closes that gap by inspecting the visualization's actual {@code
 * lengths} LUT via the package-private {@link
 * NecronomataVisualization#lengths()} test-accessor.
 */
public class NecronomataVisualizationLengthsTest {

    private static final double DELTA      = 0.0;
    private static final int    RESOLUTION = 360;

    @Test
    public void lengthsMatchLengthTableAtEveryStep() {
        Necronomata automata = new Necronomata(new Point3i(4, 4, 4));
        NecronomataVisualization visualization = new NecronomataVisualization(RESOLUTION,
                                                                               0.015f,
                                                                               automata,
                                                                               Colors.materials);
        LengthTable table = new LengthTable(RESOLUTION);
        Transform[] lengths = visualization.lengths();

        assertEquals(RESOLUTION, lengths.length);
        for (int step = 0; step < RESOLUTION; step++) {
            double expected = table.lengthAt(step);
            double actual = ((Scale) lengths[step]).getY();
            assertEquals("step " + step, expected, actual, DELTA);
        }
    }

    /**
     * Explicitly pins the four previously-overridden indices (1, 2, 3
     * forced to index 0's value; 5 forced to index 4's value): asserts
     * both that each equals {@code LengthTable}'s true value AND that it
     * is NOT equal to the old override target, so a reintroduced override
     * fails this test loudly rather than only failing the general sweep
     * above.
     */
    @Test
    public void previouslyOverriddenStepsAreNoLongerFlattened() {
        Necronomata automata = new Necronomata(new Point3i(4, 4, 4));
        NecronomataVisualization visualization = new NecronomataVisualization(RESOLUTION,
                                                                               0.015f,
                                                                               automata,
                                                                               Colors.materials);
        LengthTable table = new LengthTable(RESOLUTION);
        Transform[] lengths = visualization.lengths();

        double step0 = ((Scale) lengths[0]).getY();
        double step4 = ((Scale) lengths[4]).getY();

        for (int step : new int[] { 1, 2, 3 }) {
            double actual = ((Scale) lengths[step]).getY();
            assertEquals("step " + step, table.lengthAt(step), actual, DELTA);
            assertNotEquals("step " + step
                            + " must not be flattened to step 0's value",
                            step0, actual, DELTA);
        }

        double actual5 = ((Scale) lengths[5]).getY();
        assertEquals("step 5", table.lengthAt(5), actual5, DELTA);
        assertNotEquals("step 5 must not be flattened to step 4's value",
                        step4, actual5, DELTA);
    }
}
