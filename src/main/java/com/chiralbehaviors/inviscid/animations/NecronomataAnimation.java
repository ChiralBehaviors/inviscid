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

import static com.chiralbehaviors.inviscid.animations.Colors.blackMaterial;
import static com.chiralbehaviors.inviscid.animations.Colors.blueMaterial;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.CubicGrid;
import com.chiralbehaviors.inviscid.CubicGrid.Neighborhood;
import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.KeyFrame;
import javafx.animation.KeyValue;
import javafx.animation.Timeline;
import javafx.beans.value.WritableValue;
import javafx.scene.Group;
import javafx.scene.paint.Material;
import javafx.scene.paint.PhongMaterial;
import javafx.util.Duration;

/**
 * @author halhildebrand
 *
 */
public class NecronomataAnimation extends PolyView {
    public static class Launcher {

        public static void main(String[] argv) {
            NecronomataAnimation.main(argv);
        }
    }

    private static final Material[] edgeMaterials  = new PhongMaterial[] { blueMaterial, blueMaterial, blueMaterial,
                                                                            blueMaterial, blueMaterial, blueMaterial };

    /**
     * Alternating quanta pattern cycling over the 6 members of a cube, used
     * ONLY for the one-time initial kick before the timeline starts.
     */
    private static final float[]    POSITIVE_QUANTA = { 1, 1, -1, -1, 1, 1 };

    /** Negation of {@link #POSITIVE_QUANTA}, one-time initial kick only. */
    private static final float[]    NEGATIVE_QUANTA = { -1, -1, 1, 1, -1, -1 };

    /**
     * Uniform quanta magnitude 1 across all six members: the recurring
     * per-tick drive rotates every strut together (matches the pre-change
     * behavior, where the timeline's recurring drive fed a uniform
     * angular-resolution rate into every deltaA slot). Necronomata.QUANTUM_RATE
     * is EQUAL BY CONSTRUCTION to this visualization's angular
     * resolution: both derive from Necronomata.PHASE_RESOLUTION (the LUT
     * below is constructed with it), so a frequency of 1 advances a
     * member by exactly one LUT step per tick.
     */
    private static final float[]    UNIFORM_POSITIVE_QUANTA = { 1, 1, 1, 1, 1, 1 };

    /** Negation of {@link #UNIFORM_POSITIVE_QUANTA}, recurring per-tick drive only. */
    private static final float[]    UNIFORM_NEGATIVE_QUANTA = { -1, -1, -1, -1, -1, -1 };

    public static void main(String[] args) {
        launch(args);
    }

    /**
     * Seed frequency (the conserved quanta count) across every member,
     * cycling {@code pattern} over the 6 members of a cube. Replaces the
     * removed {@code Necronomata} drive(float[]) method, which erroneously
     * wrote directly into deltaA (an angular rate) rather than the
     * frequency quanta that now drive it via {@link Necronomata#QUANTUM_RATE}.
     */
    private static void seedFrequency(Necronomata automata, float[] pattern) {
        automata.process((angle, frequency, deltaA, deltaF) -> {
            for (int i = 0; i < frequency.length; i++) {
                frequency[i] = pattern[i % 6];
            }
        });
    }

    @Override
    protected void initializeContentModel() {
        CubicGrid grid = new CubicGrid(Neighborhood.SIX, PhiCoordinates.Cubes[3], 0);

        ContentModel content = getContentModel();
        Group group = new Group();
        content.setContent(group);
        // (4, 4, 4): smallest legal even extent on every axis, chosen to
        // satisfy com.chiralbehaviors.inviscid.lga.FccNeighborhood's
        // even-parity/PBC preconditions (all-even, >= 4 per axis —
        // inviscid-0nx.3) up front. Necronomata.step() does not yet wire
        // collision rules into process(Point3i) (reserved to
        // inviscid-0nx.14), so this visualization does not depend on
        // FccNeighborhood/wrap-around neighbors today — but keeping the
        // extent PBC-legal now means the animation stays valid without
        // rework once collision dynamics land (inviscid-7c4).
        Necronomata automata = new Necronomata(new Point3i(4, 4, 4));
        NecronomataVisualization visualization = new NecronomataVisualization(Necronomata.PHASE_RESOLUTION,
                (float) 0.015, automata,
                edgeMaterials);
        group.getChildren().add(visualization);

        group.getChildren().add(grid.construct(blackMaterial, blackMaterial, blackMaterial));

        seedFrequency(automata, POSITIVE_QUANTA);
        automata.step();
        visualization.update();
        final Timeline timeline = new Timeline();
        KeyValue keyValue = new KeyValue(new WritableValue<Float>() {
            volatile int lastIndex = 0;

            @Override
            public Float getValue() {
                return 0f;
            }

            @Override
            public void setValue(Float value) {
                int nextIndex = (int) (Math.toRadians(value) / visualization.getAngularResolution());
                int currentIndex = lastIndex;
                if (currentIndex != nextIndex) {
                    int deltaIndex = nextIndex - currentIndex;
                    lastIndex = nextIndex;
                    float[] apply = deltaIndex < 0 ? UNIFORM_NEGATIVE_QUANTA : UNIFORM_POSITIVE_QUANTA;
                    for (int step = 0; step < Math.abs(deltaIndex); step++) {
                        seedFrequency(automata, apply);
                        automata.step();
                        visualization.update();
                    }
                }
            }
        }, 45f);
        timeline.getKeyFrames().add(new KeyFrame(Duration.millis(1_000), keyValue));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
