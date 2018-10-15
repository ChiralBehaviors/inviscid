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
    private static final Material[] edgeMaterials = new PhongMaterial[] { blueMaterial,
                                                                          blueMaterial,
                                                                          blueMaterial,
                                                                          blueMaterial,
                                                                          blueMaterial,
                                                                          blueMaterial };

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    protected void initializeContentModel() {
        CubicGrid grid = new CubicGrid(Neighborhood.SIX,
                                       PhiCoordinates.Cubes[3], 0);

        ContentModel content = getContentModel();
        Group group = new Group();
        content.setContent(group);
        Necronomata automata = new Necronomata(new Point3i(5, 5, 5));
        NecronomataVisualization visualization = new NecronomataVisualization(3600,
                                                                              (float) 0.015,
                                                                              automata,
                                                                              edgeMaterials);
        group.getChildren()
             .add(visualization);

        group.getChildren()
             .add(grid.construct(blackMaterial, blackMaterial, blackMaterial));

        float[] delta = new float[6];
        float[] minusDelta = new float[6];
        for (int i = 0; i < 6; i++) {
            delta[i] = visualization.getAngularResolution();
            minusDelta[i] = -visualization.getAngularResolution();
        }
        automata.drive(NecronomataVisualization.getPositiveTet());
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
                int nextIndex = (int) (Math.toRadians(value)
                                       / visualization.getAngularResolution());
                int currentIndex = lastIndex;
                if (currentIndex != nextIndex) {
                    int deltaIndex = nextIndex - currentIndex;
                    lastIndex = nextIndex;
                    float[] apply = deltaIndex < 0 ? minusDelta : delta;
                    for (int step = 0; step < Math.abs(deltaIndex); step++) {
                        automata.drive(apply);
                        automata.step();
                        visualization.update();
                    }
                }
            }
        }, 45f);
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(1_000), keyValue));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
