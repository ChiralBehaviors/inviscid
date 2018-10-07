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
public class AutomataAnimation extends PolyView {
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
                                       PhiCoordinates.Cubes[3], 2);
        
        ContentModel content = getContentModel();
        Group group = new Group();
        content.setContent(group);
        Necronomata automata = new Necronomata(new Point3i(5, 5, 5));
        AutomataVisualization a = new AutomataVisualization(360, (float) 0.015,
                                                            automata,
                                                            edgeMaterials);
        group.getChildren()
             .addAll(a.getNodes());

        float[] plus = AutomataVisualization.getPositiveTet();
        float[] newPlus = new float[30];
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 6; j++) {
                newPlus[i * 6 + j] = plus[j];
            }
        }
        for (int i = 0; i < a.getNodes().length; i++) {
            a.getNodes()[i].setState(newPlus);
        }


        group.getChildren()
             .add(grid.construct(blackMaterial, blackMaterial, blackMaterial));
        final Timeline timeline = new Timeline();
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(2_500),
                                  new KeyValue(new WritableValue<Float>() {
                                      @Override
                                      public Float getValue() {
                                          return 0f;
                                      }

                                      @Override
                                      public void setValue(Float value) {
                                          for (int c = 0; c < 5; c++) {
                                              for (int i = 0; i < 6; i++) {
                                                  newPlus[(c * 6)
                                                          + i] = (float) (plus[i]
                                                                          + (i
                                                                             % 2 == 0 ? -Math.toRadians(value)
                                                                                      : Math.toRadians(value)));
                                              }
                                          }
                                          for (int i = 0; i < a.getNodes().length; i++) {
                                              a.getNodes()[i].setState(newPlus);
                                          }
                                      }
                                  }, 90f)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
