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

import com.chiralbehaviors.inviscid.Automata;
import com.chiralbehaviors.inviscid.CellNode;
import com.chiralbehaviors.inviscid.CubicGrid;
import com.chiralbehaviors.inviscid.CubicGrid.Neighborhood;
import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.KeyFrame;
import javafx.animation.KeyValue;
import javafx.animation.Timeline;
import javafx.beans.value.WritableValue;
import javafx.scene.Group;
import javafx.util.Duration;

/**
 * @author halhildebrand
 *
 */
public class AutomataAnimation extends PolyView {
    public static void main(String[] args) {
        launch(args);
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group group = new Group();
        CubicGrid grid = new CubicGrid(Neighborhood.SIX,
                                       PhiCoordinates.Cubes[3], 0);
        group.getChildren()
             .add(grid.construct(blackMaterial, blackMaterial, blackMaterial));
        Automata a = new Automata(360, grid, 0.1);
        double[] initialState = Automata.getPositiveTet();
        CellNode cellGroup = a.cellNode(initialState);
        grid.postition(0, 0, 0, cellGroup);
        group.getChildren()
             .add(cellGroup);
        content.setContent(group);

        final Timeline timeline = new Timeline();
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(2_500),
                                  new KeyValue(new WritableValue<Double>() {
                                      @Override
                                      public Double getValue() {
                                          return 0D;
                                      }

                                      @Override
                                      public void setValue(Double value) {
                                          double[] newState = new double[initialState.length];
                                          for (int i = 0; i < initialState.length; i++) {
                                              newState[i] = initialState[i]
                                                            + (i
                                                               % 2 == 0 ? -Math.toRadians(value)
                                                                        : Math.toRadians(value));
                                          }
                                          cellGroup.setState(newState);
                                      }
                                  }, 90D)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
