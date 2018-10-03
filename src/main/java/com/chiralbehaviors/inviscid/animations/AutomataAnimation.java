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

import java.util.ArrayList;
import java.util.List;
import static com.chiralbehaviors.inviscid.animations.Colors.*;
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
import javafx.scene.paint.Material;
import javafx.scene.paint.PhongMaterial;
import javafx.util.Duration;
import javafx.util.Pair;

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
        ContentModel content = getContentModel();
        Group group = new Group();
        content.setContent(group);
        CubicGrid grid = new CubicGrid(Neighborhood.SIX,
                                       PhiCoordinates.Cubes[3], 1);
        group.getChildren()
             .add(grid.construct(blackMaterial, blackMaterial, blackMaterial));
        double[] plus = Automata.getPositiveTet();
        double[] minus = Automata.getNegativeTet();
        Pair<List<CellNode>, List<CellNode>> cells = cellArray(group, grid,
                                                               plus, minus);

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
                                          double[] newPlus = new double[plus.length];
                                          double[] newMinus = new double[minus.length];
                                          for (int i = 0; i < plus.length; i++) {
                                              newPlus[i] = plus[i]
                                                           + (i
                                                              % 2 == 0 ? -Math.toRadians(value)
                                                                       : Math.toRadians(value));
                                              newMinus[i] = minus[i]
                                                            + (i
                                                               % 2 == 0 ? -Math.toRadians(value)
                                                                        : Math.toRadians(value));
                                          }
                                          for (int i = 0; i < cells.getKey()
                                                                   .size(); i++) {
                                              cells.getKey()
                                                   .get(i)
                                                   .setState(newPlus);
                                          }
                                          for (int i = 0; i < cells.getValue()
                                                                   .size(); i++) {
                                              cells.getValue()
                                                   .get(i)
                                                   .setState(newMinus);
                                          }
                                      }
                                  }, 90D)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }

    private Pair<List<CellNode>, List<CellNode>> cellArray(Group group,
                                                           CubicGrid grid,
                                                           double[] plus,
                                                           double[] minus) {
        Automata a = new Automata(360, grid, 0.1);
        List<CellNode> positive = new ArrayList<>();
        List<CellNode> negative = new ArrayList<>();
        CellNode cellGroup;

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(-1, -1, 0, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(-1, 0, -1, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(-1, 0, 1, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(-1, 1, 0, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(0, 0, 0, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(0, -1, -1, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(0, -1, 1, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(0, 1, -1, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(0, 1, 1, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(1, 0, -1, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(1, 0, 1, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(1, -1, 0, cellGroup);
        group.getChildren()
             .add(cellGroup);

        cellGroup = a.cellNode(plus, edgeMaterials);
        positive.add(cellGroup);
        grid.postition(1, 1, 0, cellGroup);
        group.getChildren()
             .add(cellGroup);
        return new Pair<>(positive, negative);
    }
}
