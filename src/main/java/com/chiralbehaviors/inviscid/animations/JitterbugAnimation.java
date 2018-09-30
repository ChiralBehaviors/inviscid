/**
 * Copyright (C) 2018 Chiral Behaviors, LLC. All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package com.chiralbehaviors.inviscid.animations;

import static com.chiralbehaviors.inviscid.animations.Colors.blackMaterial;
import static com.chiralbehaviors.inviscid.animations.Colors.materials;
import static com.chiralbehaviors.inviscid.animations.Colors.*;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;

import com.chiralbehaviors.inviscid.CubicGrid;
import com.chiralbehaviors.inviscid.CubicGrid.Neighborhood;
import com.chiralbehaviors.inviscid.Jitterbug;
import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.KeyFrame;
import javafx.animation.KeyValue;
import javafx.animation.Timeline;
import javafx.beans.value.WritableValue;
import javafx.scene.Group;
import javafx.util.Duration;
import mesh.Ellipse;
import mesh.polyhedra.plato.Octahedron;

/**
 * @author halhildebrand
 *
 */
public class JitterbugAnimation extends PolyView {
    public static void main(String[] args) {
        launch(args);
    }

    public void jitterbugArray(Group group, Octahedron[] octahedrons,
                               List<Jitterbug> jitterbugs, CubicGrid grid,
                               double initialAngle) {
        for (int x : new int[] { -1, 1 }) {
            for (int y : new int[] { -1, 1 }) {
                for (int z : new int[] { -1, 1 }) {
                    Jitterbug j = new Jitterbug(octahedrons[4], materials);
                    j.rotateTo(initialAngle);
                    Group jGroup = j.getGroup();
                    grid.postition(x, y, z, jGroup);
                    group.getChildren()
                         .add(jGroup);
                    jitterbugs.add(j);
                }
            }
        }
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group group = new Group();
        List<Jitterbug> jitterbugs = new ArrayList<>();

        CubicGrid grid = new CubicGrid(Neighborhood.EIGHT,
                                       PhiCoordinates.Cubes[3], 1);
        group.getChildren()
             .add(grid.construct(blackMaterial, blackMaterial, blackMaterial));

        Jitterbug j = new Jitterbug(PhiCoordinates.Octahedrons[4], materials,
                                    new boolean[] { true, false, false, false,
                                                    false, false, false,
                                                    false });
        j.rotateTo(00);
        Group jGroup = j.getGroup();
        group.getChildren()
             .add(jGroup);
        jitterbugs.add(j);

        group.getChildren()
             .add(new Ellipse(0, PhiCoordinates.Octahedrons[4],
                              0).construct(40, redMaterial, 0.015));
        group.getChildren()
             .add(new Ellipse(0, PhiCoordinates.Octahedrons[4],
                              1).construct(40, blueMaterial, 0.015));
        group.getChildren()
             .add(new Ellipse(0, PhiCoordinates.Octahedrons[4],
                              2).construct(40, greenMaterial, 0.015));

        //        jitterbugArray(group, PhiCoordinates.Octahedrons, jitterbugs, grid, (double) 0);

        content.setContent(group);
        final Timeline timeline = new Timeline();
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(5_000),
                                  new KeyValue(new WritableValue<Double>() {
                                      @Override
                                      public Double getValue() {
                                          return 0d;
                                      }

                                      @Override
                                      public void setValue(Double value) {
                                          jitterbugs.forEach(j -> j.rotateTo(value));
                                      }
                                  }, 60d)));
                timeline.setCycleCount(9000);
                timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
