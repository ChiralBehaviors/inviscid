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
import static com.chiralbehaviors.inviscid.animations.Colors.blueMaterial;
import static com.chiralbehaviors.inviscid.animations.Colors.greenMaterial;
import static com.chiralbehaviors.inviscid.animations.Colors.materials;
import static com.chiralbehaviors.inviscid.animations.Colors.redMaterial;

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

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group group = new Group();
        List<Jitterbug> jitterbugs = new ArrayList<>();
        Octahedron[] octahedrons = PhiCoordinates.Octahedrons;

        group.getChildren()
             .add(new CubicGrid(Neighborhood.EIGHT, PhiCoordinates.Cubes[3],
                                2).construct(blackMaterial, blackMaterial,
                                             blackMaterial));

        Jitterbug j = new Jitterbug(octahedrons[4], materials,
                                    new boolean[] { true, false, false, false,
                                                    false, false, false,
                                                    false });
        jitterbugs.add(j);
        j.rotateTo(0);
        group.getChildren()
             .add(j.getGroup());
        Octahedron oct = PhiCoordinates.Octahedrons[4];
        Ellipse ellipse = new Ellipse(0, oct, 0);
        group.getChildren()
             .addAll(ellipse.construct(40, redMaterial, 0.1));
        ellipse = new Ellipse(0, oct, 1);
        group.getChildren()
             .addAll(ellipse.construct(40, blueMaterial, 0.1));
        ellipse = new Ellipse(0, oct, 2);
        group.getChildren()
             .addAll(ellipse.construct(40, greenMaterial, 0.1));
        content.setContent(group);
        final Timeline timeline = new Timeline();
        AtomicReference<Double> angle = new AtomicReference<>(0d);
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(10_000),
                                  new KeyValue(new WritableValue<Double>() {
                                      @Override
                                      public Double getValue() {
                                          return angle.get();
                                      }

                                      @Override
                                      public void setValue(Double value) {
                                          angle.set(value);
                                          for (Jitterbug j : jitterbugs) {
                                              j.rotateTo(value);
                                          }
                                      }
                                  }, (double) 360)));
        timeline.setCycleCount(9000);
        content.setTimeline(timeline);
    }
}
