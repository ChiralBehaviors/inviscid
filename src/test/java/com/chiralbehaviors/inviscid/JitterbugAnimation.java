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

package com.chiralbehaviors.inviscid;

import static com.chiralbehaviors.inviscid.Colors.blackMaterial;
import static com.chiralbehaviors.inviscid.Colors.materials;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;

import com.javafx.experiments.jfx3dviewer.ContentModel;
import com.javafx.experiments.jfx3dviewer.Jfx3dViewerApp;

import javafx.animation.KeyFrame;
import javafx.animation.KeyValue;
import javafx.animation.Timeline;
import javafx.beans.value.WritableValue;
import javafx.scene.Group;
import javafx.scene.paint.Material;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.Sphere;
import javafx.util.Duration;
import mesh.Face;
import mesh.Line;
import mesh.polyhedra.Polyhedron;
import mesh.polyhedra.plato.Octahedron;

/**
 * @author halhildebrand
 *
 */
public class JitterbugAnimation extends Jfx3dViewerApp {
    public final double TAU          = (1.0 + Math.sqrt(5.0)) / 2.0;

    public final double ICOSA_ANGLE  = Math.acos((TAU * TAU)
                                                 / (2.0 * Math.sqrt(2)));
    public final double DODECA_ANGLE = Math.acos((3.0 * (TAU + (1.0 / 3.0)))
                                                 / ((TAU * TAU)
                                                    * (2.0 * Math.sqrt(2))));

    public static void main(String[] args) {
        launch(args);
    }

    public void addPolyhedron(Group group, Polyhedron poly,
                              Material meshMaterial,
                              Material[] vertexMaterials) {
        MeshView view = new MeshView(poly.toTriangleMesh()
                                         .constructRegularTexturedMesh());
        view.setMaterial(meshMaterial);
        group.getChildren()
             .add(view);
        int i = 0;
        for (Sphere s : poly.getVertices(3.5)) {
            s.setMaterial(vertexMaterials[i]);
            i = (i + 1) % vertexMaterials.length;
            group.getChildren()
                 .add(s);
        }
        poly.getEdges()
            .forEach(e -> {
                Line line = e.createLine(.015);
                line.setMaterial(blackMaterial);
                group.getChildren()
                     .addAll(line);
            });
    }

    public void addPlainPolyhedron(Group group, Polyhedron poly,
                                   Material meshMaterial) {
        MeshView view = new MeshView(poly.toTriangleMesh()
                                         .constructRegularTexturedMesh());
        view.setMaterial(meshMaterial);
        group.getChildren()
             .add(view);
        poly.getEdges()
            .forEach(e -> {
                Line line = e.createLine(.015);
                line.setMaterial(blackMaterial);
                group.getChildren()
                     .addAll(line);
            });
    }

    public void addFaces(Group group, Polyhedron poly,
                         Material[] faceMaterials) {
        int i = 0;
        for (Face face : poly.getFaces()) {
            MeshView mesh = face.constructMeshView();
            mesh.setMaterial(faceMaterials[i++]);
            group.getChildren()
                 .add(mesh);
        }
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group group = new Group();
        List<Jitterbug> jitterbugs = new ArrayList<>();
        Octahedron[] octahedrons = PhiCoordinates.octahedrons();
        Arrays.asList(0, 1, 2, 3, 4)
              .forEach(i -> {
                  Octahedron oct = octahedrons[i];
                  Jitterbug j = new Jitterbug(oct, materials);
                  jitterbugs.add(j);
                  j.rotateTo(90);
                  group.getChildren()
                       .add(j.getGroup());
              });
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
