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

import java.util.concurrent.atomic.AtomicReference;

import com.javafx.experiments.jfx3dviewer.ContentModel;
import com.javafx.experiments.jfx3dviewer.Jfx3dViewerApp;

import javafx.animation.KeyFrame;
import javafx.animation.KeyValue;
import javafx.animation.Timeline;
import javafx.beans.value.WritableValue;
import javafx.scene.Group;
import javafx.scene.paint.Material;
import javafx.scene.shape.Cylinder;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.Sphere;
import javafx.util.Duration;
import mesh.Face;
import mesh.polyhedra.Polyhedron;
import mesh.polyhedra.plato.Octahedron;

/**
 * @author halhildebrand
 *
 */
public class JitterbugAnimation extends Jfx3dViewerApp {

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
                Cylinder line = e.createLine(.015);
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
                Cylinder line = e.createLine(.015);
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
        Octahedron oct = PhiCoordinates.octahedrons()[0];
        Jitterbug j = new Jitterbug(oct, materials); 
        j.rotateTo(90);
        content.setContent(j.getGroup());
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
                                                          j.rotateTo(value);
                                                      }
                                                  }, (double) 360))); 
        timeline.setCycleCount(9000);
        content.setTimeline(timeline);
    }
}
