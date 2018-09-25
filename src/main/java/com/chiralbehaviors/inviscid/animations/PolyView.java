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

import com.javafx.experiments.jfx3dviewer.Jfx3dViewerApp;

import javafx.scene.Group;
import javafx.scene.paint.Material;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.Sphere;
import mesh.Edge;
import mesh.Face;
import mesh.Line;
import mesh.polyhedra.Polyhedron;

/**
 * @author halhildebrand
 *
 */
public class PolyView extends Jfx3dViewerApp {

    public final double __TAU        = (1.0 + Math.sqrt(5.0)) / 2.0;
    public final double DODECA_ANGLE = Math.acos((3.0 * (__TAU + (1.0 / 3.0)))
                                                 / ((__TAU * __TAU)
                                                    * (2.0 * Math.sqrt(2))));
    public final double ICOSA_ANGLE  = Math.acos((__TAU * __TAU)
                                                 / (2.0 * Math.sqrt(2)));

    public void addEdgesOf(Group group, Polyhedron poly, double radius,
                           Material material) {
        for (Edge edge : poly.getEdges()) {
            Line line = edge.createLine(radius);
            line.setMaterial(material);
            group.getChildren()
                 .add(line);
        }
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

    public void addVerticesOf(Group group, Polyhedron poly, double radius,
                              Material[] materials) {
        int i = 0;
        for (Sphere s : poly.getVertices(radius)) {
            s.setMaterial(materials[i++]);
            group.getChildren()
                 .add(s);
        }
    }

}
