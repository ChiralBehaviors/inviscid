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

import static javafx.scene.paint.Color.BLACK;
import static javafx.scene.paint.Color.BLUE;
import static javafx.scene.paint.Color.CYAN;
import static javafx.scene.paint.Color.GREEN;
import static javafx.scene.paint.Color.LAVENDER;
import static javafx.scene.paint.Color.LIME;
import static javafx.scene.paint.Color.MAGENTA;
import static javafx.scene.paint.Color.OLIVE;
import static javafx.scene.paint.Color.ORANGE;
import static javafx.scene.paint.Color.PURPLE;
import static javafx.scene.paint.Color.RED;
import static javafx.scene.paint.Color.VIOLET;
import static javafx.scene.paint.Color.YELLOW;

import javax.vecmath.Vector3d;

import com.javafx.experiments.jfx3dviewer.ContentModel;
import com.javafx.experiments.jfx3dviewer.Jfx3dViewerApp;

import javafx.geometry.Point3D;
import javafx.scene.Group;
import javafx.scene.paint.Material;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.Cylinder;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.Sphere;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Translate;
import mesh.Face;
import mesh.polyhedra.Polyhedron;
import mesh.polyhedra.plato.Octahedron;
import mesh.polyhedra.plato.Tetrahedron;

/**
 * @author halhildebrand
 *
 */
public class TestHarness extends Jfx3dViewerApp {
    protected static final PhongMaterial   blueMaterial;
    protected static final PhongMaterial   greenMaterial;
    protected static final PhongMaterial[] materials;
    protected static final PhongMaterial[] blackMaterials;
    protected static final PhongMaterial   redMaterial;
    protected static final PhongMaterial   yellowMaterial;
    protected static final PhongMaterial   blackMaterial;
    protected static final PhongMaterial   violetMaterial;
    protected static final PhongMaterial   orangeMaterial;
    protected static final PhongMaterial   cyanMaterial;
    protected static final PhongMaterial   purpleMaterial;
    protected static final PhongMaterial   magentaMaterial;
    protected static final PhongMaterial   lavenderMaterial;
    protected static final PhongMaterial   oliveMaterial;
    protected static final PhongMaterial   limeMaterial;
    protected static final PhongMaterial[] eightMaterials;
    protected static final PhongMaterial[] eight4Materials;

    static {
        redMaterial = new PhongMaterial(RED);
        blueMaterial = new PhongMaterial(BLUE);
        greenMaterial = new PhongMaterial(GREEN);
        yellowMaterial = new PhongMaterial(YELLOW);
        violetMaterial = new PhongMaterial(VIOLET);
        orangeMaterial = new PhongMaterial(ORANGE);
        cyanMaterial = new PhongMaterial(CYAN);
        purpleMaterial = new PhongMaterial(PURPLE);
        magentaMaterial = new PhongMaterial(MAGENTA);
        lavenderMaterial = new PhongMaterial(LAVENDER);
        oliveMaterial = new PhongMaterial(OLIVE);
        limeMaterial = new PhongMaterial(LIME);

        blackMaterial = new PhongMaterial(BLACK);

        materials = new PhongMaterial[] { redMaterial, blueMaterial,
                                          greenMaterial, yellowMaterial,
                                          cyanMaterial, purpleMaterial,
                                          orangeMaterial, violetMaterial,
                                          magentaMaterial, lavenderMaterial,
                                          oliveMaterial, limeMaterial };
        blackMaterials = new PhongMaterial[] { blackMaterial, blackMaterial };
        eightMaterials = new PhongMaterial[] { redMaterial, blueMaterial,
                                               greenMaterial, yellowMaterial,
                                               redMaterial, blueMaterial,
                                               greenMaterial, yellowMaterial };
        eight4Materials = new PhongMaterial[] { redMaterial, blueMaterial,
                                                greenMaterial, yellowMaterial,
                                                blackMaterial, blackMaterial,
                                                blackMaterial, blackMaterial };
    }

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
        for (Sphere s : poly.getVertices(0.15)) {
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
        Group group = new Group();

        Octahedron oct = PhiCoordinates.octahedrons()[0];
        Face f = oct.getFaces()
                       .get(0);
        Vector3d centroid = f.centroid();

        Tetrahedron tet1 = PhiCoordinates.tetrahedrons()[0];
        
        Vector3d axis = tet1.getVertexPositions().get(0);
        
        Rotate rotation = new Rotate();
        rotation.setAngle(30);
        
        rotation.setPivotX(centroid.x);
        rotation.setPivotY(centroid.y);
        rotation.setPivotZ(centroid.z);
        
        rotation.setAxis(new Point3D(axis.x, axis.y, axis.z));
        
        Vector3d direction = new Vector3d(); 
        direction.normalize(axis);
        direction.scale(1.5);
        
        Translate translate = new Translate();
        translate.setX(direction.x);
        translate.setY(direction.y);
        translate.setZ(direction.z);
       
        
        MeshView face = f.constructMeshView();
        face.setMaterial(redMaterial);
        group.getChildren().add(face);
        
        MeshView rotated = f.constructMeshView(); 
        rotated.setMaterial(blueMaterial);
        rotated.getTransforms().addAll(rotation, translate);
        
        group.getChildren().add(rotated);
        content.setContent(group);
    }
}
