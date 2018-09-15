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

import static javafx.scene.paint.Color.*;
import static javafx.scene.paint.Color.BLUE;
import static javafx.scene.paint.Color.DEEPSKYBLUE;
import static javafx.scene.paint.Color.GREEN;
import static javafx.scene.paint.Color.ORANGE;
import static javafx.scene.paint.Color.RED;
import static javafx.scene.paint.Color.VIOLET;
import static javafx.scene.paint.Color.YELLOW;

import com.chiralbehaviors.jfx.viewer3d.ContentModel;
import com.chiralbehaviors.jfx.viewer3d.Jfx3dViewerApp;

import javafx.scene.Group;
import javafx.scene.paint.Color;
import javafx.scene.paint.Material;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.Sphere;
import mesh.Edge;
import mesh.Face;
import mesh.polyhedra.Polyhedron;
import mesh.polyhedra.plato.Dodecahedron;
import mesh.polyhedra.plato.Octahedron;

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
                                          orangeMaterial, violetMaterial,
                                          orangeMaterial, cyanMaterial,
                                          purpleMaterial, magentaMaterial,
                                          lavenderMaterial, oliveMaterial,
                                          limeMaterial };
        blackMaterials = new PhongMaterial[] { blackMaterial, blackMaterial };
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
            .forEach(e -> group.getChildren()
                               .addAll(e.createLine(.015)));
    }

    @Override
    protected ContentModel createContentModel() {

        ContentModel content = super.createContentModel();
        Dodecahedron poly = PhiCoordinates.dodecahedron();
        Group group = new Group();
        int i = 0;
        for (Sphere s : poly.getVertices(0.25)) {
            s.setMaterial(materials[i]);
            group.getChildren()
                 .add(s);
            i = (i + 1) % materials.length;
        }
        for (Face f : poly.getFaces()
                          .get(0)
                          .divideIntoTriangles()) {
            f.setMesh(poly);
            for (Edge e : f.getEdges()) {
                e.setMesh(poly);
                group.getChildren()
                     .add(e.createLine(.05));
            }
        }
//        addPolyhedron(group, poly, blueMaterial, materials);
        //        addPolyhedron(group, PhiCoordinates.cubes()[0], redMaterial, materials);
        //        addPolyhedron(group, PhiCoordinates.tetrahedrons()[8], redMaterial,
        //                      blackMaterials);
        //        addPolyhedron(group, PhiCoordinates.tetrahedrons()[9], blueMaterial,
        //                      blackMaterials);

        //        for (Cube cube : PhiCoordinates.cubes()) {
        //            addPolyhedron(group, cube, blueMaterial,
        //                          new PhongMaterial[] { redMaterial, redMaterial });
        //        }
        content.setContent(group);
        return content;
    }
}
