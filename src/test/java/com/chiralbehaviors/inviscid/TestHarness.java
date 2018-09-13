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

import static javafx.scene.paint.Color.BLUE;
import static javafx.scene.paint.Color.DARKBLUE;
import static javafx.scene.paint.Color.DARKGREEN;
import static javafx.scene.paint.Color.DARKRED;
import static javafx.scene.paint.Color.GREEN;
import static javafx.scene.paint.Color.RED;
import static javafx.scene.paint.Color.YELLOW;

import com.chiralbehaviors.jfx.viewer3d.ContentModel;
import com.chiralbehaviors.jfx.viewer3d.Jfx3dViewerApp;

import javafx.scene.Group;
import javafx.scene.paint.Color;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.CullFace;
import javafx.scene.shape.MeshView;

/**
 * @author halhildebrand
 *
 */
public class TestHarness extends Jfx3dViewerApp {
    protected static final PhongMaterial   redMaterial;
    protected static final PhongMaterial   blueMaterial;
    protected static final PhongMaterial   greenMaterial;
    protected static final PhongMaterial   yellowMaterial;

    protected static final PhongMaterial[] materials;

    static {
        redMaterial = new PhongMaterial();
        redMaterial.setDiffuseColor(new Color(DARKRED.getRed(),
                                              DARKRED.getGreen(),
                                              DARKRED.getBlue(), 0.6));
        redMaterial.setSpecularColor(new Color(RED.getRed(), RED.getGreen(),
                                               RED.getBlue(), 0.6));

        blueMaterial = new PhongMaterial();
        blueMaterial.setDiffuseColor(new Color(DARKBLUE.getRed(),
                                               DARKBLUE.getGreen(),
                                               DARKBLUE.getBlue(), 0.6));
        blueMaterial.setSpecularColor(new Color(BLUE.getRed(), BLUE.getGreen(),
                                                BLUE.getBlue(), 0.6));
        greenMaterial = new PhongMaterial();
        greenMaterial.setDiffuseColor(new Color(DARKGREEN.getRed(),
                                                DARKGREEN.getGreen(),
                                                DARKGREEN.getBlue(), 0.6));
        greenMaterial.setSpecularColor(new Color(GREEN.getRed(),
                                                 GREEN.getGreen(),
                                                 GREEN.getBlue(), 0.6));
        yellowMaterial = new PhongMaterial();
        yellowMaterial.setDiffuseColor(new Color(YELLOW.getRed(),
                                                 YELLOW.getGreen(),
                                                 YELLOW.getBlue(), 0.6));
        yellowMaterial.setSpecularColor(new Color(YELLOW.getRed(),
                                                  YELLOW.getGreen(),
                                                  YELLOW.getBlue(), 0.6));

        materials = new PhongMaterial[] { redMaterial, blueMaterial,
                                          greenMaterial, yellowMaterial };
    }

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    protected ContentModel createContentModel() {
        ContentModel content = super.createContentModel();
        Group group = new Group();
        addMesh(group, PhiCoordinates.Tetrahedrons[8]);
        addMesh(group, PhiCoordinates.Tetrahedrons[9]);
        content.setContent(group);
        return content;
    }

    private void addMesh(Group group, Tetrahedron tetrahedron) {
        MeshView view = new MeshView(tetrahedron.constructRegularTexturedMesh());
        view.setMaterial(redMaterial);
        view.setCullFace(CullFace.NONE);
        group.getChildren()
             .add(view);
    }
}
