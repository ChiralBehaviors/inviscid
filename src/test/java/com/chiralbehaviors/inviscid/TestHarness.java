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
import static javafx.scene.paint.Color.RED;

import com.chiralbehaviors.jfx.viewer3d.ContentModel;
import com.chiralbehaviors.jfx.viewer3d.Jfx3dViewerApp;

import javafx.scene.Group;
import javafx.scene.paint.Color;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.Sphere;

/**
 * @author halhildebrand
 *
 */
public class TestHarness extends Jfx3dViewerApp {
    protected static final PhongMaterial   redMaterial;
    protected static final PhongMaterial   blueMaterial;
    protected static final PhongMaterial   greenMaterial;

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
        materials = new PhongMaterial[] { greenMaterial, redMaterial,
                                          blueMaterial };
    }

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    protected ContentModel createContentModel() {
        ContentModel content = super.createContentModel();
        Group group = new Group();
        for (int i = 0 ; i < 5; i ++) {
            Tetrahedron tetrahedron = PhiCoordinates.Tetrahedrons[i];
            for (Sphere s : tetrahedron.getVertices(0.5, blueMaterial)) {
                group.getChildren()
                     .add(s);
            }
            for (Edge e : tetrahedron.constructEdges()) {
                group.getChildren()
                     .add(e.createLine(0.1));
            }
        }
        content.setContent(group);
        return content;
    }
}
