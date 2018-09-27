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

import static com.chiralbehaviors.inviscid.animations.Colors.blackMaterial;

import com.chiralbehaviors.inviscid.animations.PolyView;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.scene.Group;
import mesh.Ellipse;
import mesh.PolyLine;

/**
 * @author halhildebrand
 *
 */
public class TestHarness extends PolyView {

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();

        Group group = new Group();
        for (int i = 0; i < 3; i++) {
            PolyLine line = new Ellipse(0, PhiCoordinates.Octahedrons[0],
                                        i).construct(20, blackMaterial, 0.1);
            group.getChildren()
                 .add(line);
        }
        group.getChildren()
             .add(PhiCoordinates.Octahedrons[0].getFaces()
                                               .get(0)
                                               .constructMeshView());
        content.setContent(group);
    }
}
