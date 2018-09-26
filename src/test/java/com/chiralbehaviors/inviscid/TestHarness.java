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

import static com.chiralbehaviors.inviscid.animations.Colors.*;

import com.chiralbehaviors.inviscid.CubicGrid.Neighborhood;
import com.chiralbehaviors.inviscid.animations.PolyView;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.scene.Group;
import mesh.polyhedra.plato.Cube;

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
        Cube[] cubes = PhiCoordinates.cubes();
        group.getChildren()
             .add(new CubicGrid(Neighborhood.SIX, cubes[3], 4).construct(blackMaterial,
                                                              blackMaterial,
                                                              blackMaterial));
        group.getChildren()
             .add(new CubicGrid(Neighborhood.SIX, cubes[0], 2).construct(redMaterial,
                                                              redMaterial,
                                                              redMaterial));
        content.setContent(group);
    }
}
