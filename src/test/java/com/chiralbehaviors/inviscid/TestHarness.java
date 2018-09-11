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

import com.chiralbehaviors.jfx.viewer3d.ContentModel;
import com.chiralbehaviors.jfx.viewer3d.Jfx3dViewerApp;

import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;

/**
 * @author halhildebrand
 *
 */
public class TestHarness extends Jfx3dViewerApp {
    public static void main(String[] args) {
        launch(args);
    }

    @Override
    protected ContentModel createContentModel() {
        ContentModel content = super.createContentModel();
        TriangleMesh mesh = PhiCoordinates.Tetrahedrons[0].constructRegularTexturedMesh();
        content.setContent(new MeshView(mesh));
        return content;
    }
}
