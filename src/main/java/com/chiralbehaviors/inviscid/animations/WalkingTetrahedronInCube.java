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
import static com.chiralbehaviors.inviscid.animations.Colors.materials;

import com.chiralbehaviors.inviscid.PhiCoordinates;

import javafx.scene.Group;
import mesh.polyhedra.plato.Cube;

/**
 * @author halhildebrand
 *
 */
public class WalkingTetrahedronInCube extends PolyView {
    public static void main(String[] args) {
        launch(args);
    }

    @Override
    protected void initializeContentModel() {
        Group group = new Group();
        Cube cube = PhiCoordinates.cubes()[0];
        addEdgesOf(group, cube, 0.015, blackMaterial);
        addVerticesOf(group, cube, 0.5, materials);
        addPlainPolyhedron(group, PhiCoordinates.tetrahedrons()[0],
                           blackMaterial);
        getContentModel().setContent(group);
    }

}
