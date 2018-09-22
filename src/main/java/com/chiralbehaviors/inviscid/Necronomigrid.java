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

package com.chiralbehaviors.inviscid;

import static com.chiralbehaviors.inviscid.Colors.materials;

import java.util.List;

import javax.vecmath.Vector3d;

import javafx.geometry.Point3D;
import javafx.scene.Group;
import mesh.Face;
import mesh.Line;
import mesh.polyhedra.plato.Cube;

/**
 * @author halhildebrand
 *
 */
public class Necronomigrid {

    public Group construct() {
        Group group = new Group();
        constructGrid(PhiCoordinates.cubes()[0], group);
        return group;
    }

    private void constructGrid(Cube cube, Group group) {
        List<Face> faces = cube.getFaces();
        for (int i = 0; i < 3; i++) {
            Face a = faces.get(i);
            Face b = faces.get(i + 3);
            Vector3d Ac = a.centroid();
            Vector3d Bc = b.centroid();
            Point3D positive = new Point3D(Ac.x, Ac.y, Ac.z);
            Point3D negative = new Point3D(Bc.x, Bc.y, Bc.z);
            Line axisP = new Line(0.1, positive, new Point3D(0, 0, 0));
            axisP.setMaterial(materials[i]);
            Line axisN = new Line(0.1, negative, new Point3D(0, 0, 0));
            axisN.setMaterial(materials[i+3]);
            group.getChildren()
                 .addAll(axisP, axisN);
        }
    }

}
