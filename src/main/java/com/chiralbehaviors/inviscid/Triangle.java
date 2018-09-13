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

import java.util.List;

import javafx.scene.Group;
import javafx.scene.paint.Material;
import javafx.scene.shape.Cylinder;
import javafx.scene.shape.Sphere;

/**
 * @author halhildebrand
 *
 */
public class Triangle extends Polygon {

    private static int[] edgesOf(int[] vertices) {
        return new int[] { vertices[0], vertices[1], vertices[1], vertices[2],
                           vertices[2], vertices[0] };
    }

    protected Triangle(int[] v) {
        super(v, edgesOf(v), v);
        assert vertices.length == 3 && edges.length == 3
               && triangles.length == 3;
    }

    public Group constructFace(double vertexRadius, Material[] vertexMaterial,
                               double edgeRadius, Material edgeMaterial[]) {
        Group group = new Group();
        int i = 0;
        for (Edge edge : getEdges()) {
            Cylinder l = edge.createLine(edgeRadius);
            l.setMaterial(edgeMaterial[i++]);
            group.getChildren()
                 .add(l);
        }
        List<Sphere> spheres = getVertices(vertexRadius);
        i = 0;
        for (Sphere s : spheres) {
            s.setMaterial(vertexMaterial[i++]);
        }
        group.getChildren()
             .addAll(spheres);
        return group;
    }
}
