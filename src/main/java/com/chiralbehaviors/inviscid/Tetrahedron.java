
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

import static com.chiralbehaviors.inviscid.PhiCoordinates.*;

/**
 * @author halhildebrand
 *
 */
public class Tetrahedron extends RegularPolyhedron<Triangle> {

    protected Tetrahedron(int index, int[] vertices, int[] edges, int[] faces) {
        super(index, vertices, edges, faces);
    }

    @Override
    protected Triangle[] constructFaces(int index) {
        Triangle[] faces = new Triangle[4];
        for (int i = 0; i < 4; i++) {
            faces[i] = newFace(index, i);
        }
        return faces;
    }

    private static Triangle newFace(int tet, int index) {
        int[] vertices = new int[3];
        for (int i = index * 3; i < (index + 1) * 3; i++) {
            vertices[i] = TetrahedronVertices[tet][i];
        }
        int[] edges = new int[6];
        for (int i = 0; i < 3; i++) {
            edges[i] = new int[] {};
        }
        return new Triangle(vertices, edges, TetrahedronFaceMap[index]);
    }
}
