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
public class Cube extends RegularPolyhedron<Square> {

    protected Cube(int index, int[] vertices, int[] edges, int[] faces) {
        super(index, vertices, edges, faces);
    }

    @Override
    protected Square[] constructFaces(int cube) {
        Square[] faces = new Square[6];
        for (int i = 0; i < 4; i++) {
            faces[i] = newFace(cube, i);
        }
        return faces;
    }

    private static Square newFace(int cube, int index) {
        int[] vertices = new int[4];
        for (int i = 0; i < 4; i++) {
            vertices[i] = CubeEdgeMap[cube][i + (index * 4)];
        }

        int[] triangles = new int[6];
        for (int i = 0; i < 6; i++) {
            triangles[i] = CubeFaceMap[index][i + (index * 6)];
        }

        return new Square(vertices, triangles);
    }

}
