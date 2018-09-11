
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

import java.util.Arrays;

/**
 * @author halhildebrand
 *
 */
public class Tetrahedron extends RegularPolyhedron<Triangle> {

    protected Tetrahedron(int index, int[] vertices, int[] edges, int[] faces) {
        super(index, vertices, edges, faces);
    }

    @Override
    protected Triangle[] constructFaces(int tet) {
        Triangle[] triangles = new Triangle[4];
        for (int i = 0; i < 4; i++) {
            int offset = i * 3;
            triangles[i] = new Triangle(Arrays.copyOfRange(faces, offset,
                                                           offset + 3));
        }
        return triangles;
    }
}
