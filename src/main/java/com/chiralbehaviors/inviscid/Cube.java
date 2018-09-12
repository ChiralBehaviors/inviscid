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
public class Cube extends RegularPolyhedron<Quadrilateral> {

    protected Cube(int index, int[] vertices, int[] edges, int[] faces) {
        super(index, vertices, edges, faces);
    }

    @Override
    protected Quadrilateral[] constructFaces(int cube) {
        Quadrilateral[] squares = new Quadrilateral[6];
        for (int i = 0; i < 6; i++) {
            int offset = i * 6;
            int[] vertices = new int[4];
            squares[i] = new Quadrilateral(vertices,
                                           Arrays.copyOfRange(faces, offset,
                                                              offset + 6));
        }
        return squares;
    }

}
