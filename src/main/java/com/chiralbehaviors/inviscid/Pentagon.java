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

/**
 * @author halhildebrand
 *
 */
public class Pentagon extends Polygon {

    protected Pentagon(int[] v, int[] triangles) {
        super(v, edgesOf(v), triangles);
        assert vertices.length == 5 && edges.length == 5
               && triangles.length == 15;
    }

    private static int[] edgesOf(int[] vertices) {
        return new int[] { vertices[0], vertices[1], vertices[1], vertices[2],
                           vertices[2], vertices[3], vertices[3], vertices[4],
                           vertices[4], vertices[0] };
    }

}
