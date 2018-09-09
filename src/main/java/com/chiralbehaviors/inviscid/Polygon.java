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
abstract public class Polygon {
    protected final Edge[] edges;
    protected final int[]  triangles;
    protected final int[]  vertices;

    protected Polygon(int[] vertices, int[] edges, int[] triangles) {
        this.vertices = vertices;
        this.triangles = triangles;
        this.edges = constructEdges(edges);
    }

    public Edge[] getEdges() {
        return edges;
    }

    private Edge[] constructEdges(int[] edges) {
        Edge[] e = new Edge[edges.length / 2];
        int j = 0;
        for (int i = 0; i < edges.length; i++) {
            e[j] = new Edge(new int[] { edges[1], edges[i++] });
            j++;
        }
        return e;
    }
}
