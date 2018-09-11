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

import static com.chiralbehaviors.inviscid.PhiCoordinates.p120v;

import java.util.ArrayList;
import java.util.List;

import javafx.geometry.Point3D;
import javafx.scene.paint.Material;
import javafx.scene.shape.Sphere;

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

    public List<Sphere> getVertices(double radius, Material material) {
        List<Sphere> spheres = new ArrayList<>();
        for (int index : vertices) {
            Sphere sphere = new Sphere();
            sphere.setRadius(radius);
            Point3D vertex = p120v[index];
            sphere.setTranslateX(vertex.getX());
            sphere.setTranslateY(vertex.getY());
            sphere.setTranslateZ(vertex.getZ());
            sphere.setMaterial(material);
            spheres.add(sphere);
        }
        return spheres;
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
