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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javafx.geometry.Point3D;
import javafx.scene.paint.Material;
import javafx.scene.shape.Sphere;
import javafx.scene.shape.TriangleMesh;

/**
 * @author halhildebrand
 *
 */
abstract public class Polyhedron {
    protected final int[] edges;
    protected final int[] faces;
    protected final int[] vertices;

    protected Polyhedron(int[] vertices, int[] edges, int[] faces) {
        this.vertices = vertices;
        this.edges = edges;
        this.faces = faces;
    }

    public TriangleMesh constructMesh(int[] texIndices, float[] texCoords) {
        TriangleMesh newMesh = new TriangleMesh();
        newMesh.getPoints()
               .addAll(MeshPoints);
        newMesh.getTexCoords()
               .addAll(texCoords);
        int[] facesAndTexCoords = new int[faces.length * 2];
        int i = 0;
        for (int j = 0; j < faces.length; j++) {
            facesAndTexCoords[i] = faces[j];
            facesAndTexCoords[i++] = texIndices[j];
        }
        newMesh.getFaces()
               .addAll(facesAndTexCoords);
        return newMesh;
    }

    public TriangleMesh constructRegularTexturedMesh() {
        int[] indices = new int[faces.length];
        Arrays.fill(indices, 1);
        float[] texCoords = new float[faces.length];
        Arrays.fill(texCoords, 0.5f);
        return constructMesh(indices, texCoords);
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
}
