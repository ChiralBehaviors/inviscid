/**
 * Copyright (c) 2018 Chiral Behaviors, LLC, all rights reserved.
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

import javax.vecmath.Vector3d;

import javafx.geometry.Point3D;
import javafx.scene.Group;
import javafx.scene.paint.Material;
import javafx.scene.shape.CullFace;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Translate;
import mesh.Face;
import mesh.polyhedra.plato.Octahedron;

/**
 * @author halhildebrand
 *
 */
public class Jitterbug {

    private final Group       group        = new Group();
    private final Rotate[]    rotations    = new Rotate[8];
    private final Translate[] translations = new Translate[8];
    private final double      Z;

    public Jitterbug(Octahedron oct, Material[] materials) {
        int i = 0;
        for (Face f : oct.getFaces()) {
            MeshView face = construct(f, materials[i], i++);
            group.getChildren()
                 .add(face);
        }
        Z = (oct.getEdgeLength() * Math.sqrt(2)) / Math.sqrt(3);
    }

    public Group getGroup() {
        return group;
    }

    public void rotateTo(double angle) {
        for (int i = 0; i < rotations.length; i++) {
            Rotate rotation = rotations[i];
            rotation.setAngle(angle);
            translate(i, angle);
        }
    }

    private MeshView construct(Face face, Material material, int index) {
        Vector3d centroid = face.centroid();
        Point3D axis = new Point3D(centroid.x, centroid.y, centroid.z);

        Rotate rotation = new Rotate();
        rotation.setPivotX(centroid.x);
        rotation.setPivotY(centroid.y);
        rotation.setPivotZ(centroid.z);
        rotation.setAxis(axis.normalize());

        Translate translation = new Translate();
        translation.setX(0d - axis.getX());
        translation.setY(0d - axis.getY());
        translation.setZ(0d - axis.getZ());

        int[] texIindices = new int[6];
        Arrays.fill(texIindices, 1);
        float[] texCoords = new float[6];
        float[] meshPoints = new float[9];
        Arrays.fill(texCoords, 1f);
        int i = 0;
        for (Vector3d v : face.getVertices()) {
            rotation.setAngle(-180);
            Point3D vertex = rotation.transform(new Point3D(v.x, v.y, v.z));
            vertex = translation.transform(vertex);

            meshPoints[i++] = (float) vertex.getX();
            meshPoints[i++] = (float) vertex.getY();
            meshPoints[i++] = (float) vertex.getZ();
        }

        TriangleMesh newMesh = new TriangleMesh();
        newMesh.getPoints()
               .addAll(meshPoints);
        newMesh.getTexCoords()
               .addAll(texCoords);
        newMesh.getFaces()
               .addAll(new int[] { 0, 1, 1, 1, 2, 1 });
        MeshView view = new MeshView(newMesh);
        rotations[index] = rotation;
        translations[index] = translation;
        view.getTransforms()
            .addAll(translation, rotation);
        view.setCullFace(CullFace.NONE);
        view.setMaterial(material);
        return view;
    }

    private void translate(int i, double angle) {
        Translate translation = translations[i];
        Point3D translate = rotations[i].getAxis()
                                        .multiply(Z
                                                  * Math.cos(Math.toRadians(angle)));
        translation.setX(translate.getX());
        translation.setY(translate.getY());
        translation.setZ(translate.getZ());
    }
}
