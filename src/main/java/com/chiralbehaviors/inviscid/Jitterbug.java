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

import static com.chiralbehaviors.inviscid.animations.Colors.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Vector3d;

import javafx.geometry.Point3D;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.paint.Material;
import javafx.scene.shape.CullFace;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Translate;
import mesh.Face;
import mesh.Line;
import mesh.polyhedra.plato.Octahedron;

/**
 * @author halhildebrand
 *
 */
public class Jitterbug {
    private boolean[]         INVERSES     = new boolean[] { false, false,
                                                             false, false, true,
                                                             true, true, true };
    private final Group       group        = new Group();
    private final Rotate[]    rotations    = new Rotate[8];
    private final Translate[] translations = new Translate[8];
    private final double      Z;

    public Jitterbug(Octahedron oct, Material[] materials) {
        int i = 0;
        for (Face f : oct.getFaces()) {
            group.getChildren()
                 .add(construct(f, materials[i], i++));
        }
        Z = (oct.getEdgeLength() * Math.sqrt(2)) / Math.sqrt(3);
    }

    public Group getGroup() {
        return group;
    }

    public void rotateTo(double a) {
        for (int i = 0; i < 8; i++) {
            double angle = INVERSES[i] ? -a : a;
            Rotate rotation = rotations[i];
            rotation.setAngle(angle);
            translate(i, angle);
        }
    }

    private Node construct(Face face, Material material, int index) {
        Vector3d centroid = face.centroid();
        Point3D axis = new Point3D(centroid.x, centroid.y, centroid.z);

        Rotate rotation = new Rotate();
        rotation.setPivotX(centroid.x);
        rotation.setPivotY(centroid.y);
        rotation.setPivotZ(centroid.z);
        rotation.setAxis(axis.normalize());

        Translate translation = new Translate();
        translation.setX(-axis.getX());
        translation.setY(-axis.getY());
        translation.setZ(-axis.getZ());

        int[] texIindices = new int[6];
        Arrays.fill(texIindices, 1);
        float[] texCoords = new float[6];
        float[] meshPoints = new float[9];
        Arrays.fill(texCoords, 1f);
        int i = 0;
        List<Point3D> vertices = new ArrayList<>();
        for (Vector3d v : face.getVertices()) {
            if (INVERSES[index]) {
                rotation.setAngle(60);
            } else {
                rotation.setAngle(-60);
            }
            Point3D vertex = rotation.transform(new Point3D(v.x, v.y, v.z));
            vertex = translation.transform(vertex);

            meshPoints[i++] = (float) vertex.getX();
            meshPoints[i++] = (float) vertex.getY();
            meshPoints[i++] = (float) vertex.getZ();
            vertices.add(vertex);
        }
        vertices.add(vertices.get(0));
        Group faceGroup = new Group();
        rotation.setAngle(0);
        translation.setX(0);
        translation.setY(0);
        translation.setZ(0);
        rotations[index] = rotation;
        translations[index] = translation;
        faceGroup.getTransforms()
                 .addAll(rotation, translation);

        TriangleMesh newMesh = new TriangleMesh();
        newMesh.getPoints()
               .addAll(meshPoints);
        newMesh.getTexCoords()
               .addAll(texCoords);
        newMesh.getFaces()
               .addAll(new int[] { 0, 1, 1, 1, 2, 1 });
        MeshView view = new MeshView(newMesh);
        view.setCullFace(CullFace.NONE);
        view.setMaterial(material);
        faceGroup.getChildren()
                 .addAll(view);

        for (int j = 0; j < 3; j++) {
            Line line = new Line(0.01, vertices.get(j), vertices.get(j + 1));
            line.setMaterial(blackMaterial);
            faceGroup.getChildren()
                     .add(line);
        }

        return faceGroup;
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
