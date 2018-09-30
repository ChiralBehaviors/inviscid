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

import static com.chiralbehaviors.inviscid.Constants.QUARTER_PI;
import static com.chiralbehaviors.inviscid.Constants.THREE_QUARTERS_PI;
import static com.chiralbehaviors.inviscid.Constants.TWO_PI;
import static com.chiralbehaviors.inviscid.animations.Colors.redMaterial;

import java.util.Arrays;

import javafx.geometry.Point3D;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Transform;
import javafx.scene.transform.Translate;
import mesh.Line;

/**
 * @author halhildebrand
 *
 */
public class Automata {

    private static final double[] POSITIVE_TET = new double[] { THREE_QUARTERS_PI,
                                                                QUARTER_PI,
                                                                QUARTER_PI,
                                                                THREE_QUARTERS_PI,
                                                                THREE_QUARTERS_PI,
                                                                QUARTER_PI };

    public static double[] getPositiveTet() {
        return Arrays.copyOf(POSITIVE_TET, POSITIVE_TET.length);
    }

    private final double         angularResolution;
    private final TriangleMesh[] exemplars;
    private final CubicGrid      grid;
    private LengthTable          lengths;

    public Automata(int resolution, CubicGrid grid, double radius) {
        assert resolution
               % 8 == 0 : "Angular resolution must be divisable by 8: "
                          + resolution;
        angularResolution = TWO_PI / resolution;
        this.grid = grid;
        exemplars = new TriangleMesh[resolution / 8];
        lengths = new LengthTable(resolution,
                                  grid.getIntervalX() * Math.sqrt(2) / 2.0);
        for (int i = 0; i < exemplars.length; i++) {
            exemplars[i] = Line.createLine(lengths.length(i * angularResolution)
                                           * 2, radius);
        }
    }

    public CellNode cellNode(double[] initialState) {
        double halfInterval = grid.getIntervalX() / 2.0;
        Translate xPos = new Translate(0, 0, halfInterval);
        Translate xNeg = new Translate(0, 0, -halfInterval);
        Translate yPos = new Translate(halfInterval, 0, 0);
        Translate yNeg = new Translate(-halfInterval, 0, -0);
        MeshView[][] struts = new MeshView[6][];
        for (int i = 0; i < 6; i++) {
            struts[i] = new MeshView[exemplars.length * 8];
        }
        for (int i = 0; i < exemplars.length * 8; i++) {
            MeshView view;
            view = xExemplar(i, xPos);
            view.setMaterial(redMaterial);
            struts[0][i] = view;

            view = xExemplar(i, xNeg);
            view.setMaterial(redMaterial);
            struts[1][i] = view;

            view = yExemplar(i, yPos);
            view.setMaterial(redMaterial);
            struts[2][i] = view;

            view = yExemplar(i, yNeg);
            view.setMaterial(redMaterial);
            struts[3][i] = view;

            view = zExemplar(i, yPos);
            view.setMaterial(redMaterial);
            struts[4][i] = view;

            view = zExemplar(i, yNeg);
            view.setMaterial(redMaterial);
            struts[5][i] = view;
        }
        return new CellNode(angularResolution, struts, initialState);
    }

    private TriangleMesh x(int increment, TriangleMesh[] exemplars) {
        if (increment < exemplars.length) {
            return exemplars[increment];
        }
        if (increment < exemplars.length * 2) {
            return exemplars[exemplars.length - (increment - exemplars.length)
                             - 1];
        }
        if (increment < exemplars.length * 3) {
            return exemplars[increment - (exemplars.length * 2)];
        }
        if (increment < exemplars.length * 4) {
            return exemplars[exemplars.length
                             - (increment - exemplars.length * 3) - 1];
        }
        if (increment < exemplars.length * 5) {
            return exemplars[increment - (exemplars.length * 4)];
        }
        if (increment < exemplars.length * 6) {
            return exemplars[exemplars.length
                             - (increment - exemplars.length * 5) - 1];
        }
        if (increment < exemplars.length * 7) {
            return exemplars[increment - (exemplars.length * 6)];
        }
        return exemplars[exemplars.length - (increment - exemplars.length * 7)
                         - 1];
    }

    private MeshView xExemplar(int increment, Transform translate) {
        Point3D axisOfRotation = new Point3D(0, 0, 1);
        Rotate rotateAroundCenter = new Rotate(-Math.toDegrees(increment
                                                               * angularResolution),
                                               axisOfRotation);
        MeshView line = new MeshView(x(increment, exemplars));

        line.getTransforms()
            .add(rotateAroundCenter.createConcatenation(translate));
        return line;
    }

    private MeshView yExemplar(int increment, Translate translate) {
        Point3D axisOfRotation = new Point3D(1, 0, 0);
        Rotate rotateAroundCenter = new Rotate(-Math.toDegrees(increment
                                                               * angularResolution),
                                               axisOfRotation);
        MeshView line = new MeshView(x(increment, exemplars));

        line.getTransforms()
            .add(rotateAroundCenter.createConcatenation(translate));
        return line;
    }

    private MeshView zExemplar(int increment, Translate translate) {
        Point3D axisOfRotation = new Point3D(1, 0, 0);
        Transform rotateAroundCenter = new Rotate(-Math.toDegrees(increment
                                                                  * angularResolution),
                                                  axisOfRotation);
        rotateAroundCenter = new Rotate(90, new Point3D(0, 0,
                                                        1)).createConcatenation(rotateAroundCenter);
        MeshView line = new MeshView(x(increment, exemplars));

        line.getTransforms()
            .add(rotateAroundCenter.createConcatenation(translate));
        return line;
    }
}
