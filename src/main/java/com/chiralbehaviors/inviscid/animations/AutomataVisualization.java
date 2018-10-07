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

package com.chiralbehaviors.inviscid.animations;

import static com.chiralbehaviors.inviscid.Constants.QUARTER_PI;
import static com.chiralbehaviors.inviscid.Constants.ROOT_2_DIV_2;
import static com.chiralbehaviors.inviscid.Constants.THREE_QUARTERS_PI;
import static com.chiralbehaviors.inviscid.Constants.TWO_PI;
import static com.chiralbehaviors.inviscid.CubicGrid.yAxis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.CubicGrid;
import com.chiralbehaviors.inviscid.CubicGrid.Neighborhood;
import com.chiralbehaviors.inviscid.LengthTable;
import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.PhiCoordinates;

import javafx.geometry.Point3D;
import javafx.scene.paint.Material;
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
public class AutomataVisualization {

    private static Point3D       CANONICAL_Y_AXIS = new Point3D(0, 1, 0);

    private static final float[] NEGATIVE_TET     = new float[] { QUARTER_PI,
                                                                  THREE_QUARTERS_PI,
                                                                  THREE_QUARTERS_PI,
                                                                  QUARTER_PI,
                                                                  QUARTER_PI,
                                                                  THREE_QUARTERS_PI };

    private static final float[] POSITIVE_TET     = new float[] { THREE_QUARTERS_PI,
                                                                  QUARTER_PI,
                                                                  QUARTER_PI,
                                                                  THREE_QUARTERS_PI,
                                                                  THREE_QUARTERS_PI,
                                                                  QUARTER_PI };

    public static float[] getNegativeTet() {
        return Arrays.copyOf(NEGATIVE_TET, NEGATIVE_TET.length);
    }

    public static float[] getPositiveTet() {
        return Arrays.copyOf(POSITIVE_TET, POSITIVE_TET.length);
    }

    private static Rotate base(Point3D yAxis) {
        Point3D axisOfRotation = yAxis.crossProduct(CANONICAL_Y_AXIS);
        double angle = Math.acos(yAxis.normalize()
                                      .dotProduct(CANONICAL_Y_AXIS));
        return new Rotate(-Math.toDegrees(angle), axisOfRotation);
    }

    private final float          angularResolution;
    @SuppressWarnings("unused")
    private final Necronomata    automata;
    private final TriangleMesh[] exemplars;
    private final float          halfInterval;
    private final LengthTable    lengths;
    private final CellNode[]     nodes;

    public AutomataVisualization(int resolution, float radius,
                                 Necronomata automata, Material[] materials) {
        assert resolution
               % 8 == 0 : "Angular resolution must be divisable by 8: "
                          + resolution;
        CubicGrid grid = new CubicGrid(Neighborhood.SIX,
                                       PhiCoordinates.Cubes[3],
                                       automata.getExtent().x);
        this.automata = automata;
        Point3i extent = automata.getExtent();
        Point3i halfExtent = new Point3i(extent.x / 2, extent.y / 2,
                                         extent.z / 2);
        angularResolution = TWO_PI / resolution;
        exemplars = new TriangleMesh[resolution / 8];
        halfInterval = (float) (PhiCoordinates.Cubes[0].getEdgeLength() / 2.0);
        lengths = new LengthTable(resolution,
                                  PhiCoordinates.Cubes[0].getEdgeLength()
                                              * ROOT_2_DIV_2);
        for (int i = 0; i < exemplars.length; i++) {
            exemplars[i] = Line.createLine(lengths.length(i * angularResolution)
                                           * 2, radius);
        }

        List<CellNode> n = new ArrayList<>();
        automata.forEach(c -> {
            CellNode cell = cellNode(automata.anglesOf(c), materials);
            grid.postition(c.x - halfExtent.x, c.y - halfExtent.y,
                           c.z - halfExtent.z, cell);
            n.add(cell);
        });
        nodes = n.toArray(new CellNode[n.size()]);
    }

    public CellNode[] getNodes() {
        return nodes;
    }

    public CellNode cellNode(float[] initialState, Material[] materials) {
        Translate xPos = new Translate(0, 0, halfInterval);
        Translate xNeg = new Translate(0, 0, -halfInterval);
        Translate yPos = new Translate(halfInterval, 0, 0);
        Translate yNeg = new Translate(-halfInterval, 0, -0);
        MeshView[][] struts = new MeshView[6 * 5][];
        for (int c = 0; c < PhiCoordinates.Cubes.length; c++) {
            for (int i = 0; i < 6; i++) {
                int cube_ish = c * 6;
                struts[cube_ish + i] = new MeshView[exemplars.length * 8];
            }
        }
        for (int c = 0; c < PhiCoordinates.Cubes.length; c++) {
            Rotate base = base(yAxis(PhiCoordinates.Cubes[c]));
            for (int i = 0; i < exemplars.length * 8; i++) {
                MeshView view;
                view = xExemplar(i, xPos, base);
                view.setMaterial(materials[0]);
                int cube_ish = c * 6;
                struts[cube_ish][i] = view;

                view = xExemplar(i, xNeg, base);
                view.setMaterial(materials[1]);
                struts[cube_ish + 1][i] = view;

                view = yExemplar(i, yPos, base);
                view.setMaterial(materials[2]);
                struts[cube_ish + 2][i] = view;

                view = yExemplar(i, yNeg, base);
                view.setMaterial(materials[3]);
                struts[cube_ish + 3][i] = view;

                view = zExemplar(i, yPos, base);
                view.setMaterial(materials[4]);
                struts[cube_ish + 4][i] = view;

                view = zExemplar(i, yNeg, base);
                view.setMaterial(materials[5]);
                struts[cube_ish + 5][i] = view;
            }
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

    private MeshView xExemplar(int increment, Translate translate,
                               Rotate base) {
        Point3D axisOfRotation = new Point3D(0, 0, 1);
        Rotate rotateAroundCenter = new Rotate(-Math.toDegrees(increment
                                                               * angularResolution),
                                               axisOfRotation);
        MeshView line = new MeshView(x(increment, exemplars));

        line.getTransforms()
            .add(base.createConcatenation(rotateAroundCenter.createConcatenation(translate)));
        return line;
    }

    private MeshView yExemplar(int increment, Translate translate,
                               Rotate base) {
        Point3D axisOfRotation = new Point3D(1, 0, 0);
        Rotate rotateAroundCenter = new Rotate(-Math.toDegrees(increment
                                                               * angularResolution),
                                               axisOfRotation);
        MeshView line = new MeshView(x(increment, exemplars));

        line.getTransforms()
            .add(base.createConcatenation(rotateAroundCenter.createConcatenation(translate)));
        return line;
    }

    private MeshView zExemplar(int increment, Translate translate,
                               Rotate base) {
        Point3D axisOfRotation = new Point3D(1, 0, 0);
        Transform rotateAroundCenter = new Rotate(-Math.toDegrees(increment
                                                                  * angularResolution),
                                                  axisOfRotation);
        rotateAroundCenter = new Rotate(90, new Point3D(0, 0,
                                                        1)).createConcatenation(rotateAroundCenter);
        MeshView line = new MeshView(x(increment, exemplars));

        line.getTransforms()
            .add(base.createConcatenation(rotateAroundCenter.createConcatenation(translate)));
        return line;
    }
}
