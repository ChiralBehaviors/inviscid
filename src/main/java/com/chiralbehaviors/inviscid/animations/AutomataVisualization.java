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

import static com.chiralbehaviors.inviscid.Constants.*;
import static com.chiralbehaviors.inviscid.Constants.ROOT_2_DIV_2;
import static com.chiralbehaviors.inviscid.Constants.THREE_QUARTERS_PI;
import static com.chiralbehaviors.inviscid.Constants.TWO_PI;
import static com.chiralbehaviors.inviscid.CubicGrid.yAxis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.chiralbehaviors.inviscid.CubicGrid;
import com.chiralbehaviors.inviscid.CubicGrid.Neighborhood;
import com.chiralbehaviors.inviscid.LengthTable;
import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.PhiCoordinates;

import javafx.collections.ObservableList;
import javafx.geometry.Point3D;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.paint.Material;
import javafx.scene.shape.Mesh;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Scale;
import javafx.scene.transform.Transform;
import javafx.scene.transform.Translate;
import mesh.Line;

/**
 * @author halhildebrand
 *
 */
public class AutomataVisualization extends Group {

    class CellAnimator {
        private final MeshView[][] struts;

        public CellAnimator(MeshView[][] struts) {
            this.struts = struts;
        }

        public List<MeshView> allStuts() {
            List<MeshView> allStruts = new ArrayList<>();
            for (MeshView strut : struts[3]) {
                allStruts.add(strut);
            }
            return allStruts;
        }

        public void update(int cell) {
            automata.process((angle, frequency, deltaA, deltaF) -> {
                setState(angle, cell, struts);
            });
        }
    }

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

    private static final int     STRUTS_PER_CELL  = 6 * 5;

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

    private final Necronomata        automata;
    private final List<CellAnimator> cells             = new ArrayList<>();
    private final Transform[]        lengths;
    private float                    angularResolution = TWO_PI;
    private final Transform[][]      rotations;

    public AutomataVisualization(int resolution, float radius,
                                 Necronomata automata, Material[] materials) {
        assert resolution
               % 8 == 0 : "Angular resolution must be divisable by 8: "
                          + resolution;
        angularResolution = TWO_PI / resolution;
        CubicGrid grid = new CubicGrid(Neighborhood.SIX,
                                       PhiCoordinates.Cubes[3],
                                       automata.getExtent().x);
        this.automata = automata;
        float halfInterval = (float) (PhiCoordinates.Cubes[0].getEdgeLength()
                                      / 2.0);
        TriangleMesh exemplar = Line.createLine(PhiCoordinates.Cubes[0].getEdgeLength()
                                                * ROOT_2, radius);
        rotations = new Transform[3][];
        for (int i = 0; i < 3; i++) {
            rotations[i] = new Transform[resolution];
        }
        lengths = new Transform[resolution / 8];
        buildRotations(resolution);
        LengthTable table = new LengthTable(resolution);
        for (int i1 = 0; i1 < resolution / 8; i1++) {
            lengths[i1] = new Scale(table.lengthAt(i1), 1.0, 1.0);
        }
        createCellAnimators(grid, materials, exemplar, halfInterval);
        ObservableList<Node> children = getChildren();
        cells.forEach(cell -> children.addAll(cell.allStuts()));
    }

    public void update() {
        for (int i = 0; i < cells.size(); i++) {
            cells.get(i)
                 .update(i);
        }
        setNeedsLayout(true);
        layout();
    }

    private void buildRotations(int resolution) {
        Point3D axisOfRotation = new Point3D(0, 0, 1);
        for (int i = 0; i < resolution; i++) {
            rotations[0][i] = new Rotate(-Math.toDegrees(i * angularResolution),
                                         axisOfRotation);
        }
        axisOfRotation = new Point3D(1, 0, 0);
        for (int i = 0; i < resolution; i++) {
            rotations[1][i] = new Rotate(-Math.toDegrees(i * angularResolution),
                                         axisOfRotation);
        }
        for (int i = 0; i < resolution; i++) {
            Transform rotateAroundCenter = new Rotate(-Math.toDegrees(i
                                                                      * angularResolution),
                                                      axisOfRotation);
            rotations[2][i] = new Rotate(90, new Point3D(0, 0,
                                                         1)).createConcatenation(rotateAroundCenter);
        }
    }

    private void createCellAnimators(CubicGrid grid, Material[] materials,
                                     Mesh exemplar, float halfInterval) {
        Translate xPos = new Translate(0, 0, halfInterval);
        Translate xNeg = new Translate(0, 0, -halfInterval);
        Translate yPos = new Translate(halfInterval, 0, 0);
        Translate yNeg = new Translate(-halfInterval, 0, -0);

        automata.forEach(location -> {
            MeshView[][] cell = new MeshView[5][];
            for (int cube = 0; cube < 5; cube++) {
                cell[cube] = new MeshView[6];
                Rotate base = base(yAxis(PhiCoordinates.Cubes[cube]));
                MeshView view;
                Transform position = grid.postitionTransform(location.x,
                                                             location.y,
                                                             location.z);
                view = exemplar(0, xPos, base, exemplar, position);
                view.setMaterial(materials[0]);
                cell[cube][0] = view;

                view = exemplar(0, xNeg, base, exemplar, position);
                view.setMaterial(materials[1]);
                cell[cube][1] = view;

                view = exemplar(1, yPos, base, exemplar, position);
                view.setMaterial(materials[2]);
                cell[cube][2] = view;

                view = exemplar(1, yNeg, base, exemplar, position);
                view.setMaterial(materials[3]);
                cell[cube][3] = view;

                view = exemplar(2, yPos, base, exemplar, position);
                view.setMaterial(materials[4]);
                cell[cube][4] = view;

                view = exemplar(2, yNeg, base, exemplar, position);
                view.setMaterial(materials[5]);
                cell[cube][5] = view;
            }
            cells.add(new CellAnimator(cell));
        });

    }

    private MeshView exemplar(int axis, Translate translate, Rotate base,
                              Mesh exemplar, Transform position) {
        MeshView line = new MeshView(exemplar);
        line.getTransforms()
            .addAll(scale(0).clone(), base, rotations[axis][0].clone(),
                    translate, position);
        return line;
    }

    private Transform scale(int increment) {
        if (increment < lengths.length) {
            return lengths[increment];
        }
        if (increment < lengths.length * 2) {
            return lengths[lengths.length - (increment - lengths.length) - 1];
        }
        if (increment < lengths.length * 3) {
            return lengths[increment - (lengths.length * 2)];
        }
        if (increment < lengths.length * 4) {
            return lengths[lengths.length - (increment - lengths.length * 3)
                           - 1];
        }
        if (increment < lengths.length * 5) {
            return lengths[increment - (lengths.length * 4)];
        }
        if (increment < lengths.length * 6) {
            return lengths[lengths.length - (increment - lengths.length * 5)
                           - 1];
        }
        if (increment < lengths.length * 7) {
            return lengths[increment - (lengths.length * 6)];
        }
        return lengths[lengths.length - (increment - lengths.length * 7) - 1];
    }

    private void setState(float[] angle, int cell, MeshView[][] struts) {
        System.out.println(angle[0]);
        int index = cell * STRUTS_PER_CELL;
        for (int cube = 0; cube < 5; cube++) {
            for (int side = 0; side < 3; side++) {
                for (int face = 0; face < 2; face++) {
                    int step = (int) (angle[index++] / angularResolution);
                    MeshView strut = struts[cube][(side * 2) + face];
                    List<Transform> transforms = strut.getTransforms();
                    transforms.set(0, scale(step));
                    transforms.set(2, rotations[side][step]);
                }
            }
        }
    }

}
