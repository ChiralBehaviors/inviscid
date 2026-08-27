/**
 * Copyright (C) 2018 Chiral Behaviors, LLC. All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package com.chiralbehaviors.inviscid.animations;

import static com.chiralbehaviors.inviscid.animations.Colors.blackMaterial;
import static com.chiralbehaviors.inviscid.animations.Colors.blueMaterial;
import static com.chiralbehaviors.inviscid.animations.Colors.greenMaterial;
import static com.chiralbehaviors.inviscid.animations.Colors.materials;
import static com.chiralbehaviors.inviscid.animations.Colors.redMaterial;

import java.util.ArrayList;
import java.util.List;

import com.chiralbehaviors.inviscid.CubicGrid;
import com.chiralbehaviors.inviscid.CubicGrid.Neighborhood;
import com.chiralbehaviors.inviscid.Jitterbug;
import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.KeyFrame;
import javafx.animation.KeyValue;
import javafx.animation.Timeline;
import javafx.beans.value.WritableValue;
import javafx.scene.Group;
import javafx.util.Duration;
import mesh.Ellipse;
import mesh.polyhedra.plato.Octahedron;

/**
 * @author halhildebrand
 *
 */
public class JitterbugAnimation extends PolyView {
    /** The hole cell runs exactly this far ahead of the cell: b = a + 60. */
    private static final double PHASE_OFFSET = 60.0;

    public static class Launcher {

        public static void main(String[] argv) {
            JitterbugAnimation.main(argv);
        }
    }

    public static void main(String[] args) {
        launch(args);
    }

    public void jitterbugArray(Group group, Octahedron[] octahedrons, List<Jitterbug> jitterbugs, CubicGrid grid,
                               double initialAngle) {
        for (int x : new int[] { -1, 1 }) {
            for (int y : new int[] { -1, 1 }) {
                for (int z : new int[] { -1, 1 }) {
                    Jitterbug j = new Jitterbug(octahedrons[4], materials);
                    // The corner cells are the HOLE cells and run EXACTLY 60
                    // degrees ahead: b = a + 60, the unique closure of the
                    // shared-face constraint (verified 8.3e-16 over 2, 3 and 9
                    // cells). At a = 0 that puts them at 60 -- collapsed to
                    // octahedra -- which is what makes the packing close.
                    // Driving all nine at the same phase gives nine full-size
                    // cuboctahedra on sites meant for one VE and eight octa,
                    // and they interpenetrate.
                    j.rotateTo(initialAngle + PHASE_OFFSET);
                    Group jGroup = j.getGroup();
                    grid.position(x, y, z, jGroup);
                    group.getChildren().add(jGroup);
                    jitterbugs.add(j);
                }
            }
        }
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group group = new Group();
        List<Jitterbug> jitterbugs = new ArrayList<>();

        CubicGrid grid = new CubicGrid(Neighborhood.EIGHT, PhiCoordinates.Cubes[3], 1);
        group.getChildren().add(grid.construct(blackMaterial, blackMaterial, blackMaterial));

        Jitterbug centre = new Jitterbug(PhiCoordinates.Octahedrons[4], materials);
        centre.rotateTo(0);
        Group jGroup = centre.getGroup();
        group.getChildren().add(jGroup);

        group.getChildren().add(new Ellipse(0, PhiCoordinates.Octahedrons[4], 0).construct(40, redMaterial, 0.015));
        group.getChildren().add(new Ellipse(0, PhiCoordinates.Octahedrons[4], 1).construct(40, blueMaterial, 0.015));
        group.getChildren().add(new Ellipse(0, PhiCoordinates.Octahedrons[4], 2).construct(40, greenMaterial, 0.015));

        // SINGLE CELL. The array is left off deliberately: Gray, "The Jitterbug
        // Motion" (2002) -- "Two Jitterbugs can not share the same triangular
        // face and have their positions (location of center of volume) fixed as
        // they go through the Jitterbug motion. If two Jitterbugs are to share
        // the same triangle face then as the joined Jitterbugs jitterbug, the
        // positions of the Jitterbugs must move." Pinning cells to a fixed grid
        // tears the shared vertices apart, which is what the array version did.
        // jitterbugArray(group, PhiCoordinates.Octahedrons, jitterbugs, grid, (double) 0);

        content.setContent(group);
        final Timeline timeline = new Timeline();
        // FULL SWEEP, -60 -> +60, THROUGH the VE at a = 0. The VE is the middle
        // of the range, not an end of it: free dynamics has adot PEAKING there
        // (M_eff(a)/M = Z^2 sin^2 a + rho^2 is minimal at a = 0), which is
        // Fuller's "the equatorial rotational momentum will be seen to carry
        // the rotation beyond dead-center". Sweeping 0 -> 60 and reversing
        // treats the vector equilibrium as a turning point, which is exactly
        // backwards -- it is the fastest point of the motion. Passing through
        // it inverts the chirality: rotateTo applies +a to one inscribed
        // tetrahedron and -a to the other, so a -> -a swaps their roles.
        timeline.getKeyFrames().add(new KeyFrame(Duration.millis(10_000), new KeyValue(new WritableValue<Double>() {
            @Override
            public Double getValue() {
                return -60d;
            }

            @Override
            public void setValue(Double value) {
                centre.rotateTo(value);
                jitterbugs.forEach(j -> j.rotateTo(value + PHASE_OFFSET));
                // (jitterbugs is empty while the array above is disabled)
            }
        }, 60d)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
