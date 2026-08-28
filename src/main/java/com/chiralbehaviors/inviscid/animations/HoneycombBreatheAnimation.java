/**
 * Copyright (c) 2026 Chiral Behaviors, LLC, all rights reserved.
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

import static com.chiralbehaviors.inviscid.animations.Colors.materials;

import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry;
import com.chiralbehaviors.inviscid.jitterbug.ReducedCoordinates;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.AnimationTimer;
import javafx.scene.Group;
import javafx.scene.shape.CullFace;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;

/**
 * THE MEDIUM'S ONE MOTION: a compact patch of the honeycomb breathing, bounded
 * by its own geometry.
 *
 * <p>
 * <b>Why this and not a wave.</b> The wave programme's premise was that the
 * medium is a field — one fold angle per cell, so a disturbance has somewhere to
 * be. That was measured on a chain and on a single-shell cluster, and both are
 * <b>trees</b>. The rectified cubic honeycomb is richly cyclic, and closing the
 * cycles takes the freedom away: a compact patch of 113 cells has <b>one</b>
 * internal degree of freedom, and the minimal four-cycle locks all four cells to
 * the same fold rate. So this is not a wave medium as modelled, and the
 * animation shows what it actually does instead of what it was hoped to do.
 *
 * <p>
 * A phase kick on a <em>single</em> cell of this patch produces the motion
 * below, because there is nothing else available for it to produce.
 *
 * <p>
 * <b>What is on screen.</b> The two sublattices exchange roles: at
 * {@code a = 0} the even cells are cuboctahedra and the odd ones octahedra, and
 * at {@code a = -60} they have completely swapped. The lattice constant breathes
 * {@code 1 -> 2/sqrt(3) -> 1} on the way. Nothing tumbles — no cell rotates at
 * all, to 3e-12 — and nothing passes through itself, because {@code |gamma| = 60}
 * is where a cell folds through its own octahedron and is a hard stop rather
 * than a convention. {@code b = a + 60} forces both sublattices into that window
 * at once, which is what makes the range exactly sixty degrees.
 *
 * <p>
 * <b>It is dynamics, not a sweep, and the difference is visible.</b>
 * {@link ThreeCellAnimation} drives its coherent angle at a constant rate. Here
 * the angle obeys {@code addot = -(1/2)(M'/M) adot^2} with the patch's own
 * effective mass, which runs 272 → 120 → 264 across the window — so the breathe
 * is <b>fastest at mid-swing and slowest at the ends</b>, and it lingers where
 * the cells are most open. The bounces at the limits are elastic and exact:
 * with one degree of freedom, reversing the velocity conserves energy
 * identically. Watch the printed energy hold across them.
 *
 * <p>
 * <b>No constraint solve happens here at all.</b> The motion is closed form —
 * {@code gamma_even = a}, {@code gamma_odd = a + 60}, {@code R = I},
 * {@code centre = L(a) * site} — and that is verified against the full
 * constrained model rather than assumed, over five angles, by
 * {@code jb_rc_reduced.py}'s R5h and by {@code ReducedCoordinatesTest}.
 *
 * <p>
 * Geometry comes from {@link JitterbugGeometry} rather than from
 * {@link com.chiralbehaviors.inviscid.Jitterbug}: the JavaFX renderer uses the
 * opposite {@code sigma}, so a cell drawn through it is the mirror of the cell
 * the model computes, and mirrored cells do not close their welds.
 *
 * @author halhildebrand
 */
public class HoneycombBreatheAnimation extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            HoneycombBreatheAnimation.main(argv);
        }
    }

    /** Where the swing starts. {@code -30} is mid-window, and the one angle at
     *  which the lattice constant is stationary. */
    private static final double A0        = -30.0;

    /** Initial fold rate, radians per unit time. With V = 0 this sets the clock
     *  and nothing else: scaling it rescales time exactly and leaves the path
     *  through configuration space identical. */
    private static final double ADOT0     = 0.30;

    /** Radius of the compact patch. 2.0 gives 15 cells and 32 welds with NO
     *  dangling cells, hence exactly one internal degree of freedom. Raise it
     *  for more of the lattice; the motion does not change, because there is
     *  only ever the one. */
    private static final double RADIUS    = 2.0;

    private static final int    SUBSTEPS  = 4;

    private static final double TIME_SCALE = 0.5;

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group root = new Group();

        final int[][] sites = ReducedCoordinates.ball(RADIUS);
        final int n = sites.length;
        final double scale = PhiCoordinates.Octahedrons[4].getEdgeLength()
                             / JitterbugGeometry.L_EDGE;

        final TriangleMesh[][] mesh = new TriangleMesh[n][8];
        for (int k = 0; k < n; k++) {
            for (int f = 0; f < 8; f++) {
                TriangleMesh m = new TriangleMesh();
                m.getPoints().addAll(new float[9]);
                m.getTexCoords().addAll(0, 0);
                m.getFaces().addAll(0, 0, 1, 0, 2, 0);
                MeshView view = new MeshView(m);
                // Colour by SUBLATTICE, not by cell, so the exchange is the
                // thing you see: one family opens as the other closes.
                boolean even = Math.floorMod(sites[k][0], 2) == 0
                               && Math.floorMod(sites[k][1], 2) == 0
                               && Math.floorMod(sites[k][2], 2) == 0;
                view.setMaterial(materials[(even ? 0 : 4) + f % 4]);
                view.setCullFace(CullFace.NONE);
                mesh[k][f] = m;
                root.getChildren().add(view);
            }
        }
        content.setContent(root);

        final double[] state = { A0, ADOT0 };
        final double e0 = ReducedCoordinates.swingEnergy(state, sites);

        final Runnable draw = () -> {
            double[][][] x = ReducedCoordinates.coherentPositions(sites, state[0]);
            for (int k = 0; k < n; k++) {
                for (int f = 0; f < 8; f++) {
                    float[] p = new float[9];
                    for (int c = 0; c < 3; c++) {
                        double[] w = x[k][ReducedCoordinates.SLOT[f][c]];
                        for (int t = 0; t < 3; t++) {
                            p[3 * c + t] = (float) (w[t] * scale);
                        }
                    }
                    mesh[k][f].getPoints().set(0, p, 0, 9);
                }
            }
        };
        draw.run();

        System.out.printf("a compact honeycomb patch: %d cells, ONE internal "
                          + "degree of freedom%n"
                          + "  the sublattices exchange over a in [%.0f, 0]; "
                          + "M_eff %.1f -> %.1f -> %.1f, so the breathe is "
                          + "fastest at mid-swing%n"
                          + "  V = 0, so E is the only audit -- watch it hold "
                          + "ACROSS the elastic bounces%n",
                          n, -ReducedCoordinates.FOLD_LIMIT,
                          ReducedCoordinates.effectiveMass(sites, 0.0)[0],
                          ReducedCoordinates.effectiveMass(sites, -30.0)[0],
                          ReducedCoordinates.effectiveMass(sites, -60.0)[0]);

        new AnimationTimer() {
            private int    bounces = 0;
            private long   last    = 0;
            private long   spoke   = 0;
            private double t       = 0;
            private double worst   = 0;

            @Override
            public void handle(long now) {
                if (last == 0) {
                    last = now;
                    return;
                }
                double dt = Math.min((now - last) / 1e9, 0.05) * TIME_SCALE;
                last = now;
                bounces += ReducedCoordinates.swingStep(state, sites, dt, SUBSTEPS);
                t += dt;
                draw.run();
                double e = ReducedCoordinates.swingEnergy(state, sites);
                worst = Math.max(worst, Math.abs(e - e0) / e0);
                if (System.currentTimeMillis() - spoke > 500) {
                    spoke = System.currentTimeMillis();
                    System.out.printf("t=%7.2f  a=%+8.3f  b=%+8.3f  adot=%+8.4f  "
                                      + "L=%8.6f | E=%12.9f  worst drift %8.2e | "
                                      + "bounces %d%n",
                                      t, state[0], state[0] + 60.0, state[1],
                                      ReducedCoordinates.latticeDerivatives(state[0])[0],
                                      e, worst, bounces);
                }
            }
        }.start();
    }
}
