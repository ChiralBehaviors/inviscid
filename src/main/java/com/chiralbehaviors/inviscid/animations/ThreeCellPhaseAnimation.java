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
import com.chiralbehaviors.inviscid.jitterbug.ReducedCoordinates.Assembly;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.AnimationTimer;
import javafx.scene.Group;
import javafx.scene.shape.CullFace;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;

/**
 * THREE CELLS WITH THEIR OWN PHASES: the same V that
 * {@link ThreeCellAnimation} builds, but integrated instead of driven.
 *
 * <p>
 * {@link ThreeCellAnimation} sweeps ONE coherent angle — middle cell at
 * {@code b = a + 60}, both outer cells at {@code a}, separations from the
 * shared-face law. That is a kinematic demonstration and a good one: it is how
 * the reciprocal condition and the {@code 4 -> 2 -> 1} shared-vertex decay were
 * found. But a single angle cannot represent a state in which one cell differs
 * from another, so it has exactly one degree of freedom by construction and
 * nothing can propagate through it. Anything that looked like a wave there was a
 * property of the parameterisation.
 *
 * <p>
 * Here each cell carries its own seven generalized coordinates — centre,
 * orientation, fold angle — and the array is integrated under
 * {@link ReducedCoordinates}. The V has <b>three</b> internal degrees of
 * freedom, one fold angle per cell, so a disturbance has somewhere to be. The
 * initial condition is a phase-rate kick on ONE cell, projected onto the
 * admissible set; watch it reach the others.
 *
 * <p>
 * <b>What to look for, and what not to conclude.</b>
 *
 * <ul>
 * <li>The kicked cell does not fold alone. It cannot: one cell cannot fold
 * without its neighbours folding, so the projection that makes the kick
 * admissible loads them <em>before any time passes</em>. There is no onset lag
 * and no wavefront here — a rigid constraint has infinite signal speed, so
 * nothing on screen is a signal speed.</li>
 * <li>What DOES take time is the redistribution of energy, printed below.</li>
 * <li>The cells also turn and drift. That is not a rendering artefact: with
 * cells free to rotate, an array of N cells has N internal freedoms, and the
 * orientation is how the fold angles are allowed to differ at all.</li>
 * <li>The fold angles <b>wind without bound</b> instead of oscillating. V = 0 in
 * this model — there is no potential energy anywhere — so nothing restores a
 * cell's phase. That is the open problem of the project on screen, not a
 * defect.</li>
 * </ul>
 *
 * <p>
 * <b>{@link #TIME_SCALE} is free, and exactly free.</b> With V = 0 the motion is
 * a geodesic, so scaling the initial velocity by {@code lambda} rescales time by
 * {@code 1/lambda} and leaves the trajectory through configuration space
 * identical. Slowing the animation therefore changes the clock and nothing else.
 * It is not a tuning knob on a dynamical result, because there is no dynamical
 * result it could tune.
 *
 * <p>
 * <b>Geometry comes from {@link JitterbugGeometry}, not from
 * {@link com.chiralbehaviors.inviscid.Jitterbug}.</b> The JavaFX renderer uses
 * {@code sigma = -(sx sy sz)} where the analysis harness uses
 * {@code sigma = +(sx sy sz)}; the two families agree as SETS but not pointwise,
 * so a cell rendered through the renderer is the mirror of the cell the
 * integrator computes, and mirrored cells do not close their welds unless the
 * centres and orientations are mirrored in lock-step. Rather than bridge the
 * conventions, this writes world positions {@code c + R v(gamma)} straight into
 * a mesh per triangle — the same function the integrator uses, and no transform
 * stack to get the order wrong in.
 *
 * @author halhildebrand
 */
public class ThreeCellPhaseAnimation extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            ThreeCellPhaseAnimation.main(argv);
        }
    }

    /** Which cell gets the phase-rate kick. 0 is the middle cell, whose kick is
     *  symmetric between the two outer cells and shows the constraint response
     *  cleanly; 1 is an outer cell, which breaks the symmetry and is the one that
     *  shows a disturbance travelling through the middle to the far cell. */
    private static final int    IMPULSE_CELL = 1;

    /** Runge-Kutta substeps per rendered frame. */
    private static final int    SUBSTEPS     = 4;

    /** See the class note: free, and exactly free, because V = 0. */
    private static final double TIME_SCALE   = 0.6;

    /** Two adjacent {@code <111>} diagonals, so the outer cells sit at
     *  {@code n1 . n2 = 1/3} — the same V {@link ThreeCellAnimation} builds. */
    private static final int[][] V_SITES     = { { 1, 1, 1 }, { 1, 1, -1 } };

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group root = new Group();

        final Assembly v = ReducedCoordinates.cluster(0.0, V_SITES);
        final int n = v.cells();

        // The viewer is framed for PhiCoordinates' octahedra; the model is unit
        // circumradius. One uniform factor, applied where the mesh is written.
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
                view.setMaterial(materials[f % materials.length]);
                // The cell opens: both sides of every triangle become visible.
                view.setCullFace(CullFace.NONE);
                mesh[k][f] = m;
                root.getChildren().add(view);
            }
        }
        content.setContent(root);

        final double[] state = new double[ReducedCoordinates.NQ * n
                                          + ReducedCoordinates.NU * n];
        final double[] q0 = v.initialState();
        final double[] u0 = v.foldImpulse(IMPULSE_CELL);
        System.arraycopy(q0, 0, state, 0, ReducedCoordinates.NQ * n);
        System.arraycopy(u0, 0, state, ReducedCoordinates.NQ * n,
                         ReducedCoordinates.NU * n);
        final double e0 = total(v.energyPerCell(q0, u0));

        final Runnable draw = () -> {
            double[] q = new double[ReducedCoordinates.NQ * n];
            System.arraycopy(state, 0, q, 0, q.length);
            double[][][] x = v.positions(q);
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

        System.out.printf("three cells in reduced coordinates: %d internal DOF, "
                          + "phase kick on cell %d, E = %.9f%n"
                          + "  fold rates at t=0: %s%n"
                          + "  V = 0, so the fold angles wind rather than "
                          + "oscillate and TIME_SCALE is free.%n",
                          n, IMPULSE_CELL, e0, rates(u0, n));

        new AnimationTimer() {
            private long   last  = 0;
            private long   spoke = 0;
            private double t     = 0;
            private double worst = 0;

            @Override
            public void handle(long now) {
                if (last == 0) {
                    last = now;
                    return;
                }
                // Clamped so a stalled frame cannot take one enormous step.
                double dt = Math.min((now - last) / 1e9, 0.05) * TIME_SCALE;
                last = now;
                v.step(state, dt, SUBSTEPS);
                t += dt;
                draw.run();

                double[] q = new double[ReducedCoordinates.NQ * n];
                double[] u = new double[ReducedCoordinates.NU * n];
                System.arraycopy(state, 0, q, 0, q.length);
                System.arraycopy(state, q.length, u, 0, u.length);
                double[] e = v.energyPerCell(q, u);
                double sum = total(e);
                worst = Math.max(worst, Math.abs(sum - e0) / e0);
                if (System.currentTimeMillis() - spoke > 500) {
                    spoke = System.currentTimeMillis();
                    StringBuilder sb = new StringBuilder();
                    for (int k = 0; k < n; k++) {
                        sb.append(String.format("  cell%d g=%+9.2f share=%6.4f", k,
                                                q[ReducedCoordinates.NQ * k + 7],
                                                e[k] / sum));
                    }
                    System.out.printf("t=%7.2f%s | E drift %8.2e | weld %8.2e%n", t,
                                      sb, worst, v.weldResidual(q));
                }
            }
        }.start();
    }

    private static String rates(double[] u, int n) {
        StringBuilder sb = new StringBuilder();
        for (int k = 0; k < n; k++) {
            sb.append(String.format("%s%.6f", k == 0 ? "" : ", ",
                                    u[ReducedCoordinates.NU * k + 6]));
        }
        return sb.toString();
    }

    private static double total(double[] v) {
        double s = 0;
        for (double x : v) {
            s += x;
        }
        return s;
    }
}
