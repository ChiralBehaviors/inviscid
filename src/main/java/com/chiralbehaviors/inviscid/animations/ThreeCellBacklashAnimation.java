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
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.AnimationTimer;
import javafx.scene.Group;
import javafx.scene.shape.CullFace;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;

/**
 * THREE CELLS ON A DOWEL, WITH PLAY IN THE JOINTS — the 2026-08-28 result made
 * visible, so it can be checked against the physical rig rather than believed.
 *
 * <p>
 * <b>Why this exists.</b> {@link ThreeCellPhaseAnimation} shows the RIGID model
 * and says so plainly: "There is no onset lag and no wavefront here — a rigid
 * constraint has infinite signal speed." That was correct, and it is why the
 * wave programme kept failing. This animation shows what changes when the joints
 * are given CLEARANCE, which is the one mechanism nobody had modelled and which
 * the owner named on 2026-08-28: the physical cells hold station "due to the
 * very low tolerance of the physical build".
 *
 * <p>
 * <b>The model, which is exactly the harness's</b>
 * ({@code analysis/jitterbug-variety/jb_ct_contact_chain.py}, 6 gated rows).
 * Three cells at FIXED sites along one body diagonal, phases alternating
 * {@code a, a+60, a}. Each cell has one fold angle and its own corner inertia
 * {@code m(g) = (16/3)(1 + 2 sin^2 g)}. Gray's radial law puts a cell's shared
 * face at {@code V(g) = Z cos(g)} from its centre, so two neighbours at a fixed
 * spacing are coupled by
 *
 * <pre>
 *     c_k = V(g_k) + V(g_k+1) - sep0,     |c_k| &lt;= t
 * </pre>
 *
 * INSIDE that band the cells are free of each other; at its two edges they
 * collide elastically. There is no spring anywhere — {@code V = 0} between
 * impacts — so this is BACKLASH in a permanently engaged joint, not a unilateral
 * contact that can separate, and not a compliant strut.
 *
 * <p>
 * <b>What to watch, and what each one would mean if the rig disagrees.</b>
 *
 * <ul>
 * <li><b>The centres never move.</b> Gray p.40 says two jitterbugs cannot share
 * a face and keep their positions fixed. This model holds them fixed anyway and
 * pays for it in clearance — that is the whole hypothesis. If the owner's cells
 * visibly walk apart as they fold, the premise is wrong and this line should
 * stop.</li>
 * <li><b>The joints rattle.</b> The gauge under each joint shows {@code c_k/t}
 * riding between {@code -1} and {@code +1} and striking the stops. If the rig's
 * joints are tight enough that nothing rattles, there is no play, and with no
 * play there is no finite signal speed (the harness measures the speed
 * diverging as the clearance closes).</li>
 * <li><b>The kick takes TIME to cross.</b> Only cell 0 is kicked. Cell 1 cannot
 * move until cell 0 has physically crossed its clearance, and cell 2 not until
 * cell 1 has crossed its own. The arrival times print as they happen. This is
 * the first finite propagation speed in the project, and the clearance is what
 * buys it.</li>
 * <li><b>The fold angles WIND rather than oscillate.</b> V = 0, so nothing
 * restores a cell's phase. That is the project's open problem on screen, not a
 * defect — and {@code jb_cp_contact_potential.py} measures why the usual repair
 * is unavailable: the effective potential here is an infinite square well, which
 * has no Hessian anywhere and therefore no normal modes.</li>
 * </ul>
 *
 * <p>
 * <b>{@link #PLAY} IS EXAGGERATED, AND THAT IS A PRESENTATION CHOICE.</b> A real
 * build's clearance is around 1% of an edge — 1 mm on a 100 mm strut — at which
 * scale nothing here would be visible. The default below is far larger so the
 * rattle can be seen at all. Every SCALING the harness measures is unaffected
 * (speed goes as 1/play, so a bigger play only slows the crossing), but no
 * quantity on this screen should be read as a measurement. The gated numbers
 * live in the Python harness; this is a picture of the mechanism.
 *
 * <p>
 * <b>Geometry comes from {@link JitterbugGeometry}</b>, for the reason
 * {@link ThreeCellPhaseAnimation} records: the JavaFX renderer uses
 * {@code sigma = -(sx sy sz)} where the analysis harness uses {@code +}, so a
 * cell drawn through the renderer is the mirror of the one the integrator
 * computes, and mirrored cells do not close their welds. World positions are
 * written straight into a mesh per triangle.
 *
 * @author halhildebrand
 */
public class ThreeCellBacklashAnimation extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            ThreeCellBacklashAnimation.main(argv);
        }
    }

    /** Cells on the dowel. Three is the smallest number with a cell that the
     *  kick has to reach THROUGH another one. */
    static final int    CELLS      = 3;

    /** The reference phase: midpoint of the exchange, widest lattice, and the
     *  one phase at which the harness measures a uniform front. */
    static final double A_REF      = -30.0;

    /** Joint clearance. EXAGGERATED for visibility — see the class note. A
     *  physical build is nearer 0.01 in these units. */
    static final double PLAY       = 0.12;

    /** Fold rate given to cell 0 at t = 0. Nothing else is moving. */
    static final double KICK       = 0.9;

    /** Integrator substeps per rendered frame. */
    private static final int    SUBSTEPS   = 8;

    /** Wall-clock slowdown. With V = 0 the motion is a geodesic between
     *  impacts and the impacts are elastic, so this changes the clock and
     *  nothing else. */
    private static final double TIME_SCALE = 0.25;

    /** A cell counts as "reached" when its fold rate exceeds this fraction of
     *  the kick. Printed, never gated — the gated arrival times are the
     *  harness's. */
    private static final double ARRIVED    = 0.02;

    private static final double M0         = 16.0 / 3.0;

    public static void main(String[] args) {
        launch(args);
    }

    private static double inertia(double g) {
        double s = Math.sin(g);
        return M0 * (1.0 + 2.0 * s * s);
    }

    private static double inertiaD(double g) {
        return M0 * 4.0 * Math.sin(g) * Math.cos(g);
    }

    /** Gray eq (3): centre of volume to a triangle's face centre. */
    static double radial(double g) {
        return JitterbugGeometry.Z * Math.cos(g);
    }

    private static double radialD(double g) {
        return -JitterbugGeometry.Z * Math.sin(g);
    }

    /** The face whose outward axis points most nearly along {@code d}. Plate
     *  normals are fixed in phase, so this is a phase-independent lookup. */
    static int faceToward(double[] d) {
        int best = 0;
        double bestDot = Double.NEGATIVE_INFINITY;
        for (int f = 0; f < 8; f++) {
            double[] u = JitterbugGeometry.faceAxis(f);
            double dot = u[0] * d[0] + u[1] * d[1] + u[2] * d[2];
            if (dot > bestDot) {
                bestDot = dot;
                best = f;
            }
        }
        return best;
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group root = new Group();

        final double scale = PhiCoordinates.Octahedrons[4].getEdgeLength()
                             / JitterbugGeometry.L_EDGE;

        // The dowel: one body diagonal. Cells sit sep0 apart along it, which is
        // exactly the distance at which their shared faces coincide at A_REF.
        final double norm = Math.sqrt(3.0);
        final double[] axis = { 1.0 / norm, 1.0 / norm, 1.0 / norm };
        final double sep0 = radial(Math.toRadians(A_REF))
                            + radial(Math.toRadians(A_REF + 60.0));

        final double[] g = new double[CELLS];
        final double[] gd = new double[CELLS];
        for (int k = 0; k < CELLS; k++) {
            g[k] = Math.toRadians(A_REF + 60.0 * (k % 2));
        }
        gd[0] = KICK;

        final TriangleMesh[][] mesh = new TriangleMesh[CELLS][8];
        for (int k = 0; k < CELLS; k++) {
            for (int f = 0; f < 8; f++) {
                TriangleMesh m = new TriangleMesh();
                m.getPoints().addAll(new float[9]);
                m.getTexCoords().addAll(0, 0);
                m.getFaces().addAll(0, 0, 1, 0, 2, 0);
                MeshView view = new MeshView(m);
                view.setMaterial(materials[(k * 3 + f) % materials.length]);
                view.setCullFace(CullFace.NONE);
                mesh[k][f] = m;
                root.getChildren().add(view);
            }
        }
        content.setContent(root);

        final Runnable draw = () -> {
            for (int k = 0; k < CELLS; k++) {
                double[][][] x = JitterbugGeometry.corners(Math.toDegrees(g[k]));
                double off = k * sep0;
                for (int f = 0; f < 8; f++) {
                    float[] p = new float[9];
                    for (int c = 0; c < 3; c++) {
                        for (int t = 0; t < 3; t++) {
                            p[3 * c + t] = (float) ((x[f][c][t] + off * axis[t])
                                                    * scale);
                        }
                    }
                    mesh[k][f].getPoints().set(0, p, 0, 9);
                }
            }
        };
        draw.run();

        final double e0 = energy(g, gd);
        System.out.printf("THREE CELLS ON A DOWEL, WITH PLAY IN THE JOINTS%n"
                          + "  phases %.0f / %.0f / %.0f deg, spacing %.6f, "
                          + "play %.3f (EXAGGERATED for visibility)%n"
                          + "  kick %.3f on cell 0 only; E = %.9f%n"
                          + "  watch: centres never move, joints rattle in the "
                          + "band, the kick takes TIME to cross,%n"
                          + "         and the fold angles wind because V = 0.%n"
                          + "  Nothing here is a measurement -- the gated "
                          + "numbers are in jb_ct_contact_chain.py.%n%n",
                          A_REF, A_REF + 60.0, A_REF, sep0, PLAY, KICK, e0);

        new AnimationTimer() {
            private long          last    = 0;
            private long          spoke   = 0;
            private double        clock   = 0;
            private int           hits    = 0;
            private final double[] arrive = new double[CELLS];
            private double        worst   = 0;

            {
                for (int k = 1; k < CELLS; k++) {
                    arrive[k] = Double.NaN;
                }
            }

            @Override
            public void handle(long now) {
                if (last == 0) {
                    last = now;
                    return;
                }
                double dt = Math.min((now - last) * 1e-9, 1.0 / 30.0)
                            * TIME_SCALE;
                last = now;
                double h = dt / SUBSTEPS;
                for (int s = 0; s < SUBSTEPS; s++) {
                    hits += advance(g, gd, sep0, h);
                    clock += h;
                    for (int k = 1; k < CELLS; k++) {
                        if (Double.isNaN(arrive[k])
                            && Math.abs(gd[k]) > ARRIVED * KICK) {
                            arrive[k] = clock;
                            System.out.printf("  cell %d reached at t = %.4f%n",
                                              k, clock);
                        }
                    }
                }
                for (int k = 0; k + 1 < CELLS; k++) {
                    worst = Math.max(worst, Math.abs(gap(g, k, sep0)) - PLAY);
                }
                draw.run();

                if (now - spoke > 500_000_000L) {
                    spoke = now;
                    StringBuilder sb = new StringBuilder();
                    for (int k = 0; k + 1 < CELLS; k++) {
                        sb.append(String.format("  joint %d %s", k,
                                                bar(gap(g, k, sep0) / PLAY)));
                    }
                    System.out.printf("t=%7.2f  contacts %5d  dE/E %8.1e%s%n",
                                      clock, hits,
                                      Math.abs(energy(g, gd) - e0) / e0,
                                      sb);
                }
            }
        }.start();
    }

    /** {@code c_k}: how far this joint sits from its centred position. */
    static double gap(double[] g, int k, double sep0) {
        return radial(g[k]) + radial(g[k + 1]) - sep0;
    }

    static double energy(double[] g, double[] gd) {
        double e = 0;
        for (int k = 0; k < g.length; k++) {
            e += 0.5 * inertia(g[k]) * gd[k] * gd[k];
        }
        return e;
    }

    /** A visual gauge for one joint: where it rides between its two stops. */
    private static String bar(double frac) {
        int slots = 21;
        int at = (int) Math.round((Math.max(-1, Math.min(1, frac)) + 1) / 2
                                  * (slots - 1));
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < slots; i++) {
            sb.append(i == at ? (Math.abs(frac) > 0.999 ? '#' : 'o') : '-');
        }
        return sb.append(']').toString();
    }

    /**
     * One event-driven step. Integrate; if a joint would leave its band, bisect
     * to the crossing and apply an elastic impulse along the constraint
     * gradient. Nothing is ever projected back, so energy is conserved and the
     * readout means something — the harness records that nudging a position
     * onto the stop instead bleeds 20 to 40 percent.
     */
    static int advance(double[] g, double[] gd, double sep0, double h) {
        double[] ng = g.clone();
        double[] ngd = gd.clone();
        rk4(ng, ngd, h);
        if (!violates(ng, sep0)) {
            System.arraycopy(ng, 0, g, 0, g.length);
            System.arraycopy(ngd, 0, gd, 0, gd.length);
            return 0;
        }
        double lo = 0, hi = h;
        for (int i = 0; i < 40; i++) {
            double mid = 0.5 * (lo + hi);
            double[] tg = g.clone();
            double[] tgd = gd.clone();
            rk4(tg, tgd, mid);
            if (violates(tg, sep0)) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        rk4(g, gd, lo);
        int k = 0;
        double worst = -1;
        for (int j = 0; j + 1 < g.length; j++) {
            double a = Math.abs(gap(g, j, sep0));
            if (a > worst) {
                worst = a;
                k = j;
            }
        }
        double c = gap(g, k, sep0);
        double u = radialD(g[k]) * gd[k] + radialD(g[k + 1]) * gd[k + 1];
        if ((c > 0 && u > 0) || (c < 0 && u < 0)) {
            double lam = -2.0 * u
                         / (radialD(g[k]) * radialD(g[k]) / inertia(g[k])
                            + radialD(g[k + 1]) * radialD(g[k + 1])
                              / inertia(g[k + 1]));
            gd[k] += lam * radialD(g[k]) / inertia(g[k]);
            gd[k + 1] += lam * radialD(g[k + 1]) / inertia(g[k + 1]);
            return 1;
        }
        return 0;
    }

    private static boolean violates(double[] g, double sep0) {
        for (int k = 0; k + 1 < g.length; k++) {
            if (Math.abs(gap(g, k, sep0)) - PLAY > 0) {
                return true;
            }
        }
        return false;
    }

    /** V = 0, so the only force is the configuration-dependent inertia:
     *  {@code gddot = -(1/2)(m'/m) gdot^2}. FreeDynamics' equation, per cell. */
    private static void rk4(double[] g, double[] gd, double h) {
        int n = g.length;
        double[] k1x = new double[n], k1v = new double[n];
        double[] k2x = new double[n], k2v = new double[n];
        double[] k3x = new double[n], k3v = new double[n];
        double[] k4x = new double[n], k4v = new double[n];
        deriv(g, gd, k1x, k1v);
        double[] t = new double[n], td = new double[n];
        step(g, gd, k1x, k1v, h / 2, t, td);
        deriv(t, td, k2x, k2v);
        step(g, gd, k2x, k2v, h / 2, t, td);
        deriv(t, td, k3x, k3v);
        step(g, gd, k3x, k3v, h, t, td);
        deriv(t, td, k4x, k4v);
        for (int i = 0; i < n; i++) {
            g[i] += h / 6 * (k1x[i] + 2 * k2x[i] + 2 * k3x[i] + k4x[i]);
            gd[i] += h / 6 * (k1v[i] + 2 * k2v[i] + 2 * k3v[i] + k4v[i]);
        }
    }

    private static void deriv(double[] g, double[] gd, double[] dx, double[] dv) {
        for (int i = 0; i < g.length; i++) {
            dx[i] = gd[i];
            dv[i] = -0.5 * (inertiaD(g[i]) / inertia(g[i])) * gd[i] * gd[i];
        }
    }

    private static void step(double[] g, double[] gd, double[] dx, double[] dv,
                             double h, double[] outG, double[] outGd) {
        for (int i = 0; i < g.length; i++) {
            outG[i] = g[i] + h * dx[i];
            outGd[i] = gd[i] + h * dv[i];
        }
    }
}
