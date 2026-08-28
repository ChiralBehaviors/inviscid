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

package com.chiralbehaviors.inviscid.jitterbug;

import java.util.ArrayList;
import java.util.List;

/**
 * An array of jitterbug cells in <b>reduced coordinates</b>, where the vertex
 * ellipses cannot be left because no coordinate could leave them.
 *
 * <p>
 * <b>What this replaces.</b> {@link JitterbugLinkage} is the cell as a
 * bar-and-hinge linkage: eight free triangles, twelve shared vertices, six
 * internal degrees of freedom. That is the right object for asking what the
 * linkage <em>is</em>, and it is the wrong object to integrate an array in. Its
 * constraint set holds strut lengths and nothing else, which is strictly weaker
 * than the jitterbug and admits <b>five spurious degrees of freedom per
 * cell</b>. Every array integration in this project before this class ran in
 * that space.
 *
 * <p>
 * <b>The coordinates.</b> Cell {@code k} carries a centre {@code c_k} (3), an
 * orientation {@code R_k} (3), and a fold angle {@code g_k} (1) — seven
 * generalized coordinates — and its vertex instance {@code i} sits at
 *
 * <pre>
 *     x_ki = c_k + R_k . v_i(g_k)
 * </pre>
 *
 * with {@code v} the body-frame jitterbug vertex, which is
 * {@link JitterbugGeometry#corners(double)} collapsed onto the twelve vertex
 * identities. Two things then hold <i>by construction</i>:
 *
 * <ul>
 * <li><b>All 24 strut lengths.</b> The struts are the edges of the eight rigid
 * triangles and {@code v(g)} turns each triangle rigidly about its own axis.</li>
 * <li><b>The ellipses.</b> {@code v_i} is a function of {@code g} alone, so
 * nothing in the state can move a vertex off its ellipse. The ellipse residual
 * is not small here, it is <em>absent</em>.</li>
 * </ul>
 *
 * So the only constraints are the <b>welds</b> — nine scalar rows per shared
 * face, of which exactly six are independent, because both mating triangles are
 * equilateral of edge {@code L_EDGE} at every {@code g} and their shapes
 * therefore always agree. The multiplier solve is least squares for that
 * reason, not out of caution.
 *
 * <p>
 * <b>What the count then says.</b> For a tree of welds the array has
 * {@code 7N - 6(N-1) = N + 6} dimensions, so <b>N internal degrees of freedom
 * for N cells</b> — one fold angle each. Freeze the orientations and the same
 * Jacobian, restricted to its centre and fold columns, has rank 4 per weld
 * instead of 6 and leaves <b>one</b> coordinate at any size. That contrast is
 * why the medium is a field: with orientations fixed a disturbance has nowhere
 * to be. {@code ReducedCoordinatesTest} measures both columns.
 *
 * <p>
 * <b>Congruence is not identity</b> (Fuller 1977). A cell carries 24 corner
 * instances at mass {@code 1/3}, hence {@code 2/3} on each of its
 * <em>twelve</em> vertex identities and eight units of mass at every
 * configuration. At the octahedron those twelve identities occupy six positions;
 * deduplicating by position would leave four, which is half the cell's inertia.
 * {@link #VERTEX_MASS} is counted from the slot incidence rather than written
 * down, so the accounting cannot be lost to a constant.
 *
 * <p>
 * <b>Mass model, declared:</b> unit mass per triangle, lumped {@code m/3} to
 * each corner. It matches the Python harness's {@code jb_ic} and {@code jb_rc}
 * so the two are comparable; the uniform-lamina alternative moves a period by 7%
 * elsewhere in this project, so no number here is model-independent.
 *
 * <p>
 * <b>V = 0.</b> There is no potential energy anywhere in this model — the owner
 * decision of 2026-08-27 rejects strut compliance, and finding {@code V} for the
 * rigid linkage is still the open problem. Two consequences worth stating before
 * anyone reads a wave out of an animation: energy is the <em>only</em> audit
 * there is, and the fold angles <b>wind without bound</b> rather than
 * oscillating, because nothing restores a cell's phase. What redistributes is
 * energy, not a bounded displacement.
 *
 * <p>
 * <b>No new dependencies: plain arithmetic only.</b> commons-math3 is test
 * scope in this project and stays there — see {@code Linear}'s note on why. The
 * multiplier solve is part of how the dynamics is <em>defined</em>, so it lives
 * here as a symmetric Jacobi eigendecomposition; rank and nullity are how the
 * array is <em>measured</em>, so they live in the test.
 *
 * @author halhildebrand
 */
public final class ReducedCoordinates {

    /**
     * One shared face: cell {@code a}'s vertices {@code va[i]} coincide with cell
     * {@code b}'s vertices {@code vb[i]}, for the three corners of the face.
     */
    public record Weld(int a, int b, int[] va, int[] vb) {
    }

    /**
     * A welded array of cells, and its equations of motion.
     *
     * <p>
     * State is packed as {@code q} then {@code u}: per cell
     * {@code q = (cx, cy, cz, qw, qx, qy, qz, gammaDeg)} and
     * {@code u = (cdot, omega, gammadot)} with {@code omega} the <em>spatial</em>
     * angular velocity and {@code gammadot} in <b>radians</b> per unit time. The
     * fold angle is carried in degrees because the body function is written in
     * degrees; the conversion happens in {@link #derivatives} and nowhere else.
     */
    public static final class Assembly {
        private final double[] centres;   // n x 3, initial
        private final double[] gammas;    // n, initial, degrees
        private final int      n;
        private final Weld[]   welds;

        Assembly(double[] gammas, double[][] centres, Weld[] welds) {
            this.n = gammas.length;
            this.gammas = gammas.clone();
            this.centres = new double[3 * n];
            for (int k = 0; k < n; k++) {
                System.arraycopy(centres[k], 0, this.centres, 3 * k, 3);
            }
            this.welds = welds.clone();
        }

        /**
         * Generalized accelerations from the acceleration-level equations. V = 0,
         * so the only forces are the weld reactions.
         */
        public double[] accel(double[] q, double[] u) {
            double[][][] jac = new double[n][][];
            double[][][] mass = new double[n][][];
            double[][] gyro = new double[n][];
            double[] aFree = new double[NU * n];
            for (int k = 0; k < n; k++) {
                double[][] rot = rotationOf(q, k);
                double[][][] v = body(q[NQ * k + 7], 2);
                jac[k] = cellJacobian(rot, v);
                mass[k] = massBlock(jac[k]);
                gyro[k] = cellGyro(rot, v, u, k);
                double[] f = new double[NU];
                for (int r = 0; r < 3 * NV; r++) {
                    double m = VERTEX_MASS[r / 3] * gyro[k][r];
                    for (int c = 0; c < NU; c++) {
                        f[c] -= jac[k][r][c] * m;
                    }
                }
                double[] a = solveSpd(mass[k], f, 1e-12);
                System.arraycopy(a, 0, aFree, NU * k, NU);
            }
            if (welds.length == 0) {
                return aFree;
            }
            double[][] c = constraintJacobian(jac);
            double[] d = new double[c.length];
            for (int w = 0; w < welds.length; w++) {
                Weld wd = welds[w];
                for (int i = 0; i < 3; i++) {
                    for (int t = 0; t < 3; t++) {
                        d[9 * w + 3 * i + t] = -(gyro[wd.a][3 * wd.va[i] + t]
                                                 - gyro[wd.b][3 * wd.vb[i] + t]);
                    }
                }
            }
            return constrain(mass, c, aFree, d);
        }

        public int cells() {
            return n;
        }

        /**
         * The constraint Jacobian, {@code 9W x 7N}: the weld rows and nothing
         * else. There is no ellipse row and no bar row, and that absence is the
         * whole point of the coordinate change.
         */
        public double[][] constraintJacobian(double[] q) {
            double[][][] jac = new double[n][][];
            for (int k = 0; k < n; k++) {
                jac[k] = cellJacobian(rotationOf(q, k), body(q[NQ * k + 7], 1));
            }
            return constraintJacobian(jac);
        }

        /** Kinetic energy of each cell. Cells own their own vertex instances, so
         *  there is no shared-vertex ownership to divide and no double count. */
        public double[] energyPerCell(double[] q, double[] u) {
            double[] e = new double[n];
            for (int k = 0; k < n; k++) {
                double[][] m = massBlock(cellJacobian(rotationOf(q, k),
                                                      body(q[NQ * k + 7], 1)));
                double s = 0;
                for (int i = 0; i < NU; i++) {
                    for (int j = 0; j < NU; j++) {
                        s += u[NU * k + i] * m[i][j] * u[NU * k + j];
                    }
                }
                e[k] = 0.5 * s;
            }
            return e;
        }

        /**
         * The local disturbance the field framing asks for: {@code gammadot = 1}
         * on one cell, everything else zero, projected onto the admissible set.
         *
         * <p>
         * A cell has no other internal freedom, so this is the only local kick
         * there is. The spin-on-one-triangle impulse used before this coordinate
         * change was never on the jitterbug path at all, which is why its
         * projection had to mangle it.
         */
        public double[] foldImpulse(int cell) {
            double[] u = new double[NU * n];
            u[NU * cell + 6] = 1.0;
            return projectVelocity(initialState(), u, true);
        }

        /**
         * The six global rigid motions. Translations are {@code dc = t};
         * rotations are {@code dc = omega x c}, {@code omega_k = omega}.
         *
         * <p>
         * The {@code omega x c} term is load bearing. An earlier fold-map-rank
         * run in this project built this basis with an identity orientation block
         * and no such term, projected out a six-dimensional space that is not the
         * global motions, and reported a degree-of-freedom column two too high.
         *
         * @param defective reproduce that basis instead, so a test can show the
         *                  correct one is annihilated by the constraint Jacobian
         *                  and the wrong one is not
         */
        public double[][] globalMotions(double[] q, boolean defective) {
            double[][] g = new double[6][NU * n];
            for (int d = 0; d < 3; d++) {
                for (int k = 0; k < n; k++) {
                    g[d][NU * k + d] = 1.0;
                }
            }
            for (int d = 0; d < 3; d++) {
                double[] w = new double[3];
                w[d] = 1.0;
                for (int k = 0; k < n; k++) {
                    if (!defective) {
                        double[] cr = cross(w, new double[] { q[NQ * k], q[NQ * k + 1],
                                                              q[NQ * k + 2] });
                        System.arraycopy(cr, 0, g[3 + d], NU * k, 3);
                    }
                    g[3 + d][NU * k + 3 + d] = 1.0;
                }
            }
            return g;
        }

        /** The initial configuration: given centres and fold angles, identity
         *  orientations. */
        public double[] initialState() {
            double[] q = new double[NQ * n];
            for (int k = 0; k < n; k++) {
                System.arraycopy(centres, 3 * k, q, NQ * k, 3);
                q[NQ * k + 3] = 1.0;
                q[NQ * k + 7] = gammas[k];
            }
            return q;
        }

        /**
         * Total linear and angular momentum as linear functionals of {@code u}.
         * Used to keep an impulse internal: a fold kick on one cell does induce
         * net momentum, and without removing it part of what follows would be the
         * whole patch drifting rather than transport.
         */
        public double[][] momentumRows(double[] q) {
            double[][] p = new double[6][NU * n];
            for (int k = 0; k < n; k++) {
                double[][] rot = rotationOf(q, k);
                double[][][] v = body(q[NQ * k + 7], 1);
                double[][] jac = cellJacobian(rot, v);
                for (int i = 0; i < NV; i++) {
                    double[] x = new double[3];
                    for (int t = 0; t < 3; t++) {
                        x[t] = q[NQ * k + t] + rot[t][0] * v[0][i][0]
                               + rot[t][1] * v[0][i][1] + rot[t][2] * v[0][i][2];
                    }
                    double m = VERTEX_MASS[i];
                    for (int c = 0; c < NU; c++) {
                        double[] col = { jac[3 * i][c], jac[3 * i + 1][c],
                                         jac[3 * i + 2][c] };
                        for (int t = 0; t < 3; t++) {
                            p[t][NU * k + c] += m * col[t];
                        }
                        double[] cr = cross(x, col);
                        for (int t = 0; t < 3; t++) {
                            p[3 + t][NU * k + c] += m * cr[t];
                        }
                    }
                }
            }
            return p;
        }

        /** World positions of every cell's twelve vertex identities. */
        public double[][][] positions(double[] q) {
            double[][][] x = new double[n][NV][3];
            for (int k = 0; k < n; k++) {
                double[][] rot = rotationOf(q, k);
                double[][] v = body(q[NQ * k + 7], 0)[0];
                for (int i = 0; i < NV; i++) {
                    for (int t = 0; t < 3; t++) {
                        x[k][i][t] = q[NQ * k + t] + rot[t][0] * v[i][0]
                                     + rot[t][1] * v[i][1] + rot[t][2] * v[i][2];
                    }
                }
            }
            return x;
        }

        /**
         * The nearest admissible velocity in the kinetic-energy metric:
         * minimises {@code (u-u0)' M (u-u0)} subject to {@code C u = 0} and,
         * optionally, to vanishing total linear and angular momentum.
         */
        public double[] projectVelocity(double[] q, double[] u0, boolean momentum) {
            double[][][] jac = new double[n][][];
            double[][][] mass = new double[n][][];
            for (int k = 0; k < n; k++) {
                jac[k] = cellJacobian(rotationOf(q, k), body(q[NQ * k + 7], 1));
                mass[k] = massBlock(jac[k]);
            }
            List<double[]> rows = new ArrayList<>();
            for (double[] r : constraintJacobian(jac)) {
                rows.add(r);
            }
            if (momentum) {
                for (double[] r : momentumRows(q)) {
                    rows.add(r);
                }
            }
            double[][] c = rows.toArray(new double[0][]);
            // The target is C u = 0; constrain() subtracts C.u0 itself.
            return constrain(mass, c, u0, new double[c.length]);
        }

        /**
         * Advance {@code state} (the {@code q} block followed by the {@code u}
         * block) by {@code dt}, in {@code substeps} classical Runge-Kutta steps,
         * then put the configuration and the velocity back onto the weld manifold.
         *
         * <p>
         * The projection is what keeps a fixed-step integrator honest over a long
         * run: without it the acceleration-level formulation lets the weld
         * residual drift quadratically in time. It is reported rather than hidden
         * — {@link #weldResidual} is the audit, and the test reads it unprojected
         * within a step.
         */
        public void step(double[] state, double dt, int substeps) {
            double h = dt / substeps;
            for (int s = 0; s < substeps; s++) {
                double[] k1 = derivatives(state);
                double[] k2 = derivatives(axpy(state, k1, 0.5 * h));
                double[] k3 = derivatives(axpy(state, k2, 0.5 * h));
                double[] k4 = derivatives(axpy(state, k3, h));
                for (int i = 0; i < state.length; i++) {
                    state[i] += h / 6.0
                                * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
                }
                normalize(state);
                project(state);
            }
        }

        /** Worst absolute weld closure error over every shared corner. */
        public double weldResidual(double[] q) {
            double[] g = weldGap(q);
            double w = 0;
            for (double x : g) {
                w = Math.max(w, Math.abs(x));
            }
            return w;
        }

        double[] weldGap(double[] q) {
            double[][][] x = positions(q);
            double[] g = new double[9 * welds.length];
            for (int w = 0; w < welds.length; w++) {
                Weld wd = welds[w];
                for (int i = 0; i < 3; i++) {
                    for (int t = 0; t < 3; t++) {
                        g[9 * w + 3 * i + t] = x[wd.a][wd.va[i]][t]
                                               - x[wd.b][wd.vb[i]][t];
                    }
                }
            }
            return g;
        }

        Weld[] welds() {
            return welds.clone();
        }

        /**
         * The velocity-only part of the ambient acceleration, {@code Jdot . u}:
         *
         * <pre>
         *     w x (w x Rv) + 2 w x (R v' gdot) + R v'' gdot^2
         * </pre>
         */
        private double[] cellGyro(double[][] rot, double[][][] v, double[] u, int k) {
            double[] w = { u[NU * k + 3], u[NU * k + 4], u[NU * k + 5] };
            double gd = u[NU * k + 6];
            double[] a = new double[3 * NV];
            for (int i = 0; i < NV; i++) {
                double[] r0 = apply(rot, v[0][i]);
                double[] r1 = apply(rot, v[1][i]);
                double[] r2 = apply(rot, v[2][i]);
                double[] t1 = cross(w, cross(w, r0));
                double[] t2 = cross(w, new double[] { r1[0] * gd, r1[1] * gd,
                                                      r1[2] * gd });
                for (int t = 0; t < 3; t++) {
                    a[3 * i + t] = t1[t] + 2 * t2[t] + r2[t] * gd * gd;
                }
            }
            return a;
        }

        /** Solve the KKT system in least squares and add the reaction. */
        private double[] constrain(double[][][] mass, double[][] c, double[] free,
                                   double[] rhs) {
            int rowsN = c.length;
            double[][] miCt = new double[NU * n][rowsN];
            for (int k = 0; k < n; k++) {
                double[][] mi = pinvSpd(mass[k], 1e-12);
                for (int r = 0; r < rowsN; r++) {
                    for (int i = 0; i < NU; i++) {
                        double v = 0;
                        for (int j = 0; j < NU; j++) {
                            v += mi[i][j] * c[r][NU * k + j];
                        }
                        miCt[NU * k + i][r] = v;
                    }
                }
            }
            double[][] a = new double[rowsN][rowsN];
            for (int r = 0; r < rowsN; r++) {
                for (int s = 0; s < rowsN; s++) {
                    double v = 0;
                    for (int i = 0; i < NU * n; i++) {
                        v += c[r][i] * miCt[i][s];
                    }
                    a[r][s] = v;
                }
            }
            double[] b = new double[rowsN];
            for (int r = 0; r < rowsN; r++) {
                b[r] = rhs[r] - dot(c[r], free);
            }
            double[] lam = solveSpd(a, b, 1e-10);
            double[] out = free.clone();
            for (int i = 0; i < NU * n; i++) {
                for (int r = 0; r < rowsN; r++) {
                    out[i] += miCt[i][r] * lam[r];
                }
            }
            return out;
        }

        private double[][] constraintJacobian(double[][][] jac) {
            double[][] c = new double[9 * welds.length][NU * n];
            for (int w = 0; w < welds.length; w++) {
                Weld wd = welds[w];
                for (int i = 0; i < 3; i++) {
                    for (int t = 0; t < 3; t++) {
                        int row = 9 * w + 3 * i + t;
                        for (int col = 0; col < NU; col++) {
                            c[row][NU * wd.a + col] = jac[wd.a][3 * wd.va[i] + t][col];
                            c[row][NU * wd.b + col] -= jac[wd.b][3 * wd.vb[i] + t][col];
                        }
                    }
                }
            }
            return c;
        }

        /** State derivative: qdot from u, udot from the equations of motion. */
        private double[] derivatives(double[] state) {
            double[] q = new double[NQ * n];
            double[] u = new double[NU * n];
            System.arraycopy(state, 0, q, 0, NQ * n);
            System.arraycopy(state, NQ * n, u, 0, NU * n);
            double[] out = new double[state.length];
            for (int k = 0; k < n; k++) {
                out[NQ * k] = u[NU * k];
                out[NQ * k + 1] = u[NU * k + 1];
                out[NQ * k + 2] = u[NU * k + 2];
                double[] qq = { q[NQ * k + 3], q[NQ * k + 4], q[NQ * k + 5],
                                q[NQ * k + 6] };
                double nq = Math.sqrt(dot(qq, qq));
                for (int t = 0; t < 4; t++) {
                    qq[t] /= nq;
                }
                double[] w = { u[NU * k + 3], u[NU * k + 4], u[NU * k + 5] };
                double[] qv = { qq[1], qq[2], qq[3] };
                double[] cr = cross(w, qv);
                out[NQ * k + 3] = -0.5 * dot(w, qv);
                for (int t = 0; t < 3; t++) {
                    out[NQ * k + 4 + t] = 0.5 * (qq[0] * w[t] + cr[t]);
                }
                out[NQ * k + 7] = Math.toDegrees(u[NU * k + 6]);
            }
            double[] a = accel(q, u);
            System.arraycopy(a, 0, out, NQ * n, NU * n);
            return out;
        }

        /**
         * Gauss-Newton back onto {@code g = 0}, then remove any velocity that has
         * left the tangent space.
         */
        private void project(double[] state) {
            if (welds.length == 0) {
                return;
            }
            for (int it = 0; it < 3; it++) {
                double[] q = new double[NQ * n];
                System.arraycopy(state, 0, q, 0, NQ * n);
                double[] g = weldGap(q);
                double worst = 0;
                for (double x : g) {
                    worst = Math.max(worst, Math.abs(x));
                }
                if (worst < 1e-14) {
                    break;
                }
                double[][] c = constraintJacobian(q);
                double[][] cct = new double[c.length][c.length];
                for (int r = 0; r < c.length; r++) {
                    for (int s = 0; s < c.length; s++) {
                        cct[r][s] = dot(c[r], c[s]);
                    }
                }
                double[] y = solveSpd(cct, g, 1e-10);
                double[] delta = new double[NU * n];
                for (int i = 0; i < NU * n; i++) {
                    double v = 0;
                    for (int r = 0; r < c.length; r++) {
                        v += c[r][i] * y[r];
                    }
                    delta[i] = -v;
                }
                applyIncrement(state, delta);
            }
            double[] q = new double[NQ * n];
            double[] u = new double[NU * n];
            System.arraycopy(state, 0, q, 0, NQ * n);
            System.arraycopy(state, NQ * n, u, 0, NU * n);
            double[] pu = projectVelocity(q, u, false);
            System.arraycopy(pu, 0, state, NQ * n, NU * n);
        }

        /** Apply a generalized displacement, rotating the quaternion by the
         *  exponential map so the orientation stays on SO(3). */
        private void applyIncrement(double[] state, double[] delta) {
            for (int k = 0; k < n; k++) {
                for (int t = 0; t < 3; t++) {
                    state[NQ * k + t] += delta[NU * k + t];
                }
                double[] w = { delta[NU * k + 3], delta[NU * k + 4],
                               delta[NU * k + 5] };
                double th = Math.sqrt(dot(w, w));
                double[] dq = { 1, 0, 0, 0 };
                if (th > 0) {
                    double s = Math.sin(0.5 * th) / th;
                    dq = new double[] { Math.cos(0.5 * th), s * w[0], s * w[1],
                                        s * w[2] };
                }
                double[] qq = { state[NQ * k + 3], state[NQ * k + 4],
                                state[NQ * k + 5], state[NQ * k + 6] };
                double[] pr = quatMultiply(dq, qq);
                System.arraycopy(pr, 0, state, NQ * k + 3, 4);
                state[NQ * k + 7] += Math.toDegrees(delta[NU * k + 6]);
            }
            normalize(state);
        }

        private void normalize(double[] state) {
            for (int k = 0; k < n; k++) {
                double s = 0;
                for (int t = 0; t < 4; t++) {
                    s += state[NQ * k + 3 + t] * state[NQ * k + 3 + t];
                }
                s = Math.sqrt(s);
                for (int t = 0; t < 4; t++) {
                    state[NQ * k + 3 + t] /= s;
                }
            }
        }

        private double[][] rotationOf(double[] q, int k) {
            return quatToMatrix(q[NQ * k + 3], q[NQ * k + 4], q[NQ * k + 5],
                                q[NQ * k + 6]);
        }
    }

    /** Where the slot-to-vertex map is read. Any angle whose configuration is
     *  12 vertices of multiplicity 2 will do; 30 degrees is generic. Read at the
     *  octahedron it would see six vertices and invent a different mechanism. */
    public static final double GENERIC_DEG = 30.0;

    /** Coordinates in the packed configuration block, per cell. */
    public static final int NQ = 8;

    /** Vertex identities per cell. Twelve at <em>every</em> fold angle. */
    public static final int NV = 12;

    /** Generalized velocities per cell: centre, spatial angular velocity, fold. */
    public static final int NU = 7;

    /** The eight {@code <111>} directions, in the face order of
     *  {@link JitterbugGeometry}. */
    public static final int[][] DIAGONALS;

    /** {@code [face][corner] -> vertex identity 0..11}, read at
     *  {@link #GENERIC_DEG}. */
    public static final int[][] SLOT;

    /** Mass on each vertex identity: the number of corner instances that land on
     *  it, times {@code 1/3}. Counted, not written down. */
    public static final double[] VERTEX_MASS;

    private static final double[][] FACE_DIR;

    static {
        double[][][] generic = JitterbugGeometry.corners(GENERIC_DEG);
        JitterbugGeometry.Clustering cl = JitterbugGeometry.cluster(generic, 1e-9);
        SLOT = cl.labels();
        VERTEX_MASS = new double[NV];
        for (int f = 0; f < 8; f++) {
            for (int c = 0; c < 3; c++) {
                VERTEX_MASS[SLOT[f][c]] += 1.0 / 3.0;
            }
        }
        double[][] cen = JitterbugGeometry.faceCentroids(generic);
        FACE_DIR = new double[8][3];
        DIAGONALS = new int[8][3];
        for (int f = 0; f < 8; f++) {
            double nn = Math.sqrt(cen[f][0] * cen[f][0] + cen[f][1] * cen[f][1]
                                  + cen[f][2] * cen[f][2]);
            for (int t = 0; t < 3; t++) {
                FACE_DIR[f][t] = cen[f][t] / nn;
                DIAGONALS[f][t] = FACE_DIR[f][t] > 0 ? 1 : -1;
            }
        }
    }

    /**
     * The body-frame vertices {@code v(g)} and, on request, {@code dv/dg} and
     * {@code d2v/dg2} in <b>radians</b>.
     *
     * <p>
     * With {@code x(a) = R(u, sigma(a-60)) (v-c) + u Z cos a} and {@code K} the
     * cross-product matrix of the face axis,
     *
     * <pre>
     *     dx/da   = sigma   K   R (v-c) - u Z sin(a)
     *     d2x/da2 = sigma^2 K^2 R (v-c) - u Z cos(a)
     * </pre>
     *
     * The per-face {@code sigma} is invisible in the positions and wrong in every
     * velocity if it is dropped, which is why the test checks these against
     * central differences rather than reading them.
     *
     * @return {@code [nder+1][12][3]}
     */
    public static double[][][] body(double gammaDeg, int nder) {
        double[][][] out = new double[nder + 1][NV][3];
        boolean[] seen = new boolean[NV];
        double a = Math.toRadians(gammaDeg);
        double[][][] x = JitterbugGeometry.corners(gammaDeg);
        for (int f = 0; f < 8; f++) {
            double[] u = JitterbugGeometry.faceAxis(f);
            double sigma = u[0] * u[1] * u[2] > 0 ? 1.0 : -1.0;
            for (int c = 0; c < 3; c++) {
                int i = SLOT[f][c];
                if (seen[i]) {
                    continue;
                }
                seen[i] = true;
                System.arraycopy(x[f][c], 0, out[0][i], 0, 3);
                if (nder == 0) {
                    continue;
                }
                // base = R (v-c) = x - u Z cos a, recovered rather than rebuilt
                // so the position and its derivatives can never disagree about
                // which parameterisation they belong to.
                double[] base = new double[3];
                for (int t = 0; t < 3; t++) {
                    base[t] = x[f][c][t]
                              - u[t] * JitterbugGeometry.Z * Math.cos(a);
                }
                double[] kb = cross(u, base);
                for (int t = 0; t < 3; t++) {
                    out[1][i][t] = sigma * kb[t]
                                   - u[t] * JitterbugGeometry.Z * Math.sin(a);
                }
                if (nder >= 2) {
                    double[] kkb = cross(u, kb);
                    for (int t = 0; t < 3; t++) {
                        out[2][i][t] = kkb[t]
                                       - u[t] * JitterbugGeometry.Z * Math.cos(a);
                    }
                }
            }
        }
        return out;
    }

    /**
     * A straight line of {@code n} cells along the {@code (1,1,1)} diagonal,
     * fold angles alternating {@code +-gamma0}, separations from the shared-face
     * law {@code sep = Z (cos a + cos b)}.
     */
    public static Assembly chain(int n, double gamma0) {
        double[] gam = new double[n];
        for (int k = 0; k < n; k++) {
            gam[k] = k % 2 == 0 ? gamma0 : -gamma0;
        }
        double[] axis = FACE_DIR[faceAlong(new double[] { 1, 1, 1 })];
        double zc = JitterbugGeometry.L_EDGE * Math.sqrt(2.0 / 3.0);
        double[][] ctr = new double[n][3];
        for (int k = 1; k < n; k++) {
            double sep = zc * (Math.cos(Math.toRadians(gam[k - 1]))
                               + Math.cos(Math.toRadians(gam[k])));
            for (int t = 0; t < 3; t++) {
                ctr[k][t] = ctr[k - 1][t] + sep * axis[t];
            }
        }
        Weld[] welds = new Weld[n - 1];
        for (int k = 0; k < n - 1; k++) {
            welds[k] = weld(k, k + 1, gam[k], gam[k + 1], ctr[k], ctr[k + 1]);
        }
        return new Assembly(gam, ctr, welds);
    }

    /**
     * A centre cell at {@code gammaCentre} with a neighbour at
     * {@code gammaCentre + 60} down each of the given {@code <111>} diagonals.
     *
     * <p>
     * All eight is the smallest patch with an <em>interior</em>: a line never has
     * one, because only 2 of each cell's 8 faces are shared along it, so every
     * vertex stays on an outward-facing face however long the line gets. Two
     * adjacent diagonals give the three-cell V, whose outer cells sit at
     * {@code n1 . n2 = 1/3}.
     */
    public static Assembly cluster(double gammaCentre, int[][] sites) {
        double gn = gammaCentre + 60.0;
        double zc = JitterbugGeometry.L_EDGE * Math.sqrt(2.0 / 3.0);
        double sep = zc * (Math.cos(Math.toRadians(gammaCentre))
                           + Math.cos(Math.toRadians(gn)));
        double side = sep / Math.sqrt(3.0);
        int n = sites.length + 1;
        double[] gam = new double[n];
        double[][] ctr = new double[n][3];
        gam[0] = gammaCentre;
        for (int m = 0; m < sites.length; m++) {
            gam[m + 1] = gn;
            for (int t = 0; t < 3; t++) {
                ctr[m + 1][t] = side * sites[m][t];
            }
        }
        Weld[] welds = new Weld[sites.length];
        for (int m = 0; m < sites.length; m++) {
            welds[m] = weld(0, m + 1, gammaCentre, gn, ctr[0], ctr[m + 1]);
        }
        return new Assembly(gam, ctr, welds);
    }

    /** The index of the face pointing most nearly along {@code dir}. */
    public static int faceAlong(double[] dir) {
        double nn = Math.sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
        int best = 0;
        double bv = -2;
        for (int f = 0; f < 8; f++) {
            double d = (FACE_DIR[f][0] * dir[0] + FACE_DIR[f][1] * dir[1]
                        + FACE_DIR[f][2] * dir[2]) / nn;
            if (d > bv) {
                bv = d;
                best = f;
            }
        }
        return best;
    }

    /**
     * The corner correspondence across the face shared by two cells.
     *
     * <p>
     * The search is restricted to the <b>mating face's three corners</b>. Letting
     * it range over all twelve of the neighbour's vertices — which is what the
     * Python harness did first — lets a nearer non-corner win, and the resulting
     * weld joins three vertices that form no triangle. Nothing catches that by
     * looking at positions: the doubly-written shared vertices still agree to
     * 3e-15, because a redundancy check on POSITIONS does not validate IDENTITY.
     * The assertion at the end is on the identities.
     */
    public static Weld weld(int a, int b, double ga, double gb, double[] ca,
                            double[] cb) {
        double[] d = { cb[0] - ca[0], cb[1] - ca[1], cb[2] - ca[2] };
        int fa = faceAlong(d);
        int fb = faceAlong(new double[] { -d[0], -d[1], -d[2] });
        double[][] va = body(ga, 0)[0];
        double[][] vb = body(gb, 0)[0];
        int[] target = { SLOT[fb][0], SLOT[fb][1], SLOT[fb][2] };
        int[] outA = new int[3];
        int[] outB = new int[3];
        boolean[] used = new boolean[3];
        for (int c = 0; c < 3; c++) {
            int ia = SLOT[fa][c];
            int pick = -1;
            double best = Double.MAX_VALUE;
            for (int j = 0; j < 3; j++) {
                if (used[j]) {
                    continue;
                }
                double s = 0;
                for (int t = 0; t < 3; t++) {
                    double e = (ca[t] + va[ia][t]) - (cb[t] + vb[target[j]][t]);
                    s += e * e;
                }
                if (s < best) {
                    best = s;
                    pick = j;
                }
            }
            used[pick] = true;
            outA[c] = ia;
            outB[c] = target[pick];
        }
        return new Weld(a, b, outA, outB);
    }

    /**
     * Solve {@code A x = b} for symmetric positive semi-definite {@code A}, in
     * least squares, by cyclic Jacobi eigendecomposition. Directions whose
     * eigenvalue falls below {@code relTol} times the largest are dropped, which
     * is what makes the redundant weld rows harmless: nine rows per face carry
     * six independent constraints, and the other three are exactly such
     * directions.
     */
    /**
     * The Moore-Penrose pseudo-inverse of a symmetric positive semi-definite
     * matrix, by cyclic Jacobi eigendecomposition. Directions whose eigenvalue
     * falls below {@code relTol} times the largest are dropped, which is what
     * makes the redundant weld rows harmless: nine rows per shared face carry six
     * independent constraints and the other three are exactly such directions.
     *
     * <p>
     * Jacobi rather than a library SVD because commons-math3 is test scope in
     * this project and the multiplier solve is part of how the dynamics is
     * <em>defined</em>. It is decomposed once per matrix and reused, because
     * doing it per right-hand side is what made the first version of this class
     * unusably slow.
     */
    public static double[][] pinvSpd(double[][] a, double relTol) {
        int m = a.length;
        double[][] s = new double[m][m];
        double[][] v = new double[m][m];
        double scale = 0;
        for (int i = 0; i < m; i++) {
            s[i] = a[i].clone();
            v[i][i] = 1.0;
            scale += a[i][i] * a[i][i];
        }
        double tol = 1e-26 * Math.max(scale, 1e-300);
        for (int sweep = 0; sweep < 40; sweep++) {
            double off = 0;
            for (int p = 0; p < m; p++) {
                for (int q = p + 1; q < m; q++) {
                    off += s[p][q] * s[p][q];
                }
            }
            if (off <= tol) {
                break;
            }
            for (int p = 0; p < m; p++) {
                for (int q = p + 1; q < m; q++) {
                    if (s[p][q] == 0) {
                        continue;
                    }
                    double theta = 0.5 * (s[q][q] - s[p][p]) / s[p][q];
                    double t = theta == 0 ? 1.0
                                          : Math.signum(theta)
                                            / (Math.abs(theta)
                                               + Math.sqrt(theta * theta + 1.0));
                    double cs = 1.0 / Math.sqrt(t * t + 1.0);
                    double sn = t * cs;
                    for (int k = 0; k < m; k++) {
                        double skp = s[k][p], skq = s[k][q];
                        s[k][p] = cs * skp - sn * skq;
                        s[k][q] = sn * skp + cs * skq;
                    }
                    for (int k = 0; k < m; k++) {
                        double spk = s[p][k], sqk = s[q][k];
                        s[p][k] = cs * spk - sn * sqk;
                        s[q][k] = sn * spk + cs * sqk;
                    }
                    for (int k = 0; k < m; k++) {
                        double vkp = v[k][p], vkq = v[k][q];
                        v[k][p] = cs * vkp - sn * vkq;
                        v[k][q] = sn * vkp + cs * vkq;
                    }
                }
            }
        }
        double max = 0;
        for (int i = 0; i < m; i++) {
            max = Math.max(max, Math.abs(s[i][i]));
        }
        double cut = relTol * max;
        double[][] out = new double[m][m];
        for (int i = 0; i < m; i++) {
            if (Math.abs(s[i][i]) <= cut) {
                continue;
            }
            double inv = 1.0 / s[i][i];
            for (int r = 0; r < m; r++) {
                if (v[r][i] == 0) {
                    continue;
                }
                double w = v[r][i] * inv;
                for (int c = 0; c < m; c++) {
                    out[r][c] += w * v[c][i];
                }
            }
        }
        return out;
    }

    /** Least-squares solve of {@code A x = b} for symmetric PSD {@code A}. */
    public static double[] solveSpd(double[][] a, double[] b, double relTol) {
        double[][] pinv = pinvSpd(a, relTol);
        double[] x = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            x[i] = dot(pinv[i], b);
        }
        return x;
    }

    static double[] apply(double[][] m, double[] v) {
        return new double[] {
                              m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
                              m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
                              m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2] };
    }

    static double[] axpy(double[] a, double[] b, double s) {
        double[] out = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            out[i] = a[i] + s * b[i];
        }
        return out;
    }

    /**
     * {@code J_k}, the {@code 36 x 7} map from {@code (cdot, omega, gdot)} to the
     * cell's vertex velocities: {@code I}, {@code -[Rv]x}, and {@code R v'}.
     */
    static double[][] cellJacobian(double[][] rot, double[][][] v) {
        double[][] j = new double[3 * NV][NU];
        for (int i = 0; i < NV; i++) {
            double[] r0 = apply(rot, v[0][i]);
            double[] r1 = apply(rot, v[1][i]);
            for (int t = 0; t < 3; t++) {
                j[3 * i + t][t] = 1.0;
                j[3 * i + t][6] = r1[t];
            }
            j[3 * i][4] = r0[2];
            j[3 * i][5] = -r0[1];
            j[3 * i + 1][3] = -r0[2];
            j[3 * i + 1][5] = r0[0];
            j[3 * i + 2][3] = r0[1];
            j[3 * i + 2][4] = -r0[0];
        }
        return j;
    }

    static double[] cross(double[] a, double[] b) {
        return new double[] { a[1] * b[2] - a[2] * b[1],
                              a[2] * b[0] - a[0] * b[2],
                              a[0] * b[1] - a[1] * b[0] };
    }

    static double dot(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < a.length; i++) {
            s += a[i] * b[i];
        }
        return s;
    }

    /** {@code M_k = J' diag(m) J}. Block diagonal over cells, because a cell's
     *  vertices depend on its own seven coordinates and nothing else. */
    static double[][] massBlock(double[][] j) {
        double[][] m = new double[NU][NU];
        for (int r = 0; r < 3 * NV; r++) {
            double w = VERTEX_MASS[r / 3];
            for (int i = 0; i < NU; i++) {
                if (j[r][i] == 0) {
                    continue;
                }
                for (int k = 0; k < NU; k++) {
                    m[i][k] += w * j[r][i] * j[r][k];
                }
            }
        }
        return m;
    }

    static double[] quatMultiply(double[] a, double[] b) {
        return new double[] {
                              a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
                              a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
                              a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
                              a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0] };
    }

    static double[][] quatToMatrix(double w, double x, double y, double z) {
        double nn = Math.sqrt(w * w + x * x + y * y + z * z);
        w /= nn;
        x /= nn;
        y /= nn;
        z /= nn;
        return new double[][] {
                                { 1 - 2 * (y * y + z * z), 2 * (x * y - w * z),
                                  2 * (x * z + w * y) },
                                { 2 * (x * y + w * z), 1 - 2 * (x * x + z * z),
                                  2 * (y * z - w * x) },
                                { 2 * (x * z - w * y), 2 * (y * z + w * x),
                                  1 - 2 * (x * x + y * y) } };
    }

    private ReducedCoordinates() {
    }
}
