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

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/**
 * Marching from one configuration of the linkage to another while staying on the
 * constraint variety.
 *
 * <p>
 * The Python reference drove {@code scipy.optimize.least_squares} over a global
 * 48-vector of body motions with a finite-difference Jacobian. This does the
 * same march incrementally instead: at each waypoint it linearises about the
 * <em>current</em> configuration, where the analytic Jacobians
 * {@link JitterbugLinkage#jacobian(double[][][])} and
 * {@link JitterbugLinkage#positionJacobian(double[][][])} are exact, solves for a
 * small body-motion increment, and applies it. Different solver, same variety —
 * which is the point of having two implementations.
 *
 * <p>
 * Each waypoint ends with an exact Newton projection back onto
 * {@code residual = 0}, using the minimum-norm increment, so the demonstrated
 * path provably lies <em>on</em> the variety rather than near it.
 *
 * @author halhildebrand
 */
final class VarietyWalk {

    /**
     * @param configuration    where the walk ended up
     * @param rmsToTarget      RMS corner distance to the target after optimal
     *                         rigid alignment; below 1e-6 counts as reached
     * @param worstResidual    largest post-projection hinge residual over all
     *                         waypoints — the evidence the path is on the
     *                         variety
     * @param interpenetrating how many waypoints self-intersected
     * @param waypoints        how many waypoints there were
     */
    record Result(double[][][] configuration, double rmsToTarget,
                  double worstResidual, int interpenetrating, int waypoints) {

        boolean reached() {
            return rmsToTarget < 1e-6;
        }
    }

    /** Weight on the hinge rows relative to the attraction rows. */
    private static final double HINGE_WEIGHT = 1e4;

    private VarietyWalk() {
    }

    /**
     * Predictor-corrector continuation along an internal null direction.
     *
     * <p>
     * At each step the internal null space is recomputed at the <em>current</em>
     * configuration and the branch closest to the previous heading is taken, so
     * this tracks a curve rather than repeatedly re-projecting a stale
     * direction. Returns the accumulated arclength in corner-position space; a
     * first-order flex that does not extend to a finite motion dies almost
     * immediately, while a genuine motion runs the full step budget.
     *
     * @param seed initial 48-vector heading; need not be exactly null, it is
     *             projected onto the internal null space before the first step
     */
    static double continueAlong(double[][][] start, double[] seed, double step,
                                int steps) {
        double[][][] x = start;
        double[] heading = null;
        double arclength = 0;
        for (int s = 0; s < steps; s++) {
            double[][] basis = internalNullSpace(x);
            if (basis.length == 0) {
                break;
            }
            double[] direction = closestTo(basis, heading == null ? seed
                                                                  : heading);
            if (direction == null) {
                break;
            }
            heading = direction;
            double[] increment = new double[48];
            for (int k = 0; k < 48; k++) {
                increment[k] = step * direction[k];
            }
            double[][][] predicted = JitterbugLinkage.applyBodyMotions(x,
                                                                       increment);
            double[][][] corrected = project(predicted);
            if (Linear.norm(JitterbugLinkage.residual(corrected)) > 1e-10) {
                break;
            }
            arclength += frobenius(x, corrected);
            x = corrected;
        }
        return arclength;
    }

    /**
     * An orthonormal basis for the null space of the constraint Jacobian with
     * the six global rigid motions removed — the internal degrees of freedom.
     */
    static double[][] internalNullSpace(double[][][] x) {
        double[][] j = JitterbugLinkage.jacobian(x);
        RealMatrix jm = new Array2DRowRealMatrix(j, true);
        RealMatrix jtj = jm.transpose().multiply(jm);
        EigenDecomposition eig = new EigenDecomposition(jtj);

        List<double[]> raw = new ArrayList<>();
        for (int k = 0; k < 48; k++) {
            if (Math.abs(eig.getRealEigenvalue(k)) < 1e-9) {
                raw.add(eig.getEigenvector(k).toArray());
            }
        }
        double[][] global = JitterbugLinkage.globalRigidMotions(x);
        List<double[]> globalBasis = orthonormalise(List.of(global), 1e-9);

        List<double[]> stripped = new ArrayList<>();
        for (double[] v : raw) {
            double[] w = v.clone();
            for (double[] g : globalBasis) {
                double d = dot(w, g);
                for (int k = 0; k < 48; k++) {
                    w[k] -= d * g[k];
                }
            }
            stripped.add(w);
        }
        List<double[]> basis = orthonormalise(stripped, 1e-6);
        return basis.toArray(new double[0][]);
    }

    /**
     * {@code heading} projected onto the span of {@code basis} and renormalised
     * — the tangent direction the curve should continue in. Returns null when
     * the heading has essentially no component in the tangent space, which is
     * how a dying branch, or a direction transverse to the variety, is
     * detected.
     */
    private static double[] closestTo(double[][] basis, double[] heading) {
        double[] projected = new double[48];
        for (double[] b : basis) {
            double d = dot(b, heading);
            for (int k = 0; k < 48; k++) {
                projected[k] += d * b[k];
            }
        }
        double n = Linear.norm(projected);
        if (n < 1e-6) {
            return null;
        }
        for (int k = 0; k < 48; k++) {
            projected[k] /= n;
        }
        return projected;
    }

    private static double dot(double[] a, double[] b) {
        double acc = 0;
        for (int k = 0; k < a.length; k++) {
            acc += a[k] * b[k];
        }
        return acc;
    }

    private static double frobenius(double[][][] a, double[][][] b) {
        double acc = 0;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                double d = JitterbugGeometry.distance(a[i][j], b[i][j]);
                acc += d * d;
            }
        }
        return Math.sqrt(acc);
    }

    /** Modified Gram-Schmidt, dropping vectors that collapse. */
    private static List<double[]> orthonormalise(List<double[]> vectors,
                                                 double tol) {
        List<double[]> out = new ArrayList<>();
        for (double[] v : vectors) {
            double[] w = v.clone();
            for (double[] u : out) {
                double d = dot(w, u);
                for (int k = 0; k < w.length; k++) {
                    w[k] -= d * u[k];
                }
            }
            double n = Linear.norm(w);
            if (n > tol) {
                for (int k = 0; k < w.length; k++) {
                    w[k] /= n;
                }
                out.add(w);
            }
        }
        return out;
    }

    /**
     * Rigidly align {@code p} onto {@code q} (Kabsch) and return the moved
     * {@code p}. Reflections are excluded, so this is a rotation plus a
     * translation and never a mirror.
     */
    static double[][][] align(double[][][] p, double[][][] q) {
        double[] pc = centroid(p);
        double[] qc = centroid(q);
        double[][] h = new double[3][3];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        h[a][b] += (p[i][j][a] - pc[a]) * (q[i][j][b] - qc[b]);
                    }
                }
            }
        }
        SingularValueDecomposition svd = new SingularValueDecomposition(new Array2DRowRealMatrix(h,
                                                                                                true));
        RealMatrix u = svd.getU();
        RealMatrix vt = svd.getVT();
        RealMatrix candidate = vt.transpose().multiply(u.transpose());
        double det = new org.apache.commons.math3.linear.LUDecomposition(candidate).getDeterminant();
        RealMatrix flip = new Array2DRowRealMatrix(new double[][] { { 1, 0, 0 },
                                                                    { 0, 1, 0 },
                                                                    { 0, 0,
                                                                      Math.signum(det) } });
        RealMatrix r = vt.transpose().multiply(flip).multiply(u.transpose());

        double[][][] out = new double[8][3][3];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                for (int a = 0; a < 3; a++) {
                    double acc = 0;
                    for (int b = 0; b < 3; b++) {
                        acc += r.getEntry(a, b) * (p[i][j][b] - pc[b]);
                    }
                    out[i][j][a] = acc + qc[a];
                }
            }
        }
        return out;
    }

    /**
     * Minimum-norm Newton projection of {@code x} onto {@code residual = 0}.
     * Iterated to convergence; returns the projected configuration.
     */
    static double[][][] project(double[][][] x) {
        double[][][] current = x;
        for (int it = 0; it < 40; it++) {
            double[] r = JitterbugLinkage.residual(current);
            if (Linear.norm(r) < 1e-14) {
                break;
            }
            double[][] j = JitterbugLinkage.jacobian(current);
            RealMatrix jm = new Array2DRowRealMatrix(j, true);
            RealMatrix jjt = jm.multiply(jm.transpose());
            RealVector y = new SingularValueDecomposition(jjt).getSolver()
                                                              .solve(new ArrayRealVector(r,
                                                                                         false)
                                                                                             .mapMultiply(-1));
            double[] delta = jm.transpose().operate(y).toArray();
            current = JitterbugLinkage.applyBodyMotions(current, delta);
        }
        return current;
    }

    /**
     * March from {@code start} toward {@code target} in {@code waypoints} steps,
     * projecting exactly onto the variety at every one.
     */
    static Result walk(double[][][] start, double[][][] target, int waypoints) {
        double[][][] x = start;
        double worst = 0;
        int collided = 0;
        for (int step = 1; step <= waypoints; step++) {
            double s = (double) step / waypoints;
            double[][][] aligned = align(target, x);
            double[][][] way = blend(x, aligned, s);
            for (int inner = 0; inner < 3; inner++) {
                x = attract(x, way);
            }
            x = project(x);
            worst = Math.max(worst,
                             Linear.norm(JitterbugLinkage.residual(x)));
            if (Interpenetration.count(x) > 0) {
                collided++;
            }
        }
        double[][][] aligned = align(target, x);
        double rms = 0;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                double d = JitterbugGeometry.distance(x[i][j], aligned[i][j]);
                rms += d * d;
            }
        }
        return new Result(x, Math.sqrt(rms / 24.0), worst, collided, waypoints);
    }

    /**
     * One linearised step toward {@code way}, heavily weighted to stay on the
     * variety.
     */
    private static double[][][] attract(double[][][] x, double[][][] way) {
        double[] r = JitterbugLinkage.residual(x);
        double[][] jr = JitterbugLinkage.jacobian(x);
        double[][] jp = JitterbugLinkage.positionJacobian(x);

        double[][] a = new double[108][48];
        double[] b = new double[108];
        for (int row = 0; row < 36; row++) {
            for (int c = 0; c < 48; c++) {
                a[row][c] = HINGE_WEIGHT * jr[row][c];
            }
            b[row] = -HINGE_WEIGHT * r[row];
        }
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                for (int d = 0; d < 3; d++) {
                    int src = 9 * i + 3 * j + d;
                    int row = 36 + src;
                    System.arraycopy(jp[src], 0, a[row], 0, 48);
                    b[row] = way[i][j][d] - x[i][j][d];
                }
            }
        }
        RealVector delta = new QRDecomposition(new Array2DRowRealMatrix(a, true),
                                               1e-12).getSolver()
                                                     .solve(new ArrayRealVector(b,
                                                                                false));
        return JitterbugLinkage.applyBodyMotions(x, delta.toArray());
    }

    private static double[][][] blend(double[][][] from, double[][][] to,
                                      double s) {
        double[][][] out = new double[8][3][3];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    out[i][j][k] = (1 - s) * from[i][j][k] + s * to[i][j][k];
                }
            }
        }
        return out;
    }

    private static double[] centroid(double[][][] x) {
        double[] c = new double[3];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    c[k] += x[i][j][k] / 24.0;
                }
            }
        }
        return c;
    }
}
