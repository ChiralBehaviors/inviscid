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

import static com.chiralbehaviors.inviscid.jitterbug.ReducedCoordinates.NQ;
import static com.chiralbehaviors.inviscid.jitterbug.ReducedCoordinates.NU;
import static com.chiralbehaviors.inviscid.jitterbug.ReducedCoordinates.NV;
import static com.chiralbehaviors.inviscid.jitterbug.ReducedCoordinates.SLOT;
import static com.chiralbehaviors.inviscid.jitterbug.ReducedCoordinates.VERTEX_MASS;
import static com.chiralbehaviors.inviscid.jitterbug.ReducedCoordinates.body;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.chiralbehaviors.inviscid.jitterbug.ReducedCoordinates.Assembly;

/**
 * The reduced-coordinate array, measured — and checked against the Python
 * harness that gated it first.
 *
 * <p>
 * <b>Why a second implementation is worth its cost.</b>
 * {@code analysis/jitterbug-variety/jb_rc_reduced.py} carries fifteen gate rows
 * on this same model. Nothing forces the two to agree, so when they do it means
 * something, and the load-bearing agreements here are the ones that land on
 * <em>exact rationals</em>: the nine-cell cluster splits its energy
 * {@code 5/37} to the centre and {@code 4/37} to each of eight neighbours, and
 * the three-cell V takes fold rates {@code 57/137} and {@code 12/137} at
 * {@code E = 152/137}. A wrong port does not land near a different rational; it
 * lands on a different one.
 *
 * <p>
 * <b>Rank and nullity live here, not in {@code src/main}.</b> They are how the
 * array is measured rather than how it is defined, which is the same boundary
 * {@link Linear} states for commons-math3.
 *
 * @author halhildebrand
 */
public class ReducedCoordinatesTest {

    private static final int[][] EIGHT_DIAGONALS = { { -1, -1, -1 }, { -1, -1, 1 },
                                                     { -1, 1, -1 }, { -1, 1, 1 },
                                                     { 1, -1, -1 }, { 1, -1, 1 },
                                                     { 1, 1, -1 }, { 1, 1, 1 } };

    /** The three-cell V of {@code ThreeCellAnimation}: two adjacent diagonals. */
    private static final int[][] V_SITES       = { { 1, 1, 1 }, { 1, 1, -1 } };

    /**
     * The body function is {@link JitterbugGeometry}'s, and its two derivatives
     * are the analytic ones.
     *
     * <p>
     * The per-face {@code sigma} is invisible in the positions and wrong in every
     * velocity if it is dropped, so checking the positions alone would prove
     * nothing about the dynamics. The second difference is roundoff limited at
     * {@code eps/h^2}, so agreement there is 1e-6, not 1e-12.
     */
    @Test
    public void bodyIsTheGatedGeometryAndItsDerivativesAreAnalytic() {
        for (double g : new double[] { -47.0, 0.0, 23.7, 60.0, 91.3 }) {
            double[][][] x = JitterbugGeometry.corners(g);
            double[][] v = body(g, 0)[0];
            for (int f = 0; f < 8; f++) {
                for (int c = 0; c < 3; c++) {
                    for (int t = 0; t < 3; t++) {
                        assertEquals("body != JitterbugGeometry.corners at " + g,
                                     x[f][c][t], v[SLOT[f][c]][t], 1e-15);
                    }
                }
            }
        }
        double g0 = 23.7;
        double h = 1e-5;
        double[][][] d = body(g0, 2);
        double[][] plus = body(g0 + Math.toDegrees(h), 0)[0];
        double[][] minus = body(g0 - Math.toDegrees(h), 0)[0];
        double e1 = 0;
        double e2 = 0;
        for (int i = 0; i < NV; i++) {
            for (int t = 0; t < 3; t++) {
                e1 = Math.max(e1, Math.abs(d[1][i][t]
                                           - (plus[i][t] - minus[i][t]) / (2 * h)));
                e2 = Math.max(e2,
                              Math.abs(d[2][i][t]
                                       - (plus[i][t] - 2 * d[0][i][t] + minus[i][t])
                                         / (h * h)));
            }
        }
        assertTrue("dv/dg vs central difference " + e1, e1 < 1e-8);
        assertTrue("d2v/dg2 vs central difference " + e2, e2 < 1e-4);
    }

    /**
     * Congruence is not identity, and it survives into the mass.
     *
     * <p>
     * Twelve vertex identities at every angle, {@code 2/3} each, eight units per
     * cell. At the octahedron those twelve occupy six positions — a build that
     * deduplicated by position would carry four, which is half the cell's
     * inertia and the accounting Fuller objects to.
     */
    @Test
    public void congruentQuantaAreCarriedIntoTheMass() {
        assertEquals(12, NV);
        double total = 0;
        for (int i = 0; i < NV; i++) {
            assertEquals(2.0 / 3.0, VERTEX_MASS[i], 1e-15);
            total += VERTEX_MASS[i];
        }
        assertEquals("eight triangles of unit mass", 8.0, total, 1e-12);

        double[][] octa = body(60.0, 0)[0];
        List<double[]> distinct = new ArrayList<>();
        for (double[] p : octa) {
            boolean seen = false;
            for (double[] r : distinct) {
                if (dist(p, r) < 1e-9) {
                    seen = true;
                    break;
                }
            }
            if (!seen) {
                distinct.add(p);
            }
        }
        assertEquals("twelve identities congruent as six at the octahedron", 6,
                     distinct.size());
    }

    /**
     * With V = 0 the motion is a geodesic, so energy is the only audit there is
     * — and the fixed-step integrator has no structural protection, which is what
     * makes this a test rather than a tautology.
     */
    @Test
    public void energyIsConservedAndTheWeldsStayClosed() {
        Assembly v = ReducedCoordinates.cluster(0.0, V_SITES);
        double[] u = v.foldImpulse(0);
        double[] state = new double[NQ * 3 + NU * 3];
        System.arraycopy(v.initialState(), 0, state, 0, NQ * 3);
        System.arraycopy(u, 0, state, NQ * 3, NU * 3);
        double e0 = sum(v.energyPerCell(v.initialState(), u));
        double drift = 0;
        double weld = 0;
        for (int s = 0; s < 100; s++) {
            v.step(state, 0.1, 4);
            double[] q = new double[NQ * 3];
            double[] uu = new double[NU * 3];
            System.arraycopy(state, 0, q, 0, NQ * 3);
            System.arraycopy(state, NQ * 3, uu, 0, NU * 3);
            drift = Math.max(drift, Math.abs(sum(v.energyPerCell(q, uu)) - e0) / e0);
            weld = Math.max(weld, v.weldResidual(q));
        }
        assertTrue("relative energy drift over t=10 was " + drift, drift < 1e-7);
        assertTrue("weld residual over t=10 was " + weld, weld < 1e-9);
    }

    /**
     * The ellipses are absent from the constraint set, not small in it.
     *
     * <p>
     * A vertex's body-frame position is a function of the fold angle alone, so
     * every vertex stays at perpendicular distance {@code EL/sqrt(3)} from the
     * axis of <em>both</em> faces it belongs to and keeps one Cartesian
     * coordinate identically zero. That is a property of the fixed body function,
     * checked once here over a sweep, and no run can violate it — which is why
     * the only constraint rows are the nine per shared face.
     */
    @Test
    public void theEllipsesAreStructuralAndTheConstraintSetIsWeldsOnly() {
        int[][] incident = new int[NV][2];
        int[] fill = new int[NV];
        for (int f = 0; f < 8; f++) {
            for (int c = 0; c < 3; c++) {
                incident[SLOT[f][c]][fill[SLOT[f][c]]++] = f;
            }
        }
        for (int i = 0; i < NV; i++) {
            assertEquals("each vertex belongs to exactly two faces", 2, fill[i]);
        }
        double radius = JitterbugGeometry.L_EDGE / Math.sqrt(3.0);
        double worstAxis = 0;
        double worstPlane = 0;
        for (int s = 0; s <= 60; s++) {
            double g = -95.0 + s * 190.0 / 60.0;
            double[][] v = body(g, 0)[0];
            for (int i = 0; i < NV; i++) {
                worstPlane = Math.max(worstPlane,
                                      Math.min(Math.abs(v[i][0]),
                                               Math.min(Math.abs(v[i][1]),
                                                        Math.abs(v[i][2]))));
                for (int k = 0; k < 2; k++) {
                    double[] u = JitterbugGeometry.faceAxis(incident[i][k]);
                    double along = v[i][0] * u[0] + v[i][1] * u[1] + v[i][2] * u[2];
                    double perp = 0;
                    for (int t = 0; t < 3; t++) {
                        double e = v[i][t] - along * u[t];
                        perp += e * e;
                    }
                    worstAxis = Math.max(worstAxis,
                                         Math.abs(Math.sqrt(perp) - radius));
                }
            }
        }
        assertTrue("distance to the incident face axes drifted " + worstAxis,
                   worstAxis < 1e-14);
        assertTrue("a vertex left its coordinate plane by " + worstPlane,
                   worstPlane < 1e-14);

        assertEquals("welds and nothing else, six cells", 9 * 5,
                     ReducedCoordinates.chain(6, 30.0).constraintJacobian(
                             ReducedCoordinates.chain(6, 30.0).initialState()).length);
        Assembly nine = ReducedCoordinates.cluster(0.0, EIGHT_DIAGONALS);
        assertEquals("welds and nothing else, nine cells", 9 * 8,
                     nine.constraintJacobian(nine.initialState()).length);
    }

    /**
     * The DOF table, from the reduced coordinates, with a correct
     * global-rigid-motion basis.
     *
     * <p>
     * Cells free to turn: one internal degree of freedom per cell at every size,
     * which is what makes the medium a field. Orientations frozen: the
     * <em>same</em> Jacobian restricted to its centre and fold columns has rank 4
     * per weld instead of 6 and leaves a single coordinate at every size, so a
     * disturbance has nowhere to be. Both rows are column subsets of one matrix.
     */
    @Test
    public void theDofTableIsOnePerCellFreeAndOneInTotalFrozen() {
        for (int n = 1; n <= 7; n++) {
            Assembly a = ReducedCoordinates.chain(n, 30.0);
            double[][] c = a.constraintJacobian(a.initialState());
            int free = 7 * n - rank(c) - 6;
            int frozen = 4 * n - rank(dropOrientationColumns(c, n)) - 3;
            assertEquals("chain of " + n + ", cells may rotate", n, free);
            assertEquals("chain of " + n + ", orientation fixed", 1, frozen);
            if (n > 1) {
                assertEquals("independent rows per weld", 6 * (n - 1), rank(c));
            }
        }
        Assembly nine = ReducedCoordinates.cluster(0.0, EIGHT_DIAGONALS);
        double[][] c9 = nine.constraintJacobian(nine.initialState());
        assertEquals("nine-cell cluster, cells may rotate", 9, 63 - rank(c9) - 6);
        assertEquals("nine-cell cluster, orientation fixed", 1,
                     36 - rank(dropOrientationColumns(c9, 9)) - 3);

        Assembly v = ReducedCoordinates.cluster(0.0, V_SITES);
        double[][] cv = v.constraintJacobian(v.initialState());
        assertEquals("three-cell V, cells may rotate", 3, 21 - rank(cv) - 6);
        assertEquals("three-cell V, orientation fixed", 1,
                     12 - rank(dropOrientationColumns(cv, 3)) - 3);
    }

    /**
     * The global-rigid-motion basis is the correct one, and the defective one is
     * demonstrably not.
     *
     * <p>
     * Translations are {@code dc = t}; rotations are {@code dc = omega x c},
     * {@code omega_k = omega}. An earlier fold-map-rank run in this project used
     * an identity orientation block and omitted the {@code omega x c} term,
     * projected out a six-dimensional space that is not the global motions, and
     * reported a degree-of-freedom column two too high. Two-sided: the wrong
     * basis is not annihilated, so this cannot pass by measuring nothing.
     */
    @Test
    public void theGlobalMotionBasisIsTheCorrectedOne() {
        Assembly a = ReducedCoordinates.chain(6, 30.0);
        double[] q = a.initialState();
        double[][] c = a.constraintJacobian(q);
        double good = worstProduct(c, a.globalMotions(q, false));
        double bad = worstProduct(c, a.globalMotions(q, true));
        assertTrue("correct basis is not in the nullspace: " + good, good < 1e-12);
        assertTrue("defective basis was annihilated too, so this proves nothing: "
                   + bad, bad > 1e-3);
    }

    /**
     * The nine-cell cluster's t = 0 split is exactly {@code 5/37 : 32/37}, which
     * is the strongest single agreement with the Python harness.
     *
     * <p>
     * Before any time has passed, the projection that makes a phase kick on the
     * centre admissible has already put {@code 32/37} of the energy in the shell.
     * One cell cannot fold without its eight neighbours folding, and a rigid
     * constraint has infinite signal speed — so there is no onset lag to find and
     * no wavefront to time.
     */
    @Test
    public void theNineCellSplitIsExactlyFiveThirtySevenths() {
        Assembly nine = ReducedCoordinates.cluster(0.0, EIGHT_DIAGONALS);
        double[] u = nine.foldImpulse(0);
        assertEquals("centre fold rate", 5.0 / 37.0, u[6], 1e-12);
        for (int k = 1; k < 9; k++) {
            assertEquals("shell fold rate, cell " + k, 1.0 / 37.0, u[NU * k + 6],
                         1e-12);
        }
        double[] e = nine.energyPerCell(nine.initialState(), u);
        double total = sum(e);
        assertEquals("centre share", 5.0 / 37.0, e[0] / total, 1e-12);
        for (int k = 1; k < 9; k++) {
            assertEquals("shell share, cell " + k, 4.0 / 37.0, e[k] / total, 1e-12);
        }
    }

    /**
     * The three-cell V — the configuration {@code ThreeCellAnimation} builds — in
     * reduced coordinates, against the harness.
     *
     * <p>
     * Its outer cells sit on two adjacent diagonals at {@code n1 . n2 = 1/3}, and
     * a phase kick on the middle cell projects to fold rates {@code 57/137} and
     * {@code 12/137} at total energy {@code 152/137}. The P/Q symmetry is exact,
     * which is a check on the weld construction as much as on the dynamics.
     */
    @Test
    public void theThreeCellVAgreesWithTheHarness() {
        Assembly v = ReducedCoordinates.cluster(0.0, V_SITES);
        double[] q = v.initialState();
        double n1n2 = 0;
        double na = 0;
        double nb = 0;
        for (int t = 0; t < 3; t++) {
            n1n2 += q[NQ + t] * q[2 * NQ + t];
            na += q[NQ + t] * q[NQ + t];
            nb += q[2 * NQ + t] * q[2 * NQ + t];
        }
        assertEquals("adjacent diagonals", 1.0 / 3.0,
                     n1n2 / Math.sqrt(na * nb), 1e-15);
        assertTrue("welds closed: " + v.weldResidual(q), v.weldResidual(q) < 1e-14);

        double[] u = v.foldImpulse(0);
        assertEquals("middle fold rate", 57.0 / 137.0, u[6], 1e-12);
        assertEquals("outer fold rate P", 12.0 / 137.0, u[NU + 6], 1e-12);
        assertEquals("outer fold rate Q", 12.0 / 137.0, u[2 * NU + 6], 1e-12);

        double[] e = v.energyPerCell(q, u);
        assertEquals("total energy", 152.0 / 137.0, sum(e), 1e-12);
        assertEquals("middle share", 0.4562682, e[0] / sum(e), 1e-6);
        assertEquals("P and Q are symmetric", e[1], e[2], 1e-14);
    }

    /**
     * The two implementations agree on the MOTION, not just on the initial
     * condition — which is the agreement that actually costs something to get.
     *
     * <p>
     * A phase kick on one outer cell of the V, integrated to t = 30 through
     * roughly three full turns of the driven cell's fold angle. The harness, at
     * {@code rtol = 1e-11} with an adaptive DOP853, puts the middle cell's peak
     * share at {@code 0.31264} near {@code t = 21.6} and the far cell's at
     * {@code 0.15751} near {@code t = 26.4}. This is fixed-step Runge-Kutta with
     * a projection back onto the weld manifold — a different integrator, a
     * different language, an independently written model — and it lands in the
     * same place.
     *
     * <p>
     * The transport itself is the point: the disturbance starts on one outer
     * cell, and the far cell ends up with several times the share it was given at
     * {@code t = 0}. It is <em>not</em> a wavefront. The onset is simultaneous,
     * so no time on this trajectory is an arrival time.
     */
    @Test
    public void theTransportAgreesWithTheHarnessAndNotJustTheInitialCondition() {
        Assembly v = ReducedCoordinates.cluster(0.0, V_SITES);
        double[] u = v.foldImpulse(1);
        assertEquals("driven outer fold rate", 0.630967466, u[NU + 6], 1e-9);
        double[] state = new double[NQ * 3 + NU * 3];
        System.arraycopy(v.initialState(), 0, state, 0, NQ * 3);
        System.arraycopy(u, 0, state, NQ * 3, NU * 3);
        double e0 = sum(v.energyPerCell(v.initialState(), u));

        double h = 0.05;
        double middlePeak = 0;
        double middleAt = 0;
        double farPeak = 0;
        double farAt = 0;
        double farStart = 0;
        double drift = 0;
        for (int s = 0; s * h <= 30.0; s++) {
            double[] q = new double[NQ * 3];
            double[] uu = new double[NU * 3];
            System.arraycopy(state, 0, q, 0, NQ * 3);
            System.arraycopy(state, NQ * 3, uu, 0, NU * 3);
            double[] e = v.energyPerCell(q, uu);
            double total = sum(e);
            drift = Math.max(drift, Math.abs(total - e0) / e0);
            if (s == 0) {
                farStart = e[2] / total;
            }
            if (e[0] / total > middlePeak) {
                middlePeak = e[0] / total;
                middleAt = s * h;
            }
            if (e[2] / total > farPeak) {
                farPeak = e[2] / total;
                farAt = s * h;
            }
            v.step(state, h, 4);
        }
        assertEquals("middle cell's peak share", 0.31264, middlePeak, 2e-3);
        assertEquals("and when it peaks", 21.6, middleAt, 0.5);
        assertEquals("far cell's peak share", 0.15751, farPeak, 2e-3);
        assertEquals("and when it peaks", 26.4, farAt, 0.5);
        assertTrue("the far cell must actually gain: " + farStart + " -> " + farPeak,
                   farPeak > 3 * farStart);
        assertTrue("energy drift over t=30 was " + drift, drift < 1e-8);
    }

    private static double dist(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < a.length; i++) {
            s += (a[i] - b[i]) * (a[i] - b[i]);
        }
        return Math.sqrt(s);
    }

    /** The orientation columns removed — the frozen-orientation row of the DOF
     *  table is this same matrix, not a separate measurement. */
    private static double[][] dropOrientationColumns(double[][] c, int n) {
        double[][] out = new double[c.length][4 * n];
        for (int r = 0; r < c.length; r++) {
            int w = 0;
            for (int i = 0; i < NU * n; i++) {
                int within = i % NU;
                if (within >= 3 && within <= 5) {
                    continue;
                }
                out[r][w++] = c[r][i];
            }
        }
        return out;
    }

    /** A lone cell has no welds at all, so the constraint matrix has no rows and
     *  a rank of zero. That is the n = 1 column of the table, not an edge case to
     *  skip: one cell, one fold angle, one degree of freedom either way. */
    private static int rank(double[][] m) {
        return m.length == 0 ? 0 : Linear.rank(m, 1e-8);
    }

    private static double sum(double[] v) {
        double s = 0;
        for (double x : v) {
            s += x;
        }
        return s;
    }

    private static double worstProduct(double[][] c, double[][] g) {
        double worst = 0;
        for (double[] row : c) {
            for (double[] motion : g) {
                double s = 0;
                for (int i = 0; i < row.length; i++) {
                    s += row[i] * motion[i];
                }
                worst = Math.max(worst, Math.abs(s));
            }
        }
        return worst;
    }
}
