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

import static com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.corners;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * Six internal degrees of freedom, and all six extend to finite motions.
 *
 * <p>
 * This is the retraction of the "one degree of freedom" premise. Twelve shared
 * vertices take 48 body DOF down to 12, not to 1 — the recorded inference was a
 * non-sequitur, and the 1-DOF claim had never actually been measured. Six of the
 * twelve are the global rigid motions, leaving six internal.
 *
 * <p>
 * A null-space count alone would only establish first-order flexes. The
 * continuation below is what upgrades it: each internal direction is tracked as
 * a curve, with an exact Newton projection back onto the constraint set at every
 * step, so the six directions are shown to be genuine finite motions rather than
 * infinitesimal mechanisms.
 *
 * <p>
 * The consequence: the symmetric jitterbug path is <b>one curve inside a
 * six-dimensional variety</b>, and its transverse stability is an open question,
 * not an assumption.
 *
 * @author halhildebrand
 */
public class JitterbugDegreesOfFreedomTest {

    /**
     * Six internal directions at every generic angle, after the six global rigid
     * motions are removed.
     */
    @Test
    public void thereAreSixInternalDegreesOfFreedom() {
        for (double a : new double[] { 0, 17, 30, 45, 60, 137 }) {
            assertEquals("internal null directions at a=" + a, 6,
                         VarietyWalk.internalNullSpace(corners(a)).length);
        }
        // And seven at the two rank drops, which is what "branch point" means
        // here: the tangent space gains a dimension.
        for (double a : new double[] { 90, 270 }) {
            assertEquals("internal null directions at the rank drop a=" + a, 7,
                         VarietyWalk.internalNullSpace(corners(a)).length);
        }
    }

    /**
     * The internal basis really is internal: orthonormal, annihilated by the
     * Jacobian, and orthogonal to every global rigid motion. Without the last
     * clause the count of six would be meaningless.
     */
    @Test
    public void theInternalBasisIsOrthonormalNullAndGlobalFree() {
        double[][][] x = corners(30.0);
        double[][] basis = VarietyWalk.internalNullSpace(x);
        double[][] j = JitterbugLinkage.jacobian(x);
        double[][] global = JitterbugLinkage.globalRigidMotions(x);

        for (int p = 0; p < basis.length; p++) {
            assertEquals("unit length", 1.0, Linear.norm(basis[p]), 1e-9);
            assertTrue("must be a null direction",
                       Linear.norm(Linear.apply(j, basis[p])) < 1e-8);
            for (double[] g : global) {
                double d = 0;
                for (int k = 0; k < 48; k++) {
                    d += basis[p][k] * g[k] / Linear.norm(g);
                }
                assertTrue("must carry no global component, saw " + d,
                           Math.abs(d) < 1e-7);
            }
            for (int q = p + 1; q < basis.length; q++) {
                double d = 0;
                for (int k = 0; k < 48; k++) {
                    d += basis[p][k] * basis[q][k];
                }
                assertEquals("mutually orthogonal", 0.0, d, 1e-9);
            }
        }
    }

    /**
     * All six internal directions extend to finite motions: 40 continuation
     * steps of 0.02 each, with an exact projection every step, and none of them
     * dies.
     */
    @Test
    public void allSixInternalModesExtendToFiniteMotions() {
        double[][][] x = corners(30.0);
        double[][] basis = VarietyWalk.internalNullSpace(x);
        assertEquals(6, basis.length);
        for (int m = 0; m < 6; m++) {
            double arc = VarietyWalk.continueAlong(x, basis[m], 0.02, 40);
            assertTrue("internal mode " + m
                       + " must extend to a finite motion; arclength " + arc,
                       arc > 0.5);
        }
    }

    /**
     * The continuation must be able to report a dead end, or "arclength &gt;
     * 0.5" means nothing. A direction transverse to the variety is annihilated
     * by the corrector and produces essentially no motion.
     */
    @Test
    public void continuationReportsZeroForATransverseDirection() {
        double[][][] x = corners(30.0);
        double[][] basis = VarietyWalk.internalNullSpace(x);
        double[][] global = JitterbugLinkage.globalRigidMotions(x);

        // Build a direction with no component in the null space at all.
        java.util.Random rng = new java.util.Random(20260811L);
        double[] transverse = new double[48];
        for (int k = 0; k < 48; k++) {
            transverse[k] = rng.nextDouble() - 0.5;
        }
        for (double[] b : basis) {
            subtractProjection(transverse, b);
        }
        for (double[] g : global) {
            subtractProjection(transverse, g);
        }
        assertTrue("the transverse seed must survive the stripping",
                   Linear.norm(transverse) > 1e-3);

        assertEquals("a transverse direction must go nowhere", 0.0,
                     VarietyWalk.continueAlong(x, transverse, 0.02, 40), 1e-9);
    }

    /**
     * The one-parameter symmetric family is one of those finite motions, and
     * therefore genuinely a curve in the variety rather than an artefact of the
     * parameterisation: the whole path satisfies the constraints and its tangent
     * lies in the internal null space.
     */
    @Test
    public void theSymmetricPathIsOneCurveInTheSixDimensionalVariety() {
        for (double a : new double[] { 5, 30, 55, 100 }) {
            double h = 1e-6;
            double[][][] plus = corners(a + h);
            double[][][] minus = corners(a - h);
            double[][][] x = corners(a);

            // Fit the path's corner velocity to a rigid body motion per face and
            // check the fit is exact, then check the fitted 48-vector is null.
            double[] velocity = new double[72];
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int d = 0; d < 3; d++) {
                        velocity[9 * i + 3 * j
                                 + d] = (plus[i][j][d] - minus[i][j][d])
                                        / (2 * h);
                    }
                }
            }
            double[] fitted = leastSquares(JitterbugLinkage.positionJacobian(x),
                                           velocity);
            double residual = 0;
            double[] predicted = Linear.apply(JitterbugLinkage
                                                              .positionJacobian(x),
                                              fitted);
            for (int k = 0; k < 72; k++) {
                residual = Math.max(residual,
                                    Math.abs(predicted[k] - velocity[k]));
            }
            assertTrue("the path's motion must be rigid per face at a=" + a
                       + ", residual " + residual, residual < 1e-7);
            assertTrue("the path tangent must lie in the null space at a=" + a,
                       Linear.norm(Linear.apply(JitterbugLinkage.jacobian(x),
                                                fitted)) < 1e-7);
            assertTrue("and must be a nonzero direction",
                       Linear.norm(fitted) > 1e-3);
        }
    }

    private static double[] leastSquares(double[][] a, double[] b) {
        return new org.apache.commons.math3.linear.QRDecomposition(new org.apache.commons.math3.linear.Array2DRowRealMatrix(a,
                                                                                                                            true)).getSolver()
                                                                                                                                  .solve(new org.apache.commons.math3.linear.ArrayRealVector(b,
                                                                                                                                                                                            false))
                                                                                                                                  .toArray();
    }

    private static void subtractProjection(double[] v, double[] direction) {
        double n = Linear.norm(direction);
        double d = 0;
        for (int k = 0; k < v.length; k++) {
            d += v[k] * direction[k] / n;
        }
        for (int k = 0; k < v.length; k++) {
            v[k] -= d * direction[k] / n;
        }
    }
}
