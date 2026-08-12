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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.junit.Test;

/**
 * {@link InternalMassMetric}: {@code M(x) = P^T diag(m) P} restricted to the
 * internal tangent space (the 36x48 constraint Jacobian's null space, minus the
 * six global rigid motions -- see {@link VarietyWalk#internalNullSpace}).
 *
 * <p>
 * CHOICE-FREE SCAFFOLDING (bead inviscid-6dp): no potential energy V is chosen
 * anywhere in this suite. The centrepiece is the validation anchor: restricted
 * to the symmetric 1-DOF direction, the metric built here from general 6-DOF
 * machinery must reproduce the independently pinned
 * {@link FreeDynamics#effectiveMass(double)} under the point-mass model, and
 * the corresponding 9:1 (not 3:1) swing under uniform triangular laminae. Both
 * numbers were established by a completely different route in
 * {@code FreeDynamicsTest} (central-differencing {@link JitterbugGeometry}
 * directly); reproducing them from the new tangent-space machinery is a real
 * cross-check, not a restatement.
 *
 * @author halhildebrand
 */
public class InternalMassMetricTest {

    /**
     * Generic angles: rank is 36 (internal DOF exactly 6) at every one of
     * these, per {@code JitterbugLinkageTest.theOnlyRankDropsOnTheCircleAreAtNinetyAndTwoSeventy}.
     * Deliberately skips 90 and 270, the two isolated rank-drop points.
     */
    private static final double[] GENERIC_ANGLES = { 5, 17, 30, 45, 60, 75, 85,
                                                      95, 120, 150, 179, 200,
                                                      230, 260, 265, 275, 300,
                                                      330, 355 };

    /**
     * Per-corner point mass reproducing {@link FreeDynamics#effectiveMass}
     * exactly: {@code FreeDynamicsTest.effectiveMassIsTheKineticFormOfTheCornerParameterisation}
     * measured {@code sum|dx_ij/da_rad|^2 = 48 * M_eff(a)} over the 24 corners,
     * so a uniform corner weight of {@code 1/48} makes
     * {@code sum w*|dx_ij/da|^2} equal {@code M_eff(a)} directly.
     */
    private static InternalMassMetric.PointMasses pointMassModel() {
        double[] m = new double[24];
        java.util.Arrays.fill(m, 1.0 / 48.0);
        return new InternalMassMetric.PointMasses(m);
    }

    /**
     * Per-face lamina mass. Derived (see {@link InternalMassMetric} javadoc) so
     * that the same total-mass convention as {@link #pointMassModel()} —
     * {@code faceMass = 1/16} per face — is what the corner+centroid
     * realisation needs to hit the pinned numbers: with {@code k=1/12},
     * {@code w_corner = k*faceMass = 1/192} and
     * {@code w_centroid = faceMass*(1-3k) = 3/64}, giving spin coefficient
     * {@code 16*w_corner = 1/12} exactly.
     */
    private static InternalMassMetric.UniformLaminae laminaModel() {
        double[] m = new double[8];
        java.util.Arrays.fill(m, 1.0 / 16.0);
        return new InternalMassMetric.UniformLaminae(m);
    }

    /**
     * The exact body-motion velocity of the symmetric jitterbug family at
     * {@code aDeg}, for unit {@code adot} (radians per unit time). Face
     * {@code i} spins about its own fixed axis {@code u_i} at rate
     * {@code sigma_i} (read off {@link JitterbugGeometry#corners(double)}:
     * rotation angle is {@code sigma_i*(a-60deg)}) and its centroid travels
     * along {@code u_i} as {@code Z*cos(a_rad)}, so
     * {@code d(centroid)/da_rad = -Z*sin(a_rad)*u_i}. Both are exact closed
     * forms, not finite differences -- {@link #symmetricGeneratorMatchesFiniteDifferenceCornerVelocity()}
     * is the check that they are right.
     */
    private static double[] symmetricGenerator(double aDeg) {
        double aRad = Math.toRadians(aDeg);
        double[] z = new double[48];
        for (int i = 0; i < 8; i++) {
            double[] u = JitterbugGeometry.faceAxis(i);
            double sigma = Math.signum(u[0] * u[1] * u[2]);
            for (int k = 0; k < 3; k++) {
                z[3 * i + k] = sigma * u[k];
                z[24 + 3 * i + k] = -JitterbugGeometry.Z * Math.sin(aRad) * u[k];
            }
        }
        return z;
    }

    /**
     * The symmetric generator, expressed in the internal-null-space basis, and
     * the reduced metric's value along it: {@code coeff^T M coeff}, which is
     * {@code M_eff(a)} exactly when {@code coeff} really is
     * {@code d(internal coordinates)/da} -- i.e. when {@code symmetricGenerator}
     * lies (almost) entirely inside {@code span(basis)}.
     */
    private static double measuredEffectiveMass(double aDeg,
                                                InternalMassMetric.MassModel model) {
        double[][][] x = JitterbugGeometry.corners(aDeg);
        double[] zSym = symmetricGenerator(aDeg);
        double[][] basis = VarietyWalk.internalNullSpace(x);
        double[] coeff = new double[basis.length];
        for (int k = 0; k < basis.length; k++) {
            coeff[k] = dot(zSym, basis[k]);
        }
        double[][] m = InternalMassMetric.metric(x, model);
        double acc = 0;
        for (int p = 0; p < basis.length; p++) {
            for (int q = 0; q < basis.length; q++) {
                acc += coeff[p] * m[p][q] * coeff[q];
            }
        }
        return acc;
    }

    private static double dot(double[] a, double[] b) {
        double acc = 0;
        for (int k = 0; k < a.length; k++) {
            acc += a[k] * b[k];
        }
        return acc;
    }

    /**
     * {@code symmetricGenerator} must actually generate the corner motion it
     * claims to: {@code positionJacobian(x) * symmetricGenerator(a)} has to
     * agree with a direct central difference of
     * {@link JitterbugGeometry#corners(double)}. Without this check a sign or
     * axis error in {@code symmetricGenerator} could silently corrupt every
     * other test in this class.
     */
    @Test
    public void symmetricGeneratorMatchesFiniteDifferenceCornerVelocity() {
        double hDeg = 1e-3;
        double worst = 0;
        for (double aDeg : GENERIC_ANGLES) {
            double[][][] plus = JitterbugGeometry.corners(aDeg + hDeg);
            double[][][] minus = JitterbugGeometry.corners(aDeg - hDeg);
            double den = 2 * Math.toRadians(hDeg);
            double[][][] x = JitterbugGeometry.corners(aDeg);
            double[] predicted = Linear.apply(JitterbugLinkage.positionJacobian(x),
                                              symmetricGenerator(aDeg));
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int d = 0; d < 3; d++) {
                        double fd = (plus[i][j][d] - minus[i][j][d]) / den;
                        double got = predicted[9 * i + 3 * j + d];
                        worst = Math.max(worst, Math.abs(fd - got));
                    }
                }
            }
        }
        assertTrue("symmetricGenerator vs finite difference, worst " + worst,
                   worst < 1e-5);
    }

    /**
     * Item 1 of the validation anchor: restricted to the symmetric direction,
     * the point-mass metric built from the general 6-DOF machinery here must
     * reproduce {@link FreeDynamics#effectiveMass(double)} exactly, at every
     * generic angle -- not just the ones the original derivation happened to
     * sample.
     */
    @Test
    public void symmetricPointMassMetricReproducesTheMeasuredEffectiveMass() {
        InternalMassMetric.PointMasses model = pointMassModel();
        double worst = 0;
        for (double aDeg : GENERIC_ANGLES) {
            double measured = measuredEffectiveMass(aDeg, model);
            double target = FreeDynamics.effectiveMass(Math.toRadians(aDeg));
            worst = Math.max(worst, Math.abs(measured - target));
            assertEquals("M_eff at a=" + aDeg, target, measured, 1e-6);
        }
        assertTrue("agreement should be tight, worst " + worst, worst < 1e-6);
    }

    /**
     * Item 2 of the validation anchor: uniform triangular laminae must swing
     * 9:1, not 3:1, and the pointwise formula is
     * {@code AXIAL_COEFFICIENT*sin^2(a) + 1/12} -- the same axial term as the
     * point-mass model (centroid motion, untouched by how mass is distributed
     * within a rigid face) with the corner point-mass spin coefficient
     * {@code 1/3} replaced by the lamina's {@code 1/12}.
     */
    @Test
    public void symmetricLaminaMetricGivesTheNineToOneSwing() {
        InternalMassMetric.UniformLaminae model = laminaModel();
        double lo = Double.MAX_VALUE;
        double hi = 0;
        for (int k = 0; k <= 360; k++) {
            // Unlike JitterbugLinkageTest's rank sweep, 90 and 270 are not
            // skipped here: the internal space is 7-dimensional there (one
            // extra direction from the local rank drop), but the symmetric
            // generator's coefficients in that basis still reproduce the
            // exact peak (measured: M_eff(90) = 0.7500000000000007, matching
            // AXIAL_COEFFICIENT + 1/12 to 1e-15). Skipping them would miss
            // the true maximum and turn the 9:1 ratio below into a
            // discretisation artefact instead of a measurement.
            double aDeg = k;
            double measured = measuredEffectiveMass(aDeg, model);
            // Axial term via effectiveMass(a) - effectiveMass(0), not by
            // reading AXIAL_COEFFICIENT directly: AXIAL_COEFFICIENT is a JLS
            // 15.28 compile-time constant, so a direct read here would be
            // inlined into this test's class file and could go stale under
            // incremental compilation if the constant changed (bead
            // inviscid-hep). effectiveMass(0) is exactly SPIN_COEFFICIENT
            // (sin(0) == 0.0 exactly, so the axial term vanishes identically),
            // which isolates the axial contribution through the method call
            // alone -- the same way every other target in this file routes
            // through effectiveMass() rather than the raw constants.
            double axial = FreeDynamics.effectiveMass(Math.toRadians(aDeg))
                            - FreeDynamics.effectiveMass(0.0);
            double target = axial + 1.0 / 12.0;
            assertEquals("lamina M_eff at a=" + aDeg, target, measured, 1e-6);
            lo = Math.min(lo, measured);
            hi = Math.max(hi, measured);
        }
        assertEquals("uniform laminae swing 9:1", 9.0, hi / lo, 1e-9);
        assertEquals("lamina minimum is 1/12", 1.0 / 12.0, lo, 1e-6);
    }

    /**
     * The two mass models must disagree -- otherwise the parameterisation is a
     * no-op and the 3:1 vs 9:1 distinction the bead exists to make explicit
     * would be a decoration, not a real choice.
     */
    @Test
    public void theTwoMassModelsGiveDifferentMetrics() {
        double[][][] x = JitterbugGeometry.corners(30.0);
        double[][] point = InternalMassMetric.metric(x, pointMassModel());
        double[][] lamina = InternalMassMetric.metric(x, laminaModel());
        double worst = 0;
        for (int i = 0; i < point.length; i++) {
            for (int j = 0; j < point[i].length; j++) {
                worst = Math.max(worst, Math.abs(point[i][j] - lamina[i][j]));
            }
        }
        assertTrue("the two mass models must give measurably different metrics, "
                   + "worst difference " + worst, worst > 1e-3);
    }

    /**
     * Item 3 of the deliverable: the internal tangent space has dimension 6 at
     * generic angles, measured (not assumed) at every angle in
     * {@link #GENERIC_ANGLES}.
     */
    @Test
    public void internalTangentSpaceHasDimensionSixAtGenericAngles() {
        for (double aDeg : GENERIC_ANGLES) {
            double[][] basis = VarietyWalk.internalNullSpace(JitterbugGeometry.corners(aDeg));
            assertEquals("internal DOF at a=" + aDeg, 6, basis.length);
        }
    }

    /**
     * The internal basis must be orthogonal to all six global rigid motions --
     * "the null space MINUS the global rigid motions" is a real subtraction,
     * not a restatement of the null space. A bug in
     * {@code VarietyWalk.internalNullSpace}'s Gram-Schmidt strip would show up
     * here as a nonzero overlap.
     */
    @Test
    public void globalRigidMotionsAreExcludedFromTheInternalBasis() {
        double worst = 0;
        for (double aDeg : GENERIC_ANGLES) {
            double[][][] x = JitterbugGeometry.corners(aDeg);
            double[][] basis = VarietyWalk.internalNullSpace(x);
            double[][] global = JitterbugLinkage.globalRigidMotions(x);
            for (double[] b : basis) {
                for (double[] g : global) {
                    worst = Math.max(worst, Math.abs(dot(b, g)));
                }
            }
        }
        assertTrue("internal basis vectors must be orthogonal to the six global "
                   + "rigid motions, worst overlap " + worst, worst < 1e-8);
    }

    /**
     * The metric is positive-definite at a generic configuration, for both
     * mass models -- checked via commons-math3's own Cholesky decomposition
     * rather than hand-rolled eigenvalue signs, so this reuses the same
     * positive-definiteness notion {@link GeneralizedEigensolver} relies on.
     *
     * <p>
     * The substance of this test is the Cholesky decomposition at the bottom
     * of the loop. The elementwise comparison above it is labelled as a
     * guard, not as a result.
     */
    @Test
    public void metricIsSymmetricAndPositiveDefiniteAtGenericAngles() {
        for (double aDeg : GENERIC_ANGLES) {
            double[][][] x = JitterbugGeometry.corners(aDeg);
            for (InternalMassMetric.MassModel model : new InternalMassMetric.MassModel[] {
                                                                                          pointMassModel(),
                                                                                          laminaModel() }) {
                double[][] m = InternalMassMetric.metric(x, model);
                // GUARD, not a verification: InternalMassMetric.metric()
                // writes the identical double into m[a][b] and m[b][a] (its
                // accumulation loop computes acc once per pair and assigns
                // it to both slots), so this elementwise comparison cannot
                // fail while that write pattern holds, regardless of whether
                // acc's value is physically correct. It earns its place only
                // as a regression guard against a future refactor that stops
                // populating both halves identically -- the real check is
                // the Cholesky decomposition below.
                for (int i = 0; i < m.length; i++) {
                    for (int j = 0; j < m.length; j++) {
                        assertEquals("metric must be symmetric at a=" + aDeg,
                                     m[i][j], m[j][i], 1e-12);
                    }
                }
                // Throws NonPositiveDefiniteMatrixException if not PD -- that
                // is the assertion.
                new CholeskyDecomposition(new Array2DRowRealMatrix(m, true));
            }
        }
    }

    /**
     * {@code cornerMass} and {@code faceMass} must be validated eagerly, not
     * accepted and misbehave later on a mismatched-length array.
     */
    @Test(expected = IllegalArgumentException.class)
    public void pointMassesRejectsTheWrongCornerCount() {
        new InternalMassMetric.PointMasses(new double[23]);
    }

    @Test(expected = IllegalArgumentException.class)
    public void uniformLaminaeRejectsTheWrongFaceCount() {
        new InternalMassMetric.UniformLaminae(new double[7]);
    }
}
