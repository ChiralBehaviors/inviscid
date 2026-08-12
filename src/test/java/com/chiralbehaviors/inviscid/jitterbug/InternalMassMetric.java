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

/**
 * The internal mass metric {@code M(x) = P^T diag(m) P}, restricted to the
 * internal tangent space: the null space of the 36x48 constraint Jacobian
 * ({@link JitterbugLinkage#jacobian}), minus the six global rigid motions
 * ({@link JitterbugLinkage#globalRigidMotions}) -- see
 * {@link VarietyWalk#internalNullSpace}, which this class reuses rather than
 * duplicating.
 *
 * <p>
 * <b>CHOICE-FREE SCAFFOLDING</b> (bead inviscid-6dp). This class picks no
 * potential energy {@code V}; it only builds the metric that a Hessian needs
 * in order to become frequencies via {@link GeneralizedEigensolver}. The
 * mass distribution is the only choice made here, and it is made twice, both
 * ways, as first-class alternatives -- see {@link PointMasses} and
 * {@link UniformLaminae} -- rather than one default with an option bolted on.
 *
 * <p>
 * <b>Point masses.</b> {@link PointMasses} places independent masses at the
 * 24 real corners and is {@code M(x) = P^T diag(m) P} literally, with
 * {@code P} the existing {@link JitterbugLinkage#positionJacobian}.
 *
 * <p>
 * <b>Uniform triangular laminae.</b> A rigid plate's mass is not concentrated
 * at its three vertices, and no choice of (possibly unequal) point masses at
 * just those three vertices can reproduce one: for an equilateral triangle
 * centred at its own centroid, {@code m0 r0 + m1 r1 + m2 r2 = 0} (mass
 * balances at the centroid, which every rigid motion here pivots about) forces
 * {@code m0 = m1 = m2}, and equal vertex masses always carry polar second
 * moment {@code L^2/3} per unit mass -- the point-mass value, never the
 * plate's {@code L^2/12}. {@link UniformLaminae} instead adds a fourth,
 * <em>virtual</em> mass point at the face's own centroid. Its velocity is not
 * approximated: for any rigid motion of a triangle, the centroid velocity is
 * exactly the average of the three vertices' velocities (their position
 * vectors about the centroid sum to zero, so the cross term in the average
 * drops out identically). Writing {@code k = 1/12} for the lamina's spin
 * coefficient (vs. {@code 1/3} for point masses -- see
 * {@link FreeDynamics}), matching a face's total mass and polar second moment
 * exactly requires
 *
 * <pre>
 *     w_corner   = k * faceMass       (equal at all three corners)
 *     w_centroid = faceMass * (1 - 3k)
 * </pre>
 *
 * (Substituting {@code k = 1/3} recovers {@code w_centroid = 0} -- the
 * point-mass model is the {@code k = 1/3} special case of this family, though
 * {@link PointMasses} implements it directly rather than through this
 * derivation, per the deliverable's literal formula.) This is derived, not
 * fitted: {@code InternalMassMetricTest} reproduces the independently-measured
 * {@link FreeDynamics#effectiveMass} and its lamina counterpart from it.
 *
 * @author halhildebrand
 */
final class InternalMassMetric {

    /** A mass distribution over the eight rigid triangular faces. */
    sealed interface MassModel permits PointMasses, UniformLaminae {
    }

    /**
     * Independent point masses at the 24 real corners, indexed
     * {@code 3*face + corner}, matching {@link JitterbugLinkage#positionJacobian}'s
     * row layout.
     */
    record PointMasses(double[] cornerMass) implements MassModel {
        PointMasses {
            if (cornerMass.length != 24) {
                throw new IllegalArgumentException("cornerMass must have 24 entries, has "
                                                    + cornerMass.length);
            }
        }
    }

    /**
     * A uniform triangular lamina per face, indexed by face. See the class
     * javadoc for the dynamically-equivalent corner+centroid realisation.
     */
    record UniformLaminae(double[] faceMass) implements MassModel {
        UniformLaminae {
            if (faceMass.length != 8) {
                throw new IllegalArgumentException("faceMass must have 8 entries, has "
                                                    + faceMass.length);
            }
        }

        /**
         * A uniform triangle's polar second moment about its own centroid,
         * per unit mass, is {@code L^2/12} ({@code L^2/3} for equal point
         * masses at the three corners) -- see {@link FreeDynamics}.
         */
        static final double K = 1.0 / 12.0;
    }

    private InternalMassMetric() {
    }

    /**
     * The metric on the internal tangent space at configuration {@code x}:
     * {@code n x n}, where {@code n = internalNullSpace(x).length} -- six at a
     * generic configuration, seven at the two rank-drop angles
     * ({@code JitterbugLinkageTest.theOnlyRankDropsOnTheCircleAreAtNinetyAndTwoSeventy}).
     * This method does not special-case {@code n}; the dimension is whatever
     * the null space measures.
     */
    static double[][] metric(double[][][] x, MassModel model) {
        double[] cornerWeight = new double[24];
        double[] centroidWeight = new double[8];
        weights(model, cornerWeight, centroidWeight);

        double[][] basis = VarietyWalk.internalNullSpace(x);
        double[][] p = JitterbugLinkage.positionJacobian(x);
        int n = basis.length;

        // cornerVelocity[k] is P * basis[k]: the 72 corner-velocity components
        // (row = 9*face + 3*corner + xyz) induced by the k'th internal basis
        // direction. centroidVelocity[k] is the corresponding face-centroid
        // velocity, exactly the average of its three corners' velocities --
        // see the class javadoc.
        double[][] cornerVelocity = new double[n][];
        double[][] centroidVelocity = new double[n][24];
        for (int k = 0; k < n; k++) {
            cornerVelocity[k] = Linear.apply(p, basis[k]);
            for (int f = 0; f < 8; f++) {
                for (int d = 0; d < 3; d++) {
                    double sum = 0;
                    for (int j = 0; j < 3; j++) {
                        sum += cornerVelocity[k][9 * f + 3 * j + d];
                    }
                    centroidVelocity[k][3 * f + d] = sum / 3.0;
                }
            }
        }

        double[][] m = new double[n][n];
        for (int a = 0; a < n; a++) {
            for (int b = a; b < n; b++) {
                double acc = 0;
                for (int c = 0; c < 24; c++) {
                    int face = c / 3;
                    int corner = c % 3;
                    for (int d = 0; d < 3; d++) {
                        int row = 9 * face + 3 * corner + d;
                        acc += cornerWeight[c] * cornerVelocity[a][row]
                               * cornerVelocity[b][row];
                    }
                }
                for (int f = 0; f < 8; f++) {
                    for (int d = 0; d < 3; d++) {
                        acc += centroidWeight[f] * centroidVelocity[a][3 * f + d]
                               * centroidVelocity[b][3 * f + d];
                    }
                }
                m[a][b] = acc;
                m[b][a] = acc;
            }
        }
        return m;
    }

    private static void weights(MassModel model, double[] cornerWeight,
                                double[] centroidWeight) {
        if (model instanceof PointMasses pm) {
            System.arraycopy(pm.cornerMass(), 0, cornerWeight, 0, 24);
        } else if (model instanceof UniformLaminae ul) {
            double[] faceMass = ul.faceMass();
            for (int f = 0; f < 8; f++) {
                double w = UniformLaminae.K * faceMass[f];
                cornerWeight[3 * f] = w;
                cornerWeight[3 * f + 1] = w;
                cornerWeight[3 * f + 2] = w;
                centroidWeight[f] = faceMass[f] * (1 - 3 * UniformLaminae.K);
            }
        } else {
            throw new IllegalArgumentException("unknown mass model " + model);
        }
    }
}
