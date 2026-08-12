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
 * Free motion of the symmetric jitterbug along its one-parameter path.
 *
 * <p>
 * <b>Mass model: equal point masses at the 24 triangle corners.</b> Stated first
 * because every number below is contingent on it. "Eight rigid triangles" does
 * not by itself fix a mass distribution, and a different one gives a different
 * answer — see the 9:1 paragraph.
 *
 * <p>
 * The kinetic energy of the eight triangles reduces to
 * {@code T = 0.5 * M_eff(a) * adot^2} with
 *
 * <pre>
 *     M_eff(a) / M = C_axial sin^2(a) + C_spin,   C_axial = 2/3,  C_spin = 1/3
 * </pre>
 *
 * so the effective mass swings 3:1 over the path.
 *
 * <p>
 * <b>Where the two coefficients come from.</b> They are not free parameters:
 * they are fixed by {@link JitterbugGeometry}'s corner parameterisation, and
 * {@code FreeDynamicsTest.effectiveMassIsTheKineticFormOfTheCornerParameterisation}
 * measures the tie rather than restating the literal:
 *
 * <pre>
 *     sum over the 24 corners of |dx_ij/da|^2  =  32 sin^2(a) + 16
 *                                              =  48 * [ (2/3) sin^2(a) + 1/3 ]
 * </pre>
 *
 * The 48 is a fixed normalisation and cancels out of every ratio reported here.
 * Term by term: {@code C_axial = Z*Z/2} for the axial-travel amplitude
 * {@link JitterbugGeometry#Z}, and {@code C_spin = L*L/6} for the strut length
 * {@link JitterbugGeometry#L_EDGE}, because the corners of an equilateral
 * triangle of side {@code L} sit at {@code L/sqrt(3)} from its centroid and so
 * carry polar second moment {@code L^2/3} per unit mass.
 *
 * <p>
 * <b>The 3:1 swing belongs to the mass model, not to the linkage.</b> Replace
 * the corner point masses with <b>uniform triangular laminae</b> and the polar
 * second moment per unit mass falls from {@code L^2/3} to {@code L^2/12}, a
 * factor of 4. The axial term is untouched — it is centroid motion, which no
 * redistribution of mass within a rigid triangle can change — so
 * {@code C_spin -> L*L/24 = 1/12} and the swing becomes
 * {@code (2/3 + 1/12) / (1/12) = 9:1}. Measured, by area-averaging the very same
 * corner derivatives over each triangle, in
 * {@code FreeDynamicsTest.theSwingIsAPropertyOfTheMassModelNotOfTheLinkage}.
 * Never quote the 3:1 without "point masses at the 24 corners" attached.
 *
 * <p>
 * There is <b>no potential energy in this model</b> — none has ever been
 * specified — so the Lagrangian is {@code L = T} and the equation of motion is
 *
 * <pre>
 *     addot = -0.5 * (M_eff'(a) / M_eff(a)) * adot^2
 * </pre>
 *
 * <p>
 * <b>What this retracts.</b> A previously recorded and widely propagated claim
 * held that free motion conserves the momentum {@code p = M_eff * adot}. It does
 * not, and it cannot: {@code dL/da = 0.5 * M_eff'(a) * adot^2} is nonzero, so
 * {@code a} is not a cyclic coordinate. With {@code M_eff} varying 3:1, {@code p}
 * and {@code E} cannot both be constant, and the one that is constant is
 * {@code E}. Measured: {@code p} spreads by sqrt(3), {@code E} by 1.0000000,
 * {@code adot} by sqrt(3) — not by 3, as was also recorded. The correct relation
 * is {@code adot = sqrt(2E / M_eff)}, hence
 * {@code p = sqrt(2 * E * M_eff) ~ sqrt(M_eff)}.
 *
 * <p>
 * Angles here are in <b>radians</b>, unlike {@link JitterbugGeometry}, because
 * this is a dynamical equation rather than a geometric construction.
 *
 * @author halhildebrand
 */
public final class FreeDynamics {

    /**
     * A sampled solution. Arrays are parallel and of length {@code steps + 1}.
     */
    public record Trajectory(double[] time, double[] angle,
                             double[] angularVelocity) {

        /** {@code E = 0.5 * M_eff(a) * adot^2}, per unit M. */
        public double[] energies() {
            double[] out = new double[angle.length];
            for (int i = 0; i < out.length; i++) {
                out[i] = energy(angle[i], angularVelocity[i]);
            }
            return out;
        }

        /** {@code M_eff(a)}, per unit M. */
        public double[] effectiveMasses() {
            double[] out = new double[angle.length];
            for (int i = 0; i < out.length; i++) {
                out[i] = effectiveMass(angle[i]);
            }
            return out;
        }

        /** {@code p = M_eff(a) * adot}, per unit M. Not conserved. */
        public double[] momenta() {
            double[] out = new double[angle.length];
            for (int i = 0; i < out.length; i++) {
                out[i] = momentum(angle[i], angularVelocity[i]);
            }
            return out;
        }
    }

    /**
     * The axial-travel contribution to the effective mass.
     *
     * <p>
     * <b>Not</b> the square of {@link JitterbugGeometry#Z}, despite the earlier
     * name {@code Z_SQUARED} in this slot: {@code Z = 2R/sqrt(3)} so
     * {@code Z*Z = 4/3}, twice this. The factor of two is the {@code 48 = 2*24}
     * normalisation of {@code M_eff} described in the class javadoc. The exact
     * relation {@code AXIAL_COEFFICIENT == Z*Z/2} is asserted in
     * {@code FreeDynamicsTest.theCoefficientsAreTheGeometrysAndNotIndependentLiterals}
     * so the two constants can never drift apart unnoticed.
     */
    public static final double AXIAL_COEFFICIENT = 2.0 / 3.0;

    /**
     * The in-plane spin contribution to the effective mass, for point masses at
     * the corners: {@code L_EDGE*L_EDGE/6}. Uniform triangular laminae would put
     * {@code L_EDGE*L_EDGE/24 = 1/12} here instead — see the class javadoc.
     */
    public static final double SPIN_COEFFICIENT  = 1.0 / 3.0;

    private FreeDynamics() {
    }

    /**
     * {@code addot = -0.5 * (M_eff'/M_eff) * adot^2}. Euler-Lagrange for
     * {@code L = 0.5 M_eff(a) adot^2} with no potential.
     */
    public static double angularAcceleration(double aRad, double adot) {
        return -0.5 * effectiveMassDerivative(aRad) / effectiveMass(aRad) * adot
               * adot;
    }

    /**
     * {@code M_eff(a)/M = C_axial sin^2(a) + C_spin}. Ranges over [1/3, 1] for
     * the corner point-mass model this class documents.
     */
    public static double effectiveMass(double aRad) {
        double s = Math.sin(aRad);
        return AXIAL_COEFFICIENT * s * s + SPIN_COEFFICIENT;
    }

    /** {@code d/da} of {@link #effectiveMass(double)}. */
    public static double effectiveMassDerivative(double aRad) {
        return 2.0 * AXIAL_COEFFICIENT * Math.sin(aRad) * Math.cos(aRad);
    }

    /** {@code E = 0.5 * M_eff(a) * adot^2}. This is the conserved quantity. */
    public static double energy(double aRad, double adot) {
        return 0.5 * effectiveMass(aRad) * adot * adot;
    }

    /**
     * {@code p = M_eff(a) * adot}. Present because the retracted claim was about
     * it; it is <b>not</b> conserved.
     */
    public static double momentum(double aRad, double adot) {
        return effectiveMass(aRad) * adot;
    }

    /**
     * Fixed-step classical RK4 on {@code (a, adot)}. Deterministic: no adaptive
     * stepping, no wall clock, so a given call always returns bit-identical
     * output.
     *
     * @param steps number of steps; the returned arrays have {@code steps + 1}
     *              samples
     */
    public static Trajectory rk4(double a0, double adot0, double tEnd,
                                 int steps) {
        if (steps < 1) {
            throw new IllegalArgumentException("steps must be positive");
        }
        double h = tEnd / steps;
        double[] t = new double[steps + 1];
        double[] a = new double[steps + 1];
        double[] v = new double[steps + 1];
        a[0] = a0;
        v[0] = adot0;
        for (int i = 0; i < steps; i++) {
            double ai = a[i];
            double vi = v[i];

            double k1a = vi;
            double k1v = angularAcceleration(ai, vi);

            double k2a = vi + 0.5 * h * k1v;
            double k2v = angularAcceleration(ai + 0.5 * h * k1a,
                                             vi + 0.5 * h * k1v);

            double k3a = vi + 0.5 * h * k2v;
            double k3v = angularAcceleration(ai + 0.5 * h * k2a,
                                             vi + 0.5 * h * k2v);

            double k4a = vi + h * k3v;
            double k4v = angularAcceleration(ai + h * k3a, vi + h * k3v);

            a[i + 1] = ai + h / 6.0 * (k1a + 2 * k2a + 2 * k3a + k4a);
            v[i + 1] = vi + h / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v);
            t[i + 1] = (i + 1) * h;
        }
        return new Trajectory(t, a, v);
    }
}
