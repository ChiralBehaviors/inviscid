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

import static com.chiralbehaviors.inviscid.jitterbug.FreeDynamics.AXIAL_COEFFICIENT;
import static com.chiralbehaviors.inviscid.jitterbug.FreeDynamics.SPIN_COEFFICIENT;
import static com.chiralbehaviors.inviscid.jitterbug.FreeDynamics.effectiveMass;
import static com.chiralbehaviors.inviscid.jitterbug.FreeDynamics.effectiveMassDerivative;
import static com.chiralbehaviors.inviscid.jitterbug.FreeDynamics.rk4;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * The free-dynamics retraction, pinned.
 *
 * <p>
 * A conservation law recorded as "computed and verified" said that free motion
 * conserves {@code p = M_eff * adot}. It is backwards. With no potential,
 * {@code dL/da = 0.5 M_eff'(a) adot^2} is nonzero, so {@code a} is not cyclic and
 * {@code p} is not conserved; {@code E} is. Since {@code M_eff} is not constant
 * the two claims are not even jointly possible, which is the part that should
 * have been caught before anything was built on it.
 *
 * <p>
 * <b>Mass model.</b> The 3:1 swing quoted throughout is for equal point masses
 * at the 24 triangle corners. That model is measured against
 * {@link JitterbugGeometry} here, not assumed, and the alternative reading —
 * uniform triangular laminae, which swings 9:1 — is measured alongside it. The
 * retraction above does not depend on which model is chosen: any non-constant
 * {@code M_eff} kills the momentum claim.
 *
 * <p>
 * Every number below is reproduced independently by
 * {@code analysis/jitterbug-variety/verify_critic.py} at rtol 1e-12 using an
 * adaptive integrator; this suite reaches the same values with fixed-step RK4,
 * which is a second, differently-conditioned confirmation rather than a rerun.
 *
 * @author halhildebrand
 */
public class FreeDynamicsTest {

    private static final double SQRT3 = Math.sqrt(3.0);

    /** The reference initial condition: the vector equilibrium, unit rate. */
    private static FreeDynamics.Trajectory reference() {
        return rk4(0.0, 1.0, 60.0, 300_000);
    }

    /**
     * Item 4, the headline. Over a long free run the momentum spreads by exactly
     * sqrt(3) while the energy does not move at all.
     */
    @Test
    public void energyIsConservedAndMomentumIsNot() {
        FreeDynamics.Trajectory tr = reference();

        double pSpread = spread(tr.momenta());
        double eSpread = spread(tr.energies());

        assertEquals("p = M_eff*adot spreads by sqrt(3): NOT CONSERVED", SQRT3,
                     pSpread, 1e-7);
        assertTrue("p must be measurably non-constant", pSpread > 1.7);
        assertEquals("E is CONSERVED", 1.0, eSpread, 1e-9);

        // The absolute values, not just the ratio. E = 0.5*(1/3)*1 = 1/6, and
        // p runs from sqrt(2*E*M_min)=1/3 to sqrt(2*E*M_max)=1/sqrt(3).
        assertEquals(1.0 / 6.0, min(tr.energies()), 1e-12);
        assertEquals(1.0 / 3.0, min(tr.momenta()), 1e-8);
        assertEquals(1.0 / SQRT3, max(tr.momenta()), 1e-8);
    }

    /**
     * Item 4, the second retraction. {@code adot} spreads by sqrt(3), not by 3.
     * The 3 was the {@code M_eff} spread misattributed to the rate.
     */
    @Test
    public void angularRateSpreadsBySqrtThreeNotByThree() {
        FreeDynamics.Trajectory tr = reference();
        double v = spread(tr.angularVelocity());
        assertEquals("adot spread", SQRT3, v, 1e-7);
        assertTrue("adot spread is emphatically not 3", Math.abs(v - 3.0) > 1.2);
        assertEquals("adot max", 1.0, max(tr.angularVelocity()), 1e-9);
        assertEquals("adot min", 1.0 / SQRT3, min(tr.angularVelocity()), 1e-8);
    }

    /**
     * Item 4. The 3:1 effective-mass swing survives the retraction — it is the
     * part of the old record that was right, <b>for point masses at the 24
     * corners</b>. {@link #theSwingIsAPropertyOfTheMassModelNotOfTheLinkage()}
     * shows what happens to it under the other natural reading of "eight rigid
     * triangles".
     */
    @Test
    public void effectiveMassSwingsThreeToOne() {
        FreeDynamics.Trajectory tr = reference();
        assertEquals("M_eff spread along the trajectory", 3.0,
                     spread(tr.effectiveMasses()), 1e-6);

        // Analytically, over the whole range rather than what the trajectory
        // happened to sample.
        double lo = Double.MAX_VALUE;
        double hi = 0;
        for (int k = 0; k <= 360_000; k++) {
            double m = effectiveMass(Math.toRadians(k / 1000.0));
            lo = Math.min(lo, m);
            hi = Math.max(hi, m);
        }
        assertEquals(SPIN_COEFFICIENT, lo, 1e-12);
        assertEquals(AXIAL_COEFFICIENT + SPIN_COEFFICIENT, hi, 1e-12);
        assertEquals(1.0, hi, 1e-15);
        assertEquals(3.0, hi / lo, 1e-9);
    }

    /**
     * Item 4. {@code adot} peaks at the vector equilibrium, {@code a = 0 mod pi},
     * where the effective mass is least.
     *
     * <p>
     * Started away from the peak so the maximum is genuinely located by the
     * dynamics rather than handed over as the initial condition.
     */
    @Test
    public void angularRatePeaksAtTheVectorEquilibrium() {
        FreeDynamics.Trajectory tr = rk4(1.0, 1.0, 40.0, 200_000);
        int best = 0;
        for (int i = 1; i < tr.angularVelocity().length; i++) {
            if (tr.angularVelocity()[i] > tr.angularVelocity()[best]) {
                best = i;
            }
        }
        assertTrue("the peak must not be the initial sample", best > 0);
        assertEquals("adot peaks where sin(a) = 0, i.e. at the VE", 0.0,
                     Math.abs(Math.sin(tr.angle()[best])), 1e-3);
    }

    /**
     * The equation of motion must actually be the Euler-Lagrange equation of
     * {@code L = 0.5 M_eff adot^2}: {@code d/dt(M_eff adot) = 0.5 M_eff' adot^2}.
     * Checked by finite-differencing {@code M_eff'} out of {@code M_eff}, so a
     * transcription error in either one is caught.
     */
    @Test
    public void effectiveMassDerivativeMatchesFiniteDifferences() {
        double h = 1e-6;
        double worst = 0;
        for (int k = 0; k < 1000; k++) {
            double a = k * 2 * Math.PI / 1000.0;
            double fd = (effectiveMass(a + h) - effectiveMass(a - h)) / (2 * h);
            worst = Math.max(worst,
                             Math.abs(effectiveMassDerivative(a) - fd));
        }
        assertTrue("M_eff' vs finite differences, worst " + worst,
                   worst < 1e-9);
    }

    /**
     * The retracted claim was not merely unsupported, it was self-contradictory:
     * conserving {@code p} and conserving {@code E} simultaneously forces
     * {@code M_eff} to be constant, and it is not.
     *
     * <p>
     * The substance of this test is the last line. The identity above it is
     * labelled as a guard, not as a result.
     */
    @Test
    public void bothConservationLawsCannotHoldAtOnce() {
        FreeDynamics.Trajectory tr = reference();
        double[] p = tr.momenta();
        double[] e = tr.energies();
        double[] m = tr.effectiveMasses();
        for (int i = 0; i < p.length; i += 5000) {
            // GUARD, not a verification: p^2 = 2*E*M_eff is an algebraic
            // identity of momentum(), energy() and effectiveMass() as written,
            // so this cannot fail while those three agree with each other. It
            // earns its place only as the premise of the real claim below —
            // constant p and constant E would demand constant M_eff.
            assertEquals(p[i] * p[i], 2 * e[i] * m[i], 1e-12);
        }
        assertTrue("M_eff is not constant", spread(m) > 2.9);
    }

    /**
     * {@code M_eff} must be the corner parameterisation's kinetic form and not
     * an independent literal. {@link FreeDynamics} never imports
     * {@link JitterbugGeometry}, so without this test a divergence between the
     * dynamics and the geometry could not be caught by anything.
     *
     * <p>
     * For equal point masses at the 24 corners the kinetic energy is
     * {@code 0.5 * m * sum|dx_ij/da|^2 * adot^2}, so the shape of {@code M_eff}
     * is the shape of that sum. Measured by central-differencing
     * {@link JitterbugGeometry#corners(double)} itself, which is the whole
     * point — nothing here is recomputed from the formula under test.
     */
    @Test
    public void effectiveMassIsTheKineticFormOfTheCornerParameterisation() {
        double worst = 0;
        for (int k = 0; k <= 360; k++) {
            double aDeg = k;
            double aRad = Math.toRadians(aDeg);
            double sum = 0;
            for (double[][] face : cornerDerivatives(aDeg)) {
                for (double[] d : face) {
                    sum += dot(d, d);
                }
            }
            assertEquals("closed form of the sum at a=" + aDeg,
                         32 * Math.sin(aRad) * Math.sin(aRad) + 16, sum, 1e-6);
            worst = Math.max(worst,
                             Math.abs(sum - 48 * effectiveMass(aRad)));
        }
        assertTrue("sum|dx/da|^2 must equal 48*M_eff, worst deviation " + worst,
                   worst < 1e-6);
    }

    /**
     * The coefficients are the geometry's, term by term, and this pins the
     * arithmetic that relates them so the two classes cannot drift apart.
     *
     * <p>
     * Also disambiguates a genuine name trap: the constant now called
     * {@link FreeDynamics#AXIAL_COEFFICIENT} used to be called
     * {@code Z_SQUARED} while {@link JitterbugGeometry#Z} squares to twice it.
     * Nothing compared them, so nothing would have noticed.
     *
     * <p>
     * <b>The coefficients are read out of {@link FreeDynamics#effectiveMass},
     * not off the named constants, and that is deliberate.</b> Both are
     * {@code static final double}, so the compiler inlines their values into
     * this class at <em>its</em> compile time; under incremental compilation a
     * changed {@code FreeDynamics} and an unrecompiled test can disagree while
     * every constant-to-constant assertion still passes. Verified, not
     * theorised: mutating the axial coefficient and rebuilding only the main
     * class left the constant form of this test green, and it failed only after
     * a clean rebuild. A method call cannot be inlined that way, so the
     * function form holds under either.
     */
    @Test
    public void theCoefficientsAreTheGeometrysAndNotIndependentLiterals() {
        double z = JitterbugGeometry.Z;
        double l = JitterbugGeometry.L_EDGE;

        // M_eff(0) = C_spin and M_eff(pi/2) = C_axial + C_spin, by evaluation.
        double spin = effectiveMass(0.0);
        double axial = effectiveMass(Math.PI / 2) - spin;

        assertEquals("axial coefficient is Z*Z/2, NOT Z*Z", z * z / 2, axial,
                     1e-15);
        assertTrue("and Z*Z is emphatically not it",
                   Math.abs(z * z - axial) > 0.6);

        // Corners of an equilateral triangle of side L lie L/sqrt(3) from the
        // centroid, so their polar second moment per unit mass is L^2/3.
        assertEquals("spin coefficient is L*L/6", l * l / 6, spin, 1e-15);

        // The published constants must agree with the function they name.
        assertEquals("AXIAL_COEFFICIENT names what effectiveMass uses",
                     AXIAL_COEFFICIENT, axial, 1e-15);
        assertEquals("SPIN_COEFFICIENT names what effectiveMass uses",
                     SPIN_COEFFICIENT, spin, 1e-15);
    }

    /**
     * The 3:1 swing is contingent on the mass model, and this measures the
     * alternative instead of asserting it away.
     *
     * <p>
     * Swap the corner point masses for <b>uniform triangular laminae</b> and the
     * swing is <b>9:1</b>, not 3:1. A point with fixed barycentric coordinates
     * moves affinely with the corners under a rigid motion, so
     * {@code |dp/da|^2} is quadratic in the barycentrics and the three edge
     * midpoints integrate it exactly over the triangle — degree-2 exact
     * quadrature, so this is a measurement and not a discretisation. The area
     * average obeys {@code 32 sin^2(a) + 4}: the axial term is unchanged, being
     * centroid motion, and the spin term drops by exactly the 4 that separates
     * a triangle's {@code L^2/12} from its corners' {@code L^2/3}.
     */
    @Test
    public void theSwingIsAPropertyOfTheMassModelNotOfTheLinkage() {
        double lo = Double.MAX_VALUE;
        double hi = 0;
        for (int k = 0; k <= 360; k++) {
            double aRad = Math.toRadians(k);
            double lamina = 0;
            for (double[][] d : cornerDerivatives(k)) {
                for (int j = 0; j < 3; j++) {
                    double[] mid = { 0.5 * (d[j][0] + d[(j + 1) % 3][0]),
                                     0.5 * (d[j][1] + d[(j + 1) % 3][1]),
                                     0.5 * (d[j][2] + d[(j + 1) % 3][2]) };
                    lamina += dot(mid, mid);
                }
            }
            assertEquals("area-averaged |dp/da|^2 at a=" + k,
                         32 * Math.sin(aRad) * Math.sin(aRad) + 4, lamina, 1e-6);
            lo = Math.min(lo, lamina);
            hi = Math.max(hi, lamina);
        }
        assertEquals("uniform laminae swing 9:1", 9.0, hi / lo, 1e-6);
        assertTrue("which is emphatically not the corner point-mass 3:1",
                   Math.abs(hi / lo - 3.0) > 5.0);
    }

    /**
     * {@code d(corner)/da} in radians, {@code [face][corner][xyz]}, by central
     * difference on {@link JitterbugGeometry#corners(double)}. The step is in
     * degrees because that is the parameterisation's unit; the derivative is
     * converted to radians to match {@link FreeDynamics}.
     */
    private static double[][][] cornerDerivatives(double aDeg) {
        double hDeg = 1e-3;
        double[][][] plus = JitterbugGeometry.corners(aDeg + hDeg);
        double[][][] minus = JitterbugGeometry.corners(aDeg - hDeg);
        double den = 2 * Math.toRadians(hDeg);
        double[][][] d = new double[8][3][3];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                for (int c = 0; c < 3; c++) {
                    d[i][j][c] = (plus[i][j][c] - minus[i][j][c]) / den;
                }
            }
        }
        return d;
    }

    private static double dot(double[] a, double[] b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    /**
     * The circuit-period consequence of the retraction, at the matched initial
     * condition a=0, adot=1. The recorded period was the p-conserved value,
     * 4*pi = 12.5663706; the true, E-conserved value is 8.7377526. The recorded
     * figure is wrong by a factor of 1.4381697.
     *
     * <p>
     * The period values are normalisation-dependent, so the ratio is the durable
     * number here, not the two absolutes.
     */
    @Test
    public void theRecordedCircuitPeriodIsWrongByAFactorOfOnePointFour() {
        int n = 200_000;
        double f0 = effectiveMass(0.0);
        double ifP = 0;
        double ifE = 0;
        for (int k = 0; k <= n; k++) {
            double a = 2 * Math.PI * k / n;
            double w = (k == 0 || k == n) ? 0.5 : 1.0;
            ifP += w * effectiveMass(a) / f0;
            ifE += w * Math.sqrt(effectiveMass(a) / f0);
        }
        double da = 2 * Math.PI / n;
        double periodIfMomentumConserved = ifP * da;
        double periodTrue = ifE * da;

        assertEquals("p-conserved period is 4*pi", 4 * Math.PI,
                     periodIfMomentumConserved, 1e-9);
        assertEquals("p-conserved period", 12.5663706,
                     periodIfMomentumConserved, 1e-7);
        assertEquals("true, E-conserved period", 8.7377526, periodTrue, 1e-7);
        assertEquals("the recorded period was too long by this factor",
                     1.4381697, periodIfMomentumConserved / periodTrue, 1e-7);
    }

    /**
     * The integrator must be deterministic and its accuracy must come from the
     * step count, not from luck. Halving the step size has to leave the answers
     * unchanged at the tolerances the tests above use; if the reported spreads
     * were integration artefacts, they would move.
     */
    @Test
    public void resultsAreStableUnderStepRefinementAndDeterministic() {
        FreeDynamics.Trajectory coarse = rk4(0.0, 1.0, 60.0, 150_000);
        FreeDynamics.Trajectory fine = rk4(0.0, 1.0, 60.0, 600_000);
        assertEquals(spread(coarse.momenta()), spread(fine.momenta()), 1e-9);
        assertEquals(spread(coarse.energies()), spread(fine.energies()), 1e-9);
        assertEquals(spread(coarse.angularVelocity()),
                     spread(fine.angularVelocity()), 1e-9);

        FreeDynamics.Trajectory again = rk4(0.0, 1.0, 60.0, 150_000);
        for (int i = 0; i < coarse.angle().length; i += 977) {
            assertEquals(coarse.angle()[i], again.angle()[i], 0.0);
            assertEquals(coarse.angularVelocity()[i],
                         again.angularVelocity()[i], 0.0);
        }
    }

    private static double max(double[] v) {
        double m = -Double.MAX_VALUE;
        for (double x : v) {
            m = Math.max(m, x);
        }
        return m;
    }

    private static double min(double[] v) {
        double m = Double.MAX_VALUE;
        for (double x : v) {
            m = Math.min(m, x);
        }
        return m;
    }

    private static double spread(double[] v) {
        return max(v) / min(v);
    }
}
