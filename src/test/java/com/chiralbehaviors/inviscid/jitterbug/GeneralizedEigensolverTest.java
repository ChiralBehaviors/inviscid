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

import java.util.Random;

import org.junit.Test;

/**
 * {@link GeneralizedEigensolver}: {@code H v = lambda M v} for symmetric
 * {@code H} and symmetric positive-definite {@code M}, via Cholesky-of-{@code M}.
 *
 * <p>
 * CHOICE-FREE SCAFFOLDING (bead inviscid-6dp): {@code H} throughout this suite
 * is either a hand-picked diagonal matrix or seeded random noise -- never a
 * potential's Hessian. No potential energy is chosen here.
 *
 * @author halhildebrand
 */
public class GeneralizedEigensolverTest {

    /**
     * Hand-computable case: {@code H = diag(8,3)}, {@code M = diag(4,1)}.
     * {@code det(H - lambda M) = (8-4 lambda)(3-lambda) = 0} gives
     * {@code lambda = 2} (eigenvector along e0) and {@code lambda = 3}
     * (eigenvector along e1); M-orthonormalising fixes the eigenvector
     * magnitudes to {@code (0.5, 0)} and {@code (0, 1)}.
     */
    @Test
    public void handComputedTwoByTwoCase() {
        double[][] h = { { 8, 0 }, { 0, 3 } };
        double[][] m = { { 4, 0 }, { 0, 1 } };
        GeneralizedEigensolver.Result r = GeneralizedEigensolver.solve(h, m);

        assertEquals(2, r.eigenvalues().length);
        assertEquals("ascending order, lambda_0", 2.0, r.eigenvalues()[0], 1e-9);
        assertEquals("ascending order, lambda_1", 3.0, r.eigenvalues()[1], 1e-9);

        assertEquals(0.5, Math.abs(r.eigenvectors()[0][0]), 1e-9);
        assertEquals(0.0, Math.abs(r.eigenvectors()[0][1]), 1e-9);
        assertEquals(0.0, Math.abs(r.eigenvectors()[1][0]), 1e-9);
        assertEquals(1.0, Math.abs(r.eigenvectors()[1][1]), 1e-9);
    }

    /**
     * {@code solve(M, M)} is a fixed point: every eigenvalue must be exactly
     * 1, for any positive-definite {@code M}, since {@code M v = 1 * M v}
     * trivially. A real check on the whole pipeline's wiring -- if the
     * Cholesky reduction or the eigenvector transform-back had a sign or
     * transpose bug, this would not land on exactly 1.
     */
    @Test
    public void identityCaseGivesEigenvalueOne() {
        double[][] m = { { 4, 1, 0 }, { 1, 3, 0.5 }, { 0, 0.5, 2 } };
        GeneralizedEigensolver.Result r = GeneralizedEigensolver.solve(m, m);
        for (double lambda : r.eigenvalues()) {
            assertEquals(1.0, lambda, 1e-9);
        }
    }

    /**
     * The returned eigenvectors must be M-orthonormal
     * ({@code v_i . (M v_j) = delta_ij}), the convention documented on
     * {@link GeneralizedEigensolver}, not Euclidean-orthonormal.
     */
    @Test
    public void eigenvectorsAreMOrthonormal() {
        double[][] h = { { 8, 0 }, { 0, 3 } };
        double[][] m = { { 4, 0 }, { 0, 1 } };
        GeneralizedEigensolver.Result r = GeneralizedEigensolver.solve(h, m);
        double[][] v = r.eigenvectors();
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                double acc = 0;
                for (int a = 0; a < 2; a++) {
                    for (int b = 0; b < 2; b++) {
                        acc += v[i][a] * m[a][b] * v[j][b];
                    }
                }
                assertEquals("v_" + i + " . M v_" + j, i == j ? 1.0 : 0.0, acc,
                             1e-9);
            }
        }
    }

    /**
     * On the real 6x6 internal mass metric (point-mass model, a generic
     * angle) with a seeded-random symmetric {@code H}, every returned
     * {@code (lambda_k, v_k)} pair must satisfy {@code H v_k = lambda_k M v_k}
     * to close numerical tolerance. This is the "arbitrary symmetric 6x6 H"
     * requirement, exercised against machinery {@link InternalMassMetricTest}
     * independently validates, not a synthetic toy.
     */
    @Test
    public void reconstructsHOnTheRealInternalMetricWithRandomH() {
        double[] cornerMass = new double[24];
        java.util.Arrays.fill(cornerMass, 1.0 / 48.0);
        double[][][] x = JitterbugGeometry.corners(37.0);
        double[][] m = InternalMassMetric.metric(x,
                                                 new InternalMassMetric.PointMasses(cornerMass));
        int n = m.length;
        assertEquals(6, n);

        Random rng = new Random(20260812L);
        double[][] h = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                double v = rng.nextDouble() * 2 - 1;
                h[i][j] = v;
                h[j][i] = v;
            }
        }

        GeneralizedEigensolver.Result r = GeneralizedEigensolver.solve(h, m);
        double worst = 0;
        for (int k = 0; k < n; k++) {
            double[] v = r.eigenvectors()[k];
            double lambda = r.eigenvalues()[k];
            for (int row = 0; row < n; row++) {
                double hv = 0;
                double mv = 0;
                for (int col = 0; col < n; col++) {
                    hv += h[row][col] * v[col];
                    mv += m[row][col] * v[col];
                }
                worst = Math.max(worst, Math.abs(hv - lambda * mv));
            }
        }
        assertTrue("H v = lambda M v, worst residual " + worst, worst < 1e-8);
    }

    /**
     * Eigenvalues must actually be sorted -- ascending is the documented
     * convention, and this uses a case where the natural (unsorted) order out
     * of the underlying symmetric eigensolver need not already be ascending,
     * so the sort is doing real work.
     */
    @Test
    public void eigenvaluesAreSortedAscending() {
        double[][] h = { { 1, 0, 0 }, { 0, 9, 0 }, { 0, 0, 4 } };
        double[][] m = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        GeneralizedEigensolver.Result r = GeneralizedEigensolver.solve(h, m);
        assertEquals(1.0, r.eigenvalues()[0], 1e-9);
        assertEquals(4.0, r.eigenvalues()[1], 1e-9);
        assertEquals(9.0, r.eigenvalues()[2], 1e-9);
    }

    /**
     * {@code M} with a genuinely negative eigenvalue (indefinite, not merely
     * singular) must be rejected outright: there is no metric there at all,
     * so the documented convention is to throw rather than return a nonsense
     * eigenvalue.
     */
    @Test(expected = IllegalArgumentException.class)
    public void throwsOnIndefiniteM() {
        double[][] h = { { 1, 0 }, { 0, 1 } };
        double[][] m = { { 1, 0 }, { 0, -1e-6 } };
        GeneralizedEigensolver.solve(h, m);
    }

    /**
     * {@code M} exactly singular (a genuine zero row) must also be rejected --
     * the boundary case of the indefinite check above, and the one most
     * likely to occur in practice at a rank-drop configuration.
     */
    @Test(expected = IllegalArgumentException.class)
    public void throwsOnExactlySingularM() {
        double[][] h = { { 1, 0 }, { 0, 1 } };
        double[][] m = { { 1, 0 }, { 0, 0 } };
        GeneralizedEigensolver.solve(h, m);
    }

    /**
     * {@code M} positive-definite but with a very small eigenvalue (well
     * above commons-math3's default absolute positivity threshold of
     * {@code 1e-10}) must still solve -- the documented convention is that
     * near-zero mass along a direction produces a correspondingly large
     * eigenvalue there, not a thrown exception. The small-mass direction's
     * eigenvalue must track {@code H_ii / M_ii} even though it is large.
     */
    @Test
    public void nearSingularMStillSolvesButAmplifies() {
        double small = 1e-6;
        double[][] h = { { 1, 0 }, { 0, 1 } };
        double[][] m = { { small, 0 }, { 0, 1 } };
        GeneralizedEigensolver.Result r = GeneralizedEigensolver.solve(h, m);
        assertEquals("large eigenvalue from the near-null mass direction",
                     1.0 / small, r.eigenvalues()[1], 1.0 / small * 1e-3);
        assertEquals("ordinary direction unaffected", 1.0, r.eigenvalues()[0],
                     1e-6);
    }

    /**
     * {@code solve} symmetrizes both {@code H} and {@code M} before reducing
     * to the standard eigenproblem -- both go through
     * {@code GeneralizedEigensolver#symmetrize} at the top of {@code solve},
     * so a genuinely asymmetric input must be treated as its symmetric part
     * {@code 0.5*(A + A^T)}, not as the raw matrix. Here {@code h[0][1] = 3}
     * while {@code h[1][0] = 1}, and {@code m[0][1] = 0.6} while
     * {@code m[1][0] = 0.2} -- not roundoff-level differences, so a returned
     * {@code (lambda_k, v_k)} pair genuinely distinguishes "symmetrized
     * first" from "used raw".
     *
     * <p>
     * This would fail if {@code symmetrize} were a no-op: with {@code M} left
     * genuinely asymmetric, commons-math3's {@code CholeskyDecomposition}
     * rejects non-symmetric input outright (its default relative-symmetry
     * threshold is far tighter than the 0.4 difference between
     * {@code m[0][1]} and {@code m[1][0]}), so {@code solve} would throw
     * rather than return a result. Passing the test requires both that no
     * exception is thrown and that the returned eigenpairs satisfy the
     * generalised eigenproblem for the hand-symmetrized matrices -- proving
     * {@code solve} actually used {@code 0.5*(H+H^T)} and
     * {@code 0.5*(M+M^T)}, not {@code H} and {@code M} as given.
     */
    @Test
    public void symmetrizesGenuinelyAsymmetricInputsBeforeSolving() {
        double[][] h = { { 4, 3 }, { 1, 2 } };
        double[][] m = { { 5, 0.6 }, { 0.2, 3 } };
        double[][] hSym = { { 4, 2 }, { 2, 2 } };
        double[][] mSym = { { 5, 0.4 }, { 0.4, 3 } };

        GeneralizedEigensolver.Result r = GeneralizedEigensolver.solve(h, m);
        int n = 2;
        double worst = 0;
        for (int k = 0; k < n; k++) {
            double[] v = r.eigenvectors()[k];
            double lambda = r.eigenvalues()[k];
            for (int row = 0; row < n; row++) {
                double hv = 0;
                double mv = 0;
                for (int col = 0; col < n; col++) {
                    hv += hSym[row][col] * v[col];
                    mv += mSym[row][col] * v[col];
                }
                worst = Math.max(worst, Math.abs(hv - lambda * mv));
            }
        }
        assertTrue("H_sym v = lambda M_sym v against the hand-symmetrized "
                   + "inputs, worst residual " + worst, worst < 1e-9);
    }
}
