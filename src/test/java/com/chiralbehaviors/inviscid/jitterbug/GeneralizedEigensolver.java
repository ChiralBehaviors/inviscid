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

import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * A generalised symmetric eigensolve of {@code (H, M)}: {@code H v = lambda M v}
 * for a symmetric {@code H} and a symmetric positive-definite {@code M}, via
 * the standard Cholesky-of-{@code M} reduction to a standard symmetric
 * eigenproblem.
 *
 * <p>
 * <b>CHOICE-FREE SCAFFOLDING</b> (bead inviscid-6dp). {@code H} is supplied by
 * the caller -- in the intended use a potential's Hessian, restricted to the
 * same internal tangent space as {@link InternalMassMetric#metric}. This class
 * chooses no potential and no mass model: it accepts an arbitrary symmetric
 * {@code H} and turns it, together with a metric, into frequencies.
 *
 * <p>
 * <b>STATUS, corrected 2026-08-28 (bead inviscid-qvf.30).</b> This javadoc used
 * to say "once bead inviscid-qvf.2 chooses a potential". qvf.2 closed on
 * 2026-08-15 having chosen one -- DECISION 17, the raw all-pairs kernel, six
 * real frequencies at the vector equilibrium -- so that clause was stale for a
 * fortnight and read as though no potential existed. Two things are true
 * instead. Those frequencies were computed in Python
 * ({@code analysis/jitterbug-variety/attic/jb_s_frequency_spectrum.py} and the
 * jb_o / jb_t / jb_u family), never through this class, which still has no
 * caller outside its own tests. And DECISION 17's potential is a SINGLE-UNIT
 * one: the ARRAY's coupling is contact, whose effective potential is measured
 * in {@code jb_cp_contact_potential.py} to be an infinite square well -- flat
 * inside, infinite at the stop -- which has no Hessian anywhere and therefore
 * nothing this class can consume. That is a result about the medium, not a gap
 * in this scaffolding.
 *
 * <p>
 * <b>Method.</b> {@code M = L L^T} (Cholesky); {@code A = L^-1 H L^-T} is
 * symmetric ({@code H} symmetric implies {@code A} symmetric) and shares
 * {@code (H,M)}'s generalised eigenvalues, with generalised eigenvectors
 * recovered from {@code A}'s ordinary eigenvectors {@code y} as
 * {@code v = L^-T y}. This makes the returned eigenvectors
 * <b>M-orthonormal</b> ({@code v_i . (M v_j) = delta_ij}) rather than
 * Euclidean-orthonormal: {@code v_i . (M v_j) = y_i . (L^-1 M L^-T y_j) = y_i . y_j}
 * since {@code L^-1 M L^-T = L^-1 (L L^T) L^-T = I}, and the {@code y}'s are
 * Euclidean-orthonormal out of the underlying symmetric eigensolver. A caller
 * that Euclidean-renormalises these eigenvectors will get the wrong answer.
 *
 * <p>
 * <b>Convention when {@code M} is near-singular.</b> If {@code M} is not
 * numerically positive-definite -- commons-math3's
 * {@link CholeskyDecomposition} (default absolute positivity threshold
 * {@code 1e-10}) fails, which covers both a genuinely indefinite {@code M}
 * (e.g. evaluated somewhere it should never be positive-definite at all) and
 * an {@code M} whose smallest pivot has collapsed to numerical zero -- this
 * method throws {@link IllegalArgumentException}. There is no frequency along
 * a direction with no mass, and this solver refuses to fabricate one rather
 * than return an eigenvalue no one asked for. An {@code M} that is
 * positive-definite but merely <em>close</em> to singular (a small positive
 * eigenvalue, above the threshold) does not throw: Cholesky still succeeds,
 * but the reduction amplifies the corresponding eigenvalue without bound as
 * that eigenvalue of {@code M} shrinks -- the numerically correct limit
 * (near-zero mass in a direction implies a near-infinite frequency for any
 * finite stiffness), though increasingly unreliable in floating point as the
 * amplification grows. See {@code GeneralizedEigensolverTest} for both sides
 * of this boundary.
 *
 * @author halhildebrand
 */
final class GeneralizedEigensolver {

    /**
     * @param eigenvalues  ascending order
     * @param eigenvectors {@code eigenvectors[k]} is the eigenvector for
     *                     {@code eigenvalues[k]}; see the class javadoc for
     *                     the M-orthonormal convention.
     */
    record Result(double[] eigenvalues, double[][] eigenvectors) {
    }

    private GeneralizedEigensolver() {
    }

    static Result solve(double[][] h, double[][] m) {
        int n = h.length;
        if (m.length != n) {
            throw new IllegalArgumentException("H is " + n + "x" + n + " but M is "
                                                + m.length + "x" + m.length);
        }
        RealMatrix hm = new Array2DRowRealMatrix(symmetrize(h), true);
        RealMatrix mm = new Array2DRowRealMatrix(symmetrize(m), true);

        CholeskyDecomposition chol;
        try {
            chol = new CholeskyDecomposition(mm);
        } catch (NonPositiveDefiniteMatrixException e) {
            throw new IllegalArgumentException("M is not numerically positive-definite "
                                                + "(commons-math3 CholeskyDecomposition, "
                                                + "default absolute positivity threshold "
                                                + CholeskyDecomposition.DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD
                                                + ", rejected it); the generalised "
                                                + "eigenproblem is ill-posed here -- see the "
                                                + "class javadoc's near-singular-M convention",
                                                e);
        }
        RealMatrix l = chol.getL();
        RealMatrix lInv = new LUDecomposition(l).getSolver().getInverse();
        RealMatrix lInvT = lInv.transpose();
        RealMatrix a = lInv.multiply(hm).multiply(lInvT);
        // Symmetrise away roundoff before the symmetric eigensolver: A is
        // symmetric in exact arithmetic (H is symmetric) but floating-point
        // multiplication does not guarantee it bit-for-bit.
        RealMatrix aSym = a.add(a.transpose()).scalarMultiply(0.5);

        EigenDecomposition eig = new EigenDecomposition(aSym);

        Integer[] order = new Integer[n];
        for (int k = 0; k < n; k++) {
            order[k] = k;
        }
        Arrays.sort(order, (x, y) -> Double.compare(eig.getRealEigenvalue(x),
                                                    eig.getRealEigenvalue(y)));

        double[] eigenvalues = new double[n];
        double[][] eigenvectors = new double[n][];
        for (int k = 0; k < n; k++) {
            int src = order[k];
            eigenvalues[k] = eig.getRealEigenvalue(src);
            RealVector v = lInvT.operate(eig.getEigenvector(src));
            eigenvectors[k] = v.toArray();
        }
        return new Result(eigenvalues, eigenvectors);
    }

    private static double[][] symmetrize(double[][] a) {
        int n = a.length;
        double[][] out = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                out[i][j] = 0.5 * (a[i][j] + a[j][i]);
            }
        }
        return out;
    }
}
