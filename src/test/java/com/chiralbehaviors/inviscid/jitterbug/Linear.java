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

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/**
 * Linear-algebra helpers the measurements need but production does not.
 *
 * <p>
 * commons-math3 is a <b>test-scope</b> dependency: {@code src/main} stays on
 * {@code javax.vecmath} and plain arithmetic. That boundary is deliberate —
 * rank, nullity and least squares are how the linkage is <em>measured</em>, not
 * how it is <em>defined</em>.
 *
 * @author halhildebrand
 */
final class Linear {

    private Linear() {
    }

    /**
     * Singular values in descending order. The 36x48 constraint Jacobian is
     * transposed first so the decomposition always sees rows &gt;= columns.
     */
    static double[] singularValues(double[][] m) {
        RealMatrix a = new Array2DRowRealMatrix(m, true);
        if (a.getRowDimension() < a.getColumnDimension()) {
            a = a.transpose();
        }
        return new SingularValueDecomposition(a).getSingularValues();
    }

    /**
     * Numerical rank, thresholding at {@code tol * max(1, sigma_0)} — the same
     * rule the Python reference uses, so the two agree on borderline cases.
     */
    static int rank(double[][] m, double tol) {
        double[] s = singularValues(m);
        double cut = tol * Math.max(1.0, s[0]);
        int r = 0;
        for (double v : s) {
            if (v > cut) {
                r++;
            }
        }
        return r;
    }

    /** Matrix-vector product. */
    static double[] apply(double[][] m, double[] v) {
        double[] out = new double[m.length];
        for (int i = 0; i < m.length; i++) {
            double acc = 0;
            for (int j = 0; j < v.length; j++) {
                acc += m[i][j] * v[j];
            }
            out[i] = acc;
        }
        return out;
    }

    static double norm(double[] v) {
        double acc = 0;
        for (double x : v) {
            acc += x * x;
        }
        return Math.sqrt(acc);
    }
}
