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

import com.chiralbehaviors.inviscid.jitterbug.JitterbugGeometry.Clustering;

/**
 * The jitterbug as a linkage rather than as a one-parameter family.
 *
 * <p>
 * Eight rigid triangles free in space carry {@code 8 x 6 = 48} degrees of
 * freedom. Twelve shared vertices, each identifying one corner of one triangle
 * with one corner of another, impose {@code 12 x 3 = 36} scalar constraints. Six
 * of the surviving directions are global rigid motions. What remains — measured,
 * not assumed — is <b>6 internal degrees of freedom</b>, and the symmetric
 * jitterbug path is one curve inside that variety.
 *
 * <p>
 * <b>Local coordinates.</b> Body {@code i} gets an infinitesimal rotation
 * {@code dw_i} about its own centroid plus a translation {@code dt_i}, so corner
 * {@code j} moves by {@code dw_i x r_ij + dt_i} with
 * {@code r_ij = X[i][j] - centroid_i}. The 48-vector packs all eight
 * {@code dw_i} into components 0..23 and all eight {@code dt_i} into 24..47.
 * {@link #jacobian(double[][][])} is the derivative of
 * {@link #residual(double[][][])} composed with
 * {@link #applyBodyMotions(double[][][], double[])} at zero motion, and the test
 * suite checks exactly that by finite differences.
 *
 * <p>
 * <b>The pairing is read once and held.</b> It is derived at a generic angle
 * where the configuration is 12 shared vertices of multiplicity 2. At
 * {@code a = 60, 120, 240, 300} those twelve merge into six of multiplicity 4 —
 * but the linkage does not re-wire itself, so re-deriving the pairing at a merge
 * angle would invent a different mechanism. The merge is a coincidence of point
 * positions; the extra coincidences are not constraints.
 *
 * <p>
 * No new dependencies: plain arithmetic only.
 *
 * @author halhildebrand
 */
public final class JitterbugLinkage {

    /**
     * Where the hinge combinatorics are read off. Any angle whose configuration
     * is 12x2 will do; 30 degrees is generic.
     */
    public static final double PAIRING_PROBE_DEG = 30.0;

    /**
     * The 12 hinges, each {@code {faceA, cornerA, faceB, cornerB}}.
     */
    private static final int[][] PAIRS;

    static {
        Clustering c = JitterbugGeometry.cluster(JitterbugGeometry
                                                                 .corners(PAIRING_PROBE_DEG),
                                                 1e-7);
        if (c.representatives().length != 12) {
            throw new IllegalStateException("hinge pairing must be read at a "
                                            + "generic angle; a="
                                            + PAIRING_PROBE_DEG + " gave "
                                            + c.representatives().length
                                            + " distinct corners");
        }
        for (int m : c.multiplicities()) {
            if (m != 2) {
                throw new IllegalStateException("expected 12 shared vertices of "
                                                + "multiplicity 2, saw " + m);
            }
        }
        PAIRS = new int[12][4];
        int[] filled = new int[12];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                int v = c.labels()[i][j];
                PAIRS[v][2 * filled[v]] = i;
                PAIRS[v][2 * filled[v] + 1] = j;
                filled[v]++;
            }
        }
    }

    private JitterbugLinkage() {
    }

    /**
     * Apply body motions {@code z} to base configuration {@code X0}, each body
     * rotating about its own centroid. Rotation vectors are exponentiated
     * (Rodrigues), so this is the exact nonlinear map whose derivative at
     * {@code z = 0} is {@link #jacobian(double[][][])}.
     */
    public static double[][][] applyBodyMotions(double[][][] x0, double[] z) {
        double[][] cen = JitterbugGeometry.faceCentroids(x0);
        double[][][] out = new double[8][3][3];
        for (int i = 0; i < 8; i++) {
            double[] w = { z[3 * i], z[3 * i + 1], z[3 * i + 2] };
            double theta = Math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
            double[][] r;
            if (theta < 1e-14) {
                r = new double[][] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            } else {
                r = JitterbugGeometry.rotation(new double[] { w[0] / theta,
                                                              w[1] / theta,
                                                              w[2] / theta },
                                               Math.toDegrees(theta));
            }
            for (int j = 0; j < 3; j++) {
                double[] d = { x0[i][j][0] - cen[i][0], x0[i][j][1] - cen[i][1],
                               x0[i][j][2] - cen[i][2] };
                for (int k = 0; k < 3; k++) {
                    out[i][j][k] = r[k][0] * d[0] + r[k][1] * d[1]
                                   + r[k][2] * d[2] + cen[i][k]
                                   + z[24 + 3 * i + k];
                }
            }
        }
        return out;
    }

    /**
     * The six global rigid motions as 48-vectors: three pure translations, then
     * three rotations about the origin. A global rotation {@code omega} spins
     * every body by {@code omega} <em>and</em> translates its centroid by
     * {@code omega x c_i} — omit that coupling and these stop being null
     * directions.
     */
    public static double[][] globalRigidMotions(double[][][] x) {
        double[][] cen = JitterbugGeometry.faceCentroids(x);
        double[][] g = new double[6][48];
        for (int d = 0; d < 3; d++) {
            for (int i = 0; i < 8; i++) {
                g[d][24 + 3 * i + d] = 1.0;
            }
        }
        for (int d = 0; d < 3; d++) {
            double[] e = new double[3];
            e[d] = 1.0;
            for (int i = 0; i < 8; i++) {
                g[3 + d][3 * i + d] = 1.0;
                double[] cross = { e[1] * cen[i][2] - e[2] * cen[i][1],
                                   e[2] * cen[i][0] - e[0] * cen[i][2],
                                   e[0] * cen[i][1] - e[1] * cen[i][0] };
                for (int k = 0; k < 3; k++) {
                    g[3 + d][24 + 3 * i + k] = cross[k];
                }
            }
        }
        return g;
    }

    /**
     * The 12 hinges, each {@code {faceA, cornerA, faceB, cornerB}}. Defensive
     * copy: this is the linkage's identity and must not drift.
     */
    public static int[][] hingePairs() {
        int[][] copy = new int[12][];
        for (int i = 0; i < 12; i++) {
            copy[i] = PAIRS[i].clone();
        }
        return copy;
    }

    /**
     * The 36x48 constraint Jacobian at corner configuration {@code X}.
     *
     * <p>
     * Row {@code 3c+d} is the {@code d}th component of hinge {@code c}'s
     * residual. Columns 0..23 are the angular blocks, 24..47 the translational
     * ones.
     */
    public static double[][] jacobian(double[][][] x) {
        double[][] cen = JitterbugGeometry.faceCentroids(x);
        double[][] j = new double[36][48];
        for (int c = 0; c < 12; c++) {
            int fa = PAIRS[c][0];
            int ca = PAIRS[c][1];
            int fb = PAIRS[c][2];
            int cb = PAIRS[c][3];
            double[] ra = { x[fa][ca][0] - cen[fa][0], x[fa][ca][1] - cen[fa][1],
                            x[fa][ca][2] - cen[fa][2] };
            double[] rb = { x[fb][cb][0] - cen[fb][0], x[fb][cb][1] - cen[fb][1],
                            x[fb][cb][2] - cen[fb][2] };
            for (int d = 0; d < 3; d++) {
                int row = 3 * c + d;
                double[] rowA = crossRow(ra, d);
                double[] rowB = crossRow(rb, d);
                for (int k = 0; k < 3; k++) {
                    j[row][3 * fa + k] += rowA[k];
                    j[row][3 * fb + k] -= rowB[k];
                }
                j[row][24 + 3 * fa + d] += 1.0;
                j[row][24 + 3 * fb + d] -= 1.0;
            }
        }
        return j;
    }

    /**
     * The 72x48 derivative of every corner position with respect to the body
     * motions, row-major over {@code [face][corner][xyz]}.
     *
     * <p>
     * {@link #jacobian(double[][][])} is the row-difference of this over the
     * hinge pairs, and the test suite checks that identity; the two are written
     * independently so that neither validates itself.
     */
    public static double[][] positionJacobian(double[][][] x) {
        double[][] cen = JitterbugGeometry.faceCentroids(x);
        double[][] p = new double[72][48];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                double[] r = { x[i][j][0] - cen[i][0], x[i][j][1] - cen[i][1],
                               x[i][j][2] - cen[i][2] };
                for (int d = 0; d < 3; d++) {
                    int row = 9 * i + 3 * j + d;
                    double[] cr = crossRow(r, d);
                    for (int k = 0; k < 3; k++) {
                        p[row][3 * i + k] = cr[k];
                    }
                    p[row][24 + 3 * i + d] = 1.0;
                }
            }
        }
        return p;
    }

    /**
     * The 36-component hinge residual: for each hinge, the vector from one
     * paired corner to the other. Zero exactly on the constraint variety.
     */
    public static double[] residual(double[][][] x) {
        double[] r = new double[36];
        for (int c = 0; c < 12; c++) {
            int fa = PAIRS[c][0];
            int ca = PAIRS[c][1];
            int fb = PAIRS[c][2];
            int cb = PAIRS[c][3];
            for (int d = 0; d < 3; d++) {
                r[3 * c + d] = x[fa][ca][d] - x[fb][cb][d];
            }
        }
        return r;
    }

    /**
     * Row {@code d} of the matrix {@code M} satisfying {@code M w == w x r},
     * i.e. of {@code -[r]_x}.
     */
    private static double[] crossRow(double[] r, int d) {
        switch (d) {
            case 0:
                return new double[] { 0, r[2], -r[1] };
            case 1:
                return new double[] { -r[2], 0, r[0] };
            default:
                return new double[] { r[1], -r[0], 0 };
        }
    }
}
