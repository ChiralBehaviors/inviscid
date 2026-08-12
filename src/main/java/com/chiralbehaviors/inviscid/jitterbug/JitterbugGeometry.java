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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * The symmetric jitterbug family, headless.
 *
 * <p>
 * Eight rigid equilateral triangles — the faces of an octahedron of circumradius
 * {@link #R_CIRC} — each spun about its own face axis and slid along it. The
 * corner positions are
 *
 * <pre>
 *     x_i(a) = R(u_i, sigma_i * (a - 60 deg)) * (v - c_i)  +  u_i * Z * cos(a)
 * </pre>
 *
 * with {@code c_i = s/3} the face centroid, {@code u_i = c_i/|c_i| = s/sqrt(3)}
 * the face axis, {@code sigma_i = sx*sy*sz} the face's handedness and
 * {@code Z = 2R/sqrt(3)}.
 *
 * <p>
 * <b>Where this comes from.</b> It is
 * {@code com.chiralbehaviors.inviscid.Jitterbug} with the JavaFX scene graph
 * removed, and it is <em>not</em> derivable from that class's live
 * {@code rotateTo}/{@code translate} alone. {@code Jitterbug.construct()} bakes
 * its mesh vertices at an angle of {@code -sigma*60} about the face axis with
 * the pivot at the centroid, applies a {@code -centroid} translation, and then
 * resets both transforms to identity before the animation ever runs. The baked
 * pre-rotation is why {@code a - 60} appears above. The pivot term cancels: the
 * pivot {@code c} is parallel to the rotation axis {@code u}, so
 * {@code R(u,t)*c == c} and the usual {@code c - R*c} contribution vanishes
 * identically.
 *
 * <p>
 * <b>Sign convention.</b> This class uses {@code sigma = +(sx*sy*sz)}, following
 * the Python reference harness in {@code analysis/jitterbug-variety/}. The
 * JavaFX side ({@code PhiCoordinates.JITTERBUG_INVERSES}) uses
 * {@code sigma = -(sx*sy*sz)}. The two are the same family traversed in opposite
 * senses — the negated-sigma configuration at {@code a} is this one at
 * {@code -a}, to machine precision — so no measured quantity changes. Do not
 * "fix" it; see
 * {@code JitterbugGeometryTest#negatedSigmaIsTheSameFamilyTraversedBackwards()},
 * which measures the agreement.
 *
 * <p>
 * No JavaFX, no scene graph, no new dependencies: plain arithmetic only.
 *
 * @author halhildebrand
 */
public final class JitterbugGeometry {

    /**
     * A grouping of near-identical points.
     *
     * @param representatives one point per cluster, in first-seen order
     * @param multiplicities  how many of the 24 corners landed in each cluster
     * @param labels          {@code [face][corner]} to cluster index
     */
    public record Clustering(double[][] representatives, int[] multiplicities,
                             int[][] labels) {
    }

    /** Octahedron circumradius. Everything is scaled to this. */
    public static final double     R_CIRC            = 1.0;

    /** Octahedron edge length, and therefore the rigid strut length. */
    public static final double     L_EDGE            = R_CIRC * Math.sqrt(2.0);

    /**
     * Axial travel amplitude, {@code 2R/sqrt(3)}. Matches the {@code Z} computed
     * in {@code Jitterbug}'s constructor from the octahedron edge length.
     */
    public static final double     Z                 = L_EDGE * Math.sqrt(2.0)
                                                       / Math.sqrt(3.0);

    /** The octahedron. Also, measurably, the interpenetration onset. */
    public static final double     A_OCTAHEDRON_DEG  = 60.0;

    /**
     * The icosahedron: the unique angle in (0,40) at which the open-square
     * diagonals reach strut length, giving 30 equal short edges. Located by
     * bisection in the test, not asserted from here.
     */
    public static final double     A_ICOSAHEDRON_DEG = 22.238756092965;

    /**
     * The eight sign triples, in the order
     * {@code itertools.product((1,-1), repeat=3)} produces them. Face index is
     * this index throughout the package, so the hinge pairings and the Jacobian
     * block layout depend on this order.
     */
    private static final double[][] SIGNS;

    static {
        SIGNS = new double[8][3];
        int n = 0;
        for (int x = 0; x < 2; x++) {
            for (int y = 0; y < 2; y++) {
                for (int z = 0; z < 2; z++) {
                    SIGNS[n][0] = x == 0 ? 1 : -1;
                    SIGNS[n][1] = y == 0 ? 1 : -1;
                    SIGNS[n][2] = z == 0 ? 1 : -1;
                    n++;
                }
            }
        }
    }

    private JitterbugGeometry() {
    }

    /**
     * Group corners that coincide to within {@code tol}. First-fit against the
     * representatives already seen, scanning corners in {@code [face][corner]}
     * row-major order — the same order as the reference implementation, so
     * cluster indices agree between the two.
     */
    public static Clustering cluster(double[][][] x, double tol) {
        List<double[]> reps = new ArrayList<>();
        int[][] labels = new int[8][3];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                int found = -1;
                for (int k = 0; k < reps.size(); k++) {
                    if (distance(x[i][j], reps.get(k)) < tol) {
                        found = k;
                        break;
                    }
                }
                if (found < 0) {
                    found = reps.size();
                    reps.add(x[i][j].clone());
                }
                labels[i][j] = found;
            }
        }
        int[] mult = new int[reps.size()];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                mult[labels[i][j]]++;
            }
        }
        return new Clustering(reps.toArray(new double[0][]), mult, labels);
    }

    /**
     * The 24 corner positions at angle {@code aDeg}, as
     * {@code [face][corner][xyz]}.
     */
    public static double[][][] corners(double aDeg) {
        return corners(aDeg, false, 1.0);
    }

    /**
     * @param uniformSigma when true, force sigma = +1 on every face. This is the
     *                     <em>wrong</em> choice geometrically — it breaks vertex
     *                     sharing — and exists only so that the test can
     *                     demonstrate chirality is forced rather than chosen.
     */
    public static double[][][] corners(double aDeg, boolean uniformSigma) {
        return corners(aDeg, uniformSigma, 1.0);
    }

    /**
     * @param signConvention {@code +1} for the harness convention
     *                       {@code sigma = +(sx*sy*sz)}, {@code -1} for the
     *                       JavaFX one. Provided so the equivalence of the two
     *                       can be measured; production callers want
     *                       {@code +1}.
     */
    public static double[][][] corners(double aDeg, boolean uniformSigma,
                                       double signConvention) {
        double[][][] out = new double[8][3][3];
        double axial = Z * Math.cos(Math.toRadians(aDeg));
        for (int i = 0; i < 8; i++) {
            double[] s = SIGNS[i];
            double[] c = { s[0] * R_CIRC / 3.0, s[1] * R_CIRC / 3.0,
                           s[2] * R_CIRC / 3.0 };
            double norm = Math.sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
            double[] u = { c[0] / norm, c[1] / norm, c[2] / norm };
            double sigma = signConvention
                           * (uniformSigma ? 1.0 : s[0] * s[1] * s[2]);
            double[][] rot = rotation(u, sigma * (aDeg - A_OCTAHEDRON_DEG));
            for (int j = 0; j < 3; j++) {
                double[] v = new double[3];
                v[j] = s[j] * R_CIRC;
                double[] d = { v[0] - c[0], v[1] - c[1], v[2] - c[2] };
                for (int k = 0; k < 3; k++) {
                    out[i][j][k] = rot[k][0] * d[0] + rot[k][1] * d[1]
                                   + rot[k][2] * d[2] + u[k] * axial;
                }
            }
        }
        return out;
    }

    /**
     * The unit face axis {@code u_i = s_i/sqrt(3)} about which face {@code i}
     * spins. Its coordinate signs are the face's sign triple, which is how face
     * indices map onto the vertices of the cube graph Q3.
     */
    public static double[] faceAxis(int face) {
        double[] s = SIGNS[face];
        double r = 1.0 / Math.sqrt(3.0);
        return new double[] { s[0] * r, s[1] * r, s[2] * r };
    }

    /** The centroid of each of the eight triangles, as {@code [face][xyz]}. */
    public static double[][] faceCentroids(double[][][] x) {
        double[][] cen = new double[8][3];
        for (int i = 0; i < 8; i++) {
            for (int k = 0; k < 3; k++) {
                cen[i][k] = (x[i][0][k] + x[i][1][k] + x[i][2][k]) / 3.0;
            }
        }
        return cen;
    }

    /**
     * Rodrigues rotation matrix about a unit axis, angle in degrees. Written as
     * {@code I + sin(t) K + (1-cos(t)) K^2} rather than in closed form so it is
     * term-for-term the reference implementation's.
     */
    public static double[][] rotation(double[] axis, double deg) {
        double t = Math.toRadians(deg);
        double sin = Math.sin(t);
        double omc = 1.0 - Math.cos(t);
        double[][] k = { { 0, -axis[2], axis[1] }, { axis[2], 0, -axis[0] },
                         { -axis[1], axis[0], 0 } };
        double[][] kk = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double acc = 0;
                for (int m = 0; m < 3; m++) {
                    acc += k[i][m] * k[m][j];
                }
                kk[i][j] = acc;
            }
        }
        double[][] r = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                r[i][j] = (i == j ? 1.0 : 0.0) + sin * k[i][j] + omc * kk[i][j];
            }
        }
        return r;
    }

    /**
     * The shortest distance between two shared vertices that are <em>not</em>
     * joined by a strut.
     *
     * <p>
     * Every edge of the cuboctahedron is a strut, so the nearest non-strut pair
     * is the short diagonal of a folding square. The icosahedron is the angle at
     * which that diagonal reaches strut length.
     */
    public static double shortestNonStrutDistance(double aDeg) {
        Clustering c = cluster(corners(aDeg), 1e-7);
        double[][] reps = c.representatives();
        Set<Long> struts = new HashSet<>();
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                int p = c.labels()[i][j];
                int q = c.labels()[i][(j + 1) % 3];
                struts.add(key(Math.min(p, q), Math.max(p, q)));
            }
        }
        double best = Double.MAX_VALUE;
        for (int p = 0; p < reps.length; p++) {
            for (int q = p + 1; q < reps.length; q++) {
                if (struts.contains(key(p, q))) {
                    continue;
                }
                best = Math.min(best, distance(reps[p], reps[q]));
            }
        }
        return best;
    }

    /**
     * All 24 strut lengths, three per triangle. Constant by construction — see
     * the guard comment on the corresponding test.
     */
    public static double[] strutLengths(double aDeg) {
        double[][][] x = corners(aDeg);
        double[] out = new double[24];
        int n = 0;
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                out[n++] = distance(x[i][j], x[i][(j + 1) % 3]);
            }
        }
        return out;
    }

    static double distance(double[] p, double[] q) {
        double dx = p[0] - q[0];
        double dy = p[1] - q[1];
        double dz = p[2] - q[2];
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    private static long key(int p, int q) {
        return ((long) p << 32) | q;
    }
}
