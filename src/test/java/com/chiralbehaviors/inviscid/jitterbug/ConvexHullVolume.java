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
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Convex-hull volume of a small point set.
 *
 * <p>
 * The Python reference used {@code scipy.spatial.ConvexHull}; commons-math3 has
 * no 3D hull, so this replaces it. Brute-force supporting-plane enumeration
 * rather than an incremental hull: for a dozen points that is 220 triples and it
 * handles coplanar faces — the cuboctahedron's six squares — without any of the
 * degeneracy handling an incremental algorithm needs. It is O(n^4) and is
 * intended for the twelve-point configurations here, nothing larger.
 *
 * <p>
 * Test scope. Volume is a measurement of the geometry, not part of its
 * definition, and this implementation is only justified at this scale.
 *
 * @author halhildebrand
 */
final class ConvexHullVolume {

    private ConvexHullVolume() {
    }

    /**
     * Volume of the convex hull of {@code points}.
     *
     * @param tol coplanarity tolerance; points within {@code tol} of a
     *            supporting plane are taken to lie on it
     */
    static double volume(double[][] points, double tol) {
        int n = points.length;
        if (n < 4) {
            return 0;
        }
        double[] centre = new double[3];
        for (double[] p : points) {
            for (int k = 0; k < 3; k++) {
                centre[k] += p[k] / n;
            }
        }

        Set<String> seen = new HashSet<>();
        double total = 0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                for (int k = j + 1; k < n; k++) {
                    double[] e1 = sub(points[j], points[i]);
                    double[] e2 = sub(points[k], points[i]);
                    double[] normal = cross(e1, e2);
                    double len = norm(normal);
                    // Whether a triple is degenerate is a question about the
                    // ANGLE between its two edges, not about a length:
                    // |e1 x e2| = |e1||e2| sin(theta). Comparing the raw
                    // magnitude against tol -- which is a distance tolerance,
                    // used below as one -- mixes units, and would change verdict
                    // under a uniform rescaling of the point set. Scale-relative
                    // here; the two forms agree on this domain, where a
                    // cross-check against Qhull at 601 angles found no false
                    // facet under either.
                    if (len < tol * norm(e1) * norm(e2)) {
                        continue;
                    }
                    for (int c = 0; c < 3; c++) {
                        normal[c] /= len;
                    }
                    List<Integer> on = new ArrayList<>();
                    boolean allBelow = true;
                    boolean allAbove = true;
                    for (int m = 0; m < n; m++) {
                        double d = dot(normal, sub(points[m], points[i]));
                        if (Math.abs(d) <= tol) {
                            on.add(m);
                        } else if (d > 0) {
                            allBelow = false;
                        } else {
                            allAbove = false;
                        }
                    }
                    if (!allBelow && !allAbove) {
                        continue; // not a supporting plane
                    }
                    if (!seen.add(on.toString())) {
                        continue; // this facet already contributed
                    }
                    total += facetVolume(points, on, normal, centre);
                }
            }
        }
        return total;
    }

    /**
     * Sum of the tetrahedra spanned by the interior point {@code centre} and a
     * fan triangulation of one planar facet.
     */
    private static double facetVolume(double[][] points, List<Integer> on,
                                      double[] normal, double[] centre) {
        double[] fc = new double[3];
        for (int idx : on) {
            for (int k = 0; k < 3; k++) {
                fc[k] += points[idx][k] / on.size();
            }
        }
        double[] u = orthogonal(normal);
        double[] v = cross(normal, u);
        List<Integer> ordered = new ArrayList<>(on);
        ordered.sort(Comparator.comparingDouble(idx -> {
            double[] d = sub(points[idx], fc);
            return Math.atan2(dot(d, v), dot(d, u));
        }));

        double vol = 0;
        for (int t = 1; t + 1 < ordered.size(); t++) {
            double[] a = sub(points[ordered.get(0)], centre);
            double[] b = sub(points[ordered.get(t)], centre);
            double[] c = sub(points[ordered.get(t + 1)], centre);
            vol += Math.abs(dot(a, cross(b, c))) / 6.0;
        }
        return vol;
    }

    private static double[] cross(double[] a, double[] b) {
        return new double[] { a[1] * b[2] - a[2] * b[1],
                              a[2] * b[0] - a[0] * b[2],
                              a[0] * b[1] - a[1] * b[0] };
    }

    private static double dot(double[] a, double[] b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    private static double norm(double[] a) {
        return Math.sqrt(dot(a, a));
    }

    /** Any unit vector orthogonal to the unit vector {@code n}. */
    private static double[] orthogonal(double[] n) {
        double[] seed = Math.abs(n[0]) < 0.9 ? new double[] { 1, 0, 0 }
                                             : new double[] { 0, 1, 0 };
        double[] u = cross(n, seed);
        double len = norm(u);
        return new double[] { u[0] / len, u[1] / len, u[2] / len };
    }

    private static double[] sub(double[] a, double[] b) {
        return new double[] { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
    }
}
