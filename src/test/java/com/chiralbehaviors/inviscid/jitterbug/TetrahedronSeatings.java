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
import java.util.List;

/**
 * Every way to send each triangle's three corners to three distinct tetrahedron
 * vertices such that all twelve hinge pairings agree.
 *
 * <p>
 * An equilateral triangle of edge L maps rigidly onto a tetrahedron face of edge
 * L under any of the six bijections — flipping it over is a rigid motion in 3D —
 * so combinatorial consistency <em>is</em> realisability here, and the count is
 * a count of genuine configurations rather than of labellings.
 *
 * <p>
 * Test scope: this is a measurement of the linkage, not part of it.
 * Enumeration order matches the Python reference
 * ({@code itertools.permutations(range(4), 3)}, faces recursed 0..7) so seating
 * indices are comparable between the two implementations.
 *
 * @author halhildebrand
 */
final class TetrahedronSeatings {

    /** A regular tetrahedron of edge sqrt(2) — the jitterbug's strut length. */
    static final double[][] TET = { { 0.5, 0.5, 0.5 }, { 0.5, -0.5, -0.5 },
                                    { -0.5, 0.5, -0.5 }, { -0.5, -0.5, 0.5 } };

    /** The 24 injections {0,1,2} -> {0,1,2,3}, in lexicographic order. */
    private static final int[][] OPTIONS;

    static {
        List<int[]> opts = new ArrayList<>();
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                for (int c = 0; c < 4; c++) {
                    if (a != b && b != c && a != c) {
                        opts.add(new int[] { a, b, c });
                    }
                }
            }
        }
        OPTIONS = opts.toArray(new int[0][]);
    }

    private TetrahedronSeatings() {
    }

    /**
     * Realise a seating as an 8x3x3 corner configuration.
     */
    static double[][][] build(int[][] seating) {
        double[][][] x = new double[8][3][3];
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                x[i][j] = TET[seating[i][j]].clone();
            }
        }
        return x;
    }

    /**
     * All hinge-consistent seatings, with no occupancy rule applied. Each is an
     * {@code int[8][3]} mapping (face, corner) to a tetrahedron vertex index.
     */
    static List<int[][]> all() {
        List<int[][]> out = new ArrayList<>();
        int[][] pairs = JitterbugLinkage.hingePairs();
        int[][] assign = new int[8][];
        recurse(0, assign, pairs, out);
        return out;
    }

    /**
     * The sorted multiset of tetrahedron-face occupancies for a seating, e.g.
     * {@code [2,2,2,2]} for the tetrahedron proper.
     */
    static int[] occupancyProfile(int[][] seating) {
        List<Integer> faces = new ArrayList<>();
        for (int i = 0; i < 8; i++) {
            int mask = 0;
            for (int j = 0; j < 3; j++) {
                mask |= 1 << seating[i][j];
            }
            faces.add(mask);
        }
        List<Integer> distinct = faces.stream().distinct().toList();
        int[] counts = new int[distinct.size()];
        for (int k = 0; k < distinct.size(); k++) {
            int f = distinct.get(k);
            counts[k] = (int) faces.stream().filter(v -> v == f).count();
        }
        java.util.Arrays.sort(counts);
        for (int i = 0, j = counts.length - 1; i < j; i++, j--) {
            int t = counts[i];
            counts[i] = counts[j];
            counts[j] = t;
        }
        return counts;
    }

    /**
     * The seatings that put exactly two triangles on each of the four
     * tetrahedron faces — the tetrahedron configuration proper.
     */
    static List<int[][]> tetrahedra() {
        List<int[][]> out = new ArrayList<>();
        for (int[][] s : all()) {
            int[] p = occupancyProfile(s);
            if (p.length == 4 && p[0] == 2 && p[3] == 2) {
                out.add(s);
            }
        }
        return out;
    }

    private static boolean consistent(int face, int[] option, int[][] assign,
                                      int[][] pairs) {
        for (int[] p : pairs) {
            if (p[0] == face && assign[p[2]] != null
                && option[p[1]] != assign[p[2]][p[3]]) {
                return false;
            }
            if (p[2] == face && assign[p[0]] != null
                && option[p[3]] != assign[p[0]][p[1]]) {
                return false;
            }
        }
        return true;
    }

    private static void recurse(int face, int[][] assign, int[][] pairs,
                                List<int[][]> out) {
        if (face == 8) {
            int[][] copy = new int[8][];
            for (int i = 0; i < 8; i++) {
                copy[i] = assign[i].clone();
            }
            out.add(copy);
            return;
        }
        for (int[] option : OPTIONS) {
            if (consistent(face, option, assign, pairs)) {
                assign[face] = option;
                recurse(face + 1, assign, pairs, out);
                assign[face] = null;
            }
        }
    }
}
