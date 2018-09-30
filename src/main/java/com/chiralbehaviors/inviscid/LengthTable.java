/**
 * Copyright (c) 2016 Chiral Behaviors, LLC, all rights reserved.
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

package com.chiralbehaviors.inviscid;

import javax.vecmath.Vector2d;

/**
 * @author halhildebrand
 *
 */
public class LengthTable {

    private static final double QUARTER_PI        = Math.PI / 4.0;
    private static final double HALF_PI           = Math.PI / 2.0;
    private static final double ROOT_2            = Math.sqrt(2.0);
    private static final double THREE_QUARTERS_PI = Math.PI * .75;
    private static final double TWO_PI            = 2 * Math.PI;

    private static Vector2d intersect(Vector2d a, Vector2d b, Vector2d c,
                                      Vector2d d) {
        // Line AB represented as a1x + b1y = c1 
        double a1 = b.y - a.y;
        double b1 = a.x - b.x;
        double c1 = a1 * a.x + b1 * a.y;

        // Line CD represented as a2x + b2y = c2 
        double a2 = d.y - c.y;
        double b2 = c.x - d.x;
        double c2 = a2 * c.x + b2 * c.y;

        double determinant = a1 * b2 - a2 * b1;

        if (determinant == 0) {
            return b;
        }

        double x = (b2 * c1 - b1 * c2) / determinant;
        double y = (a1 * c2 - a2 * c1) / determinant;
        return new Vector2d(x, y);
    }

    private final double[] lengths;

    private final double   resolution;

    public LengthTable(int subdivisions, double radius) {
        if (subdivisions % 8 != 0) {
            throw new IllegalArgumentException("Subdivsions must be divisible by 8: "
                                               + subdivisions);
        }

        Vector2d center = new Vector2d();
        double halfEdge = radius / ROOT_2;
        Vector2d left = new Vector2d(halfEdge, -halfEdge);
        Vector2d right = new Vector2d(halfEdge, 2 * halfEdge);
        lengths = new double[(subdivisions / 8) + 1];
        resolution = TWO_PI / subdivisions;
        double angle = 0;
        for (int i = 0; i < lengths.length; i++) {
            Vector2d pointing = new Vector2d(Math.cos(angle) * radius,
                                             Math.sin(angle) * radius);
            Vector2d intersection = intersect(center, pointing, left, right);
            lengths[i] = intersection == null ? radius
                                              : Math.min(radius,
                                                         intersection.length());
            angle += resolution;
        }
    }

    public double length(double angle) {
        assert angle <= TWO_PI;

        if (angle <= QUARTER_PI) {
            return lengths[(int) (angle / resolution)];
        }
        if (angle <= HALF_PI) {
            return lengths[lengths.length
                           - (int) ((angle - QUARTER_PI) / resolution)];
        }
        if (angle <= THREE_QUARTERS_PI) {
            return lengths[(int) ((angle - HALF_PI) / resolution)];
        }
        if (angle <= Math.PI) {
            return lengths[lengths.length
                           - (int) ((angle - THREE_QUARTERS_PI) / resolution)];
        }
        if (angle <= Math.PI + QUARTER_PI) {
            return lengths[(int) ((angle - Math.PI) / resolution)];
        }
        if (angle <= Math.PI + HALF_PI) {
            return lengths[lengths.length
                           - (int) ((angle - (Math.PI + QUARTER_PI))
                                    / resolution)];
        }
        if (angle <= Math.PI + THREE_QUARTERS_PI) {
            return lengths[(int) ((angle - (Math.PI + HALF_PI)) / resolution)];
        }
        return lengths[lengths.length
                       - (int) ((angle - (Math.PI + THREE_QUARTERS_PI))
                                / resolution)];
    }
}
