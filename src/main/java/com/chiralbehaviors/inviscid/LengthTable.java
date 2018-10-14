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

import static com.chiralbehaviors.inviscid.Constants.ROOT_2;
import static com.chiralbehaviors.inviscid.Constants.TWO_PI;

import javax.vecmath.Vector2d;

/**
 * @author halhildebrand
 *
 */
public class LengthTable {

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

    public final double[] lengths; 

    public LengthTable(int subdivisions) {
        if (subdivisions % 8 != 0) {
            throw new IllegalArgumentException("Subdivsions must be divisible by 8: "
                                               + subdivisions);
        }

        double halfEdge = 1.0 / ROOT_2;
        Vector2d center = new Vector2d();
        Vector2d left = new Vector2d(halfEdge, -halfEdge);
        Vector2d right = new Vector2d(halfEdge, 2 * halfEdge);
        lengths = new double[subdivisions / 8];
        float resolution = TWO_PI / subdivisions;
        double angle = 0;
        for (int i = 0; i < lengths.length; i++) {
            Vector2d pointing = new Vector2d(Math.cos(angle), Math.sin(angle));
            Vector2d intersection = intersect(center, pointing, left, right);
            lengths[i] = Math.min(1, intersection.length());
            angle += resolution;
        }
    } 

    public double lengthAt(int step) {
        if (step < lengths.length) {
            return lengths[step];
        }
        if (step < lengths.length * 2) {
            return lengths[lengths.length - (step - lengths.length) - 1];
        }
        if (step < lengths.length * 3) {
            return lengths[step - (lengths.length * 2)];
        }
        if (step < lengths.length * 4) {
            return lengths[lengths.length - (step - lengths.length * 3) - 1];
        }
        if (step < lengths.length * 5) {
            return lengths[step - (lengths.length * 4)];
        }
        if (step < lengths.length * 6) {
            return lengths[lengths.length - (step - lengths.length * 5) - 1];
        }
        if (step < lengths.length * 7) {
            return lengths[step - (lengths.length * 6)];
        }
        return lengths[lengths.length - (step - lengths.length * 7) - 1];
    }
}
