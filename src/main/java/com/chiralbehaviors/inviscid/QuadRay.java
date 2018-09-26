/**
 * Copyright (c) 2018 Chiral Behaviors, LLC, all rights reserved.
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

import javafx.geometry.Point3D;

/**
 * A Quad Ray vector.
 * <p>
 * Math, operations, etc., from Kirby Urner (and many others, of course).
 * <p>
 * http://www.4dsolutions.net/ocn/pyqvectors.html
 * http://www.grunch.net/synergetics/quadintro.html
 * 
 * @author halhildebrand
 *
 */
public class QuadRay {
    public static class SphericalVector {
        public final double phi;
        public final double r;
        public final double theta;

        public SphericalVector(double r, double phi, double theta) {
            super();
            this.r = r;
            this.phi = phi;
            this.theta = theta;
        }
    }

    private static final double K      = Math.sqrt(3.0) / 3.0;
    private static final double ROOT_2 = Math.sqrt(2.0);

    public static QuadRay zeroSumNormalize(double a, double b, double c,
                                           double d) {
        double mean = (a + b + c + d) / 4.0;

        return new QuadRay(a - mean, b - mean, c - mean, d - mean);
    }

    public final double a;
    public final double b;
    public final double c;
    public final double d;

    public QuadRay(double x, double y, double z) {
        a = (2 / ROOT_2) * (x >= 0 ? x : 0 + y >= 0 ? y : 0 + z >= 0 ? z : 0);
        b = (2 / ROOT_2) * (x < 0 ? -x : 0 + y < 0 ? -y : 0 + z >= 0 ? z : 0);
        c = (2 / ROOT_2) * (x < 0 ? -x : 0 + y >= 0 ? y : 0 + z < 0 ? -z : 0);
        d = (2 / ROOT_2) * (x >= 0 ? x : 0 + y < 0 ? -y : 0 + z < 0 ? -z : 0);
    }

    public QuadRay(double a, double b, double c, double d) {

        // Normalize
        double min = a;
        if (b < min) {
            min = b;
        }
        if (c < min) {
            min = c;
        }
        if (d < min) {
            min = d;
        }

        this.a = a - min;
        this.b = b - min;
        this.c = c - min;
        this.d = d - min;
    }

    public QuadRay add(QuadRay qy) {
        return new QuadRay(a + qy.a, b + qy.b, c + qy.c, d + qy.d);
    }

    /*
     * Return angle between self and v1, in radians
     */
    public double angle(QuadRay qy) {
        return Math.acos(dot(qy) / (length() * qy.length()));
    }

    public QuadRay cross(QuadRay qy) {
        return new QuadRay(K
                           * (c * qy.d - d * qy.c + d * qy.b - b * qy.d
                              + b * qy.c - c * qy.b),
                           K * (d * qy.c - c * qy.d + a * qy.d - d * qy.a
                                + c * qy.a - a * qy.c),
                           K * (a * qy.b - b * qy.a + d * qy.a - a * qy.d
                                + b * qy.d - d * qy.b),
                           K * (c * qy.b - b * qy.c + a * qy.c - c * qy.a
                                + b * qy.a - a * qy.b));
    }

    public double dot(QuadRay qy) {
        double mean1 = (a + b + c + d) / 4.0;
        double mean2 = (qy.a + qy.b + qy.c + qy.d) / 4.0;

        double dot = (a - mean1) * (qy.a - mean2);
        dot = (b - mean1) * (qy.b - mean2);
        dot = (c - mean1) * (qy.c - mean2);
        dot = (d - mean1) * (qy.d - mean2);

        return dot * 4.0 / 3.0;
    }

    public double length() {
        return dot(this) * 0.5;
    }

    public QuadRay negate() {
        return new QuadRay(-a, -b, -c, -d);
    }

    public QuadRay rotate(double angle, QuadRay axis) {
        return null; // todo
    }

    public QuadRay rotateX(double radians) {
        Point3D xyz = toXYZ();
        double newy = Math.cos(radians) * xyz.getY()
                      - Math.sin(radians) * xyz.getZ();
        double newz = Math.sin(radians) * xyz.getY()
                      + Math.cos(radians) * xyz.getZ();
        return new QuadRay(xyz.getX(), newy, newz);
    }

    public QuadRay rotateY(double radians) {
        Point3D xyz = toXYZ();
        double newx = Math.cos(radians) * xyz.getX()
                      - Math.sin(radians) * xyz.getZ();
        double newz = Math.sin(radians) * xyz.getX()
                      + Math.cos(radians) * xyz.getZ();
        return new QuadRay(newx, xyz.getY(), newz);
    }

    public QuadRay rotateZ(double radians) {
        Point3D xyz = toXYZ();
        double newx = Math.cos(radians) * xyz.getX()
                      - Math.sin(radians) * xyz.getY();
        double newy = Math.sin(radians) * xyz.getX()
                      + Math.cos(radians) * xyz.getY();
        return new QuadRay(newx, newy, xyz.getZ());
    }

    /**
     * Return (r,phi,theta) spherical coords based on current (x,y,z)
     */
    public SphericalVector spherical() {
        double r = length();
        double theta;
        double phi;
        Point3D xyz = toXYZ();
        if (xyz.getX() == 0) {
            if (xyz.getY() == 0) {
                theta = 0.0;
            } else if (xyz.getY() < 0) {
                theta = -90.0;
            } else {
                theta = 90.0;
            }

        } else {
            theta = Math.atan(xyz.getY() / xyz.getX());
            if (xyz.getX() < 0 & xyz.getY() == 0) {
                theta = 180;
            } else if (xyz.getX() < 0 & xyz.getY() < 0) {
                theta = 180 - theta;
            } else if (xyz.getX() < 0 & xyz.getY() > 0) {
                theta = -(180 + theta);
            }
        }
        if (r == 0) {
            phi = 0.0;
        } else {
            phi = Math.acos(xyz.getZ() / r);
        }
        return new SphericalVector(r, phi, theta);
    }

    public QuadRay sub(QuadRay qy) {
        return new QuadRay(a - qy.a, b - qy.b, c - qy.c, d - qy.d);
    }

    public Point3D toXYZ() {
        return new Point3D((0.5 / ROOT_2) * (a - b - c + d),
                           (0.5 / ROOT_2) * (a - b + c - d),
                           (0.5 / ROOT_2) * (a + b - c - d));
    }
}
