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

import static com.chiralbehaviors.inviscid.PhiCoordinates.p120v;

import javafx.geometry.Point3D;
import javafx.scene.shape.Cylinder;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Translate;

/**
 * @author halhildebrand
 *
 */
public class Edge {
    private final int[] endpoints;

    protected Edge(int[] endpoints) {
        assert endpoints != null && endpoints.length == 2;
        this.endpoints = endpoints;
    }

    public Cylinder createLine(double radius) {
        Point3D diff = p120v[endpoints[1]].subtract(p120v[endpoints[0]]);

        Point3D mid = p120v[endpoints[1]].midpoint(p120v[endpoints[0]]);
        Translate moveToMidpoint = new Translate(mid.getX(), mid.getY(),
                                                 mid.getZ());

        Point3D yAxis = new Point3D(0, 1, 0);
        Point3D axisOfRotation = diff.crossProduct(yAxis);
        double angle = Math.acos(diff.normalize()
                                     .dotProduct(yAxis));
        Rotate rotateAroundCenter = new Rotate(-Math.toDegrees(angle),
                                               axisOfRotation);

        Cylinder line = new Cylinder(radius, diff.magnitude());

        line.getTransforms()
            .addAll(moveToMidpoint, rotateAroundCenter);

        return line;
    }
}
