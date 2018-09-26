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

package mesh;

import javafx.geometry.Point3D;
import javafx.scene.shape.Cylinder;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Translate;

/**
 * @author halhildebrand
 *
 */
public class Line extends Cylinder {

    public Line(double radius, Point3D a, Point3D b) {
        Point3D diff = b.subtract(a);

        Point3D mid = b.midpoint(a);
        Translate moveToMidpoint = new Translate(mid.getX(), mid.getY(),
                                                 mid.getZ());

        Point3D yAxis = new Point3D(0, 1, 0);
        Point3D axisOfRotation = diff.crossProduct(yAxis);
        double angle = Math.acos(diff.normalize()
                                     .dotProduct(yAxis));
        Rotate rotateAroundCenter = new Rotate(-Math.toDegrees(angle),
                                               axisOfRotation);

        setRadius(radius);
        setHeight(diff.magnitude());

        getTransforms().addAll(moveToMidpoint, rotateAroundCenter);
    }

}
