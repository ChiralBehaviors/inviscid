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

import static com.chiralbehaviors.inviscid.animations.Colors.blueMaterial;
import static com.chiralbehaviors.inviscid.animations.Colors.greenMaterial;
import static com.chiralbehaviors.inviscid.animations.Colors.redMaterial;

import java.util.function.BiFunction;
import java.util.function.Function;

import javafx.geometry.Point3D;
import javafx.scene.Group;
import javafx.scene.paint.Material;
import javafx.util.Pair;
import mesh.Line;

/**
 * @author halhildebrand
 *
 */
public class CubicGrid {
    private final double                 intervalX;
    private final double                 intervalY;
    private final double                 intervalZ;
    private final Point3D                origin;
    private final Point3D                xAxis;
    private final Pair<Integer, Integer> xExtent;
    private final Point3D                yAxis;
    private final Pair<Integer, Integer> yExtent;
    private final Point3D                zAxis;
    private final Pair<Integer, Integer> zExtent;

    public CubicGrid() {
        this(new Point3D(0, 0, 0), new Pair<>(5, 5), new Point3D(1, 0, 0), 1,
             new Pair<>(5, 5), new Point3D(0, 1, 0), 1, new Pair<>(5, 5),
             new Point3D(0, 0, 1), 1);
    }

    public CubicGrid(Point3D origin, Pair<Integer, Integer> xExtent,
                                 Point3D xAxis, double intervalX,
                                 Pair<Integer, Integer> yExtent, Point3D yAxis,
                                 double intervalY,
                                 Pair<Integer, Integer> zExtent, Point3D zAxis,
                                 double intervalZ) {
        this.origin = origin;
        this.xExtent = xExtent;
        this.xAxis = xAxis;
        this.intervalX = intervalX;
        this.yExtent = yExtent;
        this.yAxis = yAxis;
        this.intervalY = intervalY;
        this.zExtent = zExtent;
        this.zAxis = zAxis;
        this.intervalZ = intervalZ;
    }

    public Group construct(Material xaxis, Material yaxis, Material zaxis) {
        Group grid = new Group();
        Point3D pos;
        Point3D neg;

        neg = xAxis.normalize()
                   .multiply(-intervalX * xExtent.getKey())
                   .subtract(0, intervalY * yExtent.getKey(),
                             intervalZ * zExtent.getKey());
        pos = xAxis.normalize()
                   .multiply(intervalX * xExtent.getValue())
                   .subtract(0, intervalY * yExtent.getKey(),
                             intervalZ * zExtent.getKey());

        construct(grid, neg, pos, yExtent.getKey() + yExtent.getValue(),
                  zExtent.getKey() + zExtent.getValue(), xaxis,
                  (i, p) -> p.add(0, i * intervalY, 0),
                  p -> p.add(0, 0, intervalZ));

        neg = yAxis.normalize()
                   .multiply(-intervalY * yExtent.getKey())
                   .subtract(intervalX * xExtent.getKey(), 0,
                             intervalZ * zExtent.getKey());
        pos = yAxis.normalize()
                   .multiply(intervalY * yExtent.getValue())
                   .subtract(intervalX * xExtent.getKey(), 0,
                             intervalZ * zExtent.getKey());

        construct(grid, neg, pos, xExtent.getKey() + xExtent.getValue(),
                  zExtent.getKey() + zExtent.getValue(), yaxis,
                  (i, p) -> p.add(i * intervalX, 0, 0),
                  p -> p.add(0, 0, intervalZ));

        neg = zAxis.normalize()
                   .multiply(-intervalZ * zExtent.getKey())
                   .subtract(intervalX * xExtent.getKey(),
                             intervalY * yExtent.getKey(), 0);
        pos = zAxis.normalize()
                   .multiply(intervalZ * zExtent.getValue())
                   .subtract(intervalX * xExtent.getKey(),
                             intervalY * yExtent.getKey(), 0);

        construct(grid, neg, pos, xExtent.getKey() + xExtent.getValue(),
                  yExtent.getKey() + yExtent.getValue(), zaxis,
                  (i, p) -> p.add(i * intervalX, 0, 0),
                  p -> p.add(0, intervalY, 0));

        Line axis = new Line(0.025, xAxis.normalize()
                                         .multiply(-intervalX
                                                   * xExtent.getKey()),
                             xAxis.normalize()
                                  .multiply(intervalX * xExtent.getKey()));
        axis.setMaterial(redMaterial);
        grid.getChildren()
            .addAll(axis);

        axis = new Line(0.025, yAxis.normalize()
                                    .multiply(-intervalY * yExtent.getKey()),
                        yAxis.normalize()
                             .multiply(intervalY * yExtent.getKey()));
        axis.setMaterial(blueMaterial);
        grid.getChildren()
            .addAll(axis);

        axis = new Line(0.025, zAxis.normalize()
                                    .multiply(-intervalZ * zExtent.getKey()),
                        zAxis.normalize()
                             .multiply(intervalZ * zExtent.getKey()));
        axis.setMaterial(greenMaterial);
        grid.getChildren()
            .addAll(axis);
        return grid;
    }

    public double getIntervalX() {
        return intervalX;
    }

    public double getIntervalY() {
        return intervalY;
    }

    public double getIntervalZ() {
        return intervalZ;
    }

    public Point3D getOrigin() {
        return origin;
    }

    public Point3D getxAxis() {
        return xAxis;
    }

    public Pair<Integer, Integer> getxExtent() {
        return xExtent;
    }

    public Point3D getyAxis() {
        return yAxis;
    }

    public Pair<Integer, Integer> getyExtent() {
        return yExtent;
    }

    public Point3D getzAxis() {
        return zAxis;
    }

    public Pair<Integer, Integer> getzExtent() {
        return zExtent;
    }

    private void construct(Group grid, Point3D neg, Point3D pos, Integer a,
                           Integer b, Material material,
                           BiFunction<Integer, Point3D, Point3D> advanceA,
                           Function<Point3D, Point3D> advanceB) {
        Point3D start = neg;
        Point3D end = pos;
        Line axis = new Line(0.015, start, end);
        axis.setMaterial(material);
        grid.getChildren()
            .addAll(axis);
        for (int x = 0; x <= a; x++) {
            start = advanceA.apply(x, neg);
            end = advanceA.apply(x, pos);
            axis = new Line(0.015, start, end);
            axis.setMaterial(material);
            grid.getChildren()
                .addAll(axis);
            for (int z = 0; z < b; z++) {
                start = advanceB.apply(start);
                end = advanceB.apply(end);
                axis = new Line(0.015, start, end);
                axis.setMaterial(material);
                grid.getChildren()
                    .addAll(axis);
            }
        }
    }
}
