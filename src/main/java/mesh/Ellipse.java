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

package mesh;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Vector3d;

import javafx.geometry.Point3D;
import javafx.scene.paint.Material;
import mesh.polyhedra.plato.Octahedron;

/**
 * @author halhildebrand
 *
 */
public class Ellipse {
    private static final double TWO_PI = 2 * Math.PI;

    private final Point3D       center;
    private final Point3D       u;
    private final Point3D       v;

    public Ellipse(int f, Octahedron oct, int vt) {
        Face face = oct.getFaces()
                       .get(f);
        Vector3d c = face.centroid();
        Vector3d vx = face.getVertices()
                          .get(vt);

        Point3D centroid = new Point3D(c.x, c.y, c.z);
        Point3D vertex = new Point3D(vx.x, vx.y, vx.z);
        center = new Point3D(0, 0, 0);
        u = vertex.subtract(center);
        v = vertex.subtract(centroid)
                  .crossProduct(centroid);
    }

    public Ellipse(Point3D center, Point3D a, Point3D b) {
        u = a.subtract(center);
        v = b.subtract(center);
        this.center = center;
    }

    public PolyLine construct(int segments, Material material, double radius) {
        List<Point3D> points = new ArrayList<>();
        double increment = TWO_PI / segments;
        double theta = 0;
        for (int i = 0; i <= segments; i++) {
            points.add(u.multiply(Math.cos(theta))
                        .add(v.multiply(Math.sin(theta))));
            theta += increment;
        }
        return new PolyLine(points, radius, material);
    }

    public Point3D getCenter() {
        return center;
    }

    public Point3D getU() {
        return u;
    }

    public Point3D getV() {
        return v;
    }
}
