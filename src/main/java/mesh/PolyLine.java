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

import java.util.Arrays;
import java.util.List;

import javafx.geometry.Point3D;
import javafx.scene.shape.DrawMode;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;

/**
 * @author halhildebrand
 *
 */
public class PolyLine extends MeshView {
    public PolyLine(List<Point3D> points, double width) {
        //        setDepthTest(DepthTest.ENABLE);
        TriangleMesh mesh = new TriangleMesh();
        setMesh(mesh);
        //add each point. For each point add another point shifted on Z axis by width
        //This extra point allows us to build triangles later
        for (Point3D point : points) {
            mesh.getPoints()
                .addAll((float) point.getX(), (float) point.getY(),
                        (float) point.getZ());
            mesh.getPoints()
                .addAll((float) point.getX(), (float) point.getY(),
                        (float) point.getZ() + (float) width);
        }
        //add dummy Texture Coordinate
        mesh.getTexCoords()
            .addAll(0, 0);
        //Now generate trianglestrips for each line segment
        for (int i = 2; i < points.size() * 2; i += 2) { //add each segment
            //Vertices wound counter-clockwise which is the default front face of any Triange
            //These triangles live on the frontside of the line facing the camera
            mesh.getFaces()
                .addAll(i, 0, i - 2, 0, i + 1, 0); //add primary face
            mesh.getFaces()
                .addAll(i + 1, 0, i - 2, 0, i - 1, 0); //add secondary Width face
            //Add the same faces but wind them clockwise so that the color looks correct when camera is rotated
            //These triangles live on the backside of the line facing away from initial the camera
            mesh.getFaces()
                .addAll(i + 1, 0, i - 2, 0, i, 0); //add primary face
            mesh.getFaces()
                .addAll(i - 1, 0, i - 2, 0, i + 1, 0); //add secondary Width face
        }
        setDrawMode(DrawMode.FILL); //Fill so that the line shows width
//        setCullFace(CullFace.BACK);
    }

    public PolyLine(Point3D a, Point3D b, float width) {
        this(Arrays.asList(a, b), width);
    }
}
