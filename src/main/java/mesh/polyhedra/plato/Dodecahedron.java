package mesh.polyhedra.plato;

import javafx.geometry.Point3D;

import mesh.Face;

/**
 * An implementation of a regular dodecahedron mesh.
 *
 * @author Brian Yao
 */
public class Dodecahedron extends PlatonicSolid {

    private static final double RADIUS_FACTOR = Math.sqrt(3.0) / 4.0
                                                * (1.0 + Math.sqrt(5));

    /**
     * Construct a dodecahedron mesh centered at the origin with the specified
     * edge length.
     * 
     * @param edgeLength
     *            The length of each edge of this mesh.
     */
    public Dodecahedron(double edgeLength) {
        super(edgeLength);

        // Construct vertices
        double goldenRatio = (1.0 + Math.sqrt(5.0)) / 2.0;
        double goldenRatioInv = 1.0 / goldenRatio;
        double edgeScaleFactor = edgeLength / (Math.sqrt(5.0) - 1.0);
        Point3D[] cubePoints = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            Point3D vcube = new Point3D((i & 1) == 1 ? -1.0 : 1.0,
                                        ((i >> 1) & 1) == 1 ? -1.0 : 1.0,
                                        ((i >> 2) & 1) == 1 ? -1.0 : 1.0)
                                                                         .multiply(edgeScaleFactor);
            cubePoints[i] = vcube;
        }

        Point3D[] greenVertices = new Point3D[4];
        Point3D[] pinkVertices = new Point3D[4];
        Point3D[] blueVertices = new Point3D[4];
        for (int i = 0; i < 4; i++) {
            Point3D vgreen = new Point3D((i & 1) == 1 ? -goldenRatio
                                                      : goldenRatio,
                                         ((i >> 1) & 1) == 1 ? -goldenRatioInv
                                                             : goldenRatioInv,
                                         0).multiply(edgeScaleFactor);
            greenVertices[i] = vgreen;

            Point3D vpink = new Point3D(((i >> 1) & 1) == 1 ? -goldenRatioInv
                                                            : goldenRatioInv,
                                        0,
                                        (i
                                         & 1) == 1 ? -goldenRatio
                                                   : goldenRatio).multiply(edgeScaleFactor);
            pinkVertices[i] = vpink;

            Point3D vblue = new Point3D(0,
                                        (i & 1) == 1 ? -goldenRatio
                                                     : goldenRatio,
                                        ((i >> 1)
                                         & 1) == 1 ? -goldenRatioInv
                                                   : goldenRatioInv).multiply(edgeScaleFactor);
            blueVertices[i] = vblue;
        }

        // Cube points: 0-7, green: 8-11, pink: 12-15, blue: 16-19
        addVertexPositions(cubePoints);
        addVertexPositions(greenVertices);
        addVertexPositions(pinkVertices);
        addVertexPositions(blueVertices);

        Face[] faces = new Face[12];
        for (int i = 0; i < faces.length; i++) {
            faces[i] = new Face(5);
        }

        // Construct faces
        faces[0].setAllVertexIndices(0, 16, 2, 14, 12);
        faces[1].setAllVertexIndices(1, 13, 15, 3, 18);
        faces[2].setAllVertexIndices(4, 12, 14, 6, 17);
        faces[3].setAllVertexIndices(5, 19, 7, 15, 13);
        faces[4].setAllVertexIndices(0, 12, 4, 10, 8);
        faces[5].setAllVertexIndices(2, 9, 11, 6, 14);
        faces[6].setAllVertexIndices(1, 8, 10, 5, 13);
        faces[7].setAllVertexIndices(3, 15, 7, 11, 9);
        faces[8].setAllVertexIndices(0, 8, 1, 18, 16);
        faces[9].setAllVertexIndices(4, 17, 19, 5, 10);
        faces[10].setAllVertexIndices(2, 16, 18, 3, 9);
        faces[11].setAllVertexIndices(6, 11, 7, 19, 17);

        addFaces(faces);
        setVertexNormalsToFaceNormals();
    }

    /**
     * Construct a dodecahedron mesh with the specified circumradius.
     * 
     * @param circumradius
     *            The circumradius this polyhedron will have.
     * @param dummy
     *            A dummy variable.
     */
    public Dodecahedron(double circumradius, boolean dummy) {
        this(circumradius / RADIUS_FACTOR);
    }

}
