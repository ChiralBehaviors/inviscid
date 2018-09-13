package mesh.polyhedra.plato;

import javafx.geometry.Point3D;

import mesh.Face;

/**
 * An implementation of a regular tetrahedron mesh.
 *
 * @author Brian Yao
 */
public class Tetrahedron extends PlatonicSolid {

    private static final double RADIUS_FACTOR = Math.sqrt(6.0) / 4.0;

    /**
     * Construct a tetrahedron mesh centered at the origin with the specified
     * edge length.
     * 
     * @param edgeLength
     *            The length of each edge of this mesh.
     */
    public Tetrahedron(double edgeLength) {
        super(edgeLength);

        // Vertex positions
        double edgeScaleFactor = edgeLength / (2.0 * Math.sqrt(2.0));
        Point3D v0 = new Point3D(1.0, 1.0, 1.0).multiply(edgeScaleFactor);
        Point3D v1 = new Point3D(1.0, -1.0, -1.0).multiply(edgeScaleFactor);
        Point3D v2 = new Point3D(-1.0, 1.0, -1.0).multiply(edgeScaleFactor);
        Point3D v3 = new Point3D(-1.0, -1.0, 1.0).multiply(edgeScaleFactor);
        addVertexPositions(v0, v1, v2, v3);

        // Create faces
        Face f0 = new Face(3);
        f0.setAllVertexIndices(0, 3, 1);

        Face f1 = new Face(3);
        f1.setAllVertexIndices(0, 1, 2);

        Face f2 = new Face(3);
        f2.setAllVertexIndices(2, 1, 3);

        Face f3 = new Face(3);
        f3.setAllVertexIndices(0, 2, 3);

        addFaces(f0, f1, f2, f3);
        setVertexNormalsToFaceNormals();
    }

    /**
     * Construct a tetrahedron mesh with the specified circumradius.
     * 
     * @param circumradius
     *            The circumradius this polyhedron will have.
     * @param dummy
     *            A dummy variable.
     */
    public Tetrahedron(double circumradius, boolean dummy) {
        this(circumradius / RADIUS_FACTOR);
    }

}
