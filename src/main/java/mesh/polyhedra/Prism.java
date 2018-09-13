package mesh.polyhedra;

import javafx.geometry.Point3D;

import mesh.Face;

/**
 * A class for generating right prism meshes.
 *
 * @author Brian Yao
 */
public class Prism extends Polyhedron {

    private static double RADIUS = 0.5;
    private static double HEIGHT = 1.0;

    /**
     * Constructs a prism whose base has the given circumradius and has the
     * given height. The base is a regular n-gon, where n is given as a
     * parameter. With high enough n, the geometry is an approximate cylinder.
     * 
     * @param numSides
     *            The number of sides the base has.
     * @param radius
     *            The circumradius of the base.
     * @param height
     *            The height of the prism.
     */
    public Prism(int numSides, double radius, double height) {
        Point3D toTop = new Point3D(0.0, 0.0, height);
        double centralAngle = 2.0 * Math.PI / numSides;

        Face bottomFace = new Face(numSides);
        Face topFace = new Face(numSides);
        int vertexIndex = 0;
        for (int i = 0; i < numSides; i++) {
            double totalAngle = centralAngle * i;
            double x = radius * Math.cos(totalAngle);
            double y = radius * Math.sin(totalAngle);
            Point3D bottomVertex = new Point3D(x, y, -height / 2.0);

            Point3D topVertex = bottomVertex.add(toTop);

            addVertexPositions(bottomVertex, topVertex);

            bottomFace.setVertexIndex(numSides - i - 1, vertexIndex);
            topFace.setVertexIndex(i, vertexIndex + 1);

            Face quad = new Face(4);
            int nextBottom = (vertexIndex + 2) % (2 * numSides);
            int nextTop = nextBottom + 1;
            quad.setAllVertexIndices(vertexIndex, nextBottom, nextTop,
                                     vertexIndex + 1);

            addFace(quad);

            vertexIndex += 2;
        }

        addFaces(bottomFace, topFace);

        setVertexNormalsToFaceNormals();
    }

    /**
     * Constructs a prism with a default radius (0.5) and default height (1.0).
     * The number of sides of the base is the only parameter.
     * 
     * @param numSides
     *            The number of sides the base has.
     */
    public Prism(int numSides) {
        this(numSides, RADIUS, HEIGHT);
    }

}
