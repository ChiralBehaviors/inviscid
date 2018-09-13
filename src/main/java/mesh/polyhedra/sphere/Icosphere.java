package mesh.polyhedra.sphere;

import javafx.geometry.Point3D;

import mesh.polyhedra.Polyhedron;
import mesh.polyhedra.plato.Icosahedron;

/**
 * An approximate sphere created by using a seed icosahedron and subdividing its
 * faces repeatedly, then projecting all vertices to the sphere of the desired
 * radius. Each subdivision divides each triangular face into four.
 *
 * The dual of an icosphere is a Goldberg polyhedron.
 *
 * @author Brian Yao
 */
public class Icosphere extends Polyhedron {

    /**
     * Construct an icosphere with the desired circumradius.
     * 
     * @param circumradius
     *            The circumradius of the resulting icosphere.
     * @param numSubdivisions
     *            The number of times to subdivide faces.
     */
    public Icosphere(double circumradius, int numSubdivisions) {
        Polyhedron icosphere = new Icosahedron(circumradius, true);
        for (int i = 0; i < numSubdivisions; i++) {
            icosphere = icosphere.subdivide();
        }
        addFaces(icosphere.getFaces());

        // Scale all vertices onto the sphere.
        for (Point3D vertexPos : icosphere.getVertexPositions()) {
            addVertexPosition(vertexPos.multiply(circumradius
                                                 / vertexPos.magnitude()));
        }

        setVertexNormalsToFaceNormals();
    }

}
