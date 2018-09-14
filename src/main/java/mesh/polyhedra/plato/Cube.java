package mesh.polyhedra.plato;

import javax.vecmath.Vector3d;

import math.VectorMath;
import mesh.Face;

/**
 * An implementation of a regular cube mesh.
 *
 * @author Brian Yao
 */
public class Cube extends PlatonicSolid {

    public static Vector3d[] standardPositions(double edgeLength) {
        // Generate vertex positions
        double halfEdgeLength = edgeLength / 2.0;
        Vector3d[] vs = new Vector3d[8];
        for (int i = 0; i < 8; i++) {
            Vector3d current = new Vector3d();
            current.x = (i & 1) == 1 ? halfEdgeLength : -halfEdgeLength;
            current.y = ((i >> 1) & 1) == 1 ? halfEdgeLength : -halfEdgeLength;
            current.z = ((i >> 2) & 1) == 1 ? halfEdgeLength : -halfEdgeLength;
            vs[i] = current;
        }
        return vs;
    }

    private static double edgeLength(Vector3d[] vs) {
        Vector3d diff = VectorMath.diff(vs[1], vs[0]);
        return diff.length();
    }

	/**
	 * Construct a cube mesh centered at the origin with the specified
	 * edge length.
	 * 
	 * @param edgeLength The length of each edge of this mesh.
	 */
	public Cube(double edgeLength) {
		this(edgeLength, standardPositions(edgeLength));
	}

    public Cube(double edgeLength, Vector3d[] vs) {
        super(edgeLength);

        addVertexPositions(vs);
        Face f0 = new Face(4);
        f0.setAllVertexIndices(0, 2, 3, 1);
        
        Face f1 = new Face(4);
        f1.setAllVertexIndices(0, 4, 6, 2);
        
        Face f2 = new Face(4);
        f2.setAllVertexIndices(0, 1, 5, 4);
        
        Face f3 = new Face(4);
        f3.setAllVertexIndices(7, 6, 4, 5);
        
        Face f4 = new Face(4);
        f4.setAllVertexIndices(7, 5, 1, 3);
        
        Face f5 = new Face(4);
        f5.setAllVertexIndices(7, 3, 2, 6);

        addFaces(f0, f1, f2, f3, f4, f5);
        setVertexNormalsToFaceNormals();
    }

    public Cube(Vector3d[] vs) {
        super(edgeLength(vs));
        addVertexPositions(vs);

        Face f0 = new Face(4);
        Face f1 = new Face(4);
        Face f2 = new Face(4);
        Face f3 = new Face(4);
        Face f4 = new Face(4);
        Face f5 = new Face(4);

        f0.setAllVertexIndices(0, 3, 2, 1);
        f1.setAllVertexIndices(1, 2, 6, 5);
        f2.setAllVertexIndices(5, 6, 7, 4);
        f3.setAllVertexIndices(4, 7, 3, 0);
        f4.setAllVertexIndices(0, 1, 5, 4);
        f5.setAllVertexIndices(2, 3, 7, 6);

        addFaces(f0, f1, f2, f3, f4, f5);
        setVertexNormalsToFaceNormals();
    }
}
