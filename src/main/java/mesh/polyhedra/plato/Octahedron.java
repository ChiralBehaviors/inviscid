package mesh.polyhedra.plato;

import javax.vecmath.Vector3d;

import mesh.Face;

/**
 * An implementation of a regular octahedron mesh.
 * 
 * @author Brian Yao
 */
public class Octahedron extends PlatonicSolid {

	private static final double RADIUS_FACTOR = 1.0 / Math.sqrt(2.0);
	
	/**
	 * Construct a octahedron mesh centered at the origin with the specified
	 * edge length.
	 * 
	 * @param edgeLength The length of each edge of this mesh.
	 */
	public Octahedron(double edgeLength) {
		super(edgeLength);
		
		// Generate vertices
		double distOrigin = edgeLength / Math.sqrt(2.0);
		Vector3d v0 = new Vector3d(0.0, distOrigin, 0.0);
		Vector3d v5 = new Vector3d(0.0, -distOrigin, 0.0);
		Vector3d[] middleVertices = {new Vector3d(distOrigin, 0.0, 0.0), new Vector3d(0.0, 0.0, distOrigin),
				new Vector3d(-distOrigin, 0.0, 0.0), new Vector3d(0.0, 0.0, -distOrigin)};
		addVertexPositions(middleVertices);
		addVertexPosition(v0);
		addVertexPosition(v5);
		
		// Generate faces
		for (int i = 0 ; i < middleVertices.length ; i++) {
			int next = (i + 1) % middleVertices.length;
			
			Face topFace = new Face(3);
			topFace.setAllVertexIndices(4, next, i);
			
			Face bottomFace = new Face(3);
			bottomFace.setAllVertexIndices(5, i, next);
			
			addFaces(topFace, bottomFace);
		}
		setVertexNormalsToFaceNormals();	
	}
	
	/**
	 * Construct a octahedron mesh with the specified circumradius.
	 * 
	 * @param circumradius The circumradius this polyhedron will have.
	 * @param dummy        A dummy variable.
	 */
	public Octahedron(double circumradius, boolean dummy) {
		this(circumradius / RADIUS_FACTOR);
	}

}
