package mesh.polyhedra.plato;

import javax.vecmath.Vector3d;

import mesh.Edge;
import mesh.polyhedra.Polyhedron;

/**
 * An implementation of a regular icosahedron mesh.
 * 
 * @author Brian Yao
 */
public class Icosahedron extends PlatonicSolid {
	
	private static final double RADIUS_FACTOR = Math.sin(2.0 * Math.PI / 5.0);
	
	/**
	 * Construct an icosahedron mesh centered at the origin with the specified
	 * edge length.
	 * 
	 * @param edgeLength The length of each edge of this mesh.
	 */
	public Icosahedron(double edgeLength) {
		super(edgeLength);
		
		// An icosahedron is the dual polyhedron of a dodecahedron
		Dodecahedron dodec = new Dodecahedron(1.0);
		Polyhedron dualDodecahedron = dodec.dual();
		addVertexPositions(dualDodecahedron.getVertexPositions());
		addVertexNormals(dualDodecahedron.getVertexNormals());
		addFaces(dualDodecahedron.getFaces());
		
		// Scale vertex positions
		Edge sampleEdge = new Edge(faces.get(0).getVertexIndex(0), faces.get(0).getVertexIndex(1));
		sampleEdge.setMesh(this);
		double scaleFactor = edgeLength / sampleEdge.length();
		for (Vector3d vertexPos : vertexPositions) {
			vertexPos.scale(scaleFactor);
		}
	}
	
	/**
	 * Construct a icosahedron mesh with the specified circumradius.
	 * 
	 * @param circumradius The circumradius this polyhedron will have.
	 * @param dummy        A dummy variable.
	 */
	public Icosahedron(double circumradius, boolean dummy) {
		this(circumradius / RADIUS_FACTOR);
	}
	
}
