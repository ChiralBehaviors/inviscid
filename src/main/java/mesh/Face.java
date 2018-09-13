package mesh;

import java.util.Arrays;

import javafx.geometry.Point3D;

/**
 * A class for polygon faces. These faces store vertex data in the form of
 * indices; these indices are defined by a mesh, so an instance of Face is not
 * meaningful without a mesh for it to belong to. The design of this class is
 * suited for the OBJ file format, just as the Mesh class is.
 *
 * The face stores a pointer to the mesh it belongs to, but this is null until
 * this face is added to a mesh using one of the methods in Mesh.
 *
 * By convention, the indexes of the vertex positions and normals are assumed to
 * be specified in counterclockwise order if we are looking at the face from the
 * "outside". Thus, several of the methods will rely on this convention.
 *
 * Currently, the only "data" stored by the faces are the per-vertex position
 * and normal. Both are stored as integer indices, where the index refers to the
 * index of the actual value in the Mesh this face belongs to.
 *
 * @author Brian Yao
 */
public class Face {

    private Mesh  mesh;

    private int   numVertices;
    private int[] vertexIndices;
    private int[] normalIndices;

    /**
     * Create a face with enough space to store data for the specified number of
     * vertices.
     * 
     * @param numVertices
     *            The number of vertices this face has.
     */
    public Face(int numVertices) {
        if (numVertices < 3) {
            throw new IllegalArgumentException("A polygon face cannot have fewer than 3 vertices.");
        }
        this.numVertices = numVertices;
        vertexIndices = new int[numVertices];
        normalIndices = new int[numVertices];
    }

    /**
     * Copy constructor.
     * 
     * @param face
     *            The face to copy.
     */
    public Face(Face face) {
        this(face.numVertices);
        System.arraycopy(face.vertexIndices, 0, vertexIndices, 0,
                         face.numVertices);
        System.arraycopy(face.normalIndices, 0, normalIndices, 0,
                         face.numVertices);
    }

    /**
     * Compute the vector average of all vertices of this face.
     * 
     * @return The centroid of this face.
     */
    public Point3D vertexAverage() {
        Point3D avg = new Point3D(0, 0, 0);
        for (int v : vertexIndices) {
            avg = avg.add(mesh.vertexPositions.get(v));
        }
        return avg.multiply(1.0 / numVertices);
    }

    /**
     * Compute the centroid of this face. This method will not succeed unless
     * setMesh() has been called on this face; this can be done implicitly by
     * adding this face to a mesh (see addFace() in Mesh).
     * 
     * @return The centroid location.
     */
    public Point3D centroid() {
        Face[] triangles = divideIntoTriangles();
        double totalArea = 0.0;
        double[] triangleAreas = new double[triangles.length];
        Point3D[] triangleCentroids = new Point3D[triangles.length];
        for (int i = 0; i < triangles.length; i++) {
            Point3D v0 = mesh.vertexPositions.get(triangles[i].getVertexIndex(0));
            Point3D v1 = mesh.vertexPositions.get(triangles[i].getVertexIndex(1));
            Point3D v2 = mesh.vertexPositions.get(triangles[i].getVertexIndex(2));

            // Compute centroid of the triangle
            Point3D triangleCentroid = new Point3D(0, 0, 0);
            triangleCentroid = triangleCentroid.add(v0);
            triangleCentroid = triangleCentroid.add(v1);
            triangleCentroid = triangleCentroid.add(v2);
            triangleCentroid = triangleCentroid.multiply(1.0 / 3.0);
            triangleCentroids[i] = triangleCentroid;

            // Compute area of the triangle
            Point3D diff01 = v1.subtract(v0);
            Point3D diff02 = v2.subtract(v0);
            Point3D crossProd = diff01.crossProduct(diff02);
            double triangleArea = 0.5 * crossProd.magnitude();

            totalArea += triangleArea;
            triangleAreas[i] = triangleArea;
        }

        // Compute centroid as a weighted sum of triangle centroids
        Point3D centroid = new Point3D(0, 0, 0);
        for (int i = 0; i < triangles.length; i++) {
            centroid.add(triangleCentroids[i].multiply(triangleAreas[i]
                                                       / totalArea));
        }

        return centroid;
    }

    /**
     * Converts this polygon face into an array of triangular faces whose union
     * has the same geometry as this face. If this face is a triangle, the array
     * will only have one face. In general, the number of triangular faces in
     * the returned array is numVertices - 2.
     * 
     * @return The triangular faces whose union is this one.
     */
    public Face[] divideIntoTriangles() {
        Face[] triangles = new Face[numVertices - 2];
        for (int i = 0; i < numVertices - 2; i++) {
            int v0 = i;
            int v1 = i + 1;
            int v2 = numVertices - 1;
            Face triangle = new Face(3);
            triangle.setAllVertexIndices(vertexIndices[v0], vertexIndices[v1],
                                         vertexIndices[v2]);
            triangle.setAllNormalIndices(normalIndices[v0], normalIndices[v1],
                                         normalIndices[v2]);
            triangles[i] = triangle;
        }
        return triangles;
    }

    /**
     * Compute the edges which bound this face. The edges will be specified in
     * the order that the vertices are specified. If convention is followed,
     * such that the face's vertices are specified in counterclockwise order,
     * then the edges will also be in counterclockwise order.
     * 
     * Each returned edge will have setMesh() called on it, with the same mesh
     * this face points to (which could be null).
     * 
     * @return An array of edges bounding this face.
     */
    public Edge[] getEdges() {
        Edge[] edges = new Edge[numVertices];
        for (int i = 0; i < numVertices; i++) {
            int next = (i + 1) % numVertices;
            Edge ithEdge = new Edge(vertexIndices[i], vertexIndices[next]);
            ithEdge.setMesh(mesh);
            edges[i] = ithEdge;
        }
        return edges;
    }

    /**
     * Compute the unit normal vector perpendicular to the plane this face lies
     * in. This requires this face to be assigned to a mesh (see addFace() in
     * Mesh), and assumes that the vertices of this face are specified in
     * counterclockwise order.
     * 
     * @return The 3D vector containing the normal to the face.
     */
    public Point3D getFaceNormal() {
        Point3D v0 = mesh.vertexPositions.get(vertexIndices[0]);
        Point3D diff01 = mesh.vertexPositions.get(vertexIndices[1])
                                             .subtract(v0);
        Point3D diff0n = mesh.vertexPositions.get(vertexIndices[numVertices
                                                                - 1])
                                             .subtract(v0);

        return diff01.crossProduct(diff0n).normalize();
    }

    /**
     * @param vertex
     *            The index 0 <= vertex < numVertices of the vertex whose normal
     *            (index) to retrieve.
     * @return The index of the specified vertex's normal.
     */
    public int getNormalIndex(int vertex) {
        return normalIndices[vertex];
    }

    /**
     * @return An array of all normal indices for this face.
     */
    public int[] getNormalIndices() {
        return normalIndices;
    }

    /**
     * @param vertex
     *            The index 0 <= vertex < numVertices of the vertex whose
     *            position (index) to retrieve.
     * @return The index of the specified vertex's position.
     */
    public int getVertexIndex(int vertex) {
        return vertexIndices[vertex];
    }

    /**
     * @return An array of all position indices for this face.
     */
    public int[] getVertexIndices() {
        return vertexIndices;
    }

    /**
     * @return The number of vertices this face has.
     */
    public int numVertices() {
        return numVertices;
    }

    /**
     * Sets the vertex normal indices for this face, in the order given. The
     * first normal specified corresponds to the vertex at the first position
     * stored in this face.
     * 
     * @param normalIndices
     *            Set the normal indices to the contents of the provided array.
     */
    public void setAllNormalIndices(int... normalIndices) {
        System.arraycopy(normalIndices, 0, normalIndices, 0,
                         normalIndices.length);
    }

    /**
     * Set all vertex normals to the specified index.
     * 
     * @param normalIndex
     *            The index of the normal all vertices will use.
     */
    public void setAllVertexIndicesTo(int normalIndex) {
        for (int i = 0; i < numVertices; i++) {
            normalIndices[i] = normalIndex;
        }
    }

    /**
     * Sets the vertex position indices for this face, in the order given.
     * 
     * @param positionIndices
     *            Set the position indices to the contents of the provided
     *            array.
     */
    public void setAllVertexIndices(int... positionIndices) {
        System.arraycopy(positionIndices, 0, vertexIndices, 0,
                         vertexIndices.length);
    }

    /**
     * Set the mesh that this face points to. By default, this is called when
     * this face is added to a mesh using one of addFace() or similar methods in
     * Mesh.
     * 
     * This method needs to be called at some point for certain methods in this
     * class and others to work properly, since Face objects do not store any
     * geometry; all geometry is stored in the Mesh.
     * 
     * @param mesh
     *            The mesh this face will now point to.
     */
    public void setMesh(Mesh mesh) {
        this.mesh = mesh;
    }

    /**
     * Set the position and normal for the vertex of this face stored at index
     * "vertex" (the data of this Face is stored in arrays; the index
     * corresponds to the index of the vertex in these arrays).
     * 
     * @param vertex
     *            The vertex whose data we want to set.
     * @param positionIndex
     *            The index of the vertex's position in the mesh.
     * @param normalIndex
     *            The index of the vertex's normal in the mesh.
     */
    public void setVertexData(int vertex, int positionIndex, int normalIndex) {
        setVertexIndex(vertex, positionIndex);
        setVertexNormal(vertex, normalIndex);
    }

    /**
     * Set the normal for the vertex of this face stored at index "vertex" (the
     * data of this Face is stored in arrays; the index corresponds to the index
     * of the vertex in these arrays).
     * 
     * @param vertex
     *            The vertex whose normal we want to set.
     * @param normalIndex
     *            The index of the vertex's normal in the mesh.
     */
    public void setVertexNormal(int vertex, int normalIndex) {
        normalIndices[vertex] = normalIndex;
    }

    /**
     * Set the position for the vertex of this face stored at index "vertex"
     * (the data of this Face is stored in arrays; the "vertex" corresponds to
     * the index of the vertex in these arrays).
     * 
     * @param vertex
     *            The vertex whose position we want to set.
     * @param positionIndex
     *            The index of the vertex's position in the mesh.
     */
    public void setVertexIndex(int vertex, int positionIndex) {
        vertexIndices[vertex] = positionIndex;
    }

    /**
     * Represents this face as a string, as it would be in the OBJ file format.
     *
     * @return The string with this face's data in OBJ format.
     */
    public String toOBJString() {
        StringBuilder builder = new StringBuilder("f ");
        for (int i = 0; i < vertexIndices.length; i++) {
            builder.append(String.format("%d//%d ", vertexIndices[i] + 1,
                                         normalIndices[i] + 1));
        }
        builder.deleteCharAt(builder.length() - 1);
        return builder.toString();
    }

    @Override
    public String toString() {
        return "Face:" + Arrays.toString(vertexIndices);
    }

}
