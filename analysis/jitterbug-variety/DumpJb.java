import com.chiralbehaviors.inviscid.PhiCoordinates;
import javax.vecmath.Vector3d;
import mesh.Face;
import mesh.polyhedra.plato.Octahedron;

public class DumpJb {
    public static void main(String[] args) {
        Octahedron oct = PhiCoordinates.Octahedrons[4];
        System.out.println("edgeLength = " + oct.getEdgeLength());
        int i = 0;
        for (Face f : oct.getFaces()) {
            Vector3d c = f.centroid();
            StringBuilder sb = new StringBuilder();
            for (Vector3d v : f.getVertices()) {
                sb.append(String.format("(%.6f,%.6f,%.6f) ", v.x, v.y, v.z));
            }
            double prod = Math.signum(c.x) * Math.signum(c.y) * Math.signum(c.z);
            System.out.printf(
                "face %d  inverse=%-5s  sigmaJava=%+d  centroid=(%.6f,%.6f,%.6f) |c|=%.6f  octantProduct=%+.0f  verts=%s%n",
                i, PhiCoordinates.JITTERBUG_INVERSES[i],
                PhiCoordinates.JITTERBUG_INVERSES[i] ? -1 : 1,
                c.x, c.y, c.z, c.length(), prod, sb.toString());
            i++;
        }
        // vertex positions
        System.out.println("\nvertices of the octahedron:");
        int k = 0;
        for (Vector3d v : oct.getVertexPositions()) {
            System.out.printf("  v%d = (%.6f,%.6f,%.6f)  |v|=%.6f%n", k++, v.x, v.y, v.z, v.length());
        }
    }
}
