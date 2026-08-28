package com.chiralbehaviors.inviscid.animations;

import static com.chiralbehaviors.inviscid.animations.Colors.materials;

import com.chiralbehaviors.inviscid.Jitterbug;
import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import java.util.ArrayList;
import java.util.List;

import javafx.animation.KeyFrame;
import javafx.animation.KeyValue;
import javafx.animation.Timeline;
import javafx.beans.value.WritableValue;
import javafx.geometry.Point3D;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.paint.Color;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.Cylinder;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.Sphere;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Rotate;
import javafx.scene.transform.Transform;
import javafx.scene.transform.Translate;
import javafx.util.Duration;
import javax.vecmath.Vector3d;
import mesh.Face;

/**
 * TWO CELLS SHARING ONE TRIANGULAR FACE.
 *
 * Gray, "The Jitterbug Motion" (2002): "Two Jitterbugs can not share the same
 * triangular face and have their positions (location of center of volume) fixed
 * as they go through the Jitterbug motion. If two Jitterbugs are to share the
 * same triangle face then as the joined Jitterbugs jitterbug, the positions of
 * the Jitterbugs must move."
 *
 * So the second cell's centre is DRIVEN, not pinned. Solving the shared-face
 * condition gives, verified to 1e-15 over the whole range:
 *
 *   b = a + 60                    the second cell runs exactly 60 degrees ahead
 *   separation = Z(cos a + cos b) the centres, along the shared face normal
 *   Z = EL * sqrt(2/3)            Gray's radial law, V = EL cos(gamma) sqrt(2/3)
 *
 * The separation BREATHES: 9.069136 at a = -60, a maximum 10.472136 at a = -30
 * where the two cells are congruent, and back to 9.069136 at a = 0. The ratio
 * 10.472136 / 9.069136 = 1.154701 = 2/sqrt(3) is exactly the honeycomb lattice
 * parameter's range -- which therefore is not a property of the array at all.
 * It falls out of a single shared face.
 *
 * <p>
 * <b>What is instrumented, and why it is the point.</b> The two white spheres
 * are the cells' CENTRES OF VOLUME and the rod between them is their separation
 * — so Gray's "the positions of the Jitterbugs must move" is a thing you watch
 * happen rather than infer. The console prints that separation live, against its
 * own minimum, and alongside it the worst gap between cell A's shared face and
 * cell B's. That second number is what makes the demonstration a demonstration:
 * <b>the faces stay joined to ~1e-15 WHILE the centres travel by 15%.</b> Either
 * number alone proves nothing. Together they are exactly Gray's constraint —
 * share a face and the centres are not yours to fix.
 */
public class SharedFaceAnimation extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            SharedFaceAnimation.main(argv);
        }
    }

    private static final double PHASE_OFFSET = 60.0;

    public static void main(String[] args) {
        launch(args);
    }

    /** The transformed corners of one of a cell's 8 triangles, in world. Read
     *  off the rendered mesh rather than recomputed, so the number reported is
     *  the geometry actually on screen. */
    private static List<Point3D> corners(Jitterbug j, int face, Point3D shift) {
        Group fg = (Group) j.getGroup()
                            .getChildren()
                            .get(face);
        MeshView mv = (MeshView) fg.getChildren()
                                   .get(0);
        TriangleMesh tm = (TriangleMesh) mv.getMesh();
        float[] pts = new float[tm.getPoints()
                                  .size()];
        tm.getPoints()
          .toArray(pts);
        List<Point3D> out = new ArrayList<>();
        for (int i = 0; i < pts.length; i += 3) {
            Point3D q = new Point3D(pts[i], pts[i + 1], pts[i + 2]);
            for (Transform t : fg.getTransforms()) {
                q = t.transform(q);
            }
            out.add(shift == null ? q : q.add(shift));
        }
        return out;
    }

    /** Which of the cell's 8 faces points most nearly along {@code dir}. */
    private static int faceAlong(Jitterbug j, Point3D dir) {
        int best = 0;
        double bv = -2;
        for (int f = 0; f < 8; f++) {
            List<Point3D> c = corners(j, f, null);
            Point3D cen = c.get(0)
                           .add(c.get(1))
                           .add(c.get(2))
                           .multiply(1 / 3.0);
            double d = cen.normalize()
                          .dotProduct(dir);
            if (d > bv) {
                bv = d;
                best = f;
            }
        }
        return best;
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group group = new Group();

        var oct = PhiCoordinates.Octahedrons[4];
        final double Z = oct.getEdgeLength() * Math.sqrt(2) / Math.sqrt(3);

        // The shared face's axis, measured from the octahedron itself rather
        // than assumed: an earlier version hard-coded (1,1,1) and the two cells
        // missed each other by the full separation.
        Face f0 = oct.getFaces()
                     .get(0);
        Vector3d c = f0.centroid();
        final Point3D n = new Point3D(c.x, c.y, c.z).normalize();

        Jitterbug a = new Jitterbug(oct, materials);
        Jitterbug b = new Jitterbug(oct, materials);
        final Translate bShift = new Translate();
        final long[] last = { 0 };

        // The centres of volume, and the span between them: what Gray says
        // cannot be held still.
        final PhongMaterial chalk = new PhongMaterial(Color.WHITESMOKE);
        final Sphere centreA = new Sphere(0.35);
        final Sphere centreB = new Sphere(0.35);
        centreA.setMaterial(chalk);
        centreB.setMaterial(chalk);
        final Cylinder span = new Cylinder(0.10, 1);
        span.setMaterial(new PhongMaterial(Color.GOLD));
        // the rod's direction never changes -- only its length -- so the tilt
        // off the default Y axis is computed once
        final Rotate spanTilt = new Rotate(-Math.toDegrees(Math.acos(n.dotProduct(new Point3D(0, 1, 0)))),
                                           n.crossProduct(new Point3D(0, 1, 0)));
        b.getGroup()
         .getTransforms()
         .add(bShift);

        group.getChildren()
             .addAll(a.getGroup(), b.getGroup(), centreA, centreB, span);
        content.setContent(group);

        // INITIAL STATE. Jitterbug's constructor bakes every face rotated and
        // translated to the ORIGIN, so with zero transforms all 8 triangles of a
        // cell sit piled at its centre and every cell sits at zero separation.
        // rotateTo and the placement translations are what put them where they
        // belong, and until this runs the scene is meaningless -- cells heaped
        // on top of each other. The driver is therefore captured and invoked
        // ONCE here, so the view is correct before the timeline is ever started.
        final WritableValue<Double> driver = new WritableValue<Double>() {
                    @Override
                    public Double getValue() {
                        return -60d;
                    }

                    @Override
                    public void setValue(Double value) {
                        double ga = value;
                        double gb = ga + PHASE_OFFSET;
                        a.rotateTo(ga);
                        b.rotateTo(gb);
                        double sep = Z * (Math.cos(Math.toRadians(ga)) + Math.cos(Math.toRadians(gb)));
                        Point3D d = n.multiply(sep);
                        bShift.setX(d.getX());
                        bShift.setY(d.getY());
                        bShift.setZ(d.getZ());

                        // the centres, and the span between them
                        centreB.getTransforms()
                               .setAll(new Translate(d.getX(), d.getY(), d.getZ()));
                        span.setHeight(sep);
                        span.getTransforms()
                            .setAll(new Translate(d.getX() / 2, d.getY() / 2,
                                                  d.getZ() / 2),
                                    spanTilt);

                        if (System.currentTimeMillis() - last[0] < 400) {
                            return;
                        }
                        last[0] = System.currentTimeMillis();
                        // the shared face, read off BOTH cells' rendered geometry
                        List<Point3D> fa = corners(a, faceAlong(a, n), null);
                        List<Point3D> fb = corners(b, faceAlong(b, n.multiply(-1)),
                                                   d);
                        double gap = 0;
                        for (Point3D p : fa) {
                            double best = Double.MAX_VALUE;
                            for (Point3D q : fb) {
                                best = Math.min(best, p.distance(q));
                            }
                            gap = Math.max(gap, best);
                        }
                        System.out.printf("a=%+7.2f  b=%+7.2f | centres at 0 and "
                                          + "%8.5f  (%.4f x closest approach) | "
                                          + "shared face still joined to %8.2e%n",
                                          ga, gb, sep, sep / (1.5 * Z), gap);
                    }
                };
        driver.setValue(-60d);

        final Timeline timeline = new Timeline();
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(10_000), new KeyValue(driver, 0d)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
