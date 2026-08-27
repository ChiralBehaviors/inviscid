package com.chiralbehaviors.inviscid.animations;

import static com.chiralbehaviors.inviscid.animations.Colors.materials;

import java.util.ArrayList;
import java.util.List;

import com.chiralbehaviors.inviscid.Jitterbug;
import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.KeyFrame;
import javafx.animation.KeyValue;
import javafx.animation.Timeline;
import javafx.beans.value.WritableValue;
import javafx.geometry.Point3D;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Transform;
import javafx.scene.transform.Translate;
import javafx.util.Duration;
import mesh.polyhedra.plato.Octahedron;

/**
 * THREE CELLS IN A LINE, sequential along one axis.
 *
 * The middle cell uses two OPPOSITE faces. A jitterbug's 8 triangles sit on 4
 * axes, two faces per axis (Gray), so the two chosen faces are exactly
 * antipodal -- measured dot = -1.000000 -- and the chain is genuinely straight.
 *
 *     A ---- M ---- C          along n, the shared-face axis
 *   p+60     p     p+60
 *
 * PHASES ALTERNATE; THEY DO NOT RAMP. The shared-face law admits b = a +- 60
 * (the match is mod 120 because the face is an equilateral triangle), so a
 * chain could in principle be built as a phase RAMP p-60, p, p+60. It cannot:
 * every cell's fold angle must lie in [-60, 60], and a ramp pins p = 0, a
 * single frozen configuration with no motion left. The alternating choice is
 * the one that moves, and it is what the honeycomb does along a body diagonal.
 *
 * Both spacings stay EQUAL at every angle -- Z(cos p + cos(p+60)) on each side
 * -- so the line breathes uniformly, 9.069136 at the ends of the sweep out to
 * 10.472136 at p = -30 where all three cells are congruent. Verified: both
 * shared faces mate to <= 3.4e-15 across the range.
 *
 * WHAT THIS SHOWS ABOUT THE MEDIUM. A straight chain is the simplest thing that
 * could carry a longitudinal disturbance, and this one has exactly one degree
 * of freedom: p. Fix p and the entire chain is determined -- both neighbours,
 * both spacings. There is no way to give the far cell a different phase from
 * the near one while keeping the faces shared. The chain moves as a unit.
 */
public class LineChainAnimation extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            LineChainAnimation.main(argv);
        }
    }

    private static final double PHASE_OFFSET = 60.0;

    public static void main(String[] args) {
        launch(args);
    }

    private static List<List<Point3D>> faces(Jitterbug j, Point3D extra) {
        List<List<Point3D>> out = new ArrayList<>();
        for (Node ch : j.getGroup()
                        .getChildren()) {
            Group fg = (Group) ch;
            MeshView mv = (MeshView) fg.getChildren()
                                       .get(0);
            TriangleMesh tm = (TriangleMesh) mv.getMesh();
            float[] p = new float[tm.getPoints()
                                    .size()];
            tm.getPoints()
              .toArray(p);
            List<Point3D> f = new ArrayList<>();
            for (int i = 0; i < p.length; i += 3) {
                Point3D q = new Point3D(p[i], p[i + 1], p[i + 2]);
                for (Transform t : fg.getTransforms()) {
                    q = t.transform(q);
                }
                if (extra != null) {
                    q = q.add(extra);
                }
                f.add(q);
            }
            out.add(f);
        }
        return out;
    }

    private static Point3D cen(List<Point3D> f) {
        double x = 0, y = 0, z = 0;
        for (Point3D q : f) {
            x += q.getX();
            y += q.getY();
            z += q.getZ();
        }
        return new Point3D(x / f.size(), y / f.size(), z / f.size());
    }

    private static int along(Jitterbug j, Point3D dir) {
        List<List<Point3D>> f = faces(j, null);
        int bi = 0;
        double best = -2;
        for (int k = 0; k < 8; k++) {
            double d = cen(f.get(k)).normalize()
                                    .dotProduct(dir);
            if (d > best) {
                best = d;
                bi = k;
            }
        }
        return bi;
    }

    private static double gap(List<Point3D> u, List<Point3D> v) {
        double w = 0;
        for (Point3D p : u) {
            double d = Double.MAX_VALUE;
            for (Point3D q : v) {
                d = Math.min(d, p.distance(q));
            }
            w = Math.max(w, d);
        }
        return w;
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group root = new Group();

        Octahedron oct = PhiCoordinates.Octahedrons[4];
        final double Z = oct.getEdgeLength() * Math.sqrt(2) / Math.sqrt(3);

        final Jitterbug M = new Jitterbug(oct, materials);
        final Jitterbug A = new Jitterbug(oct, materials);
        final Jitterbug C = new Jitterbug(oct, materials);
        M.rotateTo(0);
        final Point3D n = cen(faces(M, null).get(0)).normalize();

        final Translate sa = new Translate();
        final Translate sc = new Translate();
        A.getGroup()
         .getTransforms()
         .add(sa);
        C.getGroup()
         .getTransforms()
         .add(sc);
        root.getChildren()
            .addAll(A.getGroup(), M.getGroup(), C.getGroup());
        content.setContent(root);

        final long[] last = { 0 };
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
                        double p = value;
                        double q = p + PHASE_OFFSET;
                        double sep = Z * (Math.cos(Math.toRadians(p)) + Math.cos(Math.toRadians(q)));
                        M.rotateTo(p);
                        A.rotateTo(q);
                        C.rotateTo(q);
                        Point3D da = n.multiply(-sep), dc = n.multiply(sep);
                        sa.setX(da.getX());
                        sa.setY(da.getY());
                        sa.setZ(da.getZ());
                        sc.setX(dc.getX());
                        sc.setY(dc.getY());
                        sc.setZ(dc.getZ());

                        long now = System.currentTimeMillis();
                        if (now - last[0] < 500) {
                            return;
                        }
                        last[0] = now;
                        int fp = along(M, n), fm = along(M, n.multiply(-1));
                        int am = along(A, n), cm = along(C, n.multiply(-1));
                        System.out.printf("p=%+7.2f  neighbours=%+7.2f | spacing=%9.6f (both) "
                                          + "| A-M gap=%9.3e  M-C gap=%9.3e  | chain length=%9.6f%n",
                                          p, q, sep, gap(faces(M, null).get(fm), faces(A, da).get(am)),
                                          gap(faces(M, null).get(fp), faces(C, dc).get(cm)), 2 * sep);
                    }
                };
        driver.setValue(-60d);

        final Timeline timeline = new Timeline();
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(12_000), new KeyValue(driver, 0d)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
