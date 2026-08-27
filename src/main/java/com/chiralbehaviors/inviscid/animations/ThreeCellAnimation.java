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
 * THREE CELLS: VE - hole - VE, the reciprocal condition.
 *
 * The middle cell M (phase b = a + 60) shares a DIFFERENT triangular face with
 * each of two outer cells P and Q (both phase a), along two adjacent cube
 * diagonals n1, n2 with n1.n2 = 1/3. Each outer cell is positioned ONLY by its
 * own shared face with M:
 *
 *     centre(P) = n1 * Z(cos a + cos b)      centre(Q) = n2 * Z(cos a + cos b)
 *
 * Two cells were solvable by construction. Three are not: the middle cell has
 * to satisfy two shared-face constraints at once, and nothing was left free to
 * make that happen. Measured, both mate to <= 3.4e-15 at every angle.
 *
 * AND THE OUTER CELLS TOUCH EACH OTHER WITHOUT BEING ASKED TO. Nothing in the
 * construction relates P to Q -- each is placed solely by M. Yet they share
 * vertices, and the count runs
 *
 *     a = 0    ->  4 shared vertices     (the full VE-VE "square face" contact)
 *     interior ->  2
 *     a = -60  ->  1
 *
 * which is exactly the 4 -> 2 -> 1 decay this project had recorded as a
 * separate property of the honeycomb packing. It is not separate. It is a
 * consequence of two triangular-face constraints on a common neighbour.
 */
public class ThreeCellAnimation extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            ThreeCellAnimation.main(argv);
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

    private static int faceAlong(Jitterbug j, Point3D dir) {
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

    private static double faceGap(List<Point3D> u, List<Point3D> v) {
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

    private static List<Point3D> verts(Jitterbug j, Point3D sh) {
        List<Point3D> out = new ArrayList<>();
        for (List<Point3D> f : faces(j, sh)) {
            for (Point3D q : f) {
                boolean dup = false;
                for (Point3D r : out) {
                    if (r.distance(q) < 1e-6) {
                        dup = true;
                        break;
                    }
                }
                if (!dup) {
                    out.add(q);
                }
            }
        }
        return out;
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group root = new Group();

        Octahedron oct = PhiCoordinates.Octahedrons[4];
        final double Z = oct.getEdgeLength() * Math.sqrt(2) / Math.sqrt(3);

        final Jitterbug M = new Jitterbug(oct, materials);
        final Jitterbug P = new Jitterbug(oct, materials);
        final Jitterbug Q = new Jitterbug(oct, materials);
        M.rotateTo(0);

        List<List<Point3D>> mf = faces(M, null);
        final Point3D n1 = cen(mf.get(0)).normalize();
        int i2 = 1;
        double bd = 9;
        for (int k = 1; k < 8; k++) {
            double d = Math.abs(cen(mf.get(k)).normalize()
                                              .dotProduct(n1)
                                - 1.0 / 3);
            if (d < bd) {
                bd = d;
                i2 = k;
            }
        }
        final Point3D n2 = cen(mf.get(i2)).normalize();

        final Translate sp = new Translate();
        final Translate sq = new Translate();
        P.getGroup()
         .getTransforms()
         .add(sp);
        Q.getGroup()
         .getTransforms()
         .add(sq);
        root.getChildren()
            .addAll(M.getGroup(), P.getGroup(), Q.getGroup());
        content.setContent(root);

        final long[] last = { 0 };
        final Timeline timeline = new Timeline();
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(12_000), new KeyValue(new WritableValue<Double>() {
                    @Override
                    public Double getValue() {
                        return -60d;
                    }

                    @Override
                    public void setValue(Double value) {
                        double a = value;
                        double b = a + PHASE_OFFSET;
                        double sep = Z * (Math.cos(Math.toRadians(a)) + Math.cos(Math.toRadians(b)));
                        M.rotateTo(b);
                        P.rotateTo(a);
                        Q.rotateTo(a);
                        Point3D dp = n1.multiply(sep), dq = n2.multiply(sep);
                        sp.setX(dp.getX());
                        sp.setY(dp.getY());
                        sp.setZ(dp.getZ());
                        sq.setX(dq.getX());
                        sq.setY(dq.getY());
                        sq.setZ(dq.getZ());

                        long now = System.currentTimeMillis();
                        if (now - last[0] < 500) {
                            return;
                        }
                        last[0] = now;
                        int mp = faceAlong(M, n1), mq = faceAlong(M, n2);
                        int pm = faceAlong(P, n1.multiply(-1)), qm = faceAlong(Q, n2.multiply(-1));
                        double g1 = faceGap(faces(M, null).get(mp), faces(P, dp).get(pm));
                        double g2 = faceGap(faces(M, null).get(mq), faces(Q, dq).get(qm));
                        List<Point3D> vp = verts(P, dp), vq = verts(Q, dq);
                        int shared = 0;
                        for (Point3D u : vp) {
                            for (Point3D v : vq) {
                                if (u.distance(v) < 1e-6) {
                                    shared++;
                                }
                            }
                        }
                        System.out.printf("a=%+7.2f b=%+7.2f | M-P gap=%9.3e  M-Q gap=%9.3e "
                                          + "| P^Q shared verts=%d | separation=%9.6f  |P-Q|=%9.6f%n",
                                          a, b, g1, g2, shared, sep, dp.distance(dq));
                    }
                }, 0d)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
