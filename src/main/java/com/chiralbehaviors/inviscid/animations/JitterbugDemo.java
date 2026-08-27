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
import javax.vecmath.Vector3d;
import mesh.polyhedra.plato.Octahedron;

/**
 * EVERY LOAD-BEARING CLAIM, RUNNING SIDE BY SIDE WITH ITS FAILURE MODES.
 *
 * Four groups, one timeline, and a measured readout each half second. The point
 * is that the correct construction is not asserted -- it is the only one of the
 * four whose shared-face vertex gap stays at machine zero.
 *
 *   [1] SINGLE CELL. Sweeps gamma -60 -> +60 THROUGH the VE, which is Gray's
 *       actual range ("gamma = -60 being the first Octahedron position,
 *       gamma = 0 being the VE position, and gamma = +60 being the second").
 *       Readout: distinct vertex count (12, or 6 at the octahedra), the worst
 *       coincidence error among vertices that should be shared, the centre of
 *       volume (fixed), and the circumradius (breathing by sqrt(2)).
 *
 *   [2] TWO CELLS, CENTRES PINNED -- WRONG. Separation frozen at its a=0 value.
 *       This is what the nine-cell array animation did. Gray: "Two Jitterbugs
 *       can not share the same triangular face and have their positions
 *       (location of center of volume) fixed as they go through the Jitterbug
 *       motion." The gap opens and closes as the cells breathe against a fixed
 *       spacing.
 *
 *   [3] TWO CELLS, SAME FOLD ANGLE -- WRONG. Centres driven correctly but
 *       b = a instead of a + 60. Leaves a CONSTANT gap of one triangle
 *       circumradius EL/sqrt(3) = 4.275 at every angle: two concentric
 *       equilateral triangles whose corresponding vertices are one circumradius
 *       apart are rotated 60 degrees. A constant error that does not vary with
 *       the parameter is a fixed rotational offset, and 60 is the answer.
 *
 *   [4] TWO CELLS, CORRECT. b = a + 60, separation = Z(cos a + cos b),
 *       Z = EL sqrt(2/3). Gap stays at ~1e-15 for the whole sweep.
 */
public class JitterbugDemo extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            JitterbugDemo.main(argv);
        }
    }

    private static final double PHASE_OFFSET = 60.0;
    private static final double LANE         = 26.0;

    public static void main(String[] args) {
        launch(args);
    }

    /** World positions of one cell's 8 faces x 3 corners, transforms applied. */
    private static List<List<Point3D>> faces(Jitterbug j, Translate extra) {
        List<List<Point3D>> out = new ArrayList<>();
        for (Node child : j.getGroup()
                           .getChildren()) {
            Group fg = (Group) child;
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
                    q = q.add(extra.getX(), extra.getY(), extra.getZ());
                }
                f.add(q);
            }
            out.add(f);
        }
        return out;
    }

    private static Point3D centroid(List<Point3D> f) {
        double x = 0, y = 0, z = 0;
        for (Point3D q : f) {
            x += q.getX();
            y += q.getY();
            z += q.getZ();
        }
        return new Point3D(x / f.size(), y / f.size(), z / f.size());
    }

    /** Worst distance from a corner of A's shared face to the nearest corner of B's. */
    private static double sharedFaceGap(Jitterbug a, Jitterbug b, Translate bShift, Point3D n) {
        List<List<Point3D>> fa = faces(a, null);
        List<List<Point3D>> fb = faces(b, bShift);
        // Select B's shared face by ITS OWN face direction, not by the A->B
        // centre difference. When the faces coincide -- which is the whole
        // point -- that difference is ~0, and normalize() of a zero vector
        // returns garbage, so the first version of this measurement picked a
        // neighbouring face and reported a gap of EL at exactly the angles
        // where the construction is exact. The harness was broken where the
        // thing it measures works best.
        List<List<Point3D>> fbLocal = faces(b, null);
        int bi = 0;
        double best = Double.MAX_VALUE;
        for (int k = 0; k < 8; k++) {
            double d = centroid(fbLocal.get(k))
                                       .normalize()
                                       .dotProduct(n);
            if (d < best) {
                best = d;
                bi = k;
            }
        }
        double worst = 0;
        for (Point3D p : fa.get(0)) {
            double d = Double.MAX_VALUE;
            for (Point3D q : fb.get(bi)) {
                d = Math.min(d, p.distance(q));
            }
            worst = Math.max(worst, d);
        }
        return worst;
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group root = new Group();

        Octahedron oct = PhiCoordinates.Octahedrons[4];
        final double EL = oct.getEdgeLength();
        final double Z = EL * Math.sqrt(2) / Math.sqrt(3);
        Vector3d c0 = oct.getFaces()
                         .get(0)
                         .centroid();
        final Point3D n = new Point3D(c0.x, c0.y, c0.z).normalize();
        final double sep0 = Z * (Math.cos(0) + Math.cos(Math.toRadians(PHASE_OFFSET)));

        final Jitterbug solo = new Jitterbug(oct, materials);
        final Jitterbug[] A = new Jitterbug[3];
        final Jitterbug[] B = new Jitterbug[3];
        final Translate[] shift = new Translate[3];
        for (int k = 0; k < 3; k++) {
            A[k] = new Jitterbug(oct, materials);
            B[k] = new Jitterbug(oct, materials);
            shift[k] = new Translate();
            B[k].getGroup()
                .getTransforms()
                .add(shift[k]);
        }

        Group[] lanes = new Group[4];
        for (int k = 0; k < 4; k++) {
            lanes[k] = new Group();
            lanes[k].getTransforms()
                    .add(new Translate((k - 1.5) * LANE, 0, 0));
            root.getChildren()
                .add(lanes[k]);
        }
        lanes[0].getChildren()
                .add(solo.getGroup());
        for (int k = 0; k < 3; k++) {
            lanes[k + 1].getChildren()
                        .addAll(A[k].getGroup(), B[k].getGroup());
        }
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
                        double fa = Z * Math.cos(Math.toRadians(a));
                        double fb = Z * Math.cos(Math.toRadians(b));

                        solo.rotateTo(a);

                        // [2] pinned centres, correct phases
                        A[0].rotateTo(a);
                        B[0].rotateTo(b);
                        set(shift[0], n.multiply(sep0));

                        // [3] moving centres, WRONG phase (b = a)
                        A[1].rotateTo(a);
                        B[1].rotateTo(a);
                        set(shift[1], n.multiply(2 * fa));

                        // [4] correct
                        A[2].rotateTo(a);
                        B[2].rotateTo(b);
                        set(shift[2], n.multiply(fa + fb));

                        long now = System.currentTimeMillis();
                        if (now - last[0] < 500) {
                            return;
                        }
                        last[0] = now;
                        report(a, b, fa, fb, solo, A, B, shift, n);
                    }
                }, 0d)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }

    private static void set(Translate t, Point3D d) {
        t.setX(d.getX());
        t.setY(d.getY());
        t.setZ(d.getZ());
    }

    private static void report(double a, double b, double fa, double fb, Jitterbug solo, Jitterbug[] A, Jitterbug[] B,
                               Translate[] shift, Point3D n) {
        List<List<Point3D>> sf = faces(solo, null);
        List<Point3D> all = new ArrayList<>();
        for (List<Point3D> f : sf) {
            all.addAll(f);
        }
        List<Point3D> reps = new ArrayList<>();
        double worstShare = 0;
        for (Point3D q : all) {
            double best = Double.MAX_VALUE;
            for (Point3D r : reps) {
                best = Math.min(best, r.distance(q));
            }
            if (best < 1e-3) {
                worstShare = Math.max(worstShare, best);
            } else {
                reps.add(q);
            }
        }
        double cx = 0, cy = 0, cz = 0, rad = 0;
        for (Point3D q : all) {
            cx += q.getX();
            cy += q.getY();
            cz += q.getZ();
            rad = Math.max(rad, q.magnitude());
        }
        int m = all.size();
        double cen = new Point3D(cx / m, cy / m, cz / m).magnitude();

        System.out.printf("a=%+7.2f b=%+7.2f | SOLO verts=%2d share=%.1e centre=%.1e R=%8.6f "
                          + "| PINNED gap=%9.3e | SAME-GAMMA gap=%9.3e | CORRECT gap=%9.3e%n",
                          a, b, reps.size(), worstShare, cen, rad, sharedFaceGap(A[0], B[0], shift[0], n),
                          sharedFaceGap(A[1], B[1], shift[1], n), sharedFaceGap(A[2], B[2], shift[2], n));
    }
}
