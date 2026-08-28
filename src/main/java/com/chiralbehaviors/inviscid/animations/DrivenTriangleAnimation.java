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
 * DRIVE ONE TRIANGLE AT THE END OF THE CHAIN. NOTHING ELSE IS DRIVEN.
 *
 * Three cells in a line, A - M - C. The ONLY input is the own-axis rotation
 * angle theta of a single triangle -- face 1, the BLUE one, on the end cell A.
 * Every other quantity in the scene is solved from it:
 *
 *     theta          the driven blue triangle (JITTERBUG_INVERSES[1] = false,
 *                    so its own-axis angle IS its cell's fold angle)
 *     q = theta      cell A's fold angle, and cell C's
 *     p = theta - 60 cell M's fold angle
 *     spacing        Z(cos p + cos q), both sides
 *
 * THE CONSTRAINTS ARE IMMEDIATE. There is no propagation and no lag anywhere in
 * this scene: the far cell C's phase and position are computed in the SAME
 * frame as the driven triangle's angle, because they are algebraic consequences
 * of it, not the result of anything travelling. Turn the blue triangle and the
 * far end of the chain has already moved. That is what "shared vertex, rigid ->
 * constraint propagates instantly" means, and the chain's single degree of
 * freedom is why: fix theta and the whole configuration is determined.
 *
 * AND THAT IS EXACTLY WHY THIS IS NOT YET A WAVE. A wave is a disturbance that
 * takes TIME to cross the medium, and time cannot enter through the
 * constraints -- they are instantaneous by construction. It can only enter
 * through INERTIA. Nothing in this file has mass; the readout's dC/dtheta is a
 * constraint response, a ratio of displacements, not a speed. What a real
 * disturbance does here is a dynamical question this animation does not answer
 * and does not pretend to: give the triangles mass, integrate, and the far end
 * lags the near one by however long the inertia takes to be pushed. The
 * kinematics fixes WHAT moves together; only the dynamics fixes WHEN.
 */
public class DrivenTriangleAnimation extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            DrivenTriangleAnimation.main(argv);
        }
    }

    /** materials[1] is blueMaterial, and faces are built in materials order. */
    private static final int    BLUE         = 1;
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

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group root = new Group();

        Octahedron oct = PhiCoordinates.Octahedrons[4];
        final double Z = oct.getEdgeLength() * Math.sqrt(2) / Math.sqrt(3);

        final Jitterbug A = new Jitterbug(oct, materials);
        final Jitterbug M = new Jitterbug(oct, materials);
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

        final double[] prev = { Double.NaN, Double.NaN };
        final long[] last = { 0 };

        final WritableValue<Double> driver = new WritableValue<Double>() {
            @Override
            public Double getValue() {
                return 0d;
            }

            @Override
            public void setValue(Double value) {
                // theta IS the driven blue triangle's own-axis angle. Everything
                // below is solved from it -- no other input exists.
                double theta = value;
                double q = theta;
                double p = theta - PHASE_OFFSET;
                double sep = Z * (Math.cos(Math.toRadians(p)) + Math.cos(Math.toRadians(q)));
                A.rotateTo(q);
                M.rotateTo(p);
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
                // the FAR cell's blue triangle, in world coordinates, this frame
                Point3D far = cen(faces(C, dc).get(BLUE));
                double along = far.dotProduct(n);
                double resp = Double.isNaN(prev[0]) ? Double.NaN : (along - prev[1]) / (theta - prev[0]);
                prev[0] = theta;
                prev[1] = along;
                System.out.printf("DRIVEN blue triangle theta=%+7.2f | solved: A=%+6.2f M=%+6.2f "
                                  + "C=%+6.2f  spacing=%8.5f | far cell along n = %9.5f  "
                                  + "dC/dtheta = %s  lag = 0 (constraint, not propagation)%n",
                                  theta, q, p, q, sep, along,
                                  Double.isNaN(resp) ? "   --   " : String.format("%8.5f", resp));
            }
        };
        driver.setValue(0d);

        final Timeline timeline = new Timeline();
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(12_000), new KeyValue(driver, 60d)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
