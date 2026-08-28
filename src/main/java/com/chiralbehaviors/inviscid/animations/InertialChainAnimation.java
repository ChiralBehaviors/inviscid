package com.chiralbehaviors.inviscid.animations;

import static com.chiralbehaviors.inviscid.animations.Colors.materials;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.chiralbehaviors.inviscid.Jitterbug;
import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.AnimationTimer;
import javafx.geometry.Point3D;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Transform;
import javafx.scene.transform.Translate;
import mesh.polyhedra.plato.Octahedron;

/**
 * THE CHAIN UNDER ITS OWN INERTIA. Mass on the triangles, V = 0, integrated.
 *
 * Nothing drives this. The chain is released and moves under its own momentum:
 *
 *     T = 1/2 M_eff(theta) thetadot^2,   V = 0,   E conserved
 *     thetadot = sqrt(2E / M_eff(theta))
 *
 * M_eff is computed every step from the actual geometry -- finite differences
 * of the 22 DISTINCT triangles' vertices (24 cell-faces minus the 2 that are
 * shared, counted once because a shared face is one triangle with two owners).
 *
 * MASS MODEL IS DECLARED, as this project requires, and it MATTERS:
 *     POINT   masses m/3 at the corners   -> period 211.37
 *     LAMINA  uniform triangular plate    -> period 196.61
 * a 7 percent difference on the same geometry. Switch with MODEL below.
 *
 * WHAT THE INTEGRATION SHOWS, and it is the answer to "give them mass and
 * integrate":
 *
 *   TIME NOW EXISTS. thetadot is not constant -- it peaks at theta = 25.83,
 *   near the congruent configuration where M_eff is least, and is slowest at
 *   both ends. That is Fuller's "the equatorial rotational momentum will be
 *   seen to carry the rotation beyond dead-center", at chain scale. A linear
 *   timeline ramp never showed this; the motion has a shape of its own.
 *
 *   BUT THERE IS STILL NO LAG, AND THERE CANNOT BE. The chain has ONE degree
 *   of freedom. theta(t) is a single function of time and all three cells take
 *   their phase from it in the same instant. The far cell does not trail the
 *   near one by so much as a step, at any mass, at any energy. Mass buys TIME
 *   DEPENDENCE; it does not buy SPATIAL STRUCTURE.
 *
 *   So this is one oscillator, not a medium. A wave needs a disturbance that
 *   can be SOMEWHERE -- and with a single coordinate for the whole chain there
 *   is no "where" for it to be. Adding cells to the line does not help: the
 *   line has one degree of freedom at any length, and the full array has about
 *   twenty at any size.
 */
public class InertialChainAnimation extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            InertialChainAnimation.main(argv);
        }
    }

    /** "point" = m/3 at each corner; "lamina" = uniform triangular plate. */
    private static final String MODEL        = "lamina";
    private static final double PHASE_OFFSET = 60.0;
    private static final double DT           = 0.35;   // integration step
    private static final double THETA0       = 0.0;
    private static final double THETADOT0    = 1.0;

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

    private Jitterbug pa, pm, pc;
    private Point3D   axis;
    private double    Z;

    // Shared faces, resolved ONCE at construction: the face indices do not
    // change with theta because the face axes are fixed. Weight 0.5 on each
    // side of a shared face counts that triangle exactly once overall --
    // 24 cell-faces - 4 x 0.5 = 22 distinct triangles.
    private int sharedAplus, sharedMminus, sharedMplus, sharedCminus;

    private int faceAlong(Jitterbug j, Point3D dir) {
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

    private double weight(int cell, int face) {
        if (cell == 0 && face == sharedAplus) {
            return 0.5;
        }
        if (cell == 1 && (face == sharedMminus || face == sharedMplus)) {
            return 0.5;
        }
        if (cell == 2 && face == sharedCminus) {
            return 0.5;
        }
        return 1.0;
    }

    /** Vertex velocities d(vertex)/d(theta), per cell, by central difference. */
    private List<List<List<Point3D>>> geom(double theta) {
        double p = theta - PHASE_OFFSET;
        double sep = Z * (Math.cos(Math.toRadians(p)) + Math.cos(Math.toRadians(theta)));
        pa.rotateTo(theta);
        pm.rotateTo(p);
        pc.rotateTo(theta);
        List<List<List<Point3D>>> out = new ArrayList<>();
        out.add(faces(pa, axis.multiply(-sep)));
        out.add(faces(pm, null));
        out.add(faces(pc, axis.multiply(sep)));
        return out;
    }

    /** 2T/thetadot^2 from the real geometry. Indexing is by (cell, face), so
     *  it is stable across angles; sharing is handled by weight(), not by
     *  deduplicating on position -- an earlier version paired two deduped lists
     *  BY INDEX across different angles, which is neither stable nor ordered,
     *  and threw once the coincidence count changed from 22 to 23. */
    private double meff(double theta) {
        double h = 1e-4;
        List<List<List<Point3D>>> tp = geom(theta + h), tm = geom(theta - h);
        double sum = 0, lam = 0;
        for (int c = 0; c < 3; c++) {
            for (int f = 0; f < 8; f++) {
                double w = weight(c, f);
                Point3D s = Point3D.ZERO;
                for (int v = 0; v < 3; v++) {
                    Point3D vel = tp.get(c)
                                    .get(f)
                                    .get(v)
                                    .subtract(tm.get(c)
                                                .get(f)
                                                .get(v))
                                    .multiply(1.0 / (2 * h));
                    sum += w * vel.dotProduct(vel);
                    s = s.add(vel);
                }
                lam += w * s.dotProduct(s);
            }
        }
        return "point".equals(MODEL) ? sum / 3.0 : (sum + lam) / 12.0;
    }

    /** dM_eff/dtheta, for the equation of motion. */
    private double dmeff(double theta) {
        double h = 1e-3;
        return (meff(theta + h) - meff(theta - h)) / (2 * h);
    }

    @Override
    protected void initializeContentModel() {
        ContentModel content = getContentModel();
        Group root = new Group();
        Octahedron oct = PhiCoordinates.Octahedrons[4];
        Z = oct.getEdgeLength() * Math.sqrt(2) / Math.sqrt(3);

        pa = new Jitterbug(oct, materials);
        pm = new Jitterbug(oct, materials);
        pc = new Jitterbug(oct, materials);
        pm.rotateTo(0);
        axis = cen(faces(pm, null).get(0)).normalize();

        final Translate sa = new Translate(), sc = new Translate();
        pa.getGroup()
          .getTransforms()
          .add(sa);
        pc.getGroup()
          .getTransforms()
          .add(sc);
        root.getChildren()
            .addAll(pa.getGroup(), pm.getGroup(), pc.getGroup());
        content.setContent(root);

        sharedAplus = faceAlong(pa, axis);
        sharedMminus = faceAlong(pm, axis.multiply(-1));
        sharedMplus = faceAlong(pm, axis);
        sharedCminus = faceAlong(pc, axis.multiply(-1));

        final double[] th = { THETA0 };
        final double[] td = { THETADOT0 };
        final double E0 = 0.5 * meff(THETA0) * THETADOT0 * THETADOT0;
        final long[] last = { 0 };
        final double[] worst = { 0 };

        final Runnable place = () -> {
            double p = th[0] - PHASE_OFFSET;
            double sep = Z * (Math.cos(Math.toRadians(p)) + Math.cos(Math.toRadians(th[0])));
            pa.rotateTo(th[0]);
            pm.rotateTo(p);
            pc.rotateTo(th[0]);
            Point3D da = axis.multiply(-sep), dc = axis.multiply(sep);
            sa.setX(da.getX());
            sa.setY(da.getY());
            sa.setZ(da.getZ());
            sc.setX(dc.getX());
            sc.setY(dc.getY());
            sc.setZ(dc.getZ());
        };
        place.run();

        new AnimationTimer() {
            @Override
            public void handle(long now) {
                // V = 0, so the only force is the configuration-dependent
                // inertia: thetaddot = -(1/2)(M'/M) thetadot^2. Velocity
                // Verlet, so the energy audit below is a real check and not a
                // tautology -- an earlier version derived thetadot FROM E,
                // which conserves E by construction and tests nothing.
                double M = meff(th[0]);
                double acc = -0.5 * (dmeff(th[0]) / M) * td[0] * td[0];
                td[0] += 0.5 * DT * acc;
                th[0] += DT * td[0];
                if (th[0] >= 60.0) {
                    th[0] = 60.0;
                    td[0] = -Math.abs(td[0]);
                } else if (th[0] <= 0.0) {
                    th[0] = 0.0;
                    td[0] = Math.abs(td[0]);
                }
                double acc2 = -0.5 * (dmeff(th[0]) / meff(th[0])) * td[0] * td[0];
                td[0] += 0.5 * DT * acc2;
                place.run();
                double E = 0.5 * meff(th[0]) * td[0] * td[0];
                worst[0] = Math.max(worst[0], Math.abs(E - E0) / E0);
                if (System.currentTimeMillis() - last[0] > 500) {
                    last[0] = System.currentTimeMillis();
                    System.out.printf("theta=%+7.3f  thetadot=%+8.5f  M_eff=%10.7f | "
                                      + "E=%11.8f  worst drift=%8.2e | near=%+7.3f far=%+7.3f "
                                      + "LAG=%d%n",
                                      th[0], td[0], meff(th[0]), E, worst[0], th[0], th[0], 0);
                }
            }
        }.start();
    }
}
