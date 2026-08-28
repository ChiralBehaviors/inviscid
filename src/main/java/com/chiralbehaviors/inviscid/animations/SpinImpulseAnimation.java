package com.chiralbehaviors.inviscid.animations;

import static com.chiralbehaviors.inviscid.animations.Colors.materials;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.chiralbehaviors.inviscid.Jitterbug;
import com.chiralbehaviors.inviscid.PhiCoordinates;

import javafx.animation.AnimationTimer;
import javafx.geometry.Point3D;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.paint.Color;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.CullFace;
import javafx.scene.shape.MeshView;
import javafx.scene.shape.TriangleMesh;
import javafx.scene.transform.Transform;
import mesh.polyhedra.plato.Octahedron;

/**
 * SPIN ONE TRIANGLE OF THE RIGHT OCTAHEDRON, THEN STOP.  octa - VE - octa.
 *
 * The port of jb_ic_inertial_chain.py, and it exists because the previous
 * animation was built on machinery that is wrong in a way Fuller named in 1977.
 *
 * CONGRUENCE IS NOT IDENTITY. "Deceptiveness of Topology -- Quanta Lost by
 * Congruence": the jitterbug has 24 edges and 12 vertices at EVERY
 * configuration; at the octahedron they are "24 EDGES CONGRUENT AS 12, 12
 * VERTICES CONGRUENT AS 6". Deduplicating by position -- which every earlier
 * builder here did -- is exactly the accounting he objects to, and it costs an
 * octahedral cell HALF ITS INERTIA while welding its congruent vertices into a
 * solid that can never come apart. So all 24 struts and 8 triangles per cell are
 * carried as identity objects; only genuine shared-face vertices are welded.
 *
 * TWO CONSTRUCTION TRAPS, both of which produced wrong answers before being
 * caught:
 *   * the within-cell slot->vertex map must be read at a GENERIC angle; taken
 *     at gamma = 60 it sees 6 vertices per cell instead of 12.
 *   * the shared-face corner correspondence FLIPS with the direction of the 60
 *     degree offset. A chain alternates, so every second weld needs the other
 *     one; a single correspondence throughout gives triangles whose struts are
 *     not sqrt(2).
 *
 * THE IMPULSE is a SPIN -- angular velocity about the blue triangle's own axis,
 * which is the jitterbug's actual degree of freedom -- projected onto the
 * constraint tangent space so it is admissible, then released. V = 0 throughout;
 * energy is the only audit and RATTLE has to earn it.
 *
 * WHAT TO WATCH. Every cell lights up at t = 0. A rigid constraint has infinite
 * signal speed, so there is no onset lag anywhere: the projection that makes the
 * impulse admissible reaches the whole chain at once. What takes time is the
 * AMPLITUDE -- the far cell's share of the kinetic energy grows by an order of
 * magnitude. This medium transports, but it has no wavefront, and any apparent
 * front is mode superposition.
 */
public class SpinImpulseAnimation extends PolyView {
    public static class Launcher {
        public static void main(String[] argv) {
            SpinImpulseAnimation.main(argv);
        }
    }

    private static final double[] CFG  = { 60.0, 0.0, 60.0 };   // octa - VE - octa
    private static final double   GEN  = 30.0;                  // generic angle
    private static final double   H    = 0.005;
    private static final int      SUB  = 4;
    private static final double   TOL  = 1e-12;

    public static void main(String[] args) {
        launch(args);
    }

    private int      nv, ncell;
    private int[]    slot;                 // (f*3+c) -> within-cell vertex id
    private double[] q, v, mass;
    private int[]    bi, bj;
    private double[] L2;
    private int[][]  tri;
    private int[]    triCell;
    private int[][]  cellVerts;
    private Point3D  nhat;
    private double   ZC, EL;

    private static List<List<Point3D>> faces(Jitterbug j) {
        List<List<Point3D>> out = new ArrayList<>();
        for (Node ch : j.getGroup()
                        .getChildren()) {
            Group fg = (Group) ch;
            TriangleMesh tm = (TriangleMesh) ((MeshView) fg.getChildren()
                                                           .get(0)).getMesh();
            float[] p = new float[tm.getPoints()
                                    .size()];
            tm.getPoints()
              .toArray(p);
            List<Point3D> f = new ArrayList<>();
            for (int i = 0; i < p.length; i += 3) {
                Point3D t = new Point3D(p[i], p[i + 1], p[i + 2]);
                for (Transform x : fg.getTransforms()) {
                    t = x.transform(t);
                }
                f.add(t);
            }
            out.add(f);
        }
        return out;
    }

    private static String key(Point3D p) {
        return String.format("%.6f,%.6f,%.6f", p.getX(), p.getY(), p.getZ());
    }

    /** 12 vertex positions of one cell at angle g, using the FIXED topology. */
    private Point3D[] cellVertsAt(double g, Point3D origin, Octahedron oct) {
        Jitterbug j = new Jitterbug(oct, materials);
        j.rotateTo(g);
        List<List<Point3D>> fs = faces(j);
        Point3D[] V = new Point3D[nv];
        for (int f = 0; f < 8; f++) {
            for (int c = 0; c < 3; c++) {
                int i = slot[f * 3 + c];
                if (V[i] == null) {
                    V[i] = fs.get(f)
                             .get(c)
                             .add(origin);
                }
            }
        }
        return V;
    }

    /** Corner correspondence across a shared face. DIRECTION-DEPENDENT. */
    private int[][] weldFor(double ga, double gb, int fp, Octahedron oct) {
        double sep = ZC * (Math.cos(Math.toRadians(ga)) + Math.cos(Math.toRadians(gb)));
        Point3D[] A = cellVertsAt(ga, Point3D.ZERO, oct);
        Point3D[] B = cellVertsAt(gb, nhat.multiply(sep), oct);
        int[][] w = new int[3][2];
        for (int c = 0; c < 3; c++) {
            int a = slot[fp * 3 + c];
            int best = 0;
            double bd = Double.MAX_VALUE;
            for (int k = 0; k < nv; k++) {
                double d = B[k].distance(A[a]);
                if (d < bd) {
                    bd = d;
                    best = k;
                }
            }
            w[c] = new int[] { a, best };
        }
        return w;
    }

    private void build() {
        Octahedron oct = PhiCoordinates.Octahedrons[4];
        EL = oct.getEdgeLength();
        ZC = EL * Math.sqrt(2) / Math.sqrt(3);
        ncell = CFG.length;

        // within-cell topology from a GENERIC angle
        Jitterbug g0 = new Jitterbug(oct, materials);
        g0.rotateTo(GEN);
        List<List<Point3D>> gf = faces(g0);
        // Match shared corners by TOLERANCE, not by an exact key. In Java each
        // face's corner is produced through its own transform chain, so two
        // faces' copies of the same vertex differ in the last bits; a string
        // key at 1e-6 fails to match them and the cell comes out with 24
        // vertices instead of 12, leaving the chain unconnected. (Python's
        // Z.corners returns them bit-identical, which is why the port exposed
        // this and the original never did.)
        List<Point3D> reps = new ArrayList<>();
        slot = new int[24];
        for (int f = 0; f < 8; f++) {
            for (int c = 0; c < 3; c++) {
                Point3D p = gf.get(f)
                              .get(c);
                int hit = -1;
                for (int k = 0; k < reps.size(); k++) {
                    if (reps.get(k)
                            .distance(p) < 1e-6) {
                        hit = k;
                        break;
                    }
                }
                if (hit < 0) {
                    hit = reps.size();
                    reps.add(p);
                }
                slot[f * 3 + c] = hit;
            }
        }
        nv = reps.size();
        if (nv != 12) {
            throw new IllegalStateException("within-cell topology gave " + nv
                                            + " vertices, want 12 -- read it at a GENERIC angle");
        }
        nhat = cen(gf.get(0)).normalize();
        int fp = 0, fm = 0;
        double bp = -2, bm = 2;
        for (int f = 0; f < 8; f++) {
            double d = cen(gf.get(f)).normalize()
                                     .dotProduct(nhat);
            if (d > bp) {
                bp = d;
                fp = f;
            }
            if (d < bm) {
                bm = d;
                fm = f;
            }
        }

        double[] origAlong = new double[ncell];
        for (int k = 1; k < ncell; k++) {
            origAlong[k] = origAlong[k - 1]
                           + ZC * (Math.cos(Math.toRadians(CFG[k - 1])) + Math.cos(Math.toRadians(CFG[k])));
        }
        int[] par = new int[ncell * nv];
        for (int i = 0; i < par.length; i++) {
            par[i] = i;
        }
        for (int k = 0; k + 1 < ncell; k++) {
            for (int[] ab : weldFor(CFG[k], CFG[k + 1], fp, oct)) {
                union(par, k * nv + ab[0], (k + 1) * nv + ab[1]);
            }
        }
        Map<Integer, Integer> uniq = new HashMap<>();
        int[][] gid = new int[ncell][nv];
        for (int k = 0; k < ncell; k++) {
            for (int i = 0; i < nv; i++) {
                gid[k][i] = uniq.computeIfAbsent(find(par, k * nv + i), x -> uniq.size());
            }
        }
        int n = uniq.size();
        q = new double[3 * n];
        v = new double[3 * n];
        mass = new double[n];
        double worst = 0;
        boolean[] wr = new boolean[n];
        for (int k = 0; k < ncell; k++) {
            Point3D[] V = cellVertsAt(CFG[k], nhat.multiply(origAlong[k]), oct);
            for (int i = 0; i < nv; i++) {
                int g = gid[k][i];
                if (wr[g]) {
                    worst = Math.max(worst,
                                     new Point3D(q[3 * g], q[3 * g + 1], q[3 * g + 2]).distance(V[i]));
                }
                q[3 * g] = V[i].getX();
                q[3 * g + 1] = V[i].getY();
                q[3 * g + 2] = V[i].getZ();
                wr[g] = true;
            }
        }
        // ALL 8 triangles and 24 struts per cell, kept even where congruent
        tri = new int[8 * ncell][3];
        triCell = new int[8 * ncell];
        int t = 0;
        for (int k = 0; k < ncell; k++) {
            for (int f = 0; f < 8; f++) {
                for (int c = 0; c < 3; c++) {
                    tri[t][c] = gid[k][slot[f * 3 + c]];
                }
                triCell[t] = k;
                t++;
            }
        }
        int nb = 3 * tri.length;
        bi = new int[nb];
        bj = new int[nb];
        L2 = new double[nb];
        int r = 0;
        for (int[] tt : tri) {
            int[][] e = { { tt[0], tt[1] }, { tt[1], tt[2] }, { tt[0], tt[2] } };
            for (int[] ee : e) {
                bi[r] = ee[0];
                bj[r] = ee[1];
                L2[r] = d2(ee[0], ee[1]);
                r++;
            }
            for (int i : tt) {
                mass[i] += 1.0 / 3.0;
            }
        }
        cellVerts = new int[ncell][];
        for (int k = 0; k < ncell; k++) {
            java.util.TreeSet<Integer> s = new java.util.TreeSet<>();
            for (int i = 0; i < nv; i++) {
                s.add(gid[k][i]);
            }
            int[] a = new int[s.size()];
            int w = 0;
            for (int x : s) {
                a[w++] = x;
            }
            cellVerts[k] = a;
        }
        double lmin = Double.MAX_VALUE, lmax = 0;
        for (int b = 0; b < nb; b++) {
            double L = Math.sqrt(d2(bi[b], bj[b]));
            lmin = Math.min(lmin, L);
            lmax = Math.max(lmax, L);
        }
        System.out.printf("octa-VE-octa: %d verts, %d struts (24 x %d), %d triangles (8 x %d)%n",
                          n, nb, ncell, tri.length, ncell);
        System.out.printf("   struts %.6f..%.6f (want %.6f) | shared-write disagreement %.1e%n",
                          lmin, lmax, EL, worst);
    }

    private static int find(int[] p, int x) {
        while (p[x] != x) {
            p[x] = p[p[x]];
            x = p[x];
        }
        return x;
    }

    private static void union(int[] p, int a, int b) {
        int ra = find(p, a), rb = find(p, b);
        if (ra != rb) {
            p[ra] = rb;
        }
    }

    private static Point3D cen(List<Point3D> f) {
        double x = 0, y = 0, z = 0;
        for (Point3D p : f) {
            x += p.getX();
            y += p.getY();
            z += p.getZ();
        }
        return new Point3D(x / f.size(), y / f.size(), z / f.size());
    }

    private double d2(int i, int j) {
        double x = q[3 * i] - q[3 * j], y = q[3 * i + 1] - q[3 * j + 1], z = q[3 * i + 2] - q[3 * j + 2];
        return x * x + y * y + z * z;
    }

    // ---- RATTLE. Gauss-Seidel projection, which tolerates the redundant
    // ---- constraints the congruent struts introduce.
    private void shake(double[] q0) {
        for (int it = 0; it < 400; it++) {
            double worst = 0;
            for (int b = 0; b < bi.length; b++) {
                int i = bi[b], j = bj[b];
                double rx = q[3 * i] - q[3 * j], ry = q[3 * i + 1] - q[3 * j + 1],
                        rz = q[3 * i + 2] - q[3 * j + 2];
                double g = rx * rx + ry * ry + rz * rz - L2[b];
                worst = Math.max(worst, Math.abs(g));
                double ox = q0[3 * i] - q0[3 * j], oy = q0[3 * i + 1] - q0[3 * j + 1],
                        oz = q0[3 * i + 2] - q0[3 * j + 2];
                double dot = rx * ox + ry * oy + rz * oz;
                if (Math.abs(dot) < 1e-14) {
                    continue;
                }
                double im = 1.0 / mass[i] + 1.0 / mass[j];
                double lam = -g / (4.0 * dot * im);
                apply(i, j, 2 * lam * ox, 2 * lam * oy, 2 * lam * oz, true);
            }
            if (worst < TOL) {
                break;
            }
        }
    }

    private void apply(int i, int j, double dx, double dy, double dz, boolean alsoV) {
        q[3 * i] += dx / mass[i];
        q[3 * i + 1] += dy / mass[i];
        q[3 * i + 2] += dz / mass[i];
        q[3 * j] -= dx / mass[j];
        q[3 * j + 1] -= dy / mass[j];
        q[3 * j + 2] -= dz / mass[j];
        if (alsoV) {
            v[3 * i] += dx / mass[i] / H;
            v[3 * i + 1] += dy / mass[i] / H;
            v[3 * i + 2] += dz / mass[i] / H;
            v[3 * j] -= dx / mass[j] / H;
            v[3 * j + 1] -= dy / mass[j] / H;
            v[3 * j + 2] -= dz / mass[j] / H;
        }
    }

    private void projectV() {
        for (int it = 0; it < 400; it++) {
            double worst = 0;
            for (int b = 0; b < bi.length; b++) {
                int i = bi[b], j = bj[b];
                double rx = q[3 * i] - q[3 * j], ry = q[3 * i + 1] - q[3 * j + 1],
                        rz = q[3 * i + 2] - q[3 * j + 2];
                double dv = rx * (v[3 * i] - v[3 * j]) + ry * (v[3 * i + 1] - v[3 * j + 1])
                            + rz * (v[3 * i + 2] - v[3 * j + 2]);
                worst = Math.max(worst, Math.abs(dv));
                double im = 1.0 / mass[i] + 1.0 / mass[j];
                double mu = -dv / (2.0 * (rx * rx + ry * ry + rz * rz) * im);
                v[3 * i] += 2 * mu * rx / mass[i];
                v[3 * i + 1] += 2 * mu * ry / mass[i];
                v[3 * i + 2] += 2 * mu * rz / mass[i];
                v[3 * j] -= 2 * mu * rx / mass[j];
                v[3 * j + 1] -= 2 * mu * ry / mass[j];
                v[3 * j + 2] -= 2 * mu * rz / mass[j];
            }
            if (worst < TOL) {
                break;
            }
        }
    }

    private void step() {
        double[] q0 = q.clone();
        for (int i = 0; i < q.length; i++) {
            q[i] += H * v[i];
        }
        shake(q0);
        projectV();
    }

    private double[] cellKE() {
        double[] ke = new double[ncell];
        for (int k = 0; k < ncell; k++) {
            for (int i : cellVerts[k]) {
                ke[k] += 0.5 * mass[i]
                         * (v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] + v[3 * i + 2] * v[3 * i + 2]);
            }
        }
        return ke;
    }

    private double energy() {
        double e = 0;
        for (int i = 0; i < mass.length; i++) {
            e += 0.5 * mass[i]
                 * (v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] + v[3 * i + 2] * v[3 * i + 2]);
        }
        return e;
    }

    @Override
    protected void initializeContentModel() {
        build();

        // SPIN the blue triangle (face 1) of the RIGHT octahedron about its own axis
        int blue = 8 * (ncell - 1) + 1;
        double cx = 0, cy = 0, cz = 0;
        for (int c = 0; c < 3; c++) {
            cx += q[3 * tri[blue][c]] / 3;
            cy += q[3 * tri[blue][c] + 1] / 3;
            cz += q[3 * tri[blue][c] + 2] / 3;
        }
        double kx = 0, ky = 0, kz = 0;
        for (int i : cellVerts[ncell - 1]) {
            kx += q[3 * i];
            ky += q[3 * i + 1];
            kz += q[3 * i + 2];
        }
        int m = cellVerts[ncell - 1].length;
        Point3D ax = new Point3D(cx - kx / m, cy - ky / m, cz - kz / m).normalize();
        for (int c = 0; c < 3; c++) {
            int i = tri[blue][c];
            Point3D rel = new Point3D(q[3 * i] - cx, q[3 * i + 1] - cy, q[3 * i + 2] - cz);
            Point3D w = ax.crossProduct(rel);
            v[3 * i] = w.getX();
            v[3 * i + 1] = w.getY();
            v[3 * i + 2] = w.getZ();
        }
        projectV();
        final double E0 = energy();
        System.out.printf("   SPIN impulse on triangle %d (blue) of the right octahedron; "
                          + "E = %.6f, V = 0%n%n", blue, E0);

        Group root = new Group();
        final MeshView[] mv = new MeshView[tri.length];
        final PhongMaterial[] mat = new PhongMaterial[ncell];
        for (int k = 0; k < ncell; k++) {
            mat[k] = new PhongMaterial(Color.DARKSLATEBLUE);
        }
        for (int t = 0; t < tri.length; t++) {
            TriangleMesh msh = new TriangleMesh();
            msh.getPoints()
               .addAll(0, 0, 0, 0, 0, 0, 0, 0, 0);
            msh.getTexCoords()
               .addAll(0, 0);
            msh.getFaces()
               .addAll(0, 0, 1, 0, 2, 0);
            mv[t] = new MeshView(msh);
            mv[t].setCullFace(CullFace.NONE);
            mv[t].setMaterial(mat[triCell[t]]);
            root.getChildren()
                .add(mv[t]);
        }
        getContentModel().setContent(root);

        final Runnable draw = () -> {
            for (int t = 0; t < tri.length; t++) {
                float[] p = new float[9];
                for (int c = 0; c < 3; c++) {
                    p[3 * c] = (float) q[3 * tri[t][c]];
                    p[3 * c + 1] = (float) q[3 * tri[t][c] + 1];
                    p[3 * c + 2] = (float) q[3 * tri[t][c] + 2];
                }
                ((TriangleMesh) mv[t].getMesh()).getPoints()
                                                .setAll(p);
            }
            double[] ke = cellKE();
            double tot = 0;
            for (double x : ke) {
                tot += x;
            }
            for (int k = 0; k < ncell; k++) {
                double f = Math.min(1.0, ke[k] / tot * ncell / 1.5);
                mat[k].setDiffuseColor(Color.color(0.15 + 0.85 * f, 0.15 + 0.2 * f, 0.6 - 0.5 * f));
            }
        };
        draw.run();

        final long[] last = { 0 };
        final double[] ts = { 0 };
        new AnimationTimer() {
            @Override
            public void handle(long now) {
                for (int s = 0; s < SUB; s++) {
                    step();
                    ts[0] += H;
                }
                draw.run();
                if (System.currentTimeMillis() - last[0] > 500) {
                    last[0] = System.currentTimeMillis();
                    double[] ke = cellKE();
                    double tot = 0;
                    for (double x : ke) {
                        tot += x;
                    }
                    System.out.printf("t=%7.2f | cell0 %7.5f  cell1 %7.5f  cell2 %7.5f "
                                      + "| E=%.6f  drift=%8.2e%n",
                                      ts[0], ke[0] / tot, ke[1] / tot, ke[2] / tot, energy(),
                                      Math.abs(energy() - E0) / E0);
                }
            }
        }.start();
    }
}
