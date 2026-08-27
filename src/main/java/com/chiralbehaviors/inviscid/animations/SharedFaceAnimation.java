package com.chiralbehaviors.inviscid.animations;

import static com.chiralbehaviors.inviscid.animations.Colors.materials;

import com.chiralbehaviors.inviscid.Jitterbug;
import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.javafx.experiments.jfx3dviewer.ContentModel;

import javafx.animation.KeyFrame;
import javafx.animation.KeyValue;
import javafx.animation.Timeline;
import javafx.beans.value.WritableValue;
import javafx.geometry.Point3D;
import javafx.scene.Group;
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
        b.getGroup()
         .getTransforms()
         .add(bShift);

        group.getChildren()
             .addAll(a.getGroup(), b.getGroup());
        content.setContent(group);

        final Timeline timeline = new Timeline();
        timeline.getKeyFrames()
                .add(new KeyFrame(Duration.millis(10_000), new KeyValue(new WritableValue<Double>() {
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
                    }
                }, 0d)));
        timeline.setCycleCount(9000);
        timeline.setAutoReverse(true);
        content.setTimeline(timeline);
    }
}
