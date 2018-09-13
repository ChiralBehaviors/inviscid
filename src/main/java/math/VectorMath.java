package math;

import javafx.geometry.Point3D;

/**
 * A class with helper vector functions.
 *
 * @author Brian Yao
 */
public class VectorMath {

    private static final double EPSILON = 1e-15;

    /**
     * Returns true if the given vector is the zero vector, within a small
     * margin of error (component-wise).
     *
     * @param vector
     *            The vector we are comparing to zero.
     * @return True if the vector is approximately zero.
     */
    public static boolean isZero(Point3D vector) {
        return Math.abs(vector.getX()) < EPSILON
               && Math.abs(vector.getY()) < EPSILON
               && Math.abs(vector.getZ()) < EPSILON;
    }

    /**
     * Returns true if the given vector has at least one component which is NaN.
     *
     * @param vector
     *            The vector whose components we are comparing to NaN.
     * @return True if at least one component is NaN.
     */
    public static boolean isNaN(Point3D vector) {
        return Double.isNaN(vector.getX()) || Double.isNaN(vector.getY())
               || Double.isNaN(vector.getZ());
    }

    /**
     * Returns the parametric interpolation of two vectors for parameter 0 <= t
     * <= 1. At t = 0 this returns the start, and at t = 1 returns the end. In
     * general, this returns start + (end - start) * t.
     *
     * @param start
     *            The start vector (t = 0).
     * @param end
     *            The end vector (t = 1).
     * @param t
     *            The parameter (0 <= t <= 1).
     * @return The interpolation at the parameter t.
     */
    public static Point3D interpolate(Point3D start, Point3D end, double t) {
        return start.add(start.subtract(end)
                              .multiply(t));
    }

}
