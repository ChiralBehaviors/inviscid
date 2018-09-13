package util;

import javafx.geometry.Point3D;

/**
 * A utility class involving geometric computations.
 *
 * @author Brian Yao
 */
public class Geometry {

    /**
     * Computes the square of the distance from a point to a line in 3D.
     * 
     * @param point
     *            The point we are calculating distance from.
     * @param linePoint1
     *            The first of two points defining the line.
     * @param linePoint2
     *            The second of two points defining the line.
     * @return The square of the distance of point from the line defined by the
     *         two line points.
     */
    public static double pointLineDistSq(Point3D point, Point3D linePoint1,
                                         Point3D linePoint2) {
        Point3D lineVector = linePoint2.subtract(linePoint1);
        Point3D linePointVector = linePoint1.subtract(point);
        
        double magnitude = lineVector.magnitude();
        double lineVectorLengthSq = magnitude * magnitude;

        if (lineVectorLengthSq == 0) {
            double magnitude2 = linePointVector.magnitude();
            return magnitude2 * magnitude2;
        } else {
            Point3D scaledLineVector = lineVector.multiply(linePointVector.dotProduct(lineVector));
            Point3D distanceVector = linePointVector.subtract(scaledLineVector);
            double magnitude2 = distanceVector.magnitude();
            return magnitude2 * magnitude2;
        }
    }

    /**
     * Computes the distance from a point to a line in 3D.
     * 
     * @param point
     *            The point we are calculating distance from.
     * @param linePoint1
     *            The first of two points defining the line.
     * @param linePoint2
     *            The second of two points defining the line.
     * @return The distance of point from the line defined by the two line
     *         points.
     */
    public static double pointLineDist(Point3D point, Point3D linePoint1,
                                       Point3D linePoint2) {
        return Math.sqrt(pointLineDistSq(point, linePoint1, linePoint2));
    }

}
