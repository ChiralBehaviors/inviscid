package util;

import java.util.List;

import javafx.geometry.Point3D;

/**
 * Utility class for data structure related operations.
 *
 * @author Brian Yao
 */
public class Struct {

    /**
     * Deep copy of a list of 3D vectors.
     * 
     * @param vectors
     *            The list of vectors to be copied.
     * @return A copy of the list of vectors.
     */
    public static List<Point3D> copyVectorList(List<Point3D> vectors) {
        return vectors; // Point3D is immutable
    }

}
