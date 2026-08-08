/**
 * Copyright (c) 2016 Chiral Behaviors, LLC, all rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.chiralbehaviors.inviscid.lga;

import javax.vecmath.Vector3d;

/**
 * An immutable line segment between two points, in cell-local coordinates.
 * Endpoints are defensively copied on construction and on access, so
 * instances - and repeated calls against the same instance - are safe to
 * share across threads / inner loops without aliasing surprises.
 *
 * @author halhildebrand
 *
 */
public final class Segment {

    private final Vector3d a;
    private final Vector3d b;

    public Segment(Vector3d a, Vector3d b) {
        this.a = new Vector3d(a);
        this.b = new Vector3d(b);
    }

    public Vector3d getA() {
        return new Vector3d(a);
    }

    public Vector3d getB() {
        return new Vector3d(b);
    }

    /**
     * @return the Euclidean distance between the two endpoints.
     */
    public double length() {
        Vector3d diff = new Vector3d(a);
        diff.sub(b);
        return diff.length();
    }
}
