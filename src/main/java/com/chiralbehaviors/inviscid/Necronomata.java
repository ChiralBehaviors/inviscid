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

package com.chiralbehaviors.inviscid;

import java.util.Arrays;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;

import javax.vecmath.Point3i;

/**
 * @author halhildebrand
 *
 */
public class Necronomata implements Iterable<Point3i> {

    @FunctionalInterface
    public interface Processor {
        void process(float[] angle, float[] frequency, float[] deltaA,
                     float[] deltaF);
    }

    private final float[] angle;
    @SuppressWarnings("unused")
    private final float[] deltaA;
    @SuppressWarnings("unused")
    private final float[] deltaF;
    private final Point3i extent;
    @SuppressWarnings("unused")
    private final float[] frequency;

    public Necronomata(float[] angles, Point3i extent, float[] frequency) {
        assert angles.length == 6 * extent.x * extent.y
                                * extent.z : "angles are not correct in extent";
        assert frequency.length == 6 * extent.x * extent.y
                                   * extent.z : "frequencies are not correct in extent";
        this.extent = extent;
        this.angle = angles;
        this.frequency = frequency;
        this.deltaA = angles.clone();
        this.deltaF = frequency.clone();
    }

    public Necronomata(int i, int j, int k) {
        this(new float[30 * i * j * k], new Point3i(i, j, k),
             new float[30 * i * j * k]);
    }

    public Necronomata(Point3i extent) {
        this(extent.x, extent.y, extent.z);
    }

    public float[] anglesOf(Point3i c) {
        int index = indexOfCell(c);
        return Arrays.copyOfRange(angle, index, index + 30);
    }

    public int cellCount() {
        AtomicInteger i = new AtomicInteger();
        forEach(cell -> i.incrementAndGet());
        return i.get();
    }

    @Override
    public void forEach(Consumer<? super Point3i> action) {
        for (int i = 0; i < extent.x; i++) {
            for (int j = 0; j < extent.y; j++) {
                for (int k = 0; k < extent.z; k++) {
                    if ((i + j + k) % 2 == 0) {
                        action.accept(new Point3i(i, j, k));
                    }
                }
            }
        }
    }

    public Point3i getExtent() {
        return new Point3i(extent);
    }

    public int indexOfCell(int i, int j, int k) {
        return i * extent.x + j * extent.y + k * extent.z;
    }

    public int indexOfCell(Point3i cell) {
        return indexOfCell(cell.x, cell.y, cell.z);
    }

    @Override
    public Iterator<Point3i> iterator() {
        return new Iterator<Point3i>() {
            private int i = 0;
            private int j = 0;
            private int k = 0;

            @Override
            public boolean hasNext() {
                return i < extent.x;
            }

            @Override
            public Point3i next() {
                if (k > extent.z) {
                    throw new IllegalStateException("after end of cells");
                }
                Point3i cell = new Point3i(i, j, k);
                increment();
                while (i < extent.x && (i + j + k) % 2 != 0) {
                    increment();
                }
                return cell;
            }

            private void increment() {
                k += 1;
                if (k < extent.z) {
                    return;
                }
                k = 0;
                j += 1;
                if (j < extent.y) {
                    return;
                }
                j = 0;
                i += 1;
            }
        };
    }

    public void process(Point3i cell) {
        /*
        [+1] = { i+1, j  , k-1 }
        [-1] = { i-1, j  , k+1 }
        [+2] = { i  , j-1, k+1 }
        [-2] = { i  , j+1, k-1 }
        [+3] = { i+1, j  , k+1 }
        [-3] = { i-1, j  , k-1 }
        [+4] = { i+1, j+1, k   }
        [-4] = { i-1, j-1, k   }
        [+5] = { i  , j+1, k+1 }
        [-5] = { i  , j-1, k-1 }
        [+6] = { i-1, j+1, k   }
        [-6] = { i+1, j-1, k   }
        */
    }

    public void step() {
        for (int i = 0; i < angle.length; i++) {
            angle[i] = angle[i] + deltaA[i];
            frequency[i] = frequency[i] + deltaF[i];
        }
    }

    public void drive(float[] delta) {
        assert delta.length == 6;
        for (int i = 0; i < angle.length; i++) {
            deltaA[i] = delta[i % 6];
        }
    }

    public void process(Processor action) {
        action.process(angle, frequency, deltaA, deltaF);
    }
}
