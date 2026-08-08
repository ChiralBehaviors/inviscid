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

    /**
     * The phase LUT resolution QUANTUM_RATE is derived from: one frequency
     * quantum advances a member by exactly one step of a LUT with this many
     * steps per revolution. NecronomataVisualization's rotation LUT must be
     * constructed with this resolution (see NecronomataAnimation) so that
     * automaton stepping and rendered rotation stay in lock-step.
     */
    public static final int   PHASE_RESOLUTION = 3600;

    /**
     * The single coupling constant between the conserved {@code frequency}
     * quanta count and angular rate: {@code deltaA[m] == QUANTUM_RATE *
     * frequency[m]} holds immediately AFTER every {@link #step()} call
     * (not "always" unconditionally — before the first step deltaA is
     * still its zero-initialized value, even if frequency was pre-seeded
     * by the caller). Radians of angle advanced per step, per unit of
     * frequency. Derived from {@link #PHASE_RESOLUTION}: 2*pi/3600 ==
     * pi/1800.
     */
    public static final float QUANTUM_RATE = (float) (2.0 * Math.PI
                                                      / PHASE_RESOLUTION);

    /**
     * Member phase, in radians. 30 values per cell (5 cubes x 6 members).
     */
    private final float[] angle;

    /**
     * Angular rate, in radians per step. DERIVED, never independent:
     * recomputed every {@link #step()} as {@code QUANTUM_RATE *
     * frequency[m]}, using frequency AFTER that same step's deltaF has
     * been applied — so a quantum absorbed this tick moves its member on
     * this same tick. Never seeded or written directly.
     */
    private final float[] deltaA;

    /**
     * Transient collision-delta accumulator ONLY. A collision rule writes
     * the quanta it transfers into a member's slot here; {@link #step()}
     * applies it to {@code frequency} and zeroes it, so it never survives
     * past the tick that produced it.
     */
    private final float[] deltaF;
    private final Point3i extent;

    /**
     * Signed integer quanta count, held in a float slot, per member. This
     * is the conserved quantity: members are channels, frequency quanta are
     * the particles that hop between them at collisions.
     */
    private final float[] frequency;

    public Necronomata(float[] angles, Point3i extent, float[] frequency) {
        assert angles.length == 30 * extent.x * extent.y
                                * extent.z : "angles are not correct :"
                                             + angles.length + " in extent: "
                                             + 30 * extent.x * extent.y
                                               * extent.z;
        assert frequency.length == 30 * extent.x * extent.y
                                   * extent.z : "frequencies are not correct in extent :"
                                                + frequency.length
                                                + " in extent: "
                                                + 30 * extent.x * extent.y
                                                  * extent.z;
        this.extent = extent;
        this.angle = angles;
        this.frequency = frequency;
        this.deltaA = new float[angles.length];
        this.deltaF = new float[frequency.length];
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
        return 30 * ((i * extent.y + j) * extent.z + k);
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

    /**
     * Collision-rule visitation site (reserved to bead inviscid-0nx.14).
     * The 12 FCC even-parity neighbor offsets this method will visit are
     * now data, not comment: see
     * {@link com.chiralbehaviors.inviscid.lga.FccNeighborhood} for the
     * offset table, direction indexing, opposite-direction lookup, and
     * wrap-aware {@code neighbor(Point3i, int)} lookup. Note that
     * {@code Necronomata} itself still accepts any extent (odd included)
     * — the even-extent precondition for periodic-wrap parity closure
     * lives on {@code FccNeighborhood}, not here; a caller that wants
     * wrap-safe FCC neighbor lookups must construct
     * {@code FccNeighborhood} with an all-even extent.
     */
    public void process(Point3i cell) {
    }

    /**
     * Advance the automaton one tick. Per member, in order: frequency
     * absorbs whatever a collision rule accumulated in deltaF and deltaF
     * is zeroed (so it enters the next tick empty); deltaA is then
     * recomputed from the now-current frequency (never independent);
     * angle is advanced by deltaA. This ordering is same-tick: a quantum
     * transferred this tick moves its member's angle on this same
     * step(), not the next one.
     */
    public void step() {
        for (int i = 0; i < angle.length; i++) {
            frequency[i] = frequency[i] + deltaF[i];
            deltaF[i] = 0f;
            deltaA[i] = QUANTUM_RATE * frequency[i];
            angle[i] = angle[i] + deltaA[i];
        }
    }

    /**
     * Raw-array escape hatch. Two sanctioned uses: collision-rule
     * visitation writing {@code deltaF} (the conserved-transfer path,
     * bead inviscid-0nx.14) and initial-condition seeding writing
     * {@code frequency} directly (NecronomataAnimation.seedFrequency).
     * Writers of {@code angle}/{@code deltaA} are outside the contract.
     */
    public void process(Processor action) {
        action.process(angle, frequency, deltaA, deltaF);
    }
}
