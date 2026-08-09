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

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.QuantaField;
import com.chiralbehaviors.inviscid.measure.CollisionStatistics;

/**
 * The formal lattice-gas automaton (bead inviscid-0nx.21, C.4): synchronous
 * update on the even-parity sublattice using {@link FineStepContactTable} +
 * {@link CollisionTable}. No run-time geometry, no run-time trig, no
 * floating-point contact test IN THE TICK PATH - contact firing is a
 * bitset lookup (built once, at construction, from geometry - see "Contact
 * firing" below), quanta transfer is a frozen table lookup, and the only
 * floating-point value anywhere near a live tick is {@link #phaseAt(int)}'s
 * MEASUREMENT read, which lives entirely outside the tick path (see that
 * method's Javadoc).
 *
 * <h2>Contact firing (USER DECISION 2026-08-08, FINAL, revised)</h2>
 * Contact FIRING is decided by {@link FineStepContactTable}: the atlas's
 * own fine {@code geometryResolution} (360) geometry grid, replicating the
 * hybrid's live {@code ContactPredicate} point-evaluation semantics at that
 * resolution - built once at construction (reusing {@link
 * ContactComboCache} + {@link ContactPredicate}, the exact machinery that
 * built the atlas), not per-tick. This SUPERSEDES an earlier revision of
 * this class that fired on {@link ContactTable}'s bin-level ANY-OVERLAP
 * bitset directly; that approach was measured (bead inviscid-0nx.21's test
 * 8) to produce a ~6.6x collision-rate / ~3.9x effective-ratio gap against
 * the hybrid, root-caused to ANY-OVERLAP's deliberate over-inclusiveness
 * (bead inviscid-gyt) - a property of the BIN table's transcription
 * semantics, not a bug. gyt and bead inviscid-0nx.20's collision-table
 * decision (D1) are NOT overturned: {@link ContactTable} (bead
 * inviscid-0nx.19) and {@link CollisionTable} (bead inviscid-0nx.20) stand
 * UNCHANGED as the key/outcome layer - {@code ContactTable} is still
 * constructed and held (see "Header pairing" below) but no longer consulted
 * for contact-firing; a bin-level pre-filter was considered and rejected
 * (see {@link FineStepContactTable}'s own Javadoc, "No bin-level
 * pre-filter"). {@code accumulator / 10} ({@link #fineStepDivisor}, guarded
 * loudly at construction) converts the phase accumulator to a fine step,
 * exact integer arithmetic since {@code phaseResolution (3600) ==
 * geometryResolution (360) * 10} exactly.
 *
 * <h2>Phase state (user's 2A cadence decision, FINAL, 2026-08-08)</h2>
 * Per member: an integer accumulator {@code phase} in {@code [0, 3600)}
 * (never a bare {@code N_lga}-resolution counter) plus the conserved
 * {@code long quanta}. {@code M} sub-steps ({@link ContactAtlas.Header#subBinSteps()},
 * sourced from exactly one place - {@link ContactAtlasGenerator#SUB_BIN_STEPS}
 * - never {@code Necronomata.PHASE_RESOLUTION}) span one contact bin:
 * {@code bin = phase / M}. One LGA tick advances phase by exactly one
 * member's {@code quanta}, mod 3600 - the same rate a hybrid member's
 * angle advances under {@code Necronomata.step()} ({@code
 * deltaA == QUANTUM_RATE * frequency}, and {@code QUANTUM_RATE == 2*pi /
 * 3600}) - so one LGA tick == one hybrid tick (T2
 * analysis-73v-spectral-conversion-and-cadence.md §3.2, option 2A).
 * {@code N_lga} ({@link ContactTable#nLga()}) keys contact/collision table
 * lookups ONLY - it is never the phase state alphabet (73v's "Reading K",
 * not "Reading P" - see that memo's §3.1 fork).
 *
 * <h2>Layout contract</h2>
 * Mirrors {@code Necronomata}/{@link QuantaField}'s layout exactly: 30
 * slots per cell (5 cubes x 6 members), linearised as {@code 30 *
 * ((i*extent.y+j)*extent.z+k)}, even-parity cells only ({@code
 * (i+j+k)%2==0}) - this is what lets every Phase B instrument
 * ({@code ConservationAudit}, {@code StructureFactor.coarseGrainedField})
 * address this substrate identically to {@code Necronomata}.
 *
 * <h2>Synchronous update</h2>
 * {@link #tick(int)} is a strict scan/apply split: the scan phase reads
 * ONLY the pre-tick {@code phase}/{@code quanta} arrays (never mutating
 * them) and accumulates every resolved transfer into a fresh {@code
 * deltaQuanta} array; the apply phase, which runs only after the FULL
 * scan completes, adds {@code deltaQuanta} into {@code quanta} and then
 * advances {@code phase} by the POST-collision {@code quanta} (same-tick
 * absorption, matching {@code Necronomata.step()}'s own convention - a
 * quantum absorbed this tick moves its member this same tick). No cell
 * ever sees a partially-updated lattice, and the result is therefore
 * independent of cell-visitation order - see {@link #cellVisitOrder()}.
 *
 * <h2>Header pairing (bead inviscid-0nx.19's contract)</h2>
 * The constructor receives one {@link ContactTable} (built from a
 * {@link ContactAtlas.Header}) and a second {@link ContactAtlas.Header}
 * (today, in practice, the literal same instance) named {@code
 * quantizerHeader}, and asserts {@code table.header().equals(quantizerHeader)}
 * as a defence-in-depth guard (see {@code ContactTable#header()}'s own
 * Javadoc for why this cross-check exists) - a future refactor that
 * decouples the two construction paths fails loudly here, not silently
 * downstream. <b>This guard compares {@link ContactAtlas.Header}
 * equality ONLY</b> - it does not construct or hold a {@link
 * PhaseQuantizer} instance (bead inviscid-0nx.23's fix-round removed the
 * unused {@code quantizer} field and its package-private accessor, which
 * had zero callers anywhere including tests; {@code quantizerHeader}'s
 * name is retained for its OTHER live uses below - {@code subBinSteps},
 * {@code geometryResolution}, {@code memberRadius}).
 *
 * @author halhildebrand
 */
public class LatticeGasAutomaton implements QuantaField, TickDriver {

    @FunctionalInterface
    public interface Processor {
        void process(int[] phase, long[] quanta);
    }

    /**
     * One tick's outcome: {@link #signedTransferTotal} is the sum, over
     * every resolved collision this tick, of {@code delta.deltaA() +
     * delta.deltaB()} - always zero by {@link CollisionRule.Delta}'s own
     * construction-time invariant.
     */
    public record TickResult(int tick, long signedTransferTotal)
        implements TickReport {
    }

    private static final int CUBES_PER_CELL   = 5;
    private static final int MEMBERS_PER_CUBE = 6;

    private final Point3i             extent;
    private final int[]               phase;
    private final long[]              quanta;
    private List<Integer>             lastTouched = List.of();
    private final ContactTable        table;
    private final FineStepContactTable fineContacts;
    private final CollisionTable      collisions;
    private final CollisionStatistics statistics;
    private final FccNeighborhood     neighborhood;

    private final int nLga;
    private final int subBinSteps;
    private final int phaseResolution;
    private final int geometryResolution;
    private final int fineStepDivisor;

    /**
     * @param extent     the periodic-wrap extent; delegated to {@link
     *                   FccNeighborhood}'s constructor for validation
     * @param atlas      the transcribed contact atlas; {@link
     *                   ContactTable} and {@link PhaseQuantizer} are both
     *                   built from {@code atlas.header()} - see class
     *                   Javadoc, "Header pairing"
     * @param collisions the frozen, conservation-exact collision table
     *                   (bead inviscid-0nx.20) - typically {@link
     *                   CollisionTable#buildFromPhaseARule}
     * @param statistics the collision-statistics recorder this driver
     *                   populates every tick, exactly like {@link
     *                   CollisionSweep} does (bead inviscid-ckn's R6:
     *                   verified substrate-agnostic, zero {@code
     *                   Necronomata}/{@link QuantaField} type coupling)
     */
    public LatticeGasAutomaton(Point3i extent, ContactAtlas atlas,
                                CollisionTable collisions,
                                CollisionStatistics statistics) {
        this(extent, ContactTable.of(atlas), atlas.header(), collisions,
             statistics);
    }

    /**
     * Package-private seam (mirrors {@code CollisionSweep.magnitudeToRecord}'s
     * established package-visible-for-testing pattern): exposes the
     * {@link ContactTable} / {@link ContactAtlas.Header} pairing directly
     * so a same-package test can inject a deliberate mismatch and prove
     * the guard fires, without which the check in the public constructor
     * (which always pairs from the literal same instance) would be
     * untestable, vacuous defensive code.
     */
    LatticeGasAutomaton(Point3i extent, ContactTable table,
                         ContactAtlas.Header quantizerHeader,
                         CollisionTable collisions,
                         CollisionStatistics statistics) {
        // Reused, not re-implemented: validates even-and-at-least-4-per-axis.
        this.neighborhood = new FccNeighborhood(extent);
        this.extent = new Point3i(extent);
        if (!table.header().equals(quantizerHeader)) {
            throw new IllegalStateException("ContactTable/PhaseQuantizer header pairing violated (bead inviscid-0nx.19's contract): "
                                             + "table.header()=" + table.header()
                                             + " but quantizerHeader="
                                             + quantizerHeader);
        }
        this.table = table;
        this.collisions = collisions;
        this.statistics = statistics;
        this.nLga = table.nLga();
        this.subBinSteps = quantizerHeader.subBinSteps();
        // Computed, never imported: Necronomata.PHASE_RESOLUTION's own
        // javadoc forbids reaching for it as an LGA parameter.
        this.phaseResolution = nLga * subBinSteps;

        // USER DECISION 2026-08-08, FINAL (test-8 finding, see class
        // Javadoc "Contact firing" section): the fine-step contact
        // structure, built ONCE from the SAME atlas header's geometry
        // (reused, not reimplemented -- ContactComboCache +
        // ContactPredicate, the exact machinery that built the atlas).
        this.geometryResolution = quantizerHeader.geometryResolution();
        if (phaseResolution % geometryResolution != 0) {
            throw new IllegalStateException("phaseResolution (" + phaseResolution
                                             + ") must be evenly divisible by geometryResolution ("
                                             + geometryResolution
                                             + ") for the fine-step accumulator/step conversion to be exact -- "
                                             + "same divisibility discipline as PhaseQuantizer's nLga/geometryResolution guard");
        }
        this.fineStepDivisor = phaseResolution / geometryResolution;
        ContactPredicate predicate = new ContactPredicate(new MemberGeometry(geometryResolution,
                                                                               quantizerHeader.memberRadius()));
        this.fineContacts = FineStepContactTable.buildFor(predicate,
                                                            geometryResolution,
                                                            quantizerHeader.memberRadius());

        int length = 30 * extent.x * extent.y * extent.z;
        this.phase = new int[length];
        this.quanta = new long[length];
    }

    /** @return M, the sub-bin accumulator steps per contact bin. */
    public int subBinSteps() {
        return subBinSteps;
    }

    /** @return the {@link CollisionStatistics} this driver records into. */
    public CollisionStatistics statistics() {
        return statistics;
    }

    /**
     * Raw-array escape hatch for deterministic seeding (mirrors {@code
     * Necronomata#process(Necronomata.Processor)}): writes {@code phase}
     * MUST stay within {@code [0, phaseResolution)} - unlike {@code
     * Necronomata}'s wrap-on-step forgiveness, this class does not
     * re-validate a seeded out-of-range phase, since seeding is a test/
     * initial-condition concern, not a tick-path one.
     */
    public void process(Processor action) {
        action.process(phase, quanta);
    }

    // --- QuantaField ---

    @Override
    public Point3i extent() {
        return new Point3i(extent);
    }

    @Override
    public int slotCount() {
        return quanta.length;
    }

    @Override
    public long quantaAt(int slot) {
        return quanta[slot];
    }

    /**
     * Constant {@code true}: this substrate's quanta are {@code long}-
     * backed, structurally incapable of the float32 representation
     * corruption {@code Necronomata}'s {@code isExactAt} checks for (T2
     * design-ckn-lattice-seam.md §2, risk R1). A caller reporting this
     * MUST say "structurally exact (integer storage)", never present
     * this as a passing corruption check - there is no corruption class
     * this substrate can exhibit here.
     */
    @Override
    public boolean isExactAt(int slot) {
        return true;
    }

    /**
     * The FINE accumulator phase, {@code 2*pi*phase[slot]/phaseResolution}
     * - NEVER {@link PhaseQuantizer#centre(int)} (T2
     * analysis-73v-spectral-conversion-and-cadence.md §5: a bin-centre
     * read injects a deterministic 24-level staircase quantisation
     * artifact into any angle-spectrum instrument that samples this
     * accessor - see pin test {@code phaseAtIsTheFineAccumulatorNotTheContactBinCentre}).
     * A MEASUREMENT read outside the tick path - {@link #phase} itself
     * stays {@code int} throughout {@link #tick(int)}.
     */
    @Override
    public float phaseAt(int slot) {
        return (float) (2.0 * Math.PI * phase[slot] / phaseResolution);
    }

    /**
     * {@code N_lga * subBinSteps} - 3600 under the user's 2A cadence
     * decision, matching {@code Necronomata.PHASE_RESOLUTION} NUMERICALLY
     * without being SOURCED from it (pin test {@code
     * bothSubstratesReportPhaseResolution3600UnderCadence2A}).
     */
    @Override
    public int phaseResolution() {
        return phaseResolution;
    }

    @Override
    public void forEachCell(Consumer<? super Point3i> action) {
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

    @Override
    public int indexOfCell(Point3i cell) {
        return 30 * ((cell.x * extent.y + cell.y) * extent.z + cell.z);
    }

    // --- TickDriver ---

    @Override
    public QuantaField field() {
        return this;
    }

    /**
     * Table-driven, synchronous even-parity update - see class Javadoc.
     */
    @Override
    public TickResult tick(int tickNumber) {
        int length = quanta.length;
        long[] deltaQuanta = new long[length];
        long signedTotal = 0L;

        // Bead inviscid-10d's exactness-ceiling guard applies here too;
        // touched slots are tracked so the check runs ONCE per member on
        // the tick's FINAL post-delta total, mirroring CollisionSweep's
        // own "Post-throw failure contract" (never an intra-tick partial
        // sum).
        List<Integer> touched = new ArrayList<>();

        int[] stepsA = new int[30];
        int[] stepsB = new int[30];
        for (Point3i cellA : cellVisitOrder()) {
            int baseA = indexOfCell(cellA);
            for (int direction = 1; direction <= 6; direction++) {
                Point3i cellB = neighborhood.neighbor(cellA, direction);
                int baseB = indexOfCell(cellB);
                for (int s = 0; s < 30; s++) {
                    stepsA[s] = phase[baseA + s] / fineStepDivisor;
                    stepsB[s] = phase[baseB + s] / fineStepDivisor;
                }
                for (int cubeA = 0; cubeA < CUBES_PER_CELL; cubeA++) {
                    for (int memberA = 0; memberA < MEMBERS_PER_CUBE; memberA++) {
                        int slotA = cubeA * MEMBERS_PER_CUBE + memberA;
                        int indexA = baseA + slotA;
                        int stepA = stepsA[slotA];
                        for (int cubeB = 0; cubeB < CUBES_PER_CELL; cubeB++) {
                            for (int memberB = 0; memberB < MEMBERS_PER_CUBE; memberB++) {
                                int slotB = cubeB * MEMBERS_PER_CUBE
                                            + memberB;
                                if (!fineContacts.contacts(direction, cubeA,
                                                            memberA, stepA,
                                                            cubeB, memberB,
                                                            stepsB[slotB])) {
                                    continue;
                                }
                                int indexB = baseB + slotB;
                                long quantaA = quanta[indexA];
                                long quantaB = quanta[indexB];
                                CollisionRule.Delta delta = collisions.lookup(quantaA,
                                                                               quantaB);

                                deltaQuanta[indexA] += delta.deltaA();
                                deltaQuanta[indexB] += delta.deltaB();
                                signedTotal += delta.deltaA()
                                               + delta.deltaB();

                                touched.add(indexA);
                                touched.add(indexB);

                                statistics.recordCollision(cellA, cubeA,
                                                            memberA, cellB,
                                                            cubeB, memberB,
                                                            direction,
                                                            Math.abs(delta.deltaA()),
                                                            tickNumber);
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < length; i++) {
            if (deltaQuanta[i] != 0L) {
                quanta[i] += deltaQuanta[i];
            }
            // Unconditional, every slot, every tick - matching
            // Necronomata.step()'s own convention of unconditionally
            // recomputing deltaA/advancing angle regardless of whether
            // this tick changed frequency. quanta[i]==0 is still a true
            // no-op (phase + 0 == phase), so this costs nothing extra
            // for a static lattice (test zeroQuantaLatticeIsStatic).
            phase[i] = (int) Math.floorMod((long) phase[i] + quanta[i],
                                            (long) phaseResolution);
        }

        for (int index : touched) {
            checkExactnessCeiling(index, quanta[index]);
        }

        this.lastTouched = List.copyOf(touched);

        return new TickResult(tickNumber,
                               signedTransferTotalToReport(signedTotal));
    }

    /**
     * Package-private test observability seam (critic finding on the
     * .21 automaton-arc final review): the flat slot indices touched by
     * a resolved collision during the MOST RECENT {@link #tick(int)}
     * call, in scan order, with a repeated index whenever a member was
     * touched by more than one collision this same tick. Lets a test
     * assert non-vacuity for the traversal-order-invariance proof (
     * {@code LatticeGasAutomatonTest.updateIsSynchronous}) -- without
     * this, a test window with zero same-tick multi-touches cannot
     * distinguish a synchronous update from a naive sequential one, and
     * silently regressing to such a window is undetectable by inspection
     * alone.
     */
    List<Integer> lastTouchedIndices() {
        return lastTouched;
    }

    /**
     * The cell-visitation order the scan phase of {@link #tick(int)}
     * uses. Package-private, overridable ONLY for the traversal-order-
     * invariance proof ({@code LatticeGasAutomatonTest.updateIsSynchronous}):
     * because the scan phase never mutates {@link #phase}/{@link
     * #quanta} (only reads them and writes to a separate {@code
     * deltaQuanta} accumulator), the tick's result is provably
     * independent of this order - overriding it and observing identical
     * output is the empirical half of that proof.
     */
    List<Point3i> cellVisitOrder() {
        List<Point3i> cells = new ArrayList<>();
        forEachCell(cells::add);
        return cells;
    }

    /**
     * The reconciliation-negative-control seam (mirrors {@code
     * CollisionSweep.magnitudeToRecord}'s established package-visible-
     * for-testing pattern): production behaviour is the identity
     * function - {@code computed} is already provably zero by {@link
     * CollisionRule.Delta}'s own construction-time invariant. Exists
     * ONLY so a same-package negative-control test can override it to
     * lie, proving {@code AuditedRun}'s reconciliation wiring catches an
     * LGA-specific mismatch, not just a generic fake driver's.
     */
    long signedTransferTotalToReport(long computed) {
        return computed;
    }

    /**
     * Bead inviscid-10d's guard, reusing {@link
     * CollisionSweep#QUANTA_EXACTNESS_SAFETY_MARGIN} rather than
     * hardcoding the ceiling a second time.
     */
    private static void checkExactnessCeiling(int index, long quanta) {
        if (Math.abs(quanta) >= CollisionSweep.QUANTA_EXACTNESS_SAFETY_MARGIN) {
            throw new IllegalStateException("Slot " + index
                                             + " quanta magnitude " + quanta
                                             + " has reached the float32-exactness safety margin ("
                                             + CollisionSweep.QUANTA_EXACTNESS_SAFETY_MARGIN
                                             + ") - refusing to let this member's quanta random-walk further.");
        }
    }
}
