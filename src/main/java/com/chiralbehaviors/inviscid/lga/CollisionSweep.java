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
import java.util.Collections;
import java.util.List;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.measure.CollisionStatistics;

/**
 * The per-tick collision-rule loop (bead inviscid-0nx.14 / inviscid-1yk):
 * scans {@link ContactScan#scan()} for a tick's contacts, resolves each
 * with a {@link CollisionRule}, applies the resulting {@link
 * CollisionRule.Delta} through {@code Necronomata}'s sanctioned {@code
 * deltaF} accumulator (never writing {@code frequency}/{@code angle}
 * directly - see {@code Necronomata.process(Necronomata.Processor)}'s
 * Javadoc), and records each resolved contact with {@link
 * CollisionStatistics#recordCollision}.
 *
 * <p>Driving order (owned by the caller, not this class): {@code
 * ContactScan} implicitly reads the CURRENT lattice state each time
 * {@link #tick(int)} calls it, so a caller must run {@code tick(t)}, then
 * {@code Necronomata.step()}, then (if auditing) {@code
 * ConservationAudit.auditTick(t)}, in that order, once per tick - {@code
 * tick} does not call {@code step()} itself.
 *
 * <h2>Snapshot resolution within a tick (bead inviscid-72s)</h2>
 * Every contact resolved by one {@link #tick(int)} call reads quanta from
 * the FROZEN pre-tick {@code frequency} array - never {@code frequency[i]
 * + deltaF[i]} - so the outcome of resolving a contact never depends on
 * which other contacts touching the same member were resolved earlier in
 * this same tick. Deltas still accumulate additively into {@code deltaF}
 * ({@code deltaF[i] += delta.deltaX()}, so multiple contacts touching the
 * same member in one tick still all contribute), and every individual
 * {@link CollisionRule.Delta} is still exactly zero-sum by construction -
 * but the RESOLUTION DECISION for a given contact is a pure function of
 * the state as of tick start, full stop. Scan order remains the fixed,
 * deterministic order in which contacts are RECORDED (see {@link
 * TickResult#applied()}) - it just no longer affects which decision gets
 * made.
 *
 * <p>This matches two pieces of the locked design of record directly:
 * bead inviscid-0nx.15 (A.4, the double-buffered hybrid tick) already
 * mandates it in its own text - "all contacts are detected against the
 * pre-tick state... a single-buffer implementation makes the result
 * depend on scan order" - and Phase C's formal LGA (bead inviscid-0nx.20)
 * is a synchronous update on the even-parity sublattice, which is
 * exactly snapshot semantics generalized to the whole lattice, not just
 * one member. A rule loop that let earlier-in-tick contacts influence
 * later ones within the SAME tick would be a hybrid-only deviation from
 * both.
 *
 * <p><b>Rejected alternative: sequential resolution.</b> An earlier
 * version of this class resolved each contact against {@code
 * frequency[i] + deltaF[i]} (the running, same-tick-accumulated value),
 * reasoning that it was cheap and that "overdraw" (driving a member's
 * quanta negative) needed guarding against. Both were wrong: (1) scan
 * order then determined outcomes whenever a member appeared in more than
 * one contact in the same tick - not a corner case, ~18.8% of ticks
 * exhibit a same-member multi-contact at typical seed densities - and a
 * fixed scan order keyed to {@code ContactScan}'s canonical direction
 * list would then risk masquerading as genuine directional anisotropy
 * once bead B.5's directional statistics are compared, an artifact of
 * iteration order rather than physics; (2) quanta are signed {@code
 * long}s with no floor - "overdrawing" a member below zero is legal by
 * design, so there was never anything to protect against. See {@code
 * CollisionSweepTest.sameTickSameMemberMultiContactResolvesAgainstThePreTickSnapshot}
 * for the regression that pins this down (bead inviscid-72s).
 *
 * <h2>Recording convention: every resolved contact, including no-ops</h2>
 * Every contact {@link ContactScan#scan()} finds is "resolved" - even
 * when the rule decides {@link CollisionRule.Delta#noop()} - and every
 * one is recorded via {@link CollisionStatistics#recordCollision} (with
 * {@code transferMagnitude == 0} for no-ops) and retained in the returned
 * {@link TickResult#applied()} list. This is a deliberate choice, not the
 * "skip no-ops" alternative the bead notes offered: bead inviscid-1yk's
 * per-violation direction attribution needs to say which collision
 * DIRECTION(S) touched a member that tick, and a member that only
 * participated in no-op contacts is still "touched" in the sense a
 * debugging caller cares about - dropping no-ops would silently narrow
 * what attribution can ever report. See {@code
 * CollisionStatistics}'s class Javadoc ("What counts as a collision") for
 * how that class's own accessors distinguish the two.
 *
 * <h2>Reconciliation (bead inviscid-ce3)</h2>
 * Two independent cross-checks, both cheap and both assertable:
 * <ol>
 * <li><b>Recording integrity.</b> {@link
 * #magnitudeToRecord(CollisionRule.Delta)} is the single seam through
 * which a resolved delta's magnitude reaches {@link
 * CollisionStatistics#recordCollision}; {@link #tick(int)} also
 * independently sums {@code Math.abs(delta.deltaA())} straight from the
 * {@link CollisionRule.Delta} objects it applied, and throws {@link
 * ReconciliationException} if that independently-computed total ever
 * disagrees with what was actually recorded. A subclass that overrides
 * {@link #magnitudeToRecord(CollisionRule.Delta)} to lie is exactly the
 * "deliberately mis-recorded transfer" case this catches (see {@code
 * CollisionSweepTest.deliberatelyMisRecordedTransferIsCaught}).</li>
 * <li><b>Ledger reconciliation.</b> {@link
 * TickResult#signedTransferTotal()} is the sum, over every resolved
 * contact, of {@code delta.deltaA() + delta.deltaB()} - always exactly
 * zero, since {@link CollisionRule.Delta} itself enforces zero-sum at
 * construction. {@link #reconcileWithLedger(TickResult, long)} compares
 * this value against the tick's ACTUAL observed lattice-wide delta and
 * throws if they disagree - see that method's own Javadoc for the
 * caller contract and the specific failure mode this guards against.</li>
 * </ol>
 * This class intentionally does not depend on {@code ConservationAudit} -
 * the caller drives {@code scan -> tick() -> Necronomata.step() ->
 * ConservationAudit.auditTick()} in that order and passes the resulting
 * ledger delta back into {@link #reconcileWithLedger(TickResult, long)};
 * keeping the dependency one-directional (this class takes a primitive
 * {@code long}, not a {@code ConservationAudit} reference) avoids growing
 * the {@code lga}<->{@code measure} package coupling beyond the {@code
 * CollisionStatistics} dependency this loop already needs.
 *
 * <h2>Direction attribution (bead inviscid-1yk)</h2>
 * {@link TickResult#applied()} retains every resolved contact for the
 * tick; {@link #directionsTouching(TickResult, Point3i, int, int)}
 * cross-references it against a {@code (cell, cube, member)} triple (the
 * same addressing {@code ConservationAudit.Violation} reports) to answer
 * "which collision direction(s) touched this member this tick" - the
 * seam a caller uses to annotate a reported violation, without this class
 * needing any dependency on {@code ConservationAudit.Violation} itself.
 *
 * @author halhildebrand
 */
public class CollisionSweep {

    private static final int MEMBERS_PER_CUBE = 6;

    /**
     * One resolved contact for a tick: the {@link Contact} itself and the
     * {@link CollisionRule.Delta} the rule decided for it (possibly
     * {@link CollisionRule.Delta#noop()}).
     */
    public record AppliedCollision(Contact contact, CollisionRule.Delta delta) {
    }

    /**
     * The outcome of one {@link #tick(int)} call.
     *
     * @param tick                    the tick number
     * @param applied                 every resolved contact this tick, in
     *                                scan (== recording) order
     * @param appliedMagnitudeTotal   the independently-computed sum of
     *                                {@code Math.abs(delta.deltaA())}
     *                                across {@code applied}, computed
     *                                directly from the {@link
     *                                CollisionRule.Delta} objects
     * @param recordedMagnitudeTotal  the sum of the magnitudes actually
     *                                passed to {@link
     *                                CollisionStatistics#recordCollision}
     * @param signedTransferTotal     the sum, over {@code applied}, of
     *                                {@code delta.deltaA() +
     *                                delta.deltaB()} - always zero by
     *                                {@link CollisionRule.Delta}'s own
     *                                construction-time invariant
     */
    public record TickResult(int tick, List<AppliedCollision> applied,
                              long appliedMagnitudeTotal,
                              long recordedMagnitudeTotal,
                              long signedTransferTotal) {

        public TickResult {
            applied = List.copyOf(applied);
        }
    }

    /**
     * Thrown by {@link #tick(int)} when the independently-computed
     * applied-magnitude total disagrees with what was actually recorded
     * via {@link CollisionStatistics#recordCollision}, or by {@link
     * #reconcileWithLedger(TickResult, long)} when a tick's provably-zero
     * signed transfer total disagrees with the ledger's observed delta.
     */
    public static class ReconciliationException extends RuntimeException {
        private static final long serialVersionUID = 1L;

        public ReconciliationException(String message) {
            super(message);
        }
    }

    private final Necronomata         automaton;
    private final ContactScan         scan;
    private final CollisionRule       rule;
    private final CollisionStatistics statistics;

    public CollisionSweep(Necronomata automaton, ContactScan scan,
                           CollisionRule rule, CollisionStatistics statistics) {
        this.automaton = automaton;
        this.scan = scan;
        this.rule = rule;
        this.statistics = statistics;
    }

    /**
     * The recording seam: the magnitude actually passed to {@link
     * CollisionStatistics#recordCollision} for a resolved {@code delta}.
     *
     * <p><b>Visibility contract.</b> Deliberately package-private, not
     * {@code protected} or {@code public}: this is NOT a genuine
     * extension point for real subclasses (there is exactly one
     * production behavior - record the delta's true magnitude - and it
     * is not meant to vary). It exists ONLY so the same-package
     * reconciliation negative-control test ({@code
     * CollisionSweepTest.deliberatelyMisRecordedTransferIsCaught}) can
     * override it to lie about a magnitude and prove {@link #tick(int)}'s
     * recording-integrity cross-check catches the mismatch. Widening this
     * to {@code protected} would let any subclass anywhere silently
     * change what gets recorded without tripping that check's intended
     * purpose.
     */
    long magnitudeToRecord(CollisionRule.Delta delta) {
        return Math.abs(delta.deltaA());
    }

    /**
     * Resolves every contact {@link ContactScan#scan()} currently finds,
     * each against the frozen pre-tick {@code frequency} snapshot (see
     * class Javadoc, "Snapshot resolution within a tick"), applying every
     * decision via {@code Necronomata}'s {@code deltaF} accumulator and
     * recording it in scan order.
     *
     * @throws ReconciliationException if the recording-integrity
     *                                 cross-check fails (see class
     *                                 Javadoc)
     */
    public TickResult tick(int tickNumber) {
        List<Contact> contacts = scan.scan();
        List<AppliedCollision> applied = new ArrayList<>(contacts.size());
        long[] appliedTotal = { 0L };
        long[] recordedTotal = { 0L };
        long[] signedTotal = { 0L };

        automaton.process((angle, frequency, deltaA, deltaF) -> {
            for (Contact contact : contacts) {
                int indexA = automaton.indexOfCell(contact.cellA())
                             + contact.cubeA() * MEMBERS_PER_CUBE
                             + contact.memberA();
                int indexB = automaton.indexOfCell(contact.cellB())
                             + contact.cubeB() * MEMBERS_PER_CUBE
                             + contact.memberB();

                // Frozen pre-tick snapshot: read ONLY frequency[], never
                // frequency[i] + deltaF[i]. frequency[] is not written by
                // this loop (only deltaF is), so every contact this tick
                // sees the identical value here regardless of resolution
                // order - see class Javadoc, "Snapshot resolution within
                // a tick" (bead inviscid-72s).
                float currentAngleA = angle[indexA];
                float currentAngleB = angle[indexB];
                long snapshotQuantaA = Math.round((double) frequency[indexA]);
                long snapshotQuantaB = Math.round((double) frequency[indexB]);

                CollisionRule.Delta delta = rule.resolve(contact,
                                                          currentAngleA,
                                                          snapshotQuantaA,
                                                          currentAngleB,
                                                          snapshotQuantaB);

                deltaF[indexA] += delta.deltaA();
                deltaF[indexB] += delta.deltaB();

                long appliedMagnitude = Math.abs(delta.deltaA());
                long recordedMagnitude = magnitudeToRecord(delta);

                statistics.recordCollision(contact.cellA(), contact.cubeA(),
                                            contact.memberA(),
                                            contact.cellB(), contact.cubeB(),
                                            contact.memberB(),
                                            contact.direction(),
                                            recordedMagnitude, tickNumber);

                applied.add(new AppliedCollision(contact, delta));
                appliedTotal[0] += appliedMagnitude;
                recordedTotal[0] += recordedMagnitude;
                signedTotal[0] += delta.deltaA() + delta.deltaB();
            }
        });

        if (appliedTotal[0] != recordedTotal[0]) {
            throw new ReconciliationException("Recorded transfer magnitude total ("
                                               + recordedTotal[0]
                                               + ") disagrees with the applied transfer magnitude total ("
                                               + appliedTotal[0]
                                               + ") for tick " + tickNumber);
        }

        return new TickResult(tickNumber, applied, appliedTotal[0],
                               recordedTotal[0], signedTotal[0]);
    }

    /**
     * Cross-checks a tick's provably-zero {@link
     * TickResult#signedTransferTotal()} against the tick's ACTUAL
     * observed lattice-wide delta (typically {@code
     * ledgerEntry.totalAfter() - ledgerEntry.totalBefore()} from {@code
     * ConservationAudit}'s ledger, computed independently from the real
     * post-{@code step()} {@code frequency} array). See class Javadoc,
     * "Ledger reconciliation".
     *
     * <p><b>Failure mode this guards against.</b> {@code frequency} is a
     * {@code float32} slot; every legitimate value is an exact integer
     * (see {@code ConservationAudit}'s own {@code
     * Math.rint}-based corruption check), but a bug that let a
     * non-integral or out-of-precision value reach {@code frequency}
     * (e.g. a stray write bypassing {@code deltaF}, or quanta growing
     * past {@code 2^24} and silently losing exactness - see {@code
     * QuantaExchangeRuleTest.quantaStayWithinRepresentableRange}'s
     * boundary) would make the ledger's real {@code totalAfter -
   * totalBefore} disagree with this class's provably-zero {@code
     * signedTransferTotal} even though every individual {@code Delta}
     * this tick was exactly zero-sum. This method is the seam that turns
     * that silent float32-precision drop into a thrown exception instead
     * of a quietly-wrong lattice.
     *
     * <p><b>Caller contract.</b> Must be invoked once per tick, alongside
     * (immediately after) {@code ConservationAudit.auditTick(tick)}, with
     * that same tick's ledger delta - not batched, not skipped, not
     * called against a stale or mismatched tick's ledger entry.
     *
     * @throws ReconciliationException if the two disagree
     */
    public static void reconcileWithLedger(TickResult tickResult,
                                            long ledgerDelta) {
        if (tickResult.signedTransferTotal() != ledgerDelta) {
            throw new ReconciliationException("Tick " + tickResult.tick()
                                               + "'s signed transfer total ("
                                               + tickResult.signedTransferTotal()
                                               + ") disagrees with the ledger's observed delta ("
                                               + ledgerDelta + ")");
        }
    }

    /**
     * Bead inviscid-1yk: which collision direction(s) touched {@code
     * (cell, cube, member)} in {@code tickResult}'s tick.
     *
     * @return the directions (see {@link FccNeighborhood#DIRECTIONS}) of
     *         every contact in {@code tickResult} that had {@code (cell,
     *         cube, member)} as either endpoint, in {@code
     *         tickResult.applied()} order; empty if none did
     */
    public static List<Integer> directionsTouching(TickResult tickResult,
                                                     Point3i cell, int cube,
                                                     int member) {
        List<Integer> directions = new ArrayList<>();
        for (AppliedCollision applied : tickResult.applied()) {
            Contact contact = applied.contact();
            boolean touchesA = contact.cellA().equals(cell)
                                && contact.cubeA() == cube
                                && contact.memberA() == member;
            boolean touchesB = contact.cellB().equals(cell)
                                && contact.cubeB() == cube
                                && contact.memberB() == member;
            if (touchesA || touchesB) {
                directions.add(contact.direction());
            }
        }
        return Collections.unmodifiableList(directions);
    }
}
