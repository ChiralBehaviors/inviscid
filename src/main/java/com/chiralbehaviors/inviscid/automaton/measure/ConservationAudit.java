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

package com.chiralbehaviors.inviscid.automaton.measure;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.automaton.QuantaField;

/**
 * The acceptance instrument for the collision rules (bead inviscid-0nx.14 /
 * .20), built before those rules exist so it cannot be fitted to what it is
 * meant to judge.
 *
 * <p>Invariant: the total {@code frequency} quanta across the whole
 * lattice is conserved -- members are channels, frequency quanta are the
 * particles that hop between them at collisions (see the design memo,
 * inviscid/design-jitterbug-lga.md). A legitimate collision rule
 * redistributes quanta between members but never changes the lattice-wide
 * total; any tick where the total differs from the baseline captured at
 * construction is, by definition, a conservation violation.
 *
 * <p>Sums are computed with exact {@code long} arithmetic over integer
 * quanta counts -- never float summation, which would itself be an
 * epsilon trap across many entries. Every {@code frequency} slot is also
 * checked against {@link Math#rint(double)}: a slot that has drifted off
 * an integer value is a representation-corruption violation in its own
 * right, reported rather than silently rounded.
 *
 * <p>Dual-mode: {@link #auditTick(int)} always returns an
 * {@link AuditResult} usable as a test assertion
 * ({@code assertTrue(result.isClean())}); when {@link #isStrict()} is
 * enabled the same call additionally throws
 * {@link ConservationViolationException} the instant a violation is
 * detected, so the audit can be wired into the simulation loop itself as
 * an in-run assertion.
 *
 * <p><b>Construction-time corruption is not retained.</b> Any
 * representation corruption present in the lattice at construction time
 * is folded into the baseline (rounded via {@link Math#rint(double)}) and
 * discarded, not reported -- there is no prior tick to report it against.
 * If the corruption persists, it is re-flagged on the very first
 * {@link #auditTick(int)} call, since every slot is re-checked every
 * tick.
 *
 * <p><b>Scope: conservation-only, not motion-correctness.</b> This audit
 * verifies that the lattice-wide total is preserved and localizes which
 * slots changed since the last audited tick. It has a fundamental blind
 * spot: a symmetric transfer applied twice (e.g. a collision rule bug
 * that double-fires, moving +1/-1 between the same two members twice
 * instead of once) is zero-sum both times and therefore invisible to
 * conservation-only auditing. A.3/A.4 collision-rule authors must not
 * over-trust this audit for that error class -- it catches quanta
 * appearing or disappearing, not quanta moving to the wrong place or the
 * wrong number of times.
 *
 * <p><b>Frequency-array re-fetch contract (bead inviscid-36g).</b> This
 * audit does not hold a cached copy of the lattice's quanta across
 * calls. Every reader ({@link #auditTick(int)}, {@link
 * #currentTotalQuanta()}, and the internal snapshot machinery they
 * share) re-reads {@link QuantaField#quantaAt(int)} /
 * {@link QuantaField#isExactAt(int)} fresh, per slot, at the start of
 * the call -- it audits whatever {@code automaton} reports AT CALL
 * TIME, so it stays attached even if a future tick implementation
 * swaps the substrate's internal storage between buffers. Only the
 * constructor performs a one-time read, to establish
 * {@link #baselineTotal}.
 *
 * <p><b>Substrate-agnostic read seam (bead inviscid-ckn / inviscid-0nx.21).</b>
 * This class is constructed against a {@link QuantaField}, not a
 * concrete {@code Necronomata} -- any substrate that implements the
 * seam (today: {@code Necronomata}; the formal LGA once it lands) is
 * auditable unchanged. {@link QuantaField#isExactAt(int)} is constant
 * {@code true} for an integer-backed substrate; such a substrate's
 * {@link Violation.Kind#REPRESENTATION_CORRUPTION} check is
 * structurally impossible to trip, not merely one that happens to
 * pass -- a caller reporting the result should say so, not present the
 * tautology as evidence.
 *
 * @author halhildebrand
 */
public class ConservationAudit {

    /**
     * One violation: exact localization of a single {@code frequency}
     * slot whose value is unaccounted for -- naming the offending
     * {@code cell}, {@code cube}, {@code member}, and {@code tick}.
     *
     * <p><b>Direction attribution is out of scope here (bead
     * inviscid-17p).</b> This record originally carried a nullable
     * {@code direction} field, but nothing in this class ever populated
     * it -- {@code violationAt(...)} was called with {@code direction =
     * null} at every call site (both the representation-corruption path
     * in {@code snapshotExact} and the conservation-violation path in
     * {@code diff}), so the field was unconditionally null, not merely
     * "null when unknown". It has been removed rather than kept
     * structurally dead. Diffing the raw {@code frequency} array alone
     * never reveals which collision direction produced a change --
     * direction attribution requires cross-referencing
     * {@code CollisionStatistics}'s caller-supplied direction for the
     * same {@code (cell, cube, member, tick)}, a reconciliation this
     * class does not perform. NOTE the aggregate audit-stats cross-check
     * registered on bead inviscid-0nx.14 (per-tick totals reconciliation,
     * from inviscid-ce3) does NOT provide per-violation direction
     * attribution -- that narrower deliverable is tracked separately as
     * bead inviscid-1yk (blocks .14); the original bead .8 deliverable text
     * ("violation report naming cell, cube, member, direction, and
     * tick") is only partially met by this class alone until that
     * reconciliation lands.
     *
     * <p><b>Check {@link #kind()} before reading {@code cell()} /
     * {@code cube()} / {@code member()}.</b> A
     * {@link Kind#RESIDUAL_DRIFT} entry does not localize any single
     * slot -- it carries sentinel coordinates ({@code cell = (-1,-1,-1)},
     * {@code cube = -1}, {@code member = -1}) meaning "the lattice-wide
     * total is still off, but nothing in particular changed this tick".
     * Treating those sentinels as a real cell/cube/member is a bug in
     * the consumer, not in this class.
     */
    public record Violation(Point3i cell, int cube, int member, int tick,
                             Kind kind, long previousValue, long newValue) {

        public enum Kind {
            /**
             * The lattice-wide total quanta changed relative to the
             * baseline; this slot is one of the (possibly several)
             * indices whose value changed since the previous audited
             * tick and is therefore a candidate cause.
             */
            CONSERVATION_VIOLATION,
            /**
             * The raw float slot no longer holds an exact integer value
             * ({@code frequency[i] != Math.rint(frequency[i])}) --
             * representation corruption, reported unconditionally and
             * never absorbed by a tolerance.
             */
            REPRESENTATION_CORRUPTION,
            /**
             * The lattice-wide total still differs from the baseline
             * this tick, but no individual slot changed since the
             * previous audited tick -- a drift that was introduced on an
             * earlier tick (and localized there) has simply persisted
             * unchanged. {@code cell}/{@code cube}/{@code member} are
             * sentinel values ({@code (-1,-1,-1)} / {@code -1} /
             * {@code -1}); this entry exists so a caller relying solely
             * on {@code isClean()} never sees a false "clean" once a
             * drift has stabilized. See the earlier tick's ledger /
             * violation report for the original localization.
             */
            RESIDUAL_DRIFT
        }
    }

    /**
     * One tick's ledger entry: the lattice-wide total before and after
     * this audited tick, and the cumulative drift from the original
     * baseline.
     */
    public record LedgerEntry(int tick, long totalBefore, long totalAfter,
                               long cumulativeDrift) {
    }

    /**
     * The outcome of one {@link #auditTick(int)} call.
     */
    public record AuditResult(int tick, long totalQuanta, long baselineTotal,
                               List<Violation> violations) {

        public AuditResult {
            violations = List.copyOf(violations);
        }

        public boolean isClean() {
            return violations.isEmpty();
        }
    }

    /**
     * Thrown by {@link #auditTick(int)} when {@link #isStrict()} is
     * {@code true} and a violation is detected -- the in-run assertion
     * mode.
     */
    public static class ConservationViolationException extends RuntimeException {
        private static final long serialVersionUID = 1L;

        private final AuditResult result;

        public ConservationViolationException(AuditResult result) {
            super("Conservation violated at tick " + result.tick() + ": "
                  + result.violations());
            this.result = result;
        }

        public AuditResult result() {
            return result;
        }
    }

    /** 30 floats per cell: 5 cubes x 6 members (see Necronomata javadoc). */
    private static final int MEMBERS_PER_CUBE = 6;

    private final QuantaField automaton;
    private final Point3i     extent;

    private final long           baselineTotal;
    private long[]               previousSnapshot;
    private final List<LedgerEntry> ledger = new ArrayList<>();

    private boolean strict;

    public ConservationAudit(QuantaField automaton) {
        this(automaton, false);
    }

    public ConservationAudit(QuantaField automaton, boolean strict) {
        this.automaton = automaton;
        this.extent = automaton.extent();
        this.strict = strict;

        // Constructor-time read, ONLY to establish the baseline snapshot
        // (bead inviscid-36g). Every subsequent reader re-reads its own
        // fresh copy via QuantaField#quantaAt / #isExactAt -- see the
        // class Javadoc's "Frequency-array re-fetch contract" section.
        this.previousSnapshot = snapshotExact(new ArrayList<>(), 0);
        this.baselineTotal = sum(previousSnapshot);
    }

    public boolean isStrict() {
        return strict;
    }

    public void setStrict(boolean strict) {
        this.strict = strict;
    }

    /**
     * @return the current lattice-wide total quanta, as exact
     *         {@code long} arithmetic. Does not itself flag
     *         representation corruption -- corrupted slots are rounded
     *         via {@link Math#rint(double)} for the purposes of this sum;
     *         use {@link #auditTick(int)} to have corruption reported as
     *         a violation.
     */
    public long currentTotalQuanta() {
        long total = 0L;
        int n = automaton.slotCount();
        for (int i = 0; i < n; i++) {
            total += automaton.quantaAt(i);
        }
        return total;
    }

    public List<LedgerEntry> ledger() {
        return Collections.unmodifiableList(ledger);
    }

    /**
     * Re-synchronizes the per-slot localization reference
     * ({@code previousSnapshot}) to the current lattice state. This does
     * <b>not</b> reset {@link #baselineTotal} (the global conservation
     * invariant established at construction) or clear {@link #ledger()}
     * -- those are permanent for the life of this audit. Concretely: if
     * the lattice-wide total is already off from the baseline when
     * {@code checkpoint()} is called, it stays off, and the next
     * {@link #auditTick(int)} will still report a
     * {@link Violation.Kind#RESIDUAL_DRIFT} violation for it (with no
     * fresh per-slot delta, since nothing has changed since the
     * checkpoint). What {@code checkpoint()} changes is only which prior
     * state fresh per-slot diffs are computed against -- useful for a
     * caller that wants subsequent violation reports to localize only
     * NEW changes, not replay an already-reported drift as if it were
     * new. Not needed by the standard construct-then-audit workflow.
     */
    public void checkpoint() {
        this.previousSnapshot = snapshotExact(new ArrayList<>(), 0);
    }

    /**
     * Audits the lattice at {@code tick}. Always returns a result;
     * additionally throws {@link ConservationViolationException} when
     * {@link #isStrict()} and the result is not clean.
     */
    public AuditResult auditTick(int tick) {
        List<Violation> corruption = new ArrayList<>();
        long[] current = snapshotExact(corruption, tick);
        long currentTotal = sum(current);
        long totalBefore = sum(previousSnapshot);
        long tickOwnDelta = currentTotal - totalBefore;

        List<Violation> violations = new ArrayList<>(corruption);
        if (tickOwnDelta != 0) {
            // This tick's own per-slot changes do not net to zero --
            // a fresh, this-tick violation. Localize every changed slot
            // (today's behavior; not narrowed to "the" culprit slot,
            // since more than one may have moved).
            violations.addAll(diff(previousSnapshot, current, tick));
        } else if (currentTotal != baselineTotal) {
            // This tick is internally balanced (its own delta sums to
            // zero -- e.g. a legitimate two-member transfer, or no
            // activity at all) but the lattice-wide total still differs
            // from the original baseline because of a PRIOR tick's
            // unhealed drift. Do NOT mislabel this tick's balanced
            // per-slot changes as fresh CONSERVATION_VIOLATION entries
            // (that was the bug: a balanced transfer occurring after an
            // earlier leak was reported as two new violations even
            // though it moved nothing net). Report the persisting
            // divergence via the sentinel instead, so isClean() still
            // correctly reports false without falsely blaming this
            // tick's balanced activity.
            violations.add(residualDriftViolation(tick, currentTotal));
        }

        ledger.add(new LedgerEntry(tick, totalBefore, currentTotal,
                                    currentTotal - baselineTotal));
        previousSnapshot = current;

        AuditResult result = new AuditResult(tick, currentTotal,
                                              baselineTotal, violations);
        if (strict && !result.isClean()) {
            throw new ConservationViolationException(result);
        }
        return result;
    }

    /**
     * Exact per-slot snapshot as {@code long}s. Any slot that is not an
     * exact integer value (per {@link Math#rint(double)}) is both rounded
     * for the snapshot AND recorded into {@code corruptionSink} as a
     * {@link Violation.Kind#REPRESENTATION_CORRUPTION} violation -- never
     * silently absorbed.
     */
    private long[] snapshotExact(List<Violation> corruptionSink, int tick) {
        int n = automaton.slotCount();
        long[] snapshot = new long[n];
        for (int i = 0; i < n; i++) {
            long rounded = automaton.quantaAt(i);
            snapshot[i] = rounded;
            if (!automaton.isExactAt(i)) {
                corruptionSink.add(violationAt(i, tick, 0L, rounded,
                                                Violation.Kind.REPRESENTATION_CORRUPTION));
            }
        }
        return snapshot;
    }

    private List<Violation> diff(long[] before, long[] after, int tick) {
        List<Violation> changed = new ArrayList<>();
        for (int i = 0; i < after.length; i++) {
            if (before[i] != after[i]) {
                changed.add(violationAt(i, tick, before[i], after[i],
                                         Violation.Kind.CONSERVATION_VIOLATION));
            }
        }
        return changed;
    }

    private Violation violationAt(int index, int tick, long previousValue,
                                   long newValue, Violation.Kind kind) {
        int cellLinear = index / 30;
        int localIndex = index % 30;
        int k = cellLinear % extent.z;
        int rem = cellLinear / extent.z;
        int j = rem % extent.y;
        int i = rem / extent.y;
        int cube = localIndex / MEMBERS_PER_CUBE;
        int member = localIndex % MEMBERS_PER_CUBE;
        return new Violation(new Point3i(i, j, k), cube, member, tick, kind,
                              previousValue, newValue);
    }

    private Violation residualDriftViolation(int tick, long currentTotal) {
        return new Violation(new Point3i(-1, -1, -1), -1, -1, tick,
                              Violation.Kind.RESIDUAL_DRIFT, baselineTotal,
                              currentTotal);
    }

    private static long sum(long[] values) {
        long total = 0L;
        for (long v : values) {
            total += v;
        }
        return total;
    }
}
