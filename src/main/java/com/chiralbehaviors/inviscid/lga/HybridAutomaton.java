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

import com.chiralbehaviors.inviscid.Necronomata;

/**
 * The hybrid automaton tick (bead inviscid-0nx.15, Phase A.4): continuous
 * angle evolution driven by per-tick geometric contact detection and
 * discrete collision-rule resolution.
 *
 * <h2>Decision: a composed class, not a {@code Necronomata} method
 * (recorded here, as the bead requires)</h2>
 * {@code Necronomata.process(Point3i)} is the reserved visitation site the
 * original design memo pointed at, but it stays empty/unused for Phase A --
 * see that method's own Javadoc for the closure of this decision. A
 * {@code HybridAutomaton} class was chosen instead because:
 * <ul>
 * <li>{@code Necronomata} is the state substrate (flat float arrays plus
 * {@code step()}'s continuous-evolution contract) and stays that way --
 * folding contact-scan orchestration, a {@link CollisionRule}, and a
 * {@link CollisionStatistics}-driven recording loop into it would make it
 * simultaneously "the arrays" and "the physics driver", the exact
 * responsibility mixing the design memo's layering (state substrate vs.
 * measurement/driver layers) argues against.</li>
 * <li>{@link CollisionSweep} (bead inviscid-0nx.14) already owns the
 * per-tick scan-resolve-record loop and its double-buffered snapshot
 * contract; a {@code HybridAutomaton} class composes it with {@code
 * Necronomata.step()} rather than duplicating or re-deciding that
 * contract on a {@code Necronomata} instance method.</li>
 * <li>It matches the plan's package layout (bead inviscid-0nx.15 files:
 * {@code lga/HybridAutomaton.java}) and keeps {@code Necronomata} testable
 * and reusable independent of any particular collision-rule wiring -- a
 * future Phase C formal LGA (bead inviscid-0nx.20) swaps in a different
 * driver composing the same {@code Necronomata} without touching it.</li>
 * </ul>
 *
 * <h2>Tick order</h2>
 * {@link #tick(int)} does exactly two things, in order:
 * <ol>
 * <li>{@link CollisionSweep#tick(int)} -- scans {@link ContactScan#scan()}
 * against the CURRENT (pre-tick) lattice state, resolves every found
 * contact via the configured {@link CollisionRule}, and accumulates every
 * decision into {@code Necronomata}'s {@code deltaF} accumulator. This
 * class does not reimplement double-buffering: {@code CollisionSweep}
 * already resolves every contact this tick against the frozen pre-tick
 * {@code frequency} snapshot, never a same-tick partially-accumulated
 * value (see that class's Javadoc, "Snapshot resolution within a tick",
 * bead inviscid-72s) -- REUSED DIRECTLY here, not reimplemented, per the
 * bead's explicit instruction.</li>
 * <li>{@link Necronomata#step()} -- applies the accumulated {@code
 * deltaF} to {@code frequency}, zeroes {@code deltaF} for the next tick,
 * recomputes {@code deltaA} from the now-current {@code frequency}, and
 * advances {@code angle}. A quantum absorbed THIS tick therefore moves its
 * receiving member's angle on this SAME tick (see {@code Necronomata}'s
 * own Javadoc on {@code deltaA} and {@code step()}).</li>
 * </ol>
 *
 * <h2>Layering: no dependency on {@code measure.ConservationAudit}</h2>
 * This class depends only on {@code Necronomata} and {@code CollisionSweep}
 * (the latter's existing, pre-.15 dependency on {@code
 * measure.CollisionStatistics} is unchanged, not expanded). It does NOT
 * import or reference {@code ConservationAudit} -- doing so would create a
 * {@code lga}&harr;{@code measure} import cycle, since {@code measure}
 * already imports {@code lga} (e.g. {@code CollisionStatistics} imports
 * {@code FccNeighborhood}; {@code ConservationAudit} lives alongside it in
 * {@code measure}). A caller that wants {@code ConservationAudit.auditTick}
 * and {@link CollisionSweep#reconcileWithLedger(CollisionSweep.TickResult,
 * long)} invoked together, per tick, alongside this class's {@link
 * #tick(int)}, should use {@code com.chiralbehaviors.inviscid.measure.
 * AuditedRun} -- the small measure-side driver that composes this class
 * with {@code ConservationAudit} without pulling that dependency down into
 * {@code lga}.
 *
 * @author halhildebrand
 */
public class HybridAutomaton {

    private final Necronomata   automaton;
    private final CollisionSweep sweep;

    /**
     * @param automaton the lattice state this tick advances; MUST be the
     *                  same instance {@code sweep} was constructed against
     *                  (this class does not itself verify that -- {@code
     *                  CollisionSweep} already owns its {@code automaton}
     *                  reference internally and applies every resolved
     *                  delta to it via {@code deltaF}, so a mismatched
     *                  pair would silently step a different lattice than
     *                  the one collisions were resolved against)
     * @param sweep     the per-tick collision-resolution loop (bead
     *                  inviscid-0nx.14), reused directly for its
     *                  double-buffered snapshot semantics
     */
    public HybridAutomaton(Necronomata automaton, CollisionSweep sweep) {
        this.automaton = automaton;
        this.sweep = sweep;
    }

    /**
     * @return the lattice state this instance advances.
     */
    public Necronomata automaton() {
        return automaton;
    }

    /**
     * @return the per-tick collision-resolution loop this instance drives.
     */
    public CollisionSweep sweep() {
        return sweep;
    }

    /**
     * Advances the lattice exactly one tick: resolve this tick's contacts
     * against the pre-tick state (via {@link CollisionSweep#tick(int)}),
     * then apply the accumulated transfer and advance angles (via {@link
     * Necronomata#step()}). See class Javadoc, "Tick order".
     *
     * @param tickNumber the tick number, passed through to {@link
     *                   CollisionSweep#tick(int)} for its own recording
     *                   and reconciliation bookkeeping
     * @return the resolved-contact result for this tick, straight from
     *         {@link CollisionSweep#tick(int)} -- callers that also want
     *         ledger reconciliation against a {@code ConservationAudit}
     *         should feed this into {@link
     *         CollisionSweep#reconcileWithLedger(CollisionSweep.TickResult,
     *         long)} themselves (see {@code measure.AuditedRun})
     * @throws CollisionSweep.ReconciliationException if {@link
     *                                                 CollisionSweep#tick(int)}'s
     *                                                 own recording-integrity
     *                                                 cross-check fails
     */
    public CollisionSweep.TickResult tick(int tickNumber) {
        CollisionSweep.TickResult result = sweep.tick(tickNumber);
        automaton.step();
        return result;
    }
}
