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

package com.chiralbehaviors.inviscid.measure;

import java.util.List;

import com.chiralbehaviors.inviscid.lga.CollisionSweep;
import com.chiralbehaviors.inviscid.lga.HybridAutomaton;

/**
 * The audited hybrid-automaton driver (bead inviscid-0nx.15's pre-close
 * requirement, registered on that bead's NOTES from the .14 review): wires
 * {@link ConservationAudit#auditTick(int)} and {@link
 * CollisionSweep#reconcileWithLedger(CollisionSweep.TickResult, long)}
 * together, invoked once per tick, alongside {@link HybridAutomaton#tick(int)}.
 *
 * <h2>Why this class exists (the audit seam, chosen over the alternatives)</h2>
 * {@code CollisionSweep} deliberately does not depend on {@code
 * ConservationAudit} (see that class's own Javadoc, "Reconciliation") and
 * {@link HybridAutomaton} deliberately does not either (see that class's
 * own Javadoc, "Layering") -- both to avoid a {@code lga}&harr;{@code
 * measure} import cycle, since {@code measure} already imports {@code lga}
 * (this class is itself an instance of that existing, one-directional
 * {@code measure -> lga} dependency, not a new cycle). This class is the
 * "small measure-side driver" alternative the bead's NOTES offered
 * (composing {@code HybridAutomaton} + {@code ConservationAudit} +
 * {@code reconcileWithLedger}), chosen over threading a per-tick
 * hook/callback through {@code HybridAutomaton} itself: a driver class
 * keeps the wiring decision -- and the dependency on {@code
 * ConservationAudit} it requires -- entirely on the {@code measure} side,
 * with zero surface added to {@code HybridAutomaton} or {@code
 * CollisionSweep} for a concern (auditing) that only some callers need.
 *
 * <h2>What "wired per tick" means here</h2>
 * {@link #tick(int)} performs, in order, every tick:
 * <ol>
 * <li>{@link HybridAutomaton#tick(int)} -- resolve this tick's contacts
 * against the pre-tick snapshot and advance the lattice ({@code
 * Necronomata.step()}).</li>
 * <li>{@link ConservationAudit#auditTick(int)} -- audit the resulting
 * lattice-wide total against the running baseline, appending this tick's
 * {@link ConservationAudit.LedgerEntry} to the audit's ledger.</li>
 * <li>{@link CollisionSweep#reconcileWithLedger(CollisionSweep.TickResult, long)}
 * -- cross-check the tick's provably-zero {@code signedTransferTotal}
 * against that SAME ledger entry's {@code totalAfter - totalBefore}, per
 * that method's own caller contract ("once per tick... with that same
 * tick's ledger delta").</li>
 * </ol>
 * A caller using this class therefore cannot forget the reconciliation
 * call, skip it, or batch it -- it happens inside {@link #tick(int)}
 * itself, every time.
 *
 * @author halhildebrand
 */
public class AuditedRun {

    /**
     * One tick's combined outcome: the resolved-contact result and the
     * conservation-audit result for the same tick.
     */
    public record TickOutcome(CollisionSweep.TickResult collisionResult,
                               ConservationAudit.AuditResult auditResult) {
    }

    private final HybridAutomaton   automaton;
    private final ConservationAudit audit;

    /**
     * @param automaton the hybrid automaton to advance each tick
     * @param audit     the conservation audit tracking {@code automaton}'s
     *                  lattice; MUST have been constructed against the
     *                  same {@code Necronomata} instance {@code automaton}
     *                  drives (this class does not itself verify that --
     *                  a mismatched pair would silently audit a different
     *                  lattice than the one being advanced)
     */
    public AuditedRun(HybridAutomaton automaton, ConservationAudit audit) {
        this.automaton = automaton;
        this.audit = audit;
    }

    /**
     * @return the hybrid automaton this run advances.
     */
    public HybridAutomaton automaton() {
        return automaton;
    }

    /**
     * @return the conservation audit this run checks each tick.
     */
    public ConservationAudit audit() {
        return audit;
    }

    /**
     * Advances {@link #automaton()} one tick and audits/reconciles the
     * result -- see class Javadoc, "What 'wired per tick' means here".
     *
     * @param tickNumber the tick number
     * @return this tick's combined collision and audit outcome
     * @throws CollisionSweep.ReconciliationException if either {@link
     *                                                 HybridAutomaton#tick(int)}'s
     *                                                 own recording-integrity
     *                                                 cross-check fails, or
     *                                                 the ledger
     *                                                 reconciliation
     *                                                 disagrees
     * @throws ConservationAudit.ConservationViolationException if {@link
     *                                                           #audit()}
     *                                                           is in
     *                                                           strict mode
     *                                                           and a
     *                                                           conservation
     *                                                           violation
     *                                                           is detected
     */
    public TickOutcome tick(int tickNumber) {
        CollisionSweep.TickResult collisionResult = automaton.tick(tickNumber);
        ConservationAudit.AuditResult auditResult = audit.auditTick(tickNumber);

        List<ConservationAudit.LedgerEntry> ledger = audit.ledger();
        ConservationAudit.LedgerEntry entry = ledger.get(ledger.size() - 1);
        CollisionSweep.reconcileWithLedger(collisionResult,
                                            entry.totalAfter()
                                            - entry.totalBefore());

        return new TickOutcome(collisionResult, auditResult);
    }
}
