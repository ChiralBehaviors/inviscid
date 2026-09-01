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

import java.util.List;

import com.chiralbehaviors.inviscid.automaton.lga.CollisionSweep;
import com.chiralbehaviors.inviscid.automaton.lga.HybridAutomaton;
import com.chiralbehaviors.inviscid.automaton.lga.TickDriver;
import com.chiralbehaviors.inviscid.automaton.lga.TickReport;

/**
 * The audited tick driver (bead inviscid-0nx.15's pre-close requirement,
 * registered on that bead's NOTES from the .14 review; genericized over
 * {@link TickDriver} by bead inviscid-ckn / inviscid-0nx.21 -- see
 * {@link TickOutcome}'s Javadoc): wires
 * {@link ConservationAudit#auditTick(int)} and {@link
 * CollisionSweep#reconcileWithLedger(TickReport, long)}
 * together, invoked once per tick, alongside {@link TickDriver#tick(int)}.
 * Per the locked design of record, this is THE way to run an audited
 * simulation -- any {@link TickDriver}, including a future formal LGA's,
 * MUST be driven through this class, not through a second reconciliation
 * path.
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
 * <li>{@link CollisionSweep#reconcileWithLedger(TickReport, long)}
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
     *
     * <p><b>Naming debt (bead inviscid-ckn / inviscid-0nx.21, T2
     * design-ckn-lattice-seam.md §3.3).</b> {@code collisionResult}'s
     * TYPE was widened from {@code CollisionSweep.TickResult} to
     * {@link TickReport} so any {@link TickDriver} (not just the Phase A
     * hybrid) can drive an {@code AuditedRun}; the component NAME was
     * deliberately kept -- a {@link TickReport} need not describe
     * collisions, so the name now reads slightly wider than its type.
     * The single production reader ({@code ContactAtlasGenerator}, which
     * is intrinsically Phase-A-specific) narrows back to
     * {@code CollisionSweep.TickResult} via {@code instanceof}.
     */
    public record TickOutcome(TickReport collisionResult,
                               ConservationAudit.AuditResult auditResult) {
    }

    private final TickDriver        automaton;
    private final ConservationAudit audit;

    /**
     * @param automaton the tick driver to advance each tick -- any
     *                  {@link TickDriver} (today: {@link HybridAutomaton};
     *                  the formal LGA once it lands)
     * @param audit     the conservation audit tracking {@code automaton}'s
     *                  lattice; MUST have been constructed against the
     *                  same {@link com.chiralbehaviors.inviscid.automaton.QuantaField}
     *                  {@code automaton} drives (this class does not
     *                  itself verify that -- a mismatched pair would
     *                  silently audit a different lattice than the one
     *                  being advanced)
     */
    public AuditedRun(TickDriver automaton, ConservationAudit audit) {
        this.automaton = automaton;
        this.audit = audit;
    }

    /**
     * @return the tick driver this run advances.
     */
    public TickDriver automaton() {
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
        TickReport collisionResult = automaton.tick(tickNumber);
        ConservationAudit.AuditResult auditResult = audit.auditTick(tickNumber);

        List<ConservationAudit.LedgerEntry> ledger = audit.ledger();
        ConservationAudit.LedgerEntry entry = ledger.get(ledger.size() - 1);
        CollisionSweep.reconcileWithLedger(collisionResult,
                                            entry.totalAfter()
                                            - entry.totalBefore());

        return new TickOutcome(collisionResult, auditResult);
    }
}
