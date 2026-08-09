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

/**
 * v1 collision rule (bead inviscid-0nx.14, USER-decided 2026-08-08; see T2
 * {@code inviscid/decision-14-collision-rule-form.md}): TRANSFER ONE
 * QUANTUM from the higher-quanta member to the lower. {@code quantaA >
 * quantaB} yields {@code (-1, +1)}; {@code quantaA < quantaB} yields the
 * mirrored {@code (+1, -1)}; equal quanta is a defined no-op {@code (0,
 * 0)}.
 *
 * <h2>Candidate forms considered (locked decision)</h2>
 * <ul>
 * <li><b>(a) swap quanta between the two members</b> - REJECTED. Swap
 * freezes the global multiset of quanta values for all time (every
 * collision only permutes which member holds which value, never changes
 * the value set itself) - a spurious-conservation family, the classic LGA
 * "frozen invariant" disease, that would cripple any hydrodynamic
 * relaxation the coarse-grained limit is meant to exhibit.</li>
 * <li><b>(b) transfer one quantum from higher to lower</b> - CHOSEN (this
 * class). Exact signed-total conservation, participant-swap symmetric,
 * deterministic tie-break (equal quanta is defined, not arbitrary),
 * breaks the swap-rule's spurious multiset invariant, a gentle +/-1 kick
 * matched to the sparse contact density (~1.1e-4, measured in {@code
 * ContactPredicateTest}), an H-theorem-like relaxation toward a diffusive
 * coarse limit of the frequency field, and collision broadening that is
 * directly measurable as one-{@link
 * com.chiralbehaviors.inviscid.Necronomata#QUANTUM_RATE} line shifts
 * against the K=0 golden baseline (bead inviscid-0nx.7).</li>
 * <li><b>(c) negate-and-share</b> - REJECTED. Breaks exact signed-total
 * conservation in every natural reading (negating a member's quanta flips
 * the sign of its own contribution to the lattice-wide signed total, so
 * that total is not preserved); only "viable" by redefining the conserved
 * quantity itself, which would re-litigate the locked design (T2 {@code
 * inviscid/design-jitterbug-lga.md}).</li>
 * </ul>
 *
 * <h2>Reservoir model (locked decision)</h2>
 * SINGLE RESERVOIR: this rule reads and writes one quanta count per
 * member (the existing 30-floats/cell layout, {@code
 * Necronomata.frequency}), not two per-end counts. This literally
 * implements the design memo's "two-ended instant nonlocality" by
 * construction - a quantum absorbed anywhere on a member is immediately
 * available for a collision at either end of that same member, because
 * there is only one count to read. A two-ended model was rejected: it
 * would double the per-cell layout (60 floats), break every piece of
 * pinned index arithmetic that assumes 30 floats/cell ({@code
 * Necronomata#indexOfCell}, {@code ConservationAudit}, {@code
 * CollisionStatistics}), and introduce a delayed-availability semantics
 * that contradicts the memo's explicit "instant" framing.
 *
 * <h2>No long-overflow exposure</h2>
 * This rule never computes {@code quantaA - quantaB} or any other
 * arithmetic combination of the two inputs - only a three-way comparison
 * ({@code >}, {@code <}, {@code ==}) followed by one of three fixed
 * {@link Delta} constants. It is therefore immune to {@code long}
 * overflow at any input value, including {@link Long#MAX_VALUE} /
 * {@link Long#MIN_VALUE}, by construction. Applying the resulting
 * {@link Delta} is equally overflow-immune: whichever member holds the
 * extremal value is, by the comparison that selected the {@link Delta},
 * necessarily the one that LOSES a quantum - moving it away from its
 * boundary, never past it.
 *
 * @author halhildebrand
 */
public final class QuantaExchangeRule implements CollisionRule {

    private static final Delta TRANSFER_TO_B = new Delta(-1L, 1L);
    private static final Delta TRANSFER_TO_A = new Delta(1L, -1L);
    private static final Delta NO_OP         = Delta.noop();

    @Override
    public Delta resolve(Contact contact, float angleA, long quantaA,
                          float angleB, long quantaB) {
        if (quantaA > quantaB) {
            return TRANSFER_TO_B;
        }
        if (quantaA < quantaB) {
            return TRANSFER_TO_A;
        }
        return NO_OP;
    }
}
