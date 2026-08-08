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
 * A discrete, conservation-exact collision rule (bead inviscid-0nx.14):
 * given a {@link Contact} and the two participating members' (angle,
 * quanta) state, decides an exact zero-sum quanta transfer between them.
 *
 * <p>Pure function of discrete inputs by design: {@link #resolve} takes
 * only value types (a {@code Contact} plus two (angle, quanta) pairs) and
 * returns a value type ({@link Delta}) with no side effects on the
 * automaton - callers apply the delta themselves, exclusively through
 * {@code Necronomata}'s {@code deltaF} accumulator (never by writing
 * {@code frequency}/{@code angle} directly - see {@code
 * Necronomata.process(Necronomata.Processor)}'s Javadoc for that
 * contract). This shape is deliberate: Phase C (bead inviscid-0nx.20)
 * replaces the hybrid stage's live geometric contact scan with a
 * precomputed lookup table keyed on quantized inputs, and a pure function
 * of (contact identity, angleA, quantaA, angleB, quantaB) is exactly what
 * table-izes - no {@code CollisionRule} implementation needs to change
 * shape between the hybrid and formal stages, only how its inputs are
 * produced (live geometry now, LUT lookup later).
 *
 * <p><b>{@code angle} is part of the signature but unused by v1.</b>
 * {@link QuantaExchangeRule} decides purely from the two quanta counts;
 * angle is carried here so a future angle-dependent rule (e.g. one that
 * weights transfer by approach geometry) does not require an interface
 * change.
 *
 * @author halhildebrand
 */
public interface CollisionRule {

    /**
     * An exact zero-sum quanta transfer: {@code deltaA} is added to
     * member A's quanta count, {@code deltaB} to member B's. The
     * canonical constructor rejects any pair that does not sum to
     * exactly zero - conservation is enforced at the type level, not
     * merely by rule-author convention.
     */
    record Delta(long deltaA, long deltaB) {

        public Delta {
            if (deltaA + deltaB != 0) {
                throw new IllegalArgumentException("deltaA + deltaB must be exactly zero (conservation), was: "
                                                     + deltaA + " + " + deltaB
                                                     + " = " + (deltaA
                                                                 + deltaB));
            }
        }

        /** The no-op delta: both members unchanged. */
        public static Delta noop() {
            return new Delta(0L, 0L);
        }
    }

    /**
     * Resolves one contact into an exact zero-sum quanta transfer.
     *
     * <p><b>Rule-level symmetry combines with sweep-level
     * order-independence for full tick determinism, but is not a
     * substitute for it.</b> An implementation's own participant-swap
     * symmetry (swapping which member is "A" mirrors the decision - see
     * {@code QuantaExchangeRuleTest.ruleIsSymmetricUnderParticipantSwap})
     * only guarantees that a SINGLE contact's outcome does not depend on
     * which endpoint is labeled A vs B. It says nothing about a tick with
     * MULTIPLE contacts touching the same member - that is {@link
     * CollisionSweep}'s responsibility, via its snapshot-resolution
     * contract (see that class's Javadoc, "Snapshot resolution within a
     * tick"). Both are required for full determinism: a perfectly
     * symmetric rule fed through a scan-order-dependent sweep would still
     * exhibit scan-order artifacts, and a snapshot-resolving sweep fed an
     * asymmetric rule would still exhibit endpoint-labeling artifacts.
     *
     * @param contact the geometric contact being resolved (carries
     *                cell/cube/member/direction identity; this method
     *                does not itself read quanta or angle from it - the
     *                caller supplies those explicitly so the rule remains
     *                a pure function of its arguments)
     * @param angleA  member A's current angle, radians
     * @param quantaA member A's current quanta count
     * @param angleB  member B's current angle, radians
     * @param quantaB member B's current quanta count
     * @return an exact zero-sum {@link Delta}
     */
    Delta resolve(Contact contact, float angleA, long quantaA, float angleB,
                  long quantaB);
}
