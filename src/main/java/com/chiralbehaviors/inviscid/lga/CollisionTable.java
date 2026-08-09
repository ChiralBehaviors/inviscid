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

import java.util.Collections;
import java.util.EnumMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import com.chiralbehaviors.inviscid.lga.CollisionRule.Delta;

/**
 * The deterministic, conservation-exact collision TABLE (bead
 * inviscid-0nx.20, Phase C.3) - the formal-LGA counterpart of the Phase A
 * {@link QuantaExchangeRule}: input state -> output {@link Delta} by table
 * lookup, exactly conserving quanta by construction (every {@link Delta} is
 * validated zero-sum at construction time - see {@link
 * CollisionRule.Delta}'s canonical constructor).
 *
 * <h2>State-space census (why the key domain is exactly 3 elements)</h2>
 * {@link CollisionRule#resolve} is declared to take a {@code Contact} plus
 * two {@code (angle, quanta)} pairs, but a direct reading of {@link
 * QuantaExchangeRule#resolve QuantaExchangeRule.resolve(...)} (the ONLY
 * Phase A rule this bead transcribes) shows it reads none of that: the
 * {@code contact}, {@code angleA}, and {@code angleB} parameters are
 * declared but never referenced anywhere in the method body, which is a
 * pure 3-way comparison of {@code quantaA} vs {@code quantaB} ({@code >},
 * {@code <}, {@code ==}) and nothing else. (Deliberately cited by method
 * name only, not line numbers - line references drift as either file's
 * Javadoc grows; see that method's own source for the current text.)
 * {@link QuantaExchangeRuleTest} confirms this empirically (varying
 * angle/contact never changes the decision). Consequently the table's key
 * domain is NOT the raw, unbounded {@code long} quanta counts, and does NOT
 * include phase bins or contact identity (both of which {@link
 * ContactTable} / {@link PhaseQuantizer} would supply) - it is exactly the
 * ternary classification {@link QuantaOrdering}, domain size 3. (This
 * mirrors, but is not itself authorized by, this project's general
 * "minimal abstraction - extract only when repetition is proven" stance;
 * the actual justification is the dead-parameter reading above, not any
 * external mandate.)
 *
 * <p><b>Load-bearing consequence for bead inviscid-0nx.21:</b> this table
 * alone cannot decide WHETHER two members collide - only WHAT HAPPENS given
 * that they do. The synchronous update must consult {@link ContactTable}
 * (or an equivalent contact predicate) separately to decide eligibility,
 * then this table to decide the outcome.
 *
 * <h2>Inherited completeness floor (acknowledged per inviscid-0nx.19
 * substantive-critique note)</h2>
 * This class's own conservation-exactness is unconditional - every {@link
 * Delta} this table can return is validated zero-sum at construction, for
 * every state it does return. It never calls {@link ContactTable}. But the
 * larger formal-LGA pipeline's implicit claim - "every physical contact
 * gets a corresponding collision opportunity" - inherits {@link
 * ContactTable}'s 1-degree (360-step) angular-grid completeness floor: both
 * generation-time sweeps behind {@link ContactTable} ({@code
 * ContactComboCache.sweepExhaustively}, {@code
 * ContactAtlasGenerator.sweepOverlapAndCenter}) are exhaustive AT that
 * resolution, not literally continuous, so a contact ribbon narrower than
 * 1 degree could in principle be missed upstream of wherever this table is
 * consulted. This does not weaken THIS class's own conservation guarantee
 * (a missed contact simply means no collision opportunity arises there, not
 * that quanta are created or destroyed), but it bounds what the larger
 * pipeline can honestly claim.
 *
 * <h2>Semi-detailed balance: ACCEPTED VIOLATION (user decision 2026-08-08,
 * recorded verbatim in bead inviscid-0nx.20's notes)</h2>
 * {@link #checkSemiDetailedBalance(long)} computes, in closed form and
 * verified empirically (see that method's Javadoc), the exact preimage-count
 * profile of this rule family's induced map on {@code (quantaA, quantaB)}
 * pairs: {@code 3} preimages at output-diff {@code 0}; {@code 2} at
 * {@code |output-diff| == 1}; {@code 1} (the balanced value) everywhere
 * {@code |output-diff| >= 2}. A deterministic map satisfies semi-detailed
 * balance iff every target has exactly one preimage, so this profile is a
 * genuine, exact VIOLATION, concentrated at and near exact quanta equality.
 * <p>
 * The user's 2026-08-08 decision on this finding: ACCEPT the violation as
 * intended dissipative-rule behavior - it is the expected signature of a
 * relaxational, H-theorem-style rule (as {@link QuantaExchangeRule}'s own
 * Javadoc already claims: "relaxation toward a diffusive coarse limit"),
 * not a reversible/measure-preserving classical-LGA collision built to keep
 * an equilibrium distribution sane. The v1 rule form ({@link
 * QuantaExchangeRule}, locked bead inviscid-0nx.14) stays AS-IS - no
 * redesign was authorized or attempted. A mandatory flag on this finding is
 * separately registered on bead inviscid-0nx.23's gate (not this class's
 * concern; noted here only for cross-reference). Any future bead computing
 * equilibrium statistics off this collision process (bead inviscid-0nx.22)
 * must account for this non-uniform preimage profile - an equilibrium
 * distribution derived assuming semi-detailed balance would be wrong here.
 * <p>
 * {@code CollisionTableTest.semiDetailedBalanceIsCheckedAndReported} pins
 * this exact profile (not merely "balance holds/fails") so that ANY future
 * change to the rule's dissipative asymmetry - in EITHER direction,
 * including an accidental fix toward balance - fails loudly rather than
 * passing silently.
 *
 * @author halhildebrand
 */
public final class CollisionTable {

    /**
     * The ternary classification of a {@code (quantaA, quantaB)} pair that
     * is, by census (see class Javadoc), the ENTIRE input {@link
     * QuantaExchangeRule} actually reads. {@link #classify(long, long)}
     * freezes the exact same comparison, in the exact same order, as
     * {@link QuantaExchangeRule#resolve} - the tie-break convention
     * "equal quanta is a defined no-op" is pinned from that class, not
     * reinvented here.
     */
    public enum QuantaOrdering {
        A_GREATER, B_GREATER, EQUAL;

        /**
         * @return the {@link QuantaOrdering} for {@code (quantaA, quantaB)}
         *         - mirrors {@link QuantaExchangeRule#resolve}'s own
         *         three-way comparison exactly.
         */
        public static QuantaOrdering classify(long quantaA, long quantaB) {
            if (quantaA > quantaB) {
                return A_GREATER;
            }
            if (quantaA < quantaB) {
                return B_GREATER;
            }
            return EQUAL;
        }
    }

    /**
     * The recorded result of a semi-detailed-balance check (bead
     * inviscid-0nx.20 locked decision: compute and record, never
     * auto-repair). {@code preimageCountByOutputDiff} tallies, over the
     * INTERIOR of the scanned window (see {@link
     * #checkSemiDetailedBalance(long)}), how many distinct pre-collision
     * {@code (quantaA, quantaB)} pairs mapped to a post-collision pair
     * whose {@code quantaA - quantaB} equals the map's key - a deterministic
     * map satisfies semi-detailed balance iff every entry equals exactly 1
     * (i.e. the map is a bijection: every reachable state has exactly one
     * preimage). {@code balanced} is that condition, checked exhaustively
     * over the scanned interior.
     */
    public record SemiDetailedBalanceReport(long window, long statesChecked,
                                             SortedMap<Long, Long> preimageCountByOutputDiff,
                                             boolean balanced, String summary) {

        public SemiDetailedBalanceReport {
            preimageCountByOutputDiff = Collections.unmodifiableSortedMap(new TreeMap<>(preimageCountByOutputDiff));
        }
    }

    /**
     * Enforced ceiling on {@link #checkSemiDetailedBalance(long)}'s
     * {@code window} - see that method's {@code @param window} Javadoc.
     * Deliberately small enough that {@code
     * checkSemiDetailedBalanceRejectsInvalidWindow}'s boundary test
     * (exercising {@code window == MAX_WINDOW} itself, which runs the
     * full {@code O(window^2)} scan rather than short-circuiting) stays
     * fast.
     */
    public static final long MAX_WINDOW = 1_000L;

    private final Map<QuantaOrdering, Delta> entries;

    CollisionTable(Map<QuantaOrdering, Delta> entries) {
        EnumMap<QuantaOrdering, Delta> defensive = new EnumMap<>(QuantaOrdering.class);
        defensive.putAll(entries);
        this.entries = Collections.unmodifiableMap(defensive);
    }

    /**
     * Convenience: builds the canonical table by enumerating the full
     * {@link QuantaOrdering} input space and applying {@code rule} - see
     * {@link CollisionTableBuilder#fromPhaseARule(CollisionRule)} for the
     * representative-input convention.
     */
    public static CollisionTable buildFromPhaseARule(CollisionRule rule) {
        return CollisionTableBuilder.fromPhaseARule(rule).build();
    }

    /**
     * @return the frozen {@link Delta} for {@code (quantaA, quantaB)}'s
     *         {@link QuantaOrdering} - never {@code null}, throws rather
     *         than silently no-op-ing if this table was built without an
     *         entry for that class (see {@link #lookup(QuantaOrdering)}).
     */
    public Delta lookup(long quantaA, long quantaB) {
        return lookup(QuantaOrdering.classify(quantaA, quantaB));
    }

    /**
     * @return the frozen {@link Delta} for {@code key}
     * @throws IllegalStateException
     *             if this table has no entry for {@code key} - a missing
     *             entry for a legal {@link QuantaOrdering} must THROW,
     *             never silently no-op (bead inviscid-0nx.20 test 2).
     */
    public Delta lookup(QuantaOrdering key) {
        Delta delta = entries.get(key);
        if (delta == null) {
            throw new IllegalStateException("CollisionTable has no entry for " + key
                                             + " - a frozen table must be total over its domain");
        }
        return delta;
    }

    /**
     * @return an unmodifiable snapshot of every {@link QuantaOrdering} this
     *         table has an entry for
     */
    public Map<QuantaOrdering, Delta> entries() {
        return entries;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof CollisionTable other)) {
            return false;
        }
        return entries.equals(other.entries);
    }

    @Override
    public int hashCode() {
        return entries.hashCode();
    }

    /**
     * Computes and RECORDS (never auto-repairs) a semi-detailed-balance
     * check over this table's induced dynamics on the real {@code
     * (quantaA, quantaB)} pair state space - NOT over the abstracted
     * 3-element {@link QuantaOrdering} key domain, which is not closed
     * under the rule (e.g. {@code A_GREATER} with a gap of exactly 1
     * transitions to {@code B_GREATER}, with a gap of exactly 2 transitions
     * to {@code EQUAL}, and with a gap of 3+ stays {@code A_GREATER} - the
     * successor class depends on the actual gap magnitude, not just its
     * sign, so only the real pair space admits a well-defined state->state
     * transition function).
     *
     * <h2>Method</h2>
     * Scans every {@code (quantaA, quantaB)} with both coordinates in
     * {@code [-window, window]}, applies this table's own lookup to obtain
     * each pair's image, and tallies - for every image whose coordinates
     * both fall in the INTERIOR {@code [-(window-1), window-1]} - how many
     * distinct source pairs produced it. Restricting the tally to the
     * interior is not a sampling shortcut: because every {@link Delta} this
     * rule family produces has {@code |deltaA|, |deltaB| <= 1}, any true
     * preimage of an interior image differs from it by at most 1 in each
     * coordinate and is therefore GUARANTEED to already be inside the
     * scanned {@code [-window, window]} square - so interior tallies are
     * exact, not edge-truncated. {@code window} is a representative
     * sampling bound (this rule's true domain is all of {@code long x
     * long} and cannot be swept exhaustively), matching the {@code
     * QUANTA_RANGE}-style convention already established in {@code
     * QuantaExchangeRuleTest}.
     *
     * <h2>Closed-form cross-check (recorded in bead notes)</h2>
     * Writing {@code s = quantaA + quantaB} (conserved by every {@link
     * Delta}) and {@code d = quantaA - quantaB}, the rule's action on
     * {@code d} alone is: {@code d > 0 -> d - 2}; {@code d < 0 -> d + 2};
     * {@code d == 0 -> d} (fixed). Solving each branch's inverse and
     * checking its own branch-validity condition gives, for every target
     * {@code d'}: exactly 3 preimages at {@code d' == 0}; exactly 2 at
     * {@code |d'| == 1}; exactly 1 at {@code |d'| >= 2} - independent of
     * {@code s} by translation symmetry, and independent of {@code window}
     * once {@code window} is large enough to avoid interior/edge collision
     * (any {@code window >= 3} suffices). This method's empirical scan is
     * expected to reproduce that distribution exactly.
     *
     * @param window
     *            the representative scan half-width, {@code > 0} and
     *            {@code <= MAX_WINDOW} ({@code 1,000}) - the scan is
     *            {@code O(window^2)} pairs, each allocating an {@link
     *            ImagePair} tally-map entry, so an unbounded caller-
     *            supplied window is an easy way to request an accidental
     *            multi-billion-entry scan; {@code 1,000} (~4M pairs, well
     *            over the {@code window >= 3} the closed-form proof needs
     *            to avoid interior/edge collision, yet still fast enough
     *            that a boundary test can exercise {@code window ==
     *            MAX_WINDOW} directly) is a generous ceiling enforced
     *            rather than merely documented, per both reviewers'
     *            follow-up
     * @return the recorded report; {@link SemiDetailedBalanceReport#balanced()}
     *         is {@code true} only if EVERY interior target has exactly one
     *         preimage - per the closed-form above, this rule family does
     *         NOT satisfy that. This is an ACCEPTED finding (user decision
     *         2026-08-08, recorded in bead inviscid-0nx.20's notes; see
     *         class Javadoc), not a defect in this method.
     */
    public SemiDetailedBalanceReport checkSemiDetailedBalance(long window) {
        if (window <= 0) {
            throw new IllegalArgumentException("window must be > 0: " + window);
        }
        if (window > MAX_WINDOW) {
            throw new IllegalArgumentException("window must be <= " + MAX_WINDOW + " (O(window^2) scan cost): "
                                                + window);
        }
        Map<ImagePair, Long> imageCoordCounts = new LinkedHashMap<>();
        Map<Long, Long> preimageCountByDiff = new TreeMap<>();
        long interior = window - 1;
        long statesChecked = 0;
        for (long a = -window; a <= window; a++) {
            for (long b = -window; b <= window; b++) {
                Delta delta = lookup(a, b);
                long aPrime = a + delta.deltaA();
                long bPrime = b + delta.deltaB();
                if (aPrime < -interior || aPrime > interior || bPrime < -interior
                    || bPrime > interior) {
                    continue;
                }
                statesChecked++;
                imageCoordCounts.merge(new ImagePair(aPrime, bPrime), 1L, Long::sum);
            }
        }
        boolean balanced = statesChecked > 0;
        for (Map.Entry<ImagePair, Long> e : imageCoordCounts.entrySet()) {
            long dPrime = e.getKey().a() - e.getKey().b();
            long count = e.getValue();
            Long existing = preimageCountByDiff.get(dPrime);
            if (existing != null && existing.longValue() != count) {
                // Internal-consistency guard: the closed-form proof claims
                // preimage count depends only on d', not on s' - if that
                // were ever violated the check itself would be broken, so
                // fail loudly rather than silently averaging/overwriting.
                throw new IllegalStateException("preimage count for output diff " + dPrime
                                                 + " is not uniform across s' - closed-form proof violated ("
                                                 + existing + " vs " + count + ")");
            }
            preimageCountByDiff.put(dPrime, count);
            if (count != 1L) {
                balanced = false;
            }
        }
        String summary = "window=" + window + " statesChecked=" + statesChecked
                          + " preimageCountByOutputDiff=" + preimageCountByDiff
                          + " balanced=" + balanced;
        return new SemiDetailedBalanceReport(window, statesChecked, new TreeMap<>(preimageCountByDiff),
                                              balanced, summary);
    }

    /** A post-collision {@code (a, b)} pair, used only as a tally-map key. */
    private record ImagePair(long a, long b) {
    }
}
