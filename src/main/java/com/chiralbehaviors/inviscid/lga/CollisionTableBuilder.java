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

import java.util.EnumMap;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.lga.CollisionRule.Delta;
import com.chiralbehaviors.inviscid.lga.CollisionTable.QuantaOrdering;

/**
 * Mutable builder for a {@link CollisionTable} (bead inviscid-0nx.20, Phase
 * C.3): enumerate the input state space, apply a {@link CollisionRule} (or
 * hand-supply entries directly), then {@link #build()} freezes an
 * immutable, total-over-its-domain {@link CollisionTable}.
 *
 * <p>{@link #build()} does NOT require every {@link QuantaOrdering} to have
 * an entry - an intentionally incomplete builder is a legitimate way to
 * construct a table that throws on the missing key (see {@link
 * CollisionTable#lookup(QuantaOrdering)}'s Javadoc and {@code
 * CollisionTableTest.tableIsTotalOverItsDomain}). {@link
 * #fromPhaseARule(CollisionRule)} is the canonical path that always
 * produces a total table.
 *
 * @author halhildebrand
 */
public final class CollisionTableBuilder {

    /**
     * First-probe placeholder {@link Contact} used only to satisfy {@link
     * CollisionRule#resolve}'s signature when transcribing {@link
     * QuantaExchangeRule} - unused by v1's decision logic (see {@link
     * CollisionTable}'s class Javadoc, "State-space census"). Mirrors
     * {@code QuantaExchangeRuleTest.fixtureContact()}'s convention. See
     * {@link #probedResolve} for why a SECOND, deliberately different
     * {@link Contact} is also probed.
     */
    private static final Contact PLACEHOLDER_CONTACT = new Contact(new Point3i(0, 0, 0), 3, 1,
                                                                     new Point3i(1, 0, -1), 3, 0, 1,
                                                                     0.0);

    /**
     * Second-probe {@link Contact} (code-review follow-up, inviscid-
     * 0nx.20): deliberately different cells, cubes, members, and direction
     * from {@link #PLACEHOLDER_CONTACT}, so a rule whose decision actually
     * depends on contact identity diverges between the two {@link
     * #probedResolve} calls and gets caught rather than silently baked
     * into the table - the same failure mode the angle probe exists for,
     * previously unguarded because both calls passed the identical
     * {@link Contact} instance.
     *
     * <p>Two acknowledged limits of the probe surface: {@code minDistance}
     * is 0.0 in BOTH probe constants, so a rule reading only
     * {@code contact.minDistance()} would pass undetected; and
     * {@code direction=-1} here is deliberately outside the {@code +1..+6}
     * domain {@link ContactScan} emits - a synthetic maximum-divergence
     * value, never a production one.
     */
    private static final Contact PROBE_CONTACT_2 = new Contact(new Point3i(2, 1, 0), 4, 2,
                                                                 new Point3i(3, 1, -1), 4, 3, -1,
                                                                 0.0);

    /** First-probe placeholder angle - unused by v1 rules (see {@link
     * CollisionTable} class Javadoc), but see {@link #probedResolve} for
     * why a SECOND, deliberately different angle is also probed. */
    private static final float PLACEHOLDER_ANGLE = 0f;

    /**
     * Second-probe angles (code-review follow-up, inviscid-0nx.20):
     * deliberately different from {@link #PLACEHOLDER_ANGLE} on both
     * members, and from each other, so a rule whose decision actually
     * depends on angle diverges between the two {@link #probedResolve}
     * calls and gets caught rather than silently baked into the table.
     */
    private static final float PROBE_ANGLE_A_2 = 3f;
    private static final float PROBE_ANGLE_B_2 = 1f;

    private final EnumMap<QuantaOrdering, Delta> entries = new EnumMap<>(QuantaOrdering.class);

    /**
     * Builds a builder pre-populated by enumerating every {@link
     * QuantaOrdering} and applying {@code rule} with a representative
     * {@code (quantaA, quantaB)} pair for that class: {@code A_GREATER ->
     * (1, 0)}, {@code B_GREATER -> (0, 1)}, {@code EQUAL -> (0, 0)}.
     *
     * <p><b>Soundness precondition, runtime-enforced (code-review follow-
     * up):</b> this sampling is sound ONLY for a rule whose {@link Delta}
     * decision depends purely on the {@link QuantaOrdering} class - never
     * on the specific quanta magnitude, never on angle, and never on
     * contact identity. A rule that does NOT have that property (e.g. a
     * future angle-dependent or contact-dependent {@link CollisionRule})
     * would otherwise be silently baked into a WRONG table - whichever
     * answer the single representative happened to produce would be
     * frozen for the whole class. {@link #probedResolve} enforces this
     * precondition for every class by resolving a SECOND, distinct
     * representative (different quanta magnitude AND different angle AND
     * a different {@link Contact}) and requiring an identical {@link
     * Delta}; a mismatch throws {@link IllegalArgumentException} naming
     * the offending rule class rather than building a wrong table.
     *
     * @param rule
     *            the Phase A rule form to transcribe (v1: {@link
     *            QuantaExchangeRule})
     * @return a builder with all 3 entries populated, ready for {@link
     *         #build()}
     * @throws IllegalArgumentException
     *             if {@code rule} is not quanta-comparison-pure (see
     *             {@link #probedResolve})
     */
    public static CollisionTableBuilder fromPhaseARule(CollisionRule rule) {
        CollisionTableBuilder builder = new CollisionTableBuilder();
        builder.put(QuantaOrdering.A_GREATER, probedResolve(rule, QuantaOrdering.A_GREATER, 1L, 0L, 5L, 2L));
        builder.put(QuantaOrdering.B_GREATER, probedResolve(rule, QuantaOrdering.B_GREATER, 0L, 1L, 2L, 5L));
        builder.put(QuantaOrdering.EQUAL, probedResolve(rule, QuantaOrdering.EQUAL, 0L, 0L, 7L, 7L));
        return builder;
    }

    /**
     * Resolves {@code rule} against TWO distinct representatives of
     * {@code key} - {@code (quantaA1, quantaB1)} at {@link
     * #PLACEHOLDER_ANGLE} for both members under {@link
     * #PLACEHOLDER_CONTACT}, and {@code (quantaA2, quantaB2)} at {@link
     * #PROBE_ANGLE_A_2} / {@link #PROBE_ANGLE_B_2} under {@link
     * #PROBE_CONTACT_2} - and requires an identical {@link Delta} from
     * both. See {@link #fromPhaseARule(CollisionRule)}'s Javadoc for why:
     * this is the runtime enforcement of the "class-only, never magnitude,
     * never angle, never contact identity" soundness precondition that
     * sampling one representative per class otherwise silently assumes.
     *
     * @throws IllegalArgumentException
     *             naming {@code rule}'s class and both probed inputs, if
     *             the two probes disagree
     */
    private static Delta probedResolve(CollisionRule rule, QuantaOrdering key, long quantaA1, long quantaB1,
                                        long quantaA2, long quantaB2) {
        Delta first = rule.resolve(PLACEHOLDER_CONTACT, PLACEHOLDER_ANGLE, quantaA1, PLACEHOLDER_ANGLE,
                                    quantaB1);
        Delta second = rule.resolve(PROBE_CONTACT_2, PROBE_ANGLE_A_2, quantaA2, PROBE_ANGLE_B_2,
                                     quantaB2);
        if (!first.equals(second)) {
            throw new IllegalArgumentException("CollisionTableBuilder.fromPhaseARule requires a"
                                                + " quanta-comparison-pure rule (its Delta must depend ONLY on"
                                                + " the QuantaOrdering class - never on the specific quanta"
                                                + " magnitude, never on angle, never on contact identity), but "
                                                + rule.getClass().getName()
                                                + " produced different Delta values for two representatives of "
                                                + key + ": (contact=" + PLACEHOLDER_CONTACT + ", angleA="
                                                + PLACEHOLDER_ANGLE + ", quantaA=" + quantaA1 + ", angleB="
                                                + PLACEHOLDER_ANGLE + ", quantaB=" + quantaB1 + ") -> " + first
                                                + " vs (contact=" + PROBE_CONTACT_2 + ", angleA=" + PROBE_ANGLE_A_2
                                                + ", quantaA=" + quantaA2 + ", angleB=" + PROBE_ANGLE_B_2
                                                + ", quantaB=" + quantaB2 + ") -> " + second);
        }
        return first;
    }

    /**
     * @param key
     *            the {@link QuantaOrdering} this entry is for
     * @param delta
     *            the frozen {@link Delta} for {@code key}
     * @return {@code this}, for chaining
     * @throws IllegalStateException
     *             if {@code key} already has an entry - catches accidental
     *             double-population rather than silently overwriting
     */
    public CollisionTableBuilder put(QuantaOrdering key, Delta delta) {
        if (entries.containsKey(key)) {
            throw new IllegalStateException("duplicate entry for " + key);
        }
        entries.put(key, delta);
        return this;
    }

    /**
     * @return an immutable {@link CollisionTable} snapshotting this
     *         builder's current entries - subsequent {@link #put} calls on
     *         this builder do not affect an already-built table
     */
    public CollisionTable build() {
        return new CollisionTable(entries);
    }
}
