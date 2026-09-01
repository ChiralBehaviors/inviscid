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

package com.chiralbehaviors.inviscid.automaton.lga;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.SortedMap;
import java.util.TreeMap;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.automaton.lga.CollisionRule.Delta;
import com.chiralbehaviors.inviscid.automaton.lga.CollisionTable.QuantaOrdering;
import com.chiralbehaviors.inviscid.automaton.lga.CollisionTable.SemiDetailedBalanceReport;

/**
 * Behavioral tests for {@link CollisionTable} / {@link CollisionTableBuilder}
 * (bead inviscid-0nx.20, Phase C.3).
 *
 * @author halhildebrand
 */
public class CollisionTableTest {

    private static final long    QUANTA_RANGE = 8L;
    private static final Point3i CELL_A       = new Point3i(0, 0, 0);
    private static final Point3i CELL_B       = new Point3i(1, 0, -1);

    private static Contact fixtureContact() {
        return new Contact(CELL_A, 3, 1, CELL_B, 3, 0, 1, 0.0);
    }

    private static CollisionTable canonicalTable() {
        return CollisionTable.buildFromPhaseARule(new QuantaExchangeRule());
    }

    /**
     * Test 1: exhaustive over the whole table (domain size 3), exact
     * integer conservation. Note (per relay protocol, honest self-flag):
     * {@link Delta}'s own canonical constructor already rejects any
     * non-zero-summing pair at construction time, so no entry in ANY
     * {@link CollisionTable} could ever fail this - this test is a pin
     * guarding against a builder bug substituting the wrong (but
     * individually valid) {@link Delta} for a key, not a test that could
     * observe a genuine conservation violation.
     */
    @Test
    public void everyEntryConservesQuantaExactly() {
        CollisionTable table = canonicalTable();
        assertEquals("canonical table must be total over its 3-element domain", 3,
                     table.entries().size());
        for (QuantaOrdering key : QuantaOrdering.values()) {
            Delta delta = table.lookup(key);
            assertEquals("delta for " + key + " must conserve quanta exactly", 0L,
                         delta.deltaA() + delta.deltaB());
        }
    }

    /**
     * Test 2: a missing entry for a legal state must THROW, never silently
     * no-op. Also pins the positive case: the canonical (complete) table
     * never throws for any legal {@link QuantaOrdering}.
     */
    @Test
    public void tableIsTotalOverItsDomain() {
        CollisionTable incomplete = new CollisionTableBuilder().put(QuantaOrdering.A_GREATER,
                                                                      new Delta(-1L, 1L))
                                                                 .put(QuantaOrdering.B_GREATER,
                                                                      new Delta(1L, -1L))
                                                                 .build();
        assertThrows(IllegalStateException.class, () -> incomplete.lookup(QuantaOrdering.EQUAL));

        CollisionTable canonical = canonicalTable();
        for (QuantaOrdering key : QuantaOrdering.values()) {
            // must not throw
            canonical.lookup(key);
        }
    }

    /**
     * Test 3 (the A -> C consistency anchor): {@link CollisionTable}'s
     * output must equal {@link QuantaExchangeRule}'s output for every
     * shared-representable input. The census in {@link CollisionTable}'s
     * class Javadoc establishes {@link QuantaExchangeRule}'s decision
     * depends ONLY on {@code sign(quantaA - quantaB)} - so this equality
     * provably holds for ALL of {@code long x long}, not just what is
     * sampled below. The domain (all pairs of {@code long}) does not admit
     * a literal exhaustive sweep, so this test samples: the full {@code
     * [-8, 8]} grid (matching {@code QuantaExchangeRuleTest}'s own
     * QUANTA_RANGE convention, 289 pairs, genuinely exhaustive over that
     * sub-grid) plus the extremes {@code QuantaExchangeRuleTest} itself
     * uses to pin overflow-freedom ({@code Long.MIN_VALUE}, {@code
     * Long.MAX_VALUE}, the float32-exact-integer boundary {@code 2^24}).
     */
    @Test
    public void tableReproducesPhaseARuleOnSharedInputs() {
        QuantaExchangeRule rule = new QuantaExchangeRule();
        CollisionTable table = CollisionTable.buildFromPhaseARule(rule);
        Contact contact = fixtureContact();

        for (long quantaA = -QUANTA_RANGE; quantaA <= QUANTA_RANGE; quantaA++) {
            for (long quantaB = -QUANTA_RANGE; quantaB <= QUANTA_RANGE; quantaB++) {
                Delta expected = rule.resolve(contact, 0.3f, quantaA, 1.1f, quantaB);
                Delta actual = table.lookup(quantaA, quantaB);
                assertEquals("table must reproduce QuantaExchangeRule for quantaA=" + quantaA
                             + " quantaB=" + quantaB, expected, actual);
            }
        }

        long floatExactLimit = 1L << 24;
        long[][] extremes = { { floatExactLimit, floatExactLimit - 2L },
                               { Long.MAX_VALUE, Long.MIN_VALUE },
                               { Long.MIN_VALUE, Long.MAX_VALUE },
                               { Long.MAX_VALUE, Long.MAX_VALUE },
                               { Long.MIN_VALUE, Long.MIN_VALUE } };
        for (long[] pair : extremes) {
            Delta expected = rule.resolve(contact, 0f, pair[0], 0f, pair[1]);
            Delta actual = table.lookup(pair[0], pair[1]);
            assertEquals("table must reproduce QuantaExchangeRule at extreme quantaA=" + pair[0]
                         + " quantaB=" + pair[1], expected, actual);
        }
    }

    /**
     * Non-vacuity: an identity table conserves everything perfectly and
     * does nothing, which would trivially pass every other test in this
     * class.
     */
    @Test
    public void tableIsNotTheIdentity() {
        CollisionTable table = canonicalTable();
        Delta aGreater = table.lookup(QuantaOrdering.A_GREATER);
        Delta bGreater = table.lookup(QuantaOrdering.B_GREATER);
        assertNotEquals("A_GREATER must produce a nonzero transfer", 0L, aGreater.deltaA());
        assertNotEquals(0L, aGreater.deltaB());
        assertNotEquals("B_GREATER must produce a nonzero transfer", 0L, bGreater.deltaA());
        assertNotEquals(0L, bGreater.deltaB());
    }

    /**
     * Test 5, RESHAPED per the 2026-08-08 user decision recorded in bead
     * inviscid-0nx.20's notes: the semi-detailed-balance violation is
     * ACCEPTED as intended dissipative-rule behavior (v1's {@code
     * QuantaExchangeRule} form stays as-is, no redesign; a mandatory flag
     * is separately registered on bead inviscid-0nx.23's gate). This test
     * therefore no longer asserts balance holds - it asserts three things,
     * ALL of which must remain true: (a) the check RAN and produced a
     * non-vacuous, recorded report; (b) the recorded preimage-count profile
     * matches the {@link CollisionTable#checkSemiDetailedBalance(long)}
     * closed-form EXACTLY, bucket-for-bucket, so any future change to
     * either the rule or the check - in EITHER direction, including an
     * accidental fix that makes it balanced - fails this test loudly; and
     * (c) the violation is explicitly, deliberately asserted as a known,
     * accepted property ({@code balanced() == false}), not silently
     * ignored.
     */
    @Test
    public void semiDetailedBalanceIsCheckedAndReported() {
        CollisionTable table = canonicalTable();
        long window = 12L;
        SemiDetailedBalanceReport report = table.checkSemiDetailedBalance(window);

        // (a) the check ran and recorded a non-vacuous result.
        assertNotNull("semi-detailed balance check must produce a report", report);
        assertEquals(window, report.window());
        assertTrue("check must have scanned a non-vacuous number of states",
                   report.statesChecked() > 0);
        assertNotNull(report.summary());
        assertTrue("summary must be non-empty", !report.summary().isEmpty());

        // (b) exact profile match against the closed form (see
        // CollisionTable#checkSemiDetailedBalance Javadoc): 3 preimages at
        // output-diff 0, 2 at |output-diff|==1, 1 everywhere else in the
        // scanned interior [-(window-1), window-1]. Built from the formula,
        // not hand-transcribed, so it stays correct if `window` changes -
        // but any CHANGE IN SHAPE of the actual profile (e.g. the rule
        // becoming reversible/bijective, or a different asymmetry) still
        // fails this assertion, which is the point: this test must fail
        // loudly if the balance profile ever changes, in either direction.
        // NOTE: aPrime and bPrime are each independently bounded to the
        // interior [-(window-1), window-1], so their DIFFERENCE ranges over
        // twice that: [-2*(window-1), 2*(window-1)].
        long interior = window - 1;
        SortedMap<Long, Long> expected = new TreeMap<>();
        for (long dPrime = -2 * interior; dPrime <= 2 * interior; dPrime++) {
            long expectedCount = dPrime == 0 ? 3L : Math.abs(dPrime) == 1 ? 2L : 1L;
            expected.put(dPrime, expectedCount);
        }
        assertEquals("recorded preimage-count profile must match the closed form exactly - bucket for"
                     + " bucket - so any change to the rule's dissipative asymmetry (including an"
                     + " accidental fix toward balance) is caught: " + report.summary(), expected,
                     report.preimageCountByOutputDiff());

        // (c) the violation is EXPLICITLY asserted as a known, ACCEPTED
        // property - user decision 2026-08-08 (bead inviscid-0nx.20 notes):
        // ACCEPT the semi-detailed-balance violation as intended
        // dissipative-rule behavior; v1's QuantaExchangeRule form stays
        // as-is, not redesigned; a mandatory flag is separately registered
        // on bead inviscid-0nx.23's gate. This assertion is deliberately
        // the OPPOSITE polarity of a plain "balance must hold" check - if
        // this rule ever became balanced, THIS line would fail, which is
        // the intended catch for an unnoticed rule-form change.
        assertTrue("semi-detailed balance is expected (and accepted, per the 2026-08-08 user"
                   + " decision) to be VIOLATED - v1's QuantaExchangeRule is a dissipative"
                   + " relaxation rule, not a reversible/measure-preserving one; if this ever"
                   + " becomes true, the rule's dissipative character changed and that is itself"
                   + " reportable: " + report.summary(), !report.balanced());
    }

    /**
     * Test 6: no RNG in lookup (repeated lookups return equal values); the
     * builder itself is also deterministic (build twice, compare).
     */
    @Test
    public void tableIsDeterministic() {
        CollisionTable table = canonicalTable();
        for (QuantaOrdering key : QuantaOrdering.values()) {
            Delta first = table.lookup(key);
            Delta second = table.lookup(key);
            assertEquals(first, second);
        }

        CollisionTable first = CollisionTable.buildFromPhaseARule(new QuantaExchangeRule());
        CollisionTable second = CollisionTable.buildFromPhaseARule(new QuantaExchangeRule());
        assertEquals("independently built tables must be structurally equal (same entries map)", first,
                     second);
        assertEquals(first.entries(), second.entries());
    }

    /**
     * Direct coverage for {@link CollisionTable#lookup(long, long)}'s
     * missing-entry throw path (previously only exercised transitively via
     * {@link CollisionTable#lookup(CollisionTable.QuantaOrdering)} in
     * {@link #tableIsTotalOverItsDomain()}).
     */
    @Test
    public void lookupByQuantaThrowsForMissingEntry() {
        CollisionTable incomplete = new CollisionTableBuilder().put(QuantaOrdering.A_GREATER,
                                                                      new Delta(-1L, 1L))
                                                                 .build();
        // (0, 0) classifies as EQUAL, which this deliberately-incomplete
        // table has no entry for.
        assertThrows(IllegalStateException.class, () -> incomplete.lookup(0L, 0L));
    }

    /**
     * Code-review follow-up (inviscid-0nx.20): {@link
     * CollisionTableBuilder#fromPhaseARule(CollisionRule)}'s one-
     * representative-per-{@link QuantaOrdering} sampling is sound ONLY for
     * a rule whose {@link Delta} decision depends purely on the quanta
     * comparison class - never on the specific magnitudes, and never on
     * angle. A rule that reads angle (like this fixture, which ignores
     * quanta entirely and decides from {@code angleA <= angleB}) would, in
     * the absence of a soundness probe, be silently baked into a WRONG
     * table (whichever angle the single representative happened to carry
     * would freeze that arbitrary answer for the whole class). {@link
     * CollisionTableBuilder#fromPhaseARule(CollisionRule)} must instead
     * reject such a rule loudly, naming the offending rule class.
     *
     * <p>This fixture is angle-dependent for EVERY {@link QuantaOrdering}
     * class, and {@link CollisionTableBuilder#fromPhaseARule} probes
     * classes in a fixed, fail-fast sequential order (A_GREATER first) -
     * so this test specifically exercises the A_GREATER probe. See {@link
     * #probeRejectsRuleDependentOnAngleOnlyAtBGreater()} and {@link
     * #probeRejectsRuleDependentOnAngleOnlyAtEqual()} for fixtures proving
     * the B_GREATER and EQUAL probes independently fire too (substantive-
     * critique/code-review coverage-gap follow-up).
     */
    @Test
    public void probeRejectsAngleDependentRule() {
        CollisionRule angleSensitive = new CollisionRule() {
            @Override
            public Delta resolve(Contact contact, float angleA, long quantaA, float angleB,
                                  long quantaB) {
                // Deliberately NOT quanta-comparison-pure: ignores quanta
                // entirely, decides from angle alone.
                return angleA <= angleB ? new Delta(-1L, 1L) : new Delta(1L, -1L);
            }
        };
        IllegalArgumentException thrown = assertThrows(IllegalArgumentException.class,
                                                         () -> CollisionTableBuilder.fromPhaseARule(angleSensitive));
        assertTrue("exception must name the offending rule class: " + thrown.getMessage(),
                   thrown.getMessage().contains(angleSensitive.getClass().getName()));
    }

    /**
     * Coverage-gap follow-up: a fixture that is quanta-comparison-pure at
     * A_GREATER (so that probe passes and the builder proceeds), but
     * angle-dependent ONLY at B_GREATER - proving the B_GREATER probe
     * itself independently fires, not merely that SOME probe fires
     * somewhere in the sequence.
     */
    @Test
    public void probeRejectsRuleDependentOnAngleOnlyAtBGreater() {
        CollisionRule sensitiveOnlyAtBGreater = new CollisionRule() {
            @Override
            public Delta resolve(Contact contact, float angleA, long quantaA, float angleB,
                                  long quantaB) {
                if (quantaA > quantaB) {
                    return new Delta(-1L, 1L);
                }
                if (quantaA == quantaB) {
                    return Delta.noop();
                }
                // quantaA < quantaB (B_GREATER): angle-dependent here only.
                return angleA <= angleB ? new Delta(1L, -1L) : new Delta(-1L, 1L);
            }
        };
        IllegalArgumentException thrown = assertThrows(IllegalArgumentException.class,
                                                         () -> CollisionTableBuilder.fromPhaseARule(sensitiveOnlyAtBGreater));
        assertTrue("exception must name the B_GREATER class specifically: " + thrown.getMessage(),
                   thrown.getMessage().contains(QuantaOrdering.B_GREATER.toString()));
    }

    /**
     * Coverage-gap follow-up: a fixture that is quanta-comparison-pure at
     * A_GREATER and B_GREATER (both probes pass), but angle-dependent
     * ONLY at EQUAL - proving the EQUAL probe (last in sequence) itself
     * independently fires.
     */
    @Test
    public void probeRejectsRuleDependentOnAngleOnlyAtEqual() {
        CollisionRule sensitiveOnlyAtEqual = new CollisionRule() {
            @Override
            public Delta resolve(Contact contact, float angleA, long quantaA, float angleB,
                                  long quantaB) {
                if (quantaA > quantaB) {
                    return new Delta(-1L, 1L);
                }
                if (quantaA < quantaB) {
                    return new Delta(1L, -1L);
                }
                // tied (EQUAL): angle-dependent here only, instead of the
                // defined no-op.
                return angleA <= angleB ? new Delta(-1L, 1L) : new Delta(1L, -1L);
            }
        };
        IllegalArgumentException thrown = assertThrows(IllegalArgumentException.class,
                                                         () -> CollisionTableBuilder.fromPhaseARule(sensitiveOnlyAtEqual));
        assertTrue("exception must name the EQUAL class specifically: " + thrown.getMessage(),
                   thrown.getMessage().contains(QuantaOrdering.EQUAL.toString()));
    }

    /**
     * Code-review + critic follow-up (inviscid-0nx.20): the soundness
     * probe previously passed the IDENTICAL {@link Contact} instance to
     * both probe calls, so a contact-dependent rule went undetected - the
     * same failure mode the probe exists to catch for angle. This fixture
     * ignores quanta and angle entirely, deciding purely from {@code
     * contact.direction() > 0}; it must be rejected once the second probe
     * uses a genuinely different {@link Contact}.
     */
    @Test
    public void probeRejectsContactDependentRule() {
        CollisionRule contactSensitive = new CollisionRule() {
            @Override
            public Delta resolve(Contact contact, float angleA, long quantaA, float angleB,
                                  long quantaB) {
                return contact.direction() > 0 ? new Delta(-1L, 1L) : new Delta(1L, -1L);
            }
        };
        IllegalArgumentException thrown = assertThrows(IllegalArgumentException.class,
                                                         () -> CollisionTableBuilder.fromPhaseARule(contactSensitive));
        assertTrue("exception must name the offending rule class: " + thrown.getMessage(),
                   thrown.getMessage().contains(contactSensitive.getClass().getName()));
    }

    /**
     * Critic follow-up: direct boundary coverage for {@link
     * CollisionTable#checkSemiDetailedBalance(long)}'s {@code window}
     * validation - {@code window <= 0} and {@code window >
     * CollisionTable#MAX_WINDOW} must both throw.
     */
    @Test
    public void checkSemiDetailedBalanceRejectsInvalidWindow() {
        CollisionTable table = canonicalTable();
        assertThrows(IllegalArgumentException.class, () -> table.checkSemiDetailedBalance(0L));
        assertThrows(IllegalArgumentException.class, () -> table.checkSemiDetailedBalance(-1L));
        assertThrows(IllegalArgumentException.class,
                     () -> table.checkSemiDetailedBalance(CollisionTable.MAX_WINDOW + 1));
        // must NOT throw at the boundary itself
        table.checkSemiDetailedBalance(CollisionTable.MAX_WINDOW);
    }
}
