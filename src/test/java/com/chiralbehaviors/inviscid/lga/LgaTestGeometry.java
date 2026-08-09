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
 * The ONE test-scope definition of the LGA campaign's baseline member
 * radius (bead inviscid-0nx.30, E.3).
 *
 * <h2>Why this exists</h2>
 * {@code design-seeding-radius.md} §D-B makes {@code r} a PHYSICAL LEVER,
 * not a constant: {@code r} sets the contact-set measure, hence the
 * collision rate, hence the mean free path, hence the emergent viscosity.
 * The campaign's anchor is expected to MOVE off {@code 0.015} toward the
 * collision-dominated (inviscid) end of the sweep. Before this class, a
 * dozen test files each carried their own {@code private static final
 * double RADIUS = 0.015}, so retargeting the anchor was a twelve-file grep
 * whose blast radius was invisible at the point of edit. Consolidating
 * them here makes the test-scope half of that a one-line change.
 *
 * <h2>What "a one-line change" does and does not mean</h2>
 * Two qualifications, both load-bearing, because the phrase oversells
 * itself otherwise.
 * <ul>
 * <li>{@link #BASELINE_RADIUS} is an ALIAS of the PRODUCTION constant
 * {@code ContactAtlasGenerator.RADIUS}, not an independent definition. The
 * one line that moves the anchor is therefore an edit to production code,
 * and it moves more than these tests: it also retargets {@code
 * ContactAtlasGenerator.generate}'s 4-arg default overload, which every
 * caller that omits geometry parameters silently rides on. That coupling
 * is deliberate - it is what stops this class from claiming a radius the
 * generator does not actually generate at - but it is a production edit,
 * not a test-only one. Tests that reference {@code
 * ContactAtlasGenerator.RADIUS} DIRECTLY (e.g. {@code ContactComboCacheTest})
 * are dependents too, and do not appear among this constant's references.</li>
 * <li>{@code static final double} initialized from another {@code static
 * final double} is a compile-time {@code ConstantValue}: {@code javap}
 * shows consumers emitting {@code ldc2_w 0.015d} inline, with no read of
 * this field at runtime. The consolidation is a SOURCE-level one. A
 * partial recompilation that rebuilds this class but not its consumers
 * leaves the old value inlined in them; {@code mvn} does not do that, and
 * the project's regeneration discipline ({@code rm -rf target/classes
 * target/test-classes} before {@code mvn compile}) forecloses it, but an
 * IDE incremental build is not owed the same guarantee. If an anchor move
 * ever produces inexplicably mixed results, clean-build first.</li>
 * </ul>
 *
 * <h2>What deliberately does NOT use this constant</h2>
 * Two call sites keep LOCAL literals BY DOCUMENTED EXCEPTION, because for
 * them {@code 0.015} is not "the current anchor" but a FIXED PART OF A
 * PINNED CONTRACT that must not move when the anchor does:
 * <ul>
 * <li>{@code measure.SeamGoldenCompatTest} - pins Phase A golden numerics
 * at {@code r=0.015} by contract. Its {@code RADIUS} is read at exactly
 * ONE site ({@code newGeometry()}), feeding exactly one test, whose
 * radius-sensitive quantities are HARD LITERALS ({@code
 * LEDGER_TOTAL_COLLISIONS = 40L}, {@code LEDGER_EFFECTIVE_COLLISIONS =
 * 26L}). So consolidating it onto the anchor and then moving the anchor
 * makes that test go RED, not silently green - it stops comparing against
 * the golden it exists to defend, and says so. MEASURED, not assumed:
 * consolidating that file and retargeting the anchor to {@code 0.05}
 * produces {@code expected:<40> but was:<162>} (plus a second red in
 * {@code runOneSeed...}'s transport ratio). LOUD but ILLEGIBLE is the
 * actual failure mode, and it is why the exception is mechanized rather
 * than merely written down: nothing in either message points at the radius
 * consolidation as the cause. See
 * {@link #GOLDEN_PHASE_A_RADIUS}, which names the frozen value here so the
 * exception is discoverable from this side; that file itself is a PINNED
 * seam and is not edited to point back.</li>
 * <li>{@code ContactAtlasTest}'s header-string fixture - a verbatim
 * transcription of a v2 atlas header used to exercise the parser. Its
 * {@code 0.015} is TEXT under test, not a geometry parameter.</li>
 * </ul>
 * A third site, {@code animations.NecronomataVisualizationLengthsTest},
 * also keeps its local {@code 0.015f}: it is a RENDERING stroke radius on
 * the JavaFX visualization path (see {@code MemberGeometry}'s constructor
 * Javadoc, which notes the LGA radius merely "happens to" coincide with
 * what the visualization renders with), not an LGA collision parameter.
 * Moving the LGA anchor must not silently rescale the rendered struts, and
 * that test lives in {@code animations/}, outside the headless
 * {@code lga/} + {@code measure/} constraint that governs this constant's
 * consumers.
 *
 * @author halhildebrand
 */
public final class LgaTestGeometry {

    /**
     * The Phase A..C campaign's member radius, {@code 0.015} - sourced
     * from the production constant rather than restated, so this class
     * cannot drift from what {@link ContactAtlasGenerator} actually
     * generates atlases at.
     *
     * <p>
     * Per {@code design-seeding-radius.md} §D-B this value sits in the
     * BALLISTIC regime (~1.1e-4 contact fraction) - the opposite end of
     * the sweep from the epic's inviscid target - so it is a continuity
     * anchor, not a physically preferred operating point.
     */
    public static final double BASELINE_RADIUS = ContactAtlasGenerator.RADIUS;

    /**
     * The geometry LUT resolution every LGA fixture shares, {@code 360} -
     * likewise sourced, not restated. Public for the same reason as
     * {@link #BASELINE_RADIUS}: {@code measure}-package tests consume it
     * too.
     */
    public static final int BASELINE_GEOMETRY_RESOLUTION = ContactAtlasGenerator.GEOMETRY_RESOLUTION;

    /**
     * The FROZEN radius the Phase A golden numerics were captured at, as a
     * hard {@code 0.015} LITERAL - deliberately NOT {@link
     * #BASELINE_RADIUS}, and deliberately not an alias of anything that can
     * move.
     *
     * <p>
     * This constant exists to make an exception DISCOVERABLE that is
     * otherwise invisible. {@code measure.SeamGoldenCompatTest} declares
     * its own {@code private static final double RADIUS = 0.015} and
     * mentions this class nowhere, so a reader at that call site sees a
     * bare literal indistinguishable from the ones this bead consolidated,
     * and the natural tidy-up - point it at {@link #BASELINE_RADIUS} - is
     * the most damaging edit available in the suite.
     *
     * <p>
     * BE PRECISE ABOUT WHY. Not because it would pass silently: that file
     * reads {@code RADIUS} at exactly one site ({@code newGeometry()}) and
     * compares against HARD LITERALS ({@code LEDGER_TOTAL_COLLISIONS =
     * 40L}, {@code LEDGER_EFFECTIVE_COLLISIONS = 26L}), regenerating
     * nothing as a comparison input, so an anchor move makes it FAIL.
     * Measured by doing it - consolidate that file, retarget the anchor to
     * {@code 0.05}, CLEAN-build, run: test 8 reports {@code expected:<40>
     * but was:<162>}, and test 9's transport ratio reds alongside it - the
     * latter for a DIFFERENT reason, see the half-protection note below.
     * (Without the clean build the ledger test stays GREEN, because {@link
     * #BASELINE_RADIUS} is a {@code ConstantValue} still inlined in the
     * unrecompiled consumer - a live demonstration of the inlining hazard
     * this class's second qualification describes, and worth knowing
     * before trusting any single anchor-move experiment.)
     *
     * <p>
     * HALF-PROTECTION, STATED PLAINLY. This exception governs only the
     * pins that actually READ that file's local literal, which is not all
     * of them. {@code RADIUS} is read at exactly one site, {@code
     * newGeometry()}, feeding test 8's ledger pins ({@code 40L} / {@code
     * 26L}) - those ARE protected by keeping the literal. Test 9's pins
     * ({@code TRANSPORT_RATIO_PIN}, {@code
     * RUN_ONE_SEED_TOTAL_COLLISIONS = 246L}, {@code
     * RUN_ONE_SEED_EFFECTIVE_COLLISIONS = 6L}) run through {@code
     * AnisotropyProbe.phaseAHybridSubstrate}, which builds its {@code
     * MemberGeometry} from {@link ContactAtlasGenerator#RADIUS} DIRECTLY,
     * so they ride the production anchor whatever the seam file's literal
     * says. It follows that test 9's transport-ratio red in the experiment
     * above was caused by the ANCHOR MOVE ALONE and not by the
     * consolidation: it would have arrived with that file left entirely
     * untouched. Keeping the local literal is still correct - it is just
     * not sufficient to pin Phase A. Closing the other half is bead {@code
     * inviscid-yj3}.
     *
     * <p>
     * The damage is that the failure is UNATTRIBUTABLE - a collision-count
     * mismatch inside a seam-compat test, arriving alongside whatever else
     * the anchor move broke, with nothing naming the radius consolidation
     * as the cause - and that the pinned Phase A golden has by then been
     * silently retargeted whether or not anybody reads the failure
     * correctly. {@code
     * LgaTestGeometryTest.phaseAGoldenSeamKeepsItsOwnFrozenRadius} turns
     * that into a named, self-explaining red at the moment of the edit,
     * rather than an obscure one at the moment of the anchor move.
     *
     * <p>
     * CONTRAST, so the rule is legible rather than arbitrary: {@code
     * measure.AuditedRunTest} WAS correctly consolidated onto {@link
     * #BASELINE_RADIUS}. Its collision assertions are BANDS ({@code > 0}),
     * not pinned counts, so it survives an anchor move by design. What
     * decides the exception is whether the assertions are frozen numerics,
     * not whether the file is about collisions.
     *
     * <p>
     * NOT for use as a geometry parameter by anything else. Test code that
     * wants "the current campaign anchor" wants {@link #BASELINE_RADIUS};
     * this is "the value Phase A's goldens are frozen at", and the two are
     * equal today only by coincidence of the campaign not having moved yet.
     */
    public static final double GOLDEN_PHASE_A_RADIUS = 0.015;

    private LgaTestGeometry() {
    }
}
