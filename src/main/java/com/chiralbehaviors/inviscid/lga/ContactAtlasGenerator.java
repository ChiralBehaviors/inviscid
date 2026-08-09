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

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.Necronomata;
import com.chiralbehaviors.inviscid.PhiCoordinates;
import com.chiralbehaviors.inviscid.measure.AuditedRun;
import com.chiralbehaviors.inviscid.measure.CollisionStatistics;
import com.chiralbehaviors.inviscid.measure.ConservationAudit;

/**
 * Generates a {@link ContactAtlas} (bead inviscid-0nx.16, A.5): reproducible,
 * parameterized by {@code nLga} (the AUTOMATON phase resolution N_lga - a
 * USER-RESERVED decision from {@code {8, 12, 16, 24}}, made downstream of
 * this class from the measurement campaign it feeds; this class never picks
 * a final value itself).
 *
 * <h2>N_lga=24 (the committed artifact's value, USER DECISION 2026-08-08)</h2>
 * Recorded verbatim on bead inviscid-0nx.16's {@code --design} field, full
 * data in T2 {@code inviscid/analysis-nlga-candidates.md} (post
 * inviscid-0nx.16.1 correction). Summary: over the exhaustive 446-combo
 * campaign (native 360-step discovery, the provable ceiling), transcription
 * error was lowest at 24 (0.19%, 2x better than 16); 12 was empirically
 * DOMINATED by 8 (0.49% vs 0.42% - a genuine bin-boundary aliasing effect,
 * not sampling noise); the 20.18% floor risk (90/446 combos with a
 * near-point contact region narrower than even the coarsest 45-degree bin)
 * is N_lga-INVARIANT across the whole {8,12,16,24} candidate range - a
 * radius/geometry property no candidate resolution fixes, so it did not
 * discriminate between them; and every candidate's cost (rows, generation
 * wall time) was trivial. The one real cost of choosing 24 - dynamic
 * coverage dilution at finer bins - was addressed by scaling up {@code
 * ticksObserved} for the committed generation run (see {@link
 * CommittedContactAtlasTest} for the achieved-vs-targeted coverage result,
 * which came in well below the crude Poisson-saturation extrapolation - an
 * honestly-reported miss, not a silently-accepted one).
 *
 * <h2>Three data sources, merged (v2, bead inviscid-gyt)</h2>
 * <ol>
 * <li><b>Overlap + bin-center ground truth</b> ({@link
 * #sweepOverlapAndCenter}) - for every combo in the exhaustive ever-
 * contacting universe ({@link ContactComboCache}), a fine {@code
 * geometryResolution x geometryResolution} angle sweep is aggregated into
 * each {@code nLga x nLga} bin CELL (not just its center): {@code
 * overlapFraction = (fine samples in this cell that contact) / (fine
 * samples in this cell)} - the ANY-OVERLAP transcription signal (USER
 * DECISION 2026-08-08, bead inviscid-gyt): a cell is "fired" iff {@code
 * overlapFraction > 0}. {@code contact} (the original v1 signal) is
 * retained, computed at the bin center only, for comparability - see
 * {@link #sweepOverlapAndCenter}'s own Javadoc for the proof that {@code
 * contact=true} always implies {@code overlapFraction > 0}.</li>
 * <li><b>Dynamic reachability</b> ({@link #runDynamicReachability}) - drives
 * a real Phase A {@link AuditedRun} for {@code ticksObserved} ticks,
 * snapshotting {@code angle} PRE-TICK (the .15 planning note recorded on
 * this bead: no ready-made seam exists, so the harness reads the full
 * {@code angle} array via {@link Necronomata#process(Necronomata.Processor)}
 * immediately before each {@link AuditedRun#tick(int)} call) and correlating
 * {@code TickResult.applied()}'s resolved contacts against those frozen
 * pre-tick angles, quantized to {@code nLga} bins. {@link ContactScan}
 * (the source of every dynamic observation) only ever scans {@link
 * FccNeighborhood}'s 6 CANONICAL (positive) directions - so a dynamic
 * observation's {@code direction} is always positive; see "Negative-
 * direction observedCount mirroring" below for how the other 6 directions
 * get their {@code observedCount}.</li>
 * <li><b>Negative-direction observedCount mirroring</b> ({@link
 * #mirrorNegativeDirectionObservedCounts}, bead inviscid-gyt Phase A gate
 * finding): because {@link ContactScan} canonicalizes to the 6 positive
 * directions, a row keyed by a NEGATIVE direction can never receive a
 * dynamic observation directly - its {@code observedCount} would always
 * read 0, silently understating real dynamic contact rates for exactly
 * half of {@link FccNeighborhood#DIRECTIONS}. Fixed by mirroring AT
 * GENERATION TIME (option (a) of the bead's two choices - keeps the
 * artifact self-contained, so a Phase C consumer reading the atlas alone
 * sees every direction's real reachability without also having to
 * replicate this mirror logic itself): for every row with a positive
 * {@code direction} and {@code observedCount > 0}, the row keyed by
 * {@code (oppositeDirection, cubeB, memberB, cubeA, memberA, phaseBinB,
 * phaseBinA)} - the same physical contact, described from the other
 * cell's side - has its {@code observedCount} SET (not added) to match.
 * This is sound, not merely convenient: {@link ContactPredicate#minDistance}
 * is provably symmetric under this exact transform (translation-invariance
 * of Euclidean distance - {@code distance(A, B+offset) == distance(B,
 * A-offset)}), so the mirror row's {@code overlapFraction} independently
 * computes to the identical value via {@link #sweepOverlapAndCenter} alone
 * (which already sweeps all 12 directions, not just the canonical 6) -
 * mirroring only needs to carry {@code observedCount}, which {@code
 * sweepOverlapAndCenter} cannot supply.</li>
 * </ol>
 * Merge key: {@code (direction, cubeA, memberA, cubeB, memberB, phaseBinA,
 * phaseBinB)}.
 *
 * <h2>{@code gitCommit} provenance (bead inviscid-0nx.16.2)</h2>
 * {@link #resolveGitCommit()} runs {@code git rev-parse HEAD} against the
 * process's working directory at generation time - never a fixed literal
 * (a fixed literal makes {@code ContactAtlas.readValidated}'s {@code
 * gitCommit} check vacuous: it would always match). A dirty working tree
 * (per {@code git status --porcelain}) appends the {@code "-dirty"}
 * suffix; an environment without a {@code git} executable, or a directory
 * that is not a git repository, falls back to the literal {@code
 * "UNKNOWN"} rather than throwing - atlas generation must not hard-depend
 * on git being present. The 7-arg {@link #generate(int, Point3i, long,
 * int, int, double, String)} overload accepts an explicit {@code
 * gitCommit} instead, for callers (tests, reproducibility harnesses) that
 * need a fixed value independent of the invoking environment's git state.
 *
 * @author halhildebrand
 */
public final class ContactAtlasGenerator {

    public static final double  RADIUS               = 0.015;
    /** Matches the resolution {@code ContactPredicateTest}/{@code ContactScanTest} use for {@link MemberGeometry}'s internal LUT. */
    public static final int     GEOMETRY_RESOLUTION   = 360;
    public static final long    DEFAULT_SEED          = 42L;
    public static final Point3i DEFAULT_EXTENT        = new Point3i(4, 4, 4);
    public static final int     DEFAULT_TICKS         = 2000;
    /** Matches {@code AuditedRunTest}'s seeded-quanta bound convention. */
    public static final int     QUANTA_BOUND          = 6;

    private static final int CUBES_PER_CELL   = 5;
    private static final int MEMBERS_PER_CUBE = 6;
    private static final double TWO_PI        = 2 * Math.PI;

    private ContactAtlasGenerator() {
    }

    /**
     * @return the bin-center angle (radians) of bin {@code bin} out of
     *         {@code nLga} equal bins spanning {@code [0, 2*pi)}.
     */
    static double binCenter(int bin, int nLga) {
        return (bin + 0.5) * (TWO_PI / nLga);
    }

    /**
     * @return the {@code MemberGeometry}-quantized LUT step (in {@code [0,
     *         geometryResolution)}) that {@code angle} falls into - an
     *         EXACT replica of {@code MemberGeometry.stepOf}'s private
     *         arithmetic ({@code Constants.TWO_PI}, float precision), so a
     *         dynamically-observed contact's step is derived identically
     *         to the step {@link ContactPredicate} actually evaluated
     *         geometry at (see {@link #binOfStep} and {@link
     *         ContactComboCache#angleOf}'s Javadoc for why this exact
     *         replication, not the more natural {@code
     *         ContactAtlasGenerator}-local {@code TWO_PI} (double), is
     *         required for {@code recordObservedContact}'s bin to align
     *         with {@link #sweepOverlapAndCenter}'s fine-grid bins - bead
     *         inviscid-gyt Phase A gate rework).
     */
    static int stepOf(float angle, int geometryResolution) {
        float twoPi = com.chiralbehaviors.inviscid.Constants.TWO_PI;
        float normalized = angle % twoPi;
        double widened = normalized;
        if (widened < 0) {
            widened += twoPi;
        }
        double angularResolution = twoPi / geometryResolution;
        return ((int) (widened / angularResolution)) % geometryResolution;
    }

    /**
     * @return the {@code nLga} bin index a {@code geometryResolution}-step
     *         LUT {@code step} falls into - PURE INTEGER arithmetic
     *         ({@code step * nLga / geometryResolution}, Java {@code int}
     *         division truncates toward zero, equivalent to {@code floor}
     *         for non-negative operands), deliberately never routed
     *         through a reconstructed floating-point angle: that round
     *         trip is exactly what {@link ContactComboCache#angleOf}'s
     *         Javadoc documents as unreliable. Well-defined even when
     *         {@code geometryResolution} is not an exact multiple of
     *         {@code nLga} (e.g. {@code nLga=16}, {@code
     *         geometryResolution=360}): the resulting bins are merely
     *         unequal-width in step-count terms, which {@link
     *         #sweepOverlapAndCenter}'s {@code countPerBin}-based
     *         denominator already accounts for.
     */
    static int binOfStep(int step, int nLga, int geometryResolution) {
        return step * nLga / geometryResolution;
    }

    private record RowKey(int direction, int cubeA, int memberA, int cubeB,
                           int memberB, int phaseBinA, int phaseBinB) {
    }

    private static final class MutableRow {
        final RowKey key;
        boolean      contact;
        double       overlapFraction;
        double       minDistance;
        long         observedCount;

        MutableRow(RowKey key, boolean contact, double minDistance) {
            this.key = key;
            this.contact = contact;
            this.minDistance = minDistance;
        }
    }

    /**
     * Generates a {@link ContactAtlas} at the default geometry parameters
     * ({@link #GEOMETRY_RESOLUTION}, {@link #RADIUS}), resolving {@code
     * gitCommit} for real (see class Javadoc).
     */
    public static ContactAtlas generate(int nLga, Point3i extent, long seed,
                                         int ticksObserved) {
        return generate(nLga, extent, seed, ticksObserved,
                         GEOMETRY_RESOLUTION, RADIUS, resolveGitCommit());
    }

    /**
     * Resolves {@code gitCommit} for real (see class Javadoc) rather than
     * accepting a caller-supplied one - use the 7-arg overload for a fixed
     * value.
     */
    public static ContactAtlas generate(int nLga, Point3i extent, long seed,
                                         int ticksObserved,
                                         int geometryResolution,
                                         double memberRadius) {
        return generate(nLga, extent, seed, ticksObserved, geometryResolution,
                         memberRadius, resolveGitCommit());
    }

    public static ContactAtlas generate(int nLga, Point3i extent, long seed,
                                         int ticksObserved,
                                         int geometryResolution,
                                         double memberRadius,
                                         String gitCommit) {
        if (nLga <= 0) {
            throw new IllegalArgumentException("nLga must be positive: "
                                                + nLga);
        }
        MemberGeometry geometry = new MemberGeometry(geometryResolution,
                                                      memberRadius);
        ContactPredicate predicate = new ContactPredicate(geometry);

        Map<RowKey, MutableRow> rows = new LinkedHashMap<>();
        sweepOverlapAndCenter(predicate, nLga, geometryResolution,
                               memberRadius, rows);
        runDynamicReachability(predicate, nLga, geometryResolution, extent,
                                seed, ticksObserved, rows);
        mirrorNegativeDirectionObservedCounts(predicate, nLga, rows);

        List<ContactAtlas.Row> atlasRows = new ArrayList<>(rows.size());
        for (MutableRow row : rows.values()) {
            atlasRows.add(new ContactAtlas.Row(row.key.direction(),
                                                row.key.cubeA(),
                                                row.key.memberA(),
                                                row.key.cubeB(),
                                                row.key.memberB(),
                                                row.key.phaseBinA(),
                                                row.key.phaseBinB(),
                                                row.contact,
                                                row.overlapFraction,
                                                row.minDistance,
                                                row.observedCount));
        }

        ContactAtlas.Header header = new ContactAtlas.Header(ContactAtlas.ATLAS_VERSION,
                                                              ContactAtlasGenerator.class.getName(),
                                                              gitCommit,
                                                              memberRadius,
                                                              geometryResolution,
                                                              PhiCoordinates.Cubes[0].getEdgeLength(),
                                                              nLga, "Cubes[0]",
                                                              new Point3i(extent),
                                                              seed,
                                                              ticksObserved);
        return new ContactAtlas(header, atlasRows);
    }

    /**
     * Runs {@code git rev-parse HEAD} against the process's working
     * directory, appending {@code "-dirty"} if {@code git status
     * --porcelain} reports uncommitted changes. Falls back to the literal
     * {@code "UNKNOWN"} (never throws) if {@code git} is unavailable, the
     * working directory is not a git repository, or either command exits
     * non-zero - see class Javadoc.
     */
    static String resolveGitCommit() {
        String sha = runGit("rev-parse", "HEAD");
        if (sha == null || sha.isBlank()) {
            return "UNKNOWN";
        }
        return isDirty() ? sha + "-dirty" : sha;
    }

    private static boolean isDirty() {
        String status = runGit("status", "--porcelain");
        return status != null && !status.isBlank();
    }

    private static String runGit(String... args) {
        try {
            List<String> command = new ArrayList<>();
            command.add("git");
            command.addAll(List.of(args));
            Process process = new ProcessBuilder(command).redirectErrorStream(true)
                                                           .start();
            String output = new String(process.getInputStream().readAllBytes(),
                                        StandardCharsets.UTF_8).trim();
            int exit = process.waitFor();
            return exit == 0 ? output : null;
        } catch (IOException e) {
            return null;
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return null;
        }
    }

    /**
     * The v2 overlap sweep (bead inviscid-gyt): for every combo in the
     * exhaustive ever-contacting universe ({@link
     * ContactComboCache#combosFor}, {@code ~446} at this class's default
     * geometry - see that class's Javadoc for why the search space is
     * restricted to that provably-complete set rather than the full
     * {@code 12 * 30 * 30 == 10,800}), one fine {@code geometryResolution
     * x geometryResolution} {@link ContactPredicate#contacts} sweep is
     * aggregated per {@code nLga x nLga} bin cell into {@code
     * overlapFraction = hits / (fineSamplesPerBinA * fineSamplesPerBinB)}.
     * A row is created for every cell with {@code overlapFraction > 0}.
     *
     * <h2>Proof: {@code contact=true} (bin center) implies {@code
     * overlapFraction > 0}</h2>
     * {@link MemberGeometry} quantizes any continuous angle to {@code step
     * = floor(normalize(angle) / (2*pi/geometryResolution))}. The bin
     * center angle is {@code (bin+0.5) * (2*pi/nLga)}, so its step is
     * {@code floor((bin+0.5) * geometryResolution/nLga)} - algebraically
     * IDENTICAL to the fine sweep's own step-index-to-bin floor mapping
     * (see {@link #binOf}) evaluated at that same fine index. In other
     * words, the bin center's {@link ContactPredicate} evaluation and one
     * specific fine-grid sample (the one at that same LUT step) are
     * bit-for-bit the same computation. So whenever {@code contact=true}
     * at a bin center, that specific fine sample is counted as a hit in
     * this sweep's histogram for that same cell, guaranteeing {@code
     * overlapFraction >= 1/(fineSamplesPerBinA*fineSamplesPerBinB) > 0}.
     * {@link ContactAtlasTest} and {@link CommittedContactAtlasTest} both
     * assert this invariant directly over generated data, not just by
     * this proof.
     */
    private static void sweepOverlapAndCenter(ContactPredicate predicate,
                                                int nLga,
                                                int geometryResolution,
                                                double memberRadius,
                                                Map<RowKey, MutableRow> rows) {
        int[] fineBinOf = new int[geometryResolution];
        int[] countPerBin = new int[nLga];
        for (int step = 0; step < geometryResolution; step++) {
            int bin = binOfStep(step, nLga, geometryResolution);
            fineBinOf[step] = bin;
            countPerBin[bin]++;
        }

        List<ContactComboCache.Combo> combos = ContactComboCache.combosFor(predicate,
                                                                             geometryResolution,
                                                                             memberRadius);
        for (ContactComboCache.Combo combo : combos) {
            sweepComboOverlap(predicate, nLga, geometryResolution, combo,
                               fineBinOf, countPerBin, rows);
        }
    }

    private static void sweepComboOverlap(ContactPredicate predicate, int nLga,
                                           int geometryResolution,
                                           ContactComboCache.Combo combo,
                                           int[] fineBinOf, int[] countPerBin,
                                           Map<RowKey, MutableRow> rows) {
        long[][] contactCount = new long[nLga][nLga];
        for (int a = 0; a < geometryResolution; a++) {
            float angleA = ContactComboCache.angleOf(a, geometryResolution);
            int binA = fineBinOf[a];
            for (int b = 0; b < geometryResolution; b++) {
                float angleB = ContactComboCache.angleOf(b, geometryResolution);
                if (predicate.contacts(combo.cubeA(), combo.memberA(), angleA,
                                       combo.cubeB(), combo.memberB(), angleB,
                                       combo.direction())) {
                    contactCount[binA][fineBinOf[b]]++;
                }
            }
        }

        for (int binA = 0; binA < nLga; binA++) {
            for (int binB = 0; binB < nLga; binB++) {
                long hits = contactCount[binA][binB];
                if (hits == 0) {
                    continue;
                }
                double denominator = (double) countPerBin[binA]
                                     * countPerBin[binB];
                double overlapFraction = hits / denominator;

                float centerA = (float) binCenter(binA, nLga);
                float centerB = (float) binCenter(binB, nLga);
                boolean centerContact = predicate.contacts(combo.cubeA(),
                                                            combo.memberA(),
                                                            centerA,
                                                            combo.cubeB(),
                                                            combo.memberB(),
                                                            centerB,
                                                            combo.direction());
                double centerMinDistance = predicate.minDistance(combo.cubeA(),
                                                                  combo.memberA(),
                                                                  centerA,
                                                                  combo.cubeB(),
                                                                  combo.memberB(),
                                                                  centerB,
                                                                  combo.direction());

                RowKey key = new RowKey(combo.direction(), combo.cubeA(),
                                        combo.memberA(), combo.cubeB(),
                                        combo.memberB(), binA, binB);
                MutableRow row = new MutableRow(key, centerContact,
                                                 centerMinDistance);
                row.overlapFraction = overlapFraction;
                rows.put(key, row);
            }
        }
    }

    /**
     * Fixes the negative-direction {@code observedCount} asymmetry (Phase
     * A gate finding, bead inviscid-gyt) - see class Javadoc "Negative-
     * direction observedCount mirroring". Iterates a SNAPSHOT of {@code
     * rows}'s positive-direction, dynamically-observed entries (not the
     * live map) because this mutates {@code rows} in place while walking
     * it.
     */
    private static void mirrorNegativeDirectionObservedCounts(ContactPredicate predicate,
                                                                int nLga,
                                                                Map<RowKey, MutableRow> rows) {
        List<MutableRow> positivelyObserved = rows.values().stream()
                                                    .filter(row -> row.key.direction() > 0
                                                                   && row.observedCount > 0)
                                                    .toList();
        for (MutableRow source : positivelyObserved) {
            RowKey key = source.key;
            RowKey mirrorKey = new RowKey(FccNeighborhood.opposite(key.direction()),
                                          key.cubeB(), key.memberB(),
                                          key.cubeA(), key.memberA(),
                                          key.phaseBinB(), key.phaseBinA());
            MutableRow mirror = rows.get(mirrorKey);
            if (mirror == null) {
                // Defensive fallback only - see class Javadoc's symmetry
                // proof for why sweepOverlapAndCenter (which sweeps all 12
                // directions) is expected to have already created this row
                // whenever the source row's overlapFraction is positive.
                float angleA = (float) binCenter(mirrorKey.phaseBinA(), nLga);
                float angleB = (float) binCenter(mirrorKey.phaseBinB(), nLga);
                boolean centerContact = predicate.contacts(mirrorKey.cubeA(),
                                                            mirrorKey.memberA(),
                                                            angleA,
                                                            mirrorKey.cubeB(),
                                                            mirrorKey.memberB(),
                                                            angleB,
                                                            mirrorKey.direction());
                double minDistance = predicate.minDistance(mirrorKey.cubeA(),
                                                            mirrorKey.memberA(),
                                                            angleA,
                                                            mirrorKey.cubeB(),
                                                            mirrorKey.memberB(),
                                                            angleB,
                                                            mirrorKey.direction());
                mirror = new MutableRow(mirrorKey, centerContact, minDistance);
                rows.put(mirrorKey, mirror);
            }
            mirror.observedCount = source.observedCount;
        }
    }

    /**
     * Drives a real Phase A {@link AuditedRun} for {@code ticksObserved}
     * ticks, snapshotting {@code angle} pre-tick and correlating each
     * tick's resolved contacts against it - see class Javadoc.
     */
    private static void runDynamicReachability(ContactPredicate predicate,
                                                 int nLga,
                                                 int geometryResolution,
                                                 Point3i extent,
                                                 long seed, int ticksObserved,
                                                 Map<RowKey, MutableRow> rows) {
        Necronomata automaton = new Necronomata(extent);
        seedRandomAngles(automaton, extent, seed);
        seedRandomQuanta(automaton, extent, seed, QUANTA_BOUND);

        FccNeighborhood neighborhood = new FccNeighborhood(automaton);
        ContactScan scan = new ContactScan(automaton, neighborhood, predicate);
        CollisionStatistics statistics = new CollisionStatistics();
        CollisionSweep sweep = new CollisionSweep(automaton, scan,
                                                   new QuantaExchangeRule(),
                                                   statistics);
        HybridAutomaton hybrid = new HybridAutomaton(automaton, sweep);
        ConservationAudit audit = new ConservationAudit(automaton);
        AuditedRun run = new AuditedRun(hybrid, audit);

        int length = 30 * extent.x * extent.y * extent.z;
        float[] preTickAngles = new float[length];

        for (int tick = 0; tick < ticksObserved; tick++) {
            // .15 planning note: no ready-made seam retains the angles a
            // contact was resolved against - snapshot PRE-TICK via the
            // read-only process(Processor) escape hatch, before this
            // tick's HybridAutomaton.tick() (inside AuditedRun.tick())
            // advances angle.
            automaton.process((angle, frequency, deltaA, deltaF) -> System.arraycopy(angle,
                                                                                       0,
                                                                                       preTickAngles,
                                                                                       0,
                                                                                       length));

            AuditedRun.TickOutcome outcome = run.tick(tick);

            // (bead inviscid-ckn / inviscid-0nx.21) outcome.collisionResult()
            // is TickReport-typed since the seam widened AuditedRun to any
            // TickDriver. Atlas generation correlates geometric contacts
            // against pre-tick angles -- an operation with no meaning for a
            // table-driven LGA, which CONSUMES the atlas rather than
            // producing it -- so this narrowing is honest, not a
            // workaround: a non-hybrid driver here is a caller error.
            if (!(outcome.collisionResult() instanceof CollisionSweep.TickResult r)) {
                throw new IllegalStateException("the contact atlas is generated from Phase A geometric contacts; "
                                                 + "driver reported "
                                                 + outcome.collisionResult()
                                                           .getClass()
                                                           .getName());
            }
            for (CollisionSweep.AppliedCollision applied : r.applied()) {
                recordObservedContact(automaton, predicate, nLga,
                                       geometryResolution, preTickAngles,
                                       applied.contact(), rows);
            }
        }
    }

    /**
     * Bins a dynamically-observed contact via {@link #stepOf} + {@link
     * #binOfStep} - the SAME two-step, pure-integer-final-stage pipeline
     * {@link #sweepOverlapAndCenter} uses to bin its fine-grid samples -
     * rather than binning the raw continuous angle directly. This is
     * load-bearing, not stylistic (bead inviscid-gyt Phase A gate rework):
     * a continuous-angle bin can legitimately disagree with the bin the
     * fine sweep attributes to the SAME angle's LUT step (measured
     * ~2.8% mismatch rate on uniform random angles during development),
     * which would silently make some real dynamic contacts land in a
     * {@code (binA, binB)} cell {@link #sweepOverlapAndCenter} never
     * marked as fired - exactly the "falsifies the ribbon explanation"
     * anomaly {@code CommittedContactAtlasTest
     * .everyDynamicallyObservedCellHasPositiveOverlapFraction} checks for.
     * Binning via the same {@code step -> binOfStep} formula both paths
     * share eliminates the mismatch by construction, not by narrowing
     * floating-point tolerance.
     */
    private static void recordObservedContact(Necronomata automaton,
                                                ContactPredicate predicate,
                                                int nLga,
                                                int geometryResolution,
                                                float[] preTickAngles,
                                                Contact contact,
                                                Map<RowKey, MutableRow> rows) {
        int indexA = automaton.indexOfCell(contact.cellA())
                     + contact.cubeA() * MEMBERS_PER_CUBE + contact.memberA();
        int indexB = automaton.indexOfCell(contact.cellB())
                     + contact.cubeB() * MEMBERS_PER_CUBE + contact.memberB();
        int stepA = stepOf(preTickAngles[indexA], geometryResolution);
        int stepB = stepOf(preTickAngles[indexB], geometryResolution);
        int binA = binOfStep(stepA, nLga, geometryResolution);
        int binB = binOfStep(stepB, nLga, geometryResolution);

        RowKey key = new RowKey(contact.direction(), contact.cubeA(),
                                contact.memberA(), contact.cubeB(),
                                contact.memberB(), binA, binB);
        MutableRow row = rows.get(key);
        if (row == null) {
            float angleA = (float) binCenter(binA, nLga);
            float angleB = (float) binCenter(binB, nLga);
            boolean binCenterContact = predicate.contacts(contact.cubeA(),
                                                           contact.memberA(),
                                                           angleA,
                                                           contact.cubeB(),
                                                           contact.memberB(),
                                                           angleB,
                                                           contact.direction());
            double minDistance = predicate.minDistance(contact.cubeA(),
                                                        contact.memberA(),
                                                        angleA,
                                                        contact.cubeB(),
                                                        contact.memberB(),
                                                        angleB,
                                                        contact.direction());
            row = new MutableRow(key, binCenterContact, minDistance);
            rows.put(key, row);
        }
        row.observedCount++;
    }

    private static void seedRandomAngles(Necronomata automaton,
                                          Point3i extent, long seed) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        float[] angles = new float[length];
        for (int i = 0; i < length; i++) {
            angles[i] = random.nextFloat() * (float) TWO_PI;
        }
        automaton.process((angleArray, frequency, deltaA, deltaF) -> System.arraycopy(angles,
                                                                                        0,
                                                                                        angleArray,
                                                                                        0,
                                                                                        length));
    }

    private static void seedRandomQuanta(Necronomata automaton,
                                          Point3i extent, long seed,
                                          int bound) {
        Random random = new Random(seed);
        int length = 30 * extent.x * extent.y * extent.z;
        float[] quanta = new float[length];
        for (int i = 0; i < length; i++) {
            quanta[i] = random.nextInt(bound);
        }
        automaton.process((angleArray, frequency, deltaA, deltaF) -> System.arraycopy(quanta,
                                                                                        0,
                                                                                        frequency,
                                                                                        0,
                                                                                        length));
    }

    /**
     * Regenerates an atlas, reproducible from the header parameters alone
     * (bead inviscid-0nx.16's acceptance criterion): {@code args[0]} is
     * {@code nLga} (default 12, a provisional candidate value - NOT the
     * epic's user-reserved decision unless explicitly passed), {@code
     * args[1]} is {@code ticksObserved} (default {@link #DEFAULT_TICKS}),
     * {@code args[2]} is the output path (default {@code
     * target/contact-atlas-candidate-n<nLga>.tsv} - a caller writing the
     * FINAL, user-decided atlas must pass {@code
     * src/test/resources/lga/contact-atlas-v2.tsv} explicitly; this
     * method never defaults there itself, so an accidental bare
     * invocation cannot silently overwrite the committed artifact). Run
     * manually (IDE or classpath invocation - no exec plugin is
     * configured in this project).
     */
    public static void main(String[] args) throws IOException {
        int nLga = args.length > 0 ? Integer.parseInt(args[0]) : 12;
        int ticks = args.length > 1 ? Integer.parseInt(args[1])
                                     : DEFAULT_TICKS;
        Point3i extent = DEFAULT_EXTENT;
        long seed = DEFAULT_SEED;
        Path out = args.length > 2 ? Path.of(args[2])
                                    : Path.of("target",
                                              "contact-atlas-candidate-n"
                                              + nLga + ".tsv");

        long start = System.nanoTime();
        ContactAtlas atlas = generate(nLga, extent, seed, ticks);
        long elapsedMs = (System.nanoTime() - start) / 1_000_000;

        long contactRows = 0;
        long observedRows = 0;
        long overlapPositiveRows = 0;
        long anomalousRows = 0;
        double overlapSum = 0.0;
        double overlapMin = Double.POSITIVE_INFINITY;
        double overlapMax = Double.NEGATIVE_INFINITY;
        List<Double> nonZeroOverlaps = new ArrayList<>();
        for (ContactAtlas.Row row : atlas.rows()) {
            if (row.contact()) {
                contactRows++;
            }
            if (row.observedCount() > 0) {
                observedRows++;
                if (row.overlapFraction() == 0.0) {
                    anomalousRows++;
                }
            }
            if (row.overlapFraction() > 0.0) {
                overlapPositiveRows++;
                overlapSum += row.overlapFraction();
                overlapMin = Math.min(overlapMin, row.overlapFraction());
                overlapMax = Math.max(overlapMax, row.overlapFraction());
                nonZeroOverlaps.add(row.overlapFraction());
            }
        }
        nonZeroOverlaps.sort(Double::compareTo);
        double overlapMedian = nonZeroOverlaps.isEmpty() ? 0.0
                                                          : nonZeroOverlaps.get(nonZeroOverlaps.size()
                                                                                 / 2);
        double overlapMean = overlapPositiveRows == 0 ? 0.0
                                                       : overlapSum
                                                         / overlapPositiveRows;
        double overFireRatio = overlapSum == 0.0 ? 0.0
                                                  : overlapPositiveRows
                                                    / overlapSum;

        atlas.write(out);
        System.out.println("nLga=" + nLga + " ticksObserved=" + ticks
                            + " rows=" + atlas.rows().size() + " (contact="
                            + contactRows + ", no-contact="
                            + (atlas.rows().size() - contactRows)
                            + ") (observed>0=" + observedRows
                            + ", observed=0=" + (atlas.rows().size()
                                                  - observedRows)
                            + ") elapsedMs=" + elapsedMs + " gitCommit="
                            + atlas.header().gitCommit() + " -> " + out);
        System.out.println("overlapFraction: firedCells=" + overlapPositiveRows
                            + " min=" + overlapMin + " median=" + overlapMedian
                            + " mean=" + overlapMean + " max=" + overlapMax
                            + " overFireRatio(firedCells/sumOverlapFraction)="
                            + overFireRatio + " anomalousObservedZeroOverlap="
                            + anomalousRows);
    }
}
