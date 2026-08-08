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
 * <h2>Two data sources, merged</h2>
 * <ol>
 * <li><b>Geometric ground truth</b> ({@link #sweepGeometricGroundTruth}) -
 * {@link ContactPredicate} swept over the full bin grid: all 12 {@link
 * FccNeighborhood#DIRECTIONS}, all {@code 30x30} {@code (cubeA, memberA,
 * cubeB, memberB)} combinations, all {@code nLga^2} bin pairs, evaluated at
 * bin centers. Every combination the sweep finds contacting becomes a row
 * (bead: "one row per contacting combination").</li>
 * <li><b>Dynamic reachability</b> ({@link #runDynamicReachability}) - drives
 * a real Phase A {@link AuditedRun} for {@code ticksObserved} ticks,
 * snapshotting {@code angle} PRE-TICK (the .15 planning note recorded on
 * this bead: no ready-made seam exists, so the harness reads the full
 * {@code angle} array via {@link Necronomata#process(Necronomata.Processor)}
 * immediately before each {@link AuditedRun#tick(int)} call) and correlating
 * {@code TickResult.applied()}'s resolved contacts against those frozen
 * pre-tick angles, quantized to {@code nLga} bins. A dynamically-observed
 * combination not already present from the geometric sweep (its bin-center
 * verdict can legitimately disagree with the fine-grained angle the dynamic
 * run actually contacted at - exactly the quantization-fidelity risk the
 * bead's measurement campaign quantifies) is still added, with its {@code
 * contact} field recomputed at ITS bin centers, so no dynamically-observed
 * event is ever silently dropped.</li>
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
     * @return the bin index (in {@code [0, nLga)}) that {@code angle}
     *         (any real value, wrapped) falls into.
     */
    static int binOf(float angle, int nLga) {
        double normalized = ((angle % TWO_PI) + TWO_PI) % TWO_PI;
        int bin = (int) Math.floor(normalized / (TWO_PI / nLga));
        return Math.min(bin, nLga - 1);
    }

    private record RowKey(int direction, int cubeA, int memberA, int cubeB,
                           int memberB, int phaseBinA, int phaseBinB) {
    }

    private static final class MutableRow {
        final RowKey key;
        boolean      contact;
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
        sweepGeometricGroundTruth(predicate, nLga, rows);
        runDynamicReachability(predicate, nLga, extent, seed, ticksObserved,
                                rows);

        List<ContactAtlas.Row> atlasRows = new ArrayList<>(rows.size());
        for (MutableRow row : rows.values()) {
            atlasRows.add(new ContactAtlas.Row(row.key.direction(),
                                                row.key.cubeA(),
                                                row.key.memberA(),
                                                row.key.cubeB(),
                                                row.key.memberB(),
                                                row.key.phaseBinA(),
                                                row.key.phaseBinB(),
                                                row.contact, row.minDistance,
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
     * The full geometric sweep: {@code 12 * 30 * 30 * nLga^2} {@link
     * ContactPredicate#contacts} evaluations at bin centers. Only
     * combinations found contacting are added.
     */
    private static void sweepGeometricGroundTruth(ContactPredicate predicate,
                                                    int nLga,
                                                    Map<RowKey, MutableRow> rows) {
        for (int direction : FccNeighborhood.DIRECTIONS) {
            for (int cubeA = 0; cubeA < CUBES_PER_CELL; cubeA++) {
                for (int memberA = 0; memberA < MEMBERS_PER_CUBE; memberA++) {
                    for (int cubeB = 0; cubeB < CUBES_PER_CELL; cubeB++) {
                        for (int memberB = 0; memberB < MEMBERS_PER_CUBE; memberB++) {
                            sweepBinPairs(predicate, nLga, direction, cubeA,
                                          memberA, cubeB, memberB, rows);
                        }
                    }
                }
            }
        }
    }

    private static void sweepBinPairs(ContactPredicate predicate, int nLga,
                                       int direction, int cubeA, int memberA,
                                       int cubeB, int memberB,
                                       Map<RowKey, MutableRow> rows) {
        for (int binA = 0; binA < nLga; binA++) {
            float angleA = (float) binCenter(binA, nLga);
            for (int binB = 0; binB < nLga; binB++) {
                float angleB = (float) binCenter(binB, nLga);
                if (predicate.contacts(cubeA, memberA, angleA, cubeB, memberB,
                                       angleB, direction)) {
                    double minDistance = predicate.minDistance(cubeA, memberA,
                                                                angleA, cubeB,
                                                                memberB,
                                                                angleB,
                                                                direction);
                    RowKey key = new RowKey(direction, cubeA, memberA, cubeB,
                                            memberB, binA, binB);
                    rows.put(key, new MutableRow(key, true, minDistance));
                }
            }
        }
    }

    /**
     * Drives a real Phase A {@link AuditedRun} for {@code ticksObserved}
     * ticks, snapshotting {@code angle} pre-tick and correlating each
     * tick's resolved contacts against it - see class Javadoc.
     */
    private static void runDynamicReachability(ContactPredicate predicate,
                                                 int nLga, Point3i extent,
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

            for (CollisionSweep.AppliedCollision applied : outcome.collisionResult()
                                                                    .applied()) {
                recordObservedContact(automaton, predicate, nLga,
                                       preTickAngles, applied.contact(), rows);
            }
        }
    }

    private static void recordObservedContact(Necronomata automaton,
                                                ContactPredicate predicate,
                                                int nLga,
                                                float[] preTickAngles,
                                                Contact contact,
                                                Map<RowKey, MutableRow> rows) {
        int indexA = automaton.indexOfCell(contact.cellA())
                     + contact.cubeA() * MEMBERS_PER_CUBE + contact.memberA();
        int indexB = automaton.indexOfCell(contact.cellB())
                     + contact.cubeB() * MEMBERS_PER_CUBE + contact.memberB();
        int binA = binOf(preTickAngles[indexA], nLga);
        int binB = binOf(preTickAngles[indexB], nLga);

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
     * src/test/resources/lga/contact-atlas-v1.tsv} explicitly; this
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
        for (ContactAtlas.Row row : atlas.rows()) {
            if (row.contact()) {
                contactRows++;
            }
            if (row.observedCount() > 0) {
                observedRows++;
            }
        }

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
    }
}
