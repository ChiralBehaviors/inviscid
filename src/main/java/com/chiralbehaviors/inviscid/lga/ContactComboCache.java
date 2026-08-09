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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The EXHAUSTIVE, ever-contacting {@code (direction, cubeA, memberA, cubeB,
 * memberB)} combo universe (bead inviscid-gyt, C.2's negative-direction /
 * any-overlap follow-up), precomputed and cached rather than rediscovered
 * on every {@link ContactAtlasGenerator#generate} call.
 *
 * <h2>Why precomputed, not live</h2>
 * {@link #sweepExhaustively} is the exhaustive discovery every caller
 * relies on - including {@link NLgaCandidateCampaign}, which delegates its
 * own combo discovery to {@link #combosFor} rather than duplicating this
 * sweep (bead inviscid-gyt angle-quantization alignment follow-up; see
 * that class's own Javadoc "Angle-quantization alignment fix" section) -
 * for every one of the {@code 12 * 5*6 * 5*6 == 10,800}
 * {@code (direction, cubeA, memberA, cubeB, memberB)} combinations, sweep
 * the full {@code geometryResolution x geometryResolution} angle grid and
 * keep every combination that contacts anywhere. Measured directly against
 * this codebase's {@link ContactPredicate} (a single non-contacting
 * combo's full {@code 360x360} sweep, ~130K evaluations, ~40ms cold): the
 * ~10,354 combos that never contact (the ones that must run to
 * COMPLETION to prove that, unlike the ~446 that short-circuit on first
 * hit) put the full discovery at several minutes wall time - fine for a
 * one-time, offline {@link #main} run, but far too slow to pay on every
 * {@code generate()} call (it would turn {@code ContactAtlasTest}'s
 * currently-fast generated-in-test contract suite into a multi-minute
 * one). So discovery is run ONCE (this class's {@link #main}), the result
 * checked in as {@value #RESOURCE_PATH}, and {@link #combosFor} loads that
 * resource - milliseconds, not minutes - whenever its header's {@code
 * geometryResolution}/{@code memberRadius} match the caller's; the live
 * {@link #sweepExhaustively} fallback exists for correctness/generality
 * (any {@code geometryResolution}/{@code memberRadius} pair other than the
 * cached one - not exercised by any current caller, all of which share
 * {@link ContactAtlasGenerator#GEOMETRY_RESOLUTION} / {@link
 * ContactAtlasGenerator#RADIUS}) and is what {@link #main} itself calls to
 * (re)build the cache.
 *
 * <h2>Why a coarse pre-scan was rejected</h2>
 * An earlier design considered a cheap coarse (e.g. 60-step) pre-scan to
 * cull the 10,800-combo search space before the full fine sweep. Rejected:
 * bead inviscid-0nx.16's own N_lga campaign found ~20% of the 446 real
 * combos have a "near-point" contact region narrower than even the
 * coarsest 45-degree candidate bin - plausibly narrower than a 6-degree
 * (60-step) pre-scan sample spacing too, which would silently reintroduce
 * exactly the kind of completeness gap this whole bead (any-overlap
 * transcription semantics) exists to close. The precomputed-cache
 * approach has zero such risk: {@link #sweepExhaustively} IS the full
 * {@code geometryResolution}-step sweep, just paid for once and reused.
 *
 * @author halhildebrand
 */
final class ContactComboCache {

    /**
     * One ever-contacting {@code (direction, cubeA, memberA, cubeB,
     * memberB)} combination - independent of any phase-bin resolution;
     * {@link ContactAtlasGenerator} quantizes into bins downstream.
     */
    record Combo(int direction, int cubeA, int memberA, int cubeB,
                 int memberB) {
    }

    static final String RESOURCE_PATH = "lga/discovered-combos-cache.tsv";

    private static final int CUBES_PER_CELL   = 5;
    private static final int MEMBERS_PER_CUBE = 6;

    private record CacheKey(int geometryResolution, double memberRadius) {
    }

    private static final Map<CacheKey, List<Combo>> CACHE = new ConcurrentHashMap<>();

    private ContactComboCache() {
    }

    /**
     * @return every combination that contacts anywhere on a {@code
     *         geometryResolution x geometryResolution} angle grid, for
     *         {@code predicate}'s underlying geometry (which must itself
     *         have been constructed at {@code geometryResolution}/{@code
     *         memberRadius} - not re-validated here, the caller's
     *         responsibility). Resolved from the checked-in cache
     *         resource when its header matches {@code (geometryResolution,
     *         memberRadius)}; swept live otherwise (see class Javadoc).
     *         Cached in memory per distinct {@code (geometryResolution,
     *         memberRadius)} pair for the lifetime of the JVM, so repeated
     *         calls (e.g. across every {@code nLga} candidate in a
     *         campaign) never re-pay either cost more than once.
     */
    static List<Combo> combosFor(ContactPredicate predicate, int geometryResolution,
                                  double memberRadius) {
        return CACHE.computeIfAbsent(new CacheKey(geometryResolution, memberRadius),
                                      key -> loadOrSweep(predicate, key));
    }

    private static List<Combo> loadOrSweep(ContactPredicate predicate, CacheKey key) {
        List<Combo> cached = tryLoad(key);
        if (cached != null) {
            return cached;
        }
        return sweepExhaustively(predicate, key.geometryResolution());
    }

    /**
     * The checked-in {@value #RESOURCE_PATH} resource itself, bypassing
     * the in-JVM memo {@link #combosFor} consults.
     *
     * <p>
     * WHY THIS EXISTS AS A SEPARATE ENTRY POINT. {@code
     * ContactComboCacheTest.liveSweepMatchesTheCommittedCacheAtProduction
     * Resolution} is the project's only detector of the committed cache
     * having gone stale against a {@link ContactPredicate}/{@link
     * MemberGeometry}/{@link com.chiralbehaviors.inviscid.PhiCoordinates}
     * algorithm change, and its subject is THE FILE. Reaching the file via
     * {@link #combosFor} was safe only for as long as the memo could be
     * populated from nowhere else; {@link #rebuild} publishing its SWEEP
     * into that same memo makes {@code combosFor} a file-or-sweep oracle,
     * so a single earlier {@link PerRadiusRegeneration#regenerate} call at
     * the committed {@code (360, 0.015)} pair in the same JVM would have
     * silently turned that tripwire into sweep-vs-sweep - a tautology that
     * can never go red. Reading the resource directly makes the tripwire
     * independent of what any earlier caller did to the memo.
     *
     * @return the committed combos, or {@code null} if the resource is
     *         absent or its header does not match {@code
     *         (geometryResolution, memberRadius)} - the caller decides
     *         which of those is an error
     */
    static List<Combo> loadCommittedCache(int geometryResolution,
                                           double memberRadius) {
        return tryLoad(new CacheKey(geometryResolution, memberRadius));
    }

    /**
     * @return the cached combo list if {@value #RESOURCE_PATH} is present
     *         on the classpath AND its header's {@code geometryResolution}
     *         / {@code memberRadius} match {@code key} exactly (a stale or
     *         foreign cache is never silently reused - same discipline as
     *         {@link ContactAtlas}'s own header staleness contract);
     *         {@code null} otherwise (resource absent, or present but for
     *         a different geometry).
     */
    private static List<Combo> tryLoad(CacheKey key) {
        try (InputStream in = ContactComboCache.class.getClassLoader()
                                                       .getResourceAsStream(RESOURCE_PATH)) {
            if (in == null) {
                return null;
            }
            List<String> lines;
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(in,
                                                                                    StandardCharsets.UTF_8))) {
                lines = reader.lines().toList();
            }
            return parse(lines, key);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static List<Combo> parse(List<String> lines, CacheKey key) {
        Map<String, String> kv = new LinkedHashMap<>();
        List<Combo> combos = new ArrayList<>();
        for (String line : lines) {
            if (line.isBlank()) {
                continue;
            }
            if (line.startsWith("#")) {
                int eq = line.indexOf('=');
                if (eq > 0) {
                    kv.put(line.substring(1, eq).trim(), line.substring(eq + 1)
                                                               .trim());
                }
                continue;
            }
            String[] parts = line.split("\t", -1);
            if (parts.length != 5) {
                throw new IllegalStateException("Malformed combo cache row (expected 5 tab-separated columns, found "
                                                 + parts.length + "): " + line);
            }
            combos.add(new Combo(Integer.parseInt(parts[0]),
                                  Integer.parseInt(parts[1]),
                                  Integer.parseInt(parts[2]),
                                  Integer.parseInt(parts[3]),
                                  Integer.parseInt(parts[4])));
        }

        if (!kv.containsKey("geometryResolution") || !kv.containsKey("memberRadius")
        || !kv.containsKey("comboCount")) {
            return null;
        }
        int cachedResolution = Integer.parseInt(kv.get("geometryResolution"));
        double cachedRadius = Double.parseDouble(kv.get("memberRadius"));
        if (cachedResolution != key.geometryResolution()
        || Double.compare(cachedRadius, key.memberRadius()) != 0) {
            return null;
        }
        int declaredCount = Integer.parseInt(kv.get("comboCount"));
        if (declaredCount != combos.size()) {
            throw new IllegalStateException("Combo cache row count (" + combos.size()
                                             + ") disagrees with its own header comboCount ("
                                             + declaredCount
                                             + ") - the cache resource is corrupt, regenerate via ContactComboCache.main");
        }
        return List.copyOf(combos);
    }

    /**
     * The exhaustive {@code 12 * 5*6 * 5*6 * geometryResolution^2}-worst-
     * case discovery sweep - see class Javadoc for why this is not the
     * path {@link #combosFor} takes for the cached {@code
     * (geometryResolution, memberRadius)} pair.
     */
    static List<Combo> sweepExhaustively(ContactPredicate predicate, int geometryResolution) {
        List<Combo> combos = new ArrayList<>();
        for (int direction : FccNeighborhood.DIRECTIONS) {
            for (int cubeA = 0; cubeA < CUBES_PER_CELL; cubeA++) {
                for (int memberA = 0; memberA < MEMBERS_PER_CUBE; memberA++) {
                    for (int cubeB = 0; cubeB < CUBES_PER_CELL; cubeB++) {
                        for (int memberB = 0; memberB < MEMBERS_PER_CUBE; memberB++) {
                            if (contactsAnywhere(predicate, direction, cubeA,
                                                  memberA, cubeB, memberB,
                                                  geometryResolution)) {
                                combos.add(new Combo(direction, cubeA, memberA,
                                                      cubeB, memberB));
                            }
                        }
                    }
                }
            }
        }
        return List.copyOf(combos);
    }

    private static boolean contactsAnywhere(ContactPredicate predicate, int direction,
                                             int cubeA, int memberA, int cubeB,
                                             int memberB, int geometryResolution) {
        for (int a = 0; a < geometryResolution; a++) {
            float angleA = angleOf(a, geometryResolution);
            for (int b = 0; b < geometryResolution; b++) {
                float angleB = angleOf(b, geometryResolution);
                if (predicate.contacts(cubeA, memberA, angleA, cubeB, memberB,
                                        angleB, direction)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * @return the representative angle (radians) of LUT {@code step} out of
     *         {@code resolution} - the STEP CENTER, {@code (step + 0.5) *
     *         (2*pi/resolution)}, computed with {@link
     *         com.chiralbehaviors.inviscid.Constants#TWO_PI} (float
     *         precision), NOT {@code 2 * Math.PI} (double). Both choices
     *         are load-bearing, discovered empirically (bead inviscid-gyt
     *         Phase A gate rework): {@link MemberGeometry#memberSegment}
     *         quantizes any angle back down to a step via {@code
     *         Constants.TWO_PI}-based, FLOAT-precision arithmetic
     *         (compounded float-modulo then float/double division) - a
     *         genuinely non-invertible quantization, not merely
     *         low-precision. Reconstructing the STEP'S LEFT EDGE ({@code
     *         step * (2*pi/resolution)}) round-trips incorrectly through
     *         that same quantization for ~38% of steps at {@code
     *         resolution=360} (self-rounding lands a hair below the
     *         intended step); the STEP CENTER is robust to that jitter
     *         (verified 0/360 mismatches) because it sits a half-step away
     *         from every boundary the jitter could cross. Using {@code
     *         2 * Math.PI} (double) instead of {@code Constants.TWO_PI}
     *         (float) for the reconstruction is a SEPARATE, independent
     *         failure mode (a different modulus produces a scaled angle
     *         that does not correspond to the intended step at all) - both
     *         fixes are required together. Every caller that evaluates
     *         {@link ContactPredicate#contacts} for a specific LUT step
     *         (the combo discovery sweep, {@link
     *         ContactAtlasGenerator#sweepOverlapAndCenter}) MUST go through
     *         this method rather than reconstructing the angle itself.
     */
    static float angleOf(int step, int resolution) {
        return (float) ((step + 0.5) * (com.chiralbehaviors.inviscid.Constants.TWO_PI
                                         / resolution));
    }

    /**
     * Rebuilds {@value #RESOURCE_PATH} via a real {@link
     * #sweepExhaustively} run: {@code args[0]} is {@code geometryResolution}
     * (default {@link ContactAtlasGenerator#GEOMETRY_RESOLUTION}), {@code
     * args[1]} is {@code memberRadius} (default {@link
     * ContactAtlasGenerator#RADIUS}), {@code args[2]} is the output path
     * (default {@code src/main/resources/}{@value #RESOURCE_PATH} - the
     * live, on-classpath resource {@link #combosFor} actually reads from
     * for every current caller's default geometry). Run manually (no exec
     * plugin configured - see {@code ContactAtlasGenerator.main}'s own
     * Javadoc for the same convention); several minutes wall time (class
     * Javadoc).
     */
    public static void main(String[] args) throws IOException {
        int geometryResolution = args.length > 0 ? Integer.parseInt(args[0])
                                                   : ContactAtlasGenerator.GEOMETRY_RESOLUTION;
        double memberRadius = args.length > 1 ? Double.parseDouble(args[1])
                                                : ContactAtlasGenerator.RADIUS;
        Path out = args.length > 2 ? Path.of(args[2])
                                    : Path.of("src", "main", "resources", "lga",
                                              "discovered-combos-cache.tsv");

        long start = System.nanoTime();
        List<Combo> combos = rebuild(out, geometryResolution, memberRadius);
        long elapsedMs = (System.nanoTime() - start) / 1_000_000;
        System.out.println("Discovered " + combos.size()
                            + " ever-contacting combos (geometryResolution="
                            + geometryResolution + ", memberRadius="
                            + memberRadius + ") in " + elapsedMs + "ms -> "
                            + out);
    }

    /**
     * Runs a real {@link #sweepExhaustively} discovery at {@code
     * (geometryResolution, memberRadius)} and writes the result to {@code
     * out} in this class's cache format - the shared body behind both
     * {@link #main} (which rebuilds the checked-in, on-classpath resource)
     * and {@link PerRadiusRegeneration} (which writes {@code r}-stamped
     * caches into {@code target/} for the {@code design-seeding-radius.md}
     * §D-B radius sweep).
     *
     * <p>
     * Deliberately takes the output path rather than defaulting one: the
     * ONLY caller allowed to write {@value #RESOURCE_PATH} is {@link
     * #main}, and it passes that path explicitly. Nothing here can
     * silently overwrite the committed cache.
     *
     * <p>
     * The sweep result is PUBLISHED into the in-JVM memo {@link #combosFor}
     * reads, so a caller that rebuilds and then generates in the same JVM
     * (that is {@link PerRadiusRegeneration#regenerate}) pays the
     * exhaustive discovery once rather than twice. Note the asymmetry, it
     * is deliberate: this method seeds the memo but never CONSULTS it, and
     * never routes through {@link #combosFor}. Consulting would mean that
     * at the committed {@code (360, 0.015)} pair the classpath cache loads
     * and {@link #main} writes that loaded file straight back out - a
     * rebuild that rebuilt nothing. A rebuild always sweeps.
     *
     * <p>
     * WHAT THE SEEDING COSTS, AND WHERE IT WAS PAID. Publishing a sweep
     * into the memo changes what {@link #combosFor} MEANS: before this,
     * the memo could only ever be filled from the checked-in file, so any
     * test comparing {@code combosFor} against a live sweep was
     * structurally file-vs-sweep. After this it is file-OR-sweep-vs-sweep,
     * and one earlier {@link PerRadiusRegeneration#regenerate} at {@code
     * (360, 0.015)} in the same JVM would reduce that comparison to a
     * tautology. No current caller does that - but E.4/E.5 regenerating a
     * baseline for side-by-side comparison is exactly the code that would.
     * The affected tripwire therefore reads the resource directly via
     * {@link #loadCommittedCache} rather than through the memo; see that
     * method.
     *
     * @return the discovered combos, in the order written
     */
    static List<Combo> rebuild(Path out, int geometryResolution,
                                double memberRadius) throws IOException {
        ContactPredicate predicate = new ContactPredicate(new MemberGeometry(geometryResolution,
                                                                              memberRadius));
        List<Combo> combos = List.copyOf(sweepExhaustively(predicate,
                                                             geometryResolution));
        CACHE.put(new CacheKey(geometryResolution, memberRadius), combos);
        write(out, geometryResolution, memberRadius, combos);
        return combos;
    }

    private static void write(Path out, int geometryResolution, double memberRadius,
                               List<Combo> combos) throws IOException {
        StringBuilder sb = new StringBuilder();
        sb.append("# ContactComboCache - exhaustive ever-contacting combo universe (bead inviscid-gyt)\n");
        sb.append("# geometryResolution=").append(geometryResolution).append('\n');
        sb.append("# memberRadius=").append(Double.toString(memberRadius))
          .append('\n');
        sb.append("# comboCount=").append(combos.size()).append('\n');
        sb.append("# columns=direction\tcubeA\tmemberA\tcubeB\tmemberB\n");
        for (Combo combo : combos) {
            sb.append(combo.direction()).append('\t').append(combo.cubeA())
              .append('\t').append(combo.memberA()).append('\t')
              .append(combo.cubeB()).append('\t').append(combo.memberB())
              .append('\n');
        }
        if (out.getParent() != null) {
            Files.createDirectories(out.getParent());
        }
        Files.writeString(out, sb.toString());
    }
}
