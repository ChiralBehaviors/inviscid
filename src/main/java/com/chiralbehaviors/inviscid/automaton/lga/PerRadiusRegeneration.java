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

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import javax.vecmath.Point3i;

/**
 * The per-radius artifact regeneration harness (bead inviscid-0nx.30,
 * E.3): given {@code (geometryResolution, memberRadius, nLga, extent,
 * seed, ticksObserved)}, produce the two artifacts a run at that radius
 * needs - a discovered-combos cache and a contact atlas - both {@code
 * r}-stamped, both under {@code target/}.
 *
 * <h2>Why a harness at all</h2>
 * {@code design-seeding-radius.md} §D-B makes the member radius a PHYSICAL
 * LEVER rather than a constant: {@code r} sets the contact-set measure,
 * hence the collision rate, hence the mean free path, hence the emergent
 * viscosity, and the epic's target is the collision-dominated (inviscid)
 * limit that {@code r=0.015}'s ballistic regime sits at the opposite end
 * of. Sweeping {@code r} is therefore not a parameter tweak but a
 * re-derivation of every {@code r}-dependent artifact: {@link
 * ContactComboCache} misses its checked-in cache on ANY {@code
 * (geometryResolution, memberRadius)} mismatch and falls back to a full
 * exhaustive re-discovery, and the atlas that {@link LatticeGasAutomaton}
 * reads its radius back out of must be regenerated to match. This class is
 * the single entry point that does both, consistently, in one place.
 *
 * <h2>The pin-capture rule, mechanized</h2>
 * Committed artifacts are PINNED (design memo §Engineering-constraints:
 * "new {@code (r, rho)} physics produces atlas v3 + NEW versioned campaign
 * artifacts alongside old ... Never overwrite v2"). Two mechanisms enforce
 * that here rather than leaving it to caller discipline:
 * <ul>
 * <li>{@link #resolveOutputDirectory} REFUSES any output directory that is
 * not inside the process's {@code target/} tree - so no invocation,
 * mistaken or otherwise, can land on {@code src/test/resources/lga/
 * contact-atlas-v2.tsv} or {@code src/main/resources/lga/
 * discovered-combos-cache.tsv}. A rule enforced only by a careful human is
 * not a mechanism. The check resolves SYMLINKS <em>and</em> collapses
 * {@code ..} - each was independently demonstrated as a live escape past a
 * guard that handled only the other; see {@link #canonicalize}.</li>
 * <li>Filenames are {@code r}-stamped ({@link #atlasPathFor}, {@link
 * #combosCachePathFor}), so successive radii accumulate side by side
 * instead of overwriting one another.</li>
 * </ul>
 *
 * <h2>atlasVersion=3</h2>
 * Per the design memo, new-{@code (r, rho)}-physics atlases are stamped
 * {@code atlasVersion=3}. Note what this does and does not mean. The ROW
 * SCHEMA is UNCHANGED from v2 - same eleven columns, same order, same
 * semantics; {@code memberRadius} remains the header's physics-identity
 * field and is what {@link ContactAtlas#readValidated} actually compares.
 * The version bump is a provenance marker, not a format change.
 * CONSEQUENCE, stated plainly because it will bite the next caller:
 * {@code ContactAtlas.checkVersion} refuses to parse any {@code
 * atlasVersion} other than {@link ContactAtlas#ATLAS_VERSION} (== 2), so
 * an atlas written by this harness is WRITE-ONLY under the current reader.
 * Whichever bead first needs to LOAD a v3 atlas must widen that check;
 * {@code PerRadiusRegenerationTest.v3AtlasIsNotYetReadableByTheV2Reader}
 * pins the present state so the constraint surfaces as a named test rather
 * than a mid-campaign surprise.
 *
 * <h2>Cost</h2>
 * A regeneration at a new radius pays a full {@link
 * ContactComboCache#sweepExhaustively} - {@code 12 * 5*6 * 5*6 *
 * geometryResolution^2} worst case, several minutes at resolution 360 -
 * plus atlas generation. Run one at a time (the design memo's regeneration
 * discipline); no exec plugin is configured, so {@link #main} is invoked
 * directly from an IDE or classpath, matching {@link
 * ContactAtlasGenerator#main}'s own convention.
 *
 * @author halhildebrand
 */
public final class PerRadiusRegeneration {

    /**
     * The atlas format version stamped on regenerated artifacts. NOT
     * {@link ContactAtlas#ATLAS_VERSION}: that constant is what the READER
     * supports, and deliberately still says 2 so the committed v2 atlas
     * keeps loading. See class Javadoc.
     */
    public static final int ATLAS_VERSION_V3 = 3;

    /** The only directory tree this harness will ever write into. */
    private static final String OUTPUT_ROOT = "target";

    /**
     * Every parameter that makes a regeneration reproducible from its own
     * artifact headers alone.
     *
     * @param geometryResolution the member-geometry LUT resolution (must
     *                           be positive and divisible by 8 - {@link
     *                           MemberGeometry}'s contract)
     * @param memberRadius       {@code r}, the physical lever (design memo
     *                           §D-B); stamped into BOTH artifacts'
     *                           headers as the identity field
     * @param nLga               the atlas's phase-bin resolution
     * @param extent             the dynamic-reachability run's cell extent
     * @param seed               the dynamic-reachability run's RNG seed
     * @param ticksObserved      the dynamic-reachability run's tick count
     */
    public record Parameters(int geometryResolution, double memberRadius,
                              int nLga, Point3i extent, long seed,
                              int ticksObserved) {

        public Parameters {
            // isFinite FIRST: NaN fails every comparison, so a bare
            // `memberRadius <= 0` waves it through and the run emits
            // contact-atlas-v3-rNaN.tsv with a NaN physics-identity header.
            if (!Double.isFinite(memberRadius) || memberRadius <= 0) {
                throw new IllegalArgumentException("memberRadius must be positive and finite: "
                                                    + memberRadius);
            }
            extent = new Point3i(extent);
        }

        @Override
        public Point3i extent() {
            return new Point3i(extent);
        }
    }

    /**
     * What a regeneration produced.
     *
     * @param atlasPath        the {@code r}-stamped v3 atlas
     * @param combosCachePath  the {@code r}-stamped combos cache
     * @param comboCount       combos the exhaustive sweep discovered at
     *                         this {@code (resolution, radius)} - the
     *                         number a per-{@code r} test pin is measured
     *                         from
     * @param atlasRowCount    rows in the generated atlas
     * @param elapsedMs        wall time for the whole regeneration
     */
    public record Result(Path atlasPath, Path combosCachePath, int comboCount,
                          int atlasRowCount, long elapsedMs) {
    }

    private PerRadiusRegeneration() {
    }

    /**
     * @return the CANONICAL (absolute, symlink-resolved, {@code
     *         ..}-collapsed) form of {@code candidate}, if and only if it
     *         resolves inside the process's {@code target/} tree.
     *         Deliberately not {@code candidate} itself: returning the raw
     *         path would hand the write path a spelling containing {@code
     *         ..} and symlink components that was only ever CHECKED in
     *         canonical form, so every later {@code resolve} against it
     *         would re-traverse indirections the guard had already
     *         resolved away. The caller writes exactly what was approved.
     * @throws IllegalArgumentException naming {@code candidate} otherwise -
     *                                  the mechanism that makes it
     *                                  impossible for this harness to
     *                                  overwrite a pinned artifact (class
     *                                  Javadoc)
     * @throws IOException              if the filesystem cannot be
     *                                  consulted; see below for why it must
     *                                  be
     */
    public static Path resolveOutputDirectory(Path candidate) throws IOException {
        Path root = canonicalize(Path.of(OUTPUT_ROOT).toAbsolutePath());
        Path resolved = canonicalize(candidate.toAbsolutePath());
        if (!resolved.startsWith(root)) {
            throw new IllegalArgumentException("refusing to regenerate artifacts outside "
                                                + root
                                                + " - committed artifacts are PINNED (never overwrite contact-atlas-v2.tsv or the committed combos cache); requested: "
                                                + candidate + " (resolves to "
                                                + resolved + ")");
        }
        return resolved;
    }

    /**
     * The real, symlink-resolved location {@code absolute} denotes.
     *
     * <p>
     * WHY THIS CONSULTS THE FILESYSTEM. {@link Path#normalize} is purely
     * LEXICAL: it collapses {@code target/../src} but is blind to a symlink
     * {@code target/escape -> src/main/resources/lga}, which normalizes to
     * itself, passes a {@code startsWith("target")} test, and lands writes
     * on a PINNED artifact. That is not hypothetical - it was demonstrated
     * against the lexical-only version of this guard, which wrote both
     * artifacts into {@code src/main/resources/lga/}. A pin-capture rule
     * defeated by one {@code ln -s} is not a mechanism.
     *
     * <p>
     * The output directory usually does NOT exist yet, so {@link
     * Path#toRealPath} cannot be called on it directly. Instead: walk up to
     * the deepest ancestor that DOES exist, resolve THAT for real, and
     * re-append the not-yet-existing tail. Every symlink on the existing
     * prefix is therefore followed, while a path that is merely absent
     * stays comparable.
     *
     * <p>
     * WHY THE RESULT IS NORMALIZED, AND WHY THAT IS NOT REDUNDANT. The
     * re-appended tail is raw input, and {@code ..} inside it is a
     * component {@link Path#toRealPath} never got to see - so without the
     * final {@link Path#normalize} the tail escapes the very check the
     * symlink resolution above exists to enforce. Demonstrated against the
     * un-normalized version of this method, all three accepted: {@code
     * target/nonexistent/../../src/main/resources/lga}, {@code
     * target/a/b/../../../src/main/resources/lga}, {@code
     * target/./nope/../../src/test/resources/lga}. The first of those,
     * driven through {@link #regenerate} single-threaded with no attacker
     * at all, CREATED {@code src/main/resources/lga/pwned} inside the
     * pinned tree before dying on the file write; and with a concurrent
     * process creating {@code target/racy} during the (minutes-long) sweep,
     * {@code target/racy/../../src/main/resources/lga} succeeded and wrote
     * BOTH artifacts into the pinned tree. Normalizing here is sound
     * because {@code real} is already symlink-free, so a purely lexical
     * {@code ..} collapse now agrees with kernel semantics rather than
     * diverging from them - which is exactly the reason {@code normalize()}
     * alone, applied to the raw input, was never enough.
     *
     * <p>
     * ONE IMPRECISION, STATED RATHER THAN GLOSSED: {@link Files#exists}
     * FOLLOWS symlinks, so a DANGLING symlink in the path reports as
     * absent and its name is re-appended verbatim instead of being
     * resolved. That is not an escape - a write through a dangling link
     * fails rather than landing anywhere - but the tail is not, in
     * general, "a component that does not exist"; it is "a component whose
     * target does not exist", which is a weaker statement.
     */
    private static Path canonicalize(Path absolute) throws IOException {
        Path existing = absolute;
        int missing = 0;
        while (existing != null && !Files.exists(existing)) {
            existing = existing.getParent();
            missing++;
        }
        if (existing == null) {
            return absolute.normalize();
        }
        Path real = existing.toRealPath();
        for (int i = absolute.getNameCount() - missing; i < absolute.getNameCount();
             i++) {
            real = real.resolve(absolute.getName(i));
        }
        return real.normalize();
    }

    /**
     * @return the {@code r}-stamped v3 atlas path under {@code directory}
     */
    public static Path atlasPathFor(Path directory, double memberRadius) {
        return directory.resolve("contact-atlas-v3-r"
                                  + Double.toString(memberRadius) + ".tsv");
    }

    /**
     * @return the {@code r}-stamped combos-cache path under {@code
     *         directory} - deliberately NOT the bare {@code
     *         discovered-combos-cache.tsv} that {@link
     *         ContactComboCache#RESOURCE_PATH} names, so an {@code
     *         r}-stamped cache can never be mistaken for (or copied over)
     *         the committed one
     */
    public static Path combosCachePathFor(Path directory,
                                           double memberRadius) {
        return directory.resolve("discovered-combos-cache-r"
                                  + Double.toString(memberRadius) + ".tsv");
    }

    /**
     * Regenerates both artifacts into the default {@code target/}
     * directory.
     */
    public static Result regenerate(Parameters parameters) throws IOException {
        return regenerate(parameters, Path.of(OUTPUT_ROOT));
    }

    /**
     * Regenerates the combos cache and the v3 atlas at {@code
     * parameters}'s radius, into {@code directory} (which must pass
     * {@link #resolveOutputDirectory}).
     *
     * <p>
     * Order matters: the combos cache is swept and written FIRST, because
     * atlas generation itself calls {@link ContactComboCache#combosFor} for
     * the same {@code (resolution, radius)} pair. {@link
     * ContactComboCache#rebuild} publishes its sweep into that memo, so the
     * second call is a hit and the exhaustive discovery is paid once per
     * run rather than twice.
     *
     * <p>
     * This was NOT free, and the ordering alone never bought it. Until the
     * memo seeding was added, {@code rebuild} swept without publishing and
     * {@code generate}'s {@code combosFor} then missed the memo, missed the
     * classpath cache (whose header pins {@code (360, 0.015)}, so any
     * sweep radius mismatches it), and swept the whole space a second time
     * - measured at resolution 48, and the dominant term in the 205s that
     * an {@code r=0.05} regeneration at resolution 360 cost.
     */
    public static Result regenerate(Parameters parameters,
                                     Path directory) throws IOException {
        Path outputDirectory = resolveOutputDirectory(directory);
        Path combosCachePath = combosCachePathFor(outputDirectory,
                                                   parameters.memberRadius());
        Path atlasPath = atlasPathFor(outputDirectory,
                                       parameters.memberRadius());

        long start = System.nanoTime();
        List<ContactComboCache.Combo> combos = ContactComboCache.rebuild(combosCachePath,
                                                                          parameters.geometryResolution(),
                                                                          parameters.memberRadius());

        ContactAtlas generated = ContactAtlasGenerator.generate(parameters.nLga(),
                                                                 parameters.extent(),
                                                                 parameters.seed(),
                                                                 parameters.ticksObserved(),
                                                                 parameters.geometryResolution(),
                                                                 parameters.memberRadius(),
                                                                 ContactAtlasGenerator.resolveGitCommit());
        ContactAtlas versioned = new ContactAtlas(withVersion(generated.header(),
                                                                ATLAS_VERSION_V3),
                                                   generated.rows());
        versioned.write(atlasPath);
        long elapsedMs = (System.nanoTime() - start) / 1_000_000;

        return new Result(atlasPath, combosCachePath, combos.size(),
                           versioned.rows().size(), elapsedMs);
    }

    /**
     * Rebuilds {@code header} with a different {@code atlasVersion},
     * leaving every other field - crucially {@code memberRadius}, the
     * physics-identity field - untouched. Positional record construction
     * is deliberate: if {@link ContactAtlas.Header} ever gains a field,
     * this fails to COMPILE rather than silently dropping it.
     */
    private static ContactAtlas.Header withVersion(ContactAtlas.Header header,
                                                     int atlasVersion) {
        return new ContactAtlas.Header(atlasVersion, header.generatorClass(),
                                        header.gitCommit(),
                                        header.memberRadius(),
                                        header.geometryResolution(),
                                        header.cubeEdgeLength(),
                                        header.phaseResolutionNLga(),
                                        header.phiCoordinatesCubeSet(),
                                        header.extent(), header.seed(),
                                        header.ticksObserved(),
                                        header.subBinSteps());
    }

    /**
     * Regenerates one radius's artifacts: {@code args[0]} is {@code
     * memberRadius} (REQUIRED - there is deliberately no default, because
     * defaulting it would reintroduce exactly the hardcoded constant the
     * design memo's §D-B forbids), {@code args[1]} is {@code
     * geometryResolution} (default {@link
     * ContactAtlasGenerator#GEOMETRY_RESOLUTION}), {@code args[2]} is
     * {@code nLga} (default {@code 24}, matching the committed atlas's
     * {@code phaseResolutionNLga}), {@code args[3]} is {@code
     * ticksObserved} (default {@link
     * ContactAtlasGenerator#DEFAULT_TICKS}).
     *
     * <p>
     * Output always lands in {@code target/} (class Javadoc). Several
     * minutes at resolution 360; run ONE regeneration at a time.
     */
    public static void main(String[] args) throws IOException {
        if (args.length < 1) {
            System.err.println("usage: PerRadiusRegeneration <memberRadius> [geometryResolution] [nLga] [ticksObserved]");
            System.err.println("  memberRadius is REQUIRED - it is the physical lever, never a default");
            return;
        }
        double memberRadius = Double.parseDouble(args[0]);
        int geometryResolution = args.length > 1 ? Integer.parseInt(args[1])
                                                  : ContactAtlasGenerator.GEOMETRY_RESOLUTION;
        int nLga = args.length > 2 ? Integer.parseInt(args[2]) : 24;
        int ticks = args.length > 3 ? Integer.parseInt(args[3])
                                     : ContactAtlasGenerator.DEFAULT_TICKS;

        Result result = regenerate(new Parameters(geometryResolution,
                                                   memberRadius, nLga,
                                                   ContactAtlasGenerator.DEFAULT_EXTENT,
                                                   ContactAtlasGenerator.DEFAULT_SEED,
                                                   ticks));

        System.out.println("Regenerated at memberRadius=" + memberRadius
                            + " geometryResolution=" + geometryResolution
                            + " nLga=" + nLga + " ticksObserved=" + ticks
                            + " in " + result.elapsedMs() + "ms");
        System.out.println("  combos discovered: " + result.comboCount()
                            + " -> " + result.combosCachePath());
        System.out.println("  atlas rows: " + result.atlasRowCount()
                            + " (atlasVersion=" + ATLAS_VERSION_V3 + ") -> "
                            + result.atlasPath());
    }
}
