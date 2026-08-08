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
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import javax.vecmath.Point3i;

/**
 * The Phase A -> Phase C handoff artifact (bead inviscid-0nx.16): an
 * in-memory model plus versioned TSV serialization of the contact-predicate
 * table Phase C transcribes into its formal LGA collision table.
 *
 * <h2>Format contract (locked by the bead description)</h2>
 * A versioned TSV. HEADER BLOCK: comment lines ({@code # key=value})
 * recording every parameter the table depends on - {@link #read(Path)}
 * REFUSES to parse a file missing any of {@link #HEADER_KEYS}, and {@link
 * #readValidated(Path, Header)} additionally REFUSES a file whose header
 * values disagree with an expected {@link Header} a consumer supplies -
 * both loudly, via {@link HeaderMismatchException}, never silently. A
 * stale-parameter atlas (e.g. generated at a different {@code
 * phaseResolutionNLga} than the consumer expects) must never be silently
 * consumed - that is the single most important property of this class (see
 * the bead's own text).
 * <p>
 * {@code geometryResolution} (bead inviscid-0nx.16.2) is one parameter
 * beyond the bead's original 10-key list: {@link
 * ContactAtlasGenerator}'s 6-arg {@code generate} overload exposes {@code
 * MemberGeometry}'s angle-quantization LUT resolution as a real tunable
 * that changes which bins geometrically contact, so it belongs in the
 * staleness contract exactly like {@code memberRadius} does - its absence
 * from {@link #HEADER_KEYS} was a silent-staleness gap, not a deliberate
 * omission.
 * <p>
 * {@code gitCommit} is resolved for real at generation time (see {@code
 * ContactAtlasGenerator#resolveGitCommit()}) - never a fixed literal -
 * specifically so {@link #readValidated(Path, Header)}'s comparison is
 * load-bearing: an atlas regenerated after {@code ContactPredicate} /
 * {@code MemberGeometry} / {@code PhiCoordinates} code changes gets a
 * different SHA and a stale-header check that can actually fire. A commit
 * of {@code "-dirty"} suffix means the working tree had uncommitted
 * changes when the atlas was generated.
 * <p>
 * COLUMNS, one row per GEOMETRICALLY CONTACTING {@code (direction, cubeA,
 * memberA, cubeB, memberB, phaseBinA, phaseBinB)} combination (see {@link
 * ContactAtlasGenerator}), tab-separated: {@code direction cubeA memberA
 * cubeB memberB phaseBinA phaseBinB contact minDistance observedCount}.
 * {@code phaseBinA}/{@code phaseBinB} are quantised at {@code
 * phaseResolutionNLga} (the AUTOMATON phase resolution) - NOT {@code
 * Necronomata.PHASE_RESOLUTION} (the 3600-step visualisation LUT); see the
 * bead's own text for why those two resolutions are deliberately distinct.
 * {@code contact} is the geometric ground truth from {@link
 * ContactPredicate} evaluated at the bin centers; {@code observedCount} is
 * how often a real Phase A run actually hit that cell (0 is legal and
 * informative - it distinguishes "geometrically possible but dynamically
 * unreached").
 *
 * <h2>Serialization convention</h2>
 * Mirrors {@code measure.BaselineSpectrumHarness}'s existing golden-artifact
 * TSV convention ({@code # key=value} header lines, a trailing {@code #
 * columns=...} line, then tab-separated data rows) rather than inventing a
 * new one. Numeric header/row fields that must round-trip exactly ({@code
 * memberRadius}, {@code cubeEdgeLength}, {@code minDistance}) are written
 * via {@link Double#toString(double)}, which Java guarantees produces the
 * shortest decimal that reparses to the identical {@code double} - not
 * {@code %.9e} (BaselineSpectrumHarness's own summary-row convention),
 * which is lossy and would silently break byte-for-byte round-tripping.
 *
 * @author halhildebrand
 */
public final class ContactAtlas {

    public static final int ATLAS_VERSION = 1;

    /**
     * Every header parameter the format contract requires (bead
     * inviscid-0nx.16 DESCRIPTION), in write order.
     */
    private static final List<String> HEADER_KEYS = List.of("atlasVersion",
                                                              "generatorClass",
                                                              "gitCommit",
                                                              "memberRadius",
                                                              "geometryResolution",
                                                              "cubeEdgeLength",
                                                              "phaseResolutionNLga",
                                                              "phiCoordinatesCubeSet",
                                                              "extent", "seed",
                                                              "ticksObserved");

    private static final String COLUMNS = "direction\tcubeA\tmemberA\tcubeB\tmemberB\tphaseBinA\tphaseBinB\tcontact\tminDistance\tobservedCount";

    /**
     * Every provenance parameter this atlas depends on. A consumer (Phase
     * C) that loads an atlas whose header does not match what it expects
     * would silently transcribe a stale table without this record existing
     * to check against (see {@link #readValidated(Path, Header)}).
     */
    public record Header(int atlasVersion, String generatorClass,
                          String gitCommit, double memberRadius,
                          int geometryResolution, double cubeEdgeLength,
                          int phaseResolutionNLga,
                          String phiCoordinatesCubeSet, Point3i extent,
                          long seed, int ticksObserved) {

        public Header {
            extent = new Point3i(extent);
        }

        /**
         * @return a defensive copy of {@code extent} - {@link Point3i} is
         *         mutable, matching {@code Contact}'s documented
         *         convention (see that record's Javadoc): a caller must
         *         not be able to corrupt an already-built {@link Header}
         *         by mutating a {@code Point3i} it still holds a
         *         reference to, in either direction (construction-time
         *         copy above, accessor-time copy here).
         */
        @Override
        public Point3i extent() {
            return new Point3i(extent);
        }
    }

    /**
     * One contacting {@code (direction, cubeA, memberA, cubeB, memberB,
     * phaseBinA, phaseBinB)} combination.
     */
    public record Row(int direction, int cubeA, int memberA, int cubeB,
                       int memberB, int phaseBinA, int phaseBinB,
                       boolean contact, double minDistance,
                       long observedCount) {
    }

    /**
     * Thrown by {@link #read(Path)} when the header block is missing a
     * required parameter, or by {@link #readValidated(Path, Header)} when
     * a present parameter disagrees with the caller's expected {@link
     * Header} - see class Javadoc. Always names every offending parameter
     * at once (not just the first), so a caller sees the full picture in
     * one failure.
     */
    public static class HeaderMismatchException extends RuntimeException {
        private static final long serialVersionUID = 1L;

        public HeaderMismatchException(String message) {
            super(message);
        }
    }

    private final Header    header;
    private final List<Row> rows;

    public ContactAtlas(Header header, List<Row> rows) {
        this.header = header;
        this.rows = List.copyOf(rows);
    }

    public Header header() {
        return header;
    }

    public List<Row> rows() {
        return rows;
    }

    public String toTsv() {
        StringBuilder sb = new StringBuilder();
        sb.append("# ContactAtlas - Phase A -> Phase C contact predicate transcription source (bead inviscid-0nx.16)\n");
        sb.append("# atlasVersion=").append(header.atlasVersion()).append('\n');
        sb.append("# generatorClass=").append(header.generatorClass())
          .append('\n');
        sb.append("# gitCommit=").append(header.gitCommit()).append('\n');
        sb.append("# memberRadius=")
          .append(Double.toString(header.memberRadius())).append('\n');
        sb.append("# geometryResolution=")
          .append(header.geometryResolution()).append('\n');
        sb.append("# cubeEdgeLength=")
          .append(Double.toString(header.cubeEdgeLength())).append('\n');
        sb.append("# phaseResolutionNLga=")
          .append(header.phaseResolutionNLga()).append('\n');
        sb.append("# phiCoordinatesCubeSet=")
          .append(header.phiCoordinatesCubeSet()).append('\n');
        sb.append("# extent=").append(header.extent().x).append(',')
          .append(header.extent().y).append(',').append(header.extent().z)
          .append('\n');
        sb.append("# seed=").append(header.seed()).append('\n');
        sb.append("# ticksObserved=").append(header.ticksObserved())
          .append('\n');
        sb.append("# columns=").append(COLUMNS).append('\n');
        for (Row row : rows) {
            sb.append(row.direction()).append('\t').append(row.cubeA())
              .append('\t').append(row.memberA()).append('\t')
              .append(row.cubeB()).append('\t').append(row.memberB())
              .append('\t').append(row.phaseBinA()).append('\t')
              .append(row.phaseBinB()).append('\t').append(row.contact())
              .append('\t').append(Double.toString(row.minDistance()))
              .append('\t').append(row.observedCount()).append('\n');
        }
        return sb.toString();
    }

    public void write(Path path) throws IOException {
        if (path.getParent() != null) {
            Files.createDirectories(path.getParent());
        }
        Files.writeString(path, toTsv());
    }

    /**
     * Parses {@code path} into a {@link ContactAtlas}, refusing loudly
     * (via {@link HeaderMismatchException}) if the header block is missing
     * any of {@link #HEADER_KEYS}. Performs no comparison against a
     * caller's expectations - see {@link #readValidated(Path, Header)} for
     * that.
     *
     * @throws HeaderMismatchException if any required header parameter is
     *                                 absent, or a data row is malformed
     */
    public static ContactAtlas read(Path path) throws IOException {
        return parse(Files.readAllLines(path));
    }

    /**
     * {@link #read(Path)}, then refuses loudly (via {@link
     * HeaderMismatchException}) if any parsed header parameter disagrees
     * with {@code expected} - the seam a Phase C consumer uses to make
     * sure it never silently transcribes a stale atlas (bead
     * inviscid-0nx.16's single most important property).
     *
     * @throws HeaderMismatchException if the header is missing a required
     *                                 parameter (per {@link #read(Path)}),
     *                                 or if any present parameter
     *                                 disagrees with {@code expected}
     */
    public static ContactAtlas readValidated(Path path,
                                              Header expected) throws IOException {
        ContactAtlas atlas = read(path);
        validateHeader(atlas.header(), expected);
        return atlas;
    }

    private static ContactAtlas parse(List<String> lines) {
        Map<String, String> kv = new LinkedHashMap<>();
        List<Row> rows = new ArrayList<>();
        for (String line : lines) {
            if (line.isBlank()) {
                continue;
            }
            if (line.startsWith("#")) {
                String content = line.substring(1).trim();
                int eq = content.indexOf('=');
                if (eq > 0) {
                    kv.put(content.substring(0, eq).trim(),
                           content.substring(eq + 1).trim());
                }
                continue;
            }
            rows.add(parseRow(line));
        }

        List<String> missing = new ArrayList<>();
        for (String key : HEADER_KEYS) {
            if (!kv.containsKey(key)) {
                missing.add(key);
            }
        }
        if (!missing.isEmpty()) {
            throw new HeaderMismatchException("Atlas header missing required parameter(s): "
                                               + missing
                                               + " - refusing to load a table with unknown provenance");
        }

        Header header = new Header(Integer.parseInt(kv.get("atlasVersion")),
                                    kv.get("generatorClass"),
                                    kv.get("gitCommit"),
                                    Double.parseDouble(kv.get("memberRadius")),
                                    Integer.parseInt(kv.get("geometryResolution")),
                                    Double.parseDouble(kv.get("cubeEdgeLength")),
                                    Integer.parseInt(kv.get("phaseResolutionNLga")),
                                    kv.get("phiCoordinatesCubeSet"),
                                    parseExtent(kv.get("extent")),
                                    Long.parseLong(kv.get("seed")),
                                    Integer.parseInt(kv.get("ticksObserved")));
        return new ContactAtlas(header, rows);
    }

    private static Row parseRow(String line) {
        String[] parts = line.split("\t", -1);
        if (parts.length != 10) {
            throw new HeaderMismatchException("Malformed atlas data row (expected 10 tab-separated columns, found "
                                               + parts.length + "): " + line);
        }
        return new Row(Integer.parseInt(parts[0]), Integer.parseInt(parts[1]),
                       Integer.parseInt(parts[2]), Integer.parseInt(parts[3]),
                       Integer.parseInt(parts[4]), Integer.parseInt(parts[5]),
                       Integer.parseInt(parts[6]),
                       Boolean.parseBoolean(parts[7]),
                       Double.parseDouble(parts[8]),
                       Long.parseLong(parts[9]));
    }

    private static Point3i parseExtent(String value) {
        String[] parts = value.split(",");
        if (parts.length != 3) {
            throw new HeaderMismatchException("Malformed extent header value (expected \"x,y,z\"): "
                                               + value);
        }
        return new Point3i(Integer.parseInt(parts[0].trim()),
                           Integer.parseInt(parts[1].trim()),
                           Integer.parseInt(parts[2].trim()));
    }

    private static void validateHeader(Header actual, Header expected) {
        List<String> mismatches = new ArrayList<>();
        checkEquals(mismatches, "atlasVersion", expected.atlasVersion(),
                    actual.atlasVersion());
        checkEquals(mismatches, "generatorClass", expected.generatorClass(),
                    actual.generatorClass());
        checkEquals(mismatches, "gitCommit", expected.gitCommit(),
                    actual.gitCommit());
        checkEquals(mismatches, "memberRadius", expected.memberRadius(),
                    actual.memberRadius());
        checkEquals(mismatches, "geometryResolution",
                    expected.geometryResolution(),
                    actual.geometryResolution());
        checkEquals(mismatches, "cubeEdgeLength", expected.cubeEdgeLength(),
                    actual.cubeEdgeLength());
        checkEquals(mismatches, "phaseResolutionNLga",
                    expected.phaseResolutionNLga(),
                    actual.phaseResolutionNLga());
        checkEquals(mismatches, "phiCoordinatesCubeSet",
                    expected.phiCoordinatesCubeSet(),
                    actual.phiCoordinatesCubeSet());
        checkEquals(mismatches, "extent", expected.extent(), actual.extent());
        checkEquals(mismatches, "seed", expected.seed(), actual.seed());
        checkEquals(mismatches, "ticksObserved", expected.ticksObserved(),
                    actual.ticksObserved());
        if (!mismatches.isEmpty()) {
            throw new HeaderMismatchException(String.format(Locale.ROOT,
                                                              "Atlas header mismatch against expected parameters: %s",
                                                              String.join("; ",
                                                                          mismatches)));
        }
    }

    private static void checkEquals(List<String> mismatches, String key,
                                     Object expected, Object actual) {
        if (!expected.equals(actual)) {
            mismatches.add(key + ": expected " + expected + " but found "
                            + actual);
        }
    }
}
