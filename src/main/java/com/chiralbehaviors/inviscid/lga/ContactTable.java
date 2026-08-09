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
import java.nio.file.Path;

import com.chiralbehaviors.inviscid.PhiCoordinates;

/**
 * The contact-predicate TABLE (bead inviscid-0nx.19, Phase C.2),
 * transcribed once from a {@link ContactAtlas} into a single packed-bitset
 * lookup ({@link #contacts}).
 *
 * <h2>Role in the formal LGA (fix-round correction, bead inviscid-0nx.23):
 * NOT the firing path</h2>
 * {@link LatticeGasAutomaton}'s USER-DECISION-FINAL contact-firing
 * mechanism is {@link FineStepContactTable}, re-derived from LIVE {@link
 * ContactPredicate} geometry at construction (57.8M evaluations, see that
 * class's Javadoc and {@code LatticeGasAutomaton}'s class Javadoc,
 * "Contact firing") -- NOT this table's bin-level ANY-OVERLAP bitset. This
 * table's {@link #contacts} lookup is never consulted in the tick path; the
 * table is still constructed and HEADER-CHECKED (see "Header pairing" on
 * {@code LatticeGasAutomaton}), so a stale atlas still fails loudly at
 * construction, but it is otherwise dead in the firing path. That earlier
 * revision (bin-level ANY-OVERLAP firing) was measured to produce a ~6.6x
 * collision-rate gap against the hybrid (bead inviscid-0nx.21's test 8),
 * root-caused to ANY-OVERLAP's deliberate over-inclusiveness -- not a bug
 * in this class, but a property of what "transcribed once, no run-time
 * geometry" trades away for {@code LatticeGasAutomaton}'s specific tick-path
 * needs. This class's OWN transcription/loading contract below (ANY-OVERLAP
 * semantics, the REFUSES-to-load guard, the packed-index layout) is
 * unchanged and still exactly what it says -- only the claim that {@code
 * LatticeGasAutomaton} fires off it is corrected.
 *
 * <h2>Transcription semantics: ANY-OVERLAP (bead inviscid-gyt, USER
 * DECISION 2026-08-08)</h2>
 * A cell {@code (direction, cubeA, memberA, binA, cubeB, memberB, binB)}
 * is transcribed as contacting iff its atlas row's {@code overlapFraction
 * > 0} - NOT the atlas's {@code contact} column (bin-center-only
 * evaluation, retained on the atlas for comparability, per {@link
 * ContactAtlas}'s own Javadoc). Bin-center contact is always a subset of
 * any-overlap (the bin center is itself one of the fine sweep's own
 * samples - see {@code ContactAtlasGenerator.sweepOverlapAndCenter}'s
 * proof), so a table transcribed this way never DROPS a bin-center
 * contact; it additionally fires on cells the hybrid automaton's
 * continuous trajectory can clip even though the bin center itself does
 * not contact (bead inviscid-0nx.16 stage 2's 1,581-row finding - see bead
 * inviscid-gyt for the full quantitative rationale). This bead's own
 * originally-specified fidelity test ("table verdict == ContactPredicate
 * at the bin centre") is superseded accordingly - see {@code
 * ContactTableTest} for both the exact-transcription test against {@code
 * overlapFraction} and the one-sided bin-centre-implies-table-contact
 * check.
 *
 * <h2>Loading and the REFUSES-to-load contract</h2>
 * {@link #load(Path, ContactAtlas.Header)} is the sole public entry point
 * that carries the bead's "REFUSES to load a stale atlas" requirement: it
 * delegates to {@link ContactAtlas#readValidated(Path, ContactAtlas.Header)},
 * which throws {@link ContactAtlas.HeaderMismatchException} - naming every
 * mismatched field, not just the first - the moment the atlas file's
 * on-disk header disagrees with the {@code expected} header a caller
 * supplies. A caller pairs this table with a geometry consumer - a
 * {@link MemberGeometry}, or the {@link FineStepContactTable} that {@link
 * LatticeGasAutomaton} actually fires from - by building BOTH from the
 * exact same {@link ContactAtlas.Header} instance, e.g. {@code
 * new MemberGeometry(header.geometryResolution(), header.memberRadius())}
 * alongside {@code ContactTable.load(path, header)}. Deriving every
 * consumer from one header instance is the single-source-of-truth
 * convention this class relies on rather than independently re-deriving
 * or duplicating any resolution / divisibility checks of its own.
 *
 * <p><b>Which divisibility invariants survive, precisely (bead
 * inviscid-0nx.27, E.0 - do not paraphrase this as "the constructor
 * guards divisibility").</b> Two DIFFERENT invariants are at issue and
 * only one of them is still enforced anywhere in this tree:
 * <ul>
 * <li>{@code phaseResolution % geometryResolution == 0} - ENFORCED, live,
 * by {@link LatticeGasAutomaton}'s constructor. This is the guard that
 * makes {@code accumulator / subBinSteps} exact integer arithmetic.</li>
 * <li>{@code geometryResolution % nLga == 0} - <b>NO LONGER ENFORCED BY
 * ANYTHING.</b> It belonged to the phase-quantizer class bead
 * inviscid-0nx.27 retired, and died with it; no successor guard was
 * added. Losing it is not a runtime regression today, because {@code
 * N_lga} has no production consumer at all ({@link #contacts} is
 * test-only; {@link FineStepContactTable} is the firing path). It is a
 * live hazard the moment {@code N_lga} becomes a swept parameter - see
 * {@code NLgaCandidateCampaign} (bead inviscid-eho), and note that its
 * candidates ALREADY include one that violates the retired invariant:
 * measured, {@code 360 % 16 == 8}, so {@code nLga=16} does not divide
 * the 360-step geometry resolution ({@code 8}, {@code 12} and {@code 24}
 * do). The SURVIVING invariant excludes MORE of those candidates than the
 * retired one does: {@code ContactAtlasGenerator.generate} stamps {@code
 * subBinSteps = 150} into every header unconditionally, so {@code
 * phaseResolution = nLga * 150}, and {@code phaseResolution % 360} is
 * {@code 120} at {@code nLga=8}, {@code 0} at {@code 12}, {@code 240} at
 * {@code 16}, {@code 0} at {@code 24} - {@code 12} and {@code 24} are
 * LGA-constructible, {@code 8} and {@code 16} are not. {@code 3600} is
 * the fine-phase resolution only at {@code nLga=24}. An earlier version
 * of this paragraph asserted all four "happen to divide 360", and a later
 * one that all four divide "the 3600 fine-phase resolution"; the first
 * is arithmetically false, the second true but irrelevant - the
 * candidates are never checked against a fixed 3600.
 * {@link ContactAtlasGenerator#binOfStep}'s own Javadoc
 * names exactly that {@code nLga=16 / geometryResolution=360} case as
 * one it is well-defined for, at the cost of bins that are unequal-width
 * in step-count terms; whether that is acceptable for the campaign's
 * purpose is not settled here.</li>
 * </ul>
 * {@link #of(ContactAtlas)} is the unvalidated
 * escape hatch for a caller that has already validated (or otherwise
 * trusts) a {@link ContactAtlas} instance directly - prefer {@link #load}
 * whenever a caller is loading from disk.
 *
 * <h2>Packed-index design</h2>
 * The 7-tuple {@code (direction, cubeA, memberA, binA, cubeB, memberB,
 * binB)} is packed into a single non-negative {@code long} via standard
 * mixed-radix positional encoding (radixes, in packing order: 12
 * directions, {@link PhiCoordinates#Cubes}{@code .length} cubes, 6
 * members, {@code N_lga} bins, cubes again, members again, {@code N_lga}
 * bins again) - injective (in fact bijective onto {@code [0,
 * domainSize)}) by the standard mixed-radix-numeral-system argument,
 * PROVIDED every component is validated to stay within its declared
 * radix, which {@link #packIndex} does before combining them. This
 * packing order DELIBERATELY MATCHES {@link #contacts(int, int, int, int,
 * int, int, int)}'s public parameter order - {@code (cubeA, memberA,
 * binA)} grouped together, then {@code (cubeB, memberB, binB)} - see
 * {@link #packIndex}'s own Javadoc for the code-review history behind
 * that alignment. {@code direction} (one of {@code +/-1..+/-6}, never
 * {@code 0}) is mapped to a dense {@code [0, 12)} index via {@link
 * #directionIndex}: positive directions {@code 1..6} map to {@code 0..5};
 * negative directions {@code -1..-6} map to {@code 6..11} - matching
 * {@link FccNeighborhood#DIRECTIONS}'s own enumeration order, though this
 * class computes the mapping arithmetically rather than via that list, to
 * keep the hot lookup path allocation-free and free of any list/map
 * indexing. {@link #contacts} is exactly one {@code packIndex} call
 * (fixed arithmetic, bounds checks only) plus one array-index bit test -
 * O(1) and allocation-free by construction, verified empirically for the
 * allocation half by {@code
 * ContactTableTest.lookupIsConstantTimeAndAllocationFree}.
 * <p>
 * Storage is a single {@code long[]} bitset, one bit per domain cell:
 * {@code domainSize = 12 * cubeCount * 6 * cubeCount * 6 * N_lga *
 * N_lga}. At the committed atlas's {@code N_lga=24} and {@code cubeCount
 * =5}: {@code domainSize = 12*5*6*5*6*24*24 = 6,220,800} bits {@code =
 * 97,200} longs {@code = 777,600} bytes ({@code ~759.4 KiB}) -
 * measured directly via {@link #memoryFootprintBytes()}, not merely
 * computed, since a caller needs the ACTUAL allocated size, which is
 * {@code ceil(domainSize / 64) * 8} bytes (the ceiling rounds up to the
 * next whole {@code long} when {@code domainSize} is not itself a
 * multiple of 64 - though it happens to be exactly one at the committed
 * {@code N_lga=24}).
 *
 * @author halhildebrand
 */
public final class ContactTable {

    private static final int MEMBERS_PER_CUBE = 6;
    private static final int CUBE_COUNT       = PhiCoordinates.Cubes.length;
    private static final int DIRECTION_COUNT  = 12;

    private final int                nLga;
    private final long               domainSize;
    private final long[]             bits;
    private final long               firedCount;
    private final ContactAtlas.Header header;

    private ContactTable(int nLga, long domainSize, long[] bits,
                          long firedCount, ContactAtlas.Header header) {
        this.nLga = nLga;
        this.domainSize = domainSize;
        this.bits = bits;
        this.firedCount = firedCount;
        this.header = header;
    }

    /**
     * Loads a {@link ContactTable} from {@code atlasPath}, REFUSING (via
     * {@link ContactAtlas.HeaderMismatchException}) to build one from an
     * atlas whose on-disk header disagrees with {@code expected} - see
     * class Javadoc's "REFUSES-to-load contract".
     *
     * @param atlasPath
     *            path to a v2 {@link ContactAtlas} TSV file
     * @param expected
     *            the header a paired {@link MemberGeometry} /
     *            {@link FineStepContactTable} configuration was built from
     *            - the SAME {@link ContactAtlas.Header} instance a caller
     *            passes to those constructors, so both are provably built
     *            from identical parameters
     * @return a table transcribing {@code atlasPath}'s any-overlap cells
     * @throws IOException
     *             if {@code atlasPath} cannot be read
     * @throws ContactAtlas.HeaderMismatchException
     *             if the on-disk header is missing a required parameter,
     *             or any present parameter disagrees with {@code expected}
     *             - naming every mismatched field
     */
    public static ContactTable load(Path atlasPath,
                                     ContactAtlas.Header expected) throws IOException {
        return of(ContactAtlas.readValidated(atlasPath, expected));
    }

    /**
     * Builds a {@link ContactTable} directly from an already-in-hand
     * {@link ContactAtlas} - no header validation is performed here; a
     * caller using this entry point is responsible for having already
     * validated (or otherwise trusting) {@code atlas}'s provenance. Prefer
     * {@link #load(Path, ContactAtlas.Header)} when loading from disk.
     *
     * @param atlas
     *            the atlas to transcribe; {@code N_lga} is sourced from
     *            {@code atlas.header().phaseResolutionNLga()} - the single
     *            source of truth every {@code N_lga} consumer reads from,
     *            never a second hardcoded literal
     * @return a table with one bit set per row whose {@code
     *         overlapFraction > 0} (bead inviscid-gyt's any-overlap
     *         transcription signal)
     */
    public static ContactTable of(ContactAtlas atlas) {
        ContactAtlas.Header header = atlas.header();
        int nLga = header.phaseResolutionNLga();
        long domainSize = domainSize(nLga);
        long longsNeeded = (domainSize + 63) / 64;
        if (longsNeeded > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("nLga=" + nLga
                                                + " produces a domain ("
                                                + domainSize
                                                + " bits) too large for a long[]-backed bitset ("
                                                + longsNeeded
                                                + " longs exceeds Integer.MAX_VALUE)");
        }
        long[] bits = new long[(int) longsNeeded];
        long firedCount = 0;
        for (ContactAtlas.Row row : atlas.rows()) {
            if (row.overlapFraction() > 0.0) {
                long idx = packIndex(row.direction(), row.cubeA(),
                                     row.memberA(), row.phaseBinA(),
                                     row.cubeB(), row.memberB(),
                                     row.phaseBinB(), nLga);
                if (!getBit(bits, idx)) {
                    firedCount++;
                    setBit(bits, idx);
                }
            }
        }
        return new ContactTable(nLga, domainSize, bits, firedCount, header);
    }

    /**
     * The hot-path lookup: no allocation, no branching beyond {@link
     * #packIndex}'s bounds checks - see class Javadoc.
     *
     * @return {@code true} iff this table transcribed the {@code
     *         (direction, cubeA, memberA, binA, cubeB, memberB, binB)}
     *         cell as contacting (any-overlap: the atlas row's {@code
     *         overlapFraction > 0})
     * @throws IllegalArgumentException
     *             if any argument is outside its domain (see {@link
     *             #packIndex})
     */
    public boolean contacts(int direction, int cubeA, int memberA, int binA,
                            int cubeB, int memberB, int binB) {
        long idx = packIndex(direction, cubeA, memberA, binA, cubeB, memberB,
                             binB, nLga);
        return getBit(bits, idx);
    }

    /**
     * @return {@code N_lga}, sourced from the transcribed atlas's header
     */
    public int nLga() {
        return nLga;
    }

    /**
     * @return the {@link ContactAtlas.Header} this table was built from -
     *         the cross-verification hook a caller pairing this table with
     *         a geometry consumer should use: assert {@code
     *         table.header().equals(header)} where {@code header} is the
     *         SAME {@link ContactAtlas.Header} instance that consumer was
     *         built from (code-review follow-up, inviscid-0nx.19 Important
     *         finding 3: {@link #load(Path, ContactAtlas.Header)}'s
     *         REFUSES contract only checks the ON-DISK atlas against the
     *         caller-supplied {@code expected} header - it does nothing to
     *         verify {@code expected} itself is the header a co-existing
     *         consumer was actually built from. This accessor makes that
     *         cross-check POSSIBLE for a caller to perform; it does not
     *         perform the check itself). {@link LatticeGasAutomaton}'s
     *         constructor is where that check is actually performed at
     *         runtime today - see its "Header pairing" Javadoc.
     */
    public ContactAtlas.Header header() {
        return header;
    }

    /**
     * @return the full 7-tuple domain size this table partitions - see
     *         class Javadoc's packed-index design
     */
    public long domainSize() {
        return domainSize;
    }

    /**
     * @return how many distinct domain cells this table transcribed as
     *         contacting (i.e. how many bits are set) - always {@code <=
     *         domainSize}, and strictly less for any real atlas (the
     *         atlas is sparse: it only ever records FIRED cells - see
     *         {@code ContactAtlasGenerator.sweepComboOverlap}'s {@code if
     *         (hits == 0) continue;})
     */
    public long firedCount() {
        return firedCount;
    }

    /**
     * @return the table's actual allocated backing-array size, in bytes
     *         ({@code bits.length * Long.BYTES}) - the MEASURED footprint
     *         the bead's acceptance criterion requires recording, not
     *         merely {@code domainSize / 8} (which would omit the
     *         ceiling-to-a-whole-{@code long} rounding).
     */
    public long memoryFootprintBytes() {
        return (long) bits.length * Long.BYTES;
    }

    /**
     * Packs the 7-tuple into a single non-negative bit index - see class
     * Javadoc's packed-index design. Package-private (not private) so
     * {@code ContactTableTest.packedIndexingRoundTrips} can verify
     * injectivity directly, exhaustively, over the full domain.
     * <p>
     * Parameter order DELIBERATELY MATCHES {@link #contacts(int, int, int,
     * int, int, int, int)}'s public signature (code-review follow-up,
     * inviscid-0nx.19: the two orders originally diverged - {@code
     * contacts} groups {@code (cubeA, memberA, binA)} then {@code (cubeB,
     * memberB, binB)}, while this method grouped {@code (cubeA, memberA,
     * cubeB, memberB)} then {@code (binA, binB)} - a correct-today but
     * silent maintenance footgun for any future edit to either signature).
     * Both the bit-packing order below and every call site ({@link
     * #contacts}, {@link #of(ContactAtlas)}, {@code
     * ContactTableTest.packedIndexingRoundTrips}) were updated together;
     * the resulting bit LAYOUT differs from before (an internal-only
     * change - {@link #domainSize(int)}'s total size and every public
     * behavior are unaffected, since a mixed-radix encoding's total range
     * does not depend on radix ORDER, only injectivity does, which {@link
     * ContactTableTest#packedIndexingRoundTrips()} re-verifies exhaustively
     * for the new layout).
     *
     * @throws IllegalArgumentException
     *             if {@code direction} is {@code 0} or outside {@code
     *             [-6,6]}, if either cube index is outside {@code [0,
     *             cubeCount)}, if either member index is outside {@code
     *             [0,6)}, or if either bin index is outside {@code [0,
     *             nLga)}
     */
    static long packIndex(int direction, int cubeA, int memberA, int binA,
                          int cubeB, int memberB, int binB, int nLga) {
        int dirIdx = directionIndex(direction);
        validateCube(cubeA);
        validateCube(cubeB);
        validateMember(memberA);
        validateMember(memberB);
        validateBin(binA, nLga);
        validateBin(binB, nLga);

        long idx = dirIdx;
        idx = idx * CUBE_COUNT + cubeA;
        idx = idx * MEMBERS_PER_CUBE + memberA;
        idx = idx * nLga + binA;
        idx = idx * CUBE_COUNT + cubeB;
        idx = idx * MEMBERS_PER_CUBE + memberB;
        idx = idx * nLga + binB;
        return idx;
    }

    /**
     * @return a dense {@code [0, 12)} index for {@code direction} - {@code
     *         1..6} map to {@code 0..5}; {@code -1..-6} map to {@code
     *         6..11} - see class Javadoc
     */
    static int directionIndex(int direction) {
        if (direction == 0 || direction < -6 || direction > 6) {
            throw new IllegalArgumentException("direction must be one of +/-1..+/-6: "
                                                + direction);
        }
        return direction > 0 ? direction - 1 : 5 - direction;
    }

    private static long domainSize(int nLga) {
        return (long) DIRECTION_COUNT * CUBE_COUNT * MEMBERS_PER_CUBE
               * CUBE_COUNT * MEMBERS_PER_CUBE * nLga * nLga;
    }

    private static boolean getBit(long[] bits, long idx) {
        return (bits[(int) (idx >>> 6)] & (1L << (idx & 63))) != 0;
    }

    private static void setBit(long[] bits, long idx) {
        bits[(int) (idx >>> 6)] |= (1L << (idx & 63));
    }

    private static void validateCube(int cube) {
        if (cube < 0 || cube >= CUBE_COUNT) {
            throw new IllegalArgumentException("cube: " + cube);
        }
    }

    private static void validateMember(int member) {
        if (member < 0 || member >= MEMBERS_PER_CUBE) {
            throw new IllegalArgumentException("member: " + member);
        }
    }

    private static void validateBin(int bin, int nLga) {
        if (bin < 0 || bin >= nLga) {
            throw new IllegalArgumentException("bin: " + bin);
        }
    }
}
