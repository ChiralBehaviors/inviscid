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

import static com.chiralbehaviors.inviscid.Constants.TWO_PI;

import com.chiralbehaviors.inviscid.PhiCoordinates;

/**
 * Phase quantization for the formal lattice-gas automaton (bead
 * inviscid-0nx.18, Phase C.1): continuous rotation angles are mapped onto
 * exactly {@code N_lga} discrete phase bins partitioning {@code [0, 2*pi)},
 * and the member-segment geometry for every {@code (cube, member, bin)}
 * triple is precomputed once at construction, so the running automaton
 * never calls trigonometry.
 *
 * <h2>Distinct from the visualization LUT</h2>
 * {@code N_lga} is the small, physics-model phase count chosen and
 * justified in inviscid-A.5 (24, per the committed contact atlas' {@code
 * phaseResolutionNLga} header field) - it is NOT {@code
 * Necronomata.PHASE_RESOLUTION} (3600), the independent, much finer LUT
 * {@code NecronomataVisualization} uses purely for smooth rendering. This
 * class must never reference {@code Necronomata.PHASE_RESOLUTION}, and the
 * visualization LUT is untouched by this class - the two quantizations are
 * separate concerns that happen to share a normalize-then-divide
 * mechanism, not a shared resolution.
 *
 * <h2>Angle normalization</h2>
 * {@link #bin(float)} reuses, rather than reinvents, the exact
 * normalize-into-{@code [0, 2*pi)}-before-dividing convention {@code
 * NecronomataVisualization.setState} and {@code MemberGeometry.stepOf}
 * already establish: {@code normalized = angle % TWO_PI; if (normalized <
 * 0) normalized += TWO_PI;}, using {@link
 * com.chiralbehaviors.inviscid.Constants#TWO_PI} (float precision)
 * throughout - never {@code 2 * Math.PI} (double), which is a documented,
 * independent source of quantization mismatch (see {@code
 * ContactComboCache#angleOf}).
 *
 * <h2>Geometry LUT provenance</h2>
 * Every {@link #segment(int, int, int)} entry is exactly {@code
 * MemberGeometry.memberSegment(cube, member, centre(bin))} for a {@link
 * MemberGeometry} built at the atlas' own {@code geometryResolution} and
 * {@code memberRadius} - the LUT is tied to the single verified geometry
 * source, not to a second, independent derivation.
 *
 * @author halhildebrand
 */
public final class PhaseQuantizer {

    private static final int MEMBERS_PER_CUBE = 6;

    /**
     * @return a {@link PhaseQuantizer} sourcing {@code N_lga}, the geometry
     *         LUT resolution, and the member radius from a single already-
     *         parsed atlas header - {@link
     *         ContactAtlas.Header#phaseResolutionNLga()} is the ONE place
     *         {@code N_lga} is defined; callers must not hardcode it a
     *         second time.
     */
    public static PhaseQuantizer of(ContactAtlas.Header header) {
        return new PhaseQuantizer(header.phaseResolutionNLga(),
                                   header.geometryResolution(),
                                   header.memberRadius());
    }

    private final float         binWidth;
    private final Segment[][][] lut;
    private final int           nLga;

    /**
     * @param nLga
     *            N_lga, the number of discrete phase bins partitioning
     *            {@code [0, 2*pi)}. Must be > 0. Prefer {@link
     *            #of(ContactAtlas.Header)} over calling this constructor
     *            directly with a literal, so {@code N_lga} stays sourced
     *            from the atlas header in exactly one place.
     * @param geometryResolution
     *            the {@link MemberGeometry} LUT resolution the segment
     *            geometry is sourced from (atlas header {@code
     *            geometryResolution}; must be > 0 and divisible by 8, per
     *            {@link MemberGeometry}'s own constructor requirement).
     *            {@code nLga} MUST evenly divide this value - see the
     *            {@code IllegalArgumentException} thrown below.
     * @param memberRadius
     *            member (strut) radius (atlas header {@code
     *            memberRadius}).
     * @throws IllegalArgumentException
     *             if {@code nLga} does not evenly divide {@code
     *             geometryResolution}. {@link #bin(float)} and {@code
     *             ContactAtlasGenerator.binOfStep(step, nLga,
     *             geometryResolution)} are independent derivations
     *             (continuous-angle division vs. integer step-index
     *             division) that provably agree ONLY when {@code nLga}
     *             divides {@code geometryResolution} exactly - substantive-
     *             critique mutation testing on inviscid-0nx.18 measured
     *             40/360 bin mismatches at {@code nLga=100},
     *             {@code geometryResolution=360} (today's production
     *             {@code nLga=24} agrees only because {@code 24 | 360}).
     *             A future N_lga revision that is not a divisor would
     *             otherwise silently misalign live {@code bin()}
     *             classification against the inviscid-0nx.19-transcribed
     *             contact table, corrupting collision lookups with no
     *             test catching it - so this fails loudly at construction
     *             instead.
     */
    public PhaseQuantizer(int nLga, int geometryResolution, double memberRadius) {
        if (nLga <= 0) {
            throw new IllegalArgumentException("nLga must be > 0: " + nLga);
        }
        if (geometryResolution % nLga != 0) {
            throw new IllegalArgumentException("nLga (" + nLga
                                                + ") must evenly divide geometryResolution ("
                                                + geometryResolution + ")");
        }
        this.nLga = nLga;
        this.binWidth = TWO_PI / nLga;
        MemberGeometry geometry = new MemberGeometry(geometryResolution, memberRadius);
        this.lut = buildLut(geometry);
    }

    /**
     * <b>Deliberate departure from the {@code stepOf} siblings:</b> this
     * method CLAMPS an out-of-range division result to {@code nLga - 1},
     * whereas {@code MemberGeometry.stepOf} and {@code
     * ContactAtlasGenerator.stepOf} both WRAP via {@code % resolution}.
     * That divergence is intentional, not an oversight: a bare wrap sends
     * the near-{@code 2*pi} float-overshoot sliver (the top-of-circle bug
     * fixed on this bead) to bin 0 instead of the true last bin. See
     * follow-up bead inviscid-ann for whether the {@code stepOf} siblings
     * should adopt the same clamp.
     *
     * @param angle
     *            continuous rotation angle, radians, any sign or
     *            magnitude
     * @return the phase bin (0..{@code nLga}-1) {@code angle} falls into
     */
    public int bin(float angle) {
        float normalized = angle % TWO_PI;
        if (normalized < 0) {
            normalized += TWO_PI;
        }
        // normalized < TWO_PI is guaranteed here, so any idx >= nLga is a
        // division-rounding artifact (not a real extra bin) - clamp rather
        // than wrap via % nLga, which would send the near-2*pi sliver to
        // bin 0 instead of the true last bin nLga - 1 (code review finding,
        // inviscid-0nx.18: reviewer-verified trigger at nLga=360).
        int idx = (int) (normalized / binWidth);
        return idx >= nLga ? nLga - 1 : idx;
    }

    /**
     * Computed as float-only multiplication ({@code (bin + 0.5f) *
     * binWidth}), NOT {@code ContactComboCache#angleOf}'s double-then-
     * narrow-to-float convention. Deliberate: the two forms differ by at
     * most 1 ULP (substantive-critique measured up to ~4.8e-7 rad at
     * nLga=100/360/1000/3600/7919), immaterial against this class' 1e-6
     * geometry tolerance ({@code lutMatchesMemberGeometryAtBinCentres}),
     * so the simpler float-only form was kept rather than switched.
     *
     * @param bin
     *            phase bin index, 0..{@code nLga}-1
     * @return the bin's representative CENTER angle, radians - {@code (bin
     *         + 0.5) * (2*pi/nLga)}. Centers, not left edges, are the
     *         round-trip-safe reconstruction (see {@code
     *         ContactComboCache#angleOf}'s javadoc for the empirical
     *         reason: a left-edge reconstruction self-rounds a hair below
     *         its intended bin for a substantial fraction of bins under
     *         this same quantization scheme).
     */
    public float centre(int bin) {
        validateBin(bin);
        return (bin + 0.5f) * binWidth;
    }

    /**
     * @return {@code N_lga}, the number of phase bins this quantizer
     *         partitions {@code [0, 2*pi)} into
     */
    public int nLga() {
        return nLga;
    }

    /**
     * @return the precomputed member segment for {@code (cube, member)} at
     *         {@code bin}'s center angle - a table lookup, no
     *         trigonometry at call time
     */
    public Segment segment(int cube, int member, int bin) {
        validateCube(cube);
        validateMember(member);
        validateBin(bin);
        return lut[cube][member][bin];
    }

    private Segment[][][] buildLut(MemberGeometry geometry) {
        Segment[][][] table = new Segment[PhiCoordinates.Cubes.length][MEMBERS_PER_CUBE][nLga];
        for (int cube = 0; cube < table.length; cube++) {
            for (int member = 0; member < MEMBERS_PER_CUBE; member++) {
                for (int bin = 0; bin < nLga; bin++) {
                    table[cube][member][bin] = geometry.memberSegment(cube, member,
                                                                        centre(bin));
                }
            }
        }
        return table;
    }

    private void validateBin(int bin) {
        if (bin < 0 || bin >= nLga) {
            throw new IllegalArgumentException("bin: " + bin);
        }
    }

    private void validateCube(int cube) {
        if (cube < 0 || cube >= PhiCoordinates.Cubes.length) {
            throw new IllegalArgumentException("cube: " + cube);
        }
    }

    private void validateMember(int member) {
        if (member < 0 || member >= MEMBERS_PER_CUBE) {
            throw new IllegalArgumentException("member: " + member);
        }
    }
}
