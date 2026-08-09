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

package com.chiralbehaviors.inviscid.measure;

import javax.vecmath.Point3i;

import com.chiralbehaviors.inviscid.QuantaField;

/**
 * Campaign-runner substrate injection seam (bead inviscid-ckn /
 * inviscid-0nx.21, T2 design-ckn-lattice-seam.md §4). Lets {@link
 * AnisotropyProbe#runOneSeed} / {@link AnisotropyProbe#runCampaign} drive
 * a substrate other than the Phase A hybrid without touching their
 * existing 5-arg / 4-arg signatures, which keep delegating to {@link
 * AnisotropyProbe#phaseAHybridSubstrate} -- so Phase A reproducibility
 * is preserved BY CONSTRUCTION, not by re-measurement.
 *
 * <p><b>RNG draw-order contract.</b> A {@link Substrate} implementation
 * that draws randomness MUST match {@link
 * AnisotropyProbe#phaseAHybridSubstrate}'s order (random angles, then
 * the localized quanta packet) if it wants the same seed to reproduce
 * the same trajectory -- reordering those draws changes every seeded
 * run even at an identical seed. See {@code
 * SeamGoldenCompatTest#runOneSeedThroughTheSeamMatchesPinnedPhaseANumerics}
 * for the pin this contract protects.
 *
 * @author halhildebrand
 */
@FunctionalInterface
public interface SubstrateFactory {

    /**
     * One seed's fully-wired, ready-to-run audited substrate: the
     * {@link QuantaField} instruments read, the {@link AuditedRun}
     * driving it per tick, and the {@link CollisionStatistics} it
     * records into.
     */
    record Substrate(QuantaField field, AuditedRun run,
                      CollisionStatistics statistics) {
    }

    Substrate create(Point3i extent, long seed, int packetQuanta,
                      Point3i originCell);
}
