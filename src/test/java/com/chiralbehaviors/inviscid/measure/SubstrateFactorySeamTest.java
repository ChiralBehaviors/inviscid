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

import static org.junit.Assert.assertEquals;

import javax.vecmath.Point3i;

import org.junit.Test;

import com.chiralbehaviors.inviscid.measure.AnisotropyProbe.SeedResult;

/**
 * Conformance tests for the {@link SubstrateFactory} campaign-runner
 * injection seam (bead inviscid-ckn / inviscid-0nx.21, T2
 * design-ckn-lattice-seam.md §4). The existing 5-arg
 * {@link AnisotropyProbe#runOneSeed} signature is unchanged and now
 * delegates to the 6-arg factory overload with
 * {@link AnisotropyProbe#phaseAHybridSubstrate} -- both call shapes MUST
 * produce byte-identical results, since one is defined in terms of the
 * other.
 *
 * @author halhildebrand
 */
public class SubstrateFactorySeamTest {

    /**
     * Non-vacuity + delegation-equivalence in one: the 5-arg call and
     * the explicit 6-arg factory call over the SAME seed must agree on
     * every field, AND (non-vacuity) must actually see collisions.
     */
    @Test
    public void fiveArgRunOneSeedDelegatesToPhaseAHybridSubstrateFactoryUnchanged() {
        Point3i extent = new Point3i(6, 6, 6);
        Point3i origin = AnisotropyProbe.nearestEvenParityCenter(extent);

        SeedResult viaDefault = AnisotropyProbe.runOneSeed(extent, 42L, 16,
                                                             100, origin);
        SeedResult viaFactory = AnisotropyProbe.runOneSeed(extent, 42L, 16,
                                                             100, origin,
                                                             AnisotropyProbe::phaseAHybridSubstrate);

        assertEquals(viaDefault.seed(), viaFactory.seed());
        assertEquals(viaDefault.totalCollisions(),
                     viaFactory.totalCollisions());
        assertEquals(viaDefault.effectiveCollisions(),
                     viaFactory.effectiveCollisions());
        assertEquals(viaDefault.transport(), viaFactory.transport());
        assertEquals(viaDefault.spectral(), viaFactory.spectral());

        assertEquals("non-vacuity: this seed must actually produce collisions",
                     246L, viaDefault.totalCollisions());
    }
}
