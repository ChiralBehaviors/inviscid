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

import com.chiralbehaviors.inviscid.automaton.QuantaField;

/**
 * One tick of a lattice substrate, plus the field instruments read (bead
 * inviscid-ckn / inviscid-0nx.21, T2 design-ckn-lattice-seam.md §3).
 * {@link HybridAutomaton} adopts this by adding {@link #field()}; its
 * existing {@link HybridAutomaton#tick(int)} signature is unchanged
 * (Java covariant returns satisfy {@code TickReport tick(int)} since
 * {@link CollisionSweep.TickResult} implements {@link TickReport}).
 *
 * @author halhildebrand
 */
public interface TickDriver {

    TickReport tick(int tickNumber);

    QuantaField field();
}
