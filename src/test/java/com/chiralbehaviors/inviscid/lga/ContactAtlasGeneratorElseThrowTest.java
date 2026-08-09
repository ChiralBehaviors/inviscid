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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import com.chiralbehaviors.inviscid.measure.AuditedRun.TickOutcome;
import com.chiralbehaviors.inviscid.measure.ConservationAudit;

/**
 * Conformance tests for {@link ContactAtlasGenerator#requireHybridApplied}
 * (code-review non-blocking suggestion on the inviscid-0nx.21 seam
 * checkpoint): the else-throw path a non-hybrid {@link TickDriver} takes
 * inside atlas generation, now directly testable without paying for a
 * real generation run.
 *
 * @author halhildebrand
 */
public class ContactAtlasGeneratorElseThrowTest {

    /** A {@link TickReport} with no {@link CollisionSweep} involvement at all. */
    private record FakeTickReport(int tick, long signedTransferTotal)
        implements TickReport {
    }

    /**
     * The else-throw path: a driver that is not {@code
     * CollisionSweep.TickResult}-backed must be refused loudly, naming
     * the actual reported class.
     */
    @Test
    public void requireHybridAppliedThrowsForANonTickResultReport() {
        TickOutcome outcome = new TickOutcome(new FakeTickReport(0, 0L),
                                               dummyAuditResult());

        IllegalStateException thrown = assertThrows(IllegalStateException.class,
                                                       () -> ContactAtlasGenerator.requireHybridApplied(outcome));
        assertTrue("expected the actual reported class named in the message: "
                   + thrown.getMessage(),
                   thrown.getMessage().contains(FakeTickReport.class.getName()));
    }

    /** Positive control: a real {@code CollisionSweep.TickResult} passes through unchanged. */
    @Test
    public void requireHybridAppliedPassesThroughARealTickResult() {
        CollisionSweep.TickResult real = new CollisionSweep.TickResult(3,
                                                                          List.of(),
                                                                          0L, 0L,
                                                                          0L);
        TickOutcome outcome = new TickOutcome(real, dummyAuditResult());

        List<CollisionSweep.AppliedCollision> applied = ContactAtlasGenerator.requireHybridApplied(outcome);
        assertSame(real.applied(), applied);
        assertEquals(0, applied.size());
    }

    private static ConservationAudit.AuditResult dummyAuditResult() {
        return new ConservationAudit.AuditResult(0, 0L, 0L, List.of());
    }
}
