/**
 * Copyright (c) 2026 Chiral Behaviors, LLC, all rights reserved.
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
package com.chiralbehaviors.inviscid;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import org.junit.Test;

/**
 * The wall between the two lines that share this tree.
 *
 * <p>
 * The <b>jitterbug medium</b> (model of record: {@code analysis/}; Java
 * expression: {@code Jitterbug}, {@code CubicGrid}, the golden-ratio
 * geometry, the {@code jitterbug} package, the medium animations) and the
 * <b>cellular-automaton line</b> ({@code automaton} — Necronomata,
 * QuantaField, lga, measure, and their visualizations) are distinct
 * projects with distinct dynamics. The CA renders jitterbug geometry; it
 * is not the jitterbug medium, and decisions in one line do not bind the
 * other (the conflation was made once — see
 * {@code analysis/notes/su2_boundary_conditions.md} §6/§9, 2026-09-01 —
 * and this test exists so structure, not prose, prevents the next one).
 *
 * <p>
 * The dependency direction is one-way by construction: the automaton may
 * import medium geometry (it renders with it); nothing outside
 * {@code automaton} may reference anything inside it. This test scans the
 * raw sources, so javadoc links and reflection strings count as
 * references too — the wall has no doc-comment-sized holes. The needle is
 * assembled at runtime so this file does not flag itself.
 *
 * <p>
 * Non-vacuity: the scan must SEE both sides of the wall at realistic
 * sizes. If the tree moves and this test scans nothing, it fails loudly
 * rather than skip-passing (the gate discipline of {@code analysis/}'s
 * gates.sh, applied to the Java tree).
 */
public class ArchitectureWallTest {

    private static final String NEEDLE = "com.chiralbehaviors.inviscid." + "automaton";

    @Test
    public void nothingOutsideTheAutomatonReferencesIt() throws IOException {
        List<Path> automaton = new ArrayList<>();
        List<Path> outside = new ArrayList<>();
        for (String root : new String[] { "src/main/java", "src/test/java" }) {
            Path r = Path.of(root);
            assertTrue("source root must exist for this wall to mean anything: " + r,
                       Files.isDirectory(r));
            try (Stream<Path> walk = Files.walk(r)) {
                walk.filter(p -> p.toString().endsWith(".java"))
                    .forEach(p -> (p.toString().contains("/automaton/") ? automaton
                                                                        : outside).add(p));
            }
        }
        // non-vacuity: both sides of the wall are populated at realistic size
        assertTrue("expected the CA line under automaton/ (saw " + automaton.size()
                   + " files); if it moved, move this wall with it",
                   automaton.size() >= 20);
        assertTrue("expected the medium side outside automaton/ (saw " + outside.size()
                   + " files); if the tree moved, this wall is scanning nothing",
                   outside.size() >= 30);

        List<String> offenders = new ArrayList<>();
        for (Path p : outside) {
            if (Files.readString(p, StandardCharsets.UTF_8).contains(NEEDLE)) {
                offenders.add(p.toString());
            }
        }
        assertEquals("files outside automaton/ reference the CA line — the medium must "
                     + "not depend on, link to, or name the automaton in code: " + offenders,
                     List.of(), offenders);
    }
}
