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

package com.chiralbehaviors.inviscid;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3i;

import org.junit.Test;

/**
 * @author halhildebrand
 *
 */
public class AutomataTest {
    @Test
    public void testIteration() {
        Necronomata automata = new Necronomata(7, 6, 6);
        List<Point3i> loop = new ArrayList<>();
        
        automata.forEach(c -> loop.add(c));
        
        List<Point3i> collected = new ArrayList<>();
        for (Point3i c: automata) {
            collected.add(c);
        }
        
        assertEquals(loop, collected);
        
    }
}
