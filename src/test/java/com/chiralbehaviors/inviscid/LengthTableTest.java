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

import org.junit.Test;

/**
 * @author halhildebrand
 *
 */
public class LengthTableTest {

    @Test
    public void testIt() {
        OldLengthTable t = new OldLengthTable(360, 1);
        assertEquals(Math.sqrt(2) / 2, t.length(Math.PI / 2), 0.001);
        double degree = (2 * Math.PI) / 360;
        assertEquals(1.0, t.length(45 * degree), 0.001);
        assertEquals(Math.sqrt(2) / 2, t.length(90 * degree), 0.001);
        assertEquals(1.0, t.length(135 * degree), 0.001);
        assertEquals(Math.sqrt(2) / 2, t.length(180 * degree), 0.001);
        assertEquals(1.0, t.length(225 * degree), 0.001);
        assertEquals(Math.sqrt(2) / 2, t.length(270 * degree), 0.001);
        assertEquals(1.0, t.length(315 * degree), 0.001);
    }
}
