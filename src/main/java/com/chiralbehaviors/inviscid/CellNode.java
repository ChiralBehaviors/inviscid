/**
 * Copyright (c) 2018 Chiral Behaviors, LLC, all rights reserved.
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

import static com.chiralbehaviors.inviscid.Constants.*;
import javafx.scene.Group;
import javafx.scene.shape.MeshView;

/**
 * @author halhildebrand
 *
 */
public class CellNode extends Group {
    private final double       resolution;
    private final MeshView[][] struts;

    public CellNode(double resolution, MeshView[][] struts) {
        this.struts = struts;
        this.resolution = resolution;
    }

    public CellNode(double resolution, MeshView[][] struts,
                    double[] initialState) {
        this(resolution, struts);
        setState(initialState);
    }

    public void setState(double[] angles) {
        getChildren().clear();
        for (int i = 0; i < 6; i++) {
            double angle = angles[i] % TWO_PI;
            if (angle < 0) {
                angle = TWO_PI + angle;
            }
            int index = (int) (angle / resolution);
            getChildren().add(struts[i][index]);
        }
        requestLayout();
    }
}