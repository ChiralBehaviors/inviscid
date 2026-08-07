# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Inviscid — cellular automata based on "jitterbugs" (Buckminster Fuller's jitterbug transformation), visualized in JavaFX 3D. Single-module Maven project, JavaFX 20.

## Commands

- Build: `mvn compile`
- Test: `mvn test`
- Single test: `mvn test -Dtest=AutomataTest`

No exec/javafx Maven plugin is configured — the animation apps are launched by running their `main` methods directly (IDE or classpath invocation). Runnable entry points:
- `com.chiralbehaviors.inviscid.animations.JitterbugAnimation`
- `com.chiralbehaviors.inviscid.animations.NecronomataAnimation`
- `com.javafx.experiments.jfx3dviewer.Jfx3dViewerApp` (generic 3D viewer)

Build quirk: `pom.xml` properties declare source/target 20, but the maven-compiler-plugin's explicit `<source>1.8</source>/<target>1.8</target>` configuration overrides them. JavaFX 20 requires a JDK 17+ runtime regardless.

## Architecture

Three layers, top to bottom:

1. **`com.chiralbehaviors.inviscid`** — the original code of this project.
   - Core model: `Necronomata` is the cellular automaton state — flat float arrays (angle, frequency, deltas) over a 3D `Point3i` extent, 30 values per cell, updated via a `Processor` functional interface. `Jitterbug` renders one jitterbug: the 8 faces of an `Octahedron` given independent `Rotate`/`Translate` transforms so faces open/close with a rotation angle. `CubicGrid` lays out cells on a cubic lattice (SIX or EIGHT neighborhood). `PhiCoordinates` and `LengthTable`/`Constants` hold golden-ratio-based geometry constants; `QuadRay` is a quadray (tetrahedral) coordinate system.
   - `animations/` — JavaFX apps. They extend `Jfx3dViewerApp` (the vendored viewer below) and override content creation; `PolyView` is the shared base that adds polyhedron edges/vertices to the scene. `Colors` holds shared materials.

2. **`mesh.*`, `util.*`, `math.*`** — vendored polyhedra geometry library. `mesh.polyhedra.Polyhedron` is the base; `plato/` (Tetrahedron…Icosahedron) and `archimedes/` build specific solids; `Polyhedron.toTriangleMesh()` bridges to JavaFX meshes. Uses `javax.vecmath` types (`Vector3d`, `Point3i`), not JavaFX geometry — conversion happens at the rendering boundary.

3. **`com.javafx.experiments.*`** — vendored Oracle JavaFX 3D sample code: the `jfx3dviewer` application (camera, timeline, FXML UI in `src/main/resources`), model importers (OBJ, Maya, 3DS Max, Collada), `shape3d` polygon meshes, and `utils3d` geometry/transform classes. Treat as third-party: don't restyle or refactor it; it changes only when the core layers need something from it.

The single test (`AutomataTest`) covers `Necronomata` neighbor/iteration logic — the automaton model is testable headlessly; everything touching JavaFX scene graph is not.
