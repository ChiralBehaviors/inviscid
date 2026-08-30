# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Inviscid — cellular automata based on "jitterbugs" (Buckminster Fuller's jitterbug transformation), visualized in JavaFX 3D. Single-module Maven project, JavaFX 20.

## Commands

- Build: `mvn compile`
- Test: `mvn test`
- Single test: `mvn test -Dtest=AutomataTest`

No exec/javafx Maven plugin is configured — the animation apps are launched by running their `main` methods directly (IDE or classpath invocation). Runnable entry points:
- `com.chiralbehaviors.inviscid.animations.HoneycombBreatheAnimation` (a compact honeycomb patch doing the medium's one motion, bounded by elastic contact)
- `com.chiralbehaviors.inviscid.animations.ThreeCellPhaseAnimation` (three cells, each with its own phase, integrated in reduced coordinates) — a free tree fragment, so it tumbles and winds; see `inviscid-qvf.26`
- `com.chiralbehaviors.inviscid.animations.JitterbugAnimation`
- `com.chiralbehaviors.inviscid.animations.NecronomataAnimation`
- `com.javafx.experiments.jfx3dviewer.Jfx3dViewerApp` (generic 3D viewer)

The build compiles at source/target 20 (pom properties; the old maven-compiler-plugin 1.8 override was removed in 3ae3f77), so Java 16+ features like records are available. JavaFX 20 requires a JDK 17+ runtime.

## Architecture

Three layers, top to bottom:

1. **`com.chiralbehaviors.inviscid`** — the original code of this project.
   - Core model: `Necronomata` is the cellular automaton state — flat float arrays (angle, frequency, deltas) over a 3D `Point3i` extent, 30 values per cell, updated via a `Processor` functional interface. `Jitterbug` renders one jitterbug: the 8 faces of an `Octahedron` given independent `Rotate`/`Translate` transforms so faces open/close with a rotation angle. `CubicGrid` lays out cells on a cubic lattice (SIX or EIGHT neighborhood). `PhiCoordinates` and `LengthTable`/`Constants` hold golden-ratio-based geometry constants; `QuadRay` is a quadray (tetrahedral) coordinate system.
   - `animations/` — JavaFX apps. They extend `Jfx3dViewerApp` (the vendored viewer below) and override content creation; `PolyView` is the shared base that adds polyhedron edges/vertices to the scene. `Colors` holds shared materials.

2. **`mesh.*`, `util.*`, `math.*`** — vendored polyhedra geometry library. `mesh.polyhedra.Polyhedron` is the base; `plato/` (Tetrahedron…Icosahedron) and `archimedes/` build specific solids; `Polyhedron.toTriangleMesh()` bridges to JavaFX meshes. Uses `javax.vecmath` types (`Vector3d`, `Point3i`), not JavaFX geometry — conversion happens at the rendering boundary.

3. **`com.javafx.experiments.*`** — vendored Oracle JavaFX 3D sample code: the `jfx3dviewer` application (camera, timeline, FXML UI in `src/main/resources`), model importers (OBJ, Maya, 3DS Max, Collada), `shape3d` polygon meshes, and `utils3d` geometry/transform classes. Treat as third-party: don't restyle or refactor it; it changes only when the core layers need something from it.

The single test (`AutomataTest`) covers `Necronomata` neighbor/iteration logic — the automaton model is testable headlessly; everything touching JavaFX scene graph is not.

## Analysis (Python)

`analysis/` holds the jitterbug-medium derivation as self-gating Python scripts (numpy/scipy; interpreter `/opt/homebrew/opt/python@3.13/libexec/bin/python3`). `analysis/README.md` is the map: `model/` is the model of record (rigid struts, soft joints, voids empty), `model/double_covering/` the same line before the voids were emptied, `retired/` the closed lines (hard wall, rig lock, strut springs, attic), `history/` the long derivation record. `model/first_principles/` is the tied array built up one body at a time (one cell, face to face, the vertex point, the ring, the overdrive), with the page exporters in `pages/`. Run from the repo root: `analysis/gates.sh` (all), `analysis/gates.sh live`, or `python -m analysis.model.dispersion`. Every module's exit code is its verdict. The old two-letter names (`jb_rc`, `jb_bz`, …) that T2 and beads cite are listed beside each new name in the README.


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:6cd5cc61 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.

## Agent Context Profiles

The managed Beads block is task-tracking guidance, not permission to override repository, user, or orchestrator instructions.

- **Conservative (default)**: Use `bd` for task tracking. Do not run git commits, git pushes, or Dolt remote sync unless explicitly asked. At handoff, report changed files, validation, and suggested next commands.
- **Minimal**: Keep tool instruction files as pointers to `bd prime`; use the same conservative git policy unless active instructions say otherwise.
- **Team-maintainer**: Only when the repository explicitly opts in, agents may close beads, run quality gates, commit, and push as part of session close. A current "do not commit" or "do not push" instruction still wins.

## Session Completion

This protocol applies when ending a Beads implementation workflow. It is subordinate to explicit user, repository, and orchestrator instructions.

1. **File issues for remaining work** - Create beads for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Handle git/sync by active profile**:
   ```bash
   # Conservative/minimal/default: report status and proposed commands; wait for approval.
   git status

   # Team-maintainer opt-in only, unless current instructions forbid it:
   git pull --rebase
   git push
   git status
   ```
5. **Hand off** - Summarize changes, validation, issue status, and any blocked sync/commit/push step

**Critical rules:**
- Explicit user or orchestrator instructions override this Beads block.
- Do not commit or push without clear authority from the active profile or the current user request.
- If a required sync or push is blocked, stop and report the exact command and error.
<!-- END BEADS INTEGRATION -->
