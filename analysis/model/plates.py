"""plates -- the plate geometry the MEDIUM work needs, and nothing else.

WHY THIS FILE EXISTS. The honeycomb modules -- jb_hc, jb_bt, jb_fl -- were
importing `jb_z_quasistatic_array`, all 5164 lines and 276 KB of it, to reach
exactly two things: `corners`, which is not even jb_z's (it re-exports
`jb_a_family`'s), and `plate_normal`, which is the two lines below. The only
other reference, `_is_piercing`, appears solely inside a docstring warning
against using it.

That import was expensive, not merely untidy. jb_z carries the quasistatic
contact/crank machinery and a 152-row gate whose ARRAY results are all computed
on the single-vertex-contact topology inviscid-ia5 retired; every edit anywhere
in that line invalidated the jb_cache entries the medium work sat behind and
turned a seconds-long check into a two-minute one. The medium and the rig are
different structures (T2 inviscid/the-rig-and-the-medium-are-different-structures)
and after this file they are different dependency trees as well.

WHAT BELONGS HERE: geometry that is phase-independent or purely per-unit, used
by the honeycomb work. Nothing about arrays, contacts, cranks, wires or
topologies. If something here starts needing jb_z, it is in the wrong file.

ONE jb_z REFERENCE SURVIVES ON PURPOSE, and it is not this file's: jb_w's W2 row
imports `fold_halves` locally to check lam(a) against an INDEPENDENT
implementation. A cross-check that reaches outside its own tree is doing its job;
copying the function across would make it a self-check. That import is local to
the row, so the medium's runtime path stays free of jb_z.
"""

from __future__ import annotations

from analysis.model.jitterbug import corners, faces  # noqa: F401  (corners is re-exported)

#: The 8 octahedron faces as (v, c, u, sigma). PHASE INDEPENDENT -- jb_hc's Z0
#: row gates that a plate's outward normal does not move with the fold, and the
#: whole contact kernel rests on it.
_FACES = faces()


def plate_normal(face_idx):
    """The FIXED outward unit normal of plate `face_idx`, independent of phase."""
    return _FACES[face_idx][2]
