# Triangle-network probes (2026-08-30)

Owner decision 2026-08-30 (T2 `inviscid`, DECISION 22): the hinge between two
triangles inside one jitterbug and the tie between two jitterbugs are the same
joint. Then "the cell" is a modelling construct and the medium is a periodic
network of rigid triangles with uniform compliant four-corner junctions.

These are PROBES, not gated modules. They are committed so their numbers stay
traceable (the 2026-08-29 handoff's lesson). Run with
`/opt/homebrew/opt/python@3.13/libexec/bin/python3 -W ignore <file>`; the
export scripts take the scratch directory as their argument and write
`frames.json` there, which `patch-template.html` embeds at `__DATA__`.

- `tri_network.py` — the network builder (every triangle a rigid body, every
  junction a blob); rigid-junction mechanism counts and soft spectra on blocks.
- `tri_bloch.py` — the periodic 48-band calculation, sound speeds by direction,
  and where a finite block's mechanisms live.
- `export_frames.py`, `export_kick2.py` — the cell-model breathe and kick
  frames (jb_rc.honeycomb_single, jb_je.state at n = 2) for the viewer.
- `export_net.py` — the same kick on the network.
- `patch-template.html` — the three-tab viewer (Breathe / Kick / Kick, network),
  published as the "Patch of Record" artifact.
