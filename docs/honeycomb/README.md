# Honeycomb figures

Interactive figures for the VE/octa rectified cubic honeycomb. Each is a
self-contained HTML page — open it in a browser, no build step, no network
except Google Fonts. Everything they draw is computed in the page from the
closed form `p(a) = √2·cos a·m + √(2/3)·sin a·b`, verified against the real
`verts(a)` to 4.4e-16.

They were published as Claude Artifacts during the 2026-08-25 session and are
copied here because at that point they were the *only* durable record of the
corrected exchange — the T2 entry still carried a retracted version, and no
code existed. The code now lives in `../analysis/jitterbug-variety/jb_hc_honeycomb.py`.

| file | what it shows |
|---|---|
| `orbits.html` | one face on its three vertex ellipses; the exact constants; the 30° first-contact wall on the **symmetric path** |
| `edge.html` | two VEs at a shared vertex, four collinear strut pairs — **the wrong packing**, kept as the record of a corrected mistake |
| `exchange.html` | one VE and one hole cell: the exchange, `b = a+60`, driven to the full swap at `a = −60°` |
| `three.html` | VE — hole — VE: the reciprocal condition, and the square contact decaying 4 → 2 → 1 |
| `eight.html` | eight VEs around one octahedral cell; the whole cluster turning inside out |

`edge.html` documents a superseded configuration. Neighbours in this packing
share **faces**, not a single vertex; that page shows the vertex-touching
arrangement the codebase's `build_topologies` builds, which is not how these
cells pack. It is kept deliberately — the retraction history is part of the
record.

Canonical prose state: T2 `inviscid/qvf-epic-consolidated-state`.
