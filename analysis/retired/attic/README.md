# attic — the single-unit variety chain

Twenty modules, ~612 KB, moved here 2026-08-26. **Nothing is deleted; git has
the full history and `git log --follow` works through the move.**

## Why they are here

This project *landed* on its current model rather than starting from it, and the
build carries the path dependence. These modules are the single-unit **variety**
study — the jitterbug's configuration manifold, its curvature, the Riemannian
Hessian, the mode structure, the frequency spectrum. Real work, and the source
of results the canonical record still cites.

They are here because **nothing on the medium path imports any of them.** They
were not in the verification loop, so they were already unmaintained in
practice; keeping them in the live directory cost nothing in CPU and a great
deal in legibility. Moving them makes the live surface answer the question "what
is actually load-bearing?" without anyone having to trace imports to find out.

## What is still cited from here

The canonical entry (T2 `inviscid/qvf-epic-consolidated-state`, section (c)(3))
lists the single-unit and symmetric-path results that survived the topology
correction. Those remain true and quotable:

- the **vertex ellipses** — semi-major √2, semi-minor √(2/3), axis ratio exactly
  √3, the crank angle as eccentric anomaly. `jb_hc_honeycomb` recomputes these
  itself and gates them, so the live path does not depend on the attic for them.
- **gap(a) = 4 cos(a + 60°)** and **a\*(t) = arccos(t/4) − 60°**
- the **arm-A k-table** and the **arm-A / arm-C crossover at k_c = √(2/3)**

If you need to re-derive any of those, the module that produced them is here.

## What is NOT here

`jb_a_family`, `jb_b_variety`, `jb_g_strut_clearance`, `jb_cache` and
`jb_x_array_linkage` stayed in the live directory: the medium work or `jb_z`
still import them.

`jb_b_variety` is here-by-necessity rather than by category — it belongs to this
chain, but `jb_g_strut_clearance` takes `PAIRS` from it, `jb_z` takes
`segment_distance` from `jb_g`, and `jb_w`'s W2 cross-check imports `jb_z`. That
path runs through a *local* import inside a function, so a static scan from the
medium modules does not see it; the archive was moved, jb_w went red, and the
chain was read off the traceback. Worth knowing before moving anything else.

`jb_y_dephasing` and `jb_z_quasistatic_array` also stayed. They are the **rig /
bench** line rather than the medium line, and their array results are computed on
the single-vertex-contact topology `inviscid-ia5` retired — but their *machinery*
was corrected on 2026-08-26 (the wire span, `inviscid-l1d`, and the per-unit
phase threading), so they are the right tool if bench predictions are wanted
again. See T2 `inviscid/the-rig-and-the-medium-are-different-structures`.

## Running one

They import `jb_a_family` and friends from the parent directory, so put it on
the path:

```
cd analysis/jitterbug-variety
PYTHONPATH=. python3 -W ignore attic/jb_u_riemannian_hessian.py
```
