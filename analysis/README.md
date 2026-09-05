# analysis — the map

Everything here is a self-gating Python script: run it, read the `PASS` /
`FAIL` rows, the exit code is the verdict. Run from the repo root by module
name, which is what makes the imports resolve:

    analysis/gates.sh                    # every module, live first
    analysis/gates.sh live               # the model of record only
    python -m analysis.model.dispersion  # one module

Interpreter: `/opt/homebrew/opt/python@3.13/libexec/bin/python3` (needs numpy,
scipy). A bare `python3` in a nohup'd shell has no numpy and every gate then
reports nothing.

Reorganised 2026-08-30 from a flat directory of two-letter codes. The old code
is in brackets after every name because the record (T2, beads, commit
messages before this date) cites them. `git log --follow` works through the
moves; nothing was deleted.

## model/ — the model of record

Rigid struts. Soft rubber joints at the vertices. The octahedral voids empty
(one cell per all-even site, six axis neighbours, two tied vertices per
neighbour). One fold angle per cell.

| module | what it is | old |
|---|---|---|
| `jitterbug.py` | one cell: the symmetric jitterbug family, faces, body frame | jb_a_family |
| `plates.py` | the rigid triangular plates | jb_gp_plate_geometry |
| `cell.py` | vertex identities, strut length √2, `cell_verts` | jb_ic_inertial_chain |
| `cluster.py` | which face faces which neighbour | jb_cl_cluster |
| `assembly.py` | cells + welds in reduced coordinates; `honeycomb` (both sublattices, the old covering) and `honeycomb_single` (voids empty); lattice constant; the coherent breathe | jb_rc_reduced |
| `kinematics.py` | tied pairs, mass matrices, free acceleration, separations, band rows, the standard kicks — split out of the impact-law module so nothing live imports a retired file | (from jb_mj) |
| `dispersion.py` | the 7 bands, two sound speeds in the ratio √2, the zone corner | jb_1c_single_covering |
| `block_spectrum.py` | finite blocks against the bands: seven free modes at every size, the top under √3, the bottom acoustic | jb_1s_single_spectrum |
| `joint_exponent.py` | the speed–amplitude exponent family; the Hertz check | jb_1e_single_exponent |
| `linkage_variety.py` | the jitterbug as a linkage, early facts | jb_b_variety |
| `strut_clearance.py` | the 60–120 interference without faces | jb_g_strut_clearance |
| `interpenetration_check.py` | independent edge-through-face count that `strut_clearance` cross-checks against | verify_critic |
| `kick_response.py` | one cell driven, the medium's reaction: a scalar fold wave at speed 1/2 | — |
| `shear_response.py` | one cell shoved sideways or stirred: polarized translational waves; the longitudinal one degenerate with the fold at 1/2 | — |
| `soliton.py` | the plane-symmetric medium as an exact two-field lattice (strain, fold): both branches ω = sin(q/2), the strain exactly harmonic, the fold's on-site quartic; the envelope soliton and its collision | — |
| `envelope.py` | rung 4: the envelope equation with every coefficient the chain's own (P, Q_ring, Q_dent, Q_x); where it stops toward long carriers | — |
| `longwave.py` | rung 4b, three fields: envelope + strain field + the gapped mean fold; the umklapp exchange, the cold launch, the dent's schedule, the gap K_d(q), the π/4 instability at a 6 | — |
| `instability.py` | rung 4b, five fields: + the second-harmonic envelope C (the plane wave's four-wave decay) and the strain envelope U at the carrier (the packet's strain-mediated growth, resonant because the two branches are degenerate), + the on-site quartic; the rates at a 1.5–6 and π/3, the π/2 null, live strain's damping, buckling | — |

### model/double_covering/

The same soft-joint line as it was measured before the voids were emptied,
with a cell on every site. Kept because the live modules reproduce its
numbers as their instrument check, and because it is where the results were
first found. Its absolute numbers are properties of that covering, not of
the medium.

| module | old |
|---|---|
| `dispersion.py` — 14 bands | jb_bz_dispersion |
| `soft_joint_spectrum.py` — the soft joint clears qvf.2 | jb_sj_soft_joint |
| `joint_exponent.py` — the exponent family, first measurement | jb_je_joint_exponent |
| `asymmetric_joint.py` — an asymmetric joint costs the Hessian | jb_ja_asymmetric_joint |

### model/first_principles/

The tied array built up one body at a time, 2026-08-30, each step a gated
measurement and each with a page the owner looked at. The record of the
conversation is in T2 [23789] [23791] [23794]-[23796]; the physics in one
paragraph: one cell passes through the VE and chooses a sense there; an
octahedron with its joints as they are opens one way; an octahedron hung on a
VE face is a passenger; a void bounded by plates of different VEs in a closed
ring is a VE that expands, its fold the VE's plus sixty; so the tied array's
motion is one sixty-degree segment with an octahedron at each end. The
overdrive continues the parametrisation past both ends: joints hold by
identity, plates turn inside out, and any patch collapses onto one octahedron
at +60 and +240.

| module | what it is |
|---|---|
| `geometry.py` | the shared measuring tools: coincidence classes, the twelve joints, segment distance (scalar and vectorised crossings), a rigid fit, front signs |
| `one_cell.py` | one jitterbug −60 → +60: joints permanent, six vertices at ±60, the two octahedra pair the joints differently, the octahedron a dead end even as a free linkage, the VE the choice |
| `face_to_face.py` | a VE and an octahedron cell on one shared plate: two independent folds, turn = a − b + 60, six corners per shared joint, never touching; two passengers stay passengers |
| `vertex_point.py` | one vertex at the VE: four cells, eight corners, two permanent joints; the tied block walks through the VE on its weld manifold and is stopped by the two joints passing head-on |
| `ring.py` | the smallest closed ring: one motion, void = VE + 60, no rotation, the void expands, both ends octahedra |
| `overdrive.py` | the ring, HC15 and the 3×3×3 block through 360°: welds by identity, the spacing sinusoid, the side flips, crossings only in the passages, the collapse |

`pages/` holds the exporters and templates behind the pages ("One Jitterbug",
"Face to Face", "One Vertex, Four Cells", "Two Joints, One Point", "The
Exchange", "A Full Cycle", "Overdrive" ×3). `python -m
analysis.model.first_principles.pages.export_ring` and friends write frames to
`analysis/.pages/data/`; `python -m analysis.model.first_principles.pages.build`
inlines them into self-contained HTML in `analysis/.pages/` (git-ignored).

### model/su2/

The SU(2)/boundary-conditions gate lane from
`notes/su2_boundary_conditions.md` §8.6 (read §8 before citing §1–7).
Vocabulary is classical ℤ₂ holonomy / antiperiodic Bloch sector — no spin,
no quantum number.

| module | what it is |
|---|---|
| `lift.py` | G0: sign-continuous SU(2) lift of a sampled SO(3) path with its null controls (trivial loop +q, bare 2π −q, bare 4π +q); refuses coarse sampling and unclosed loops |
| `plate_holonomy.py` | G1: every plate of every body on the record's three patches lifts to −q over one 360° overdrive cycle — a confirmation of the record (§8.3), necessary not sufficient for any SU(2)-medium claim |
| `inside_out_cover.py` | G4: X(a+180) = −X(a) exactly on configurations, const = 0, inversion centre the lattice origin — deck-map candidate (b) pinned; sheet-touching in ambient space is patch geometry (site-symmetric patches coincide everywhere, uncentred ones only at the collapses) |
| `joint_twist.py` | G2 (subsumes G3): the relative-rotation lift across BOTH co-located vertex joints ({O,A}, {B,D}) is +q — the ℤ₂ cancels pairwise; it lives in absolute orientation histories only, so even a tether-capable joint would return untwisted each cycle; face-weld relative history is constant |
| `physical_arc.py` | G7: the physical 60° arc never exercises the holonomy — crossing-free exactly on [−60, 0] with walls within 1° of both ends, breathe lifts +q everywhere, oscillation of any amplitude lifts +q (only monotone traversal reaches −q); closed patches are one-freedom, the free block's nullity 15 identified as breathe + 8 hanging corner-VE folds |
| `screw_prerequisites.py` | τ prerequisites for G6: single applications of τ leave the coherent family (only τ¹ carries ex-cells, τ⁵ ex-voids, τ²⁻⁴ nothing); τ⁶ = pure translation exactly, τ³ = translate ∘ per-body inversion; free + properly discontinuous on the arc (spacing √3·L, L ≥ 1), degenerate at the collapses — G6 must use the mapping-torus identification, Bieberbach naming stays withdrawn |
| `doubled_block.py` | G5: the §3c double built (every half-weld closed, all cells degree 6) — it does NOT lose the seven zero modes: spectrum = free (Neumann) ⊎ pinned-boundary (Dirichlet, no zeros) exactly; torus reference = wrapped operator ≡ Bloch bands at commensurate k (4 zeros); the double runs ~3× closer to the Bloch DOS than the free block and at the torus's level from side 3, winning at matched cell count |
| `screw_sector.py` | G6, resolved in the negative: odd multiples of (1,1,1) put zero solid sites on solid sites (the screw has nothing to act on under DECISION 21); no isometry pairs cell and void bodies at any generic fold (0/48) — the glide pairing exists only at a = −30 (12+12 of 48) and belongs to the retired both-kinds covering; the BC question lands on G5's race, double vs torus |

## retired/ — lines the owner has closed

Still runnable, still gated, not maintained. Do not quote their numbers as
properties of the medium.

### retired/hard_wall/ — a clearance with a stop, and the impact law

The joint as a ball of play with a hard wall: the n → ∞ end of the
joint-law family, not the rubber joint the rig is made of. Superseded by the
soft-joint decision of 2026-08-28 (DECISION 19 reformulated the criterion
for it; the rubber decision made the original criterion reachable).

| module | old |
|---|---|
| `contact_chain.py` | jb_ct_contact_chain |
| `contact_potential.py` | jb_cp_contact_potential |
| `impact_law.py` — Moreau-Jean on the honeycomb | jb_mj_inertial_honeycomb |
| `transport.py` | jb_tr_transport |
| `sonic_vacuum.py` | jb_sv_sonic_vacuum |
| `lattice_front.py` | jb_lf_lattice_front |
| `sonic_vacuum_single.py` — the same under one covering | jb_1v_single_vacuum |
| `phase_ramp.py` — fixed centres, the rig's clearance budget | jb_pr_phase_ramp |

### retired/rig_lock/ — the qvf.11 lock and the quasi-static crank

The physical rig's one-sided lock, withdrawn as evidence about the medium
(DECISION 20, 2026-08-29): the rig held centres Gray requires to move. These
also fill the voids.

| module | old |
|---|---|
| `array_linkage.py` | jb_x_array_linkage |
| `dephasing.py` | jb_y_dephasing |
| `quasistatic_crank.py` | jb_z_quasistatic_array |
| `honeycomb_exchange.py` | jb_w_honeycomb |
| `honeycomb_topology.py` | jb_ht_honeycomb_topology |
| `held_correspondences.py` | jb_hh_held_correspondences |
| `cache.py` | jb_cache |

### retired/strut_springs/ — the void fork

Springs on the struts: rejected by the owner on 2026-08-27 ("the jitterbug
transformation *is* the wave medium, not vibration modes in the struts").
Kinematic facts in these survive; anything read as a frequency does not.

| module | old |
|---|---|
| `inertial_array.py` | jb_aa_inertial_array |
| `flexible_lattice.py` | jb_fl_flexible_lattice |
| `band_touching.py` | jb_bt_band_touching |
| `calibration.py` | jb_cal_calibration |
| `honeycomb_waves.py` | jb_hc_honeycomb |

### retired/attic/

The single-unit variety study (twenty modules, `jb_c` … `jb_v`, moved
2026-08-26) and the loose review probes (`crit_*`, `rv_*`, `verify_*`,
`DumpJb.java`). Not runnable from here without their old sibling imports;
kept as history, with their own README.

## history/

`derivation-record.md` — the 704-line record that used to be this directory's
README: what was overturned, when, and why. Read it for the path; read the
tables above for the state.
