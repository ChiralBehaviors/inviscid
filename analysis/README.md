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
