#!/bin/zsh
# Run every gated analysis module by its module name and count PASS/FAIL rows.
#
#   analysis/gates.sh                 every module, live first, then retired
#   analysis/gates.sh live            the model of record only
#   analysis/gates.sh model.dispersion retired.hard_wall.transport   named modules
#
# Each module is a self-gating script: its exit code is the verdict. Output
# goes to $OUT/<module>.txt (default: analysis/.gates/), one file per module,
# so a failure is readable rather than swallowed. Runs from the repo root via
# `python -m`, which is what makes the package imports resolve.
#
# The interpreter is pinned: a bare `python3` in a nohup'd shell resolves to a
# system interpreter without numpy and every gate then reports pass=0 fail=0.
PY=${PY:-/opt/homebrew/opt/python@3.13/libexec/bin/python3}
ROOT=${ROOT:-$(cd "$(dirname "$0")/.." && pwd)}
OUT=${OUT:-$ROOT/analysis/.gates}
mkdir -p "$OUT"

LIVE=(
  model.jitterbug model.plates model.linkage_variety model.strut_clearance model.interpenetration_check model.cell model.cluster
  model.assembly model.dispersion model.block_spectrum model.joint_exponent
  model.double_covering.dispersion model.double_covering.soft_joint_spectrum
  model.double_covering.joint_exponent model.double_covering.asymmetric_joint
  model.first_principles.one_cell model.first_principles.face_to_face model.first_principles.vertex_point
  model.first_principles.ring model.first_principles.overdrive
  model.su2.lift model.su2.plate_holonomy model.su2.inside_out_cover model.su2.joint_twist
  model.su2.physical_arc
)
RETIRED=(
  retired.hard_wall.contact_chain retired.hard_wall.contact_potential retired.hard_wall.impact_law
  retired.hard_wall.transport retired.hard_wall.sonic_vacuum retired.hard_wall.lattice_front
  retired.hard_wall.sonic_vacuum_single retired.hard_wall.phase_ramp
  retired.rig_lock.array_linkage retired.rig_lock.dephasing retired.rig_lock.quasistatic_crank
  retired.rig_lock.honeycomb_exchange retired.rig_lock.honeycomb_topology retired.rig_lock.held_correspondences
  retired.strut_springs.inertial_array retired.strut_springs.flexible_lattice retired.strut_springs.band_touching
  retired.strut_springs.calibration retired.strut_springs.honeycomb_waves
)

case "$1" in
  "")     mods=("${LIVE[@]}" "${RETIRED[@]}") ;;
  live)   mods=("${LIVE[@]}") ;;
  *)      mods=("$@") ;;
esac

tp=0; tf=0; bad=0
for m in "${mods[@]}"; do
  t0=$(date +%s)
  (cd "$ROOT" && "$PY" -W ignore -m "analysis.$m") > "$OUT/$m.txt" 2>&1
  rc=$?
  t1=$(date +%s)
  p=$(grep -c '^  PASS' "$OUT/$m.txt"); f=$(grep -c '^  FAIL' "$OUT/$m.txt")
  tp=$((tp+p)); tf=$((tf+f)); [ "$rc" -ne 0 ] && bad=$((bad+1))
  printf '%-44s rc=%d pass=%3d fail=%3d %4ds\n' "$m" "$rc" "$p" "$f" "$((t1-t0))"
done
printf 'TOTAL pass=%d fail=%d nonzero-exit=%d\n' "$tp" "$tf" "$bad"
[ "$tf" -eq 0 ] && [ "$bad" -eq 0 ]
