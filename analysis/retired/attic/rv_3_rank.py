"""REVIEW CHECK 3: rank determination robustness + completeness of the
singularity list + independent count of the tetrahedron combinatorics.
"""
import itertools
import numpy as np
from jb_a_family import corners
from jb_b_variety import jacobian, PAIRS
from jb_d_tet import build_tet_config
from jb_c_branches import tet_assignments

np.set_printoptions(precision=6, suppress=True)

print("=== A. full singular tail; is the rank cut anywhere near a real value? ===")
for a in (0.0, 30.0, 60.0, 75.0, 90.0):
    s = np.linalg.svd(jacobian(corners(a)), compute_uv=False)
    print(f"  a={a:6.2f}  s[0]={s[0]:.4f}  s[33..35]={s[33]:.4e},{s[34]:.4e},{s[35]:.4e}")

print("\n=== B. rank as a function of the tolerance (7 decades) ===")
tets = [s for s in tet_assignments()
        if len(set(tuple(sorted(d.values())) for d in s)) == 4]
Xtet = build_tet_config(tets[0])
cases = [("a=0", corners(0.0)), ("a=60", corners(60.0)), ("a=59.999999", corners(59.999999)),
         ("a=90", corners(90.0)), ("a=90+1e-6", corners(90.000001)), ("tetrahedron", Xtet)]
print(f"  {'config':14s} " + " ".join(f"{t:>8s}" for t in
      ("1e-14", "1e-12", "1e-10", "1e-9", "1e-8", "1e-6", "1e-4")))
for name, X in cases:
    s = np.linalg.svd(jacobian(X), compute_uv=False)
    row = [int(np.sum(s > t * max(1.0, s[0]))) for t in
           (1e-14, 1e-12, 1e-10, 1e-9, 1e-8, 1e-6, 1e-4)]
    print(f"  {name:14s} " + " ".join(f"{r:8d}" for r in row))

print("\n=== C. FULL sweep for rank drops: min sigma_36 over a in [0,360) ===")
A = np.linspace(0.0, 360.0, 3601)
S36 = np.array([np.linalg.svd(jacobian(corners(a)), compute_uv=False)[35] for a in A])
print(f"  global min sigma_36 = {S36.min():.4e} at a={A[np.argmin(S36)]:.3f}")
# local minima below a generous threshold
loc = [(A[i], S36[i]) for i in range(1, len(A) - 1)
       if S36[i] < S36[i-1] and S36[i] < S36[i+1] and S36[i] < 1e-2]
print(f"  local minima with sigma_36 < 1e-2: {[(round(a,3), f'{v:.2e}') for a, v in loc]}")
print(f"  sigma_36 at a=0,60,120,180,240,300: "
      f"{[round(float(np.linalg.svd(jacobian(corners(a)), compute_uv=False)[35]),4) for a in (0,60,120,180,240,300)]}")
print("  (coarse sweep step 0.1 deg; the a=90 zero is transversal with slope ~0.03/deg,")
print("   so any additional transversal zero would show as a dip well below 1e-2)")

print("\n=== D. is the a=90 rank drop a real extra MECHANISM or a parameterisation artifact? ===")
X90 = corners(90.0)
J90 = jacobian(X90)
U, S, Vt = np.linalg.svd(J90)
null = Vt[35:]                              # 13 null directions
cen = X90.mean(axis=1)
print(f"  all 8 centroids at the origin? max|cen| = {np.abs(cen).max():.3e}")
glob = np.zeros((6, 48))
for d in range(3):
    for i in range(8):
        glob[d, 24 + 3*i + d] = 1.0
for d in range(3):
    e = np.eye(3)[d]
    for i in range(8):
        glob[3 + d, 3*i:3*i+3] = e
        glob[3 + d, 24 + 3*i:24 + 3*i + 3] = np.cross(e, cen[i])
print(f"  rank of the 6 global vectors at a=90 = {np.linalg.matrix_rank(glob)}")
Q, _ = np.linalg.qr(glob.T)
internal = null - (null @ Q) @ Q.T
Si = np.linalg.svd(internal.T, compute_uv=False)
print(f"  internal directions after removing globals: {int(np.sum(Si > 1e-8))}  (a=90)")

print("\n=== E. independent recount of the tetrahedron combinatorics ===")
# Reformulate: colour each of the 12 SHARED VERTICES with a tet vertex 0..3;
# a face is admissible iff its 3 corner colours are distinct.  This is an
# algorithmically different enumeration from jb_c_branches' per-face backtracking.
vert_of = {}
for v, ((i, j), (k, l)) in enumerate(PAIRS):
    vert_of[(i, j)] = v
    vert_of[(k, l)] = v
face_verts = [[vert_of[(i, j)] for j in range(3)] for i in range(8)]
print(f"  each face's 3 shared-vertex ids: {face_verts}")
total, profiles = 0, {}
for col in itertools.product(range(4), repeat=12):
    tri = []
    for fv in face_verts:
        t = (col[fv[0]], col[fv[1]], col[fv[2]])
        if len(set(t)) != 3:
            break
        tri.append(tuple(sorted(t)))
    else:
        total += 1
        load = {}
        for t in tri:
            load[t] = load.get(t, 0) + 1
        key = tuple(sorted(load.values(), reverse=True))
        profiles[key] = profiles.get(key, 0) + 1
print(f"  total hinge-consistent colourings = {total}   (jb_c reports 9216)")
for k in sorted(profiles, key=lambda t: (-len(t), t)):
    print(f"    profile {str(k):14s} -> {profiles[k]:6d}"
          f"{'   <== tetrahedron' if k == (2,2,2,2) else ''}")
