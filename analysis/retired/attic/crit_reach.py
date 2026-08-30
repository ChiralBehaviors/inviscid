"""Critic check on item 4: the 3 'undetermined' targets.
(a) What IS the endpoint of a failed walk? Is it itself a tetrahedron?
(b) Are the 8 'sampled' targets a meaningful sample?
"""
import numpy as np
from jb_a_family import corners
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import build_tet_config, kabsch
from jb_e_tighten import is_tet, walk

tets = [s for s in tet_assignments() if is_tet(s)]
print(f"total tet assignments: {len(tets)};  jb_e samples tets[0:8] -- enumeration order, not random")
print("first 8 assignments, face->(sorted tet-face) signature:")
for t in range(8):
    print("   t%d: %s" % (t, [ "".join(map(str,sorted(d.values()))) for d in tets[t] ]))

X0 = corners(0.0)
print("\n=== endpoint characterisation of each of the 8 walks ===")
for t in range(8):
    Xt = build_tet_config(tets[t])
    X, dist, wa, ws = walk(X0, Xt, n=80)
    pts = X.reshape(-1,3)
    uniq=[]
    for p in pts:
        if not any(np.linalg.norm(p-q)<1e-7 for q in uniq): uniq.append(p)
    D=[np.linalg.norm(uniq[i]-uniq[j]) for i in range(len(uniq)) for j in range(i+1,len(uniq))]
    tag = "REACHED" if dist<1e-6 else "not reached"
    print(f"  t{t}: RMS={dist:.3e} {tag:12s} | endpoint distinct corners={len(uniq):2d} "
          f"| pairwise {min(D):.6f}..{max(D):.6f} | hinge res {np.linalg.norm(residual(X)):.1e}")
