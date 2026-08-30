"""REVIEW CHECK 7: pairing() combinatorics — correct, stable in a, not re-derived.
Plus: is "12 shared vertices at every a except 60" true over a full sweep?
"""
import numpy as np
from jb_a_family import corners, cluster
from jb_b_variety import pairing, PAIRS

print("=== A. is the hinge pairing the SAME at every generic a? ===")
base = set(PAIRS)
bad = []
for a in np.concatenate([np.linspace(0.5, 59.5, 60), np.linspace(60.5, 359.5, 60)]):
    try:
        p = set(pairing(a_probe=float(a)))
    except AssertionError:
        bad.append((round(float(a), 2), "not 12x2"))
        continue
    if p != base:
        bad.append((round(float(a), 2), "DIFFERENT pairing"))
print(f"  probed 120 values of a; disagreements: {bad if bad else 'NONE'}")

print("\n=== B. face-adjacency graph implied by the pairing ===")
adj = {i: set() for i in range(8)}
for (i, j), (k, l) in PAIRS:
    adj[i].add(k)
    adj[k].add(i)
deg = sorted(len(v) for v in adj.values())
print(f"  8 triangles, degrees = {deg}  (cuboctahedron: every triangle meets 3 others)")
edges = sorted({tuple(sorted((i, k))) for i in adj for k in adj[i]})
print(f"  {len(edges)} face-face hinges (cube graph has 12): {edges}")
# 3-regular bipartite & triangle-free => cube graph Q3
tri = sum(1 for a_ in range(8) for b in adj[a_] for c in adj[a_] if b < c and c in adj[b])
print(f"  triangles in the face-adjacency graph: {tri}  (Q3 is triangle-free => 0)")

print("\n=== C. distinct-corner count over a DENSE sweep (claim: 12 everywhere but 60) ===")
odd = []
for a in np.linspace(0.0, 360.0, 721):
    reps, mult, _ = cluster(corners(float(a)), tol=1e-7)
    if len(reps) != 12:
        odd.append((round(float(a), 2), len(reps), sorted(set(mult.tolist()))))
print(f"  values of a (0.5 deg grid) where the count is NOT 12: {odd}")

print("\n=== D. is PAIRS re-derived anywhere downstream? ===")
import subprocess
out = subprocess.run(["grep", "-n", "pairing(\\|PAIRS",
                      "jb_a_family.py", "jb_b_variety.py", "jb_c_branches.py",
                      "jb_d_tet.py", "jb_e_tighten.py"],
                     capture_output=True, text=True).stdout
print(out or "  (no matches)")
