"""Critic check: does the symmetric path pass through SELF-INTERSECTING configurations?
At a=90 every centroid is at Z*cos(90)=0, i.e. all 8 triangles are centred on the origin.
The constraint model has no non-interpenetration term, so the 'full circuit' the hbar
periodicity argument and the period formula both rely on may pass through total collision.
"""
import numpy as np, itertools
from jb_a_family import corners

def tri_tri_intersect(A, B, eps=1e-9):
    """Separating-axis test for two triangles in 3D (coplanar cases ignored -> treated as hit)."""
    def axes(P, Q):
        n1 = np.cross(P[1]-P[0], P[2]-P[0]); n2 = np.cross(Q[1]-Q[0], Q[2]-Q[0])
        ax = [n1, n2]
        for i in range(3):
            for j in range(3):
                c = np.cross(P[(i+1)%3]-P[i], Q[(j+1)%3]-Q[j])
                ax.append(c)
        return ax
    for ax in axes(A, B):
        n = np.linalg.norm(ax)
        if n < eps: continue
        ax = ax/n
        a0,a1 = (A@ax).min(), (A@ax).max()
        b0,b1 = (B@ax).min(), (B@ax).max()
        if a1 < b0 - eps or b1 < a0 - eps:
            return False          # separated on this axis
    return True                   # no separating axis found -> intersecting

def count_hits(a):
    X = corners(a)
    hits = 0
    for i, j in itertools.combinations(range(8), 2):
        # faces sharing a vertex always "touch"; only count genuine overlap
        shared = sum(1 for p in X[i] for q in X[j] if np.linalg.norm(p-q) < 1e-7)
        if tri_tri_intersect(X[i], X[j]) and shared == 0:
            hits += 1
    return hits

print(" a(deg)  intersecting non-adjacent face pairs (of 28 total pairs)")
for a in [0,10,20,22.24,30,40,50,55,60,62,65,70,75,80,85,89,90,91,95,100,110,120,150,180]:
    print(f"  {a:7.2f}   {count_hits(a)}")
# bisect the onset
lo, hi = 60.0, 90.0
while hi-lo > 1e-6:
    mid=(lo+hi)/2
    if count_hits(mid)==0: lo=mid
    else: hi=mid
print(f"\nfirst self-intersection on the symmetric path at a = {hi:.6f} deg")
