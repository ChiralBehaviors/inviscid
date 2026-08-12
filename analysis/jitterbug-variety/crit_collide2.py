"""Independent confirmation by a different method: Moller-Trumbore edge-vs-triangle
with STRICT interior tolerance. Reports actual penetration, not mere contact."""
import numpy as np, itertools
from jb_a_family import corners

def edge_hits_tri(p0, p1, T, tol=1e-6):
    """Does segment p0->p1 cross the STRICT interior of triangle T?"""
    e1, e2 = T[1]-T[0], T[2]-T[0]
    d = p1-p0
    h = np.cross(d, e2); a = e1 @ h
    if abs(a) < 1e-12: return None
    f = 1.0/a; s = p0 - T[0]
    u = f*(s@h)
    q = np.cross(s, e1)
    v = f*(d@q)
    t = f*(e2@q)
    if u > tol and v > tol and u+v < 1-tol and t > tol and t < 1-tol:
        return t
    return None

def penetration(a):
    X = corners(a); out=[]
    for i, j in itertools.combinations(range(8), 2):
        shared = sum(1 for p in X[i] for q in X[j] if np.linalg.norm(p-q) < 1e-7)
        if shared: continue
        for (A,B) in ((i,j),(j,i)):
            for k in range(3):
                t = edge_hits_tri(X[A][k], X[A][(k+1)%3], X[B])
                if t is not None:
                    out.append((i,j,t))
    return out

for a in [55,59,59.9,60,60.1,61,62,70,90,110,119,119.9,120,121,125]:
    h = penetration(a)
    depth = f"min param {min(x[2] for x in h):.4f}" if h else "-"
    print(f"  a={a:7.2f}  strict interior edge-through-face crossings: {len(h):3d}   {depth}")

lo,hi=59.0,62.0
while hi-lo>1e-7:
    m=(lo+hi)/2
    if len(penetration(m))==0: lo=m
    else: hi=m
print(f"\n  collision onset (strict interior): a = {hi:.7f} deg")
lo,hi=118.0,121.0
while hi-lo>1e-7:
    m=(lo+hi)/2
    if len(penetration(m))>0: lo=m
    else: hi=m
print(f"  collision release (strict interior): a = {hi:.7f} deg")
