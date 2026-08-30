"""Is the SUCCESSFUL VE->tetrahedron path collision-free? The variety has no
non-interpenetration term, so 'reachable' may mean 'reachable through states where
the rigid plates pass through each other' -- not a motion a mechanism can perform."""
import numpy as np, itertools
from jb_a_family import corners
from jb_c_branches import tet_assignments, place, residual
from jb_d_tet import build_tet_config, kabsch
from jb_e_tighten import is_tet, project
from scipy.optimize import least_squares
from crit_collide2 import penetration as _pen

def pen_X(X, tol=1e-6):
    import crit_collide2 as cc
    out=0
    for i,j in itertools.combinations(range(8),2):
        shared = sum(1 for p in X[i] for q in X[j] if np.linalg.norm(p-q)<1e-7)
        if shared: continue
        for (A,B) in ((i,j),(j,i)):
            for k in range(3):
                if cc.edge_hits_tri(X[A][k], X[A][(k+1)%3], X[B], tol) is not None: out+=1
    return out

def walk_traced(X0, Xt, n=80, pull=1.0, w=1e4):
    z=np.zeros(48); trace=[]
    for s in np.linspace(1.0/n,1.0,n):
        Xa = kabsch(Xt.reshape(-1,3), place(X0,z).reshape(-1,3)).reshape(8,3,3)
        way=(1-s)*place(X0,z)+s*Xa
        sol=least_squares(lambda zz: np.concatenate([w*residual(place(X0,zz)),
             pull*(place(X0,zz)-way).reshape(-1)]), z, xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        zp,_=project(X0,sol.x); z=zp; trace.append(place(X0,z))
    return trace

tets=[s for s in tet_assignments() if is_tet(s)]
X0=corners(0.0)
tr=walk_traced(X0, build_tet_config(tets[1]), n=80)
pens=[pen_X(X) for X in tr]
print("collision profile along the SUCCESSFUL VE -> tetrahedron path (target t1, 80 waypoints)")
print("  waypoint : strict-interior edge-through-face crossings")
for k in range(0,80,5):
    print(f"   {k:3d} : {pens[k]}")
print(f"   79 : {pens[79]}  (endpoint = the tetrahedron)")
nz=[k for k,v in enumerate(pens) if v>0]
print(f"\n  waypoints with self-intersection: {len(nz)}/80   first at {nz[0] if nz else '-'}")
print(f"  max crossings on the path: {max(pens)}")
