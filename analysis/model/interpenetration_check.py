"""Independent verification of the critic's two criticals.

Written from scratch rather than by running the critic's scripts: the project's
standing rule is reproduce, do not accept.

A. Which quantity does the free jitterbug conserve, p = M_eff*adot or E?
B. Does the model self-intersect, and if so where?
"""
import numpy as np
from scipy.integrate import solve_ivp
from analysis.model.jitterbug import corners
from analysis.model.linkage_variety import PAIRS

Z2, RHO2 = 2.0 / 3.0, 1.0 / 3.0          # from the spec: M_eff/M = Z^2 sin^2 a + rho^2


def f(a):
    return Z2 * np.sin(a) ** 2 + RHO2


def fp(a):
    return Z2 * 2 * np.sin(a) * np.cos(a)


# ----------------------------------------------------------------- A. dynamics
def eom(t, y):
    """Euler-Lagrange for L = 0.5 f(a) adot^2 (M L^2 = 1, no potential):
       d/dt(f adot) = 0.5 f' adot^2  =>  addot = -0.5 (f'/f) adot^2"""
    a, ad = y
    return [ad, -0.5 * fp(a) / f(a) * ad ** 2]


def run_dynamics():
    print("=== A. WHAT DOES THE FREE JITTERBUG CONSERVE? ===")
    print("   L = 0.5 M_eff(a) adot^2, no potential.")
    print(f"   dL/da = 0.5 M_eff'(a) adot^2, which is NONZERO => a is NOT cyclic => p NOT conserved.")
    sol = solve_ivp(eom, [0, 60], [0.0, 1.0], rtol=1e-12, atol=1e-14, dense_output=True)
    a, ad = sol.y
    p = f(a) * ad
    E = 0.5 * f(a) * ad ** 2
    print(f"   integrated to t=60, rtol 1e-12")
    print(f"   p  = M_eff*adot : min {p.min():.10f}  max {p.max():.10f}  "
          f"spread {p.max()/p.min():.7f}   -> {'CONSERVED' if p.max()/p.min()<1+1e-6 else 'NOT conserved'}")
    print(f"   E               : min {E.min():.10f}  max {E.max():.10f}  "
          f"spread {E.max()/E.min():.7f}   -> {'CONSERVED' if E.max()/E.min()<1+1e-6 else 'NOT conserved'}")
    print(f"   adot            : min {ad.min():.10f}  max {ad.max():.10f}  "
          f"spread {ad.max()/ad.min():.7f}   (spec says 3; sqrt(3) = {np.sqrt(3):.7f})")
    print(f"   M_eff           : spread {f(a).max()/f(a).min():.7f}  (3:1 claim SURVIVES)")
    i = np.argmax(ad)
    print(f"   adot peaks at a = {np.degrees(a[i]) % 360:.4f} deg  "
          f"(VE is a=0 mod 180 -> claim SURVIVES)")

    # period, both hypotheses, matched at the SAME initial condition a=0, adot=1
    grid = np.linspace(0, 2 * np.pi, 200001)
    trap = getattr(np, "trapezoid", None) or np.trapz
    Tp = trap(f(grid) / f(0.0), grid)          # p conserved: T = int M_eff da / p
    TE = trap(np.sqrt(f(grid) / f(0.0)), grid)  # E conserved: T = int da / adot
    print(f"   circuit period at a=0,adot=1:  if p conserved (WRONG) {Tp/1.0:.7f}"
          f"   truth, E conserved {TE:.7f}   ratio {Tp/TE:.7f}")
    print("   NOTE: period VALUES are normalisation-dependent; the sqrt(3) spread is not.")


# --------------------------------------------------------------- B. collisions
def tri_edge_hits(P, Q, shared, eps=1e-9):
    """Moller-Trumbore: does any edge of triangle P pierce triangle Q's interior?
    Intersections at a shared hinge vertex are contact, not interpenetration."""
    hits = 0
    q0, q1, q2 = Q
    e1, e2 = q1 - q0, q2 - q0
    n = np.cross(e1, e2)
    if np.linalg.norm(n) < 1e-12:
        return 0
    for k in range(3):
        o, d = P[k], P[(k + 1) % 3] - P[k]
        pv = np.cross(d, e2)
        det = e1 @ pv
        if abs(det) < 1e-12:
            continue                      # parallel; coplanar handled as no-hit
        inv = 1.0 / det
        tv = o - q0
        u = (tv @ pv) * inv
        if u < eps or u > 1 - eps:
            continue
        qv = np.cross(tv, e1)
        v = (d @ qv) * inv
        if v < eps or u + v > 1 - eps:
            continue
        t = (e2 @ qv) * inv
        if t < eps or t > 1 - eps:
            continue
        X = o + t * d
        if any(np.linalg.norm(X - s) < 1e-7 for s in shared):
            continue                      # touching at a hinge is legal
        hits += 1
    return hits


def interpenetrations(X, eps=1e-9):
    """Count edge-through-face events over all triangle pairs."""
    total = 0
    for i in range(8):
        for j in range(8):
            if i == j:
                continue
            shared = [X[i, c] for c in range(3)
                      for ((p, cc), (q, dd)) in PAIRS
                      if (p, cc) == (i, c) and q == j or (q, dd) == (i, c) and p == j]
            total += tri_edge_hits(X[i], X[j], shared, eps)
    return total


def run_collisions():
    print("\n=== B. DOES THE MODEL SELF-INTERSECT? ===")
    print("   edge-through-face events, hinge contacts excluded")
    for a in (0, 15, 30, 45, 55, 59, 59.9, 60, 60.1, 61, 70, 80, 90,
              100, 110, 119, 119.9, 120, 120.1, 130, 150, 180):
        print(f"   a={a:7.2f}  interpenetrations = {interpenetrations(corners(a)):3d}")


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    run_dynamics()
    run_collisions()
