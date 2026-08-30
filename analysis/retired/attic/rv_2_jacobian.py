"""REVIEW CHECK 2: verify jacobian() against a finite difference of the EXACT
nonlinear constraint map used by the continuation code.

jb_c_branches.place()/unpack() is the genuine nonlinear map:
    X(z) = [ exp([w_i]x) (X0_i - cen_i) + cen_i + t_i ]
    residual(X) = concat over the 12 hinge pairs of (X[i,j] - X[k,l])

So d residual / d z at z = 0 IS the object jb_b_variety.jacobian(X0) claims to be.
Compare them entry by entry.  This settles the sign/convention question outright.
"""
import numpy as np
from jb_a_family import corners
from jb_b_variety import jacobian, PAIRS, cross_row
from jb_c_branches import place, residual


def fd_jacobian(X0, h=1e-6):
    z0 = np.zeros(48)
    J = np.zeros((36, 48))
    for k in range(48):
        zp, zm = z0.copy(), z0.copy()
        zp[k] += h
        zm[k] -= h
        J[:, k] = (residual(place(X0, zp)) - residual(place(X0, zm))) / (2 * h)
    return J


print("=== A. analytic J vs central-difference of residual(place(X0, z)) at z=0 ===")
for a in (0.0, 22.238756093, 30.0, 60.0, 90.0):
    X0 = corners(a)
    Ja, Jf = jacobian(X0), fd_jacobian(X0)
    print(f"  a={a:12.6f}   max|Ja-Jf| = {np.abs(Ja - Jf).max():.3e}"
          f"   ||Ja-Jf||_F/||Ja||_F = {np.linalg.norm(Ja-Jf)/np.linalg.norm(Ja):.3e}")

print("\n=== A2. sign-error control: does the WRONG omega sign show up in this test? ===")


def jacobian_wrongsign(X):
    cen = X.mean(axis=1)
    J = np.zeros((36, 48))
    for c, ((i, j), (k, l)) in enumerate(PAIRS):
        rij, rkl = X[i, j] - cen[i], X[k, l] - cen[k]
        for d in range(3):
            row = 3 * c + d
            J[row, 3*i:3*i+3] -= cross_row(rij, d)      # deliberately flipped
            J[row, 3*k:3*k+3] += cross_row(rkl, d)
            J[row, 24+3*i+d] += 1.0
            J[row, 24+3*k+d] -= 1.0
    return J


X0 = corners(30.0)
Jw, Jf = jacobian_wrongsign(X0), fd_jacobian(X0)
sw = np.linalg.svd(Jw, compute_uv=False)
print(f"  wrong-sign J: max|Jw-Jf| = {np.abs(Jw-Jf).max():.3e}  (test IS sensitive)")
print(f"  wrong-sign J rank at a=30 = {int(np.sum(sw > 1e-9*max(1,sw[0]))):d}"
      f"  <- would have given internal DOF {48-int(np.sum(sw>1e-9*max(1,sw[0])))-6}")

print("\n=== B. cross_row identity:  cross_row(r,d) . w  ==  (w x r)_d ? ===")
rng = np.random.default_rng(0)
err = 0.0
for _ in range(200):
    r, w = rng.normal(size=3), rng.normal(size=3)
    for d in range(3):
        err = max(err, abs(cross_row(r, d) @ w - np.cross(w, r)[d]))
print(f"  max error over 200 random (r,w) = {err:.3e}")

print("\n=== C. does place()/unpack() agree with the (dw,dt) linearisation? ===")
X0 = corners(30.0)
cen = X0.mean(axis=1)
for eps in (1e-3, 1e-5, 1e-7):
    z = rng.normal(size=48) * eps
    Xlin = np.array([X0[i] + np.array([np.cross(z[3*i:3*i+3], X0[i, j] - cen[i])
                                       + z[24+3*i:24+3*i+3] for j in range(3)])
                     for i in range(8)])
    print(f"  eps={eps:.0e}  ||place - linearisation|| = "
          f"{np.linalg.norm(place(X0, z) - Xlin):.3e}   (should be O(eps^2))")
