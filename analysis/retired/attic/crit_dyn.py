"""Critic check: for L = 1/2 M_eff(a) adot^2 with M_eff depending on a,
is p = M_eff*adot conserved, or is E = 1/2 M_eff adot^2 conserved?
The spec of record asserts BOTH. Those are incompatible when M_eff varies."""
import numpy as np
from scipy.integrate import solve_ivp

Z2, R2 = 2.0/3.0, 1.0/3.0
f  = lambda a: Z2*np.sin(a)**2 + R2          # M_eff / (M L^2)
fp = lambda a: Z2*2*np.sin(a)*np.cos(a)      # d f / d a

# Euler-Lagrange for L = 1/2 f(a) adot^2:  d/dt(f adot) = 1/2 f'(a) adot^2
#   => f addot + f' adot^2 = 1/2 f' adot^2  => addot = -f'/(2f) * adot^2
def rhs(t, y):
    a, ad = y
    return [ad, -fp(a)/(2*f(a)) * ad**2]

sol = solve_ivp(rhs, [0, 12], [0.0, 1.0], rtol=1e-12, atol=1e-14, dense_output=True, max_step=0.01)
a, ad = sol.y
p = f(a)*ad
E = 0.5*f(a)*ad**2

print("integrating the free jitterbug EOM from a=0, adot=1")
print(f"  p = M_eff*adot :  min {p.min():.10f}  max {p.max():.10f}  spread {p.max()/p.min():.6f}x")
print(f"  E = p^2/(2Meff):  min {E.min():.10f}  max {E.max():.10f}  spread {E.max()/E.min():.6f}x")
print(f"  adot           :  min {ad.min():.10f}  max {ad.max():.10f}  spread {ad.max()/ad.min():.6f}x")
print(f"  f = Meff       :  min {f(a).min():.10f}  max {f(a).max():.10f}  ratio {f(a).max()/f(a).min():.6f}x")
print()
print("predicted if E conserved: adot = sqrt(2E/f), so adot spread = sqrt(3) =", np.sqrt(3))
print("predicted if p conserved: adot = p/f,        so adot spread = 3")
print()
# period of one full 2pi circuit, both hypotheses
from scipy.integrate import quad
E0 = 0.5*f(0.0)*1.0**2
T_E,_ = quad(lambda x: np.sqrt(f(x)/(2*E0)), 0, 2*np.pi, limit=400)
p0 = f(0.0)*1.0
T_p,_ = quad(lambda x: f(x)/p0, 0, 2*np.pi, limit=400)
print(f"  circuit period, E-conserving (true EOM): {T_E:.10f}")
print(f"  circuit period, p-conserving (spec fmla, oint f da / p): {T_p:.10f}   [= (4pi/3)/p0 = {4*np.pi/3/p0:.10f}]")
# measure it from the integration
ev = sol.sol
tt = np.linspace(0, 12, 2000001)
aa = ev(tt)[0]
i = np.argmax(aa >= 2*np.pi)
print(f"  circuit period, MEASURED from integration:  {tt[i]:.10f}")
