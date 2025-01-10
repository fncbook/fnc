---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-stiffness-explicit)

The `BDF` solver is good for stiff problems and needs few time steps to solve the Oregonator from {numref}`Demo {number} <demo-stiffness-oregon>`.

```{code-cell}
from scipy.integrate import solve_ivp
tspan = (0, 25)
start = timer()
sol = solve_ivp(ode, tspan, u0, method="BDF")
print(f"stiff solver took {timer() - start:.3f} seconds with {len(sol.t) - 1} time steps")
```

But if we apply {numref}`Function {number} <function-rk23>` to the problem, the step size will be made small enough to cope with the large negative eigenvalue. 

```{code-cell}
start = timer()
t, u = FNC.rk23(ode, tspan, u0, 1e-6)
print(f"rk23 solver took {timer() - start:.3f} seconds with {len(t) - 1} time steps")
```

Starting from the eigenvalues of the Jacobian matrix, we can find an effective $\zeta(t)$ by multiplying with the local time step size. The values of $\zeta(t)$ for each time level are plotted below and color coded by component of the diagonalized system.

```{code-cell}
:tags: [hide-input]
zeta = zeros([len(t)- 1, 3]) + 0j    # complex array
for i in range(len(t) - 1):
    dt = t[i+1] - t[i]
    lamb = eigvals(J(u[:, i]))
    zeta[i] = lamb * dt
plot(real(zeta), imag(zeta), ".")
axis("equal")
xlabel("Re $\\zeta$")
ylabel("Im $\\zeta$")
title("Oregonator stability");
```

Roughly speaking, the $\zeta$ values stay within or close to the RK2 stability region in {numref}`figure-stabreg_bd_rk`. Momentary departures from the region are possible, but time stepping repeatedly in that situation would cause instability. 

