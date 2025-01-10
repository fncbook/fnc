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
[**Demo %s**](#demo-methodlines-auto)

We set up the semidiscretization and initial condition in $x$ just as before.

```{code-cell}
m = 100
x, Dx, Dxx = FNC.diffper(m, [0, 1])
u0 = exp(-60 * (x - 0.5) ** 2)
```

Now, however, we apply a standard solver using `solve_ivp` to the initial-value problem $\mathbf{u}'=\mathbf{D}_{xx}\mathbf{u}$.

```{code-cell}
from scipy.integrate import solve_ivp
tfinal = 0.05
f = lambda t, u: Dxx @ u
sol = solve_ivp(f, [0, tfinal], u0, method="RK45", dense_output=True)

t = linspace(0, 0.05, 5)
plot(x, sol.sol(t))
xlabel("$x$"),  ylabel("$u(x,t)$")
legend([f"$t={tj:.4g}$" for tj in t])
title("Heat equation by RK45");
```

The solution appears to be correct. But the number of time steps that were selected automatically is surprisingly large, considering how smoothly the solution changes.

```{code-cell}
print(f"RK45 took {len(sol.t) - 1} steps")
```

Now we apply a different solver called `BDF`.

```{code-cell}
sol = solve_ivp(f, [0, tfinal], u0, method="BDF")
print(f"BDF took {len(sol.t) - 1} steps")
```

The number of steps selected was reduced by a factor of 20!
