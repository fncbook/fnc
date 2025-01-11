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
[**Demo %s**](#demo-rk-converge)

We solve the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$. We start by getting a reference solution to validate against.

```{code-cell}
from scipy.integrate import solve_ivp
du_dt = lambda t, u: sin((t + u)**2)
tspan = (0.0, 4.0)
u0 = -1.0
sol = solve_ivp(du_dt, tspan, [u0], dense_output=True, atol=1e-13, rtol=1e-13)
u_ref = sol.sol
```

Now we perform a convergence study of our two Runge–Kutta implementations.

```{code-cell}
n = array([int(2 * 10**k) for k in linspace(0, 3, 7)])
err = {"IE2" : [], "RK4" : []}
results = PrettyTable(["n", "IE2 error", "RK4 error"])
for i in range(len(n)):
    t, u = FNC.ie2(du_dt, tspan, u0, n[i])
    err["IE2"].append( abs(u_ref(4)[0] - u[0][-1]) )
    t, u = FNC.rk4(du_dt, tspan, u0, n[i])
    err["RK4"].append( abs(u_ref(4)[0] - u[0][-1]) )
    results.add_row([n[i], err["IE2"][-1], err["RK4"][-1]])

print(results)
```

The amount of computational work at each time step is assumed to be proportional to the number of stages. Let's compare on an apples-to-apples basis by using the number of $f$-evaluations on the horizontal axis.

```{code-cell}
loglog(2 * n, err["IE2"], "-o", label="IE2")
loglog(4 * n, err["RK4"], "-o", label="RK4")
plot(2 * n, 0.5 * err["IE2"][-1] * (n / n[-1])**(-2), "--", label="2nd order")
plot(4 * n, 0.5 * err["RK4"][-1] * (n / n[-1])**(-4), "--", label="4th order")

xlabel("f-evaluations"),  ylabel("inf-norm error")
legend()
title("Convergence of RK methods");
```

The fourth-order variant is more efficient in this problem over a wide range of accuracy.
