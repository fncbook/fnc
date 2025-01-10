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
[**Demo %s**](#demo-wave-boundaries)


```{code-cell}
m = 200
x, Dx, Dxx = FNC.diffmat2(m, [-1, 1])
```

The boundary values of $u$ are given to be zero, so they are not unknowns in the ODEs. Instead they are added or removed as necessary.

```{code-cell}
chop = lambda u: u[1:-1]
extend = lambda v: hstack([0, v, 0])
```

The following function computes the time derivative of the system at interior points.

```{code-cell}
def dw_dt(t, w):
    u = extend(w[:m-1])
    z = w[m-1:]
    du_dt = Dx @ z
    dz_dt = c**2 * (Dx @ u)
    return hstack([chop(du_dt), dz_dt])
```

Our initial condition is a single hump for $u$.

```{code-cell}
u_init = exp(-100 * x**2)
z_init = -u_init
w_init = hstack([chop(u_init), z_init])
```

Because the wave equation is hyperbolic, we can use a nonstiff explicit solver.

```{code-cell}
from scipy.integrate import solve_ivp
c = 2
sol = solve_ivp(dw_dt, (0, 2), w_init, dense_output=True)
u = lambda t: extend(sol.sol(t)[:m-1])   # extract the u component
```

We plot the results for the original $u$ variable only. Its interior values are at indices `1:m-1` of the composite $\mathbf{w}$ variable.

```{code-cell}
t = linspace(0, 2, 80)
U = [u(tj) for tj in t]
contour(x, t, U, levels=24, cmap="RdBu", vmin=-1, vmax=1)
xlabel("$x$"),  ylabel("$t$")
title("Wave equation with boundaries");
```

```{code-cell}
:tags: [hide-input]
from matplotlib.animation import FuncAnimation
fig, ax = subplots()
curve = ax.plot(x, u_init)[0]
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$")
ax.set_ylabel("$u(x,t)$")
ax.set_ylim(-1.05, 1.05)
ax.set_title("Wave equation with boundaries")

def snapshot(t):
    curve.set_ydata(u(t))
    time_text.set_text(f"t = {t:.2f}")

anim = FuncAnimation(fig, snapshot, frames=linspace(0, 2, 161))
anim.save("wave-boundaries.mp4", fps=30)
close()
```

![Wave equation with boundaries](wave-boundaries.mp4)

The original hump breaks into two pieces of different amplitudes, each traveling with speed $c=2$. They pass through one another without interference. When a hump encounters a boundary, it is perfectly reflected, but with inverted shape. At time $t=2$, the solution looks just like the initial condition.

