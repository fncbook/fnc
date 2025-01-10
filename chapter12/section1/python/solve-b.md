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
Now we use an initial condition with a larger bump. Note that the scale on the $y$-axis is much different for this solution.

```{code-cell}
from scipy.integrate import solve_ivp
rho_c, rho_m, q_m, ep = (1080, 380, 10000, 0.02)
Q0prime = (
    lambda rho: q_m * 4 * rho_c**2 * (rho_c - rho_m) * rho_m * (rho_m - rho)
    / (rho * (rho_c - 2 * rho_m) + rho_c * rho_m) ** 3
)
x, Dx, Dxx = FNC.diffper(800, [0, 4])
ode = lambda t, rho: -Q0prime(rho) * (Dx @ rho) + ep * (Dxx @ rho)
```

```{code-cell}
rho_init = 400 + 80 * exp(-16 * (x - 3) ** 2)
sol = solve_ivp(ode, [0, 0.5], rho_init, method="Radau", dense_output=True)
```

```{code-cell}
:tags: [hide-input]
for t in linspace(0, 0.5, 6):
    plot(x, sol.sol(t), label=f"t = {t:.1f}")
xlabel("$x$"),  ylabel("car density")
legend(),  title("Traffic jam");
```

```{code-cell}
:tags: [hide-input]
from matplotlib.animation import FuncAnimation
fig, ax = subplots()
curve = ax.plot(x, rho_init)[0]
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$")
ax.set_ylabel("density")
ax.set_ylim(400, 480)
ax.set_title("Traffic jam")
def snapshot(t):
    curve.set_ydata(sol.sol(t))
    time_text.set_text(f"t = {t:.2f}")

anim = FuncAnimation(fig, snapshot, frames=linspace(0, 0.5, 101))
anim.save("traffic-jam.mp4", fps=30)
close()
```

![Traffic jam simulation](traffic-jam.mp4)

In this case the density bump travels backward along the road. It also steepens on the side facing the incoming traffic and decreases much more slowly on the other side. A motorist would experience this as an abrupt increase in density, followed by a much more gradual decrease in density and resulting gradual increase in speed. (You also see some transient, high-frequency oscillations. These are caused by instabilities, as we discuss in simpler situations later in this chapter.)

