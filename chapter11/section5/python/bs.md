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
[**Demo %s**](#demo-boundaries-bs)


```{code-cell}
K = 3;  sigma = 0.06;  r = 0.08;  Smax = 8;
phi = lambda t, x, u, ux, uxx: sigma**2/2 * (x**2 * uxx) + r*x*ux - r*u
ga = lambda u, ux: u
gb = lambda u, ux: ux - 1
```

```{code-cell}
u0 = lambda x: maximum(0, x - K)
x, u = FNC.parabolic(phi, (0, Smax), 80, ga, gb, (0, 15), u0);
```

```{code-cell}
:tags: [hide-input]
from matplotlib.animation import FuncAnimation
fig = figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(0, Smax), ylim=(-0.5, 8))
line, = ax.plot([], [], '-')
ax.set_title("Black–Scholes equation with boundaries")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def snapshot(t):
    line.set_data(x, u(t))
    time_text.set_text(f"t = {t:.2e}")
    return line, time_text

anim = FuncAnimation(fig, animate, frames=linspace(0, 15, 151), blit=True)
anim.save("boundaries-bs.mp4", fps=30)
close()
```

![Black–Scholes equation with boundaries](boundaries-bs.mp4)

Recall that $u$ is the value of the call option, and time runs backward from the strike time. The longer the horizon, the more value the option has due to anticipated growth in the stock price.
