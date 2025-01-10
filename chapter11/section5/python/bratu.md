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
[**Demo %s**](#demo-boundaries-bratu)


```{code-cell}
phi = lambda t, x, u, ux, uxx: u**2 + uxx
ga = lambda u, ux: u
gb = lambda u, ux: ux
init = lambda x: 400 * x**4 * (1 - x)**2
x, u = FNC.parabolic(phi, (0, 1), 60, ga, gb, (0, 0.1), init);
```

```{code-cell}
:tags: [hide-input]
from matplotlib.animation import FuncAnimation
fig = figure()
ax = fig.add_subplot(autoscale_on=False, xlim=(0, 1), ylim=(0, 10))
line, = ax.plot([], [], '-')
ax.set_title("Heat equation with source")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def snapshot(t):
    line.set_data(x, u(t))
    time_text.set_text(f"t = {t:.2e}")
    return line, time_text

anim = FuncAnimation(fig, snapshot, frames=linspace(0, 0.1, 101), blit=True)
anim.save("boundaries-source.mp4", fps=30)
close()
```

![Heat equation with source](boundaries-source.mp4)
