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
[**Demo %s**](#demo-interpolation-pwise)

Let us recall the data from {numref}`Demo %s <demo-interpolation-global>`.

```{code-cell}
clf
n = 12
t = linspace(-1, 1, n + 1)
y = t**2 + t + 0.5 * sin(20 * t)
fig, ax = subplots()
scatter(t, y, label="data")
xlabel("$x$"),  ylabel("$y$");
```

Here is an interpolant that is linear between each consecutive pair of nodes, using `plinterp` from {numref}`section-localapprox-pwlin`.

```{code-cell}
from scipy.interpolate import interp1d
tt = linspace(-1, 1, 400)
p = interp1d(t, y, kind="linear")
ax.plot(tt, p(tt), label="piecewise linear")
ax.legend()
fig
```

```{index} ! Julia; Spline1D
```

We may prefer a smoother interpolant that is piecewise cubic:

```{code-cell}
scatter(t, y, label="data")
p = interp1d(t, y, kind="cubic")
tt = linspace(-1, 1, 400)
plot(tt, p(tt), label="cubic spline")
xlabel("$x$"),  ylabel("$y$")
legend();
```
