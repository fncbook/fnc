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
[**Demo %s**](#demo-interp-cond)

In {numref}`Demo %s <demo-interpolation-global>` and {numref}`Demo %s <demo-interpolation-pwise>` we saw a big difference between polynomial interpolation and piecewise polynomial interpolation of some arbitrarily chosen data. The same effects can be seen clearly in the cardinal functions, which are closely tied to the condition numbers.

```{code-cell}
from scipy.interpolate import interp1d
n = 18
t = linspace(-1, 1, n + 1)
y = zeros(n + 1)
y[9] = 1.0
p = interp1d(t, y, kind="cubic")

scatter(t, y, label="data")
tt = linspace(-1, 1, 400)
plot(tt, p(tt), label="cardinal function")
title("Cubic spline cardinal function")
legend();
```

The piecewise cubic cardinal function is nowhere greater than one in absolute value. This happens to be true for all the cardinal functions, ensuring a good condition number for any interpolation with these functions. But the story for global polynomials is very different.

```{code-cell}
from numpy.polynomial.polynomial import polyfit, poly1d
p = poly1d(polyfit(t, y, n))
scatter(t, y, label="data")
plot(tt, p(tt), label="cardinal function")
xlabel("$x$")
ylabel("$y$")
title("Polynomial cardinal function")
legend();
```

From the figure we can see that the condition number for polynomial interpolation on these nodes is at least 500.
