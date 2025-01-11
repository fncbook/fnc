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
[**Demo %s**](#demo-interpolation-global)

Here are some points that we could consider to be observations of an unknown function on $[-1,1]$.

```{code-cell}
n = 5
t = linspace(-1, 1, n + 1)
y = t**2 + t + 0.05 * sin(20 * t)
fig, ax = subplots()
plot(t, y, "o", label="data")
xlabel("$x$"),  ylabel("$y$");
```

```{index} ! Julia; fit
```

The polynomial interpolant, as computed using `fit`, looks very sensible. It's the kind of function you'd take home to meet your parents.

```{code-cell}
p = poly1d(polyfit(t, y, n))  # interpolating polynomial
tt = linspace(-1, 1, 400)
ax.plot(tt, p(tt), label="interpolant")
ax.legend()
fig
```

But now consider a different set of points generated in almost exactly the same way.

```{code-cell}
n = 18
t = linspace(-1, 1, n + 1)
y = t**2 + t + 0.05 * sin(20 * t)
fig, ax = subplots()
plot(t, y, "o", label="data")
xlabel("$x$"),  ylabel("$y$");
```

The points themselves are unremarkable. But take a look at what happens to the polynomial interpolant.

```{code-cell}
p = poly1d(polyfit(t, y, n))
ax.plot(tt, p(tt), label="interpolant")
ax.legend()
fig
```

Surely there must be functions that are more intuitively representative of those points!
