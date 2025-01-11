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
[**Demo %s**](#demo-secant-line)


```{code-cell}
f = lambda x: x * exp(x) - 2
xx = linspace(0.25, 1.25, 400)

fig, ax = subplots()
ax.plot(xx, f(xx), label="function")
ax.set_xlabel("$x$")
ax.set_ylabel("$f(x)$")
ax.grid();
```

From the graph, it's clear that there is a root near $x=1$. To be more precise, there is a root in the interval $[0.5,1]$. So let us take the endpoints of that interval as _two_ initial approximations.

```{code-cell}
x1 = 1
y1 = f(x1)
x2 = 0.5
y2 = f(x2)
ax.plot([x1, x2], [y1, y2], "ko", label="initial points")
ax.legend()
fig
```

Instead of constructing the tangent line by evaluating the derivative, we can construct a linear model function by drawing the line between the two points $\bigl(x_1,f(x_1)\bigr)$ and $\bigl(x_2,f(x_2)\bigr)$. This is called a _secant line_.

```{code-cell}
slope2 = (y2 - y1) / (x2 - x1)
secant2 = lambda x: y2 + slope2 * (x - x2)
ax.plot(xx, secant2(xx), "--", label="secant line")
ax.legend()
fig
```

As before, the next root estimate in the iteration is the root of this linear model.

```{code-cell}
x3 = x2 - y2 / slope2
ax.plot(x3, 0, "o", label="root of secant")
y3 = f(x3)
print(y3)
ax.legend()
fig
```

For the next linear model, we use the line through the two most recent points. The next iterate is the root of that secant line, and so on.

```{code-cell}
slope3 = (y3 - y2) / (x3 - x2)
x4 = x3 - y3 / slope3
print(f(x4))
```
