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
[**Demo %s**](#demo-secant-iqi)

Here we look for a root of $x+\cos(10x)$ that is close to 1.

```{code-cell}
f = lambda x: x + cos(10 * x)
xx = linspace(0.5, 1.5, 400)
fig, ax = subplots()
ax.plot(xx, f(xx), label="function")
ax.grid()
xlabel("$x$"), ylabel("$y$")
fig
```

We choose three values to get the iteration started.

```{code-cell}
x = array([0.8, 1.2, 1])
y = f(x)
ax.plot(x, y, "ko", label="initial points")
ax.legend()
fig
```

If we were using forward interpolation, we would ask for the polynomial interpolant of $y$ as a function of $x$. But that parabola has no real roots.

```{code-cell}
q = poly1d(polyfit(x, y, 2))  # interpolating polynomial
ax.plot(xx, q(xx), "--", label="interpolant")
ax.set_ylim(-0.1, 3), ax.legend()
fig
```

To do inverse interpolation, we swap the roles of $x$ and $y$ in the interpolation.

```{code-cell}
plot(xx, f(xx), label="function")
plot(x, y, "ko", label="initial points")

q = poly1d(polyfit(y, x, 2))  # inverse interpolating polynomial
yy = linspace(-0.1, 2.6, 400)
plot(q(yy), yy, "--", label="inverse interpolant")

grid(), xlabel("$x$"), ylabel("$y$")
legend();
```

We seek the value of $x$ that makes $y$ zero. This means evaluating $q$ at zero.

```{code-cell}
x = hstack([x, q(0)])
y = hstack([y, f(x[-1])])
print("x:", x, "\ny:", y)
```

We repeat the process a few more times.

```{code-cell}
for k in range(6):
    q = poly1d(polyfit(y[-3:], x[-3:], 2))
    x = hstack([x, q(0)])
    y = hstack([y, f(x[-1])])
print(f"final residual is {y[-1]:.2e}")
```

Here is the sequence of errors.

```{code-cell}
from scipy.optimize import root_scalar
r = root_scalar(f, bracket=[0.9, 1]).root
err = x - r
print(err)
```

The error seems to be superlinear, but subquadratic:

```{code-cell}
logerr = log(abs(err))
for i in range(len(err) - 1):
    print(logerr[i+1] / logerr[i])
```
