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
[**Demo %s**](#demo-newton-line)


Suppose we want to find a root of this function:

```{code-cell}
f = lambda x: x * exp(x) - 2
xx = linspace(0, 1.5, 400)

fig, ax = subplots()
ax.plot(xx, f(xx), label="function")
ax.grid()
ax.set_xlabel("$x$")
ax.set_ylabel("$y$");
```

From the graph, it is clear that there is a root near $x=1$. So we call that our initial guess, $x_1$.

```{code-cell}
x1 = 1
y1 = f(x1)
ax.plot(x1, y1, "ko", label="initial point")
ax.legend()
fig
```

Next, we can compute the tangent line at the point $\bigl(x_1,f(x_1)\bigr)$, using the derivative.

```{code-cell}
df_dx = lambda x: exp(x) * (x + 1)
slope1 = df_dx(x1)
tangent1 = lambda x: y1 + slope1 * (x - x1)

ax.plot(xx, tangent1(xx), "--", label="tangent line")
ax.set_ylim(-2, 4)
ax.legend()
fig
```

In lieu of finding the root of $f$ itself, we settle for finding the root of the tangent line approximation, which is trivial. Call this $x_2$, our next approximation to the root.

```{code-cell}
x2 = x1 - y1 / slope1
ax.plot(x2, 0, "ko", label="tangent root")
ax.legend()
fig
```

```{code-cell}
y2 = f(x2)
print(y2)
```

The residual (i.e., value of $f$) is smaller than before, but not zero. So we repeat the process with a new tangent line based on the latest point on the curve.

```{code-cell}
xx = linspace(0.83, 0.88, 200)

plot(xx, f(xx))
plot(x2, y2, "ko")
grid(), xlabel("$x$"), ylabel("$y$")

slope2 = df_dx(x2)
tangent2 = lambda x: y2 + slope2 * (x - x2)
plot(xx, tangent2(xx), "--")
x3 = x2 - y2 / slope2
plot(x3, 0, "ko")
title("Second iteration");
```

```{code-cell}
y3 = f(x3)
print(y3)
```

Judging by the residual, we appear to be getting closer to the true root each time.
