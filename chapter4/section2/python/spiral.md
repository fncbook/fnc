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
[**Demo %s**](#demo-fp-spiral)

Let's convert the roots of a quadratic polynomial $f(x)$ to a fixed point problem.

```{code-cell}
f = poly1d([1, -4, 3.5])
r = f.roots
print(r)
```

We define $g(x)=x - f(x)$. 

```{code-cell}
g = lambda x: x - f(x)
```

Intersections of $y=g(x)$ with the line $y=x$ are fixed points of $g$ and thus roots of $f$. (Only one is shown in the chosen plot range.)

```{code-cell}
fig, ax = subplots()
g = lambda x: x - f(x)
xx = linspace(2, 3, 400)
ax.plot(xx, g(xx), label="y=g(x)")
ax.plot(xx, xx, label="y=x")
axis("equal"), legend()
title("Finding a fixed point");
```

If we evaluate $g(2.1)$, we get a value of almost 2.6, so this is not a fixed point.

```{code-cell}
x = 2.1
y = g(x)
print(y)
```

However, $y=g(x)$ is considerably closer to the fixed point at around 2.7 than $x$ is. Suppose then that we adopt $y$ as our new $x$ value. Changing the $x$ coordinate in this way is the same as following a horizontal line over to the graph of $y=x$.

```{code-cell}
ax.plot([x, y], [y, y], "r:", label="")
fig
```

Now we can compute a new value for $y$. We leave $x$ alone here, so we travel along a vertical line to the graph of $g$.

```{code-cell}
x = y
y = g(x)
print("y:", y)
ax.plot([x, x], [x, y], "k:")
fig
```

You see that we are in a position to repeat these steps as often as we like. Let's apply them a few times and see the result.

```{code-cell}
for k in range(5):
    ax.plot([x, y], [y, y], "r:")
    x = y       # y --> new x
    y = g(x)    # g(x) --> new y
    ax.plot([x, x], [x, y], "k:")  
fig
```

The process spirals in beautifully toward the fixed point we seek. Our last estimate has almost 4 accurate digits.

```{code-cell} 
print(abs(y - max(r)) / max(r))
```

Now let's try to find the other fixed point $\approx 1.29$ in the same way. We'll use 1.3 as a starting approximation.

```{code-cell}
xx = linspace(1, 2, 400)
fig, ax = subplots()
ax.plot(xx, g(xx), label="y=g(x)")
ax.plot(xx, xx, label="y=x")
ax.set_aspect(1.0)
ax.legend()

x = 1.3
y = g(x)
for k in range(5):
    ax.plot([x, y], [y, y], "r:")
    x = y
    y = g(x)
    ax.plot([x, x], [x, y], "k:")
ylim(1, 2.5)
title("No convergence");
```

This time, the iteration is pushing us _away from_ the correct answer.
