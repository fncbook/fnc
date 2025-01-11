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
[**Demo %s**](#demo-pwlin-usage)

We generate a piecewise linear interpolant of $f(x)=e^{\sin 7x}$.

```{code-cell}
f = lambda x: exp(sin(7 * x))
x = linspace(0, 1, 400)
fig, ax = subplots()
plot(x, f(x), label="function")
xlabel("$x$")
ylabel("$f(x)$");
```

First we sample the function to create the data.

```{code-cell}
t = array([0, 0.075, 0.25, 0.55, 0.7, 1])  # nodes
y = f(t)  # function values

ax.plot(t, y, "o", label="nodes")
ax.legend()
fig
```

Now we create a callable function that will evaluate the piecewise linear interpolant at any $x$, and then plot it.

```{code-cell}
p = FNC.plinterp(t, y)
ax.plot(x, p(x), label="interpolant")
ax.legend()
fig
```
