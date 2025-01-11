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
[**Demo %s**](#demo-orthogonal-approx)

Let's approximate $e^x$ over the interval $[−1,1]$. We can sample it at, say, 15 points, and find the best-fitting straight line to that data.

```{code-cell}
from numpy.linalg import lstsq
t = linspace(-1, 1, 15)
y = exp(t)
plot(t, y, label="function")

V = [[ti**j for j in range(2)] for ti in t]
c = lstsq(V, y, rcond=None)[0]
print("fit coeffs:", c)

x = linspace(-1, 1, 600)
plot(x, c[1] + c[0] * x, label="fit")
xlabel("x"),  ylabel("value")
legend(),  title("Least squares fit of exp(x)");
```

There's nothing special about 15 points. Choosing more doesn't change the result much.

```{code-cell}
t = linspace(-1, 1, 150)
y = exp(t)
plot(t, y, label="function")

V = [[ti**j for j in range(2)] for ti in t]
c = lstsq(V, y, rcond=None)[0]
print("fit coeffs:", c)

x = linspace(-1, 1, 600)
plot(x, c[1] + c[0] * x, label="fit")
xlabel("x"),  ylabel("value")
legend(),  title("Least squares fit of exp(x)");
```

This situation is unlike interpolation, where the degree of the interpolant increases with the number of nodes. Here, the linear fit is apparently approaching a limit that we may think of as a continuous least-squares fit.

```{code-cell}
n = arange(40, 420, 60)
results = PrettyTable(["n", "intercept", "slope"])
slope = zeros(n.size)
intercept = zeros(n.size)

for k in range(n.size):
    t = linspace(-1, 1, n[k])
    y = exp(t)
    V = [[ti**j for j in range(2)] for ti in t]
    c = lstsq(V, y, rcond=None)[0]
    results.add_row([n[k], c[1], c[0]])

print(results)
```
