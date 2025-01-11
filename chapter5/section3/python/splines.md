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
[**Demo %s**](#demo-splines-splines)

For illustration, here is a spline interpolant using just a few nodes.

```{code-cell}
f = lambda x: exp(sin(7 * x))

x = linspace(0, 1, 500)
fig, ax = subplots()
ax.plot(x, f(x), label="function")

t = array([0, 0.075, 0.25, 0.55, 0.7, 1])  # nodes
y = f(t)  # values at nodes

xlabel("$x$")
ylabel("$y$")
ax.scatter(t, y, label="nodes")
```

```{code-cell}
S = FNC.spinterp(t, y)
ax.plot(x, S(x), label="spline")
ax.legend()
fig
```

Now we look at the convergence rate as the number of nodes increases.

```{code-cell}
N = floor(2 ** linspace(3, 8, 17)).astype(int)
err = zeros(N.size)
for i, n in enumerate(N):
    t = linspace(0, 1, n + 1)  # interpolation nodes
    p = FNC.spinterp(t, f(t))
    err[i] = max(abs(f(x) - p(x)))
print(err)
```

Since we expect convergence that is $O(h^4)=O(n^{-4})$, we use a log-log graph of error and expect a straight line of slope $-4$.

```{code-cell}
order4 = (N / N[0]) ** (-4)
loglog(N, err, "-o", label="observed error")
loglog(N, order4, "--", label="4th order")
xlabel("$n$")
ylabel("$\|f-S\|_\infty$")
legend();
```
