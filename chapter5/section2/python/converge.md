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
[**Demo %s**](#demo-pwlin-converge)

We measure the convergence rate for piecewise linear interpolation of $e^{\sin 7x}$ over $x \in [0,1]$.

```{code-cell}
f = lambda x: exp(sin(7 * x))
x = linspace(0, 1, 10000)  # sample the difference at many points
N = 2 ** arange(3, 11)
err = zeros(N.size)
for i, n in enumerate(N):
    t = linspace(0, 1, n + 1)  # interpolation nodes
    p = FNC.plinterp(t, f(t))
    err[i] = max(abs(f(x) - p(x)))
print(err)
```

As predicted, a factor of 10 in $n$ produces a factor of 100 in the error. In a convergence plot, it is traditional to have $h$ *decrease* from left to right, so we expect a straight line of slope $-2$ on a log-log plot.

```{code-cell}
order2 = 0.1 * (N / N[0]) ** (-2)
loglog(N, err, "-o", label="observed error")
loglog(N, order2, "--", label="2nd order")
xlabel("$n$")
ylabel("$\|f-p\|_\infty$")
legend();
```
