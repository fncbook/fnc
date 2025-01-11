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
[**Demo %s**](#demo-stability-rungefix)

Here again is the function from {numref}`Demo {number} <demo-stability-runge>` that provoked the Runge phenomenon when using equispaced nodes.

```{code-cell} 
f = lambda x: 1 / (x**2 + 16)
```

```{code-cell} 
:tags: [hide-input]
x = linspace(-1, 1, 1601)
labels = []
for k, n in enumerate([4, 10, 16, 40]):
    t = -cos(pi * arange(n + 1) / n)         # Chebyshev nodes
    y = f(t)                                 # interpolation data
    p = FNC.polyinterp(t, y)
    err = abs(f(x) - p(x))
    semilogy(x, err, ".", markersize=2)
    labels.append(f"degree {n}")

xlabel("$x$"),  ylabel("$|f(x)-p(x)|$"),  ylim([1e-20, 1])
legend(labels),  title("Error for Chebyshev interpolants");
```

By degree 16 the error is uniformly within machine epsilon, and, importantly, it stays there as $n$ increases. Note that as predicted by the error indicator function, the error is uniform over the interval at each value of $n$.
