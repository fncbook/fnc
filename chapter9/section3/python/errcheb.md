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
[**Demo %s**](#demo-stability-errcheb)

Now we look at the error indicator function $\Phi$ for Chebyshev node sets.

```{code-cell} 
:tags: [hide-input]
x = linspace(-1, 1, 1601)
labels = []
for n in range(10, 60, 10):
    theta = pi * arange(n + 1) / n
    t = -cos(theta)
    Phi = array([prod(xk - t) for xk in x])
    semilogy(x, abs(Phi), ".")
    labels.append(f"degree {n}")

xlabel("$x$"),  ylabel("$|\\Phi(x)|$"),  ylim([1e-18, 1e-2])
legend(labels),  title("Error indicator for Chebyshev nodes");
```

In contrast to the equispaced case, $|\Phi|$ decreases exponentially with $n$ almost uniformly across the interval.
