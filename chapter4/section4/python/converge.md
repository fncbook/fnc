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
[**Demo %s**](#demo-secant-converge)

We check the convergence of the secant method from {numref}`Demo %s <demo-secant-line>`.

```{code-cell}
f = lambda x: x * exp(x) - 2
x = FNC.secant(f, 1, 0.5)
print(x)
```

We don't know the exact root, so we use `root_scalar` to get a substitute.

```{code-cell}
from scipy.optimize import root_scalar
r = root_scalar(f, bracket=[0.5, 1]).root
print(r)
```

Here is the sequence of errors.

```{code-cell}
err = r - x
print(err)
```

It's not easy to see the convergence rate by staring at these numbers. We can use {eq}`superlinear-rate` to try to expose the superlinear convergence rate.

```{code-cell}
logerr = log(abs(err))
for i in range(len(err) - 2):
    print(logerr[i+1] / logerr[i])
```

As expected, this settles in at around 1.618.
