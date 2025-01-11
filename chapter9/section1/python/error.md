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
[**Demo %s**](#demo-polynomial-error)


```{code-cell}
from scipy.interpolate import BarycentricInterpolator as interp
t = array([1, 1.6, 1.9, 2.7, 3])
p = interp(t, log(t))
```

```{code-cell}
from scipy.interpolate import BarycentricInterpolator as interp
t = array([1, 1.6, 1.9, 2.7, 3])
p = interp(t, log(t))
x = linspace(1, 3, 500)
Phi = lambda x: prod([x - ti for ti in t])
plot(x, [Phi(xj) / 5 for xj in x], label="$\\frac{1}{5}|\\Phi(x)|$")
plot(x, abs(log(x) - p(x)), label="$|f(x)-p(x)|$")
plot(t, zeros(t.size), "ko", label="nodes")
xlabel("$x$"),  ylabel("error")
title("Interpolation error and upper bound"),  legend();
```

The error is zero at the nodes, by the definition of interpolation. The error bound, as well as the error itself, has one local maximum between each consecutive pair of nodes.
