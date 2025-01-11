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
[**Demo %s**](#demo-rootproblem-bessel)


```{code-cell}
import scipy.special as special
def J3(x):
    return special.jv(3.0, x)

xx = linspace(0, 20, 500)
fig, ax = subplots()
ax.plot(xx, J3(xx))
ax.grid()
xlabel("$x$"), ylabel("$J_3(x)$")
title("Bessel function");
```
From the graph we see roots near 6, 10, 13, 16, and 19. We use `root_scalar` from the `scipy.optimize` package to find these roots accurately.

```{code-cell}
from scipy.optimize import root_scalar

omega = []
for guess in [6.0, 10.0, 13.0, 16.0, 19.0]:
    s = root_scalar(J3, bracket=[guess - 0.5, guess + 0.5]).root
    omega.append(s)

results = PrettyTable()
results.add_column("root estimate", omega)
results.add_column("function value", [J3(ω) for ω in omega])
print(results)
```

```{code-cell}
ax.scatter(omega, J3(omega))
ax.set_title("Bessel function roots")
fig
```

If instead we seek values at which $J_3(x)=0.2$, then we must find roots of the function $J_3(x)-0.2$.

```{code-cell}
omega = []
for guess in [3., 6., 10., 13.]:
    f = lambda x: J3(x) - 0.2
    s = root_scalar(f, x0=guess).root
    omega.append(s)

ax.scatter(omega, J3(omega))
fig
```
