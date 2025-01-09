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
[**Demo %s**](#demo-nonlinear2d-mems)

All we need to define are $\phi$ from {eq}`nonlinpdepde` for the PDE, and a trivial zero function for the boundary condition.

```{code-cell}
lamb = 1.5
phi = lambda x, y, u, ux, uxx, uy, uyy: uxx + uyy - lamb / (u + 1)**2
g = lambda x, y: 0
```

Here is the solution for $m=15$, $n=8$.

```{code-cell}
u = FNC.elliptic(phi, g, 15, [0, 2.5], 8, [0, 1])

print(f"solution at (2, 0.6) is {u(2, 0.6):.7f}")
```

```{code-cell}
:tags: [hide-input]
x = linspace(0, 2.5, 90)
y = linspace(0, 1, 60)
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)
U = mtx(u)

pcolormesh(X.T, Y.T, U.T, cmap="viridis")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar()
title("Solution of the MEMS equation in 2d");
```

In the absence of an exact solution, how can we be confident that the solution is accurate? First, the Levenberg iteration converged without issuing a warning, so we should feel confident that the discrete equations were solved. We can check the boundary values easily. For example,

```{code-cell}
err = norm(u(x, 0) - g(x, 0), inf)
print(f"max error on bottom edge: {err:.2e}")
```

Assuming that we encoded the PDE correctly, the remaining source error is truncation from the discretization. We can estimate that by refining the grid a bit and seeing how much the numerical solution changes.

```{code-cell}
x_test = linspace(0, 2.5, 6)
y_test = linspace(0, 1, 6)
mtx_test, X_test, Y_test, _, _, _ = FNC.tensorgrid(x_test, y_test)

with printoptions(precision=7, suppress=True):
    print(mtx_test(u))
```

```{code-cell}
u = FNC.elliptic(phi, g, 25, [0, 2.5], 14, [0, 1])
with printoptions(precision=7, suppress=True):
    print(mtx_test(u))
```

The original solution seems to be accurate to about four digits.

